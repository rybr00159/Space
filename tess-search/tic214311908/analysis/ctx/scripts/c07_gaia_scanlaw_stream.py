"""Stream the Gaia full-mission commanded scanning law (2.5 GB zip from cosmos.esa.int, never written to disk) and keep
every 10-s sample whose FOV1 or FOV2 centre lies within 1.5 deg of TIC 214311908. Reads the zip from stdin:

  curl -sS https://anonftp.cosmos.esa.int/pub/GAIA_PUBLIC_DATA/Gaia_General/GaiaScanningLaw/FullGaiaMissionScanningLaw.zip \
     | python scripts/c07_gaia_scanlaw_stream.py out/gaia_scanlaw_near.csv

The zip uses data descriptors (sizes unknown in the local headers), so it is parsed by hand: local header ->
self-terminating deflate stream -> (gzip) -> CSV.
"""
import sys, struct, zlib, io
import numpy as np
import pandas as pd

RA0, DEC0 = 260.969674, -51.323821   # M3V (Gaia DR3, J2016); the Mira is 18" away, irrelevant at this scale
RAD = 1.5
out_path = sys.argv[1]
src = sys.stdin.buffer
buf = bytearray()


def read(n):
    global buf
    while len(buf) < n:
        chunk = src.read(1 << 20)
        if not chunk:
            break
        buf += chunk
    out = bytes(buf[:n]); del buf[:n]; return out


def entries():
    while True:
        hdr = read(30)
        if len(hdr) < 30 or struct.unpack("<I", hdr[:4])[0] != 0x04034b50:
            return
        _, ver, flag, comp, mt, md, crc, cs, us, nl, el = struct.unpack("<IHHHHHIIIHH", hdr)
        name = read(nl).decode(); read(el)
        yield name, flag, comp, cs


def deflate_stream(comp):
    """Yield decompressed chunks of the current entry; leaves the stream positioned after the compressed data."""
    global buf
    if comp == 0:
        raise RuntimeError("stored entries with data descriptor not supported")
    d = zlib.decompressobj(-15)
    while True:
        if not buf:
            chunk = src.read(1 << 20)
            if not chunk:
                return
            buf += chunk
        data = bytes(buf); buf = bytearray()
        out = d.decompress(data)
        if out:
            yield out
        if d.eof:
            buf = bytearray(d.unused_data); return


def skip_descriptor(flag):
    if flag & 8:
        sig = read(4)
        if struct.unpack("<I", sig)[0] == 0x08074b50:
            read(12)    # crc, csize, usize (32-bit); zip64 would be 20 -> handled below
        else:
            read(8)
        # zip64 descriptor check: if next 4 bytes are not a header signature, consume 8 more
        global buf
        while len(buf) < 4:
            chunk = src.read(1 << 20)
            if not chunk:
                break
            buf += chunk
        if bytes(buf[:4]) not in (b"PK\x03\x04", b"PK\x01\x02"):
            read(8)


cr0 = np.radians(RA0); cd0 = np.radians(DEC0)
v0 = np.array([np.cos(cd0) * np.cos(cr0), np.cos(cd0) * np.sin(cr0), np.sin(cd0)])
cosmax = np.cos(np.radians(RAD))
fout = open(out_path, "w"); nrows = 0; nkeep = 0; header = None


def process_lines(lines):
    global nrows, nkeep, header
    lines = [l for l in lines if not l.startswith("#")]
    if not lines:
        return
    if header is None:
        header = lines[0]; fout.write(header + "\n"); lines = lines[1:]
        print("header:", header, flush=True)
    cols = header.split(",")
    i_r1, i_d1, i_r2, i_d2 = cols.index("ra_fov1"), cols.index("dec_fov1"), cols.index("ra_fov2"), cols.index("dec_fov2")
    if not lines:
        return
    arr = pd.read_csv(io.StringIO("\n".join(lines)), header=None, usecols=[i_r1, i_d1, i_r2, i_d2]).to_numpy(float)
    keep = np.zeros(len(arr), bool)
    for a, b in ((0, 1), (2, 3)):
        r = np.radians(arr[:, a]); d = np.radians(arr[:, b])
        cosd = np.cos(d) * np.cos(r) * v0[0] + np.cos(d) * np.sin(r) * v0[1] + np.sin(d) * v0[2]
        keep |= cosd > cosmax
    nrows += len(arr)
    if keep.any():
        for i in np.where(keep)[0]:
            fout.write(lines[i] + "\n")
        nkeep += int(keep.sum())


for name, flag, comp, cs in entries():
    print("entry", name, "flag", flag, "comp", comp, flush=True)
    if not name.endswith(".csv.gz"):
        for _ in deflate_stream(comp):
            pass
        skip_descriptor(flag); continue
    gz = zlib.decompressobj(16 + zlib.MAX_WBITS); tail = ""; nchunk = 0
    for chunk in deflate_stream(comp):
        text = gz.decompress(chunk)
        while gz.eof and gz.unused_data:     # concatenated gzip members
            rest = gz.unused_data; gz = zlib.decompressobj(16 + zlib.MAX_WBITS); text += gz.decompress(rest)
        if not text:
            continue
        s = tail + text.decode(); lines = s.split("\n"); tail = lines.pop()
        process_lines([l for l in lines if l])
        nchunk += 1
        if nchunk % 200 == 0:
            print(f"rows {nrows} kept {nkeep}", flush=True); fout.flush()
    if tail:
        process_lines([tail])
    skip_descriptor(flag)
fout.close()
print("done rows", nrows, "kept", nkeep, flush=True)
