"""w06: copy the files listed in bundle_manifest.txt into a target folder laid out as the repository's tic214311908/,
rewriting the working-layout paths that the analysis scripts and READMEs hard-code.

Usage (from the working layout only): python w06_copy_bundle.py <target>/tic214311908
  * "writeup/<path>" -> <target>/<path>; "<folder>/<path>" -> <target>/analysis/<folder>/<path>.
  * In every copied text file (.py .sh .md .json .csv .txt .ecsv) the absolute working paths are rewritten:
      <analysis root>/writeup  -> tic214311908
      <analysis root>          -> tic214311908/analysis
      <scratch folder>         -> <local>      (other local folders, e.g. the TESS PRF files)
    so that scripts run from the repository root (tess-search/) find the committed files, and no absolute local path
    remains (the READMEs' 'Folder:' lines then name the repository folder).
  * Afterwards the copy is searched for absolute /tmp paths and for the analysis root; the script exits with an
    error if any trace is left.
The working paths are computed from this file's location, not written into it. Nothing outside <target> is changed.
The target must not exist yet (no overwriting)."""
import os, re, shutil, sys

W = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))     # writeup/
ROOT = os.path.dirname(W)                                            # the analysis root (holds phot/, ..., writeup/)
SCRATCH = os.path.dirname(ROOT)
TEXT = {'.py', '.sh', '.md', '.json', '.csv', '.txt', '.ecsv'}


def rewrite(s):
    s = s.replace(ROOT + '/writeup', 'tic214311908')
    s = s.replace(ROOT, 'tic214311908/analysis')
    s = s.replace(SCRATCH, '<local>')
    return s


def main(dst):
    if not os.path.isdir(os.path.join(ROOT, 'clean')):
        sys.exit('run this from the working layout (the analysis root must hold clean/)')
    if os.path.exists(dst):
        sys.exit(f'{dst} exists; refusing to overwrite')
    man = [l.strip() for l in open(os.path.join(W, 'bundle_manifest.txt')) if l.strip() and not l.startswith('#')]
    n = 0
    for rel in man:
        src = os.path.join(ROOT, rel)
        out = os.path.join(dst, rel[len('writeup/'):]) if rel.startswith('writeup/') else os.path.join(dst, 'analysis', rel)
        os.makedirs(os.path.dirname(out), exist_ok=True)
        if os.path.splitext(rel)[1].lower() in TEXT:
            open(out, 'w', encoding='utf-8').write(rewrite(open(src, encoding='utf-8').read()))
        else:
            shutil.copy2(src, out)
        n += 1
    left = []
    for dp, _, fns in os.walk(dst):
        for fn in fns:
            p = os.path.join(dp, fn)
            if os.path.splitext(fn)[1].lower() not in TEXT: continue
            s = open(p, encoding='utf-8', errors='replace').read()
            if re.search(r'/tmp/[A-Za-z0-9_.-]', s) or ROOT in s or SCRATCH in s: left.append(os.path.relpath(p, dst))
    print(n, 'files copied to', dst)
    print('files with a local path left:', len(left))
    for p in left: print('  ', p)
    if left: sys.exit(1)


if __name__ == '__main__':
    main(sys.argv[1])
