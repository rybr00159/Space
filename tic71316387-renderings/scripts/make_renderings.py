"""Artist's-impression renderings of the TIC 71316387 system (procedural, numpy only).

Writes four PNGs to ../renderings/. Sizes, temperatures, orbits and brightness ratios are the adopted values of the
campaign report (tess-search/exofop_leads/TIC71316387_campaign/README.md, 2026-10-15). Everything that has not been
measured is invented for illustration: starspot patterns, viewing angles, the orientation of both orbits, the
companion's position angle, and the background sky.

Requires numpy, scipy, pillow and matplotlib (only to locate the Inter font).  Usage:  python make_renderings.py
"""
import os

import numpy as np
from scipy.ndimage import gaussian_filter, zoom
from scipy.special import j1
from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'renderings')

# ---- adopted physical values (campaign report) ------------------------------------------------------------
T_G, R_G, M_G = 5880, 1.07, 1.03          # G star: Teff (K), radius (Rsun), mass (Msun)
T_1, R_1, M_1 = 3550, 0.48, 0.45          # eclipsing primary (M2V)
T_2, R_2, M_2 = 3160, 0.21, 0.20          # eclipsing secondary (M4V)
A_IN = 4.35                               # inner orbit, Rsun (0.0202 au); P = 1.30367 d
DKS = 2.6                                 # contrast of the pair against the G star, Ks band (mag)


# ---- fonts ------------------------------------------------------------------------------------------------
_FONTS, _CMAPS = {}, {}


def font(size, weight='regular', family='Inter'):
    key = (size, weight, family)
    if key not in _FONTS:
        from matplotlib import font_manager as fm
        f = None
        for fam in (family, 'DejaVu Sans'):
            try:
                path = fm.findfont(fm.FontProperties(family=fam, weight=weight), fallback_to_default=False)
                f = ImageFont.truetype(path, size)
                break
            except Exception:
                continue
        _FONTS[key] = f or ImageFont.load_default()
    return _FONTS[key]


def _has(f, ch):
    path = getattr(f, 'path', None)
    if path is None:
        return True
    if path not in _CMAPS:
        from fontTools.ttLib import TTFont
        _CMAPS[path] = set(TTFont(path, lazy=True).getBestCmap())
    return ord(ch) in _CMAPS[path]


def _runs(text, f):
    """Split text into runs drawable with f, or with DejaVu Sans where f lacks the glyph (e.g. the Sun symbol)."""
    fb = font(f.size, family='DejaVu Sans')
    runs = []
    for ch in text:
        ff = f if _has(f, ch) else fb
        if runs and runs[-1][1] is ff:
            runs[-1][0] += ch
        else:
            runs.append([ch, ff])
    return runs


def text_width(d, text, f):
    return sum(d.textlength(t, font=ff) for t, ff in _runs(text, f))


def draw_text(d, xy, text, f, fill, shadow=True):
    """Draw text with its top-left at xy (baseline-aligned runs, glyph fallback, soft shadow)."""
    x, y = xy
    base = y + f.getmetrics()[0]
    for t, ff in _runs(text, f):
        if shadow:
            d.text((x + 1.5, base + 1.5), t, font=ff, fill=(0, 0, 0, 170), anchor='ls')
        d.text((x, base), t, font=ff, fill=fill, anchor='ls')
        x += d.textlength(t, font=ff)


# ---- colour of a blackbody, as an RGB chromaticity ------------------------------------------------------
def _g(l, mu, s1, s2):
    s = np.where(l < mu, s1, s2)
    return np.exp(-0.5 * ((l - mu) / s) ** 2)


def bb_rgb(T, sat=1.0):
    """Linear-sRGB colour of a blackbody (CIE 1931 fit of Wyman et al. 2013), max channel = 1."""
    lam = np.linspace(380, 780, 401)
    B = lam ** -5 / np.expm1(1.4388e7 / (lam * T))
    xb = 1.056 * _g(lam, 599.8, 37.9, 31.0) + 0.362 * _g(lam, 442.0, 16.0, 26.7) - 0.065 * _g(lam, 501.1, 20.4, 26.2)
    yb = 0.821 * _g(lam, 568.8, 46.9, 40.5) + 0.286 * _g(lam, 530.9, 16.3, 31.1)
    zb = 1.217 * _g(lam, 437.0, 11.8, 36.0) + 0.681 * _g(lam, 459.0, 26.0, 13.8)
    xyz = np.array([(B * xb).sum(), (B * yb).sum(), (B * zb).sum()])
    m = np.array([[3.2406, -1.5372, -0.4986], [-0.9689, 1.8758, 0.0415], [0.0557, -0.2040, 1.0570]])
    rgb = np.clip(m @ xyz, 0, None)
    rgb /= rgb.max()
    if sat != 1.0:
        lum = rgb @ np.array([0.2126, 0.7152, 0.0722])
        rgb = np.clip(lum + sat * (rgb - lum), 0, None)
        rgb /= rgb.max()
    return rgb


# ---- value noise ------------------------------------------------------------------------------------------
def _hash(ix, iy, iz):
    h = np.sin(ix * 127.1 + iy * 311.7 + iz * 74.7) * 43758.5453
    return h - np.floor(h)


def vnoise(p):
    pi = np.floor(p)
    f = p - pi
    w = f * f * (3 - 2 * f)
    ix, iy, iz = pi[..., 0], pi[..., 1], pi[..., 2]
    wx, wy, wz = w[..., 0], w[..., 1], w[..., 2]
    c = {(a, b, d): _hash(ix + a, iy + b, iz + d) for a in (0, 1) for b in (0, 1) for d in (0, 1)}
    x00 = c[0, 0, 0] + wx * (c[1, 0, 0] - c[0, 0, 0])
    x10 = c[0, 1, 0] + wx * (c[1, 1, 0] - c[0, 1, 0])
    x01 = c[0, 0, 1] + wx * (c[1, 0, 1] - c[0, 0, 1])
    x11 = c[0, 1, 1] + wx * (c[1, 1, 1] - c[0, 1, 1])
    y0 = x00 + wy * (x10 - x00)
    y1 = x01 + wy * (x11 - x01)
    return y0 + wz * (y1 - y0)


def fbm(p, octaves=5, gain=0.5):
    tot, amp, norm = 0.0, 1.0, 0.0
    for o in range(octaves):
        tot = tot + amp * vnoise(p)
        norm += amp
        p = p * 2.03 + 17.17
        amp *= gain
    return tot / norm


def smoothstep(e0, e1, x):
    t = np.clip((x - e0) / (e1 - e0), 0, 1)
    return t * t * (3 - 2 * t)


# ---- canvas helpers ---------------------------------------------------------------------------------------
def rot(tilt, spin):
    """Rotation from view frame (x right, y up, z to viewer) to the star frame."""
    a, b = np.radians(tilt), np.radians(spin)
    rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    return ry @ rx


def sph(lon, lat):
    lo, la = np.radians(lon), np.radians(lat)
    return np.array([np.cos(la) * np.sin(lo), np.sin(la), np.cos(la) * np.cos(lo)])


def window(img, cx, cy, r):
    H, W = img.shape[:2]
    x0, x1 = max(0, int(cx - r)), min(W, int(cx + r) + 2)
    y0, y1 = max(0, int(cy - r)), min(H, int(cy + r) + 2)
    yy, xx = np.mgrid[y0:y1, x0:x1].astype(float) + 0.5
    return (slice(y0, y1), slice(x0, x1)), xx, yy


def add_glow(img, cx, cy, R, rgb, amp, scale, reach=8.0):
    sl, xx, yy = window(img, cx, cy, R + reach * scale)
    d = np.clip(np.hypot(xx - cx, yy - cy) - R, 0, None)
    img[sl] += amp * np.exp(-d / scale)[..., None] * rgb


_LUT = {}


def colour_lut(sat):
    if sat not in _LUT:
        temps = np.arange(1500, 9001, 25.0)
        _LUT[sat] = (temps, np.array([bb_rgb(t, sat) for t in temps]))
    return _LUT[sat]


def render_star(img, cx, cy, R, T, *, inten=1.0, sat=1.7, u=(0.74, 0.08), spots=(), seed=0.0, tilt=15.0,
                spin=0.0, irr=None, gran=0.22, mottle=0.11):
    """Composite a self-luminous star onto img (linear radiance): limb darkening with limb reddening, convective
    texture, soft starspot groups with faculae, colour from the local temperature, anti-aliased edge."""
    sl, xx, yy = window(img, cx, cy, R + 2)
    nx, ny = (xx - cx) / R, -(yy - cy) / R
    r = np.hypot(nx, ny)
    cov = np.clip((1 - r) * R + 0.5, 0, 1)
    nz = np.sqrt(np.clip(1 - r * r, 0, 1))
    n = np.stack([nx, ny, nz], -1)
    p = n @ rot(tilt, spin).T
    mu = nz
    limb = 1 - u[0] * (1 - mu) - u[1] * (1 - mu) ** 2
    # convection: large mottling plus ridged turbulence (dark lanes between bright cells)
    turb, amp, norm = 0.0, 1.0, 0.0
    for o, fr in enumerate((16, 32, 64, 128)):
        turb = turb + amp * np.abs(2 * vnoise(p * fr + seed + 7.1 * o) - 1)
        norm += amp
        amp *= 0.6
    turb = turb / norm
    tex = 1 + mottle * (2 * fbm(p * 3.2 + seed + 2.0, 4) - 1) + gran * (turb - 0.45)
    # starspot groups: soft dark umbrae in grey penumbrae, irregular edges, bright faculae around them
    spot_i = np.ones_like(mu)
    spot_t = np.ones_like(mu)
    fac = np.zeros_like(mu)
    edge = 2 * fbm(p * 9 + seed + 11, 4) - 1
    for lon, lat, size, depth in spots:
        ang = np.degrees(np.arccos(np.clip(p @ sph(lon, lat), -1, 1)))
        a_n = ang * (1 + 0.45 * edge)
        pen = smoothstep(size, 0.55 * size, a_n)
        umb = smoothstep(0.5 * size, 0.15 * size, a_n)
        spot_i = np.minimum(spot_i, 1 - 0.30 * pen - (0.70 - depth) * umb)
        spot_t = np.minimum(spot_t, 1 - 0.05 * pen - 0.12 * umb)
        fac = np.maximum(fac, smoothstep(2.4 * size, 1.2 * size, ang) * (1 - pen))
    fac *= smoothstep(0.85, 0.2, mu) * (turb > 0.4)               # faculae show up towards the limb
    tex *= 1 + 0.12 * fac
    if irr is not None:
        d, st = irr
        tex *= 1 + st * np.clip(n @ (np.asarray(d) / np.linalg.norm(d)), 0, None) ** 1.5
    inten_map = limb * tex * spot_i
    t_loc = T * np.clip(limb * tex, 0.05, None) ** 0.32 * spot_t
    temps, lut = colour_lut(sat)
    col = np.stack([np.interp(t_loc, temps, lut[:, k]) for k in range(3)], -1)
    rad = inten * inten_map[..., None] * col
    img[sl] = img[sl] * (1 - cov[..., None]) + rad * cov[..., None]


def star_glow(img, cx, cy, R, T, inten, sat=1.7):
    """Chromosphere-like rim glow plus a wide halo, in the star's limb colour."""
    c = bb_rgb(T * 0.9, sat)
    add_glow(img, cx, cy, R, c, 0.55 * inten, max(2.0, 0.018 * R))
    add_glow(img, cx, cy, R, c, 0.16 * inten, 0.12 * R)
    add_glow(img, cx, cy, R, c, 0.05 * inten, 0.5 * R, reach=6)


def point_star(img, cx, cy, rgb, amp, core=2.5, halo=(1.2, 22), wide=(0.25, 110), spikes=(0.9, 520), angles=(45, 135)):
    """A very bright unresolved star: core, halo, wide glow and diffraction spikes."""
    sl, xx, yy = window(img, cx, cy, max(wide[1] * 7, spikes[1]))
    dx, dy = xx - cx, yy - cy
    d = np.hypot(dx, dy)
    f = amp * np.exp(-0.5 * (d / core) ** 2) + halo[0] * np.exp(-d / halo[1]) + wide[0] * np.exp(-d / wide[1])
    for a in angles:
        t = np.radians(a)
        along = np.abs(dx * np.cos(t) + dy * np.sin(t))
        perp = np.abs(-dx * np.sin(t) + dy * np.cos(t))
        f += spikes[0] * np.exp(-0.5 * (perp / 1.3) ** 2) * np.clip(1 - along / spikes[1], 0, 1) ** 2 / (1 + along / 40)
    img[sl] += f[..., None] * rgb


def starfield(img, n, seed, density=None, maxflux=1.2):
    rng = np.random.default_rng(seed)
    H, W = img.shape[:2]
    if density is None:
        xs, ys = rng.uniform(0, W, n), rng.uniform(0, H, n)
    else:                                         # rejection-sample a density map
        xs, ys = [], []
        while len(xs) < n:
            x, y = rng.uniform(0, W, 4 * n), rng.uniform(0, H, 4 * n)
            k = rng.uniform(0, 1, x.size) < density[y.astype(int), x.astype(int)]
            xs += list(x[k]); ys += list(y[k])
        xs, ys = np.array(xs[:n]), np.array(ys[:n])
    flux = maxflux * 10 ** (-0.4 * rng.exponential(2.2, n))
    temps = np.clip(rng.lognormal(np.log(5200), 0.35, n), 2800, 14000)
    for x, y, fl, T in zip(xs, ys, flux, temps):
        s = 0.5 + 0.9 * fl ** 0.5
        sl, xx, yy = window(img, x, y, 4 * s + 1)
        img[sl] += (fl * np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / s ** 2))[..., None] * bb_rgb(T, 1.25)


def milky_way(img, p0, p1, width, amp, seed):
    """Faint diffuse band with dust lanes, computed at quarter resolution."""
    H, W = img.shape[:2]
    h, w = H // 4, W // 4
    yy, xx = np.mgrid[0:h, 0:w].astype(float) * 4
    (x0, y0), (x1, y1) = p0, p1
    L = np.hypot(x1 - x0, y1 - y0)
    d = ((x1 - x0) * (y0 - yy) - (x0 - xx) * (y1 - y0)) / L
    p = np.stack([xx / 420, yy / 420, np.full_like(xx, seed)], -1)
    band = np.exp(-0.5 * (d / width) ** 2) * (0.45 + 0.9 * fbm(p, 6))
    lanes = smoothstep(0.52, 0.72, fbm(p * 2.2 + 3.3, 5)) * np.exp(-0.5 * (d / (0.6 * width)) ** 2)
    band *= 1 - 0.75 * lanes
    warm, cool = np.array([1.0, 0.82, 0.62]), np.array([0.55, 0.68, 1.0])
    core = np.exp(-0.5 * (d / (0.5 * width)) ** 2)[..., None]
    col = warm * core + cool * (1 - core)
    lay = amp * band[..., None] * col
    img += zoom(lay, (H / h, W / w, 1), order=1)[:H, :W]
    return zoom(np.exp(-0.5 * (d / width) ** 2), (H / h, W / w), order=1)[:H, :W]


def bloom(img, levels=((4, 2.0, 0.10), (8, 4.0, 0.08), (16, 6.0, 0.06))):
    H, W = img.shape[:2]
    out = img.copy()
    for f, s, w in levels:
        h, wd = H // f, W // f
        small = img[:h * f, :wd * f].reshape(h, f, wd, f, 3).mean(axis=(1, 3))
        small = gaussian_filter(small, (s, s, 0))
        out += w * zoom(small, (H / h, W / wd, 1), order=1)[:H, :W]
    return out


def tonemap(img, exposure=1.0, white=6.0):
    x = np.clip(img * exposure, 0, None)
    L = x @ np.array([0.2126, 0.7152, 0.0722]) + 1e-9
    Lm = L * (1 + L / white ** 2) / (1 + L)
    y = x * (Lm / L)[..., None]
    over = np.clip(y.max(-1) - 1, 0, None)[..., None]           # bright cores roll off to white
    y = np.clip(y + over * 0.6, 0, None)
    y = np.clip(y, 0, 1)
    y = np.where(y <= 0.0031308, 12.92 * y, 1.055 * np.power(y, 1 / 2.4) - 0.055)
    return Image.fromarray((np.clip(y, 0, 1) * 255 + 0.5).astype(np.uint8))


def vignette(img, strength=0.35):
    H, W = img.shape[:2]
    yy, xx = np.mgrid[0:H, 0:W]
    r = np.hypot((xx - W / 2) / (W / 2), (yy - H / 2) / (H / 2))
    img *= (1 - strength * smoothstep(0.55, 1.45, r))[..., None]


def caption(im, lines, x=60, y=None, sizes=(34, 21), alpha=(235, 195)):
    d = ImageDraw.Draw(im, 'RGBA')
    fonts = [font(sizes[0], 'semibold')] + [font(sizes[1])] * (len(lines) - 1)
    if y is None:
        y = im.height - 56 - sum(f.size * 1.45 for f in fonts)
    for i, (t, f) in enumerate(zip(lines, fonts)):
        draw_text(d, (x, y), t, f, (255, 255, 255, alpha[0] if i == 0 else alpha[1]))
        y += f.size * 1.45


def label(im, xy, text, anchor_xy, size=21, align='left'):
    """Small white label at xy with a thin leader line to anchor_xy."""
    d = ImageDraw.Draw(im, 'RGBA')
    f = font(size)
    d.line([anchor_xy, xy], fill=(255, 255, 255, 110), width=2)
    tx = xy[0] + (10 if align == 'left' else -10)
    for k, line in enumerate(text.split('\n')):
        w = text_width(d, line, f)
        px = tx if align == 'left' else tx - w
        py = xy[1] - f.size * 0.72 + k * f.size * 1.32
        draw_text(d, (px, py), line, f, (255, 255, 255, 230 if k == 0 else 185))


def save(im, name):
    os.makedirs(OUT, exist_ok=True)
    im.save(os.path.join(OUT, name), optimize=True)
    print('wrote', name)


def spot_groups(groups, seed):
    """Expand (lon, lat, size, n) groups into a leader spot plus n smaller followers."""
    rng = np.random.default_rng(seed)
    out = []
    for lon, lat, size, nsub in groups:
        out.append((lon, lat, size, 0.40))
        for _ in range(nsub):
            out.append((lon + rng.normal(0, 0.9 * size) / max(0.3, np.cos(np.radians(lat))),
                        lat + rng.normal(0, 0.5 * size), size * rng.uniform(0.25, 0.55), rng.uniform(0.42, 0.56)))
    return out


# invented spot patterns (spot modulation at the orbital period shows that the stars are spotted)
SPOTS_1 = spot_groups([(20, 52, 12, 4), (-35, 18, 8, 3), (65, -28, 9, 3), (0, 84, 15, 2), (150, 25, 10, 3),
                       (-110, -35, 9, 2)], 1)
SPOTS_2 = spot_groups([(-10, 45, 14, 3), (60, -20, 11, 3), (0, 82, 17, 1), (-90, 10, 10, 2)], 2)


# ---- 1: the triple, wide view ----------------------------------------------------------------------------
def render_triple():
    W, H = 2400, 1350
    img = np.zeros((H, W, 3))
    yy, xx = np.mgrid[0:H, 0:W]
    img += (0.005 + 0.010 * np.exp(-((xx - 1950) ** 2 + (yy - 290) ** 2) / (2 * 650 ** 2)))[..., None] * \
        np.array([0.55, 0.62, 1.0])
    dens = milky_way(img, (-200, 1500), (2600, -350), 230, 0.045, 3.0)
    starfield(img, 1700, 11, maxflux=0.8)
    starfield(img, 2200, 12, density=0.1 + 0.9 * dens, maxflux=0.35)
    # the G star, ~100-200 au away: a dazzling point (magnitude about -15.5 to -17 seen from the pair)
    point_star(img, 1940, 290, bb_rgb(T_G, 1.2), amp=60, halo=(2.0, 26), wide=(0.35, 140), spikes=(1.1, 430),
               angles=(22, 112))
    # the eclipsing pair: M4 behind and to the right, M2 in front
    R1 = 245
    R2 = R1 * R_2 / R_1 * 0.92
    x2, y2, x1, y1 = 1590, 560, 820, 790
    star_glow(img, x2, y2, R2, T_2, 0.8, sat=2.2)
    render_star(img, x2, y2, R2, T_2, inten=0.8, sat=2.2, spots=SPOTS_2, seed=4.2, tilt=18, spin=-35,
                irr=((-0.7, -0.25, 0.65), 0.15))
    star_glow(img, x1, y1, R1, T_1, 1.4)
    render_star(img, x1, y1, R1, T_1, inten=1.4, spots=SPOTS_1, seed=1.7, tilt=14, spin=25)
    img = bloom(img, levels=((4, 2.0, 0.10), (8, 4.0, 0.09), (16, 7.0, 0.08)))
    vignette(img, 0.4)
    im = tonemap(img, exposure=1.0)
    label(im, (1095, 505), 'M2 red dwarf · 0.45 M☉ · 3,550 K', (1022, 615))
    label(im, (1745, 700), 'M4 red dwarf · 0.20 M☉ · 3,160 K', (1667, 628))
    label(im, (1700, 285), 'the G star (the catalogued target)\n1.03 M☉ · 5,880 K · about 100–200 au away',
          (1912, 289), align='right')
    caption(im, ['TIC 71316387: a hidden triple star system (artist’s impression)',
                 'Two red dwarfs circle each other every 1.30 days, only 0.02 au apart. The Sun-like G star that '
                 'TESS was watching shines from far away.'])
    save(im, '01_hidden_triple.png')


# ---- 2: the deeper eclipse ---------------------------------------------------------------------------------
def render_eclipse():
    W, H = 2400, 1350
    img = np.zeros((H, W, 3)) + np.array([0.004, 0.005, 0.011])
    starfield(img, 2200, 21, maxflux=0.6)
    point_star(img, 300, 235, bb_rgb(T_G, 1.2), amp=50, halo=(1.6, 22), wide=(0.25, 120), spikes=(0.9, 360),
               angles=(30, 120))
    R1 = 540
    R2 = R1 * R_2 / R_1 * 1.04
    x1, y1 = 1500, 690
    x2, y2 = x1 - 0.40 * R1, y1 - 0.06 * R1
    star_glow(img, x1, y1, R1, T_1, 1.4)
    render_star(img, x1, y1, R1, T_1, inten=1.4, spots=SPOTS_1, seed=1.7, tilt=8, spin=-60)
    render_star(img, x2, y2, R2, T_2, inten=0.6, sat=2.2, u=(0.78, 0.06), spots=SPOTS_2, seed=4.2, tilt=8,
                spin=-25)
    img = bloom(img, levels=((4, 2.0, 0.10), (8, 4.0, 0.09), (16, 7.0, 0.08)))
    vignette(img, 0.35)
    im = tonemap(img, exposure=1.0)
    label(im, (700, 545), 'M4 red dwarf, in front\n0.21 R☉ · 3,160 K: cooler and dimmer', (x2 - R2 * 0.82, y2 - R2 * 0.3),
          align='right')
    label(im, (2085, 165), 'M2 red dwarf, behind\n0.48 R☉ · 3,550 K', (x1 + R1 * 0.52, y1 - R1 * 0.64))
    label(im, (540, 345), 'the G star, far away', (318, 252))
    caption(im, ['The deeper eclipse (artist’s impression)',
                 'Every 1.30 days the smaller M4 dwarf passes in front of the M2 dwarf, as seen from Earth. Because '
                 'the pair adds only ~2.6% of the',
                 'system’s light, TESS sees just a 0.58% dip. Half an orbit later the M4 dwarf hides behind its '
                 'partner: the 0.26% dip.'])
    save(im, '02_the_deeper_eclipse.png')


# ---- 3: map of the system ----------------------------------------------------------------------------------
def ring(img, cx, cy, r, rgb, amp, w, dash=None, aspect=1.0):
    sl, xx, yy = window(img, cx, cy, r + 6 * w + 4)
    dx, dy = xx - cx, (yy - cy) / aspect
    d = np.hypot(dx, dy)
    f = amp * np.exp(-0.5 * ((d - r) / w) ** 2) + 0.25 * amp * np.exp(-0.5 * ((d - r) / (4 * w)) ** 2)
    if dash is not None:
        ang = np.arctan2(dy, dx) * r
        f *= (np.mod(ang, dash[0]) < dash[1])
    img[sl] += f[..., None] * rgb


def segment(img, p0, p1, rgb, amp, w):
    (x0, y0), (x1, y1) = p0, p1
    cx, cy, rr = (x0 + x1) / 2, (y0 + y1) / 2, np.hypot(x1 - x0, y1 - y0) / 2 + 6 * w
    sl, xx, yy = window(img, cx, cy, rr)
    vx, vy = x1 - x0, y1 - y0
    t = np.clip(((xx - x0) * vx + (yy - y0) * vy) / (vx * vx + vy * vy), 0, 1)
    d = np.hypot(xx - (x0 + t * vx), yy - (y0 + t * vy))
    img[sl] += (amp * np.exp(-0.5 * (d / w) ** 2))[..., None] * rgb


def render_map():
    W, H = 2400, 1350
    img = np.zeros((H, W, 3)) + np.array([0.004, 0.005, 0.010])
    starfield(img, 1400, 31, maxflux=0.25)
    grey = np.array([0.55, 0.62, 0.8])
    cG, c1, c2 = bb_rgb(T_G, 1.05), bb_rgb(T_1, 1.2), bb_rgb(T_2, 1.2)
    # outer orbit, drawn face-on and circular (true shape and tilt unknown); a = 150 au -> 1000 px
    bx, by, A = 760, 690, 1000
    fG = (M_1 + M_2) / (M_G + M_1 + M_2)
    rG, rB = fG * A, (1 - fG) * A
    ring(img, bx, by, rG, cG, 0.10, 1.2)
    ring(img, bx, by, rB, c1, 0.10, 1.2)
    th = np.radians(200)
    gx, gy = bx + rG * np.cos(th), by - rG * np.sin(th)
    ex, ey = bx - rB * np.cos(th), by + rB * np.sin(th)
    ring(img, gx, gy, 30 / 150 * A, grey, 0.07, 1.0, dash=(14, 7))        # Neptune's orbit for scale
    segment(img, (bx - 9, by), (bx + 9, by), grey, 0.25, 1.0)
    segment(img, (bx, by - 9), (bx, by + 9), grey, 0.25, 1.0)
    point_star(img, gx, gy, cG, amp=35, halo=(1.2, 14), wide=(0.15, 60), spikes=(0.5, 160), angles=(45, 135))
    point_star(img, ex, ey, c1, amp=6, halo=(0.5, 9), wide=(0.06, 40), spikes=(0.15, 80), angles=(45, 135))
    # inset: the inner binary to scale, 4.35 Rsun -> 480 px
    ix, iy, IR = 1930, 610, 400
    ring(img, ix, iy, IR, grey, 0.10, 1.4)
    s = 480 / A_IN
    f1 = M_2 / (M_1 + M_2)
    r1o, r2o = f1 * A_IN * s, (1 - f1) * A_IN * s
    ring(img, ix, iy, r1o, c1, 0.12, 1.1)
    ring(img, ix, iy, r2o, c2, 0.12, 1.1)
    ph = np.radians(35)
    p1 = (ix - r1o * np.cos(ph), iy + r1o * np.sin(ph))
    p2 = (ix + r2o * np.cos(ph), iy - r2o * np.sin(ph))
    star_glow(img, p1[0], p1[1], R_1 * s, T_1, 1.4)
    star_glow(img, p2[0], p2[1], R_2 * s, T_2, 0.8, sat=2.2)
    render_star(img, p1[0], p1[1], R_1 * s, T_1, inten=1.4, spots=SPOTS_1, seed=1.7, tilt=40, spin=20)
    render_star(img, p2[0], p2[1], R_2 * s, T_2, inten=0.8, sat=2.2, spots=SPOTS_2, seed=4.2, tilt=40, spin=-20)
    segment(img, (ix - 7, iy), (ix + 7, iy), grey, 0.25, 1.0)
    segment(img, (ix, iy - 7), (ix, iy + 7), grey, 0.25, 1.0)
    # leader from the pair on the outer orbit to the inset
    ang = np.arctan2(iy - ey, ix - ex)
    segment(img, (ex + 16 * np.cos(ang), ey + 16 * np.sin(ang)),
            (ix - IR * np.cos(ang), iy - IR * np.sin(ang)), grey, 0.10, 1.0)
    img = bloom(img, levels=((4, 2.0, 0.06), (8, 4.0, 0.05)))
    im = tonemap(img, exposure=1.0)
    d = ImageDraw.Draw(im, 'RGBA')
    f_t, f_s, f_x = font(40, 'semibold'), font(22), font(19)
    d.text((60, 46), 'The TIC 71316387 system from above', font=f_t, fill=(255, 255, 255, 235))
    d.text((60, 102), 'Most likely layout (outer orbit drawn circular and face-on; its true size, shape and tilt '
                      'are not yet measured)', font=f_s, fill=(255, 255, 255, 170))
    label(im, (gx + 165, gy + 228), 'G star · 1.03 M☉ · 5,880 K', (gx + 10, gy + 10))
    label(im, (gx + 40, gy - 268), 'Neptune’s orbit (30 au), for scale', (gx + 120, gy - 170))
    label(im, (ex - 40, ey + 140), 'the eclipsing pair\n0.65 M☉ in total', (ex - 6, ey + 10), align='right')
    label(im, (bx + 60, by + 60), 'centre of mass', (bx + 6, by + 6))
    d = ImageDraw.Draw(im, 'RGBA')
    # outer-orbit scale bar: 50 au
    sb = 50 / 150 * A
    y0 = H - 120
    d.line([(60, y0), (60 + sb, y0)], fill=(255, 255, 255, 200), width=3)
    d.text((60, y0 + 12), '50 au', font=f_x, fill=(255, 255, 255, 200))
    d.text((60, y0 + 42), 'Outer orbit: about 100–200 au, one lap every ~600–1,700 years (estimates; '
                          'range 30–650 au)', font=f_x, fill=(255, 255, 255, 170))
    # inset texts
    d.text((ix - IR + 20, iy + IR + 26), 'Zoom on the eclipsing pair, to scale', font=f_s, fill=(255, 255, 255, 220))
    draw_text(d, (ix - IR + 20, iy + IR + 60), 'Inner orbit 4.35 R☉ = 0.020 au across, one lap every 1.30 days',
              f_x, (255, 255, 255, 170))
    sb2 = 1.0 * s
    d.line([(ix + IR - 20 - sb2, iy + IR - 20), (ix + IR - 20, iy + IR - 20)], fill=(255, 255, 255, 200), width=3)
    draw_text(d, (ix + IR - 20 - sb2, iy + IR - 52), '1 R☉', f_x, (255, 255, 255, 200))
    label(im, (p1[0] - 40, p1[1] + 110), 'M2 dwarf · 0.45 M☉', (p1[0] - 10, p1[1] + R_1 * s), align='right')
    label(im, (p2[0] - 40, p2[1] - 95), 'M4 dwarf · 0.20 M☉', (p2[0] - 6, p2[1] - R_2 * s), align='right')
    save(im, '03_system_map.png')


# ---- 4: what telescopes see --------------------------------------------------------------------------------
def render_telescopes():
    import csv
    REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
    rows = list(csv.DictReader(open(os.path.join(REPO, 'tess-search', 'exofop_leads', 'TIC71316387_campaign',
                                                 'partial', 'localise', 'gaia', 'neighbours_120arcsec.csv'))))
    ra0, de0 = float(rows[0]['ra']), float(rows[0]['dec'])
    stars = [((float(r['ra']) - ra0) * np.cos(np.radians(de0)) * 3600, (float(r['dec']) - de0) * 3600,
              float(r['phot_g_mean_mag'])) for r in rows]
    S = 720                       # panel size, px
    pad, top = 30, 150
    W, H = 3 * S + 4 * pad, S + top + 160
    canvas = Image.new('RGB', (W, H), (8, 9, 14))

    def to_img(a, gamma=0.5):
        a = np.clip(a / a.max(), 0, 1) ** gamma
        tint = np.array([0.86, 0.92, 1.0])
        return Image.fromarray((a[..., None] * tint * 255).astype(np.uint8))

    # panel A: TESS, 9 x 9 pixels of 21"
    npx, scale = 9, 21.0
    sub = 8
    g = np.zeros((npx * sub, npx * sub))
    yy, xx = (np.mgrid[0:npx * sub, 0:npx * sub] + 0.5) / sub - npx / 2
    for dE, dN, gm in stars + [(0.4 * np.sin(np.radians(-20)), 0.4 * np.cos(np.radians(-20)), 16.4)]:
        f = 10 ** (-0.4 * (gm - 11.8))
        x, y = -dE / scale + 0.18, -dN / scale + 0.12   # East left; pixel phase arbitrary
        g += f * np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / 0.62 ** 2)
    pix = g.reshape(npx, sub, npx, sub).sum(axis=(1, 3))
    pix += 2e-4 * pix.max() * (1 + 0.15 * np.random.default_rng(5).standard_normal(pix.shape))
    imA = to_img(np.clip(pix, 0, None), 0.35).resize((S, S), Image.NEAREST)
    dA = ImageDraw.Draw(imA)
    for k in range(1, npx):
        dA.line([(k * S / npx, 0), (k * S / npx, S)], fill=(60, 66, 80), width=1)
        dA.line([(0, k * S / npx), (S, k * S / npx)], fill=(60, 66, 80), width=1)
    # panel B: ground-based, 1" seeing, 8" field
    n = 480
    fov = 8.0
    yy, xx = (np.mgrid[0:n, 0:n] + 0.5) / n * fov - fov / 2
    sig = 1.0 / 2.355
    comp = (-0.4 * np.sin(np.radians(-20)), 0.4 * np.cos(np.radians(-20)))   # illustrative: 0.4" at PA 340 deg
    fB = np.exp(-0.5 * (xx ** 2 + yy ** 2) / sig ** 2)
    fB += 10 ** (-0.4 * 4.6) * np.exp(-0.5 * ((xx - comp[0]) ** 2 + (yy + comp[1]) ** 2) / sig ** 2)
    imB = to_img(fB, 0.45).resize((S, S), Image.LANCZOS)
    # panel C: adaptive optics in Ks on an 8-m telescope, 2.4" field, Strehl ~0.5
    n = 720
    fov = 3.0
    yy, xx = (np.mgrid[0:n, 0:n] + 0.5) / n * fov - fov / 2

    def psf(x0, y0):
        r = np.hypot(xx - x0, yy - y0) + 1e-9
        lam_d = 2.15e-6 / 8.0 * 206265          # lambda / D in arcsec
        v = np.pi * r / lam_d
        airy = (2 * j1(v) / v) ** 2
        halo = 0.02 * np.exp(-r / 0.25)
        return 0.5 * airy + halo
    fC = psf(0, 0) + 10 ** (-0.4 * DKS) * psf(comp[0], -comp[1])
    imC = to_img(fC, 0.28)
    dC = ImageDraw.Draw(imC, 'RGBA')
    rr = 1.24 / fov * S
    dC.ellipse([S / 2 - rr, S / 2 - rr, S / 2 + rr, S / 2 + rr], outline=(120, 170, 255, 150), width=2)

    panels = [(imA, 'What TESS sees', '21″ pixels: the G star and the eclipsing pair fall in one pixel',
               ('63″ = 3 pixels', 3 * S / npx)),
              (imB, 'A normal telescope image', '1″ seeing: still a single star', ('1″', S / 8.0)),
              (imC, 'Adaptive optics, near-infrared (8-m)', 'a faint companion about 0.4″ away (ΔK ≈ 2.6)',
               ('0.5″', 0.5 / 3.0 * S))]
    d = ImageDraw.Draw(canvas, 'RGBA')
    f_t, f_s, f_x = font(40, 'semibold'), font(23), font(19)
    d.text((pad, 34), 'How TIC 71316387 looks through different telescopes (simulated)', font=f_t,
           fill=(255, 255, 255, 235))
    d.text((pad, 92), 'Gaia star positions and the adopted brightness of the M-dwarf pair; its exact position '
                      'angle is not known, so the one shown is illustrative.', font=f_s, fill=(255, 255, 255, 170))
    for k, (im, t1, t2, (sbl, sbw)) in enumerate(panels):
        x0 = pad + k * (S + pad)
        canvas.paste(im, (x0, top))
        d.text((x0, top + S + 18), t1, font=font(27, 'semibold'), fill=(255, 255, 255, 230))
        d.text((x0, top + S + 58), t2, font=f_x, fill=(255, 255, 255, 175))
        d.line([(x0 + S - 24 - sbw, top + S - 30), (x0 + S - 24, top + S - 30)], fill=(255, 255, 255, 220), width=3)
        bb = d.textbbox((0, 0), sbl, font=f_x)
        d.text((x0 + S - 24 - sbw / 2 - (bb[2] - bb[0]) / 2, top + S - 62), sbl, font=f_x, fill=(255, 255, 255, 220))
    x0 = pad + 2 * (S + pad)
    d.text((x0 + 18, top + 16), 'blue circle: 1.24″, the TESS 95% limit', font=f_x, fill=(150, 190, 255, 220))
    cx, cy = x0 + S / 2 + comp[0] / 3.0 * S, top + S / 2 - comp[1] / 3.0 * S
    d.line([(cx + 14, cy - 14), (cx + 90, cy - 90)], fill=(255, 255, 255, 150), width=2)
    d.text((cx + 96, cy - 112), 'the M-dwarf pair\n(its two stars are only\n0.00007″ apart)', font=f_x,
           fill=(255, 255, 255, 210))
    d.text((pad, H - 46), 'North up, east left. Panel widths: 189″, 8″ and 3″.', font=f_x,
           fill=(255, 255, 255, 140))
    save(canvas, '04_through_the_telescope.png')


if __name__ == '__main__':
    import sys
    which = sys.argv[1:] or ['1', '2', '3', '4']
    for k in which:
        {'1': render_triple, '2': render_eclipse, '3': render_map, '4': render_telescopes}[k]()
