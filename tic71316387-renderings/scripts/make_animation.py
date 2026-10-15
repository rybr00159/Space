"""Animation: how the two red dwarfs of TIC 71316387 make the dips in its light curve.

Three synchronised views: the orbit seen from above, the pair seen from Earth, and the brightness curve drawn as
the stars move. Writes ../animation/tic71316387_eclipses.mp4 (1920x1080, 30 fps), a README preview
../animation/tic71316387_eclipses.gif (960x540, 15 fps) and a still frame ../animation/tic71316387_eclipses_poster.png.

Not to scale: the stars are drawn 3.4 M2-radii apart instead of about 9, and the orbit is tilted 85 degrees instead
of about 89, so the eclipses last longer and the paths are easy to follow. The radius ratio (0.44), the mass
ratio and the 5,770 : 2,600 ppm depth ratio follow the campaign's adopted model. The brightness curve is computed from
the same geometry that is drawn.

Requires numpy, scipy, pillow, matplotlib (fonts) and ffmpeg.  Usage:  python make_animation.py
"""
import os
import shutil
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
import numpy as np                                                      # noqa: E402
from PIL import Image, ImageDraw                                         # noqa: E402

import make_renderings as mr                                            # noqa: E402

OUT = os.path.join(mr.HERE, '..', 'animation')
W, H, FPS = 1920, 1080, 30
GIF_FPS, GIF_W, GIF_COLOURS = 15, 960, 128      # README preview
P_DAYS = 1.30367228
TESS_DEEP, TESS_SHALLOW = 5770, 2600          # adopted TESS depths (ppm)

# ---- illustrative geometry (units of the M2 dwarf's radius) -------------------------------------------------
K = 0.44                                      # radius ratio M4/M2
Q = mr.M_2 / mr.M_1                           # mass ratio
A = 3.4                                       # separation (really ~9.1)
INC = np.radians(85.0)                        # tilt (really ~89 deg)
R1O, R2O = A * Q / (1 + Q), A / (1 + Q)       # orbit radii of M2 and M4 about the centre of mass
LD = (0.38, 0.28)                             # quadratic limb darkening (approximate, red-optical)

# ---- colours ---------------------------------------------------------------------------------------------------
BG, PANEL, BORDER = (10, 12, 19), (16, 19, 28), (40, 46, 62)
INK, INK2, MUTED = (240, 239, 236), (195, 194, 183), (137, 135, 129)
GRID, AXIS = (44, 48, 60), (70, 74, 86)
BLUE = (57, 135, 229)                         # categorical slot 1, dark-mode step


def lin(c):
    c = np.asarray(c, float) / 255
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


# ---- light-curve model -----------------------------------------------------------------------------------------
def _disk(n, R):
    g = (np.arange(n) + 0.5) / n * 2 - 1
    x, y = np.meshgrid(g * R, g * R)
    r = np.hypot(x, y) / R
    mu = np.sqrt(np.clip(1 - r ** 2, 0, 1))
    return x, y, np.where(r < 1, 1 - LD[0] * (1 - mu) - LD[1] * (1 - mu) ** 2, 0.0) * (2 * R / n) ** 2


X1, Y1, W1 = _disk(900, 1.0)
X2, Y2, W2 = _disk(500, K)


def positions(phase):
    """Orbit-plane positions of M2 and M4 (x right, y away from Earth). Phase 0: M4 at the left; 0.25: M4 in front
    of M2 (deep eclipse); 0.75: M4 behind M2 (shallow eclipse)."""
    th = 2 * np.pi * phase + np.pi
    return (-R1O * np.cos(th), -R1O * np.sin(th)), (R2O * np.cos(th), R2O * np.sin(th))


def flux(phase, J):
    """Brightness of the pair (1 = both stars fully visible) for surface-brightness ratio J (M4/M2)."""
    (x2, y2), (x4, y4) = positions(phase)
    dx, dy = x4 - x2, (y4 - y2) * np.cos(INC)
    f1, f2 = W1.sum(), J * W2.sum()
    if np.hypot(dx, dy) >= 1 + K:
        return 1.0
    if y4 < y2:                                                    # M4 nearer to Earth: covers part of M2
        blocked = W1[np.hypot(X1 - dx, Y1 - dy) < K].sum()
    else:                                                          # M4 behind M2
        blocked = J * W2[np.hypot(X2 + dx, Y2 + dy) < 1].sum()
    return (f1 + f2 - blocked) / (f1 + f2)


def tune_j():
    """Surface-brightness ratio that reproduces the observed deep : shallow depth ratio."""
    lo, hi = 0.2, 1.2
    for _ in range(40):
        j = 0.5 * (lo + hi)
        r = (1 - flux(0.25, j)) / (1 - flux(0.75, j))
        lo, hi = (j, hi) if r > TESS_DEEP / TESS_SHALLOW else (lo, j)
    return j


# ---- timeline: two orbits, slowed down during the eclipses ---------------------------------------------------
SLOW, T_ORBIT, HOLD_START, HOLD_END = 1 / 3, 6.0, 0.7, 4.5


def eclipse_weight(ph):
    d = np.minimum(np.abs((ph - 0.25 + 0.5) % 1 - 0.5), np.abs((ph - 0.75 + 0.5) % 1 - 0.5))
    return mr.smoothstep(0.115, 0.08, d)


def timeline():
    ph = np.linspace(0, 2, 40001)
    speed = 1 - (1 - SLOW) * eclipse_weight(ph)
    t = np.concatenate([[0], np.cumsum(np.diff(ph) / (0.5 * (speed[1:] + speed[:-1])))]) * T_ORBIT
    n_run = int(round(t[-1] * FPS))
    run = np.interp(np.arange(n_run + 1) / FPS * t[-1] / (n_run / FPS), t, ph)
    return np.concatenate([np.zeros(int(HOLD_START * FPS)), run, np.full(int(HOLD_END * FPS), 2.0)])


# ---- layout ----------------------------------------------------------------------------------------------------
PA = (40, 150, 640, 650)                      # seen from above
PB = (680, 150, 1880, 650)                    # seen from Earth
PC = (40, 690, 1880, 1060)                    # brightness
A_C, A_S = (300, 250), 76.0                              # panel-A centre and px per unit
A_R1, A_R2 = A_S, K * A_S                                  # star radii, same scale as the orbit
B_C, B_S = (600, 225), 112.0                               # panel-B centre and px per unit
X0, X1P, YT, YB = 170, 1820, 742, 975                     # plot box of panel C
YMIN, YMAX = 70.0, 110.0


def cx(t):
    return X0 + (X1P - X0) * t / (2 * P_DAYS)


def cy(pct):
    return YB - (YB - YT) * (pct - YMIN) / (YMAX - YMIN)


def text_img(lines, sizes, colours, weights=None, align='left', spacing=1.35, shadow=False):
    """Render lines of text into a tight RGBA image."""
    weights = weights or ['regular'] * len(lines)
    fonts = [mr.font(s, w) for s, w in zip(sizes, weights)]
    tmp = ImageDraw.Draw(Image.new('RGBA', (8, 8)))
    widths = [mr.text_width(tmp, t, f) for t, f in zip(lines, fonts)]
    w = int(max(widths)) + 8
    h = int(sum(f.size * spacing for f in fonts)) + 6
    im = Image.new('RGBA', (w, h), (0, 0, 0, 0))
    d = ImageDraw.Draw(im)
    y = 2
    for t, f, c, tw in zip(lines, fonts, colours, widths):
        x = 4 if align == 'left' else (w - tw) / 2 if align == 'center' else w - 4 - tw
        mr.draw_text(d, (x, y), t, f, c + (255,) if len(c) == 3 else c, shadow=shadow)
        y += f.size * spacing
    return im


def paste(frame, im, xy, alpha=1.0, anchor='lt'):
    if alpha <= 0.004:
        return
    x, y = xy
    if anchor[0] == 'm':
        x -= im.width / 2
    elif anchor[0] == 'r':
        x -= im.width
    if anchor[1] == 'm':
        y -= im.height / 2
    elif anchor[1] == 'b':
        y -= im.height
    if alpha < 1:
        im = im.copy()
        im.putalpha(im.getchannel('A').point(lambda v: int(v * alpha)))
    frame.alpha_composite(im, (int(round(x)), int(round(y))))


def supersampled(size, draw_fn, scale=3):
    """RGBA layer drawn at `scale` times the size and shrunk, for anti-aliased lines."""
    big = Image.new('RGBA', (size[0] * scale, size[1] * scale), (0, 0, 0, 0))
    draw_fn(ImageDraw.Draw(big), scale)
    return big.resize(size, Image.LANCZOS)


def mask_from(size, draw_fn, scale=3):
    big = Image.new('L', (size[0] * scale, size[1] * scale), 0)
    draw_fn(ImageDraw.Draw(big), scale)
    return np.asarray(big.resize(size, Image.LANCZOS), float) / 255


def dashed_ellipse(d, box, width, fill, n=72, duty=0.55, start=0, end=360):
    a0 = start
    step = (end - start) / n
    for k in range(n):
        d.arc(box, a0 + k * step, a0 + (k + duty) * step, fill=fill, width=width)


# ---- static layers ---------------------------------------------------------------------------------------------
def static_base():
    base = Image.new('RGBA', (W, H), BG + (255,))

    def panels(d, s):
        for p in (PA, PB, PC):
            d.rounded_rectangle([p[0] * s, p[1] * s, p[2] * s, p[3] * s], radius=14 * s, fill=PANEL + (255,))
        # panel C: grid, axes, ticks
        for pct in (80, 90, 100):
            d.line([(X0 * s, cy(pct) * s), (X1P * s, cy(pct) * s)], fill=GRID + (255,), width=s)
        d.line([(X0 * s, YB * s), (X1P * s, YB * s)], fill=AXIS + (255,), width=s)
        for tt in np.arange(0, 2.61, 0.5):
            d.line([(cx(tt) * s, YB * s), (cx(tt) * s, (YB + 7) * s)], fill=AXIS + (255,), width=s)
    base.alpha_composite(supersampled((W, H), panels))
    d = ImageDraw.Draw(base)
    mr.draw_text(d, (60, 26), 'How two stars make the dips in the light of TIC 71316387', mr.font(44, 'semibold'),
                 INK + (255,), shadow=False)
    mr.draw_text(d, (60, 88), 'Two red dwarfs orbit each other every 1.30 days. Not to scale, and slowed down '
                 'during the eclipses.', mr.font(26), INK2 + (255,), shadow=False)
    mr.draw_text(d, (64, 700), 'Brightness of the two red dwarfs', mr.font(28, 'semibold'), INK + (255,),
                 shadow=False)
    f = mr.font(22)
    for pct in (80, 90, 100):
        lab = f'{pct}%'
        mr.draw_text(d, (X0 - 14 - mr.text_width(d, lab, f), cy(pct) - 14), lab, f, INK2 + (255,), shadow=False)
    for tt in np.arange(0, 2.61, 0.5):
        lab = f'{tt:g}' + (' days' if tt == 2.5 else '')
        x = cx(tt) - (mr.text_width(d, f'{tt:g}', f) / 2)
        mr.draw_text(d, (x, YB + 8), lab, f, INK2 + (255,), shadow=False)
    mr.draw_text(d, (64, 1020), 'TESS measures only 0.58% and 0.26%, because the bright G star next to the pair '
                 'adds about 97% of the light. Dips drawn wider than in reality.', mr.font(21), MUTED + (255,),
                 shadow=False)
    return base


def legend_items():
    """(text, local text position, icon centre, temperature, icon radius, brightness, saturation) for panel B."""
    d = ImageDraw.Draw(Image.new('RGBA', (8, 8)))
    rows = [('M2 dwarf: bigger and hotter (3,550 K)', mr.T_1, 12, 1.4, 1.7),
            ('M4 dwarf: smaller and cooler (3,160 K)', mr.T_2, 6, 0.8, 2.2)]
    x = PB[2] - PB[0] - 24 - max(mr.text_width(d, r[0], mr.font(22)) for r in rows)
    return [(txt, (x, 22 + 40 * k), (x - 24, 22 + 40 * k + 14), T, r, inten, sat)
            for k, (txt, T, r, inten, sat) in enumerate(rows)]


def static_overlay(j):
    """Text and frames drawn on top of the two star panels."""
    ov = Image.new('RGBA', (W, H), (0, 0, 0, 0))

    def frames(d, s):
        for p in (PA, PB):
            d.rounded_rectangle([p[0] * s, p[1] * s, p[2] * s, p[3] * s], radius=14 * s, outline=BORDER + (255,),
                                width=s)
        # arrow to Earth in panel A
        ax, ay = PA[0] + A_C[0], PA[1] + 452
        d.line([(ax * s, ay * s), (ax * s, (ay + 30) * s)], fill=INK2 + (255,), width=2 * s)
        d.polygon([((ax - 8) * s, (ay + 26) * s), ((ax + 8) * s, (ay + 26) * s), (ax * s, (ay + 40) * s)],
                  fill=INK2 + (255,))
    ov.alpha_composite(supersampled((W, H), frames))
    d = ImageDraw.Draw(ov)
    mr.draw_text(d, (PA[0] + 22, PA[1] + 14), 'Seen from above', mr.font(28, 'semibold'), INK + (255,))
    mr.draw_text(d, (PB[0] + 22, PB[1] + 14), 'Seen from Earth', mr.font(28, 'semibold'), INK + (255,))
    mr.draw_text(d, (PB[0] + 22, PB[1] + 54), 'as TESS would see it if it could zoom in', mr.font(22),
                 INK2 + (255,))
    f = mr.font(21)
    mr.draw_text(d, (PA[0] + A_C[0] + 16, PA[1] + 452), 'towards Earth and TESS', f, INK2 + (255,))
    mr.draw_text(d, (PA[0] + 468, PA[1] + 56), 'far side', f, MUTED + (255,))
    mr.draw_text(d, (PA[0] + 468, PA[1] + 424), 'near side', f, MUTED + (255,))
    # legend in panel B (top right); the star icons are drawn into the panel itself
    for txt, (x, y), *_ in legend_items():
        mr.draw_text(d, (PB[0] + x, PB[1] + y), txt, mr.font(22), INK2 + (255,))
    return ov


def panel_a_background():
    w, h = PA[2] - PA[0], PA[3] - PA[1]
    img = np.zeros((h, w, 3)) + lin(PANEL)
    mr.starfield(img, 160, 41, maxflux=0.12)
    c1, c2 = mr.bb_rgb(mr.T_1, 1.7), mr.bb_rgb(mr.T_2, 2.2)
    for r, col, amp in ((R2O, c2, 0.16), (R1O, c1, 0.16)):
        m = mask_from((w, h), lambda d, s: d.ellipse([(A_C[0] - r * A_S) * s, (A_C[1] - r * A_S) * s,
                                                     (A_C[0] + r * A_S) * s, (A_C[1] + r * A_S) * s],
                                                    outline=255, width=2 * s))
        img += amp * m[..., None] * col
    m = mask_from((w, h), lambda d, s: (d.line([((A_C[0] - 7) * s, A_C[1] * s), ((A_C[0] + 7) * s, A_C[1] * s)],
                                               fill=255, width=2 * s),
                                        d.line([(A_C[0] * s, (A_C[1] - 7) * s), (A_C[0] * s, (A_C[1] + 7) * s)],
                                               fill=255, width=2 * s)))
    img += 0.12 * m[..., None]
    return img


def panel_b_background():
    w, h = PB[2] - PB[0], PB[3] - PB[1]
    img = np.zeros((h, w, 3)) + lin(PANEL)
    mr.starfield(img, 420, 42, maxflux=0.14)
    for _, _, (ix, iy), T, r, inten, sat in legend_items():
        mr.star_glow(img, ix, iy, r, T, inten, sat=sat)
        mr.render_star(img, ix, iy, r, T, inten=inten, sat=sat, gran=0.0, mottle=0.0)
    rx, ry = R2O * B_S, R2O * np.cos(INC) * B_S

    def half(far):
        def fn(d, s):
            box = [(B_C[0] - rx) * s, (B_C[1] - ry) * s, (B_C[0] + rx) * s, (B_C[1] + ry) * s]
            dashed_ellipse(d, box, 2 * s, 255, n=40, start=180 if far else 0, end=360 if far else 180)
        return mask_from((w, h), fn)
    return img, half(True), half(False)


def curve_layer(times, pct):
    def fn(d, s):
        pts = [(cx(t) * s, cy(p) * s) for t, p in zip(times, pct)]
        d.line(pts, fill=BLUE + (255,), width=3 * s, joint='curve')
    return supersampled((W, H), fn)


def glow(img, x, y, R, T, inten, sat):
    """Rim glow plus a modest halo (no wide halo or bloom: fewer changing pixels between frames)."""
    c = mr.bb_rgb(T * 0.9, sat)
    mr.add_glow(img, x, y, R, c, 0.55 * inten, max(2.0, 0.018 * R))
    mr.add_glow(img, x, y, R, c, 0.10 * inten, 0.10 * R, reach=6)


# ---- per-frame panels ------------------------------------------------------------------------------------------
def render_panel_a(bg, phase, w_ecl):
    img = bg.copy()
    (x2, y2), (x4, y4) = positions(phase)
    # the column of sky that M2 covers as seen from Earth: M4 inside it means an eclipse
    dx = np.abs(np.arange(img.shape[1]) + 0.5 - (A_C[0] + x2 * A_S))
    fill = mr.smoothstep(A_R1 + 1.0, A_R1 - 1.0, dx)
    edge = np.exp(-0.5 * ((dx - A_R1) / 1.1) ** 2)
    img += ((0.010 + 0.035 * w_ecl) * fill + (0.035 + 0.16 * w_ecl) * edge)[None, :, None] * np.array([0.35, 0.6, 1.0])
    for (x, y), r, T, inten, sat in (((x4, y4), A_R2, mr.T_2, 0.8, 2.2), ((x2, y2), A_R1, mr.T_1, 1.4, 1.7)):
        px, py = A_C[0] + x * A_S, A_C[1] - y * A_S
        glow(img, px, py, r, T, inten, sat)
        mr.render_star(img, px, py, r, T, inten=inten, sat=sat, spots=mr.SPOTS_1 if r == A_R1 else mr.SPOTS_2,
                       seed=1.7 if r == A_R1 else 4.2, tilt=80, spin=0)
    return mr.tonemap(img)


def render_panel_b(bg, far_path, near_path, phase):
    img = bg.copy()
    img += 0.07 * far_path[..., None]
    (x2, y2), (x4, y4) = positions(phase)
    stars = [dict(x=x2, y=y2, R=B_S, T=mr.T_1, inten=1.4, sat=1.7, spots=mr.SPOTS_1, seed=1.7),
             dict(x=x4, y=y4, R=K * B_S, T=mr.T_2, inten=0.8, sat=2.2, spots=mr.SPOTS_2, seed=4.2)]
    for s in sorted(stars, key=lambda s: -s['y']):                 # farthest first
        s['px'], s['py'] = B_C[0] + s['x'] * B_S, B_C[1] - s['y'] * np.cos(INC) * B_S
        glow(img, s['px'], s['py'], s['R'], s['T'], s['inten'], s['sat'])
        mr.render_star(img, s['px'], s['py'], s['R'], s['T'], inten=s['inten'], sat=s['sat'], spots=s['spots'],
                       seed=s['seed'], tilt=90 - np.degrees(INC), spin=-360 * phase + 40)
    # the near half of M4's path passes in front of M2, but not across M4 itself
    m4 = stars[1]
    yy, xx = np.mgrid[0:img.shape[0], 0:img.shape[1]]
    keep = np.hypot(xx - m4['px'], yy - m4['py']) > m4['R'] + 3
    img += 0.07 * (near_path * keep)[..., None]
    im = mr.tonemap(img)
    return im, stars


def main():
    j = tune_j()
    d_deep, d_shallow = 1 - flux(0.25, j), 1 - flux(0.75, j)
    print(f'surface-brightness ratio J = {j:.3f}; pair-only dips {100 * d_deep:.1f}% and {100 * d_shallow:.1f}%; '
          f'pair light fraction for the TESS depths {TESS_DEEP / 1e6 / d_deep * 100:.2f}%')
    phases = timeline()
    n = phases.size
    print(f'{n} frames ({n / FPS:.1f} s)')
    # brightness curve on a fine grid
    fine = np.linspace(0, 2, 6001)
    pct = np.array([100 * flux(p, j) for p in fine])
    times = fine * P_DAYS

    base = static_base()
    over = static_overlay(j)
    curve = curve_layer(times, pct)
    bg_a = panel_a_background()
    bg_b, far_path, near_path = panel_b_background()

    # text sprites
    cap_sizes, cap_cols = (32, 26), (INK, INK2)
    captions = {
        'out': text_img(['Both stars are in full view, so we see all of their light.'], (32,), (INK,),
                        align='center', shadow=True),
        'deep': text_img(['The small M4 dwarf passes in front of the M2 dwarf',
                          'and blocks part of its light: the deep dip.'], cap_sizes, (INK, INK), align='center',
                         shadow=True),
        'shallow': text_img(['The M4 dwarf passes behind the M2 dwarf and is hidden: the shallow dip.',
                             'Hiding the small, cooler star costs less light than covering part of the hotter one.'],
                            cap_sizes, cap_cols, align='center', shadow=True),
        'end': text_img(['Each 1.30-day orbit gives one deep and one shallow dip.',
                         'A planet would make every dip the same: these are two eclipsing stars.'],
                        (32, 32), (INK, INK), align='center', shadow=True),
    }
    slow = text_img(['slowed down during the eclipse'], (21,), (MUTED,))
    lined = text_img(['lined up: eclipse'], (21,), (INK,))
    m2_lab, m4_lab = text_img(['M2'], (24,), (INK2,), weights=['semibold']), \
        text_img(['M4'], (24,), (INK2,), weights=['semibold'])
    dips = [(0.25, text_img([f'deep dip: about {100 * d_deep:.0f}%', 'TESS measures 0.58%'], (23, 20), (INK, INK2))),
            (0.75, text_img([f'shallow dip: about {100 * d_shallow:.0f}%', 'TESS measures 0.26%'], (23, 20),
                            (INK, INK2))),
            (1.25, text_img(['deep'], (22,), (INK2,))),
            (1.75, text_img(['shallow'], (22,), (INK2,)))]
    brackets = [(0.25, 0.75, 101.8, '0.65 days: the "planet" period'),
                (0.25, 1.25, 107.4, '1.30 days: one orbit')]
    br_layers = []
    for p0, p1, yv, txt in brackets:
        lay = supersampled((W, H), lambda d, s, p0=p0, p1=p1, yv=yv: d.line(
            [(cx(p0 * P_DAYS) * s, (cy(yv) + 6) * s), (cx(p0 * P_DAYS) * s, cy(yv) * s),
             (cx(p1 * P_DAYS) * s, cy(yv) * s), (cx(p1 * P_DAYS) * s, (cy(yv) + 6) * s)],
            fill=INK2 + (255,), width=2 * s))
        br_layers.append((p1, lay, text_img([txt], (21,), (INK2,), align='center'),
                          (cx(0.5 * (p0 + p1) * P_DAYS), cy(yv) - 2)))

    # caption weights per frame, smoothed into cross-fades
    w_ecl = eclipse_weight(phases)
    state = np.zeros((n, 4))
    for i, ph in enumerate(phases):
        frac = ph % 1
        if i >= n - int(HOLD_END * FPS) + 15:
            state[i, 3] = 1
        elif abs(frac - 0.25) < 0.085 and ph < 2:
            state[i, 1] = 1
        elif abs(frac - 0.75) < 0.085:
            state[i, 2] = 1
        else:
            state[i, 0] = 1
    ker = np.ones(9) / 9
    state = np.stack([np.convolve(np.pad(state[:, k], 4, mode='edge'), ker, 'valid') for k in range(4)], 1)

    def make_frame(i):
        ph = phases[i]
        frame = base.copy()
        frame.paste(render_panel_a(bg_a, ph, float(w_ecl[i])), (PA[0], PA[1]))
        im_b, stars = render_panel_b(bg_b, far_path, near_path, ph)
        frame.paste(im_b, (PB[0], PB[1]))
        # ghost outline of M4 while it is hidden behind M2
        m2, m4 = stars
        dist = np.hypot(m4['px'] - m2['px'], m4['py'] - m2['py'])
        behind = m4['y'] > m2['y']
        if behind and dist < m2['R'] + m4['R']:
            ghost = Image.new('RGBA', (W, H), (0, 0, 0, 0))
            box = [PB[0] + m4['px'] - m4['R'], PB[1] + m4['py'] - m4['R'], PB[0] + m4['px'] + m4['R'],
                   PB[1] + m4['py'] + m4['R']]
            dashed_ellipse(ImageDraw.Draw(ghost), box, 2, (255, 255, 255, 150), n=24)
            frame.alpha_composite(ghost)
        frame.alpha_composite(over)
        # star labels in panel B (M4's fades out while it overlaps M2)
        paste(frame, m2_lab, (PB[0] + m2['px'], PB[1] + m2['py'] + m2['R'] + 14), anchor='mt')
        a4 = float(np.clip((dist - (m2['R'] + m4['R'])) / 30 + 1, 0, 1)) if dist < m2['R'] + m4['R'] + 30 else 1.0
        paste(frame, m4_lab, (PB[0] + m4['px'], PB[1] + m4['py'] + m4['R'] + 14), alpha=a4, anchor='mt')
        paste(frame, slow, (PB[0] + 18, PB[1] + 86), alpha=float(w_ecl[i]), anchor='lt')
        paste(frame, lined, (PA[0] + 18, PA[1] + 450), alpha=float(w_ecl[i]), anchor='lt')
        for k, key in enumerate(('out', 'deep', 'shallow', 'end')):
            paste(frame, captions[key], (PB[0] + B_C[0], PB[3] - 22), alpha=float(state[i, k]), anchor='mb')
        # brightness curve, revealed up to now
        t_now = min(ph, 2.0) * P_DAYS
        x_now = cx(t_now)
        frame.alpha_composite(curve.crop((0, 0, int(np.ceil(x_now)) + 2, H)), (0, 0))
        for p_dip, lab in dips:
            if ph >= p_dip + 0.02:
                a = float(np.clip((ph - p_dip - 0.02) / 0.03, 0, 1))
                ymin = 100 * flux(p_dip, j)
                paste(frame, lab, (cx(p_dip * P_DAYS) + 66, cy(ymin)), alpha=a, anchor='lm')
        for p1, lay, txt, (tx, ty) in br_layers:
            if ph >= p1 + 0.04:
                a = float(np.clip((ph - p1 - 0.04) / 0.04, 0, 1))
                paste(frame, lay, (0, 0), alpha=a)
                paste(frame, txt, (tx, ty), alpha=a, anchor='mb')
        if True:
            y_now = cy(100 * flux(min(ph, 2.0), j))
            marker = Image.new('RGBA', (W, H), (0, 0, 0, 0))
            dm = ImageDraw.Draw(marker)
            dm.line([(x_now, YT), (x_now, YB)], fill=BLUE + (70,), width=1)
            dm.ellipse([x_now - 9, y_now - 9, x_now + 9, y_now + 9], fill=(255, 255, 255, 60))
            dm.ellipse([x_now - 6, y_now - 6, x_now + 6, y_now + 6], fill=INK + (255,), outline=PANEL + (255,),
                       width=2)
            frame.alpha_composite(marker)
        return frame.convert('RGB')

    if len(sys.argv) > 2 and sys.argv[1] == '--preview':
        os.makedirs(sys.argv[2], exist_ok=True)
        for ph_want in (0.12, 0.2, 0.25, 0.75, 1.62, 2.0):
            i = int(np.argmin(np.abs(phases - ph_want))) if ph_want < 2 else n - 1
            make_frame(i).save(os.path.join(sys.argv[2], f'preview_{ph_want:.2f}.png'))
            print('preview', ph_want)
        return

    ffmpeg = shutil.which('ffmpeg')
    if ffmpeg is None:
        import imageio_ffmpeg
        ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
    gif_dir = tempfile.mkdtemp(prefix='tic71316387_gif_')
    os.makedirs(OUT, exist_ok=True)
    mp4 = os.path.join(OUT, 'tic71316387_eclipses.mp4')
    proc = subprocess.Popen([ffmpeg, '-y', '-loglevel', 'error', '-f', 'rawvideo', '-pix_fmt', 'rgb24',
                             '-s', f'{W}x{H}', '-r', str(FPS), '-i', '-', '-c:v', 'libx264', '-preset', 'slow',
                             '-crf', '18', '-pix_fmt', 'yuv420p', '-movflags', '+faststart', mp4],
                            stdin=subprocess.PIPE)
    poster_at = int(np.argmin(np.abs(phases - 0.25)))
    for i in range(n):
        rgb = make_frame(i)
        proc.stdin.write(rgb.tobytes())
        if i % (FPS // GIF_FPS) == 0:             # lossless, downscaled frames for the GIF
            rgb.resize((GIF_W, GIF_W * H // W), Image.LANCZOS).save(
                os.path.join(gif_dir, f'{i // (FPS // GIF_FPS):04d}.png'), compress_level=1)
        if i == poster_at:
            rgb.save(os.path.join(OUT, 'tic71316387_eclipses_poster.png'), optimize=True)
        if i % 60 == 0:
            print(f'  frame {i}/{n}', flush=True)
    proc.stdin.close()
    proc.wait()
    print('wrote', mp4)
    make_gif(ffmpeg, gif_dir, os.path.join(OUT, 'tic71316387_eclipses.gif'))
    shutil.rmtree(gif_dir)


def make_gif(ffmpeg, frames_dir, gif):
    """Two-pass GIF from lossless frames: one palette from whole frames (so the static blue curve keeps its
    colour), no dithering (unchanged pixels stay identical between frames) and only the changed pixels stored."""
    pattern = os.path.join(frames_dir, '%04d.png')
    palette = os.path.join(frames_dir, 'palette.png')
    subprocess.run([ffmpeg, '-y', '-loglevel', 'error', '-framerate', str(GIF_FPS), '-i', pattern, '-vf',
                    f'palettegen=max_colors={GIF_COLOURS}:stats_mode=full', palette], check=True)
    subprocess.run([ffmpeg, '-y', '-loglevel', 'error', '-framerate', str(GIF_FPS), '-i', pattern, '-i', palette,
                    '-lavfi', 'paletteuse=dither=none:diff_mode=rectangle', '-gifflags', '+transdiff',
                    '-loop', '0', gif], check=True)
    print('wrote', gif, f'({os.path.getsize(gif) / 1e6:.1f} MB)')


if __name__ == '__main__':
    main()
