"""Builds the supporting PDF report for the ExoFOP-TESS observing note on TIC 71316387 (CTOI 71316387.01).
Run with the system python (reportlab):  python3 build_report.py
All numbers are taken from the refereed campaign results (exofop_leads/TIC71316387_campaign/)."""
import json
import os
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT, TA_CENTER
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas as rl_canvas
from reportlab.platypus import (BaseDocTemplate, PageTemplate, Frame, Paragraph, Spacer, Table, TableStyle, Image,
                                KeepTogether, NextPageTemplate, CondPageBreak)
from reportlab.lib.fonts import addMapping

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, 'figs')
# ExoFOP file-name pattern: <target>O-<initials><UT date>-<info>.pdf ('xx' = uploader's initials, to be filled in)
OUTNAME = os.environ.get('REPORT_NAME', 'TIC71316387.01O-xx20261016-EB-analysis.pdf')
OUT = os.path.join(HERE, OUTNAME)
VERSION = '1.0'
DATE = '2026-10-16'
AUTHOR = os.environ.get('REPORT_AUTHOR', '[Author name]')
AUTHOR_LINE = os.environ.get('REPORT_AUTHOR_LINE', '[Author name] (ExoFOP-TESS user [username])')

# ------------------------------------------------------------------ fonts
LIB = '/usr/share/fonts/truetype/liberation'
DJV = '/usr/share/fonts/truetype/dejavu'
for name, f in (('Serif', 'LiberationSerif-Regular'), ('Serif-Bold', 'LiberationSerif-Bold'),
                ('Serif-Italic', 'LiberationSerif-Italic'), ('Serif-BoldItalic', 'LiberationSerif-BoldItalic'),
                ('Sans', 'LiberationSans-Regular'), ('Sans-Bold', 'LiberationSans-Bold'),
                ('Sans-Italic', 'LiberationSans-Italic'), ('Sans-BoldItalic', 'LiberationSans-BoldItalic')):
    pdfmetrics.registerFont(TTFont(name, f'{LIB}/{f}.ttf'))
pdfmetrics.registerFont(TTFont('DJ', f'{DJV}/DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('DJ-Bold', f'{DJV}/DejaVuSans-Bold.ttf'))
addMapping('Serif', 0, 0, 'Serif'); addMapping('Serif', 1, 0, 'Serif-Bold')
addMapping('Serif', 0, 1, 'Serif-Italic'); addMapping('Serif', 1, 1, 'Serif-BoldItalic')
addMapping('Sans', 0, 0, 'Sans'); addMapping('Sans', 1, 0, 'Sans-Bold')
addMapping('Sans', 0, 1, 'Sans-Italic'); addMapping('Sans', 1, 1, 'Sans-BoldItalic')
addMapping('DJ', 0, 0, 'DJ'); addMapping('DJ', 1, 0, 'DJ-Bold'); addMapping('DJ', 0, 1, 'DJ'); addMapping('DJ', 1, 1, 'DJ-Bold')

INK = colors.HexColor('#0b0b0b')
INK2 = colors.HexColor('#52514e')
RULE = colors.HexColor('#b9b7b1')
BLUE = colors.HexColor('#2a78d6')
NAVY = colors.HexColor('#173a66')
TINT = colors.HexColor('#f1f5fb')
GREYBG = colors.HexColor('#f4f4f2')

# glyphs missing from Liberation are drawn from DejaVu Sans
SPECIAL = {'☉': '☉', '≲': '≲', '≳': '≳', '⊕': '⊕', '★': '★'}


def g(s):
    for ch in SPECIAL:
        s = s.replace(ch, f'<font name="DJ">{ch}</font>')
    return s


def Msun():
    return 'M<sub>☉</sub>'


def Rsun():
    return 'R<sub>☉</sub>'


# ------------------------------------------------------------------ styles
BODY = ParagraphStyle('body', fontName='Serif', fontSize=9.6, leading=13.0, textColor=INK, alignment=TA_JUSTIFY,
                      spaceAfter=5.5)
BODYL = ParagraphStyle('bodyl', parent=BODY, alignment=TA_LEFT)
BUL = ParagraphStyle('bul', parent=BODY, leftIndent=12, bulletIndent=2, spaceAfter=3.2, alignment=TA_LEFT,
                     bulletFontName='Serif')
H1 = ParagraphStyle('h1', fontName='Sans-Bold', fontSize=11.2, leading=14, textColor=NAVY, spaceBefore=9, spaceAfter=4.5)
H2 = ParagraphStyle('h2', fontName='Sans-Bold', fontSize=9.4, leading=12, textColor=INK, spaceBefore=5, spaceAfter=2.5)
CAP = ParagraphStyle('cap', fontName='Sans', fontSize=7.9, leading=10.2, textColor=INK2, alignment=TA_JUSTIFY, spaceBefore=3,
                     spaceAfter=9)
TH = ParagraphStyle('th', fontName='Sans-Bold', fontSize=7.8, leading=10.8, textColor=INK)
TD = ParagraphStyle('td', fontName='Sans', fontSize=7.8, leading=10.8, textColor=INK)
TDS = ParagraphStyle('tds', parent=TD, textColor=INK2)
TDC = ParagraphStyle('tdc', parent=TD, alignment=TA_CENTER)
SUMH = ParagraphStyle('sumh', fontName='Sans-Bold', fontSize=9.0, leading=11.5, textColor=NAVY, spaceAfter=3)
SUM = ParagraphStyle('sum', fontName='Serif', fontSize=9.4, leading=12.6, textColor=INK, leftIndent=11, bulletIndent=0,
                     spaceAfter=2.6, bulletFontName='Serif')
REF = ParagraphStyle('ref', fontName='Serif', fontSize=8.4, leading=10.6, textColor=INK, leftIndent=12, firstLineIndent=-12,
                     spaceAfter=1.6, alignment=TA_LEFT)
SMALL = ParagraphStyle('small', fontName='Sans', fontSize=7.6, leading=9.8, textColor=INK2, alignment=TA_LEFT)

PAGE_W, PAGE_H = letter
LM = RM = 0.85 * inch
TM, BM = 0.95 * inch, 0.80 * inch
TW = PAGE_W - LM - RM


def P(t, st=BODY):
    return Paragraph(g(t), st)


def bullets(items, st=BUL, sym='•'):
    return [Paragraph(g(t), st, bulletText=sym) for t in items]


def table(rows, widths, header=True, zebra=False, small_first=False, valign='TOP', pad=3.2, spans=()):
    data = []
    for i, r in enumerate(rows):
        st = TH if (header and i == 0) else TD
        data.append([c if not isinstance(c, str) else Paragraph(g(c), st) for c in r])
    t = Table(data, colWidths=widths, repeatRows=1 if header else 0, hAlign='LEFT')
    cmds = [('VALIGN', (0, 0), (-1, -1), valign),
            ('LEFTPADDING', (0, 0), (-1, -1), 4), ('RIGHTPADDING', (0, 0), (-1, -1), 4),
            ('TOPPADDING', (0, 0), (-1, -1), pad), ('BOTTOMPADDING', (0, 0), (-1, -1), pad),
            ('LINEABOVE', (0, 0), (-1, 0), 0.9, INK), ('LINEBELOW', (0, -1), (-1, -1), 0.9, INK)]
    if header:
        cmds += [('LINEBELOW', (0, 0), (-1, 0), 0.5, INK), ('BACKGROUND', (0, 0), (-1, 0), GREYBG)]
    if zebra:
        for i in range(1 if header else 0, len(rows)):
            if (i % 2 == 0):
                cmds.append(('BACKGROUND', (0, i), (-1, i), colors.HexColor('#fafaf8')))
    for i in range(1 if header else 0, len(rows) - 1):
        cmds.append(('LINEBELOW', (0, i), (-1, i), 0.25, colors.HexColor('#dedcd6')))
    cmds += [('SPAN', a, b) for a, b in spans]
    t.setStyle(TableStyle(cmds))
    return t


def tcap(text):
    return Paragraph(g(text), ParagraphStyle('tcap', parent=CAP, keepWithNext=1, spaceAfter=4))


def figure(path, caption, width=TW):
    from reportlab.lib.utils import ImageReader
    ir = ImageReader(path)
    w, h = ir.getSize()
    img = Image(path, width=width, height=width * h / w)
    img.hAlign = 'CENTER'
    return FloatBlock([Spacer(1, 7), img, Paragraph(g(caption), CAP)])


class FloatBlock(KeepTogether):
    """A figure or table that may move to the top of the next page while the text keeps flowing (like a LaTeX float)."""


class FlushFloats(Spacer):
    """Marks the point after which no deferred float may be held back any longer."""
    def __init__(self):
        Spacer.__init__(self, 1, 0)


class FloatDoc(BaseDocTemplate):
    def __init__(self, *a, **k):
        BaseDocTemplate.__init__(self, *a, **k)
        self._deferred, self._newpage = [], False

    def handle_pageBegin(self):
        BaseDocTemplate.handle_pageBegin(self)
        self._newpage = True

    def handle_flowable(self, flowables):
        if self._newpage and self._deferred:
            # Top of a fresh page: the held-back floats go first. The first is placed unconditionally (so the
            # layout always advances); the others stay floats and may move on again if they do not fit. The text
            # that was waiting for this page may now land later, so clear its 'postponed' mark.
            if flowables and hasattr(flowables[0], '_postponed'):
                del flowables[0]._postponed
            first, rest = self._deferred[0], self._deferred[1:]
            flowables[0:0] = [KeepTogether(first._content)] + rest
            self._deferred = []
        self._newpage = False
        f = flowables[0]
        if isinstance(f, FlushFloats) and self._deferred:
            flowables[0:1] = [KeepTogether(d._content) for d in self._deferred] + [Spacer(1, 0)]
            self._deferred = []
            return
        if isinstance(f, FloatBlock):
            fr = self.frame
            avail = fr._y - fr._y1p
            w, h = f.wrapOn(self.canv, fr._aW, avail)
            if self._deferred or h > avail:
                self._deferred.append(flowables.pop(0))
                return
            flowables[0] = KeepTogether(f._content)
        BaseDocTemplate.handle_flowable(self, flowables)


# ------------------------------------------------------------------ page decoration
class NumberedCanvas(rl_canvas.Canvas):
    def __init__(self, *a, **k):
        rl_canvas.Canvas.__init__(self, *a, **k)
        self._saved = []

    def showPage(self):
        self._saved.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        n = len(self._saved)
        for st in self._saved:
            self.__dict__.update(st)
            self.draw_footer(n)
            rl_canvas.Canvas.showPage(self)
        rl_canvas.Canvas.save(self)

    def draw_footer(self, n):
        self.saveState()
        self.setStrokeColor(RULE); self.setLineWidth(0.5)
        self.line(LM, BM - 0.30 * inch, PAGE_W - RM, BM - 0.30 * inch)
        self.setFont('Sans', 7.2); self.setFillColor(INK2)
        self.drawString(LM, BM - 0.46 * inch, f'TIC 71316387 / CTOI 71316387.01  ·  supporting report v{VERSION}  ·  {DATE}')
        self.drawRightString(PAGE_W - RM, BM - 0.46 * inch, f'Page {self._pageNumber} of {n}')
        self.restoreState()


def first_page(c, doc):
    c.saveState()
    c.setFillColor(NAVY)
    c.rect(0, PAGE_H - 0.20 * inch, PAGE_W, 0.20 * inch, stroke=0, fill=1)
    c.restoreState()


def later_pages(c, doc):
    c.saveState()
    c.setFillColor(NAVY)
    c.rect(0, PAGE_H - 0.09 * inch, PAGE_W, 0.09 * inch, stroke=0, fill=1)
    c.setFont('Sans', 7.4); c.setFillColor(INK2)
    c.drawString(LM, PAGE_H - 0.55 * inch, RUNNING_TITLE)
    c.drawRightString(PAGE_W - RM, PAGE_H - 0.55 * inch, 'ExoFOP-TESS supporting report')
    c.setStrokeColor(RULE); c.setLineWidth(0.5)
    c.line(LM, PAGE_H - 0.62 * inch, PAGE_W - RM, PAGE_H - 0.62 * inch)
    c.restoreState()


# ------------------------------------------------------------------ content
story = []
RUNNING_TITLE = 'TIC 71316387 (CTOI 71316387.01): eclipsing binary, probably an M-dwarf pair in a hierarchical triple'

# ---- title block
story.append(Paragraph('EXOFOP-TESS OBSERVING NOTE  ·  SUPPORTING REPORT', ParagraphStyle(
    'kick', fontName='Sans-Bold', fontSize=7.8, leading=10, textColor=BLUE, spaceAfter=6)))
story.append(Paragraph(g('CTOI 71316387.01 is an eclipsing binary at twice the catalogued period, '
                         'most likely an M-dwarf pair bound to the G star TIC 71316387'), ParagraphStyle(
    'title', fontName='Sans-Bold', fontSize=15, leading=19, textColor=INK, spaceAfter=7)))
story.append(Paragraph(g('Light-curve, pixel-level, colour and Gaia analysis of public data, with predictions for '
                         'high-resolution imaging, spectroscopy and multi-colour photometry'), ParagraphStyle(
    'sub', fontName='Serif-Italic', fontSize=10.5, leading=13.5, textColor=INK2, spaceAfter=9)))

meta = [
    ['Target', g('TIC 71316387 = TYC 8001-154-1 = Gaia DR3 6571343840265900544 = 2MASS J22033305-4250076'),
     'Prepared by', g(AUTHOR_LINE)],
    ['Candidate', 'CTOI 71316387.01 (listed P = 0.6518369 d, 5,495 ppm; user disposition PC)',
     'Date / version', f'{DATE} / v{VERSION}'],
    ['Data', 'TESS Sectors 1, 28, 68, 95, 105, 106; Gaia DR3; SuperWASP; NEOWISE; VISTA VHS', 'Status',
     'Independent analysis of public data; no new observations'],
]
mrows = [[Paragraph(a, ParagraphStyle('mk', fontName='Sans-Bold', fontSize=7.6, leading=9.6, textColor=INK2)),
          Paragraph(b, ParagraphStyle('mv', fontName='Sans', fontSize=7.8, leading=9.8, textColor=INK)),
          Paragraph(c, ParagraphStyle('mk', fontName='Sans-Bold', fontSize=7.6, leading=9.6, textColor=INK2)),
          Paragraph(d, ParagraphStyle('mv', fontName='Sans', fontSize=7.8, leading=9.8, textColor=INK))] for a, b, c, d in meta]
mt = Table(mrows, colWidths=[0.72 * inch, 3.50 * inch, 0.92 * inch, TW - 5.14 * inch], hAlign='LEFT')
mt.setStyle(TableStyle([('VALIGN', (0, 0), (-1, -1), 'TOP'), ('LEFTPADDING', (0, 0), (-1, -1), 0),
                        ('RIGHTPADDING', (0, 0), (-1, -1), 6), ('TOPPADDING', (0, 0), (-1, -1), 1.6),
                        ('BOTTOMPADDING', (0, 0), (-1, -1), 1.6),
                        ('LINEABOVE', (0, 0), (-1, 0), 0.9, INK), ('LINEBELOW', (0, -1), (-1, -1), 0.5, RULE)]))
story += [mt, Spacer(1, 10)]

# ---- summary box
summary_items = [
    '<b>Not a planet.</b> The signal is an eclipsing binary (EB) with P = 1.30367228 d, twice the CTOI period; the '
    'CTOI\'s listed epoch falls on a secondary eclipse. Alternate eclipses differ by a factor of about two (5,770 and '
    '2,600 ppm), with the deeper eclipse on the same parity in all six TESS sectors from 2018 to 2026.',
    '<b>The G star is not an eclipsing component.</b> Its Gaia DR3 radial velocity is constant over 15 RVS field-of-view '
    'transits spanning 946 d (χ<super>2</super> p = 0.28), whereas any eclipsing partner of 0.1 ' + Msun() +
    ' or more in a 1.30 d orbit would move it by at least 18 km s<super>−1</super>.',
    '<b>The eclipses come from red light within 1.24 arcsec of the G star (95%).</b> TESS difference-image centroids '
    'show offsets of 2.1–3.7 arcsec, but these are a colour effect of the TESS optics. Once it is calibrated with field '
    'stars, all six sectors agree, and the size of the shift itself gives BP−RP ≈ 2.3 for the eclipsed star. Eclipse '
    'depths from SuperWASP (2006–07) and NEOWISE W1 independently support a red eclipsed star (shallower in the optical, '
    'deeper at 3.4 µm; together 3.2σ against the most favourable G/F case).',
    '<b>Probably a hierarchical triple.</b> Under main-sequence assumptions the eclipsing pair is an early-M + mid-M '
    'binary (about 0.45 + 0.20 ' + Msun() + ') supplying about 2.6% of the TESS flux. A physically bound companion is '
    'favoured over a chance-aligned background or foreground EB by estimated odds of order 10<super>4</super> '
    '(5×10<super>3</super>–1.4×10<super>4</super> for baseline assumptions; about 200 in the most pessimistic case). Gaia '
    'image diagnostics show an unresolved companion within about 1.2 arcsec; this constrains the separation, not the association.',
    '<b>Not yet confirmed.</b> There is no resolved image, common-proper-motion measurement or radial velocity of the '
    'M dwarfs. Near-infrared AO imaging is the most direct test (predicted ΔK<sub>s</sub> ≈ 2.4–2.8 mag, within about '
    '0.9 arcsec and possibly inside 0.1 arcsec). Section 7 lists further tests and upcoming eclipse times.',
]
sbox = [Paragraph('SUMMARY', SUMH)] + [Paragraph(g(t), SUM, bulletText='•') for t in summary_items]
sbox.append(Spacer(1, 2))
sbox.append(Paragraph(g('<b>Suggested disposition:</b> false positive (blended eclipsing binary). The EB lies within '
                        '1.24 arcsec of the target and is most likely a physically bound companion in a hierarchical '
                        'triple, not a background star.'),
                      ParagraphStyle('disp', parent=SUM, leftIndent=0, spaceBefore=2)))
bt = Table([[sbox]], colWidths=[TW])
bt.setStyle(TableStyle([('BACKGROUND', (0, 0), (-1, -1), TINT), ('LINEBEFORE', (0, 0), (0, -1), 2.4, BLUE),
                        ('LEFTPADDING', (0, 0), (-1, -1), 11), ('RIGHTPADDING', (0, 0), (-1, -1), 11),
                        ('TOPPADDING', (0, 0), (-1, -1), 8), ('BOTTOMPADDING', (0, 0), (-1, -1), 8)]))
story += [bt, Spacer(1, 6)]

# ---- 1. target
story.append(P('1  Target and catalogued candidate', H1))
story.append(P(
    'TIC 71316387 is a T = 11.36 G dwarf at about 297 pc (Stassun et al. 2019; Gaia Collaboration 2023). The CTOI was '
    'submitted to ExoFOP-TESS on 2026-08-27 as a high-probability planet from the RAVEN pipeline (Lafarga et al. 2026), '
    'with P = 0.6518369 d, a depth of 5,495 ± 1,338 ppm, a duration of 1.46 h and an implied radius of 8.5 R<sub>⊕</sub>. '
    'Table 1 gives the properties of the star that are used below.'))
t1 = [
    ['Quantity', 'Value', 'Source'],
    ['RA, Dec (ICRS, epoch 2016.0)', '22h03m33.11s, −42°50′07.8″ (330.887948°, −42.835511°)', 'Gaia DR3'],
    ['G, BP−RP, T', '11.775, 0.786, 11.358', 'Gaia DR3; TIC v8'],
    ['V, J, K<sub>s</sub>', '11.857, 10.728, 10.341', 'APASS; 2MASS'],
    ['Parallax, distance', '3.3622 ± 0.0170 mas; 297.4 ± 1.5 pc (inverse parallax; 293.9 pc after the −0.041 mas '
     'zero-point correction of Lindegren et al. 2021)', 'Gaia DR3'],
    ['Proper motion (μ<sub>α*</sub>, μ<sub>δ</sub>)', '(+38.503, −3.134) ± 0.013 mas yr<super>−1</super>', 'Gaia DR3'],
    ['Radial velocity', '−6.07 ± 0.86 km s<super>−1</super>; constant over 15 field-of-view transits (12 visibility '
     'periods, 946 d; χ<super>2</super> p = 0.275)', 'Gaia DR3 RVS'],
    ['T<sub>eff</sub>, R<sub>★</sub>, M<sub>★</sub>, age', '≈5,880 K, 1.08 ' + Rsun() + ', 1.03 ' + Msun() + ', ≈5 Gyr',
     'Gaia DR3 GSP-Phot (T<sub>eff</sub>); FLAME (R, M, age)'],
    ['RUWE; NSS; PM anomaly', '1.026; no NSS solution; none (Tycho-Gaia, p = 0.46)', 'Gaia DR3'],
    ['Nearest Gaia neighbour', '18.2 arcsec', 'Gaia DR3'],
]
story += [FloatBlock([Spacer(1, 4), tcap('<b>Table 1.</b> Target properties. Magnitudes include the light of the eclipsing pair '
               '(about 1.4% in G, 2.6% in T and 8.5% in K<sub>s</sub>).'),
          table(t1, [1.75 * inch, 3.55 * inch, TW - 5.30 * inch], pad=2.6), Spacer(1, 8)])]

# ---- 2. ephemeris
story.append(P('2  The signal is an eclipsing binary at twice the CTOI period', H1))
story.append(P(
    'Folded at the CTOI period, the odd and even events have clearly different depths (Figure 1a). At P = 1.30367228 d, '
    'the deeper (primary) and shallower (secondary) eclipses fall on the same parity in every sector, and each keeps a '
    'stable depth from 2018 to 2026 (Figure 1b). The depth difference is significant at more than 15σ in independent '
    'checks, and the epoch listed for the CTOI falls on a secondary eclipse (phase 0.500). The secondary eclipse is at '
    'phase 0.5 (e cos ω = 0.0006 ± 0.0005), so the orbit is consistent with circular. Spot modulation at the orbital '
    'period in five of six sectors shows that at least one component, most likely the spotted primary, rotates '
    'synchronously, as expected for a tidally locked 1.3 d binary. The deep eclipse is also recovered in SuperWASP data '
    'from 2006–07 (Pollacco et al. 2006; Butters et al. 2010) at the predicted time to within a few minutes (−3.8 min, '
    '1σ range −7.5 to +0.9 min, relative to a TESS-only ephemeris), so a linear ephemeris holds over 20 years.'))
t2 = [
    ['Parameter', 'Value', 'Notes'],
    ['Orbital period P', '1.30367228 ± 0.00000019 d', '= 2.00000 × CTOI period'],
    ['Primary epoch T<sub>0</sub>', 'BJD<sub>TDB</sub> 2460165.670416 ± 0.000180', 'BTJD 3165.670416'],
    ['CTOI epoch (BJD 2459067.3272)', 'phase 0.500', 'a secondary eclipse of the EB'],
    ['Primary / secondary depth (TESS, diluted)', '5,770 ± 100 ppm / 2,600 ± 55 ppm',
     'ratio 0.45; ±8% between sectors and data products'],
    ['Eclipse duration T<sub>14</sub>', '1.61 h (both eclipses)', 'flat secondary bottom ≈36 min in the main-sequence solution'],
    ['e cos ω; e sin ω', '0.0006 ± 0.0005; 0.01 ± 0.01', 'circular within the errors'],
    ['Out-of-eclipse variability', 'spots at P<sub>EB</sub>, 280–480 ppm; ellipsoidal < 74 ppm (95%)', 'synchronised rotation'],
    ['Eclipse-timing limits', '|dP/dt| < 0.032 s yr<super>−1</super>; curvature < 0.077 min yr<super>−2</super>', 'no third-body signal detected'],
]
story += [FloatBlock([Spacer(1, 4), tcap('<b>Table 2.</b> Eclipse ephemeris and light-curve properties. Times come from TESS-SPOC FFI light '
               'curves (S1, S28 and S68; Caldwell et al. 2020) and from our own TESSCut photometry (S95, S105, S106; '
               'Brasseur et al. 2019), with time stamps re-barycentred to the target position. QLP time stamps were not '
               'used for timing.'),
          table(t2, [2.15 * inch, 2.40 * inch, TW - 4.55 * inch], pad=2.6), Spacer(1, 8)])]

story.append(figure(f'{FIG}/fig1_lightcurve.png', width=0.93 * TW, caption=
    '<b>Figure 1.</b> TESS photometry of TIC 71316387. Panels a, c and d use Sectors 68, 95, 105 and 106 (200-s cadence) '
    'in 6-min bins; panel b uses all six sectors. '
    '<b>(a)</b> Folded at the CTOI period of 0.65184 d: odd (blue) and even (orange) events have different depths. '
    '<b>(b)</b> Diluted primary and secondary depths per sector; dashed lines and bands mark the adopted values. '
    'Sectors 1 and 28 come from the 30- and 10-min TESS-SPOC FFI light curves. Their primary depths sit about 8% '
    'higher, while their secondaries agree within the ~5% errors; this sector/cadence systematic is covered by the '
    'adopted ±8% depth-scale uncertainty. '
    '<b>(c, d)</b> Primary and secondary eclipses at P<sub>EB</sub> = 1.30367 d with the adopted main-sequence eclipse '
    'model and its residuals. In that model the secondary eclipse is a total occultation with a flat bottom lasting '
    'about 36 min; shape-only fits also allow a partial occultation.'))
# ---- 3. G star not EB
story.append(P('3  The G star is not one of the eclipsing stars', H1))
story.append(P(
    'The 2,600 ppm secondary eclipse (45% of the primary depth) requires a self-luminous companion with roughly 45% of the '
    'eclipsed star\'s surface brightness in the TESS band. If the G star (1.03 ' + Msun() + ') were the eclipsed star, '
    'its partner would be a K dwarf of about 0.7 ' + Msun() + ', giving the G star a radial-velocity semi-amplitude of '
    'about 95 km s<super>−1</super>; even a 0.1 ' + Msun() + ' partner would give K ≥ 18 km s<super>−1</super>. Gaia DR3 '
    'RVS (Katz et al. 2023) measured the star in 15 field-of-view transits (12 visibility periods) over 946 d and found '
    'no variability: the robust peak-to-peak amplitude is 9.9 km s<super>−1</super>, and the constancy test '
    '(χ<super>2</super> p = 0.275) shows that this scatter is consistent with measurement noise. The colour of the '
    'eclipsed light (Section 5) excludes the G star independently. The eclipses must therefore come from a fainter star '
    'blended with it.'))

# ---- 4. localisation
story.append(P('4  Location of the eclipsing source', H1))
story.append(P(
    'We fitted the TESS pixel response function to stacked in-eclipse minus out-of-eclipse difference images in each of the '
    'six sectors, relative to a fit to the direct image. Injection-recovery tests at 0, 1, 3 and 5 arcsec calibrated the '
    'bias and scatter. The bias-corrected offsets of the deeper eclipse are 2.1–3.7 arcsec. In every sector they point '
    'radially away from that camera\'s boresight, with tangential components consistent with zero (Figure 2a). A fixed sky '
    'offset cannot explain them (χ<super>2</super> = 217 for 22 degrees of freedom). This is the signature of the lateral '
    'colour of the TESS optics: the image of a red star sits farther from the boresight than that of a blue star at the '
    'same position. We measured the effect directly from field stars in the same six sectors: redder stars are shifted '
    'away from the boresight by 1.52 ± 0.02 arcsec per magnitude of BP−RP (2,560 fits; 1.56 ± 0.02 in a wider sample of '
    '5,560 fits; 1.45–1.64 per sector), more steeply at the red end (Figure 3a).'))
story.append(P(
    'A joint fit with one colour shift per eclipse type, projected onto each sector\'s boresight direction and scaled by '
    'that sector\'s measured coefficient, describes all 12 offsets (deep and shallow eclipses in six sectors). A source '
    'exactly on the G star with only a colour term gives χ<super>2</super> = 20.8 for 22 degrees of freedom. The '
    'colour-corrected offset is (ΔE, ΔN) = (−0.17 ± 0.18, +0.53 ± 0.28) arcsec, consistent with zero (p = 0.09 in the '
    'baseline analysis; 0.03–0.69 across variants). The source lies within 1.24 arcsec of the G star at 95% confidence and '
    'within 1.48 arcsec at 99.7%, and separations of 1.5 arcsec or more are excluded at p = 0.002 (Figure 2b). Across the '
    '15 analysis variants that include the deep eclipses, the 95% limit ranges from 1.06 to 1.54 arcsec (the shallow '
    'eclipses alone give 2.46 arcsec), and 2 arcsec or more stays excluded at p ≤ 0.005 even with all errors inflated by '
    '1.5. Every catalogued Gaia neighbour (the nearest at 18.2 arcsec) is therefore excluded as the source. VISTA VHS '
    'K<sub>s</sub> images (McMahon et al. 2013) recover injected companions with ΔK<sub>s</sub> ≤ 2.5 beyond 0.6 arcsec and '
    'ΔK<sub>s</sub> ≤ 3.0 beyond 1.0 arcsec. Given the predicted ΔK<sub>s</sub> ≈ 2.4–2.8 (Table 4), this places the pair '
    'inside about 0.8–1.0 arcsec.'))
# ---- 5. colour
story.append(P('5  The eclipsed star is red', H1))
story.append(P('Three independent measurements show that the deeper eclipse comes from a red star (Figure 3).', BODY))
story += bullets([
    '<b>TESS centroid colour shift.</b> The deep-eclipse light is shifted by 2.94 ± 0.28 arcsec relative to the G star. '
    'The field-star relation steepens for red stars, so the linear coefficient does not apply here: stars at BP−RP 2.14 '
    'and 2.42 are shifted by 2.42 and 3.29 arcsec. Read off this empirical relation (Figure 3a), the measured shift gives '
    'BP−RP ≈ 2.31 (1σ range 2.22–2.40, with a systematic of about ±0.1–0.2 mag because only a few dozen calibration stars '
    'are this red), an M1.5–M2.5 dwarf. A K7 dwarf (BP−RP ≈ 1.65) would be shifted by about 1.35 arcsec, 5.7σ '
    '(statistical) below the measured shift.',
    '<b>SuperWASP, 2006–07.</b> The deep eclipse is detected at S/N 4.8; none of 3,000 trials at false periods reaches '
    'that S/N. Its depth relative to TESS is 0.42 ± 0.15, or 0.36 ± 0.11 (statistical only) in an independent '
    're-reduction. An eclipsed G/F star (5,600–7,000 K) would give 0.82–1.31.',
    '<b>NEOWISE/AllWISE W1 (3.4 µm; Mainzer et al. 2014).</b> The deep-eclipse depth is 3.7 ± 1.1 times the TESS depth '
    '(S/N 4.3, from 8 in-eclipse points in 8 events), against 0.70–1.35 for an eclipsed G/F star; W1 alone is 2.0σ above '
    'the most favourable G/F case.',
])
story.append(P(
    'Together, SuperWASP and W1 reject the most favourable G/F eclipsed star at p = 6×10<super>−4</super> (3.2σ). '
    'Combined with the TESS colour shift, G/F and early-K eclipsed stars are excluded. A late-K (K5–K7) eclipsed star is '
    'disfavoured by the TESS colour shift, but the archival depth ratios alone still allow it: they give '
    'T<sub>eff</sub> = 4,040 (+330/−280) K for the eclipsed star.'))
story.append(figure(f'{FIG}/fig2_localisation.png', width=0.95 * TW, caption=
    '<b>Figure 2.</b> Localisation of the eclipsing light relative to the G star (star symbol). <b>(a)</b> Bias-corrected '
    'but not colour-corrected PRF difference-image offsets of the deeper eclipse in each sector. Grey arrows (of arbitrary '
    'length) show the direction away from that camera\'s boresight: the offsets follow these directions rather than a '
    'common sky position, as expected from lateral colour. <b>(b)</b> Joint colour-corrected position: 68, 95 and 99.7% '
    'regions calibrated by injection, and the best fit (×). Orange circles mark VHS K<sub>s</sub> companion limits. The '
    'diamond marks an unexplained faint VHS residual (ΔK<sub>s</sub> ≈ 4.0–4.4) that is too faint to be the EB.'))

story.append(figure(f'{FIG}/fig3_colour.png', width=0.95 * TW, caption=
    '<b>Figure 3.</b> Colour of the eclipsed light. <b>(a)</b> TESS lateral colour: centroid shift away from the camera '
    'boresight against Gaia BP−RP for field stars (grey points, 5,560 fits; black squares are binned means), relative to '
    'a star of the G star\'s colour. The relation steepens above BP−RP ≈ 1.9, where calibration stars are few (21, 21 '
    'and 4 in the three reddest bins). The horizontal band is the measured shift of the deep-eclipse light, and the '
    'vertical band is the implied colour. <b>(b)</b> Deep-eclipse depth relative to TESS against wavelength. Orange '
    'points: SuperWASP (two reductions; filled: campaign, open: independent re-reduction) and NEOWISE W1. Blue line and '
    'band: per-band predictions for the adopted early-M + mid-M pair (g′, V, r′, TESS, z′, J, H, K<sub>s</sub>, W1), '
    '±15% in the optical and ±20% in the infrared; at V the band spans the model-dependent range 0.28–0.33. Grey boxes: '
    'range expected if a G/F star (5,600–7,000 K) were eclipsed.'))

# ---- 6. components and triple
story.append(P('6  The eclipsing pair and its relation to the G star', H1))
story.append(P('6.1  Eclipse model', H2))
story.append(P(
    'We fitted a physical two-star eclipse model with third light to all six sectors, with limb darkening limited to the '
    'range of the Claret (2017) TESS coefficients for cool dwarfs. We used MCMC, likelihood profiles and '
    'injection-recovery tests; an independent re-analysis with a separately written model (one of the AI referee passes '
    'described in Section 9) reproduced the key limit. The light-curve shape alone requires the pair to supply more than '
    '2.2% of the TESS flux (Δχ<super>2</super> = 9) but sets no useful upper limit; larger fractions are disfavoured only '
    'by the main-sequence requirement, the ellipsoidal limit, the colours and the spectral energy distribution (≤5.8–11.6%). '
    'Requiring the two stars to follow main-sequence relations (Pecaut &amp; Mamajek 2013; Mann et al. 2015, 2019) gives '
    'the solution in Table 3. Its main-sequence photometric distance, about 210–330 pc, is consistent with the G star\'s '
    '297 pc. That agreement is expected if the pair is bound, but most chance-aligned M dwarfs of this colour and '
    'brightness would match too, so it adds little evidence. The component masses are not dynamical.'))
t3 = [
    ['Parameter', 'Primary (deeper eclipse)', 'Secondary', 'Notes'],
    ['Mass', '0.45 ± 0.04 ' + Msun(), '0.20 ± 0.02 ' + Msun(),
     'main-sequence relations; metallicity-dependent ([M/H] −0.35 to +0.10, not propagated); colour alone suggests '
     '0.40–0.44 ' + Msun()],
    ['Radius', '0.48 ± 0.02 ' + Rsun(), '0.21 ± 0.01 ' + Rsun(), 'primary ≈6% larger than main sequence'],
    ['T<sub>eff</sub>; spectral type', '≈3,550 K; M2V (M1–M2.5V)', '≈3,160 K; M4–M4.5V', 'archival colours give 4,040 K'],
    ['Pair share of TESS flux', '2.6 ± 0.4% (plus ±5% relative depth-scale systematic); shape-only lower limit 2.2%', '',
     'secondary ≈10% of the pair'],
    ['Undiluted eclipse depths', '≈22% (18–26%)', '≈10% (8–12%)', 'secondary total only in the main-sequence solution'],
    ['Orbit', 'a = 4.35 ± 0.09 ' + Rsun() + ' (0.020 au); i = 89.1 (+0.6/−0.65)°', '', 'i depends on the main-sequence prior'],
    ['Predicted RV semi-amplitude', 'K<sub>1</sub> = 52 ± 4 km s<super>−1</super>', 'K<sub>2</sub> = 117 ± 3 km s<super>−1</super>', ''],
    ['v sin i (synchronised)', '18.6 ± 1.0 km s<super>−1</super>', '8.2 ± 0.5 km s<super>−1</super>', ''],
    ['Apparent magnitudes of the pair', 'T ≈ 15.3, G ≈ 16.4, I<sub>c</sub> ≈ 14.9–15.3, K<sub>s</sub> ≈ 13.0', '', ''],
]
story += [FloatBlock([Spacer(1, 4), tcap('<b>Table 3.</b> Adopted parameters of the eclipsing pair (main-sequence solution). Values that span '
               'two columns describe the pair or its orbit.'),
          table(t3, [1.50 * inch, 1.95 * inch, 1.25 * inch, TW - 4.70 * inch], pad=2.6,
                spans=[((1, 4), (2, 4)), ((1, 6), (2, 6)), ((1, 9), (2, 9))]), Spacer(1, 8)])]
story.append(P('6.2  Evidence for an unresolved companion in Gaia DR3', H2))
story.append(P(
    'Gaia DR3 (Lindegren et al. 2021; Riello et al. 2021) sees a double-peaked image in 41% of the target\'s scans '
    '(ipd_frac_multi_peak = 41) and flags BP and RP blending in 47 of 47 and 48 of 48 field-of-view transits, while RUWE '
    'stays normal (1.026). We calibrated these diagnostics against about 4,850 pairs resolved by speckle or AO imaging and '
    '6,100 non-detections (SOAR speckle, Kepler follow-up imaging and Robo-AO). Among 2,650 stars with G = 9–14 in which '
    'speckle or AO imaging found no companion, only 0.19% reach ipd_frac_multi_peak ≥ 30, and none (95% upper limit '
    '0.11%) matches all four of the target\'s diagnostics together: ipd_frac_multi_peak 30–50, RUWE < 1.4, no Gaia '
    'neighbour, and BP/RP blending in at least 90% of transits. This indicates an unresolved companion at roughly '
    '0.1–1.2 arcsec. With few calibration pairs, the lower end of that range is uncertain. A bound and a chance-aligned '
    'EB would produce the same signature, so it constrains the separation, not the association.'))
story.append(P('6.3  Bound companion or chance alignment?', H2))
story.append(P(
    'Gaia star counts give about 4×10<super>−6</super> unrelated stars of the right colour (BP−RP 2.1–2.5) and brightness '
    '(supplying 2.2–6.5% of the TESS light) within a circle of radius 1.24 arcsec. Standard companion statistics for a '
    'solar-type star (Raghavan et al. 2010; Moe &amp; Di Stefano 2017 give the same within about 20%) give a 2.6% chance '
    'of a bound companion of matching mass at 0.1–1.24 arcsec. Assuming equal EB occurrence in both cases, the prior '
    'odds are about 6.5×10<super>3</super> in favour of a bound companion. Folding in the position (TESS, VHS and the Gaia '
    'signature) and the light-fraction likelihoods leaves the estimated posterior odds of order 10<super>4</super>. '
    'Baseline variants give 5×10<super>3</super> to 1.4×10<super>4</super>, the 160 tested combinations of analysis '
    'choices give about 6×10<super>2</super> to 4×10<super>4</super>, and about 200 remains when every choice goes against '
    'the triple and a further factor-3 prior systematic is applied. These numbers depend on the adopted companion '
    'statistics, equal EB occurrence and the colour-corrected TESS localisation. They are a probabilistic argument, not a '
    'demonstration.'))
story.append(P(
    'If the pair is bound, its projected separation is within about 0.9 arcsec (95% upper bound; about 270 au). If the '
    'Gaia signature comes from the pair, the separation is about 0.15–0.85 arcsec (median 0.4 arcsec). If it comes from '
    'another source, such as the unexplained VHS residual (Section 8), the separation is about 0.03–0.5 arcsec. The '
    'position angle is not constrained. The outer orbit would have a<sub>out</sub> ≈ 30–650 au and '
    'P<sub>out</sub> ≈ 10<super>2</super>–10<super>4</super> yr; P<sub>out</sub>/P<sub>in</sub> of about '
    '3×10<super>4</super> or more makes it dynamically very stable.'))

# ---- 7. predictions
story.append(P('7  Predictions for follow-up', H1))
story.append(P(
    'Table 4 lists the tests that can confirm or refute the triple interpretation. High-resolution imaging is the most '
    'direct. The predicted contrasts and depths assume a main-sequence pair at 297 pc; the contrasts carry a systematic '
    'uncertainty of about ±0.2 mag, and the depths about ±15% (optical) and ±20% (near-infrared).'))
NB = ' '
t4 = [
    ['Technique', 'Prediction', 'What it tests'],
    ['Near-IR AO (J, H, K<sub>s</sub>), 8-m class',
     'Companion within about 0.9 arcsec (possibly inside 0.1 arcsec) with ΔK<sub>s</sub> ≈ 2.4–2.8, ΔH ≈ 2.6–2.9, '
     'ΔJ ≈ 2.9–3.2 mag',
     'Presence and separation of the EB. A non-detection to ΔK<sub>s</sub> ≈ 3.5 beyond 0.1 arcsec would put the EB '
     'inside 0.1 arcsec (under 2% of the posterior if the Gaia signature is the EB\'s, 20–30% if not) or refute the '
     'model, so AO should reach inside 0.1 arcsec. A companion with ΔK<sub>s</sub> ≤ 2.5 beyond 0.6 arcsec would '
     'contradict VHS'],
    ['Optical speckle (8-m)', f'Δm(832{NB}nm) ≈ 3.7–4.0; ΔI<sub>c</sub> ≈ 3.7–4.1; Δm(562{NB}nm) ≈{NB}5.1–5.6 mag',
     'Detectable at 832 nm beyond ≈0.15–0.2 arcsec; SOAR I-band is borderline'],
    ['Second imaging epoch (1–2 yr)', 'Relative proper motion < 2.5–4.3 mas yr<super>−1</super> (95%) if bound',
     'Common proper motion. None of 77 field M dwarfs in the colour window moves this slowly (likelihood ratio '
     'about 25–340)'],
    ['Red / near-IR spectroscopy', 'Main-sequence predictions: K<sub>1</sub> ≈ 52, K<sub>2</sub> ≈ 117 km '
     's<super>−1</super> at P = 1.30367 d; systemic velocity within a few km s<super>−1</super> of '
     '−6.1 km s<super>−1</super> if bound; G-star lines constant. The pair gives about 3% of the light at I and 7–8.5% '
     'at H/K<sub>s</sub>, so expect at best an SB1 from the M2 lines in the blend; K<sub>2</sub> needs high-S/N near-IR '
     'or AO-fed spectroscopy',
     'M-dwarf nature, mass ratio and association (systemic velocity is a weak test: likelihood ratio about 6–12)'],
    ['Multi-colour eclipse photometry', 'Deep / shallow depths (ppm):<br/>'
     f'g′ 1,150{NB}/{NB}200; V 1,600–1,900{NB}/{NB}370;<br/>'
     f'i′ 5,600–6,700{NB}/{NB}2,200–3,000; z′ 7,200{NB}/{NB}3,800;<br/>'
     f'J 12,100{NB}/{NB}8,900; K<sub>s</sub> 18,000{NB}/{NB}13,800',
     'Red, M-dwarf eclipsed star (a G/F star would give a V depth of 4,600 ppm or more). This does not test whether '
     'the pair is bound'],
    ['Gaia DR4 epoch photometry', 'Eclipses of about 3.6 mmag (primary) and 1.2 mmag (secondary) in G; detection likely '
     'marginal (few in-eclipse transits)', 'Possible independent ephemeris check'],
]
story += [FloatBlock([Spacer(1, 4), tcap('<b>Table 4.</b> Testable predictions.'),
          table(t4, [1.40 * inch, 2.85 * inch, TW - 4.25 * inch], pad=2.6), Spacer(1, 8)])]
story.append(Spacer(1, 4))
ew = json.load(open(os.path.join(HERE, 'eclipse_windows.json')))
rows = [['Epoch E', 'Mid-eclipse (UTC)', 'BJD (TDB)', 'Observable from']]
SN = {'CTIO': 'Cerro Tololo (CTIO)', 'SAAO': 'Sutherland (SAAO)', 'SSO': 'Siding Spring (SSO)'}
for r in ew:
    rows.append([str(r['epoch']), r['utc_mid'].replace('T', ' '), f"{r['bjd_tdb']:.5f}", ', '.join(SN[s] for s in r['sites'])])
story.append(FloatBlock([Spacer(1, 4), tcap(
    '<b>Table 5.</b> Primary eclipses between 2026-10-17 and 2027-01-20 that are fully observable at night from at least '
    'one of three southern sites (Sun below −12° and target above 30° from 1.3 h before to 1.3 h after mid-eclipse); the '
    'site that can observe each event is listed. No later event in this window meets these criteria, because the target '
    'sets early in the evening from December. UTC times are rounded to the nearest minute; the 1σ timing uncertainty is '
    '0.4 min. Other times follow from T<sub>0</sub> + E × P (Table 2); secondary eclipses fall at E + 0.5.'),
    table(rows, [0.65 * inch, 1.45 * inch, 1.35 * inch, TW - 3.45 * inch], pad=2.2), Spacer(1, 8)]))

# ---- 8. caveats
story.append(P('8  Caveats and open issues', H1))
story += bullets([
    '<b>No direct evidence of association.</b> The bound interpretation rests on chance-alignment statistics and '
    'consistency checks. There is no resolved image, common proper motion or systemic velocity.',
    '<b>Light fraction and masses depend on priors.</b> The light fraction, masses, radius ratio, inclination and '
    'photometric distance all rest on the main-sequence assumption. The light-curve shape alone gives only a lower limit '
    '(more than 2.2%); with e sin ω free, a uniform-prior fit prefers about 7%. Larger fractions are disfavoured by the '
    'main-sequence requirement, the ellipsoidal limit, the colours and the SED, not by the eclipse shape. Without the '
    'main-sequence term, the odds fall to about 10<super>3</super> (roughly 1×10<super>3</super>–4×10<super>3</super>).',
    '<b>Temperature tension.</b> The archival depth ratios give 4,040 (+330/−280) K for the eclipsed star, against about '
    '3,550 K from the eclipse model and the TESS colour shift (a 1.6–1.9σ difference). For this reason a late-K primary is '
    'disfavoured but not excluded.',
    '<b>Shallow-eclipse colours.</b> The W1 and SuperWASP depths of the shallow eclipse disagree with the model. Both come '
    'from sparse archival data and are sensitive to detrending.',
    '<b>Separation.</b> The lower edge depends on a sparse calibration of the Gaia image diagnostics, and the VHS limits '
    'are somewhat optimistic because the injected companions used the fitted PSF.',
    '<b>Unexplained VHS feature.</b> A faint K<sub>s</sub> residual at 1.12–1.14 arcsec, PA 345–350° '
    '(ΔK<sub>s</sub> ≈ 4.0–4.4), appears in all three VHS K<sub>s</sub> frames. It would supply only about 0.5% of the TESS '
    'light, too little to be the EB. It may be a fourth body or a fitting artefact. It lies inside the separation range '
    'implied by the Gaia image diagnostics and could itself produce that signature, in which case the EB would lie at '
    'about 0.03–0.5 arcsec (Section 6.3).',
    '<b>TRICERATOPS not used.</b> Run at the CTOI period on a single fold, the false-positive probabilities of TRICERATOPS '
    '(Giacalone et al. 2021) are driven by its priors and do not describe this system, so we do not quote them.',
])

# ---- 9. data, methods
story.append(FlushFloats())
story.append(P('9  Data, methods and provenance', H1))
story.append(P(
    'Every input is public. The TESS data (Ricker et al. 2015) are TESS-SPOC FFI light curves (Sectors 1, 28 and 68, at '
    '30-min, 10-min and 200-s cadence; Jenkins et al. 2016; Caldwell et al. 2020) and our own aperture photometry from '
    'TESSCut full-frame-image cutouts (Sectors 95, 105 and 106; all six sectors for the pixel analysis; Brasseur et al. '
    '2019). QLP light curves (Huang et al. 2020) were used only as a depth cross-check. We also used Gaia DR3 astrometry, '
    'photometry, RVS velocities and image parameters; the TIC v8; APASS and 2MASS photometry (Henden et al. 2016; '
    'Skrutskie et al. 2006); the SuperWASP public archive; NEOWISE and AllWISE single-exposure photometry (Wright et al. '
    '2010; Mainzer et al. 2014); and VISTA VHS, SkyMapper (Onken et al. 2019) and DESI Legacy Survey (Dey et al. 2019) '
    'images. The work had four parts: localisation, eclipse colours, eclipse modelling and timing, and stellar and '
    'statistical analysis. An adversarial referee pass recomputed the key number of each part, and a final referee '
    'independently recomputed the bound-versus-chance odds over 160 combinations of analysis choices. The wording of this '
    'report follows that final referee\'s corrections. Scripts and machine-readable results are available from the author '
    'on request.'))
story.append(P(
    '<b>AI assistance.</b> The analysis and this report were produced with an AI agent (Claude, by Anthropic) working '
    'under the author\'s direction. The referee passes described above were also carried out by AI agents. No professional '
    'astronomer has yet reviewed the results. Independent checks are welcome.'))

# ---- references
story.append(FlushFloats())
story.append(P('References', H1))
refs = [
    'Brasseur, C. E., Phillip, C., Fleming, S. W., Mullally, S. E., & White, R. L. 2019, Astrocut: Tools for creating cutouts of TESS images, ASCL, ascl:1905.007',
    'Butters, O. W., West, R. G., Anderson, D. R., et al. 2010, A&A, 520, L10 (SuperWASP public archive)',
    'Caldwell, D. A., Tenenbaum, P., Twicken, J. D., et al. 2020, RNAAS, 4, 201 (TESS-SPOC FFI light curves)',
    'Claret, A. 2017, A&A, 600, A30 (TESS limb darkening)',
    'Dey, A., Schlegel, D. J., Lang, D., et al. 2019, AJ, 157, 168 (DESI Legacy Imaging Surveys)',
    'Gaia Collaboration, Vallenari, A., Brown, A. G. A., et al. 2023, A&A, 674, A1',
    'Giacalone, S., Dressing, C. D., Jensen, E. L. N., et al. 2021, AJ, 161, 24 (TRICERATOPS)',
    'Henden, A. A., Templeton, M., Terrell, D., et al. 2016, VizieR Online Data Catalog, II/336 (APASS DR9)',
    'Huang, C. X., Vanderburg, A., Pál, A., et al. 2020, RNAAS, 4, 204 (QLP)',
    'Jenkins, J. M., Twicken, J. D., McCauliff, S., et al. 2016, Proc. SPIE, 9913, 99133E (SPOC)',
    'Katz, D., Sartoretti, P., Guerrier, A., et al. 2023, A&A, 674, A5 (Gaia DR3 radial velocities)',
    'Lafarga, M., Armstrong, D. J., Cui, K., et al. 2026, MNRAS, 548, stag512 (RAVEN; doi:10.1093/mnras/stag512)',
    'Lindegren, L., Klioner, S. A., Hernández, J., et al. 2021, A&A, 649, A2 (Gaia EDR3 astrometric solution)',
    'Mainzer, A., Bauer, J., Cutri, R. M., et al. 2014, ApJ, 792, 30 (NEOWISE)',
    'Mann, A. W., Feiden, G. A., Gaidos, E., Boyajian, T., & von Braun, K. 2015, ApJ, 804, 64',
    'Mann, A. W., Dupuy, T., Kraus, A. L., et al. 2019, ApJ, 871, 63',
    'McMahon, R. G., Banerji, M., Gonzalez, E., et al. 2013, The Messenger, 154, 35 (VHS)',
    'Moe, M., & Di Stefano, R. 2017, ApJS, 230, 15',
    'Onken, C. A., Wolf, C., Bessell, M. S., et al. 2019, PASA, 36, e033 (SkyMapper DR2)',
    'Pecaut, M. J., & Mamajek, E. E. 2013, ApJS, 208, 9',
    'Pollacco, D. L., Skillen, I., Collier Cameron, A., et al. 2006, PASP, 118, 1407 (SuperWASP)',
    'Raghavan, D., McAlister, H. A., Henry, T. J., et al. 2010, ApJS, 190, 1',
    'Ricker, G. R., Winn, J. N., Vanderspek, R., et al. 2015, JATIS, 1, 014003 (TESS)',
    'Riello, M., De Angeli, F., Evans, D. W., et al. 2021, A&A, 649, A3 (Gaia EDR3 photometry)',
    'Skrutskie, M. F., Cutri, R. M., Stiening, R., et al. 2006, AJ, 131, 1163 (2MASS)',
    'Stassun, K. G., Oelkers, R. J., Paegert, M., et al. 2019, AJ, 158, 138 (TIC v8)',
    'Wright, E. L., Eisenhardt, P. R. M., Mainzer, A. K., et al. 2010, AJ, 140, 1868 (WISE)',
]
story += [Paragraph(g(r.replace('&', '&amp;')), REF) for r in refs]
story.append(Spacer(1, 6))
story.append(P('This research uses data from the TESS mission, obtained from the MAST data archive at the Space Telescope '
               'Science Institute, and from the ESA Gaia mission, processed by the Gaia DPAC. It also uses the SuperWASP '
               'public archive, NEOWISE, the VISTA Hemisphere Survey, and the Exoplanet Follow-up Observing Program '
               '(ExoFOP; DOI 10.26134/ExoFOP5) website, operated by the California Institute of Technology under contract '
               'with NASA.', SMALL))

# ------------------------------------------------------------------ build
doc = FloatDoc(OUT, pagesize=letter, leftMargin=LM, rightMargin=RM, topMargin=TM, bottomMargin=BM,
                      title=RUNNING_TITLE,
                      author=AUTHOR, subject='Supporting report for an ExoFOP-TESS observing note',
                      keywords='TESS, CTOI 71316387.01, TIC 71316387, eclipsing binary, hierarchical triple, false positive',
                      creator='reportlab')
f1 = Frame(LM, BM, TW, PAGE_H - 0.55 * inch - BM, id='f1', leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
f2 = Frame(LM, BM, TW, PAGE_H - TM - BM, id='f2', leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
doc.addPageTemplates([PageTemplate(id='first', frames=[f1], onPage=first_page, autoNextPageTemplate='later'),
                      PageTemplate(id='later', frames=[f2], onPage=later_pages)])
doc.build(story, canvasmaker=NumberedCanvas)
print('wrote', OUT)
