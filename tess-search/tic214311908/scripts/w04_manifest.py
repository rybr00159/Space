"""w04: writeup/bundle_manifest.txt - files to copy from tic214/<folder>/ into the repository as
tess-search/tic214311908/analysis/<folder>/, plus the write-up's own files (top level of tic214311908/).
Rules: each folder's README.md; scripts/ (.py, .sh, small .json/.txt configs); out/ files with extension
json/csv/txt/png/ecsv/md of at most 1.5 MB; figs/ or fig/ PNGs of at most 1.5 MB; small top-level json/csv deliverables.
Skipped: data/, logs/, lit/, ext/, ext_src/, explore/, lc/ (downloads, logs, literature, re-creatable inputs),
npz/npy/pkl/fits and anything larger than 1.5 MB. Also listed: the final review panels (panel_*: scripts and out/),
and the write-up's figures, figscripts (with its small data file), novelty check and review log. Checks: every
key_numbers.json source and every README link resolves to a listed file. Copy with scripts/w06_copy_bundle.py."""
import os, json, re
W = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # the write-up folder (writeup/ or tic214311908/)
B = os.path.dirname(W) if os.path.isdir(os.path.join(os.path.dirname(W), 'clean')) else os.path.join(W, 'analysis')
FOLDERS = ['phot', 'referee_phot', 'loc', 'referee_loc', 'star', 'referee_star', 'ctx', 'referee_ctx',
           'shape', 'referee2_shape', 'search', 'referee2_search', 'host', 'referee2_host', 'fpp', 'referee2_fpp',
           'clean', 'referee3_clean', 'bd', 'referee3_bd', 'followup', 'referee3_followup']
# the final review of the write-up (three AI reviewers; no README of their own, described in review_log.md):
# their top-level or scripts/ .py files and their out/ json/csv/txt files
PANELS = ['panel_systematics', 'panel_fpp', 'panel_audit']
PANEL_SKIP = {'w03_copy.py'}          # a scratch copy of writeup/scripts/w03_key_numbers.py, not part of the review
LIM = 1.5e6
OUT_EXT = {'.json', '.csv', '.txt', '.png', '.ecsv', '.md'}
SKIP_DIRS = {'data', 'logs', 'lit', 'ext', 'ext_src', 'explore', 'lc', '__pycache__'}
sel, skipped_big = [], []
for f in FOLDERS + PANELS:
    root = f'{B}/{f}'
    for dp, dns, fns in os.walk(root):
        rel_dir = os.path.relpath(dp, root)
        parts = [] if rel_dir == '.' else rel_dir.split(os.sep)
        if parts and parts[0] in SKIP_DIRS: dns[:] = []; continue
        dns[:] = [d for d in dns if d not in ('__pycache__',)]
        for fn in fns:
            p = os.path.join(dp, fn); rel = os.path.relpath(p, B); ext = os.path.splitext(fn)[1].lower(); sz = os.path.getsize(p)
            top = parts[0] if parts else ''
            ok = False
            if f in PANELS:
                ok = (fn not in PANEL_SKIP and ((not parts or top == 'scripts') and ext == '.py' or
                                                top == 'out' and ext in {'.json', '.csv', '.txt'})) and sz <= LIM
            elif not parts and fn == 'README.md': ok = True
            elif top == 'scripts' and ext in {'.py', '.sh', '.json', '.txt'}: ok = sz <= LIM
            elif top == 'out' and ext in OUT_EXT: ok = sz <= LIM
            elif top in ('figs', 'fig') and ext == '.png': ok = sz <= LIM
            elif not parts and ext in {'.json', '.csv'}: ok = sz <= LIM
            if ok: sel.append((rel, sz))
            elif (top in ('out', 'figs', 'fig', 'scripts', '') and ext in OUT_EXT | {'.py', '.sh'} and sz > LIM):
                skipped_big.append((rel, sz))
sel.sort()
# write-up files (top level of tic214311908/)
MINE = ['README.md', 'key_numbers.json', 'ephemeris_adopted.json', 'predicted_transits_2026_2028.csv',
        'predicted_transits_2027.csv', 'fpp_ladder_adopted.json', 'corrections_to_repo.md', 'review_log.md',
        'review_findings.json', 'bundle_manifest.txt']
def _walk(sub, exts):
    out = []
    for dp, dns, fns in os.walk(f'{W}/{sub}'):
        dns[:] = [d for d in dns if d != '__pycache__']
        out += [os.path.relpath(os.path.join(dp, fn), W) for fn in fns if fn.endswith(exts)]
    return sorted(out)
wfiles = MINE + _walk('observability', ('.csv', '.md')) + _walk('scripts', ('.py', '.json')) \
    + _walk('figures', ('.png', '.md')) + _walk('figscripts', ('.py', '.sh', '.npz')) \
    + _walk('novelty', ('.md', '.py', '.json', '.csv', '.txt', '.log', '.xml', '.tsv', '.html', '.dat'))
OTHER = sorted(d for d in os.listdir(W) if os.path.isdir(f'{W}/{d}') and d not in ('observability', 'scripts', 'figures', 'figscripts', 'novelty', '__pycache__'))
# ---- checks ----
listed = {r for r, _ in sel} | {f'writeup/{w}' for w in wfiles}
kn = json.load(open(f'{W}/key_numbers.json'))
src_missing = sorted({s for e in kn['entries'] for s in e['source'] if s not in listed})
T = open(f'{W}/README.md').read()
links = re.findall(r'\]\(([^)]+)\)', T)
link_missing = []
for l in links:
    if l.startswith('http') or l.startswith('../'): continue
    if l.startswith('analysis/'):
        tgt = l[len('analysis/'):]
        if tgt not in listed and not any(x.startswith(tgt.rstrip('/') + '/') for x in listed): link_missing.append(l)
    else:
        t = l.rstrip('/')
        if t not in wfiles and not any(w.startswith(t + '/') for w in wfiles): link_missing.append(l)
figs = sorted(set(re.findall(r'\]\((figures/[^)]+)\)', T)))
tot = sum(s for _, s in sel)
wsz = sum(os.path.getsize(f'{W}/{w}') for w in wfiles if os.path.exists(f'{W}/{w}'))
by_folder = {}
for r, s in sel:
    f = r.split('/')[0]; by_folder.setdefault(f, [0, 0]); by_folder[f][0] += 1; by_folder[f][1] += s
lines = []
lines.append('# Bundle manifest for tess-search/tic214311908/ (prepared 2026-10-17; nothing has been copied or committed)')
lines.append('#')
lines.append('# Source root: the tic214/ analysis folder of the session scratchpad (the folder that holds phot/, loc/, ..., writeup/).')
lines.append('# Copy each "<folder>/<path>" below to "tic214311908/analysis/<folder>/<path>", and each "writeup/<path>" to')
lines.append('# "tic214311908/<path>" (write-up, figures, figscripts, novelty check, observability tables, review log).')
lines.append('# Use scripts/w06_copy_bundle.py <dest>/tic214311908: it copies this list and rewrites the absolute working paths')
lines.append('# that the analysis scripts and READMEs hard-code (analysis root -> tic214311908/analysis, write-up -> tic214311908,')
lines.append('# other local folders -> <local>), then fails if any absolute /tmp path is left. Do not copy the analysis folders by hand:')
lines.append('# about 430 of their files contain the absolute scratch path. After copying, search the copy once more for absolute /tmp paths.')
lines.append('# The review panels (panel_*) have no README; review_log.md describes them.')
lines.append('# Figures referenced by the README: ' + ', '.join(figs) + '.')
if OTHER: lines.append('# Folders in writeup/ not listed: ' + ', '.join(d + '/' for d in OTHER) + '.')
lines.append('# Rules: README.md; scripts/ (.py .sh and small .json/.txt configs); out/ json/csv/txt/png/ecsv/md <= 1.5 MB;')
lines.append('# figs/ and fig/ PNGs <= 1.5 MB; small top-level json/csv deliverables. Skipped: data/, logs/, lit/, ext/,')
lines.append('# ext_src/, explore/, lc/, npz/npy/pkl/fits and anything > 1.5 MB (pixel cubes, injection arrays, full light curves);')
lines.append('# the scripts regenerate them from MAST/Gaia.')
lines.append('#')
lines.append(f'# Analysis files: {len(sel)} ({tot/1e6:.1f} MB); write-up files: {len(wfiles)} ({wsz/1e6:.2f} MB, with figures/);')
lines.append(f'# total {(tot + wsz)/1e6:.1f} MB (limit ~60 MB).')
lines.append(f'# Checks: key_numbers.json sources not in this list: {len(src_missing)}; README links not resolved: {len(link_missing)}.')
lines.append('#')
lines.append('# Per folder (files, MB):')
for f in FOLDERS + PANELS:
    n, s = by_folder.get(f, [0, 0]); lines.append(f'#   {f:18s} {n:4d} {s/1e6:6.2f}')
lines.append('#')
lines.append('# Larger files left out on purpose (> 1.5 MB), for reference:')
for r, s in sorted(skipped_big):
    lines.append(f'#   {r} ({s/1e6:.1f} MB)')
lines.append('')
lines.append('## write-up files (copy to tic214311908/)')
for w in wfiles: lines.append(f'writeup/{w}')
lines.append('')
lines.append('## analysis files (copy to tic214311908/analysis/)')
for r, s in sel: lines.append(r)
open(f'{W}/bundle_manifest.txt', 'w').write('\n'.join(lines) + '\n')
print(len(sel), 'files', round(tot/1e6, 2), 'MB; writeup', len(wfiles), round(wsz/1e6, 2), 'MB')
print('src missing', src_missing); print('link missing', link_missing)
print({f: (n, round(s/1e6, 2)) for f, (n, s) in by_folder.items()})
print('skipped big', len(skipped_big))
