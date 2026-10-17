import re, os, sys, json
T='tic214311908/analysis'
W=T+'/writeup'
man=[l.strip() for l in open(W+'/bundle_manifest.txt') if l.strip() and not l.startswith('#')]
# map to published paths
pub=set()
for p in man:
    if p.startswith('writeup/'): pub.add(p[len('writeup/'):])
    else: pub.add('analysis/'+p)
# expected figures
for f in ['false_positives','fold','followup_2027','host','localisation','per_transit','scenarios']:
    pub.add('figures/%s.png'%f)
print('n published', len(pub))
txt=open(W+'/README.md').read()
links=re.findall(r'\]\(([^)]+)\)', txt)
for l in links:
    if l.startswith('http'): continue
    if l.startswith('../'):
        rp='/home/user/Space/tess-search/'+l[3:]
        print('REPO', l, os.path.exists(rp)); continue
    if l.endswith('/'):
        ok=any(p.startswith(l) for p in pub)
    else:
        ok=l in pub
    print('OK ' if ok else 'MISSING', l)
# check key_numbers sources
kn=json.load(open(W+'/key_numbers.json'))
miss=set()
for e in kn['entries']:
    for s in e['source']:
        if s.startswith('writeup/'):
            q=s[len('writeup/'):]
        else: q='analysis/'+s
        if q not in pub: miss.add(s)
        if not os.path.exists(T+'/'+s): print('NOT ON DISK', s)
print('key_numbers sources not in bundle:', sorted(miss))
# corrections_to_repo links (relative to tic214311908/ for ../tic214311908 style)
c=open(W+'/corrections_to_repo.md').read()
for l in sorted(set(re.findall(r'\]\(([^)]+)\)', c))): print('CORR link', l)
# captions sources exist?
cap=open(W+'/figures/captions.md').read()
for s in sorted(set(re.findall(r'`([a-z0-9_]+/[^`]+)`', cap))):
    print('CAP src', s, os.path.exists(T+'/'+s), ('analysis/'+s) in pub)
