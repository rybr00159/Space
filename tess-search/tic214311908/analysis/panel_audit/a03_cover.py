import re, json
W='tic214311908'
txt=open(W+'/README.md').read()
kn=json.load(open(W+'/key_numbers.json'))['entries']
alltext=' '.join(e['readme_text'] for e in kn)
missing=[]
for i,line in enumerate(txt.split('\n'),1):
    if line.startswith('![') : continue
    for m in re.finditer(r'[−+-]?\d[\d,]*\.?\d*', line):
        tok=m.group(0)
        core=tok.lstrip('−+-')
        if core in alltext: continue
        missing.append((i,tok,line[max(0,m.start()-30):m.end()+20]))
for x in missing: print(x)
print(len(missing))
