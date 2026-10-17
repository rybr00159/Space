"""s09: fetch arXiv abstracts (export.arxiv.org API, read-only) for literature numbers used by the bd/ priors.
Usage: s09_fetch_abs.py id1 id2 ...  -> lit/abs_<id>.xml and printed title/journal/abstract"""
import sys, re, time, urllib.request, html
for aid in sys.argv[1:]:
    url = f'https://export.arxiv.org/api/query?id_list={aid}'
    t = urllib.request.urlopen(url, timeout=40).read().decode()
    open(f'tic214311908/analysis/referee3_bd/lit/abs_{aid}.xml', 'w').write(t)
    ent = t.split('<entry>')[1] if '<entry>' in t else ''
    ti = re.findall(r'<title>(.*?)</title>', ent, re.S); s = re.findall(r'<summary>(.*?)</summary>', ent, re.S)
    j = re.findall(r'<arxiv:journal_ref>(.*?)</arxiv:journal_ref>', ent, re.S); d = re.findall(r'<arxiv:doi>(.*?)</arxiv:doi>', ent, re.S)
    au = re.findall(r'<name>(.*?)</name>', ent)
    print('==', aid, '|', ' '.join(ti[0].split()) if ti else '', '|', au[:3], '|', j, d)
    print(' '.join(html.unescape(s[0]).split()) if s else 'NO SUMMARY'); print(flush=True)
    time.sleep(3)
