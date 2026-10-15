import subprocess, os, sys, time
D='/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/astro/data/svo/'
filters=['GALEX/GALEX.NUV','GAIA/GAIA3.G','GAIA/GAIA3.Gbp','GAIA/GAIA3.Grp','GAIA/GAIA3.Grvs','2MASS/2MASS.J','2MASS/2MASS.H','2MASS/2MASS.Ks',
'WISE/WISE.W1','WISE/WISE.W2','WISE/WISE.W3','WISE/WISE.W4','SkyMapper/SkyMapper.u','SkyMapper/SkyMapper.v','SkyMapper/SkyMapper.g','SkyMapper/SkyMapper.r',
'SkyMapper/SkyMapper.i','SkyMapper/SkyMapper.z','Generic/Johnson.B','Generic/Johnson.V','SLOAN/SDSS.gprime_filter','SLOAN/SDSS.rprime_filter','SLOAN/SDSS.iprime_filter',
'TYCHO/TYCHO.B','TYCHO/TYCHO.V','DENIS/DENIS.I','TESS/TESS.Red','Generic/Cousins.I','Generic/Cousins.R','Paranal/VISTA.Ks','Paranal/VISTA.J']
for f in filters:
    fn=D+'filters/'+f.replace('/','__')+'.xml'
    if os.path.exists(fn) and os.path.getsize(fn)>1000: continue
    for sysn in ['Vega','AB']:
        fn2=fn.replace('.xml',f'_{sysn}.xml')
        subprocess.run(['curl','-sS','--max-time','60','-o',fn2,f'http://svo2.cab.inta-csic.es/theory/fps/fps.php?ID={f}&PhotCalID={f}/{sysn}'])
    print(f,flush=True)
if len(sys.argv)>1:
    # BT-Settl CIFIST spectra: list then download
    import re
    for logg in [4.5,5.0]:
        lst=subprocess.run(['curl','-sS','--max-time','90',f'http://svo2.cab.inta-csic.es/theory/newov2/ssap.php?model=bt-settl-cifist&teffmin=2600&teffmax=7000&loggmin={logg}&loggmax={logg}'],capture_output=True,text=True).stdout
        rows=re.findall(r'<TD>(\d+)</TD><TD>([\d.]+)</TD><TD>0</TD>.*?fid=(\d+)\]\]',lst)
        print(logg,len(rows),flush=True)
        for teff,lg,fid in rows:
            fn=D+f'btsettl/bts_{int(teff):05d}_{float(lg):.1f}.xml'
            if os.path.exists(fn) and os.path.getsize(fn)>10000: continue
            subprocess.run(['curl','-sS','--max-time','120','-o',fn,f'http://svo2.cab.inta-csic.es/theory/newov2/ssap.php?model=bt-settl-cifist&fid={fid}'])
        print('done',logg,flush=True)
