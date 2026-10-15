for f in tess2018222125942-s0001-1-2-0120-s_ffic.fits tess2018215152942-s0001-1-2-0120-s_ffic.fits \
 tess2020217211915-s0028-1-1-0190-s_ffic.fits \
 tess2023223181152-s0068-1-1-0262-s_ffic.fits tess2023215212152-s0068-1-1-0262-s_ffic.fits \
 tess2025218053243-s0095-1-2-0292-s_ffic.fits tess2025212112924-s0095-1-2-0292-s_ffic.fits \
 tess2026177080230-s0105-2-3-0307-s_ffic.fits tess2026171095550-s0105-2-3-0307-s_ffic.fits \
 tess2026204174229-s0106-1-4-0308-s_ffic.fits tess2026199045549-s0106-1-4-0308-s_ffic.fits; do
  [ -s $f ] || timeout 300 curl -sS -L -o $f "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:TESS/product/$f"; echo "$f $(stat -c %s $f)"; done
