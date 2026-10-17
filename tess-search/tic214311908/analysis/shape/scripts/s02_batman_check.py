"""s02: validation of batman (quadratic LD) for radius ratios k = 0.1-5 (k > 1 = occultation of a small star by a
larger body) against a brute-force 2-D integration of the limb-darkened disk; also timing per model call.
Output: out/batman_check.json"""
import numpy as np, batman, time, json
def brute(z, k, u1, u2, n=1500):
    x = np.linspace(-1, 1, n); X, Y = np.meshgrid(x, x); r2 = X**2 + Y**2; ins = r2 < 1
    mu = np.sqrt(np.clip(1 - r2, 0, 1)); I = (1 - u1*(1-mu) - u2*(1-mu)**2) * ins
    tot = I.sum(); out = []
    for zz in z:
        occ = (X - zz)**2 + Y**2 < k**2
        out.append(1 - (I*occ).sum()/tot)
    return np.array(out)
res = {}
for k in [0.1, 0.5, 0.95, 1.0, 1.5, 3.0, 5.0]:
    p = batman.TransitParams(); p.t0=0; p.per=10.; p.rp=k; p.a=30.; p.inc=90.; p.ecc=0.; p.w=90.; p.limb_dark='quadratic'; p.u=[0.3,0.3]
    # map z(t): z = a * sin(2 pi t/P) for inc=90
    z = np.linspace(0, 1+k+0.05, 25)
    t = np.arcsin(z/p.a)*p.per/(2*np.pi)
    m = batman.TransitModel(p, t)
    fb = m.light_curve(p); fr = brute(z, k, 0.3, 0.3)
    res[str(k)] = float(np.max(np.abs(fb-fr)))
    print(k, 'max |batman - brute| =', res[str(k)], 'min flux', fb.min(), fr.min())
# timing
t = np.linspace(-0.3, 0.3, 4000)
p = batman.TransitParams(); p.t0=0; p.per=5.69; p.rp=0.125; p.a=25.; p.inc=89.; p.ecc=0.; p.w=90.; p.limb_dark='quadratic'; p.u=[0.12,0.53]
m = batman.TransitModel(p, t, supersample_factor=5, exp_time=200/86400.)
t0=time.time()
for i in range(200): m.light_curve(p)
dtc = (time.time()-t0)/200
print('time per call (4000 pts x 5 supersample):', dtc)
json.dump(dict(max_abs_diff=res, sec_per_call_4000x5=dtc), open('out/batman_check.json','w'), indent=1)
