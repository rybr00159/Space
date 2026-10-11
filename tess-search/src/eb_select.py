"""Curate the final list of newly characterised eclipsing binaries from all eb_analysis outputs.

Quality cuts (all required):
  * >= 2 TESS sectors; difference-image offset <= 2.5 sigma in every measured sector (on target)
  * primary eclipse >= 0.5 %; a coherent ephemeris (timing rms < 15 min, or primary >= 1 % with secondary >= 10 sigma)
  * no period-matched or within-1-pixel entry in any EB/variable catalogue, no variability catalogue hit in the all-VizieR
    sweep, SIMBAD type not a known variable/binary
Detached vs contact/ellipsoidal is decided by the out-of-eclipse scatter. When two candidates share a period within 3'
(light leaking from one EB into a neighbour's aperture) only the better-localised one is kept.
"""
import json, glob
import numpy as np, pandas as pd

W = "/home/user/work"
SKIP = ("IV/38", "IV/39", "J/AJ/155/39", "I/355/epphot", "I/358/varisum")
VAR = ("EB*", "SB*", "V*", "Ro*", "RR*", "dS*", "BY*", "RS*", "Pu*", "LP*", "Mi*", "sr*", "Ce*", "WU*", "Al*", "El*", "**", "gD*")

# visual review of fold plots (vet/eb_contact_*.png): sinusoidal variables or inconsistent folds, not EBs
REJECT_VISUAL = {237329004, 166225041, 60829240, 422003560, 387291911}
# periods fixed by eclipse timing (equal eclipses folded on half/quarter of the orbit)
PERIOD_FIX = {100666720: 7.341134, 49617965: 3.765198}

known = set()
for f in glob.glob(f"{W}/results/new_s0*.csv"):
    n = pd.read_csv(f); known |= set(n[n.viz_known].tic.astype(int))

rows = []
for f in glob.glob(f"{W}/vet/EB_TIC*/eb.json"):
    r = json.load(open(f))
    if "P" not in r: continue
    s = r["star"]; loc = [L for L in r.get("localisation", []) if L.get("offset_sigma") == L.get("offset_sigma")]
    viz = [v for v in r.get("vizier_variability_15", []) if not v.startswith(SKIP)]
    sim = [x["otype"] for x in r.get("simbad_30", []) if float(x["sep"]) < 10]
    eph = r.get("timing_ephemeris") or {}
    alt = r.get("alt_depths") or [np.nan, np.nan]
    alt_ratio = min(alt) / max(alt) if all(a == a and a > 0 for a in alt) else np.nan
    rows.append(dict(alt_ratio=alt_ratio, dur_h=r.get("dur_h", np.nan), tic=r["tic"], ra=s["ra"], dec=s["dec"], Tmag=s["Tmag"], Teff=s["Teff"], P=r["P"], t0=r["t0"],
                     prim=r["depth_prim"], sec=r["depth_sec"], phase_sec=r["phase_sec"], sec_sig=r["sec_sig"], oot_ppm=r["oot_std"] * 1e6,
                     nsec=len(r["sectors"]), worst_loc=max([L["offset_sigma"] for L in loc]) if loc else np.nan, nloc=len(loc),
                     timing_rms=eph.get("rms_min", np.nan), viz=";".join(viz), simbad=";".join(sim), cat_known=r["tic"] in known))
d = pd.DataFrame(rows)
d["on_target"] = (d.nloc >= 1) & (d.worst_loc <= 2.5)
d["coherent"] = (d.timing_rms < 15) | ((d.prim >= 0.01) & (d.sec_sig >= 10))
d["novel"] = ~d.cat_known & (d.viz == "") & ~d.simbad.apply(lambda s: any(t in VAR for t in s.split(";") if t))
d["ok"] = (d.nsec >= 2) & d.on_target & (d.prim >= 0.005) & d.coherent & d.novel
# morphology from the fraction of the orbit spent in eclipse (noise in faint stars makes scatter a poor guide)
# near-equal eclipses folded on half the orbit: alternate eclipse depths differ -> true period is 2P
# (only when nothing comparable already sits at phase 0.5, i.e. the second eclipse family is on the same grid)
d["P_doubled"] = (d.alt_ratio < 0.88) & (d.prim >= 0.01) & (d.sec / d.prim < 0.5)
d.loc[d.P_doubled, "P"] = 2 * d.loc[d.P_doubled, "P"]
for t, p in PERIOD_FIX.items(): d.loc[d.tic == t, "P"] = p
d["ok"] &= ~d.tic.isin(REJECT_VISUAL)
d["duty"] = d.dur_h / (d.P * 24)
d["kind"] = np.where((d.duty < 0.12) & (d.P > 0.8), "detached", "contact/ellipsoidal")
# drop period twins within 3' (keep the better-localised)
ok = d[d.ok].sort_values("worst_loc")
keep = []
for r in ok.itertuples():
    twin = [k for k in keep if abs(k.P - r.P) / r.P < 1e-4 and np.hypot((k.ra - r.ra) * np.cos(np.radians(r.dec)), k.dec - r.dec) * 3600 < 180]
    if not twin: keep.append(r)
final = pd.DataFrame(keep).drop(columns=["Index"], errors="ignore").sort_values(["kind", "prim"], ascending=[True, False])
final.to_csv(f"{W}/vet/eb_selected.csv", index=False)
print(f"{len(d)} analysed; {int(d.ok.sum())} pass cuts; {len(final)} after removing contamination twins")
print(final.groupby("kind").size().to_dict())
pd.set_option("display.width", 220)
print(final[["tic", "Tmag", "Teff", "P", "P_doubled", "alt_ratio", "prim", "sec", "phase_sec", "sec_sig", "oot_ppm", "nsec", "worst_loc", "timing_rms", "kind"]].round(4).to_string())
