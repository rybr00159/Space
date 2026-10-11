"""Novelty checks: is a (TIC, period) signal already known?"""
import os, re, time, functools
import numpy as np, pandas as pd


def retry(fn):
    """Retry network lookups with backoff (VizieR/SIMBAD occasionally reset connections)."""
    @functools.wraps(fn)
    def wrap(*a, **k):
        for i in range(5):
            try:
                return fn(*a, **k)
            except Exception:
                if i == 4: raise
                time.sleep(3 * 2 ** i)
    return wrap

K = os.path.join(os.path.dirname(__file__), "..", "known")
_cache = {}


def _tables():
    if _cache: return _cache
    toi = pd.read_csv(f"{K}/toi.csv")
    _cache["toi"] = toi.rename(columns={"TIC ID": "tic", "Period (days)": "P", "TOI": "name", "TFOPWG Disposition": "disp"})[["tic", "P", "name", "disp"]]
    ctoi = pd.read_csv(f"{K}/ctoi.csv")
    _cache["ctoi"] = ctoi.rename(columns={"TIC ID": "tic", "Period (days)": "P", "CTOI": "name", "User Disposition": "disp"})[["tic", "P", "name", "disp"]]
    ps = pd.read_csv(f"{K}/ps.csv")
    ps["tic"] = pd.to_numeric(ps.tic_id.astype(str).str.replace("TIC ", ""), errors="coerce")
    _cache["ps"] = ps.rename(columns={"pl_orbper": "P", "pl_name": "name"})[["tic", "P", "name"]].assign(disp="confirmed")
    tce = pd.read_parquet(f"{K}/tce_all.parquet")
    _cache["tce"] = tce.rename(columns={"ticid": "tic", "tce_period": "P", "file": "name"})[["tic", "P", "name"]].assign(disp="SPOC TCE")
    return _cache


def period_match(P1, P2, tol=0.01, kmax=3):
    """True when P1/P2 is within tol of k or 1/k for k <= kmax."""
    if not (np.isfinite(P1) and np.isfinite(P2)) or P2 <= 0: return False
    r = P1 / P2
    return any(abs(r - k) < tol * k or abs(r - 1 / k) < tol / k for k in range(1, kmax + 1))


def known_signals(tic, P):
    """Return (any_entry_on_star, matching_entries) across TOI, CTOI, confirmed planets and SPOC TCEs."""
    t = _tables(); on_star, match = [], []
    for k in ("toi", "ctoi", "ps", "tce"):
        rows = t[k][t[k].tic == tic]
        for _, r in rows.iterrows():
            e = f"{k}:{r['name']}:{r.P:.4f}" if np.isfinite(r.P) else f"{k}:{r['name']}:nan"
            on_star.append(e)
            if period_match(P, r.P): match.append(e)
    return on_star, match


VIZ_CATS = {"B/vsx/vsx": "VSX", "II/366/catv2021": "ASAS-SN var", "J/ApJS/258/16": "TESS EB (Prsa+22)",
            "I/358/veb": "Gaia DR3 EB", "I/358/vclassre": "Gaia DR3 var class",
            "J/ApJS/279/50/table3": "TESS 10k EB (Kostov+25) vetted new", "J/ApJS/279/50/table4": "TESS 10k EB (Kostov+25) known",
            "J/MNRAS/517/2190": "ASAS-SN EB (Rowan+22)", "J/AJ/156/241": "ATLAS var (Heinze+18)",
            "J/MNRAS/502/1299": "SuperWASP var (Thiemann+21)", "J/ApJS/276/57": "TESS periodic var (Gao+25)",
            "J/A+A/691/A242": "TESS OBAF EB (IJspeert+24)", "J/MNRAS/503/3975": "WISE var (Petrosky+21)"}


@retry
def kostov_unvetted(ra, dec, radius_arcsec=21.0):
    """Membership in the ~873k unvetted neural-network 'eclipse-like' list of Kostov+2025 (table 2, no periods)."""
    from astroquery.vizier import Vizier
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    r = Vizier(columns=["TIC", "_r"], row_limit=5).query_region(SkyCoord(ra * u.deg, dec * u.deg), radius=radius_arcsec * u.arcsec,
                                                                 catalog="J/ApJS/279/50/table2")
    return [(int(x["TIC"]), float(x["_r"])) for t in r for x in t]


@retry
def vizier_cone(ra, dec, radius_arcsec=42.0):
    """Variable-star / eclipsing-binary catalogue entries within radius (2 TESS pixels by default)."""
    from astroquery.vizier import Vizier
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    v = Vizier(columns=["**", "_r"], row_limit=50)
    res = v.query_region(SkyCoord(ra * u.deg, dec * u.deg), radius=radius_arcsec * u.arcsec, catalog=list(VIZ_CATS))
    out = []
    for tab in res:
        name = VIZ_CATS.get(tab.meta.get("name", ""), None) or next((n for c, n in VIZ_CATS.items() if tab.meta.get("name", "").startswith(c)), tab.meta.get("name"))
        for r in tab:
            d = {c: (r[c].item() if hasattr(r[c], "item") else r[c]) for c in tab.colnames}
            per = next((d[c] for c in ("Period", "Per", "period", "P") if c in d and d[c] not in (None, "") and np.isfinite(float(d[c] or np.nan))), None)
            if per is None and "Freq" in d and d["Freq"] not in (None, "") and np.isfinite(float(d["Freq"] or np.nan)) and float(d["Freq"]) > 0:
                per = 1.0 / float(d["Freq"])
            typ = next((d[c] for c in ("Type", "VarType", "Class", "Cl", "best_class_name") if c in d), None)
            out.append(dict(cat=name, sep=float(d.get("_r", np.nan)), period=per, type=typ,
                            id=next((d[c] for c in ("Name", "ID", "Source", "TIC", "OID") if c in d), None)))
    return out


KEYWORDS = ("variab", "eclips", "binar", "EB ", "EBs", "transit", "TOI", "planet", "periodic", "light curve", "lightcurve", "pulsat", "rotation")


@retry
def vizier_all(ra, dec, radius_arcsec=15.0):
    """Every VizieR catalogue with a source within radius, keeping those whose title suggests variability work."""
    from astroquery.vizier import Vizier
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    v = Vizier(columns=["*", "_r"], row_limit=5, timeout=300)
    res = v.query_region(SkyCoord(ra * u.deg, dec * u.deg), radius=radius_arcsec * u.arcsec)
    hits = []
    for name in res.keys():
        desc = res[name].meta.get("description", "") or ""
        if any(k.lower() in desc.lower() for k in KEYWORDS):
            hits.append(f"{name}: {desc[:90]}")
    return hits, len(res)


@retry
def simbad(ra, dec, radius_arcsec=30.0):
    """SIMBAD objects within radius with their main identifier and object type."""
    import requests, io
    q = f"""SELECT main_id, otype, DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {ra}, {dec}))*3600 AS sep
            FROM basic WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {radius_arcsec / 3600}))=1 ORDER BY sep"""
    r = requests.post("https://simbad.cds.unistra.fr/simbad/sim-tap/sync",
                      data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q}, timeout=120)
    r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text)).to_dict("records")
