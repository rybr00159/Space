"""One-line verdict inputs for a list of vet_report.py outputs: size, geometry, density, odd/even, secondary,
held-out depth spread, difference-image offsets and Gaia RUWE."""
import sys, json
import numpy as np

W = "/home/user/work"


def summary(tic):
    try: r = json.load(open(f"{W}/vet/TIC{tic}/report.json"))
    except Exception as e: return f"TIC {tic}: no report ({e})"
    s, f, c = r["star"], r["fit"], r["combined"]
    ho = [h["depth"] for h in r.get("held_out", []) if h.get("depth") and h.get("snr", 0) > 5]
    spread = (max(ho) / min(ho)) if len(ho) > 1 and min(ho) > 0 else np.nan
    loc = [(L["sector"], round(L["offset_arcsec"], 1), round(L["offset_sigma"], 1)) for L in r.get("localisation", [])
           if L.get("offset_sigma") == L.get("offset_sigma")]
    g = r.get("gaia") or {}
    ruwe = g.get("ruwe") if isinstance(g, dict) else None
    return (f"TIC {tic}: T={s['Tmag']:.2f} Teff={s['Teff']} R={s['rad']} | P={f['P']:.5f} depth={c['depth']*1e6:.0f}ppm "
            f"Rp={r['rp_earth']:.1f}Re b={f['b']:.2f} rho_fit/tic={r['rho_fit']:.2f}/{r['rho_tic']} | odd/even={c['oddeven_sig']:.1f}s "
            f"sec={c['sec_depth']*1e6:.0f}ppm({c['sec_sig']:.1f}s) oot_rms={c['oot_rms_ratio']:.2f} | held-out max/min={spread:.2f} n={len(ho)} "
            f"| loc={loc} | ruwe={ruwe} contratio={s.get('contratio')}")


if __name__ == "__main__":
    for t in sys.argv[1:]: print(summary(int(t)))
