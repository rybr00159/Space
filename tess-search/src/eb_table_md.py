"""Markdown table of the final EB list for the README."""
import csv

W = "/home/user/work"


def f(x, n):
    try: return f"{float(x):.{n}f}"
    except (TypeError, ValueError): return ""


rows = list(csv.DictReader(open(f"{W}/vet/new_ebs_final.csv")))
lines = ["| TIC | Other name | T mag | Period (d) | Primary | Secondary (phase) | Sectors | Type | Notes |",
         "|---|---|---|---|---|---|---|---|---|"]
for r in sorted(rows, key=lambda r: (r["kind"] != "detached", -float(r["primary_pct"]))):
    try: perr = float(r["P_err"])
    except ValueError: perr = float("nan")
    pdig = 6 if perr < 5e-6 else (5 if perr < 5e-5 else 4)
    lines.append(f"| {r['tic']} | {r['name'] or '-'} | {f(r['Tmag'], 2)} | {float(r['P']):.{pdig}f} | {f(r['primary_pct'], 1)} % "
                 f"| {f(r['secondary_pct'], 2)} % ({f(r['phase_sec'], 3)}) | {r['n_sectors']} | {r['kind']} | {r['note']} |")
open(f"{W}/vet/new_ebs_table.md", "w").write("\n".join(lines) + "\n")
print("table rows", len(lines) - 2)
