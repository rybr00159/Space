import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
H = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/lcmodel/results/'
def prof(res, keyfix='L3fixed', small_k_only=True):
    best = {}
    for r in res:
        k = r[keyfix]
        if k is None: continue
        tot = r.get('tot', r.get('chi2', 0) + r.get('prior2', 0))
        if small_k_only and r['k'] > 1.0: continue
        if k not in best or tot < best[k]: best[k] = tot
    return best
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
R = json.load(open(H + 'ml_profile_A.json')) + json.load(open(H + 'ml_profile_A_hiL3.json'))
for r in R: r['tot'] = r['chi2'] + r['prior2']
b = prof(R); x = np.array(sorted(b)); y = np.array([b[v] for v in x]); y -= y.min()
ax[0].plot(x, y, 'ko-', label='real data (small-companion branch)')
for c, col in [('A', 'C3'), ('B', 'C0'), ('C', 'C2')]:
    J = json.load(open(H + f'inj_recover_{c}.json'))
    bb = prof(J['all']); xx = np.array(sorted(bb)); yy = np.array([bb[v] for v in xx]); yy -= min(yy.min(), min(r['tot'] for r in J['all'] if r['L3fixed'] is None))
    ax[0].plot(xx, yy, 'o--', color=col, ms=4, label=f"injection {c}: true L3={J['true']['L3']}, k={J['true']['k']}, bn={J['true']['bn']}")
    ax[0].axvline(J['true']['L3'], color=col, ls=':', lw=0.8)
ax[0].set_xlabel('third light L3 (fixed)'); ax[0].set_ylabel('profile delta chi2'); ax[0].set_ylim(-1, 30); ax[0].legend(fontsize=7)
ax[0].set_title('Profile likelihood in L3 (geometric model)')
for fn, lab in [('etv_companion_prelim_M31.03_MEB0.64.json', '')]:
    D = json.load(open(H + fn))
    for mode, ls in [('circular', '-'), ('thermal', '--')]:
        g = [r for r in D['grid'] if r['ecc'] == mode]
        ax[1].semilogx([r['a_au'] for r in g], [r['frac_excluded'] for r in g], ls, label=f'{mode} outer orbit')
ax[1].axvspan(297, 535, color='0.85', label='1.0-1.8" at 297 pc')
ax[1].set_xlabel('outer semi-major axis (au), M3 = 1.03 Msun, M_EB = 0.64 Msun'); ax[1].set_ylabel('fraction of orientations excluded (99.7%)')
ax[1].set_title('ETV (light-travel-time) exclusion of a tertiary'); ax[1].legend(fontsize=8); ax[1].set_ylim(0, 1.02)
plt.tight_layout(); plt.savefig(H + '../plots/profiles_companion.png', dpi=110)
