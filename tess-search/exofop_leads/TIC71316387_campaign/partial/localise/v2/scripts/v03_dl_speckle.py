from astroquery.vizier import Vizier
import warnings; warnings.filterwarnings('ignore')
O = '/tmp/claude-0/-home-user-Space/121b4aaf-900a-5ff8-a67e-bfbb30b2a844/scratchpad/tic713/localise/v2/speckle'
tabs = ['J/AJ/159/19/table6', 'J/AJ/162/192/table1', 'J/AJ/162/192/table2', 'J/AJ/165/180/table1',
        'J/AJ/168/28/table1', 'J/AJ/168/28/table2', 'J/AJ/164/58/table2', 'J/AJ/164/58/table3',
        'J/AJ/162/41/table2', 'J/AJ/162/41/table3', 'J/AJ/160/7/table2', 'J/AJ/160/7/table3',
        'J/AJ/158/48/table2', 'J/AJ/158/48/table3', 'J/AJ/155/235/table2', 'J/AJ/155/235/table3',
        'J/AJ/166/139/table2', 'J/AJ/166/139/table3']
v = Vizier(row_limit=-1, columns=['**'], timeout=120)
for t in tabs:
    try:
        r = v.get_catalogs(t)
        if len(r) == 0:
            print(t, 'EMPTY'); continue
        tab = r[0]
        fn = O + '/' + t.replace('/', '_') + '.ecsv'
        tab.write(fn, overwrite=True, format='ascii.ecsv')
        print(t, len(tab), tab.meta.get('description', '')[:70], '|', ','.join(tab.colnames)[:200])
    except Exception as e:
        print(t, 'ERR', e)
