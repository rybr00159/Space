import sys, time
from triceratops.funcs import query_TRILEGAL, save_trilegal
url=query_TRILEGAL(330.887948,-42.835511,verbose=1)
print('url',url,flush=True)
for k in range(40):
    try:
        fn=save_trilegal(url,71316387); print('saved',fn,flush=True); break
    except Exception as e:
        print('wait',str(e)[:100],flush=True); time.sleep(30)
