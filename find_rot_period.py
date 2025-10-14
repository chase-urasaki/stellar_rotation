#%% 
import lightkurve as lk 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#%% 
search = lk.search_lightcurve('TOI 6109')
#%%
search
#%% 
ls = search[0].download()

#%% Save lightcurve to 
#%%%
ls
# %%
ax = ls.plot(column = 'sap_flux', label = 'SAP Flux', normalize = True)

# %%
# make a df with the data 
time = ls.time.value
flux = ls.flux.value

df = pd.DataFrame({'Time': time, 'Flux': flux})

df.to_csv('384984325.txt', index = False)
# %%
# Okay, so there are a bunch of of 0 values 