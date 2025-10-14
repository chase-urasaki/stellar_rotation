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
# %%
# Try gyro interp 
from gyrointerp import gyro_age_posterior

Prot = 3.1944784409461136
Prot_err = 0.5555614679906284

Teff, Teff_err = 5660, 100

# Unifrom grid between 0 and 1000 Myr
age_grid = np.linspace(0, 1000, 500)

# Calculate the age posterior at each step 
age_poster = gyro_age_posterior(
    Prot, Teff,
    Prot_err = Prot_err,
    Teff_err = Teff_err,
    age_grid = age_grid
)
# %%
from gyrointerp import get_summary_statistics

result = get_summary_statistics(age_grid, age_poster)

print(f"age = {result['median']} + {result['+1sigma']} -{result['-1sigma']} Myr.")
# %%
fig, ax = plt.subplots()
ax.plot(age_grid, 1e3*age_poster, c = 'k', lw = 0.8)
ax.update({
    'xlabel': 'Age [Myr]',
    'ylabel': r'Probablility ($10^{-3}\$Myr$^{-1}$)',
    'title': f'P_rot = {Prot} days, Teff = {Teff} K',
    'xlim': [0, 1000]
})
plt.show()
# %%
