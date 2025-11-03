#%% 
import lightkurve as lk 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#%% 
TARGET = 'TOI 6109'
#%%
search = lk.search_lightcurve(TARGET)
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

df.to_csv('lc/raw/105840719.txt', index = False)
# %%
# Okay, so there are a bunch of of 0 values 
# %%
# Try gyro interp 
from gyrointerp import gyro_age_posterior

Prot = 5.763533783432416
Prot_err = 0.61827224641155

Teff, Teff_err = 5432, 125

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
    'ylabel': r'Probablility ($10^{-3}$Myr$^{-1}$)',
    'title': f'P_rot = {Prot} days, Teff = {Teff} K',
    'xlim': [0, 1000]
})
# Plt the maximum value 
plt.axvline(result['median'], ls = '--', c = 'C1', label = 'Median Age')
plt.text(result['median'] + 20, ax.get_ylim()[1]*0.8, f"{result['median']:.0f} Myr", color = 'C1')\

# Plot the error bars on this 
plt.fill_betweenx(
    y = ax.get_ylim(),
    x1 = result['median'] - result['-1sigma'],
    x2 = result['median'] + result['+1sigma'],
    color = 'C1',
    alpha = 0.3,
    label = '1$\sigma$ Interval'
)
plt.text()
# Label the target star as a text box off to the side
plt.text(800, ax.get_ylim()[1]*0.9, TARGET, fontsize = 12, bbox = dict(facecolor = 'white', edgecolor = 'black'))
plt.show()
# %%
