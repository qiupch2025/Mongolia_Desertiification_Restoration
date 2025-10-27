import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import statsmodels.api as sm
import matplotlib.ticker as mticker
import matplotlib as mpl

# -------------------------
# Read WRF .nc files
# -------------------------
wrf_2023_file = '/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/figure_wrf/dust_shao04_2023.nc'
ds_2023 = xr.open_dataset(wrf_2023_file)
dust_2023 = ds_2023['dust'].values

wrf_2005_file = '/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/figure_wrf/dust_shao04_2005.nc'
ds_2005 = xr.open_dataset(wrf_2005_file)
dust_2005 = ds_2005['dust'].values

wrf_ideal_file = '/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/figure_wrf/dust_shao04_ideal_new.nc'
ds_ideal = xr.open_dataset(wrf_ideal_file)
dust_ideal = ds_ideal['dust'].values

lon_wrf = ds_ideal['lon'].values
lat_wrf = ds_ideal['lat'].values

# Time slice (assuming [time, lat, lon])
dust_layer_2023 = np.mean(dust_2023[6*7-2:6*99-2, :, :], axis=0)
dust_layer_2005 = np.mean(dust_2005[6*7-2:6*99-2, :, :], axis=0)
dust_layer_ideal = np.mean(dust_ideal[6*7-2:6*99-2, :, :], axis=0)

# Flatten WRF grid
points_wrf = np.column_stack((lon_wrf.ravel(), lat_wrf.ravel()))

# -------------------------
# Read stations & observations
# -------------------------
station_data = pd.read_excel('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_DATA/air_quality/zhandian.xlsx')
station_data['lon'] = pd.to_numeric(station_data['lon'], errors='coerce')
station_data['lat'] = pd.to_numeric(station_data['lat'], errors='coerce')

lon_station = station_data['lon'].values
lat_station = station_data['lat'].values
points_station = np.column_stack((lon_station, lat_station))

# Explicitly read PM10 observations (no daily averaging first)
all_data = pd.DataFrame()
dates = pd.date_range(start="2023-03-01", end="2023-05-31", freq="D")

for date in dates:
    file_path = f'/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/airquality/china_sites_{date.strftime("%Y%m%d")}.csv'
    daily_data = pd.read_csv(file_path, delimiter=',', encoding='utf-8')
    pm10_data = daily_data[daily_data['type'] == 'PM10']
    pm10_values = pm10_data.iloc[:, 3:].apply(pd.to_numeric, errors='coerce')
    all_data = pd.concat([all_data, pm10_values], axis=0)

pm10_obs = all_data.mean(axis=0, skipna=True).values

# Keep valid station entries only
valid_idx = (~np.isnan(lon_station)) & (~np.isnan(lat_station)) & (~np.isnan(pm10_obs))
lon_station = lon_station[valid_idx]
lat_station = lat_station[valid_idx]
pm10_obs = pm10_obs[valid_idx]
points_station = points_station[valid_idx, :]

# -------------------------
# Interpolate model dust to station points (convert to mg/m^2)
# -------------------------
dust_interp_2023 = griddata(points_wrf, dust_layer_2023.ravel(), points_station, method='linear') / 1000
dust_interp_2005 = griddata(points_wrf, dust_layer_2005.ravel(), points_station, method='linear') / 1000
dust_interp_ideal = griddata(points_wrf, dust_layer_ideal.ravel(), points_station, method='linear') / 1000

# Save to NetCDF
result_df = pd.DataFrame({
    'lon': lon_station,
    'lat': lat_station,
    'pm10_obs': pm10_obs,
    'pm10_2023': dust_interp_2023,
    'pm10_2005': dust_interp_2005,
    'pm10_ideal': dust_interp_ideal
})
result_ds = xr.Dataset.from_dataframe(result_df)
result_ds.to_netcdf('pm10_dust.nc')

# -------------------------
# Plot: correlation between observed PM10 and simulated dust (2023)
# -------------------------
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 1.3

fig = plt.figure(figsize=(10, 7))
axe = plt.subplot(111)

# Open the saved NetCDF
file_path = './pm10_dust.nc'
nc_file = Dataset(file_path, 'r')
x = nc_file.variables['pm10_obs'][:]
y = nc_file.variables['pm10_2023'][:]
nc_file.close()

# Clean NaNs consistently for regression & plotting
mask = np.isfinite(x) & np.isfinite(y)
x = x[mask]
y = y[mask]

# Scatter
plt.scatter(x, y, color=[64/255, 127/255, 181/255], label='Data Points', s=100, marker='o', alpha=0.5)

# OLS fit (with intercept)
x_fit = np.linspace(x.min(), x.max(), 200)
X = sm.add_constant(x)                # add intercept
model = sm.OLS(y, X, missing='drop')  # safety: drop remaining missing
results = model.fit()
y_fit = results.predict(sm.add_constant(x_fit))

# Fit line
plt.plot(x_fit, y_fit, color=[167/255, 31/255, 45/255], label='Fit Line', linewidth=3, alpha=1)

# Axes limits
plt.xlim([0, 400])
plt.ylim([0, 400])

# Labels & title
plt.title('(b) Correlation between PM10 and Simulated Dust', loc='left', fontsize=20, pad=8)
plt.xlabel('Observed PM10 [µg/m$^3$]', fontsize=20)
plt.ylabel('Simulated Dust [mg/m$^2$]', fontsize=20)

# Stats: R², p-value for slope, Pearson r (with sign)
r_squared = results.rsquared
p_value = results.pvalues[1] if results.pvalues.shape[0] > 1 else np.nan
if p_value < 0.01:
    p_value_text = 'p < 0.01'
elif p_value < 0.05:
    p_value_text = 'p < 0.05'
else:
    p_value_text = 'p ≥ 0.05'

# Pearson r
r = np.corrcoef(x, y)[0, 1]

# Text box (upper-left inside axes)
plt.text(0.55, 0.25, f'r = {r:.2f}\n{p_value_text}',
         fontsize=30, transform=plt.gca().transAxes, verticalalignment='top')

# Ticks
plt.xticks(fontsize=19)
plt.yticks(fontsize=19)
plt.tick_params(top='on', right='on', which='both')

# Major/minor ticks
axe.minorticks_on()
axe.tick_params(axis="both", which="major", direction="out", width=1.3, length=7)
axe.tick_params(axis="both", which="minor", direction="out", width=1.3, length=3.5)
axe.xaxis.set_minor_locator(mticker.MultipleLocator(50))
axe.yaxis.set_minor_locator(mticker.MultipleLocator(50))  # <- sensible for 0–400 range

# Hide top/right spines' ticks only (keep spines themselves)
axe.tick_params(top=False, right=False)
axe.tick_params(axis="x", which="both", top=False)
axe.tick_params(axis="y", which="both", right=False)

plt.subplots_adjust(left=0.12, bottom=0.12, right=0.95, top=0.93, wspace=0.2, hspace=0.1)
plt.savefig('Figure_S5b.png', dpi=500)
