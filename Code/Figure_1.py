import matplotlib.pyplot as plt
import netCDF4 as nc
import cartopy.crs as ccrs
import numpy as np
import matplotlib.ticker as mticker
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.io.shapereader as shpreader
from matplotlib.colors import ListedColormap, BoundaryNorm
import pandas as pd
import geopandas as gpd
import matplotlib.colors as colors
import cmaps
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats
import matplotlib as mpl
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors

"""
================================================================================
Comprehensive Multi-Panel Figure Script
--------------------------------------------------------------------------------
(a) PM10 time series for selected stations and mean curve
(b) Observed Dust AOD (MERRA2)
(c) Spatial PM10 averages and exceedance days (Mar–May 2023)
(d) Dust emission [µg/m²/s] + 10 m winds
(e) 700 hPa geopotential height [gpm] + winds
================================================================================
"""

# === Global font and figure settings ===
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 1.3

fig = plt.figure(figsize=(16, 16), dpi=500)
gs = gridspec.GridSpec(21, 20, figure=fig)
ax = fig.add_subplot(gs[0:6, 0:20])

# === Time series setup ===
dates = pd.date_range(start="2023-03-01", end="2023-05-31", freq="D")

light_colors = [
    (0.85, 0.34, 0.34), (0.85, 0.62, 0.34), (0.80, 0.85, 0.34),
    (0.53, 0.85, 0.34), (0.34, 0.85, 0.43), (0.34, 0.85, 0.71),
    (0.34, 0.71, 0.85), (0.34, 0.43, 0.85), (0.53, 0.34, 0.85),
    (0.80, 0.34, 0.85), (0.85, 0.34, 0.62)
]
dark_color = 'darkblue'

all_data = pd.DataFrame()

# === Load daily PM10 station data ===
for date in dates:
    file_path = f'/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/airquality/china_sites_{date.strftime("%Y%m%d")}.csv'
    daily_data = pd.read_csv(file_path, delimiter=',', encoding='utf-8')
    hourly_data = daily_data.iloc[3::15].copy()
    time_index = pd.date_range(start=date, periods=len(hourly_data), freq='h')
    hourly_data.index = time_index
    all_data = pd.concat([all_data, hourly_data], axis=0)

# === Selected monitoring stations ===
selected_stations = ["1001A", "1015A", "1029A", "1488A", "1301A",
                     "1081A", "1317A", "1094A", "1476A", "1130A", "1119A"]
selected_labels = ["Beijing", "Tianjin", "Shijiazhuang", "Yinchuan",
                   "Jinan", "Taiyuan", "Zhengzhou", "Hohhot",
                   "Lanzhou", "Harbin", "Changchun"]

# === Plot PM10 time series ===
for station, label, light_color in zip(selected_stations, selected_labels, light_colors[:len(selected_stations)]):
    if station in all_data.columns:
        ax.plot(all_data.index, all_data[station], label=label, color=light_color, alpha=0.5)

# === Compute and plot average curve ===
average_data = all_data[selected_stations].mean(axis=1, skipna=True)
ax.plot(all_data.index, average_data, color=dark_color, linewidth=3, label='Average')
ax.axhline(y=1000, color='#E3738B', linestyle='--', linewidth=1.5)

# === Axis settings ===
ax.set_yscale('linear')
ax.set_ylim(0, 3500)
ax.set_yticks(np.arange(0, 2501, 200))
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: '{:.0f}'.format(x)))
ax.yaxis.set_major_locator(ticker.FixedLocator(np.arange(0, 3501, 500)))

ax.set_title('(a) PM10 Time Series [µg/m$^3$]', loc='left', fontsize=20, pad=8)
ax.legend(fontsize=17, ncol=4, loc='upper right', frameon=False)
ax.tick_params(axis='both', labelsize=21)
plt.xticks(rotation=0)

ax.minorticks_on()
ax.tick_params(axis="both", which="major", direction="out", width=1.3, length=7)
ax.tick_params(axis="both", which="minor", direction="out", width=1.3, length=3.5)
ax.tick_params(top=False, right=False)

# ==========================================================================================
# (c) Observed PM10 spatial statistics
# ==========================================================================================
axe = plt.subplot(324, projection=ccrs.PlateCarree(), aspect="auto")
sta = pd.read_excel('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_DATA/air_quality/zhandian.xlsx',
                    sheet_name='站点列表-2022 02 13起')

levels = np.arange(0, 160.1, 10)
norm = BoundaryNorm(levels, cmaps.WhiteBlueGreenYellowRed.N, clip=True)

pm10_sum, pm10_count, pm10_days = {}, {}, {}

# === Loop through daily files to compute mean and exceedance days ===
for date in pd.date_range("2023-03-01", "2023-05-31"):
    file_path = f'/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/airquality/china_sites_{date.strftime("%Y%m%d")}.csv'
    try:
        data_daily = pd.read_csv(file_path, delimiter=',', encoding='utf-8')
        for name in data_daily.columns[3::1]:
            if not pd.isna(data_daily[name][0]) and not pd.isna(data_daily[name][3]):
                if name not in pm10_sum:
                    pm10_sum[name], pm10_count[name], pm10_days[name] = 0, 0, 0
                daily_values = data_daily[name].iloc[3::15]
                daily_avg = daily_values.mean()
                pm10_sum[name] += daily_avg
                pm10_count[name] += 1
                if daily_avg > 150:
                    pm10_days[name] += 1
    except FileNotFoundError:
        continue

pm10_avg = {name: pm10_sum[name] / pm10_count[name] if pm10_count[name] > 0 else 0 for name in pm10_sum}

# === Define size scaling for exceedance days ===
def get_marker_size(days):
    if days < 10: return 35
    elif 10 <= days < 15: return 55
    elif 15 <= days < 20: return 75
    elif 20 <= days < 25: return 95
    elif 25 <= days < 30: return 115
    else: return 140

# === Plot PM10 mean values and marker size ===
for name in pm10_avg:
    if name in sta['监测点编码'].values:
        lon = float(sta[sta['监测点编码'] == name]['lon'].values[0])
        lat = float(sta[sta['监测点编码'] == name]['lat'].values[0])
        pm10 = pm10_avg[name]
        days = pm10_days[name]
        color = cmaps.WhiteBlueGreenYellowRed(norm(pm10))
        size = get_marker_size(days)
        axe.scatter(lon, lat, s=size, color=color, transform=ccrs.PlateCarree())

# === Legend for exceedance days ===
legend_x, legend_y = 93, 50
legend_sizes = [35, 55, 75, 95, 115, 140]
legend_labels = ["<10", "10-15", "15-20", "20-25", "25-30", ">30"]

cols, spacing_x, spacing_y = 3, 8, 2
for idx, (size, label) in enumerate(zip(legend_sizes, legend_labels)):
    row, col = idx // cols, idx % cols
    x_pos = legend_x + col * spacing_x
    y_pos = legend_y - row * spacing_y
    axe.scatter(x_pos, y_pos, s=size, color='gray', edgecolors='black', transform=ccrs.PlateCarree())
    axe.text(x_pos + 1, y_pos, label, fontsize=14, verticalalignment='center', transform=ccrs.PlateCarree())

axe.text(legend_x - 1, legend_y + spacing_y, "PM10 > 150 µg/m$^3$ Days",
         fontsize=16, fontweight='bold', transform=ccrs.PlateCarree())

axe.set_title('(c) Observed PM10 [µg/m$^3$]', loc='left', fontsize=20, pad=8)
axe.add_feature(cfeat.COASTLINE.with_scale('10m'), linewidth=0, color='k')
axe.set_extent([75, 135, 30, 55], crs=ccrs.PlateCarree())

gl = axe.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0, color='gray', linestyle=':')
gl.top_labels = gl.bottom_labels = gl.right_labels = gl.left_labels = False
gl.xlocator = mticker.FixedLocator(np.arange(80, 135, 15))
gl.ylocator = mticker.FixedLocator(np.arange(30, 56, 10))

axe.set_xticks(np.arange(80, 135, 15), crs=ccrs.PlateCarree())
axe.set_yticks(np.arange(30, 56, 10), crs=ccrs.PlateCarree())
axe.xaxis.set_major_formatter(LongitudeFormatter())
axe.yaxis.set_major_formatter(LatitudeFormatter())

plt.xticks(fontsize=21)
plt.yticks(fontsize=21)

# === Add shapefiles ===
shp_china = shpreader.Reader('/home/qiupch2023/data/shp/china/china.shp').geometries()
axe.add_geometries(shp_china, ccrs.PlateCarree(), facecolor='none', edgecolor='k', linewidth=0.5, zorder=1)
gdf = gpd.read_file('/home/qiupch2023/data/shp/china/china.shp')
axe.add_geometries(gdf.geometry, ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=0.2)

axe.tick_params(axis="both", which="major", direction="out", width=1.3, length=7)
axe.tick_params(axis="both", which="minor", direction="out", width=1.3, length=3.5)
axe.xaxis.set_minor_locator(mticker.MultipleLocator(5))
axe.yaxis.set_minor_locator(mticker.MultipleLocator(5))
axe.tick_params(top=False, right=False)

norm = BoundaryNorm(levels, cmaps.WhiteBlueGreenYellowRed.N, clip=True)
cb3 = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmaps.WhiteBlueGreenYellowRed),
                   ax=axe, orientation='horizontal', pad=0.1, shrink=1, aspect=22, extend='max')
cb3.outline.set_edgecolor('black')
cb3.ax.tick_params(labelsize=17, width=0, length=0, direction='out', colors='black')
ticks = np.arange(0, 161, 20)
cb3.set_ticks(ticks)
cb3.set_ticklabels([x for x in ticks])

# (The rest of the figure: Dust, Dust emission + wind, 700hPa GH + wind panels)
# ---------------------------------------------------------------------------------

# === Final layout and save ===
plt.subplots_adjust(left=0.06, bottom=0.03, right=0.95, top=0.95, wspace=0.2, hspace=0.1)
plt.savefig('figure1_Arial.png', dpi=500)
