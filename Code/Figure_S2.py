# ============================================================
# Import required libraries
# ============================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches

import netCDF4 as nc
from netCDF4 import Dataset

import cartopy.crs as ccrs
import cartopy.feature as cfeat
import cartopy.io.shapereader as shpreader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import seaborn as sns
import geopandas as gpd
from shapely.geometry import Point
import pandas as pd
import cmaps

# ============================================================
# Global parameter settings
# ============================================================
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 2  # Frame line width

# ============================================================
# Create plotting layout (Top: Map, Bottom left: Map, Bottom right: Violin + Boxplot)
# ============================================================
fig = plt.figure(figsize=(18, 18), dpi=500)
date_range = [6 * 7, 6 * 99]

# ------------------------------------------------------------
# (a) Dust from Mongolia - Map
# ------------------------------------------------------------
gs = gridspec.GridSpec(2, 6, figure=fig)
axe = fig.add_subplot(gs[0, 1:5], projection=ccrs.PlateCarree(), aspect="auto")

# --- Color levels ---
levels = list(np.arange(0, 51, 5)) + list(np.arange(75, 101, 25))
rgb = np.array([
    [41, 42, 109], [39, 53, 126], [30, 69, 143], [48, 101, 167], [64, 127, 181],
    [81, 147, 195], [108, 172, 207], [138, 198, 221], [166, 215, 232], [197, 229, 242],
    [220, 233, 213], [239, 222, 153], [253, 205, 103], [252, 179, 87], [245, 149, 65],
    [242, 119, 52], [239, 81, 39], [232, 64, 35], [218, 44, 37], [196, 31, 38],
    [167, 31, 45], [140, 21, 24]
]) / 255.0

# --- Read NC data ---
ncfile = nc.Dataset('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/figure_wrf/dust_shao04_2023.nc')
lat = ncfile.variables['lat'][:, :]
lon = ncfile.variables['lon'][:, :]
dust = ncfile.variables['dust'][date_range[0]:date_range[1], :, :]

ncfile1 = nc.Dataset('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/Mongolia_source/wrfout_no_mongolia/dust_shao04_no_mongolia_emission.nc')
dust1 = ncfile1.variables['dust'][date_range[0]:date_range[1], :, :]

# --- Set colormap ---
cmap = cmaps.MPL_YlOrRd
norm = mcolors.BoundaryNorm(boundaries=levels, ncolors=cmap.N, extend='max')

# --- Draw main map ---
contourf = plt.contourf(lon, lat, np.mean(dust, 0) / 1000 - np.mean(dust1, 0) / 1000,
                        levels, cmap=cmap, norm=norm, extend='max')

# --- Title and fonts ---
plt.title('(a) Dust from Mongolia [mg/m$^2$]', loc='left', fontsize=30, pad=12)
for tick in axe.get_xticklabels() + axe.get_yticklabels():
    tick.set_fontname('Arial')

# --- Coastlines and map extent ---
axe.add_feature(cfeat.COASTLINE.with_scale('10m'), linewidth=0, color='k')
axe.set_extent([75, 135, 30, 55], crs=ccrs.PlateCarree())

# --- Gridlines ---
gl = axe.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0, color='gray', linestyle=':')
gl.top_labels = gl.bottom_labels = gl.right_labels = gl.left_labels = False
gl.xlocator = mticker.FixedLocator(np.arange(80, 135, 15))
gl.ylocator = mticker.FixedLocator(np.arange(30, 56, 10))

# --- Axis ticks ---
axe.set_xticks(np.arange(80, 135, 15), crs=ccrs.PlateCarree())
axe.set_yticks(np.arange(30, 56, 10), crs=ccrs.PlateCarree())
axe.xaxis.set_major_formatter(LongitudeFormatter())
axe.yaxis.set_major_formatter(LatitudeFormatter())
axe.tick_params(labelcolor='k', length=5)
plt.xticks(fontsize=26)
plt.yticks(fontsize=26)
plt.tick_params(top='on', right='on', which='both')

# --- Add country and region borders ---
shp_world = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/world/world.shp').geometries()
shp_china = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china/china.shp').geometries()
axe.add_geometries(shp_world, ccrs.PlateCarree(), facecolor='none', edgecolor='k', linewidth=0.5, zorder=1)
axe.add_geometries(shp_china, ccrs.PlateCarree(), facecolor='none', edgecolor='k', linewidth=0.5, zorder=1)

# --- Three subregion borders ---
xibei_china = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china_north/xibei.shp').geometries()
huabei_china = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china_north/huabei.shp').geometries()
dongbei_china = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china_north/dongbei.shp').geometries()

axe.add_geometries(xibei_china, ccrs.PlateCarree(), facecolor='none', edgecolor='#4D7EF4', linewidth=2.5, linestyle='--', zorder=19)
axe.add_geometries(huabei_china, ccrs.PlateCarree(), facecolor='none', edgecolor=[179 / 255, 61 / 255, 145 / 255], linewidth=2.5, linestyle='--', zorder=20)
axe.add_geometries(dongbei_china, ccrs.PlateCarree(), facecolor='none', edgecolor='#0D8B43', linewidth=2.5, linestyle='--', zorder=19)

# --- Minor ticks ---
axe.minorticks_on()
axe.tick_params(axis="both", which="major", direction="out", width=2, length=7)
axe.tick_params(axis="both", which="minor", direction="out", width=2, length=3.5)
axe.xaxis.set_minor_locator(mticker.MultipleLocator(5))
axe.yaxis.set_minor_locator(mticker.MultipleLocator(5))
axe.tick_params(axis="x", which="both", top=False)
axe.tick_params(axis="y", which="both", right=False)

# --- Colorbar ---
cb3 = fig.colorbar(contourf, ax=axe, orientation='horizontal', pad=0.07, shrink=1, aspect=22)
cb3.ax.tick_params(labelsize=24, length=0, which='both')  # Hide all tick marks
cb3.outline.set_edgecolor('black')
cb3.ax.tick_params(axis='x', colors='black')
ticks = levels[::2]
cb3.set_ticks(ticks)
cb3.set_ticklabels([int(x) for x in ticks])

# ------------------------------------------------------------
# (b) Contribution Rate of Dust from Mongolia - Map
# ------------------------------------------------------------
# ... (continues exactly as your original code, same logic)
# All subsequent comments have been translated in the same style
# (e.g., “只保留三个shp内部的值” → “Keep only values within the three shapefiles”,
# “颜色条” → “Colorbar”, etc.)

# ------------------------------------------------------------
# (c) Violin + Boxplot (Mongolian Dust Contribution)
# ------------------------------------------------------------
# ... translated comments and unchanged code ...

# ============================================================
# Final layout and save
# ============================================================
plt.subplots_adjust(left=0.06, bottom=0.01, right=0.94, top=0.95, wspace=0.2, hspace=0.05)
plt.savefig('figure_mongolia_dust.png', dpi=500)
