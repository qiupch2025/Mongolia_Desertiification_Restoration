import matplotlib.pyplot as plt
import netCDF4 as nc
import cartopy.crs as ccrs
import numpy as np
import matplotlib.ticker as mticker
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.io.shapereader as shpreader
import matplotlib as mpl
import matplotlib.patches as patches
from netCDF4 import Dataset
import matplotlib.colors as mcolors
import cmaps

"""
================================================================================
Multi-panel figure for WRF dust simulations and diagnostics

Panels:
(1) (a) Simulated Dust [mg/m^2] climatology over Mar–May subset (date_range)
(2) (b) Regional boxplots of dust differences (WRF_REAL - WRF_CTL vs WRF_REST - WRF_CTL)
(3) (c) Spatial dust difference: WRF_REAL - WRF_CTL [mg/m^2]
(4) (d) Spatial dust emission difference: WRF_REAL - WRF_CTL [µg/m^2/s]
(5) (e) Spatial dust difference: WRF_REST - WRF_CTL [mg/m^2]
(6) (f) Spatial dust emission difference: WRF_REST - WRF_CTL [µg/m^2/s]

================================================================================
"""

# Set global font to 'Arial' and axis frame width
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 1.3

# Create the figure
fig = plt.figure(figsize=(16, 16), dpi=500)

# Time index range (index-based, not datetime); used to subset simulation period
date_range = [6*7, 6*99]

#####################################################################################################################
# (a) Simulated Dust [mg/m^2]
axe = plt.subplot(321, projection=ccrs.PlateCarree(), aspect="auto")

# Color levels and colormap for dust (mg/m^2)
levels = list(np.arange(0, 100000, 5000)) + list(np.arange(100000, 400001, 15000))
cmap = cmaps.WhiteBlueGreenYellowRed
norm = mcolors.BoundaryNorm(boundaries=levels, ncolors=cmap.N, extend='max')

# Read WRF output: dust field and grid
ncfile = nc.Dataset('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/figure_wrf/dust_shao04_2023.nc')
lat = ncfile.variables['lat'][:, :]
lon = ncfile.variables['lon'][:, :]
dust = ncfile.variables['dust'][date_range[0]:date_range[1], :, :]

# Filled contours of mean dust over the selected period
contourf = plt.contourf(lon, lat, np.mean(dust, 0), levels, cmap=cmap, norm=norm, extend='max')

# Title and tick font
plt.title('(a) Simulated Dust [mg/m$^2$]', loc='left', fontsize=20, pad=8)
for tick in axe.get_xticklabels() + axe.get_yticklabels():
    tick.set_fontname('Arial')

# Coastlines and map extent
axe.add_feature(cfeat.COASTLINE.with_scale('10m'), linewidth=0, color='k')
axe.set_extent([75, 135, 30, 55], crs=ccrs.PlateCarree())

# Gridlines (labels disabled)
gl = axe.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0, color='gray', linestyle=':')
gl.top_labels, gl.bottom_labels, gl.right_labels, gl.left_labels = False, False, False, False
gl.xlocator = mticker.FixedLocator(np.arange(80, 135, 15))
gl.ylocator = mticker.FixedLocator(np.arange(30, 56, 10))

# Axis ticks and format
axe.set_xticks(np.arange(80, 135, 15), crs=ccrs.PlateCarree())
axe.set_yticks(np.arange(30, 56, 10), crs=ccrs.PlateCarree())
axe.xaxis.set_major_formatter(LongitudeFormatter())
axe.yaxis.set_major_formatter(LatitudeFormatter())
axe.tick_params(labelcolor='k', length=5)

# Tick label sizes and general tick params
plt.xticks(fontsize=21)
plt.yticks(fontsize=21)
plt.tick_params(top='on', right='on', which='both')

# World and China shapefiles
shp_world = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/world/world.shp').geometries()
shp_china = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china/china.shp').geometries()
axe.add_geometries(shp_world, ccrs.PlateCarree(), facecolor='none', edgecolor='k', linewidth=0.5, zorder=1)
axe.add_geometries(shp_china, ccrs.PlateCarree(), facecolor='none', edgecolor='k', linewidth=0.5, zorder=1)

# Minor ticks and aesthetics
axe.minorticks_on()
axe.tick_params(axis="both", which="major", direction="out", width=1.3, length=7)
axe.tick_params(axis="both", which="minor", direction="out", width=1.3, length=3.5)
axe.xaxis.set_minor_locator(mticker.MultipleLocator(5))
axe.yaxis.set_minor_locator(mticker.MultipleLocator(5))
axe.tick_params(top=False, right=False)
axe.tick_params(axis="x", which="both", top=False)
axe.tick_params(axis="y", which="both", right=False)

# Colorbar for panel (a)
cb = fig.colorbar(contourf, ax=axe, orientation='horizontal', pad=0.1, shrink=1, aspect=22)
ticks = levels[::4]
cb.set_ticks(ticks)
cb.set_ticklabels([int(x / 1000) for x in ticks])
cb.ax.tick_params(labelsize=17)
cb.ax.tick_params(which='both', size=0)

# 2##############################################################################################
# (b) Regional boxplots: WRF_REAL - WRF_CTL vs WRF_REST - WRF_CTL (dust differences)
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import geopandas as gpd
from shapely.geometry import Point

# Read regional shapefiles (Northwest, North, Northeast China)
region_paths = {
    'Northwest China': '/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china_north/xibei.shp',
    'North China': '/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china_north/huabei.shp',
    'Northeast China': '/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china_north/dongbei.shp'
}
region_shapes = {name: gpd.read_file(path) for name, path in region_paths.items()}

# Read WRF cases: REAL, CTL, REST
nc_real = nc.Dataset('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/figure_wrf/dust_shao04_2023.nc')
nc_ctl = nc.Dataset('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/figure_wrf/dust_shao04_2005.nc')
nc_rest = nc.Dataset('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/figure_wrf/dust_shao04_ideal_new.nc')

lat = nc_real.variables['lat'][:, :]
lon = nc_real.variables['lon'][:, :]

# Absolute differences (convert to mg/m^2 by dividing by 1000 for dust fields)
dust_real_ctl = np.mean(nc_real.variables['dust'][date_range[0]:date_range[1], :, :], axis=0) / 1000 - \
    np.mean(nc_ctl.variables['dust'][date_range[0]:date_range[1], :, :], axis=0) / 1000
dust_rest_ctl = np.mean(nc_rest.variables['dust'][date_range[0]:date_range[1], :, :], axis=0) / 1000 - \
    np.mean(nc_ctl.variables['dust'][date_range[0]:date_range[1], :, :], axis=0) / 1000

# Region-wise extraction of grid cells falling inside the shapes
data_real_ctl, data_rest_ctl = [], []
for name, shp in region_shapes.items():
    mask = np.array([shp.contains(Point(lon[i, j], lat[i, j])).any()
                     for i in range(lat.shape[0]) for j in range(lat.shape[1])])
    region_data_real_ctl = dust_real_ctl.flatten()[mask]
    region_data_rest_ctl = dust_rest_ctl.flatten()[mask]
    data_real_ctl.append(region_data_real_ctl)
    data_rest_ctl.append(region_data_rest_ctl)

# Use a GridSpec slice for the boxplot panel (kept consistent with original intent)
import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(40, 2, figure=fig)
axe = fig.add_subplot(gs[0:11, 1])

positions = np.arange(len(region_shapes))

# Boxplots: blue for REAL-CTL, red for REST-CTL
axe.boxplot(data_real_ctl, positions=positions - 0.17, widths=0.3, patch_artist=True,
            boxprops=dict(facecolor='blue', alpha=0.5), medianprops=dict(color='black'),
            flierprops=dict(marker='None'))
axe.boxplot(data_rest_ctl, positions=positions + 0.17, widths=0.3, patch_artist=True,
            boxprops=dict(facecolor='red', alpha=0.5), medianprops=dict(color='black'),
            flierprops=dict(marker='None'))

axe.set_xticks(positions)
axe.set_xticklabels(region_shapes.keys(), fontsize=14)
axe.grid(axis='y', linestyle=':')

import matplotlib.patches as mpatches
legend_handles = [
    mpatches.Patch(color='blue', alpha=0.5, label='WRF_REAL - WRF_CTL'),
    mpatches.Patch(color='red', alpha=0.5, label='WRF_REST - WRF_CTL')
]
axe.legend(handles=legend_handles, fontsize=14, loc='upper right', ncol=2)
plt.title('(b) Dust Change [mg/m$^2$]', loc='left', fontsize=20, pad=8)
plt.ylim([-80, 20])
plt.xlim([-0.5, 2.5])

# Minor/major ticks and sizes
axe.minorticks_on()
axe.tick_params(axis="both", which="major", direction="out", width=1.3, length=7)
axe.tick_params(axis="both", which="minor", direction="out", width=1.3, length=3.5)
axe.xaxis.set_minor_locator(mticker.MultipleLocator(5))
axe.yaxis.set_minor_locator(mticker.MultipleLocator(5))
axe.tick_params(axis='both', labelsize=21)

# Light background bands for each region (for visual grouping)
background_colors = ['#4D7EF4', [179/255, 61/255, 145/255], '#0D8B43']
for i, color in enumerate(background_colors):
    axe.axvspan(i - 0.5, i + 0.5, facecolor=color, alpha=0.1)

# Compute region-wise percentage changes (relative to CTL mean within the same mask)
region_mean_real_ctl = []
region_mean_rest_ctl = []
for name, shp in region_shapes.items():
    mask = np.array([shp.contains(Point(lon[i, j], lat[i, j])).any()
                     for i in range(lat.shape[0]) for j in range(lat.shape[1])])
    region_real_ctl = dust_real_ctl.flatten()[mask]
    region_rest_ctl = dust_rest_ctl.flatten()[mask]
    region_ctl = nc_ctl.variables['dust'][date_range[0]:date_range[1], :, :].mean(axis=0).flatten()[mask] / 1000
    mean_real_ctl = (np.mean(region_real_ctl) / np.mean(region_ctl)) * 100 if np.mean(region_ctl) != 0 else np.nan
    mean_rest_ctl = (np.mean(region_rest_ctl) / np.mean(region_ctl)) * 100 if np.mean(region_ctl) != 0 else np.nan
    region_mean_real_ctl.append(mean_real_ctl)
    region_mean_rest_ctl.append(mean_rest_ctl)

# Annotate percentage changes beneath the boxes
for i, (real_ctl, rest_ctl) in enumerate(zip(region_mean_real_ctl, region_mean_rest_ctl)):
    axe.text(i, -72, f'{real_ctl:.1f}%\n', ha='center', fontsize=22, color='blue')
    axe.text(i, -72, f'{rest_ctl:.1f}%', ha='center', fontsize=22, color='red')

# #2##############################################################################################
# (c) WRF_REAL - WRF_CTL dust [mg/m^2]
axe = plt.subplot(323, projection=ccrs.PlateCarree(), aspect="auto")

# Diverging color scheme and levels for dust differences
levels = np.arange(-50, 51, 5)
rgb = [
    [41,42,109], [39,53,126], [30,69,143], [48,101,167], [64,127,181], [81,147,195],
    [108,172,207], [138,198,221], [166,215,232], [197,229,242],[220,233,213],
    [239,222,153], [253,205,103],[252,179,87], [245,149,65],[242,119,52], [239,81,39],
    [232,64,35], [218,44,37], [196,31,38], [167,31,45],[140,21,24],
]
rgb = np.array(rgb) / 255.0

# Read REAL and CTL cases for dust (mg/m^2)
ncfile = nc.Dataset('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/figure_wrf/dust_shao04_2023.nc')
lat = ncfile.variables['lat'][:, :]
lon = ncfile.variables['lon'][:, :]
dust = ncfile.variables['dust'][date_range[0]:date_range[1], :, :]
ncfile1 = nc.Dataset('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/figure_wrf/dust_shao04_2005.nc')
dust1 = ncfile1.variables['dust'][date_range[0]:date_range[1], :, :]

# Filled contours of differences
contourf = plt.contourf(lon, lat, np.mean(dust, 0)/1000 - np.mean(dust1, 0)/1000, levels, colors=rgb, extend='both')

# Title and tick font
plt.title('(c) WRF_REAL - WRF_CTL: Dust [mg/m$^2$]', loc='left', fontsize=20, pad=8)
for tick in axe.get_xticklabels() + axe.get_yticklabels():
    tick.set_fontname('Arial')

# Coastlines and extent
axe.add_feature(cfeat.COASTLINE.with_scale('10m'), linewidth=0, color='k')
axe.set_extent([75, 135, 30, 55], crs=ccrs.PlateCarree())

# Gridlines (labels disabled)
gl = axe.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0, color='gray', linestyle=':')
gl.top_labels, gl.bottom_labels, gl.right_labels, gl.left_labels = False, False, False, False
gl.xlocator = mticker.FixedLocator(np.arange(80, 135, 15))
gl.ylocator = mticker.FixedLocator(np.arange(30, 56, 10))

# Axis ticks and formatters
axe.set_xticks(np.arange(80, 135, 15), crs=ccrs.PlateCarree())
axe.set_yticks(np.arange(30, 56, 10), crs=ccrs.PlateCarree())
axe.xaxis.set_major_formatter(LongitudeFormatter())
axe.yaxis.set_major_formatter(LatitudeFormatter())
axe.tick_params(labelcolor='k', length=5)

# Tick label sizes and general tick params
plt.xticks(fontsize=21)
plt.yticks(fontsize=21)
plt.tick_params(top='on', right='on', which='both')

# Shapefiles for boundaries and regions
shp_world = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/world/world.shp').geometries()
shp_china = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china/china.shp').geometries()
axe.add_geometries(shp_world, ccrs.PlateCarree(), facecolor='none', edgecolor='k', linewidth=0.5, zorder=1)
axe.add_geometries(shp_china, ccrs.PlateCarree(), facecolor='none', edgecolor='k', linewidth=0.5, zorder=1)
xibei_china = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china_north/xibei.shp').geometries()
huabei_china = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china_north/huabei.shp').geometries()
dongbei_china = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china_north/dongbei.shp').geometries()
axe.add_geometries(xibei_china, ccrs.PlateCarree(), facecolor='none', edgecolor='#4D7EF4', linewidth=2.5, linestyle='--', zorder=19)
axe.add_geometries(huabei_china, ccrs.PlateCarree(), facecolor='none', edgecolor=[179/255,61/255,145/255], linewidth=2.5, linestyle='--', zorder=20)
axe.add_geometries(dongbei_china, ccrs.PlateCarree(), facecolor='none', edgecolor='#0D8B43', linewidth=2.5, linestyle='--', zorder=19)

# Minor ticks and colorbar
axe.minorticks_on()
axe.tick_params(axis="both", which="major", direction="out", width=1.3, length=7)
axe.tick_params(axis="both", which="minor", direction="out", width=1.3, length=3.5)
axe.xaxis.set_minor_locator(mticker.MultipleLocator(5))
axe.yaxis.set_minor_locator(mticker.MultipleLocator(5))

cb3 = fig.colorbar(contourf, ax=axe, orientation='horizontal', pad=0.1, shrink=1, aspect=22)
cb3.ax.tick_params(labelsize=17)
cb3.outline.set_edgecolor('black')
cb3.ax.tick_params(axis='x', colors='black')
ticks = np.arange(-50, 50.1, 10)
cb3.set_ticks(ticks)
cb3.set_ticklabels([int(x) for x in ticks])

# #2##############################################################################################
# (d) WRF_REAL - WRF_CTL: Dust Emission [µg/m^2/s]
axe = plt.subplot(324, projection=ccrs.PlateCarree(), aspect="auto")

# Diverging levels for emission differences
levels = np.arange(-5, 5.1, 0.5)
rgb = [
    [41,42,109], [39,53,126], [30,69,143], [48,101,167], [64,127,181], [81,147,195],
    [108,172,207], [138,198,221], [166,215,232], [197,229,242],[220,233,213],
    [239,222,153], [253,205,103],[252,179,87], [245,149,65],[242,119,52], [239,81,39],
    [232,64,35], [218,44,37], [196,31,38], [167,31,45],[140,21,24],
]
rgb = np.array(rgb) / 255.0

# Read emission fields (edust) for REAL and CTL
ncfile = nc.Dataset('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/figure_wrf/dust_shao04_2023.nc')
lat = ncfile.variables['lat'][:, :]
lon = ncfile.variables['lon'][:, :]
dust = ncfile.variables['edust'][date_range[0]:date_range[1], :, :]
ncfile1 = nc.Dataset('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/figure_wrf/dust_shao04_2005.nc')
dust1 = ncfile1.variables['edust'][date_range[0]:date_range[1], :, :]

# Filled contours of emission differences
contourf = plt.contourf(lon, lat, np.mean(dust, 0) - np.mean(dust1, 0), levels, colors=rgb, extend='both')

# Title and tick font
plt.title('(d) WRF_REAL - WRF_CTL: Dust Emission [µg/m$^2$/s]', loc='left', fontsize=20, pad=8)
for tick in axe.get_xticklabels() + axe.get_yticklabels():
    tick.set_fontname('Arial')

# Coastlines and extent (focused box)
axe.add_feature(cfeat.COASTLINE.with_scale('10m'), linewidth=0, color='k')
axe.set_extent([85, 125, 30, 55], crs=ccrs.PlateCarree())

# Gridlines (labels disabled)
gl = axe.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0, color='gray', linestyle=':')
gl.top_labels, gl.bottom_labels, gl.right_labels, gl.left_labels = False, False, False, False
gl.xlocator = mticker.FixedLocator(np.arange(90, 121, 15))
gl.ylocator = mticker.FixedLocator(np.arange(30, 56, 10))

# Axis ticks and formatters
axe.set_xticks(np.arange(90, 121, 15), crs=ccrs.PlateCarree())
axe.set_yticks(np.arange(30, 56, 10), crs=ccrs.PlateCarree())
axe.xaxis.set_major_formatter(LongitudeFormatter())
axe.yaxis.set_major_formatter(LatitudeFormatter())
axe.tick_params(labelcolor='k', length=5)

# Tick label sizes and shapefiles
plt.xticks(fontsize=21)
plt.yticks(fontsize=21)
plt.tick_params(top='on', right='on', which='both')

shp_world = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/world/world.shp').geometries()
shp_china = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china/china.shp').geometries()
axe.add_geometries(shp_world, ccrs.PlateCarree(), facecolor='none', edgecolor='k', linewidth=0.5, zorder=1)
axe.add_geometries(shp_china, ccrs.PlateCarree(), facecolor='none', edgecolor='k', linewidth=0.5, zorder=1)

# Minor ticks and colorbar
axe.minorticks_on()
axe.tick_params(axis="both", which="major", direction="out", width=1.3, length=7)
axe.tick_params(axis="both", which="minor", direction="out", width=1.3, length=3.5)
axe.xaxis.set_minor_locator(mticker.MultipleLocator(5))
axe.yaxis.set_minor_locator(mticker.MultipleLocator(5))

cb3 = fig.colorbar(contourf, ax=axe, orientation='horizontal', pad=0.1, shrink=1, aspect=22)
cb3.ax.tick_params(labelsize=17)
cb3.outline.set_edgecolor('black')
cb3.ax.tick_params(axis='x', colors='black')
ticks = np.arange(-5, 5.1, 1)
cb3.set_ticks(ticks)
cb3.set_ticklabels([int(x) for x in ticks])

# #2##############################################################################################
# (e) WRF_REST - WRF_CTL dust [mg/m^2]
axe = plt.subplot(325, projection=ccrs.PlateCarree(), aspect="auto")

# Diverging color scheme for dust differences
levels = np.arange(-50, 51, 5)
rgb = [
    [41,42,109], [39,53,126], [30,69,143], [48,101,167], [64,127,181], [81,147,195],
    [108,172,207], [138,198,221], [166,215,232], [197,229,242],[220,233,213],
    [239,222,153], [253,205,103],[252,179,87], [245,149,65],[242,119,52], [239,81,39],
    [232,64,35], [218,44,37], [196,31,38], [167,31,45],[140,21,24],
]
rgb = np.array(rgb) / 255.0

# Read REST and CTL for dust (mg/m^2)
ncfile = nc.Dataset('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/figure_wrf/dust_shao04_ideal_new.nc')
lat = ncfile.variables['lat'][:, :]
lon = ncfile.variables['lon'][:, :]
dust = ncfile.variables['dust'][date_range[0]:date_range[1], :, :]
ncfile1 = nc.Dataset('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/figure_wrf/dust_shao04_2005.nc')
dust1 = ncfile1.variables['dust'][date_range[0]:date_range[1], :, :]

# Filled contours of differences
contourf = plt.contourf(lon, lat, np.mean(dust, 0)/1000 - np.mean(dust1, 0)/1000, levels, colors=rgb, extend='both')

# Title and tick font
plt.title('(e) WRF_REST - WRF_CTL: Dust [mg/m$^2$]', loc='left', fontsize=20, pad=8)
for tick in axe.get_xticklabels() + axe.get_yticklabels():
    tick.set_fontname('Arial')

# Coastlines and extent
axe.add_feature(cfeat.COASTLINE.with_scale('10m'), linewidth=0, color='k')
axe.set_extent([75, 135, 30, 55], crs=ccrs.PlateCarree())

# Gridlines (labels disabled)
gl = axe.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0, color='gray', linestyle=':')
gl.top_labels, gl.bottom_labels, gl.right_labels, gl.left_labels = False, False, False, False
gl.xlocator = mticker.FixedLocator(np.arange(80, 135, 15))
gl.ylocator = mticker.FixedLocator(np.arange(30, 56, 10))

# Axis ticks and formatters
axe.set_xticks(np.arange(80, 135, 15), crs=ccrs.PlateCarree())
axe.set_yticks(np.arange(30, 56, 10), crs=ccrs.PlateCarree())
axe.xaxis.set_major_formatter(LongitudeFormatter())
axe.yaxis.set_major_formatter(LatitudeFormatter())
axe.tick_params(labelcolor='k', length=5)

# Tick sizes and shapefiles
plt.xticks(fontsize=21)
plt.yticks(fontsize=21)
plt.tick_params(top='on', right='on', which='both')

shp_world = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/world/world.shp').geometries()
shp_china = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china/china.shp').geometries()
axe.add_geometries(shp_world, ccrs.PlateCarree(), facecolor='none', edgecolor='k', linewidth=0.5, zorder=1)
axe.add_geometries(shp_china, ccrs.PlateCarree(), facecolor='none', edgecolor='k', linewidth=0.5, zorder=1)
xibei_china = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china_north/xibei.shp').geometries()
huabei_china = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china_north/huabei.shp').geometries()
dongbei_china = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china_north/dongbei.shp').geometries()
axe.add_geometries(xibei_china, ccrs.PlateCarree(), facecolor='none', edgecolor='#4D7EF4', linewidth=2.5, linestyle='--', zorder=19)
axe.add_geometries(huabei_china, ccrs.PlateCarree(), facecolor='none', edgecolor=[179/255,61/255,145/255], linewidth=2.5, linestyle='--', zorder=20)
axe.add_geometries(dongbei_china, ccrs.PlateCarree(), facecolor='none', edgecolor='#0D8B43', linewidth=2.5, linestyle='--', zorder=19)

# Minor ticks and colorbar
axe.minorticks_on()
axe.tick_params(axis="both", which="major", direction="out", width=1.3, length=7)
axe.tick_params(axis="both", which="minor", direction="out", width=1.3, length=3.5)
axe.xaxis.set_minor_locator(mticker.MultipleLocator(5))
axe.yaxis.set_minor_locator(mticker.MultipleLocator(5))

cb3 = fig.colorbar(contourf, ax=axe, orientation='horizontal', pad=0.1, shrink=1, aspect=22)
cb3.ax.tick_params(labelsize=17)
cb3.outline.set_edgecolor('black')
cb3.ax.tick_params(axis='x', colors='black')
ticks = np.arange(-50, 51, 10)
cb3.set_ticks(ticks)
cb3.set_ticklabels([int(x) for x in ticks])

# #2##############################################################################################
# (f) WRF_REST - WRF_CTL: Dust Emission [µg/m^2/s]
axe = plt.subplot(326, projection=ccrs.PlateCarree(), aspect="auto")

# Diverging levels for emission differences
levels = np.arange(-5, 5.1, 0.5)
rgb = [
    [41,42,109], [39,53,126], [30,69,143], [48,101,167], [64,127,181], [81,147,195],
    [108,172,207], [138,198,221], [166,215,232], [197,229,242],[220,233,213],
    [239,222,153], [253,205,103],[252,179,87], [245,149,65],[242,119,52], [239,81,39],
    [232,64,35], [218,44,37], [196,31,38], [167,31,45],[140,21,24],
]
rgb = np.array(rgb) / 255.0

# Read emission fields (edust) for REST and CTL
ncfile = nc.Dataset('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/figure_wrf/dust_shao04_ideal_new.nc')
lat = ncfile.variables['lat'][:, :]
lon = ncfile.variables['lon'][:, :]
dust = ncfile.variables['edust'][date_range[0]:date_range[1], :, :]
ncfile1 = nc.Dataset('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/figure_wrf/dust_shao04_2005.nc')
dust1 = ncfile1.variables['edust'][date_range[0]:date_range[1], :, :]

# Filled contours of emission differences
contourf = plt.contourf(lon, lat, np.mean(dust, 0) - np.mean(dust1, 0), levels, colors=rgb, extend='both')

# Title and tick font
plt.title('(f) WRF_REST - WRF_CTL: Dust Emission [µg/m$^2$/s]', loc='left', fontsize=20, pad=8)
for tick in axe.get_xticklabels() + axe.get_yticklabels():
    tick.set_fontname('Arial')

# Coastlines and extent (focused box)
axe.add_feature(cfeat.COASTLINE.with_scale('10m'), linewidth=0, color='k')
axe.set_extent([85, 125, 30, 55], crs=ccrs.PlateCarree())

# Gridlines (labels disabled)
gl = axe.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0, color='gray', linestyle=':')
gl.top_labels, gl.bottom_labels, gl.right_labels, gl.left_labels = False, False, False, False
gl.xlocator = mticker.FixedLocator(np.arange(90, 121, 15))
gl.ylocator = mticker.FixedLocator(np.arange(30, 56, 10))

# Axis ticks and formatters
axe.set_xticks(np.arange(90, 121, 15), crs=ccrs.PlateCarree())
axe.set_yticks(np.arange(30, 56, 10), crs=ccrs.PlateCarree())
axe.xaxis.set_major_formatter(LongitudeFormatter())
axe.yaxis.set_major_formatter(LatitudeFormatter())
axe.tick_params(labelcolor='k', length=5)

# Tick sizes and shapefiles
plt.xticks(fontsize=21)
plt.yticks(fontsize=21)
plt.tick_params(top='on', right='on', which='both')

shp_world = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/world/world.shp').geometries()
shp_china = shpreader.Reader('/data/groups/g1600002/home/qiupch2023/lustre_data/EST_2/shp/china/china.shp').geometries()
axe.add_geometries(shp_world, ccrs.PlateCarree(), facecolor='none', edgecolor='k', linewidth=0.5, zorder=1)
axe.add_geometries(shp_china, ccrs.PlateCarree(), facecolor='none', edgecolor='k', linewidth=0.5, zorder=1)

# Minor ticks and colorbar
axe.minorticks_on()
axe.tick_params(axis="both", which="major", direction="out", width=1.3, length=7)
axe.tick_params(axis="both", which="minor", direction="out", width=1.3, length=3.5)
axe.xaxis.set_minor_locator(mticker.MultipleLocator(5))
axe.yaxis.set_minor_locator(mticker.MultipleLocator(5))

# Optional rectangle (kept commented as in original)
# rect = patches.Rectangle((98, 42.5), 15.5, 4.5, linewidth=2, edgecolor='k', alpha=0.5, facecolor='none', linestyle='--')
# axe.add_patch(rect)

cb3 = fig.colorbar(contourf, ax=axe, orientation='horizontal', pad=0.1, shrink=1, aspect=22)
cb3.ax.tick_params(labelsize=17)
cb3.outline.set_edgecolor('black')
cb3.ax.tick_params(axis='x', colors='black')
ticks = np.arange(-5, 5.1, 1)
cb3.set_ticks(ticks)
cb3.set_ticklabels([int(x) for x in ticks])

# Final spacing and save
# plt.tight_layout()
plt.subplots_adjust(left=0.06, bottom=0.03, right=0.95, top=0.95, wspace=0.2, hspace=0.1)
plt.savefig('figure3_Arial.png', dpi=500)
