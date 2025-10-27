import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
from scipy.stats import pearsonr
from netCDF4 import Dataset

# -----------------------------------------------------------------------------
# Global style
# -----------------------------------------------------------------------------
plt.rcParams.update({'font.size': 15, 'font.family': 'Arial'})
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 1.3

# -----------------------------------------------------------------------------
# Load .mat data
# -----------------------------------------------------------------------------
data = scipy.io.loadmat('相关性分析数据.mat')

# Original variables (flattened)
raw_data = [
    data['ndviFiltered'].flatten(),
    data['eviFiltered'].flatten(),
    data['laiFiltered'].flatten(),
    data['albedoFiltered'].flatten()
]
var_names = ['NDVI', 'EVI', 'LAI', 'Albedo']

# -----------------------------------------------------------------------------
# Replace NDVI/Albedo for the first panel using NetCDF file
# -----------------------------------------------------------------------------
point_nc = Dataset('./selected_points_data.nc')
ndvi_from_nc = point_nc.variables['ndvi'][:].flatten()
albedo_from_nc = point_nc.variables['albedo'][:].flatten()
point_nc.close()

# -----------------------------------------------------------------------------
# Build x,y lists for three panels (all UN-normalized, as-is)
# Panel (a): NDVI(from nc) vs Albedo(from nc)
# Panel (b): EVI(from .mat) vs Albedo(from .mat)
# Panel (c): LAI(from .mat) vs Albedo(from .mat)
# -----------------------------------------------------------------------------
x_data = [ndvi_from_nc, raw_data[1], raw_data[2]]
y_data = [albedo_from_nc, raw_data[3], raw_data[3]]

# Titles/labels
titles = ['(a)', '(b)', '(c)']
x_labels = ['NDVI', 'EVI', 'LAI']
y_labels = ['Albedo'] * 3

# -----------------------------------------------------------------------------
# Figure and axes
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Common scatter size per panel
sizes = [10, 10, 20]

# Placeholder for colorbar handle (last scatter plotted)
last_scatter = None

# -----------------------------------------------------------------------------
# Plot each panel
# -----------------------------------------------------------------------------
for i in range(3):
    ax = axes[i]
    x = np.asarray(x_data[i]).astype(float)
    y = np.asarray(y_data[i]).astype(float)

    # Keep only finite values and in-range [0,1] for both axes
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if x.size == 0 or y.size == 0:
        # If nothing left, just set up blank axes
        ax.set_title(titles[i], fontsize=20, loc='left')
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel(x_labels[i], fontsize=18)
        ax.set_ylabel(y_labels[i], fontsize=18)
        continue

    # Clip to [0,1] (optional: keep consistency with histogram range)
    x = np.clip(x, 0, 1)
    y = np.clip(y, 0, 1)

    # Pearson correlation and p-value
    corr_coeff, p_value = pearsonr(x, y)
    if p_value < 0.001:
        p_text = '$p < 0.001$'
    elif p_value < 0.01:
        p_text = '$p < 0.01$'
    else:
        p_text = f'$p = {p_value:.2f}$'

    # 2D histogram-based density (range fixed to [0,1]×[0,1])
    density, xedges, yedges = np.histogram2d(x, y, bins=30, range=[[0, 1], [0, 1]])
    density = density.T
    xcent = 0.5 * (xedges[:-1] + xedges[1:])
    ycent = 0.5 * (yedges[:-1] + yedges[1:])
    # Simple 1D interpolation of mean density by x-bin (kept same approach as your code)
    dens_color = np.interp(x, xcent, density.mean(axis=0))

    # Scatter colored by density
    last_scatter = ax.scatter(
        x, y,
        c=dens_color,
        s=sizes[i],
        cmap='jet',
        alpha=0.7,
        edgecolors='k',
        linewidths=0,
        marker='o'
    )

    # Linear fit (degree 1)
    coeffs = np.polyfit(x, y, deg=1)
    poly_eq = np.poly1d(coeffs)
    x_fit = np.linspace(0, 1, 200)
    y_fit = poly_eq(x_fit)
    ax.plot(x_fit, y_fit, color='red', linewidth=2)

    # Text box: regression equation, r, p
    ax.text(
        0.50, 0.78,
        f'$y = {coeffs[0]:.2f}x + {coeffs[1]:.2f}$\n'
        f'$r = {corr_coeff:.2f}$\n{p_text}',
        fontsize=17,
        transform=ax.transAxes
    )

    # Ticks & grid
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", direction="out", width=1.3, length=4)
    ax.tick_params(axis="both", which="minor", direction="out", width=1.3, length=0)

    # Minor tick spacing suited for [0,1]
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(0.1))

    # Axes cosmetics
    ax.set_title(titles[i], fontsize=20, loc='left')
    ax.set_xlabel(x_labels[i], fontsize=18)
    ax.set_ylabel(y_labels[i], fontsize=18)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(True, linestyle='--', alpha=0.6)

# Layout
plt.subplots_adjust(left=0.05, bottom=0.26, right=0.98, top=0.92, wspace=0.2, hspace=0.1)

# Horizontal colorbar (density)
cbar_ax = fig.add_axes([0.25, 0.11, 0.5, 0.03])
if last_scatter is not None:
    fig.colorbar(last_scatter, cax=cbar_ax, orientation='horizontal', label='Density')

# Save
plt.savefig('figure_v_choose.png', dpi=500)
