"""Shared basemap for the lower-48 array maps. Cartopy cache is pinned into the
workspace because the user's ~/.local is not writable from the sandbox."""
import os
import cartopy
cartopy.config['data_dir'] = os.path.abspath('work/cartopy')
import cartopy.crs as ccrs
import cartopy.feature as cfeature

AEA = ccrs.AlbersEqualArea(central_longitude=-96, standard_parallels=(29.5, 45.5))
PC = ccrs.PlateCarree()
LAND = "#f4f2ee"
LINE = "#cfc9c0"
COAST = "#9a938a"


def basemap(ax, extent=(-125, -66, 24, 49.5), lw_state=.35):
    ax.set_extent(extent, PC)
    ax.add_feature(cfeature.LAND, fc=LAND, zorder=0)
    ax.add_feature(cfeature.OCEAN, fc="white", zorder=0)
    ax.add_feature(cfeature.STATES.with_scale('50m'), ec=LINE, lw=lw_state, zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), ec=COAST, lw=.4, zorder=1)
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), ec=COAST, lw=.4, zorder=1)
    ax.spines['geo'].set_visible(False)
    return ax
