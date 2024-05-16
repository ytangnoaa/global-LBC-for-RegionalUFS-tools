#!/usr/bin/env python3
import xarray as xr
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
myvalue=[10,20,30,40,50,60,70,80]
c=xr.open_dataset('/scratch2/NCEPDEV/naqfc/Youhua.Tang/GFDL_AM4/tracer_level.201901-201912.all.nc')
plt.figure(figsize=(12,7))
for n in range(12):
  month=str(n+1).zfill(2)
  
  ax=plt.axes(projection=ccrs.PlateCarree())
  cs3=plt.pcolormesh(c.lon,c.lat,c.O3[n,-1,:,:]*1e9,
    cmap='jet',norm=colors.BoundaryNorm(boundaries=myvalue,ncolors=256))
  plt.axis([-180,-40,10,80]) # plt.axis([180,320,10,80])
  plt.title('AM4 Surface O3 (ppbv) In 2019'+month)
  cbar=plt.colorbar(cs3,fraction=0.07,orientation='vertical',ticks=myvalue)
  cbar.set_label('O$_3$ (ppbv)',size=18)
  cbar.ax.tick_params(direction='in',labelsize=10)
  ax.add_feature(cfeature.BORDERS)
  ax.add_feature(cfeature.COASTLINE)
  ax.tick_params(which='major', width=1.00, length=5)
  ax.tick_params(which='major', width=1.00, length=5)
  gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, 
    linestyle='--')
  gl.top_labels = False
  gl.right_labels = False
  gl.xformatter = LONGITUDE_FORMATTER
  gl.yformatter = LATITUDE_FORMATTER
  gl.ylabel_style = {'size':12, 'weight':'bold'}
  gl.xlabel_style = {'size':12, 'weight':'bold'}
  plt.tight_layout()
  plt.savefig('figures/am4_2019'+month+'_surface_o3.png',dpi=300)
  plt.clf()
