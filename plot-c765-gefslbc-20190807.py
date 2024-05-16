#!/usr/bin/env python
import xarray as xr
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.colors as colors
#import matplotlib.backends.backend_pdf as bpdf

iflag=[True,True]
stime=datetime.strptime('2019080712','%Y%m%d%H') # in datetime obj

myvalues=[1,2,3,5,10,15,20,30,50]

#nowstep=4
#nt=nowstep-1

plt.figure(figsize=(12,7))

for nt in range(0,49,6):
 
 step4=stime+timedelta(hours=nt)
 nowstring=step4.strftime('%HZ, %m/%d/%Y')
 ftstring=step4.strftime('%y%m%d%H')
# nowfile='INPUT/gfs_bndy.tile7.'+str(nt).zfill(3)+'.nc'
 nowfile='/scratch2/NCEPDEV/stmp3/Chan-hoo.Jeon/expt_dirs/uwm_aqm_conus13_2days/2019080712/INPUT/gfs_bndy.tile7.'+str(nt).zfill(3)+'.nc'
 print('nowfile=',nowfile)
 c=xr.open_dataset(nowfile)
# North LBC 
 if iflag[0]:
  cs3=plt.pcolormesh(c.lon,c.zh_top[:,0,:],c.aecj_top[:,0,:]+c.aorgcj_top[:,0,:],cmap='jet',norm=colors.BoundaryNorm(boundaries=myvalues,ncolors=256))
  plt.ylim(0,20000)
  plt.xlim(c.lon.max(),0)
  cbar=plt.colorbar(cs3)
  cbar.set_label('AECJ+AORGCJ ($\mu$g/kg)',fontsize=18)
  cbar.ax.tick_params(direction='in',labelsize=16)
  plt.xticks(fontsize=16)
  plt.yticks(fontsize=16) 
  plt.xlabel('C765 Model North Boundary Periphery (grid)',fontsize=18)
  plt.ylabel('Altitude (m)',fontsize=18)
  plt.title('GEFS AECJ+AORGCJ north Lateral Boundary Condition at '+nowstring+'\n',fontsize=20)
  plt.tight_layout()
  plt.savefig('figures/c765-2019080712-gefslbc-ecoc-north-'+ftstring+'.png',dpi=300)
  plt.clf()
# East LBC
 if iflag[1]:
  cs3=plt.pcolormesh(c.lat,c.zh_right[:,:,0],c.asoil_right[:,:,0],cmap='jet',norm=colors.BoundaryNorm(boundaries=myvalues,ncolors=256))
  plt.ylim(0,20000)
  plt.xlim(c.lat.max(),0)
  cbar=plt.colorbar(cs3)
  cbar.set_label('ASOIL ($\mu$g/kg)',fontsize=18)
  cbar.ax.tick_params(direction='in',labelsize=16)
  plt.xticks(fontsize=16)
  plt.yticks(fontsize=16) 
  plt.xlabel('C401 Model East Boundary Periphery (grid)',fontsize=18)
  plt.ylabel('Altitude (m)',fontsize=18)
  plt.title('GEFS ASOIL East Lateral Boundary Condition at '+nowstring+'\n',fontsize=20)
  plt.tight_layout()
  plt.savefig('figures/c765-2019080712-gefslbc-east-'+ftstring+'.png',dpi=300)
  plt.clf()


#pdf.close() 
