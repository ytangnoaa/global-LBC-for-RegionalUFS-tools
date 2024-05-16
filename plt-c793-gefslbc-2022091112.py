#!/usr/bin/env python3

import xarray as xr
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.colors as colors

iflag=[True,True]
ndatetime='2022090112'
nfolder='/scratch2/NCEPDEV/naqfc/Kai.Wang/Diag_LBC/files_fromwcoss2/lfs/h2/emc/ptmp/jianping.huang/para/com/aqm/v7.0/aqm.v7.0.c4/With_GEFS_Aero_LBCs/'
case='c3'
stime=datetime.strptime(ndatetime,'%Y%m%d%H') # in datetime obj

myvalues=[1,2,3,5,10,15,20,30,50]

plt.figure(figsize=(12,7))



for nt in range(0,30,6):
 
 step4=stime+timedelta(hours=nt)
 nowstring=step4.strftime('%HZ, %m/%d/%Y')
 ftstring=step4.strftime('%y%m%d%H')
# nowfile='INPUT/gfs_bndy.tile7.'+str(nt).zfill(3)+'.nc'
 
 nowfile=nfolder+ndatetime+'/INPUT/gfs_bndy.tile7.'+str(nt).zfill(3)+'.nc'
 print('nowfile=',nowfile)
 c=xr.open_dataset(nowfile)
# North LBC 
 if iflag[0]:
  cs3=plt.pcolormesh(c.lon,c.zh_top[:,0,:],c.aecj_top[:,0,:]+c.aorgcj_top[:,0,:],
   cmap='jet',norm=colors.BoundaryNorm(boundaries=myvalues,ncolors=256))
#  plt.ylim(0,20000)
  plt.xlim(c.lon.max(),0)
  cbar=plt.colorbar(cs3)
  cbar.set_label('AECJ+AORGCJ ($\mu$g/kg)',fontsize=18)
  cbar.ax.tick_params(direction='in',labelsize=16)
  plt.xticks(fontsize=16)
  plt.yticks(fontsize=16) 
  plt.xlabel('C793 Model North Boundary Periphery (grid)',fontsize=18)
  plt.ylabel('Altitude (m)',fontsize=18)
  plt.title(case+' GEFS AECJ+AORGCJ North Lateral Boundary Condition at '+nowstring+'\n',fontsize=20)
  plt.tight_layout()
  plt.savefig('figures/'+case+'-c793-'+ndatetime+'-gefslbc-ecoc-north-'+ftstring+'.png',dpi=300)
  plt.clf()
# South LBC
 if iflag[1]:
  cs3=plt.pcolormesh(c.lon,c.zh_bottom[:,0,:],c.aecj_bottom[:,0,:]+c.aorgcj_bottom[:,0,:],
   cmap='jet',norm=colors.BoundaryNorm(boundaries=myvalues,ncolors=256))
#  plt.ylim(0,20000)
  plt.xlim(c.lon.max(),0)
  cbar=plt.colorbar(cs3)
  cbar.set_label('AECJ+AORGCJ ($\mu$g/kg)',fontsize=18)
  cbar.ax.tick_params(direction='in',labelsize=16)
  plt.xticks(fontsize=16)
  plt.yticks(fontsize=16) 
  plt.xlabel('C793 Model South Boundary Periphery (grid)',fontsize=18)
  plt.ylabel('Altitude (m)',fontsize=18)
  plt.title(case+' GEFS AECJ+AORGCJ South Lateral Boundary Condition at '+nowstring+'\n',fontsize=20)
  plt.tight_layout()
  plt.savefig('figures/'+case+'-c793-'+ndatetime+'-gefslbc-ecoc-south-'+ftstring+'.png',dpi=300)
  plt.clf()
