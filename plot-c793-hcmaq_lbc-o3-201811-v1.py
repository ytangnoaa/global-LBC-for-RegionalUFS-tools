#!/usr/bin/env python
import xarray as xr
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.colors as colors
import matplotlib.backends.backend_pdf as bpdf

iflag=[True,True]
#stime=datetime.strptime('2013080712','%Y%m%d%H') # in datetime obj

myvalues=[20,30,40,50,60,70,80,100]

#nowstep=4
#nt=nowstep-1
pdf=bpdf.PdfPages('figures/c793-201811-hcmaq_lbc-v1-o3-1.pdf')

plt.figure(figsize=(12,7))
a=xr.open_dataset('INPUT-C793-201811/hcmaq_bndy.c793.201811.v1.nc')

# North LBC 
if iflag[0]:
  cs3=plt.pcolormesh(a.lon,a.zh_top[:,0,:],a.o3_top[:,0,:]*1000,cmap='jet',edgecolor='face',
   norm=colors.BoundaryNorm(boundaries=myvalues,ncolors=256))
  plt.ylim(0,20000)
#  plt.xlim(c.lon.max(),0)
  cbar=plt.colorbar(cs3)
  cbar.set_label('O$_3$ (ppbv)',fontsize=18)
  cbar.ax.tick_params(direction='in',labelsize=16)
  plt.xticks(fontsize=16)
  plt.yticks(fontsize=16) 
  plt.xlabel('C793 Model North Boundary Periphery (grid)',fontsize=18)
  plt.ylabel('Altitude (m)',fontsize=18)
  plt.title('HCMAQ O$_3$ north Lateral Boundary Condition in Nov',fontsize=20)
  plt.tight_layout()
#  plt.savefig('figures/c793-201811-hcmaq_lbc-v1-o3-north.png',dpi=300)
  pdf.savefig(plt.figure(1))
  plt.clf()
# West LBC
if iflag[1]:
  cs3=plt.pcolormesh(a.lat,a.zh_left[:,:,0],a.o3_left[:,:,0]*1000,cmap='jet',edgecolor='face',
   norm=colors.BoundaryNorm(boundaries=myvalues,ncolors=256))
  plt.ylim(0,20000)
#  plt.xlim(c.lat.max(),0)
  cbar=plt.colorbar(cs3)
  cbar.set_label('O$_3$ (ppbv)',fontsize=18)
  cbar.ax.tick_params(direction='in',labelsize=16)
  plt.xticks(fontsize=16)
  plt.yticks(fontsize=16) 
  plt.xlabel('C793 Model West Boundary Periphery (grid)',fontsize=18)
  plt.ylabel('Altitude (m)',fontsize=18)
  plt.title('HCMAQ O$_3$ West Lateral Boundary Condition in Nov',fontsize=20)
  plt.tight_layout()
#  plt.savefig('figures/c793-201811-hcmaq_lbc-v1-o3-west.png',dpi=300)
  pdf.savefig(plt.figure(1))
  plt.clf()


pdf.close() 
