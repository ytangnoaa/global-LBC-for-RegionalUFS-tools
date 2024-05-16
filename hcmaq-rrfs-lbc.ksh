#!/bin/ksh -x
##. /etc/profile
##module load intel/18.0.5.274  impi/2018.0.4 netcdf/4.7.0 nco

export  SALLOC_ACCOUNT=arl
export  SBATCH_ACCOUNT=arl
export  SLURM_QOS=debug

mon=2
year=2018

while [ $mon -le 12 ]; do

amon=$mon

typeset -Z2 amon

#ORIGDIR=/scratch1/NCEPDEV/stmp4/Jianping.Huang/expt_dirs/test_community$mon/$PDY$cyc/INPUT
#ORIGDIR=sample-input-C791/2019080100/INPUT/
#if [ ! -s $ORIGDIR/gfs_bndy.tile7.000.nc ]; then
# echo " no original LBC file $ORIGDIR/gfs_bndy.tile7.000.nc "
# exit 1
#fi

TDIR=INPUT-C793-$year$amon
if [ ! -s $TDIR ] ;then
 mkdir $TDIR
fi

cp -p INPUT-C793-2013$amon/geos5_bndy.c793.2013$mon.nc $TDIR/hcmaq_bndy.c793.$year$amon.v1.nc

export HCMAQ=/scratch2/NCEPDEV/naqfc/Youhua.Tang/HEMIS-CMAQ-2018/CCTM_CONC_v532_intel18.0_CMAQv53_TS_108NHEMI_2018${amon}_monthlyav.nc

cat > hcmaq-rrfs-lbc.ini <<EOF
&control
 iprint=1
 lbcfile='$TDIR/hcmaq_bndy.c793.$year$amon.v1.nc'
 topofile='/scratch2/NCEPDEV/naqfc/RRFS_CMAQ/DOMAIN_DATA/AQM_NA_13km/C793_oro_data.tile7.halo4.nc'
 bndname='no2','no','o3','no3','oh','ho2','n2o5','hno3','nh3','hono','pna','meoh','etoh',
  'h2o2','co','so2','pan','panx','mepx','mgly','cl','cl2','hcl','clo','hocl','facd','aacd','pacd',
  'open','etha','cres','form','ald2','aldx','rooh','par','ole','isop','ispd','eth','xylmn',
  'xo2','o1d','c2o3','prpa',
  'anai','acli','aso4j','ano3j','anh4j','aecj','asoil','anaj','aseacat','aothrj','apocj',
  'aclj','aclk','numatkn','numacc','numcor'
 o3cap=0.11
 bk_numatkn=9.583000e+07, 9.583000e+07, 9.583000e+07, 9.583000e+07, 9.583000e+07,
       9.583000e+07, 9.583000e+07, 9.583000e+07, 9.583000e+07, 9.583000e+07,
       9.583000e+07, 9.583000e+07, 9.583000e+07, 9.583000e+07, 9.583000e+07,
       9.583000e+07, 9.583000e+07, 9.583000e+07, 9.583000e+07, 9.583000e+07,
       9.583000e+07, 9.560066e+07, 9.526129e+07, 9.507742e+07, 9.507114e+07,
       9.524665e+07, 9.557631e+07, 9.598170e+07, 9.630892e+07, 9.632580e+07,
       9.727337e+07, 1.020295e+08, 1.072891e+08, 1.219373e+08, 1.358889e+08,
       1.469180e+08, 1.537770e+08, 1.555209e+08, 1.911703e+08, 2.729436e+08,
       3.282245e+08, 5.999171e+08, 8.832315e+08, 1.137202e+09, 1.342505e+09,
       1.518161e+09, 1.675616e+09, 1.804869e+09, 1.893318e+09, 1.923304e+09,
       1.917325e+09, 1.916367e+09, 1.917183e+09, 1.916952e+09, 1.917007e+09,
       1.916999e+09, 1.917000e+09, 1.917000e+09, 1.917000e+09, 1.917000e+09,
       1.917000e+09, 1.917000e+09, 1.917000e+09, 1.917000e+09, 1.917000e+09
 bk_numacc=1.416000e+07, 1.416000e+07, 1.416000e+07, 1.416000e+07, 1.416000e+07,
       1.416000e+07, 1.416000e+07, 1.416000e+07, 1.416000e+07, 1.416000e+07,
       1.416000e+07, 1.416000e+07, 1.416000e+07, 1.416000e+07, 1.416000e+07,
       1.416000e+07, 1.416000e+07, 1.416000e+07, 1.416000e+07, 1.416000e+07,
       1.416000e+07, 1.412615e+07, 1.407606e+07, 1.404892e+07, 1.404800e+07,
       1.407390e+07, 1.412256e+07, 1.418239e+07, 1.423069e+07, 1.423318e+07,
       1.437304e+07, 1.507502e+07, 1.585132e+07, 1.801333e+07, 2.007252e+07,
       2.170036e+07, 2.271280e+07, 2.296984e+07, 2.823377e+07, 4.030898e+07,
       4.847186e+07, 8.859494e+07, 1.304345e+08, 1.679405e+08, 1.982594e+08,
       2.242000e+08, 2.474527e+08, 2.665407e+08, 2.796027e+08, 2.840309e+08,
       2.831480e+08, 2.830065e+08, 2.831270e+08, 2.830930e+08, 2.831011e+08,
       2.830998e+08, 2.831000e+08, 2.831000e+08, 2.831000e+08, 2.831000e+08,
       2.831000e+08, 2.831000e+08, 2.831000e+08, 2.831000e+08, 2.831000e+08
 bk_numcor=1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
       1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
       1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
       1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. 
 ak_gfs = 0., 20, 64.247, 137.79, 221.958, 318.266, 428.434, 554.424, 698.457, 863.058, 
    1051.08, 1265.752, 1510.711, 1790.051, 2108.366, 2470.788, 2883.038, 
    3351.46, 3883.052, 4485.493, 5167.146, 5937.05, 6804.874, 7777.15, 
    8832.537, 9936.614, 11054.85, 12152.94, 13197.07, 14154.32, 14993.07, 
    15683.49, 16197.97, 16511.74, 16611.6, 16503.14, 16197.32, 15708.89, 
    15056.34, 14261.43, 13348.67, 12344.49, 11276.35, 10171.71, 9057.051, 
    7956.908, 6893.117, 5884.206, 4945.029, 4086.614, 3316.217, 2637.553, 
    2051.15, 1554.789, 1143.988, 812.489, 552.72, 356.223, 214.015, 116.899, 
    55.712, 21.516, 5.741, 0.575, 0, 0
 bk_gfs = 0., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    3.697e-05, 0.00043106, 0.00163591, 0.00410671, 0.00829402, 0.01463712, 
    0.02355588, 0.03544162, 0.05064684, 0.06947458, 0.09216691, 0.1188122, 
    0.1492688, 0.1832962, 0.2205702, 0.2606854, 0.3031641, 0.3474685, 
    0.3930182, 0.4392108, 0.4854433, 0.5311348, 0.5757467, 0.6187996, 
    0.659887, 0.6986829, 0.7349452, 0.7685147, 0.7993097, 0.8273188, 
    0.8525907, 0.8752236, 0.895355, 0.913151, 0.9287973, 0.9424911, 
    0.9544341, 0.9648276, 0.9738676, 0.9817423, 0.9886266, 0.9946712, 1      
&end
EOF

hcmaq-rrfs-lbc.x

if [ $(stat -Lc %s $TDIR/hcmaq_bndy.c793.$year$amon.v1.nc) -gt 180000000 ]; then
 ncks -C -O -x -v i_bottom,j_bottom,i_top,j_top,i_right,j_right,i_left,j_left,ps_bottom,ps_top,ps_right,ps_left,\
t_bottom,t_top,t_right,t_left,w_bottom,w_top,w_right,w_left,zh_bottom,zh_top,zh_right,zh_left,sphum_bottom,sphum_top,sphum_right,sphum_left,\
liq_wat_bottom,liq_wat_top,liq_wat_right,liq_wat_left,o3mr_bottom,o3mr_top,o3mr_right,o3mr_left,ice_wat_bottom,ice_wat_top,ice_wat_right,ice_wat_left,\
rainwat_bottom,rainwat_top,rainwat_right,rainwat_left,snowwat_bottom,snowwat_top,snowwat_right,snowwat_left,graupel_bottom,graupel_top,graupel_right,graupel_left,\
i_w_bottom,j_w_bottom,i_w_top,j_w_top,i_w_right,j_w_right,i_w_left,j_w_left,u_w_bottom,u_w_top,u_w_right,u_w_left,v_w_bottom,v_w_top,v_w_right,v_w_left,\
i_s_bottom,j_s_bottom,i_s_top,j_s_top,i_s_right,j_s_right,i_s_left,j_s_left,u_s_bottom,u_s_top,u_s_right,u_s_left,v_s_bottom,v_s_top,v_s_right,v_s_left\
 $TDIR/hcmaq_bndy.c793.$year$amon.v1.nc /scratch2/NCEPDEV/naqfc/RRFS_CMAQ_NA13km/LBCS/RRFS_NA13km_HCMAQ_v1/hcmaq_bndy.c793.$year$amon.v1.nc
 

else
 echo "HCMAQ LBC failed"
 exit 1
fi   

let mon=$mon+1
done
