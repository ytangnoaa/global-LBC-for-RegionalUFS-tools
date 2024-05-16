#!/bin/ksh -xa
. /etc/profile
module load intel/18.0.5.274  impi/2018.0.4 netcdf/4.7.0 nco
#PDY=20190801

export  SALLOC_ACCOUNT=arl
export  SBATCH_ACCOUNT=arl
export  SLURM_QOS=debug

year=2013
mon=1

while [ $mon -le 12 ]; do

amon=$mon

typeset -Z2 amon
PDY=2020${amon}01
cyc=00

ORIGDIR=/scratch1/NCEPDEV/stmp4/Jianping.Huang/expt_dirs/test_community$mon/$PDY$cyc/INPUT


geoscyc=12
let tstepdiff=$cyc-$geoscyc

#ORIGDIR=sample-input-C791/2019080100/INPUT/

if [ ! -s $ORIGDIR/gfs_bndy.tile7.000.nc ]; then
 echo " no original LBC file $ORIGDIR/gfs_bndy.tile7.000.nc "
 exit 1
fi

TDIR=INPUT-C793-$year$amon
if [ ! -s $TDIR ] ;then
 mkdir $TDIR
fi

cp -p $ORIGDIR/gfs_bndy.tile7.000.nc $TDIR/

NUMTS=1

cat > geos-rrfs-lbc-para.ini <<EOF
&control
 iprint=1
 tstepdiff=$tstepdiff
 dtstep=0 
 bndname='aothrj','aecj','aorgcj','asoil','numacc','numcor'
 prefix='/scratch2/NCEPDEV/naqfc/Youhua.Tang/geos5-nc/data/CGG.tavg24_3d_dac_Nv.monthly.$year${amon}',
    '/scratch2/NCEPDEV/naqfc/Youhua.Tang/geos5-nc/data/CGG.tavg24_3d_dae_Nv.monthly.$year${amon}',
    '/scratch2/NCEPDEV/naqfc/Youhua.Tang/geos5-nc/data/CGG.inst3_3d_aer_Nv.monthly.$year${amon}',
 suffix='.nc4','.nc4','.nc4'
 lbcfile='$TDIR/gfs_bndy.tile7.','.nc'
 topofile='$ORIGDIR/oro_data.tile7.halo4.nc'
 bndname='no2','no','o3','no3','oh','ho2','n2o5','hno3','nh3','hono','pna','meoh','etoh',
  'h2o2','co','so2','pan','panx','mepx','mgly','cl','cl2','hcl','clo','hocl','facd','aacd','pacd',
  'open','etha','cres','form','ald2','aldx','rooh','par','ole','isop','ispd','isopx','isopo2','eth','xyl',
  'xo2','o1d','c2o3','ntri','ntrio2','maco3','prpa',
  'anai','acli','aso4j','ano3j','anh4j','aecj','asoil','anaj','aseacat','aothrj','apocj',
  'aclj','aclk','numatkn','numacc','numcor'
 o3cap=1.0
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
&end

Species converting Factor
# Global ppbv to regional ppbv
'CH2O'    1
'form'   1.0
'CO'      1
'co'     1.0
'HCOOH'   1
'facd'   1.0
'HNO2'    1
'hono'   1.0
'HNO3'    1
'hno3'   1.0
'HNO4'    1
'pna'    1.0
'HO2'     1
'ho2'    1.0
'H2O2'    1
'h2o2'   1.0
'MO2'     1   #CH3O2 
'xo2'    1.0  # RO2 for NO->NO2
'MOH'     1  # methanol
'meoh'   1.0
'MP'      1 # methylhydroperoxide
'mepx'   1.0
#'N2O'   non-reactive
'NO2'      1
'no2'     1.0
'NO'       1
'no'      1.0
'NO3'      1
'no3'     1.0
'N2O5'     1
'n2o5'    1.0
'O1D'      1
'o1d'     1.0
'O3'       1  # unit is in kg/kg
'o3'      1.0
'OH'       1
'oh'      1.0
# 'Br'  'BrCl' 'BrO' 'BrONO2' 'HBr' 'HOBr'
'Cl'       1
'cl'      1.0
'Cl2'      1
'cl2'     1.0
'ClO'      1
'clo'     1.0
# 'Cl2O2' 
'ClONO2'   1
'clno2'   1.0
'HCl'      1
'hcl'     1.0
'HOCl'     1
'hocl'    1.0
# 'OClO' 'CH3Br' 'CH3Cl' 'CH3CCl3' 'CCl4' 'CFCl3'  'CF2Cl2' 'CFC113' 'CFC114' 'CFC115' 'HCFC22'
# 'HCFC141b'  'HCFC142b' 'CF2Br2' 'CF2ClBr' 'CF3Br' 'H2402'
'A3O2'     2   # primary RO2 from C3H8:  CH3CH2CH2OO
'par'     1.0    'xo2'  1.0
'ACTA'     1  # acetic acid
'aacd'    1.0
'ALD2'     1
'ald2'    1.0
'ALK4'     1  # >= C4 alkanes
'par'     4.0
'ATO2'     2  # RO2 from acetone  CH3C(O)CH2O2
'xo2'  1.0    'par'  2.0
'B3O2'     2  #  secondary RO2 from C3H8: CH3CH(OO)CH3
'xo2'  1.0    'par'  2.0
'C2H6'     1
'etha'    1.0
'C3H8'     1
'prpa'    1.0 
'EOH'      1
'etoh'    1.0
'ETO2'     2  # ethylperoxy radical: CH3CH2OO
'meo2' 1.0   'par' 1.0
'ETP'    2    # ethylhydroperoxide: CH3CH2OOH
'MEPX' 1.0   'PAR' 1.0 
'GCO3'   1     # HOCH2C(O)OO  hydroxy peroxyacetyl radical
'c2o3'  1.0 
'GLYX'  2  # glyoxal  CHOCHO
'form' 1.0   'par'  1.0
'GLYC'   2     # HOCH2CHO Glycolaldehyde (Hydroxyacetaldehyde)
'form'  1.0   'par' 2.0
'GP'     1     # HOCH2C(O)OOH   Peroxide from GCO3
'rooh'  1.0   # Higher organic peroxide
'GPAN'   1     # HOCH2C(O)OONO2          Peroxyacylnitrate from GCO3
'panx'  1.0  # C3 and higher peroxyacyl nitrate
'HAC'    1    # HOCH2C(O)CH3  hydroxyacetone
'par'   2.0  
'IALD'   1    # hydroxy carbonyl alkenes from isoprene  HOCH2C(CH3)=CHCHO
'isopx' 1.0
'IAO2'    1  # RO2 from isoprene oxidation products
'isopo2' 1.0
'IAP'     1   # Peroxide from IAO2, HOCH2C(CH3)(OOH)CH(OH)CHO
'rooh'   1.0
'INO2'    7  #  O2NOCH2C(OO)(CH3)CH=CH2: RO2 from ISOP+NO3
'ispd'  0.2  'ntr' 0.8 'xo2' 1.0 'ho2' 0.8 'no2' 0.2 'aldx' 0.8 'par' 2.4
'INPN'   7 # peroxide from INO2  O2NOCH2C(OOH)(CH3)CH=CH2, HO2+INO2=INPN
'ispd'  0.2  'ntr' 0.8 'rooh' 1.0 'h2o2' 0.8 'pna' 0.2 'aldx' 0.8 'par' 2.4
'ISN1'     1   # RO2 from isoprene nitrate
'ntri'    1.0  # second generation isoprene nitrate
'ISNP'     1   # Peroxide from ISN1, HO2 + ISN1 = ISNP
'ntrio2'  1.0  # hydroxyperoxy radicals from second generation isoprene nitrate
'ISOP'     1
'isop'    1.0
'KO2'      2  # RO2 from >C3 ketones
'xo2'     1.0    'par'  1.0
'MACR'     1  # CHOC(CH2)CH3            Methacrolein
'ispd'    1.0 
'MAN2'     2  # RO2 from MACR+NO3
'ho2'  0.925  'xo2' 0.075
'MAO3'   1 # Peroxyacyl from MVK and MACR
'maco3' 1.0  # peroxyacyl radical from methacrolein 
'MAOP'   1  # Peroxide from MAO3
'ispd'  1.0
'MAP'    1   #Peroxyacetic acid, CH3C(O)OOH
'pacd'  1.0
'MCO3'   1     # Peroxyacetyl radical
'c2o3'  1.0
'MEK'    1    #  >C3 ketones
'par'   4.0
'MGLY'   1   # Methylglyoxal  CH3COCHO
'mgly'  1.0
'MRO2'   2   # RO2 from MACR+OH, HOCH2C(OO)(CH3)CHO
'xo2'   0.713  'ho2'  0.503
'MRP'    1 # Peroxide from MRO2, HOCH2C(OOH)(CH3)CHO
'rooh'  1.0
'MVK'    1  # methylvinylketone, CH2=CHC(O)CH3
'ispd'  1.0
'MVN2'   2   # RO2 from MVK+NO3 = MVN2
'ho2'   0.925   'xo2' 0.075
'PAN'    1
'pan'   1.0
'PMN'    1   # CH2=C(CH3)C(O)OONO2     Peroxymethacryloyl nitrate (MPAN)
'open'  1.0
'PO2'    1   #  HOC3H6O2                RO2 from propene
'xo2'   1.0
'PP'     1  #  HOC3H6OOH               Peroxide from PO2
'rooh'  1.0
'PPN'    1   # CH3CH2C(O)OONO2         Peroxypropionylnitrate
'panx'  1.0
'PRN1'   1    #RO2 from propene+NO3
'xo2'   1.0
'PRPE'   2  # C3H6
'ole'   1.0   'par'  1.0
'PRPN'   1  # Peroxide from PRN1
'rooh'  1.0
'R4N1'   2  # RO2 from R4N2, C4,5 alkylnitrates
'rooh'  1.0   'par'  1.0 
'R4N2'   2   # C4,5 alkylnitrates
'ntr'   1.   'par'  2.
'R4O2'  1 # RO2 from ALK4
'xo2'  1.0
'R4P'   1      # Peroxide from R4O2
'rooh' 1.0
'RA3P'  1     # Peroxide from A3O2
'rooh' 1.0
'RB3P'  1    # Peroxide from B3O2
'rooh' 1.0
'RCHO'  1   # >C2 aldehydes
'aldx' 1.0
'RCO3'  1   #  Peroxypropionyl radical, CH3CH2C(O)OO 
'xo2'  1.0
'RCOOH'  1  # >C2 organic acids
'rooh'  1.0
'RIO1'    1 # RO2 from isoprene oxidation products, HOCH2C(OO)(CH3)CH=CHOH
'ispd'   1.0
'RIO2'    1  # HOCH2C(OO)(CH3)CH=CH2	RO2 from isoprene
'isopo2' 1.0
'RIP'     1  # Peroxide from RIO2
'isopx'  1.0
'ROH'     1  # >C2 alcohols
'par'   3.0
'RP'     1  # Peroxide from RCO3
'rooh'  1.0
'VRO2'    1 # RO2 from MVK+OH
'isopo2' 1.0
'VRP'     1  # Peroxide from VRO2
'rooh'   1.0
'ACET'    1
'par'    3.0
'HNO3COND' 1  # condensed HNO3
'hno3'    1.0
# need unit conversion
'SO2'   1
'so2'  1.0
'SO2v'  1
'so2'  1.0
'NH3'   1
'nh3'  1.0
## aerosolin kg/kg to ug/m3
'BCPHILIC'   2         # convert to kg/kg, to ug/std m3  12*0.0409, std air density is 1.225kg/m3
'aecj'      0.0    'numacc'  0.0
# 'aecj'      1.0    'numacc' 27205909.  # assuming mean diameter is 0.3 um (volume= 0.01414x10^-18 m3) and density is 2.6x10^3 kg/m3 or 2.6x10^12 ug/m3.so 1 particle = 0.036x10^-6 ug 
'BCPHOBIC'   2
'aecj'      0.0    'numacc'  0.0
#'aecj'      1.0    'numacc' 27205909.
'OCPHILIC'   2
'apocj'     0.0    'numacc'  0.0
#'apocj'     1.0    'numacc' 27205909.
'OCPHOBIC'   2
'apocj'     0.0    'numacc'  0.0
#'apocj'     1.0    'numacc' 27205909.
'SO4'        2
'aso4j'     1.0    'numacc' 27205909.
'SO4v'       2
'aso4j'     1.0    'numacc' 27205909.
'NH4a'       2
'anh4j'     1.0    'numacc' 27205909.
'NO3an1'     2  # mean radius 0.2695, volume is 0.01026x10^-18 m3, density is 1.8x10^3 kg/m3 or 1.8x10^12 ug/m3, so  1 particle = 0.01847x10^-6 ug
'ano3j'     1.0    'numacc' 54141851.
'NO3an2'     4  # mean radius 2.1
'ano3j'   0.8   'numacc'  5420652.   'ano3k'  0.2  'numcor'  12000.
'NO3an3'     2  # mean radius 7.57
'ano3k'   1.0   'numcor' 6000.
'DU001'      2  ## 0.2-2 um diameter:  assuming mean diameter is 0.3 um (volume= 0.01414x10^-18 m3) and density is 2.6x10^3 kg/m3 or 2.6x10^12 ug/m3.so 1 particle = 0.036x10^-6 ug 
'aothrj'  0.0   'numacc' 0.0
#'aothrj'  1.0   'numacc' 27205909.
'DU002'    4  ## 2-3.6 um
'aothrj'  0.0    'numacc'  0.0  'asoil'  0.0   'numcor'  0.0 
#'aothrj'  0.45    'numacc'  330882.  'asoil'  0.55   'numcor'  50607. 
'DU003'    2  ## 3.6-6um
'asoil'   0.0     'numcor' 0.0
#'asoil'   1.0     'numcor' 11501.
'DU004'    1  ## 6-12um
'asoil'  0.0   'numcor' 0.0
#'asoil'  0.7586   'numcor' 1437.
'DU005'    1  #12-20 um
'asoil'   0.0
## sea salt split into Na and Cl
'SS001'     3    ## 0.06-0.2um  diameter 35.45+22.99, density is 2.2x10^3 kg/m3, mean diameter 0.1um
'anai'    0.0      'acli'    0.0   'numatkn'  0.0
#'anai'    0.39      'acli'    0.61   'numatkn'  736106000.
'SS002'     3    ## 0.2-1 um, mean diameter 0.3 um
'anaj'    0.0      'aclj'    0.0   'numacc' 0.0
#'anaj'    0.39      'aclj'    0.61   'numacc' 27263185.
'SS003'     6    ## 1-3 um   8:2 split    mean diameter 1.6, 2.4
'anaj'    0.0  'aseacat' 0.0    'aclj'  0.0   'aclk'  0.0   'numacc'  0.    'numcor'  0.
#'anaj'    0.312  'aseacat' 0.078    'aclj'  0.488   'aclk'  0.122   'numacc'  169910.    'numcor'  12586.
'SSOO4'    3   # 3-10 um, mean diameter 4um
'aseacat'  0.0     'aclk'   0.0   'numcor'  0.0 
#'aseacat'    0.39     'aclk'   0.61   'numcor'  13592. 
'SS005'     1    ## 10-20um
'ASEACAT'    0.0
EOF

srun --time=30:00 -n $NUMTS geos-rrfs-lbc-para.x

if [ $(stat -Lc %s $TDIR/gfs_bndy.tile7.000.nc) -gt 180000000 ]; then
 mv $TDIR/gfs_bndy.tile7.000.nc $TDIR/geos5_bndy.c793.$year$mon.nc
 ncks -C -O -x -v i_bottom,j_bottom,i_top,j_top,i_right,j_right,i_left,j_left,ps_bottom,ps_top,ps_right,ps_left,\
t_bottom,t_top,t_right,t_left,w_bottom,w_top,w_right,w_left,zh_bottom,zh_top,zh_right,zh_left,sphum_bottom,sphum_top,sphum_right,sphum_left,\
liq_wat_bottom,liq_wat_top,liq_wat_right,liq_wat_left,o3mr_bottom,o3mr_top,o3mr_right,o3mr_left,ice_wat_bottom,ice_wat_top,ice_wat_right,ice_wat_left,\
rainwat_bottom,rainwat_top,rainwat_right,rainwat_left,snowwat_bottom,snowwat_top,snowwat_right,snowwat_left,graupel_bottom,graupel_top,graupel_right,graupel_left,\
i_w_bottom,j_w_bottom,i_w_top,j_w_top,i_w_right,j_w_right,i_w_left,j_w_left,u_w_bottom,u_w_top,u_w_right,u_w_left,v_w_bottom,v_w_top,v_w_right,v_w_left,\
i_s_bottom,j_s_bottom,i_s_top,j_s_top,i_s_right,j_s_right,i_s_left,j_s_left,u_s_bottom,u_s_top,u_s_right,u_s_left,v_s_bottom,v_s_top,v_s_right,v_s_left\
 $TDIR/geos5_bndy.c793.$year$mon.nc /scratch2/NCEPDEV/naqfc/RRFS_CMAQ_NA13km/LBCS/RRFS_NA13km_GEOS5_v2/geos5_bndy_v2.c793.$year$amon.nc
 
# /scratch2/NCEPDEV/naqfc/RRFS_CMAQ_NA13km/LBCS/RRFS_NA13km_GEOS5_v1/geos5_bndy_v1.c793.$year$mon.nc

else
 echo "GEOS5 LBC failed"
 exit 1
fi   

let mon=$mon+1
done
