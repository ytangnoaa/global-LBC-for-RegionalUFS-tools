      program am4_bnd
!-------------------------------------------------------------------------------------      
!   interpolate AM4 Concentration to be lateral boundary condition for regional 
!   air quality model,  also output a layer result for checking purpose
!
!
!   Author: Youhua Tang
!   Revisions: AM4 monthly data to RRFS-CMAQ
!-------------------------------------------------------------------------------------
!
!                             nhalo_model=3
!
!                    |----------- nxp-1 -----------| <-- east/west compute points
!                 |---------- north BC data ----------|
!                 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!       ---       ooo           ---j=1---           ooo     ---         ---
!        |        ooo                               ooo      |           |
!        |        ooo                              |ooo      |           |
!                 ooo                        i=1-->|ooo
!   west BC data  ooo|                             |ooo east BC data    nyp-1 <-- north/south compute points
!                 ooo|<--i=isd-nhalo_model          ooo
!        |        ooo|                              ooo      |           |
!        |        ooo                               ooo      |           |
!       ---       ooo    ---j=jsd-nhalo_model---    ooo     ---         ---
!                 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                 |---------- south BC data ----------|
      
      use netcdf 

      parameter(maxfile=300,nspecies=200, ngas1=32, naerosol=15, &
       ngeos=ngas1+naerosol)

      real sfact(ngeos,nspecies),val(nspecies),  &
         checkfact(ngeos,nspecies)

      double precision,allocatable  :: glon(:), glat(:)
      real,allocatable  :: pmid_am4(:,:,:),pface_am4(:,:,:),tmpa(:),vgeos(:,:,:), press(:,:,:),  &
       worka(:),workb(:),workc(:), work(:), work1(:), work2(:), work3(:), xlat(:,:), xlon(:,:), &
       bndx(:,:,:,:,:), bndy(:,:,:,:,:),bndcoordx(:,:,:,:) ,bndcoordy(:,:,:,:),  &
       checkcoord(:,:,:),checksp(:,:,:), ps_x(:,:,:),ps_y(:,:,:), pface_tmp(:), &
       pmid_x(:,:,:,:),pmid_y(:,:,:,:),ps_am4(:,:),tmpbndx(:,:,:),tmpbndy(:,:,:)
      
      real bk_numatkn(200),bk_numacc(200),bk_numcor(200), &  ! background aerosol number concentration
         pk_am4(200),bk_am4(200),ak_gfs(200),bk_gfs(200)
	 
      character bndname(nspecies)*16,geosname(ngeos)*8,ctmp*16,  &
       echar(nspecies)*16,am4file*200,checkname(nspecies)*16,     &
       aline*200,gdatatype*4,modelname*4,gtype*16,arank*2, topofile*200, &
       lbcfile*200,gas1name(ngas1)*8, aeroname(naerosol)*8
      
      integer netindex(ngeos),checklayer,modate(maxfile),         &
       mosecs(maxfile),julian,ismotime(maxfile),iemotime(maxfile),  &
       idate(7),tlmeta,iret,istart2d(3),icount2d(3),istart3d(4),icount3d(4)
      logical ingeos,lflag,extrameta,indexfind(nspecies)
      integer monthday(12),dimids(3),ifind(1)
      data monthday/31,30,31,30,31,30,31,31,30,31,30,31/

      data gas1name/'O3','C10H16','C2H4','C2H5OH','C2H6','C3H6','C3H8','C4H10','CH2O','CH3CHO', &
       'CH3COCH3','CH3OH','CH3OOH','CO','DMS','H2O2','HNO3','HO2','ISOP','N2O5', &
       'N2O','NH3','NH4NO3','NO2','NO','NOy','PAN','SO2','SO4','Cly','Bry','H2'/

      data aeroname/'bcphil','bcphob','omphil','omphob','dust1','dust2','dust3','dust4','dust5','SOA', &
       'ssalt1','ssalt2','ssalt3','ssalt4','ssalt5'/
       
      data indexfind/nspecies*.false./
      
      integer  begyear,begdate,begtime,dtstep,numts,tstepdiff      
      namelist /control/iprint,bndname,o3cap,am4file,imonth, &	  !  input file preffix and suffix
       lbcfile,topofile,bk_numatkn,bk_numacc,bk_numcor,      &
       pk_am4,bk_am4,ak_gfs,bk_gfs,inblend  ! inside blending layers
            
      call aq_blank(16*nspecies,bndname)
      call aq_blank(16*nspecies,checkname)

      sfact(1:ngeos,1:nspecies)=0.
      checkfact(1:ngeos,1:nspecies)=0.

      geosname(1:ngas1)=gas1name(1:ngas1)
      geosname(ngas1+1:ngeos)=aeroname(1:naerosol)
      
      inblend=0.
      bk_numatkn=-1.; bk_numacc=-1; bk_numcor=-1
      pk_am4=-1.; bk_am4=-1.; ak_gfs=-1. ; bk_gfs=-1.
      
! read converting information

      open(7,file='am4-rrfs-lbc.ini')
      read(7,control)
      
      ifind=findloc(pk_am4,-1.)
      kmax1_am4=ifind(1)-1
      ifind=findloc(bk_am4,-1.)
      if(kmax1_am4.ne.ifind(1)-1) then
       print*,'kmax1_am4 inconsistent ',kmax1_am4,ifind(1)-1
       stop
      endif 
      
      ifind=findloc(ak_gfs,-1.)
      kmax1_gfs=ifind(1)-1

      ifind=findloc(bk_gfs,-1.)
      if(kmax1_gfs.ne.ifind(1)-1) then
       print*,'kmax1_gfs inconsistent ',kmax1_gfs,ifind(1)-1
       stop
      endif 
      print*,'kmax1_am4,kmax1_gfs=',kmax1_am4,kmax1_gfs
      
      ifind=findloc(bk_numatkn,-1.)
      if(kmax1_gfs.ne.ifind(1)-1) then
       print*,'bk_numatkn kmax1_gfs inconsistent ',kmax1_gfs,ifind(1)-1
       stop
      endif
      ifind=findloc(bk_numacc,-1.)
      if(kmax1_gfs.ne.ifind(1)-1) then
       print*,'bk_numacc kmax1_gfs inconsistent ',kmax1_gfs,ifind(1)-1
       stop
      endif
      ifind=findloc(bk_numcor,-1.)
      if(kmax1_gfs.ne.ifind(1)-1) then
       print*,'bk_numcor kmax1_gfs inconsistent ',kmax1_gfs,ifind(1)-1
       stop
      endif
          
      call aq_find(nspecies,' ',bndname,lpsec,iflag)   ! BND species
      noutbnd=lpsec-1
      call aq_find(nspecies,' ',checkname,lpsec,iflag)   ! BND species
      ncheck=lpsec-1

      call aq_locate(7,'Species converting Factor',iflag)      
      if(iflag.ne.0) then
        print*,'can not find Species converting Factor'
       stop
      endif
      do while(.true.)
       call aq_readhd(7)
       read(7,*,end=98,err=99) ctmp,num
       call aq_find(ngeos,ctmp,geosname,lpsec1,iflag)
       if(iflag.eq.0) then
        read(7,*)(echar(i),val(i),i=1,num)
	do i=1,num
	 call aq_find(noutbnd,echar(i),bndname,lpsec2,iflag)
	 if(iflag.eq.0) then
	  sfact(lpsec1,lpsec2)=val(i)
	  indexfind(lpsec2)=.true.
	 endif 
	 call aq_find(ncheck,echar(i),checkname,lpsec2,iflag)
	 if(iflag.eq.0) checkfact(lpsec1,lpsec2)=val(i)   
	end do
       endif 
      print*,' Converting factor for ',geosname(lpsec1),' is ', &
      (sfact(lpsec1,lm),lm=1,noutbnd) 
 99   continue
      enddo 
 98   close(7)

      call check(nf90_open(trim(topofile),nf90_nowrite, ncid))
      call check(nf90_inq_dimid(ncid,'lon',iddim_lon))
      call check(nf90_inquire_dimension(ncid,iddim_lon,len=imax))
      call check(nf90_inq_dimid(ncid,'lat',iddim_lat))
      call check(nf90_inquire_dimension(ncid,iddim_lat,len=jmax1))

      allocate(xlat(imax,jmax1))
      allocate(xlon(imax,jmax1))      
      
      call check(nf90_inq_varid(ncid,'geolon',idvar_geolon))
      call check(nf90_get_var(ncid,idvar_geolon,xlon))
      call check(nf90_inq_varid(ncid,'geolat',idvar_geolat))
      call check(nf90_get_var(ncid,idvar_geolat,xlat))
      
      do i=1,imax
       do j=1,jmax1
       if(xlon(i,j).lt.0) xlon(i,j)=xlon(i,j)+360
       enddo
      enddo 
      call check(nf90_close(ncid))
      print*,'finish reading topofile' 

! open LBC file for rewrite

      call check(nf90_open(trim(lbcfile),nf90_write, ncid))
      call check(nf90_inq_dimid(ncid,'lon',iddim_lon))
      call check(nf90_inquire_dimension(ncid,iddim_lon,len=nlon))
      call check(nf90_inq_dimid(ncid,'lat',iddim_lat))
      call check(nf90_inquire_dimension(ncid,iddim_lat,len=nlat))
      call check(nf90_inq_dimid(ncid,'halo',iddim_halo))
      call check(nf90_inquire_dimension(ncid,iddim_halo,len=nhalo))
      
      nhalo_outside=nhalo-inblend ! outside halo layers
      jmax=jmax1-nhalo_outside*2
      if(nlon.ne.imax.or.nlat.ne.jmax) then
        print*,'dimension mismatch ',nlon,imax,nlat,jmax
	stop
      endif
      
! read ps
      call check(nf90_inq_dimid(ncid,'lev',iddim_lev))
      call check(nf90_inquire_dimension(ncid,iddim_lev,len=kmax))
      call check(nf90_inq_dimid(ncid,'levp',iddim_levp))
      call check(nf90_inquire_dimension(ncid,iddim_levp,len=kmax1))
      if(kmax.ne.kmax1_gfs) then
        print*,'profile for background aerosol numbers need adjustment'
	stop
      endif
      if(kmax1.ne.(kmax+1)) then
       print*,'kmax,kmax1=',kmax,kmax1
       stop
      endif 
      
      allocate(pmid_x(imax,nhalo,kmax,2),pmid_y(nhalo,jmax,kmax,2), &
         ps_x(imax,nhalo,2),ps_y(nhalo,jmax,2),pface_tmp(kmax1) )
      allocate(bndcoordx(imax,nhalo,2,2+kmax),bndcoordy(nhalo,jmax,2,2+kmax))
      allocate(bndx(imax,nhalo,kmax,2,noutbnd),bndy(nhalo,jmax,kmax,2,noutbnd))
      bndx=0.
      bndy=0.
      do jj=1,2
       do i=1,imax
        do j=1,nhalo
	 do L=1,noutbnd        
	 if(index(bndname(L),'numatkn').gt.0) bndx(i,j,1:kmax,jj,L)=bk_numatkn(1:kmax)
         if(index(bndname(L),'numacc').gt.0) bndx(i,j,1:kmax,jj,L)=bk_numacc(1:kmax)
         if(index(bndname(L),'numcor').gt.0) then
	   L_numcor=L
	   bndx(i,j,1:kmax,jj,L)=bk_numcor(1:kmax)
	 endif  
         enddo
	enddo
       enddo
       do i=1,nhalo
        do j=1,jmax
	 do L=1,noutbnd
	 if(index(bndname(L),'numatkn').gt.0) bndy(i,j,1:kmax,jj,L)=bk_numatkn(1:kmax)
         if(index(bndname(L),'numacc').gt.0) bndy(i,j,1:kmax,jj,L)=bk_numacc(1:kmax)
         if(index(bndname(L),'numcor').gt.0) bndy(i,j,1:kmax,jj,L)=bk_numcor(1:kmax)
	 enddo
        enddo
       enddo
      enddo		
      print*,'1 bndx(numcor) min,max=',minval(bndx(:,:,:,:,L_numcor)),maxval(bndx(:,:,:,:,L_numcor))
      print*,'1 bndy(numcor) min,max=',minval(bndy(:,:,:,:,L_numcor)),maxval(bndy(:,:,:,:,L_numcor))
      
      allocate(tmpbndx(imax,nhalo,kmax),tmpbndy(nhalo,jmax,kmax))

!! calculate pressure	  
!---process ps_left, ps_right
      print*,'read ps_left, ps_right'
      call check(nf90_inq_varid(ncid,'ps_right',idvar_ps_right))
      call check(nf90_get_var(ncid,idvar_ps_right,ps_y(:,:,2)))

      call check(nf90_inq_varid(ncid,'ps_left',idvar_ps_left))
      call check(nf90_get_var(ncid,idvar_ps_left,ps_y(:,:,1)))

      do i=1,nhalo
        do j=1,jmax
	  do m=1,2
	   pface_tmp(1:kmax1)=bk_gfs(1:kmax1)*ps_y(i,j,m)+ak_gfs(1:kmax1)
           do k=1,kmax       
            delp=ak_gfs(k+1)+bk_gfs(k+1)*ps_y(i,j,m)-(ak_gfs(k)+bk_gfs(k)*ps_y(i,j,m))
	    if(abs(delp).le.1e-33) then
	     pmid_y(i,j,k,m)=0.
	    else 
	     pmid_y(i,j,k,m)=delp/(log(pface_tmp(k+1))-log(pface_tmp(k)))
	    endif 
	   enddo
	  enddo
	enddo
      enddo

      print*,'read ps_bottom, ps_top'
      call check(nf90_inq_varid(ncid,'ps_bottom',idvar_ps_bottom))
      call check(nf90_get_var(ncid,idvar_ps_bottom,ps_x(:,:,1)))
      call check(nf90_inq_varid(ncid,'ps_top',idvar_ps_top))
      call check(nf90_get_var(ncid,idvar_ps_top,ps_x(:,:,2)))
      
      do i=1,imax
        do j=1,nhalo
	  do m=1,2
	   pface_tmp(1:kmax1)=bk_gfs(1:kmax1)*ps_x(i,j,m)+ak_gfs(1:kmax1)
           do k=1,kmax       
            delp=ak_gfs(k+1)+bk_gfs(k+1)*ps_x(i,j,m)-(ak_gfs(k)+bk_gfs(k)*ps_x(i,j,m))
	    if(abs(delp).le.1e-33) then
	       pmid_x(i,j,k,m)=0.
	    else   
	       pmid_x(i,j,k,m)=delp/(log(pface_tmp(k+1))-log(pface_tmp(k)))
	    endif   
	   enddo
	  enddo
	enddo
      enddo
      
      print*,'imax,jmax,kmax=',imax,jmax,kmax
      print*,'pmid_x bottom min/max=',minval(pmid_x(:,1,:,1)),maxval(pmid_x(:,1,:,1))

! -open AM4

      call check(nf90_open(trim(am4file),nf90_nowrite, id_file))
      call check(nf90_inq_dimid(id_file,'lat',id_nlat))
      call check(nf90_inquire_dimension(id_file,id_nlat,len=nlatgeos))
      call check(nf90_inq_dimid(id_file,'lon',id_nlon))
      call check(nf90_inquire_dimension(id_file,id_nlon,len=nlongeos))
      call check(nf90_inq_dimid(id_file,'pfull',id_nlev))
      call check(nf90_inquire_dimension(id_file,id_nlev,len=nlevgeos))
      
       igeos=nlongeos
       jgeos=nlatgeos
       kgeos=nlevgeos
       if(kmax1_am4.ne.kgeos+1) then
        print*,'inconsistent AM4 cooridate ',kgeos,kmax1_am4
	stop
       endif
       	
       istart2d=[1,1,imonth]; icount2d=[igeos,jgeos,1]
       istart3d=[1,1,1,imonth]; icount3d=[igeos,jgeos,kgeos,1]
       
       allocate(glon(igeos),glat(jgeos),tmpa(kgeos),ps_am4(igeos,jgeos),pmid_am4(igeos,jgeos,kgeos))
       allocate(vgeos(igeos,jgeos,kgeos))
       print*,'igeos,jgeos,kgeos=',igeos,jgeos,kgeos

!--- AM4 variables             
        call check(nf90_inq_varid(id_file,'lon',id_lon))
        call check(nf90_get_var(id_file,id_lon,glon))
	call check(nf90_inq_varid(id_file,'lat',id_lat))
        call check(nf90_get_var(id_file,id_lat,glat))
	call check(nf90_inq_varid(id_file,'ps',id_ps))
	print*,'read ps_am4'
	call check(nf90_get_var(id_file,id_ps,ps_am4,start=istart2d,count=icount2d))  
	print*,'finish reading ps_am4'
	
	do i=1,igeos
	 do j=1,jgeos
	  pface_tmp(1:kmax1_am4)=bk_am4(1:kmax1_am4)*ps_am4(i,j)+pk_am4(1:kmax1_am4)	  
	  do k=1,kgeos
	   delp=pk_am4(k+1)+bk_am4(k+1)*ps_am4(i,j)-(pk_am4(k)+bk_am4(k)*ps_am4(i,j))
	   pmid_am4(i,j,k)=delp/(log(pface_tmp(k+1))-log(pface_tmp(k)))
	  enddo
	 enddo
	enddo  

	glonint=sngl((glon(igeos)-glon(1))/(igeos-1))
	glatint=sngl((glat(jgeos)-glon(1))/(jgeos-1))
	
!	ishift=0
!	if((glon(1)+180).lt.2) then	  
!	  do i=1,igeos
!	   if(abs(glon(i)).lt.glonint/4) exit  ! shift to start from lon=0.0 for North America 
!	  enddo	  
!	  if(i.le.igeos) then
!	   ishift=i-1
!	   glon(i)=0   ! fix the precision issue for lon
!	  endif 
!	  where(glon.lt.0) glon=glon+360
!	  if(i.gt.1) then
!	   glon=cshift(glon,shift=ishift,dim=1)
!	   zgeos=cshift(zgeos,shift=ishift,dim=1)
!	   airgeos=cshift(airgeos,shift=ishift,dim=1)
!	   print*,ishift,'after shift, glon=',glon
!	  endif 
!	endif

!---calculating lateral boundary horizontal index in AM4 coordinate	

! --- top and bottom

       do i=1,imax
        ix=i
	do j=1,nhalo
 	 do n=1,2
	  
	  select case(n)
	   case(1)	 
	    jy=j      ! bottom
	   case(2) 
	    jy=jmax1-nhalo+j ! top
	  end select    
         
	  do i2=1,igeos-1
           if(xlon(ix,jy).ge.glon(i2).and.xlon(ix,jy).le.glon(i2+1)) then
	    bndcoordx(i,j,n,1)=i2+(xlon(ix,jy)-glon(i2))/     &   ! i in AM4 coordiate
      	    (glon(i2+1)-glon(i2))  
	    exit
	   endif
	  enddo
	 
          do j2=1,jgeos
           if(xlat(ix,jy).ge.glat(j2).and.xlat(ix,jy).le.glat(j2+1)) then
	    bndcoordx(i,j,n,2)=j2+(xlat(ix,jy)-glat(j2))/     &   ! j in AM4 coordiate
      	    (glat(j2+1)-glat(j2))  
	    exit
	   endif
	  enddo
	  
!---vertical index for top/bottom     
	  x=bndcoordx(i,j,n,1) 	 
	  y=bndcoordx(i,j,n,2)
	  xratio=x-int(x)
	  yratio=y-int(y)
	 
	  do kp=1,kgeos
      	   tmpa(kp)=(1-yratio)*(pmid_am4(int(x),int(y),kp)*      &    ! horizontally interpolate pressure
     	    (1-xratio)+pmid_am4(int(x)+1,int(y),kp)*xratio)+     &
     	    yratio*(pmid_am4(int(x),int(y)+1,kp)*(1-xratio)+     &
     	    pmid_am4(int(x)+1,int(y)+1,kp)*xratio)
          enddo

          do k=1,kmax
	   if(pmid_x(i,j,k,n).ge.tmpa(kgeos)) then  ! surface
	     bndcoordx(i,j,n,k+2)=real(kgeos)  ! k index start from 3
	   else  
 	     do kp=2,kgeos
 	      if(pmid_x(i,j,k,n).le.tmpa(kp).and.pmid_x(i,j,k,n).ge.tmpa(kp-1)) then  ! both use top-down coordinate
                bndcoordx(i,j,n,k+2)=kp-1+(pmid_x(i,j,k,n)-tmpa(kp-1))/  &
     	          (tmpa(kp)-tmpa(kp-1))
                 exit
 	      endif  
	     enddo
           endif
           if(pmid_x(i,j,k,n).le.tmpa(1)) bndcoordx(i,j,n,k+2)=1.   ! model top  
          enddo
        
	enddo
       enddo 	
      enddo

! --- left and right
       
      do j=1,jmax
       jy=j+nhalo_outside
           
	do i=1,nhalo
	 do n=1,2
	  if (n.eq.1) then
	   ix=i      !   left
	  else if (n.eq.2) then
	   ix=imax-nhalo+i ! right
	  endif    
         
	 do i2=1,igeos-1
          if(xlon(ix,jy).ge.glon(i2).and.xlon(ix,jy).le.glon(i2+1)) then
	   bndcoordy(i,j,n,1)=i2+(xlon(ix,jy)-glon(i2))/     &   ! i in AM4 coordiate
      	    (glon(i2+1)-glon(i2))  
	   exit
	  endif
	 enddo
	 
         do j2=1,jgeos-1
 	  if(xlat(ix,jy).ge.glat(j2).and.xlat(ix,jy).le.glat(j2+1)) then
	   bndcoordy(i,j,n,2)=j2+(xlat(ix,jy)-glat(j2))/     &   ! j in AM4 coordiate
     	    (glat(j2+1)-glat(j2))  
	   exit
	  endif
	 enddo

!---vertical index for left/right   
	  x=bndcoordy(i,j,n,1) 	 
	  y=bndcoordy(i,j,n,2)
	  xratio=x-int(x)
	  yratio=y-int(y)
	 
	  do kp=1,kgeos
      	   tmpa(kp)=(1-yratio)*(pmid_am4(int(x),int(y),kp)*      &    ! horizontally interpolate height
     	    (1-xratio)+pmid_am4(int(x)+1,int(y),kp)*xratio)+     &
     	    yratio*(pmid_am4(int(x),int(y)+1,kp)*(1-xratio)+     &
     	    pmid_am4(int(x)+1,int(y)+1,kp)*xratio)
          enddo

          do k=1,kmax
	   if(pmid_y(i,j,k,n).ge.tmpa(kgeos)) then ! surface
	     bndcoordy(i,j,n,k+2)=real(kgeos)  ! k index start from 3
	   else  
 	     do kp=2,kgeos
 	      if(pmid_y(i,j,k,n).le.tmpa(kp).and.pmid_y(i,j,k,n).ge.tmpa(kp-1)) then  ! both use top-down coordinate
                bndcoordy(i,j,n,k+2)=kp-1+(pmid_y(i,j,k,n)-tmpa(kp-1))/  &
     	          (tmpa(kp)-tmpa(kp-1))
                 exit
 	      endif  
	     enddo
           endif
           if(pmid_y(i,j,k,n).le.tmpa(1)) bndcoordy(i,j,n,k+2)=1.     
          enddo

         enddo
        enddo
       enddo
		  
	  
       if(iprint.eq.1) then
!         open(27,file='dust2.bin',form='unformatted',access='direct',recl=igeos*jgeos*4)
         open(27,file='o3-tmp.bin',form='unformatted',access='direct',recl=igeos*jgeos*4)
	 open(28,file='pmid-tmp.bin',form='unformatted',access='direct',recl=imax*kmax1*4)
	 write(28,rec=1)pmid_x(1:imax,1,1:kmax,1)
	 write(28,rec=2)pmid_x(1:imax,1,1:kmax,2)
	 print*,'pmid_x2 top =',pmid_x(1:imax,1,1,2)
	 print*,'pmid_x bottom min/max after print=',minval(pmid_x(:,1,:,1)),maxval(pmid_x(:,1,:,1))
         print*,'pmid_x top min/max after print=',minval(pmid_x(:,1,:,2)),maxval(pmid_x(:,1,:,2))
	 close(28) 
       endif 
    
    if(m.eq.1) then
       nowspc=ngas1
    else if (m.eq.2) then
       nowspc=ngas2
    else if (m.eq.3) then
       nowspc=naerosol
    endif 
       	
  ! begin species interpolation                
   do L=1,ngeos
    
    ctmp=geosname(L)
    call check(nf90_inq_varid(id_file,trim(ctmp),id_var1))
    call check(nf90_get_var(id_file,id_var1,vgeos))
!    if(ishift.gt.0) vgeos=cshift(vgeos,shift=ishift,dim=1)
     
    if(L.le.ngas1) then
       vgeos(:,:,:)=vgeos(:,:,:)*1e6   ! gaseous: mole/mole to ppmV
    else
       vgeos(:,:,:)=vgeos(:,:,:)*1e9   ! aerosol kg/kg - > ug/kg
    endif 
     
    vgeos(:,:,:)=amax1(0.,vgeos(:,:,:)) ! remove negative value	
    if(ctmp.eq.'O3'.and.o3cap.gt.0)  vgeos(:,:,:)=amin1(o3cap,vgeos(:,:,:)) ! capped o3 in ppmv	
 
    if(iprint.eq.1) then
      if(ctmp.eq.'O3') write(27,rec=1) vgeos(:,:,kgeos)
      if(ctmp.eq.'CO') write(27,rec=2) vgeos(:,:,kgeos)
    endif
	
     do i=1,imax       ! for top/bottom boundary conditions
      do j=1,nhalo
       do n=1,2
        x=bndcoordx(i,j,n,1)
	y=bndcoordx(i,j,n,2)
	xratio=x-int(x)
	yratio=y-int(y)

	tmpa(1:kgeos)=(1-yratio)*(vgeos(int(x),int(y),     & ! horizontally interpolate values
     	 1:kgeos)*(1-xratio)+vgeos(int(x)+1,int(y),       &
     	 1:kgeos)*xratio)+yratio*(vgeos(int(x),int(y)+1,  &
     	 1:kgeos)*(1-xratio)+vgeos(int(x)+1,int(y)+1,     &
     	 1:kgeos)*xratio)
	
	 do k=1,kmax
	  z=bndcoordx(i,j,n,k+2)
	  if(z.ge.kgeos) then
	   tmpvalue=tmpa(int(z))
	  else 
	   zratio=z-int(z)
	   tmpvalue=(1-zratio)*tmpa(int(z))+zratio*tmpa(int(z)+1)     ! vertically interpolate values
	  endif 
	  do L2=1,noutbnd
	   bndx(i,j,k,n,L2)=bndx(i,j,k,n,L2)+amax1(tmpvalue,0.)*sfact(L,L2)
	  enddo
	 enddo
	 
	enddo
       enddo
      enddo

     do i=1,nhalo       ! for left/right boundary conditions
      do j=1,jmax
       do n=1,2
        x=bndcoordy(i,j,n,1)
	y=bndcoordy(i,j,n,2)
	xratio=x-int(x)
	yratio=y-int(y)

	tmpa(1:kgeos)=(1-yratio)*(vgeos(int(x),int(y),     & ! horizontally interpolate values
     	 1:kgeos)*(1-xratio)+vgeos(int(x)+1,int(y),       &
     	 1:kgeos)*xratio)+yratio*(vgeos(int(x),int(y)+1,  &
     	 1:kgeos)*(1-xratio)+vgeos(int(x)+1,int(y)+1,     &
     	 1:kgeos)*xratio)
	
	 do k=1,kmax
	  z=bndcoordy(i,j,n,k+2)
	  if(z.ge.kgeos) then
	   tmpvalue=tmpa(int(z))
	  else 
	   zratio=z-int(z)
	   tmpvalue=(1-zratio)*tmpa(int(z))+zratio*tmpa(int(z)+1)     ! vertically interpolate values
	  endif 
	  do L2=1,noutbnd
	   bndy(i,j,k,n,L2)=bndy(i,j,k,n,L2)+amax1(tmpvalue,0.)*sfact(L,L2)
	  enddo
	 enddo
	 
	enddo
       enddo
      enddo
      
    enddo  ! end  geos species loop for one file

      print*,'2 bndx(numcor) min,max=',minval(bndx(:,:,:,:,L_numcor)),maxval(bndx(:,:,:,:,L_numcor))
      print*,'2 bndy(numcor) min,max=',minval(bndy(:,:,:,:,L_numcor)),maxval(bndy(:,:,:,:,L_numcor))

  
! begin output
     fillval=-9e33
     print*,'start write' 
      do L=1,noutbnd       ! check if global model supplies all species, otherwise do not overwrite the existing aerosol

        do m=1,2

!! -top/bottom	  
	  if(sum(bndx(1:imax,1:nhalo,1:kmax,m,L)).le.1e-20 ) then
	    print*, 'm=,', m,' X skip for ', bndname(L)
	    cycle  
	  endif
	  if (m.eq.1) then
	   aline=trim(bndname(L))//'_bottom'
	  else
	   aline=trim(bndname(L))//'_top'
	  endif 
          if(nf90_inq_varid(ncid,trim(aline),idvar_tmp).ne.nf90_noerr) then ! if not exit, create
	    dimids=(/iddim_lon,iddim_halo,iddim_lev/)
	    call check(nf90_redef(ncid))
	    print*,'add variable ',trim(aline)
	    call check(nf90_def_var(ncid,trim(aline),nf90_real,dimids,idvar_tmp))
	    call check(nf90_put_att(ncid,idvar_tmp,'_FillValue',fillval))
	    call check(nf90_enddef(ncid))
	    call check(nf90_inq_varid(ncid,trim(aline),idvar_tmp)) ! check exist
	  endif
	  
	  print*,'write ',trim(aline),idvar_tmp,minval(bndx(1:imax,1:nhalo,1:kmax,m,L)),maxval(bndx(1:imax,1:nhalo,1:kmax,m,L))
	  
	  tmpbndx(1:imax,1:nhalo,1:kmax)=bndx(1:imax,1:nhalo,1:kmax,m,L)
	  call check(nf90_put_var(ncid,idvar_tmp,tmpbndx))
	  
!! left/right
	  if(sum(bndy(1:nhalo,1:jmax,1:kmax,m,L)).le.1e-20 ) then
	    print*, 'm=,',m, ' Y skip for ', bndname(L)
	    cycle  
	  endif
	  if (m.eq.1) then
	   aline=trim(bndname(L))//'_left'
	  else
	   aline=trim(bndname(L))//'_right'
	  endif 

          if(nf90_inq_varid(ncid,trim(aline),idvar_tmp).ne.nf90_noerr) then ! if not exit, create
 	    call check(nf90_redef(ncid))
	    dimids=(/iddim_halo,iddim_lat,iddim_lev/)
	    call check(nf90_def_var(ncid,trim(aline),nf90_real,dimids,idvar_tmp))
	    call check(nf90_put_att(ncid,idvar_tmp,'_FillValue',fillval))
	    call check(nf90_enddef(ncid))
	    call check(nf90_inq_varid(ncid,trim(aline),idvar_tmp)) ! check exist
	  endif
          
	  print*,'write ',trim(aline),idvar_tmp,minval(bndy(1:nhalo,1:jmax,1:kmax,m,L)),maxval(bndy(1:nhalo,1:jmax,1:kmax,m,L))

	  tmpbndy(1:nhalo,1:jmax,1:kmax)=bndy(1:nhalo,1:jmax,1:kmax,m,L)
	  call check(nf90_put_var(ncid,idvar_tmp,tmpbndy))         
	 enddo
	enddo
      call check(nf90_close(ncid))


contains
    subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
    end subroutine check
      end


      subroutine aq_blank(ntot,y)

      character*1 y(ntot)
      do i=1,ntot
      y(i)=' '
      enddo
      return
      end

      subroutine aq_locate(iunit,char,iflag)
!***********************************************************************
      character*(*) char
      character*80 dum1
      nchar=len(char)
      iflag=0
      do iter=1,10000
      read(iunit,'(a)',end=98) dum1(1:nchar)
!      print*,'dum1= ',dum1(1:nchar)
      if(dum1(1:nchar).eq.char) return
      enddo
98    iflag=1
!      print*,'dum1= ',dum1(1:nchar)
      return
      end

!**********************************************************************
      subroutine aq_find(num,cdum1,sname,lpsec,iflag)
!***********************************************************************
      dimension sname(1)
      character*(*) sname,cdum1
      iflag=0
      do 15 l=1,num
      if(cdum1.ne.sname(l)) go to 15
      lpsec=l
      return
15    continue
      iflag=1
      return
      end

      subroutine aq_readhd(iunit)
!***********************************************************************
      character*1 char(80)
       do iter=1,1000
	read(iunit,100,end=8,err=8)  (char(i),i=1,80)
	if(char(1).eq.'$') then
	  write(6,200) iunit,(char(i),i=1,80)
        else if(char(1).eq.'#') then
	else
	  backspace iunit
	  return
        endif
       end do
 8    backspace iunit
100   format(80a1)
200   format(2x,'iunit=',i3,2x,80a1)
      end
