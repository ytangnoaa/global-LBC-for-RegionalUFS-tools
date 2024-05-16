      program hcmaq_bnd
!-------------------------------------------------------------------------------------      
!   interpolate Hemispheric CMAQ Concentration to be lateral boundary condition for regional 
!   air quality model,  also output a layer result for checking purpose
!
!
!   Author: Youhua Tang, Jan. 2023
!   Revisions: Hemispheric CMAQ monthly data to RRFS-CMAQ
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
      use m3utilio
      
      parameter(maxfile=300,nspecies=200)

      real val(nspecies)

      double precision,allocatable  :: glon(:), glat(:)
      real,allocatable  :: p_hcmaq(:,:,:),dens_hcmaq(:,:,:),vhcmaq(:,:,:), press(:,:,:),  &
       worka(:),workb(:),workc(:), work(:), work1(:), work2(:), work3(:), xlat(:,:), xlon(:,:), &
       bndx(:,:,:,:,:), bndy(:,:,:,:,:),bndcoordx(:,:,:,:) ,bndcoordy(:,:,:,:),tmpa(:),  &
       checkcoord(:,:,:),checksp(:,:,:), ps_x(:,:,:),ps_y(:,:,:), pface_tmp(:), &
       pmid_x(:,:,:,:),pmid_y(:,:,:,:),tmpbndx(:,:,:),tmpbndy(:,:,:)
      
      real bk_numatkn(65),bk_numacc(65),bk_numcor(65), &  ! background aerosol number concentration
         ak_gfs(200),bk_gfs(200)
	 
      character bndname(nspecies)*16,ctmp*16,  &
       echar(nspecies)*16,checkname(nspecies)*16,     &
       aline*200,gdatatype*4,modelname*4,gtype*16,arank*2, topofile*200, &
       lbcfile*200
      
      integer checklayer,modate(maxfile),         &
       mosecs(maxfile),ismotime(maxfile),iemotime(maxfile),  &
       idate(7),tlmeta,iret,istart2d(3),icount2d(3),istart3d(4),icount3d(4)
      logical lflag,extrameta,indexfind(nspecies)
      integer monthday(12),dimids(3),ifind(1)
      data monthday/31,30,31,30,31,30,31,31,30,31,30,31/
       
      data indexfind/nspecies*.false./
      
      integer  begyear,begdate,begtime,dtstep,numts,tstepdiff      
      namelist /control/iprint,bndname,o3cap, &	  !  input file preffix and suffix
       lbcfile,topofile,bk_numatkn,bk_numacc,bk_numcor,      &
       ak_gfs,bk_gfs
            
      call aq_blank(16*nspecies,bndname)
      call aq_blank(16*nspecies,checkname)


      ak_gfs=-1. ; bk_gfs=-1.
      
! read converting information

      open(7,file='hcmaq-rrfs-lbc.ini')
      read(7,control)
      
      ifind=findloc(ak_gfs,-1.)
      kmax1_gfs=ifind(1)-1
      if(kmax1_gfs.ne.66) then
        print*,'kmax1_gfs not equal to 66',kmax1_gfs
	stop
      endif	
      ifind=findloc(bk_gfs,-1.)
      if(kmax1_gfs.ne.ifind(1)-1) then
       print*,'kmax1_gfs inconsistent ',kmax1_gfs,ifind(1)-1
       stop
      endif 
      print*,'kmax1_gfs=',kmax1_gfs
      
      call aq_find(nspecies,' ',bndname,lpsec,iflag)   ! BND species
      noutbnd=lpsec-1
      call aq_find(nspecies,' ',checkname,lpsec,iflag)   ! BND species
      ncheck=lpsec-1


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
      jmax=jmax1-nhalo*2
      if(nlon.ne.imax.or.nlat.ne.jmax) then
        print*,'dimension mismatch ',nlon,imax,nlat,jmax
	stop
      endif
      
! read ps
      call check(nf90_inq_dimid(ncid,'lev',iddim_lev))
      call check(nf90_inquire_dimension(ncid,iddim_lev,len=kmax))
      call check(nf90_inq_dimid(ncid,'levp',iddim_levp))
      call check(nf90_inquire_dimension(ncid,iddim_levp,len=kmax1))
      if(kmax.ne.65) then
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

! -open HCMAQ file
      if(.not.open3('HCMAQ',FSREAD3,'my')) stop
      if(.not.desc3('HCMAQ')) stop
      ihcmaq=ncols3d
      jhcmaq=nrows3d
      khcmaq=nlays3d
      allocate(p_hcmaq(ihcmaq,jhcmaq,khcmaq),dens_hcmaq(ihcmaq,jhcmaq,khcmaq),&
        vhcmaq(ihcmaq,jhcmaq,khcmaq),tmpa(khcmaq))
      
      if(.not.read3('HCMAQ','PRES',ALLAYS3,sdate3d,stime3d,p_hcmaq)) stop
      if(.not.read3('HCMAQ','DENS',ALLAYS3,sdate3d,stime3d,dens_hcmaq)) stop
            
      iflag=setenvvar('IOAPI_ISPH','20') ! same as WRF R=6370km
      if(.not.setpol(sngl(p_alp3d),sngl(p_bet3d),sngl(p_gam3d),sngl(xcent3d),sngl(ycent3d))) stop

!---calculating lateral boundary horizontal index in HCMAQ coordinate	

! --- top and bottom

       do i=1,imax
        ix=i
	do j=1,nhalo
 	 do n=1,2
	  
	  select case(n)
	   case(1)	 
	    jy=j      ! _bottom
	   case(2) 
	    jy=jmax1-nhalo+j ! top
	  end select    
         
	  if(.not.LL2POL(xlon(ix,jy),xlat(ix,jy),x,y)) stop
	  bndcoordx(i,j,n,1)=(x-sngl(xorig3d))/sngl(xcell3d)+0.5 ! grid center
	  bndcoordx(i,j,n,2)=(y-sngl(yorig3d))/sngl(ycell3d)+0.5
	  
!---vertical index for top/bottom     
	  x=bndcoordx(i,j,n,1) 	 
	  y=bndcoordx(i,j,n,2)
	  xratio=x-int(x)
	  yratio=y-int(y)
	 
	  do kp=1,khcmaq
      	   tmpa(kp)=(1-yratio)*(p_hcmaq(int(x),int(y),kp)*      &    ! horizontally interpolate pressure
     	    (1-xratio)+p_hcmaq(int(x)+1,int(y),kp)*xratio)+     &
     	    yratio*(p_hcmaq(int(x),int(y)+1,kp)*(1-xratio)+     &
     	    p_hcmaq(int(x)+1,int(y)+1,kp)*xratio)
          enddo

          do k=1,kmax
	   if(pmid_x(i,j,k,n).ge.tmpa(1)) then  ! surface
	     bndcoordx(i,j,n,k+2)=1.  ! k index start from 3
	   else  
 	     do kp=2,khcmaq
 	      if(pmid_x(i,j,k,n).ge.tmpa(kp).and.pmid_x(i,j,k,n).le.tmpa(kp-1)) then 
                bndcoordx(i,j,n,k+2)=kp-1+(pmid_x(i,j,k,n)-tmpa(kp-1))/  &
     	          (tmpa(kp)-tmpa(kp-1))
                 exit
 	      endif  
	     enddo
           endif
           if(pmid_x(i,j,k,n).le.tmpa(khcmaq)) bndcoordx(i,j,n,k+2)=real(khcmaq)  ! model top  
          enddo
        
	enddo
       enddo 	
      enddo

! --- left and right
       
      do j=1,jmax
       jy=j+nhalo
           
	do i=1,nhalo
	 do n=1,2
	  if (n.eq.1) then
	   ix=i      !   left
	  else if (n.eq.2) then
	   ix=imax-nhalo+i ! right
	  endif    
          
	  if(.not.LL2POL(xlon(ix,jy),xlat(ix,jy),x,y)) stop
	  bndcoordy(i,j,n,1)=(x-sngl(xorig3d))/sngl(xcell3d)+0.5 ! grid center
	  bndcoordy(i,j,n,2)=(y-sngl(yorig3d))/sngl(ycell3d)+0.5
	  
!---vertical index for left/right   
	  x=bndcoordy(i,j,n,1) 	 
	  y=bndcoordy(i,j,n,2)
	  xratio=x-int(x)
	  yratio=y-int(y)
	 
	  do kp=1,khcmaq
      	   tmpa(kp)=(1-yratio)*(p_hcmaq(int(x),int(y),kp)*      &    ! horizontally interpolate height
     	    (1-xratio)+p_hcmaq(int(x)+1,int(y),kp)*xratio)+     &
     	    yratio*(p_hcmaq(int(x),int(y)+1,kp)*(1-xratio)+     &
     	    p_hcmaq(int(x)+1,int(y)+1,kp)*xratio)
          enddo

          do k=1,kmax
	   if(pmid_y(i,j,k,n).ge.tmpa(1)) then ! surface
	     bndcoordy(i,j,n,k+2)=1.  ! k index start from 3
	   else  
 	     do kp=2,khcmaq
 	      if(pmid_y(i,j,k,n).ge.tmpa(kp).and.pmid_y(i,j,k,n).le.tmpa(kp-1)) then 
                bndcoordy(i,j,n,k+2)=kp-1+(pmid_y(i,j,k,n)-tmpa(kp-1))/  &
     	          (tmpa(kp)-tmpa(kp-1))
                 exit
 	      endif  
	     enddo
           endif
           if(pmid_y(i,j,k,n).le.tmpa(khcmaq)) bndcoordy(i,j,n,k+2)=real(khcmaq)     
          enddo

         enddo
        enddo
       enddo
		  
	  
       if(iprint.eq.1) then
!         open(27,file='dust2.bin',form='unformatted',access='direct',recl=ihcmaq*jhcmaq*4)
         open(27,file='o3-tmp.bin',form='unformatted',access='direct',recl=ihcmaq*jhcmaq*4)
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
   
   do L=1,noutbnd
    
    ctmp=bndname(L)
    call upcase(ctmp)
    L1=index1(ctmp,nvars3d,vname3d)
    if(L1.le.0) then
      print*,'can not find ',ctmp
      stop 9
    endif
    if(.not.read3('HCMAQ',ctmp,ALLAYS3,sdate3d,stime3d,vhcmaq)) stop
    
    if(index(units3d(L1),'ug m-3').ge.1) vhcmaq(:,:,:)=vhcmaq(:,:,:)/dens_hcmaq(:,:,:)  ! ug/m3 -> ug/kg 

    vhcmaq(:,:,:)=amax1(0.,vhcmaq(:,:,:)) ! remove negative value	
    if(ctmp.eq.'O3'.and.o3cap.gt.0)  vhcmaq(:,:,:)=amin1(o3cap,vhcmaq(:,:,:)) ! capped o3 in ppmv	
 
    if(iprint.eq.1) then
      if(ctmp.eq.'O3') write(27,rec=1) vhcmaq(:,:,1)
      if(ctmp.eq.'CO') write(27,rec=2) vhcmaq(:,:,1)
    endif
	
     do i=1,imax       ! for top/bottom boundary conditions
      do j=1,nhalo
       do n=1,2
        x=bndcoordx(i,j,n,1)
	y=bndcoordx(i,j,n,2)
	xratio=x-int(x)
	yratio=y-int(y)

	tmpa(1:khcmaq)=(1-yratio)*(vhcmaq(int(x),int(y),     & ! horizontally interpolate values
     	 1:khcmaq)*(1-xratio)+vhcmaq(int(x)+1,int(y),       &
     	 1:khcmaq)*xratio)+yratio*(vhcmaq(int(x),int(y)+1,  &
     	 1:khcmaq)*(1-xratio)+vhcmaq(int(x)+1,int(y)+1,     &
     	 1:khcmaq)*xratio)
	
	 do k=1,kmax
	  z=bndcoordx(i,j,n,k+2)
	  if(z.ge.khcmaq) then
	   tmpvalue=tmpa(int(z))
	  else 
	   zratio=z-int(z)
	   tmpvalue=(1-zratio)*tmpa(int(z))+zratio*tmpa(int(z)+1)     ! vertically interpolate values
	  endif 

	  bndx(i,j,k,n,L)=amax1(tmpvalue,0.)
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

	tmpa(1:khcmaq)=(1-yratio)*(vhcmaq(int(x),int(y),     & ! horizontally interpolate values
     	 1:khcmaq)*(1-xratio)+vhcmaq(int(x)+1,int(y),       &
     	 1:khcmaq)*xratio)+yratio*(vhcmaq(int(x),int(y)+1,  &
     	 1:khcmaq)*(1-xratio)+vhcmaq(int(x)+1,int(y)+1,     &
     	 1:khcmaq)*xratio)
	
	 do k=1,kmax
	  z=bndcoordy(i,j,n,k+2)
	  if(z.ge.khcmaq) then
	   tmpvalue=tmpa(int(z))
	  else 
	   zratio=z-int(z)
	   tmpvalue=(1-zratio)*tmpa(int(z))+zratio*tmpa(int(z)+1)     ! vertically interpolate values
	  endif 

	  bndy(i,j,k,n,L)=amax1(tmpvalue,0.)
	 enddo
	 
	enddo
       enddo
      enddo
      
    enddo  ! end  species loop for one file

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
