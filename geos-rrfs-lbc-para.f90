      program geos_bnd
!-------------------------------------------------------------------------------------      
!   interpolate GEOS Concentration to be lateral boundary condition for regional 
!   air quality model,  also output a layer result for checking purpose
!
!
!   Author: Youhua Tang
!   Revisions: GEOS5 monthly data to RRFS-CMAQ
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
      use mpi
      parameter(maxfile=300,nspecies=200, ngas1=33, ngas2=52, naerosol=25, &
       ngeos=ngas1+ngas2+naerosol)

      real sfact(ngeos,nspecies),val(nspecies),  &
         checkfact(ngeos,nspecies)

      double precision,allocatable  :: glon(:), glat(:)
      real,allocatable  :: zgeos(:,:,:),zfull_geos(:,:,:),tmpa(:),vgeos(:,:,:), press(:,:,:),  &
       worka(:),workb(:),workc(:), work(:), work1(:), work2(:), work3(:), xlat(:,:), xlon(:,:), &
       bndx(:,:,:,:,:), bndy(:,:,:,:,:),bndcoordx(:,:,:,:) ,bndcoordy(:,:,:,:),  &
       checkcoord(:,:,:),checksp(:,:,:), topo(:,:),topox(:,:,:),topoy(:,:,:), &
       airgeos(:,:,:),zhx(:,:,:,:),zhy(:,:,:,:),tmpbndx(:,:,:),tmpbndy(:,:,:)
      
      real bk_numatkn(65),bk_numacc(65),bk_numcor(65)  ! background aerosol number concentration
      
      character bndname(nspecies)*16,geosname(ngeos)*8,ctmp*16,  &
       echar(nspecies)*16,prefix(3)*200,suffix(3)*200,checkname(nspecies)*16,     &
       aline*200,gdatatype*4,modelname*4,gtype*16,arank*2, topofile*200, &
       lbcfile(2)*200,gas1name(ngas1)*8,gas2name(ngas2)*8, aeroname(naerosol)*8
      
      integer netindex(ngeos),checklayer,modate(maxfile),         &
       mosecs(maxfile),julian,ismotime(maxfile),iemotime(maxfile),  &
       idate(7),tlmeta,iret
      logical ingeos,lflag,extrameta,indexfind(nspecies)
      integer monthday(12),dimids(3)
      data monthday/31,30,31,30,31,30,31,31,30,31,30,31/

      data gas1name/'CH2O','ACET','ALD2','ALK4','C2H6','C3H8','Cl','Cl2','ClO','ClONO2',&
       'H2O2','HCl','HNO3','HNO3COND','HNO4','HO2','HOCl','ISOP','MO2','MP', &
       'MVK','N2O5','NO','NO2','NO3','CO','O1D','O3','OH','PAN','PRPE','R4N2','RCHO'/  ! Dac species end

      data gas2name/'A3O2','ACTA','ATO2','B3O2','EOH','ETO2','ETP','GCO3','GLYC','GLYX','GPAN', & ! 11
      'HAC','HCOOH','HNO2','IALD','IAO2','IAP','INO2','INPN','ISN1','ISNP','KO2','MACR','MAN2', & ! 13
      'MAO3','MAOP','MAP','MCO3','MEK','MGLY','MOH','MRO2','MRP','MVN2','PMN','PO2','PP', &       ! 13 
      'PPN','PRN1','PRPN','R4N1','R4O2','R4P','RA3P','RB3P','RIO1','RIO2','RIP','ROH','RP', &     ! 13
      'VRO2','VRP'/

      data aeroname/'BCPHILIC','BCPHOBIC','OCPHOBIC','OCPHILIC','SO2','SO2v',    &
       'SO4','SO4v','DMS','MSA','NH3','NH4a','NO3an1','NO3an2','NO3an3',   &
       'DU001','DU002','DU003','DU004','DU005','SS001','SS002','SS003','SS004','SS005'/
       
      data indexfind/nspecies*.false./
      
      integer  begyear,begdate,begtime,dtstep,numts,tstepdiff      
      namelist /control/iprint,bndname,tstepdiff,dtstep,o3cap,prefix,suffix, &	  !  input file preffix and suffix
       lbcfile,topofile,bk_numatkn,bk_numacc,bk_numcor
      
      CALL MPI_Init(ierr)
      CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
      CALL MPI_Comm_size(MPI_COMM_WORLD, npe, ierr) 
      
      call aq_blank(16*nspecies,bndname)
      call aq_blank(16*nspecies,checkname)

      sfact(1:ngeos,1:nspecies)=0.
      checkfact(1:ngeos,1:nspecies)=0.

      geosname(1:ngas1)=gas1name(1:ngas1)
      geosname(ngas1+1:ngas1+ngas2)=gas2name(1:ngas2)
      geosname(ngas1+ngas2+1:ngeos)=aeroname(1:naerosol)

! read converting information

      open(7,file='geos-rrfs-lbc-para.ini')
      read(7,control)

      call aq_find(nspecies,' ',bndname,lpsec,iflag)   ! BND species
      noutbnd=lpsec-1
      call aq_find(nspecies,' ',checkname,lpsec,iflag)   ! BND species
      ncheck=lpsec-1

      call aq_locate(7,'Species converting Factor',iflag)      

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
      allocate(xlon(imax,jmax1),topo(imax,jmax1))      
      
      call check(nf90_inq_varid(ncid,'geolon',idvar_geolon))
      call check(nf90_get_var(ncid,idvar_geolon,xlon))
      call check(nf90_inq_varid(ncid,'geolat',idvar_geolat))
      call check(nf90_get_var(ncid,idvar_geolat,xlat))
      call check(nf90_inq_varid(ncid,'orog_raw',idvar_topo))
      call check(nf90_get_var(ncid,idvar_topo,topo))
      
      do i=1,imax
       do j=1,jmax1
       if(xlon(i,j).lt.0) xlon(i,j)=xlon(i,j)+360
       enddo
      enddo 
      call check(nf90_close(ncid))
      print*,'finish reading topofile' 

! open LBC file for rewrite
      nowstep=0+my_rank*dtstep
      jfhour=tstepdiff+nowstep
      write(aline,'(a,i3.3,a)')trim(lbcfile(1)),nowstep,trim(lbcfile(2))

      call check(nf90_open(trim(aline),nf90_write, ncid))
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
      
! read zh
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
      
      allocate(zhx(imax,nhalo,kmax1,2),zhy(nhalo,jmax,kmax1,2), &
      topox(imax,nhalo,2),topoy(nhalo,jmax,2))
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

      topox(1:imax,1:nhalo,1)=topo(1:imax,1:nhalo)         ! bottom
      topox(1:imax,1:nhalo,2)=topo(1:imax,jmax1-nhalo+1:jmax1) ! top
      topoy(1:nhalo,1:jmax,1)=topo(1:nhalo,1:jmax)  ! left
      topoy(1:nhalo,1:jmax,2)=topo(imax-nhalo+1:imax,1:jmax) ! right
	  
      print*,'read zh_bottom, zh_top'
      call check(nf90_inq_varid(ncid,'zh_bottom',idvar_zh_bottom))
      call check(nf90_get_var(ncid,idvar_zh_bottom,zhx(:,:,:,1)))
      call check(nf90_inq_varid(ncid,'zh_top',idvar_zh_top))
      call check(nf90_get_var(ncid,idvar_zh_top,zhx(:,:,:,2)))
      do k=1,kmax
        zhx(1:imax,1:nhalo,k,1:2)=0.5*(zhx(1:imax,1:nhalo,k,1:2)+zhx(1:imax,1:nhalo,k+1,1:2))  ! convert from interface ASL to layer ASL
      enddo
      print*,'imax,jmax,kmax=',imax,jmax,kmax
      print*,'zhx bottom min/max before topo=',minval(zhx(:,1,:,1)),maxval(zhx(:,1,:,1))
      print*,'zhx top min/max before topo=',minval(zhx(:,1,:,2)),maxval(zhx(:,1,:,2))
      do k=1,kmax
        zhx(:,:,k,:)=amax1(0.,zhx(:,:,k,:)-topox(:,:,:))  ! convert to AGL
      enddo
      print*,'zhx bottom min/max after topo=',minval(zhx(:,1,:,1)),maxval(zhx(:,1,:,1))
      print*,'zhx top min/max after topo=',minval(zhx(:,1,:,2)),maxval(zhx(:,1,:,2))
	
!      if(maxval(zhx(:,:,:,:)).gt.1e12) then 
!          print*,'1 zhx overflowed ', maxval(zhx(:,:,:,:))
!	  stop
!       endif 

      print*,'read zh_left, zh_right'
      call check(nf90_inq_varid(ncid,'zh_left',idvar_zh_left))
      call check(nf90_get_var(ncid,idvar_zh_left,zhy(:,:,:,1)))
      call check(nf90_inq_varid(ncid,'zh_right',idvar_zh_right))
      call check(nf90_get_var(ncid,idvar_zh_right,zhy(:,:,:,2)))
      do k=1,kmax
        zhy(:,:,k,:)=0.5*(zhy(:,:,k,:)+zhy(:,:,k+1,:))  ! convert from interface ASL to layer ASL
      enddo
      do k=1,kmax
       zhy(:,:,k,:)=amax1(0.,zhy(:,:,k,:)-topoy(:,:,:))  ! convert to AGL
      enddo
      
!       if(maxval(zhx(:,:,:,:)).gt.1e12) then 
!          print*,'2 zhx overflowed ', maxval(zhx(:,:,:,:))
!	  stop
!       endif

! -open GEOS5
     do m=1,3
      if(dtstep.eq.0) then
       aline=trim(prefix(m))//trim(suffix(m))
      else 
       write(aline,'(a,i3.3,a)')trim(prefix(m)),jfhour,trim(suffix(m))
      endif      
      
      if(m.gt.1) call check(nf90_close(id_file))
      print*,'open ',trim(aline)      
      call check(nf90_open(trim(aline),nf90_nowrite, id_file))
      call check(nf90_inq_dimid(id_file,'lat',id_nlat))
      call check(nf90_inquire_dimension(id_file,id_nlat,len=nlatgeos))
      call check(nf90_inq_dimid(id_file,'lon',id_nlon))
      call check(nf90_inquire_dimension(id_file,id_nlon,len=nlongeos))
      call check(nf90_inq_dimid(id_file,'lev',id_nlev))
      call check(nf90_inquire_dimension(id_file,id_nlev,len=nlevgeos))
      
      if(m.eq.1) then
       igeos=nlongeos
       jgeos=nlatgeos
       kgeos=nlevgeos
       
       allocate(glon(igeos),glat(jgeos),tmpa(kgeos))
       allocate(zgeos(igeos,jgeos,kgeos),zfull_geos(igeos,jgeos,kgeos+1))
       allocate(vgeos(igeos,jgeos,kgeos),airgeos(igeos,jgeos,kgeos))
      else
       if(nlongeos.ne.igeos.or.nlatgeos.ne.jgeos.or.nlevgeos.ne.kgeos) then
        print*,'inconsisent dimensions '
	print*,'igeos,jgeos,kgeos=',igeos,jgeos,kgeos
	print*,'nlongeos,nlatgeos,nlevgeos=',nlongeos,nlatgeos,nlevgeos
	stop
       endif
      endif 
             
      print*,'opened GEOS5 file ',trim(aline),igeos,jgeos,kgeos 
      if(m.eq.1) then
        call check(nf90_inq_varid(id_file,'lon',id_lon))
        call check(nf90_get_var(id_file,id_lon,glon))
	call check(nf90_inq_varid(id_file,'lat',id_lat))
        call check(nf90_get_var(id_file,id_lat,glat))
	call check(nf90_inq_varid(id_file,'DELP',id_delp))
	call check(nf90_get_var(id_file,id_delp,zgeos))   ! store delp in zgeos
	call check(nf90_inq_varid(id_file,'AIRDENS',id_airdens))
	call check(nf90_get_var(id_file,id_airdens,airgeos))
	zfull_geos(:,:,kgeos+1)=0.  ! top down
	do k=kgeos,1,-1
	 zfull_geos(:,:,k)=zfull_geos(:,:,k+1)+zgeos(:,:,k)/9.8/airgeos(:,:,k)
	enddo
	do k=1,kgeos
	 zgeos(:,:,k)=0.5*(zfull_geos(:,:,k)+zfull_geos(:,:,k+1))
	enddo  
	glonint=sngl((glon(igeos)-glon(1))/(igeos-1))
	glatint=sngl((glat(jgeos)-glon(1))/(jgeos-1))
	
	ishift=0
	if((glon(1)+180).lt.2) then	  
	  do i=1,igeos
	   if(abs(glon(i)).lt.glonint/4) exit  ! shift to start from lon=0.0 for North America 
	  enddo	  
	  if(i.le.igeos) then
	   ishift=i-1
	   glon(i)=0   ! fix the precision issue for lon
	  endif 
	  where(glon.lt.0) glon=glon+360
	  if(i.gt.1) then
	   glon=cshift(glon,shift=ishift,dim=1)
	   zgeos=cshift(zgeos,shift=ishift,dim=1)
	   airgeos=cshift(airgeos,shift=ishift,dim=1)
	   print*,ishift,'after shift, glon=',glon
	  endif 
	endif

!---calculating lateral boundary horizontal index in geos coordinate	

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
         
	  do i2=1,igeos-1
           if(xlon(ix,jy).ge.glon(i2).and.xlon(ix,jy).le.glon(i2+1)) then
	    bndcoordx(i,j,n,1)=i2+(xlon(ix,jy)-glon(i2))/     &   ! i in geos coordiate
      	    (glon(i2+1)-glon(i2))  
	    exit
	   endif
	  enddo
	 
          do j2=1,jgeos
           if(xlat(ix,jy).ge.glat(j2).and.xlat(ix,jy).le.glat(j2+1)) then
	    bndcoordx(i,j,n,2)=j2+(xlat(ix,jy)-glat(j2))/     &   ! j in geos coordiate
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
      	   tmpa(kp)=(1-yratio)*(zgeos(int(x),int(y),kp)*      &    ! horizontally interpolate height
     	    (1-xratio)+zgeos(int(x)+1,int(y),kp)*xratio)+     &
     	    yratio*(zgeos(int(x),int(y)+1,kp)*(1-xratio)+     &
     	    zgeos(int(x)+1,int(y)+1,kp)*xratio)
          enddo

          do k=1,kmax
	   if(zhx(i,j,k,n).le.tmpa(kgeos)) then
	     bndcoordx(i,j,n,k+2)=real(kgeos)  ! k index start from 3
	   else  
 	     do kp=2,kgeos
 	      if(zhx(i,j,k,n).ge.tmpa(kp).and.zhx(i,j,k,n).le.tmpa(kp-1)) then  ! both use top-down coordinate
                bndcoordx(i,j,n,k+2)=kp-1+(tmpa(kp-1)-zhx(i,j,k,n))/  &
     	          (tmpa(kp-1)-tmpa(kp))
                 exit
 	      endif  
	     enddo
           endif
           if(zhx(i,j,k,n).ge.tmpa(1)) bndcoordx(i,j,n,k+2)=1.     
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
         
	 do i2=1,igeos-1
          if(xlon(ix,jy).ge.glon(i2).and.xlon(ix,jy).le.glon(i2+1)) then
	   bndcoordy(i,j,n,1)=i2+(xlon(ix,jy)-glon(i2))/     &   ! i in geos coordiate
      	    (glon(i2+1)-glon(i2))  
	   exit
	  endif
	 enddo
	 
         do j2=1,jgeos-1
 	  if(xlat(ix,jy).ge.glat(j2).and.xlat(ix,jy).le.glat(j2+1)) then
	   bndcoordy(i,j,n,2)=j2+(xlat(ix,jy)-glat(j2))/     &   ! j, in geos coordiate
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
      	   tmpa(kp)=(1-yratio)*(zgeos(int(x),int(y),kp)*      &    ! horizontally interpolate height
     	    (1-xratio)+zgeos(int(x)+1,int(y),kp)*xratio)+     &
     	    yratio*(zgeos(int(x),int(y)+1,kp)*(1-xratio)+     &
     	    zgeos(int(x)+1,int(y)+1,kp)*xratio)
          enddo

          do k=1,kmax
	   if(zhy(i,j,k,n).le.tmpa(kgeos)) then
	     bndcoordy(i,j,n,k+2)=real(kgeos)  ! k index start from 3
	   else  
 	     do kp=2,kgeos
 	      if(zhy(i,j,k,n).ge.tmpa(kp).and.zhy(i,j,k,n).le.tmpa(kp-1)) then  ! both use top-down coordinate
                bndcoordy(i,j,n,k+2)=kp-1+(tmpa(kp-1)-zhy(i,j,k,n))/  &
     	          (tmpa(kp-1)-tmpa(kp))
                 exit
 	      endif  
	     enddo
           endif
           if(zhy(i,j,k,n).ge.tmpa(1)) bndcoordy(i,j,n,k+2)=1.     
          enddo

         enddo
        enddo
       enddo
	  
     endif ! m=1	  
	  
       if(iprint.eq.1.and.my_rank.eq.0.and.m.eq.1) then
!         open(27,file='dust2.bin',form='unformatted',access='direct',recl=igeos*jgeos*4)
         open(27,file='o3-tmp.bin',form='unformatted',access='direct',recl=igeos*jgeos*4)
	 open(28,file='zhx-tmp.bin',form='unformatted',access='direct',recl=imax*kmax1*4)
	 write(28,rec=1)zhx(1:imax,1,1:kmax1,1)
	 write(28,rec=2)zhx(1:imax,1,1:kmax1,2)
	 print*,'zhx2 top =',zhx(1:imax,1,1,2)
	 print*,'zhx bottom min/max after print=',minval(zhx(:,1,:,1)),maxval(zhx(:,1,:,1))
         print*,'zhx top min/max after print=',minval(zhx(:,1,:,2)),maxval(zhx(:,1,:,2))
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
   do L=1,nowspc
    if(m.eq.1) then 
      ctmp=gas1name(L)
      L1=L
    else if(m.eq.2) then
     ctmp=gas2name(L)
     L1=L+ngas1
    else if(m.eq.3) then
     ctmp=aeroname(L)
     L1=L+ngas1+ngas2  ! index in whole species
    endif
    call check(nf90_inq_varid(id_file,trim(ctmp),id_var1))
    call check(nf90_get_var(id_file,id_var1,vgeos))
    if(ishift.gt.0) vgeos=cshift(vgeos,shift=ishift,dim=1)
     
    if(m.le.2) then
       vgeos(:,:,:)=vgeos(:,:,:)*1e6   ! gaseous: mole/mole to ppmV
    else if(m.eq.3) then 
       if(ctmp.eq.'SO2'.or.ctmp.eq.'SO2v') then
	  vgeos(:,:,:)=vgeos(:,:,:)/64.066*28.97*1e6  ! kg/kg to ppmV
       else if(ctmp.eq.'NH3') then
	  vgeos(:,:,:)=vgeos(:,:,:)/17.031*28.97*1e6  ! to ppmV
       else if(ctmp.eq.'MSA') then
	  vgeos(:,:,:)=vgeos(:,:,:)/96.10*28.97*1e6  ! to ppmV
       else if (ctmp.eq.'DMS') then
	  vgeos(:,:,:)=vgeos(:,:,:)/62.13*28.97*1e6
       else ! aerosols kg/kg -> ug/kg
	  vgeos(:,:,:)=vgeos(:,:,:)*1e9
       endif 
    endif 
    vgeos(:,:,:)=amax1(0.,vgeos(:,:,:)) ! remove negative value	
    if(ctmp.eq.'O3'.and.o3cap.gt.0)  vgeos(:,:,:)=amin1(o3cap,vgeos(:,:,:)) ! capped o3 in ppmv	
 
    if(iprint.eq.1.and.my_rank.eq.0) then
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
	   bndx(i,j,k,n,L2)=bndx(i,j,k,n,L2)+amax1(tmpvalue,0.)*sfact(L1,L2)
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
	   bndy(i,j,k,n,L2)=bndy(i,j,k,n,L2)+amax1(tmpvalue,0.)*sfact(L1,L2)
	  enddo
	 enddo
	 
	enddo
       enddo
      enddo
      
    enddo  ! end  geos species loop for one file
  enddo ! end of m loop
      print*,'2 bndx(numcor) min,max=',minval(bndx(:,:,:,:,L_numcor)),maxval(bndx(:,:,:,:,L_numcor))
      print*,'2 bndy(numcor) min,max=',minval(bndy(:,:,:,:,L_numcor)),maxval(bndy(:,:,:,:,L_numcor))

  
! begin output
     fillval=0./0.
     print*,'start write' 
      do L=1,noutbnd       ! check if geos supplies all species, otherwise do not overwrite the existing aerosol

        do m=1,2

!! -top/bottom	  
	  if(sum(bndx(1:imax,1:nhalo,1:kmax,m,L)).le.1e-20 ) then
	    print*, 'm=,', m,' X skip for ', bndname(L), 'my_rank=',my_rank
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
	    print*,'add variable ',trim(aline),my_rank
	    call check(nf90_def_var(ncid,trim(aline),nf90_real,dimids,idvar_tmp))
	    call check(nf90_put_att(ncid,idvar_tmp,'_FillValue',fillval))
	    call check(nf90_enddef(ncid))
	    call check(nf90_inq_varid(ncid,trim(aline),idvar_tmp)) ! check exist
	  endif
	  
	  print*,my_rank,'write ',trim(aline),idvar_tmp,minval(bndx(1:imax,1:nhalo,1:kmax,m,L)),maxval(bndx(1:imax,1:nhalo,1:kmax,m,L))
	  
!	  if(index(bndname(L),'num').gt.0) then
!!	    call check(nf90_get_var(ncid,idvar_tmp,tmpbndx))
!	    tmpbndx(1:imax,1:nhalo,1:kmax)=tmpbndx(1:imax,1:nhalo,1:kmax)+ &
!	      bndx(1:imax,1:nhalo,1:kmax,m,L)
!	  else
	    tmpbndx(1:imax,1:nhalo,1:kmax)=bndx(1:imax,1:nhalo,1:kmax,m,L)
	    if(index(bndname(L),'aecj').gt.0) then
	     do k=1,4
	      tmpbndx(1:imax,1:nhalo,k)=tmpbndx(1:imax,1:nhalo,5) ! for bug in GEFS EC 
	     enddo
	    endif   
!	  endif
	  call check(nf90_put_var(ncid,idvar_tmp,tmpbndx))
	  
!! left/right
	  if(sum(bndy(1:nhalo,1:jmax,1:kmax,m,L)).le.1e-20 ) then
	    print*, 'm=,',m, ' Y skip for ', bndname(L), 'my_rank=',my_rank
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

!	  if(index(bndname(L),'num').gt.0) then
!	    call check(nf90_get_var(ncid,idvar_tmp,tmpbndy))
!	    tmpbndy(1:nhalo,1:jmax,1:kmax)=tmpbndy(1:nhalo,1:jmax,1:kmax)+ &
!	      bndy(1:nhalo,1:jmax,1:kmax,m,L)
!	  else
	    tmpbndy(1:nhalo,1:jmax,1:kmax)=bndy(1:nhalo,1:jmax,1:kmax,m,L)
	    if(index(bndname(L),'aecj').gt.0) then
	     do k=1,4
	      tmpbndy(1:nhalo,1:jmax,k)=tmpbndy(1:nhalo,1:jmax,5) ! for bug in GEFS EC
	     enddo
	    endif 
!	  endif    
	    call check(nf90_put_var(ncid,idvar_tmp,tmpbndy))         
	 enddo
	enddo
      call check(nf90_close(ncid))
      call MPI_Finalize(ierr)

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
	read(iunit,100,end=8)  (char(i),i=1,80)
	if(char(1).eq.'$') then
	  write(6,200) iunit,(char(i),i=1,80)
        else if(char(1).eq.'#') then
	else
	  backspace iunit
	  return
        endif
       end do
 8    continue
100   format(80a1)
200   format(2x,'iunit=',i3,2x,80a1)
      end
