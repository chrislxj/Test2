
#include <define.h>

SUBROUTINE aggregation_LAI( dir_rawdata,dir_model_landdata, &
                            lon_points,lat_points, &
                            nrow_start,nrow_end,ncol_start,ncol_end, &
                            nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                            sinn,sins,lonw_rad,lone_rad,sinn_i,sins_i,lonw_rad_i,lone_rad_i,&
                            READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB) 
! ----------------------------------------------------------------------
! 1. Global land cover types (updated with the specific dataset)
!
! 2. Global Plant Leaf Area Index
!    (http://globalchange.bnu.edu.cn)
!    Yuan H., et al., 2011:
!    Reprocessing the MODIS Leaf Area Index products for land surface 
!    and climate modelling. Remote Sensing of Environment, 115: 1171-1187.
!
! Created by Yongjiu Dai, 02/2014
!
! ________________
! REVISION HISTORY:
!   /07/2014, Siguang Zhu & Xiangxiang Zhang: weight average considering 
!               partial overlap between fine grid and model grid for a user
!               defined domain file.
!
!   Corrected by Nan Wei, 2020.02
! ----------------------------------------------------------------------
use precision

IMPLICIT NONE

! arguments:
#if(defined USGS_CLASSIFICATION)
      integer, parameter :: N_land_classification = 24 ! GLCC USGS number of land cover category
#endif
#if(defined IGBP_CLASSIFICATION)
      integer, parameter :: N_land_classification = 17 ! MODIS IGBP number of land cover category
#endif
      integer, parameter :: nlat = 21600  ! 180*(60*2)
      integer, parameter :: nlon = 43200  ! 360*(60*2)

      character(LEN=256), intent(in) :: dir_rawdata
      character(LEN=256), intent(in) :: dir_model_landdata

      integer, intent(in) :: lon_points ! number of model longitude grid points
      integer, intent(in) :: lat_points ! number of model latitude grid points
      integer, intent(in) :: nrow_start
      integer, intent(in) :: nrow_end
      integer, intent(in) :: ncol_start
      integer, intent(in) :: ncol_end
      integer, intent(in) :: nx_fine_gridcell
      integer, intent(in) :: ny_fine_gridcell
      
      real(r8), intent(in) :: sinn(lat_points)        ! grid cell latitude, northern edge(sin)  
      real(r8), intent(in) :: sins(lat_points)        ! grid cell latitude, southern edge(sin)
      real(r8), intent(in) :: lonw_rad(lon_points)    ! grid cell longitude, western edge (radian)
      real(r8), intent(in) :: lone_rad(lon_points)    ! grid cell longitude, eastern edge (radian)
      real(r8), intent(in) :: sinn_i(nlat)            ! fine grid cell latitude, northern edge(sin)
      real(r8), intent(in) :: sins_i(nlat)            ! fine grid cell latitude, southern edge(sin)
      real(r8), intent(in) :: lonw_rad_i(nlon)        ! fine grid cell longitude, western edge (radian)
      real(r8), intent(in) :: lone_rad_i(nlon)        ! fine grid cell longitude, eastern edge (radian)
      real(r8), intent(in) :: area_fine_gridcell(nlon,nlat)  ! rawdata fine cell area (km**2)      
      integer,  intent(in) :: READ_row_UB(lat_points) ! north boundary index in rawdata files for each model grid cell
      integer,  intent(in) :: READ_col_UB(lon_points) ! west boundary index in rawdata files for each model grid cell  
      integer,  intent(in) :: READ_row_LB(lat_points) ! south boundary index in rawdata files for each model grid cell
      integer,  intent(in) :: READ_col_LB(lon_points) ! east boundary index in rawdata files for each model grid cell

! local variables:
! ----------------------------------------------------------------------
      character(len=256) lndname
      character(len=1) land_chr1(nlon)
      character(len=2) land_chr2(nlon)
      integer(kind=1)  land_int1(nlon)
      integer(kind=2)  land_int2(nlon)

      integer iunit
      integer length
      integer i, j, L, i1, i2, j1, j2, isfirst
      integer nrow, ncol, ncol_mod
      integer LL, np
      integer n_fine_gridcell

      integer , allocatable :: num_patches(:) 
      real(r8), allocatable :: a_area_patches(:,:)

      integer,  allocatable :: landtypes(:,:) ! GLCC USGS/MODIS IGBP land cover types 
      real(r8), allocatable :: LAI(:,:)          ! plant leaf area index (m2/m2)
      real(r8), allocatable :: a_LAI_patches(:,:) 
      real(r8), allocatable :: LAI_patches(:,:,:)
      integer N8, Julian_day
      character(LEN=256) c

      real(r8), external :: find_min_area

! ........................................
! ... (1) gloabl land cover types
! ........................................
      iunit = 100
      inquire(iolength=length) land_chr1 
      allocate ( landtypes (nlon,nlat) )
 
#if(defined USE_POINT_DATA)

#if(defined USGS_CLASSIFICATION)
      landtypes(ncol_start,nrow_start) = USGS_CLASSIFICATION
#endif

#if(defined IGBP_CLASSIFICATION)
      landtypes(ncol_start,nrow_start) = IGBP_CLASSIFICATION
#endif

#else

#if(defined USGS_CLASSIFICATION)
     ! GLCC USGS classification
     ! -------------------
      lndname = trim(dir_rawdata)//'RAW_DATA_updated/landtypes_usgs_update.bin' 
      print*,lndname
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old') 
      do nrow = nrow_start, nrow_end
         read(iunit,rec=nrow,err=100) land_chr1 
         landtypes(:,nrow) = ichar(land_chr1(:)) 
      enddo 
      close (iunit)
#endif

#if(defined IGBP_CLASSIFICATION)
     ! MODIS IGBP classification
     ! -------------------
      lndname = trim(dir_rawdata)//'RAW_DATA_updated/landtypes_igbp_update.bin'
      print*,lndname
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
      do nrow = nrow_start, nrow_end
         read(iunit,rec=nrow,err=100) land_chr1
         landtypes(:,nrow) = ichar(land_chr1(:))
      enddo
      close (iunit)
#endif 

#endif

! ................................................
! ... (2) global plant leaf area index
! ................................................
      n_fine_gridcell = nx_fine_gridcell * ny_fine_gridcell

      allocate (LAI           (nlon,nlat))
      allocate (num_patches   (0:N_land_classification))
      allocate (a_area_patches(0:N_land_classification,1:n_fine_gridcell))
      allocate (a_LAI_patches (0:N_land_classification,1:n_fine_gridcell))
      allocate (LAI_patches   (0:N_land_classification,1:lon_points,1:lat_points))

      iunit = 100
      inquire(iolength=length) land_chr1

      DO N8 = 1, 46
         ! -----------------------
         ! read in leaf area index
         ! -----------------------
         Julian_day = 1 + (N8-1)*8
         write(c,'(i3.3)') Julian_day

         lndname = trim(dir_rawdata)//'lai/global_30s_10_year_avg/LAI_BNU_'//trim(c)
         print*,lndname

         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old') 

         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) land_chr1
            LAI(:,nrow) = ichar(land_chr1(:))*0.1
         enddo 
         close (iunit)
         print*, minval(LAI(:,nrow_start:nrow_end)), maxval(LAI(:,nrow_start:nrow_end))

         ! ---------------------------------------------------------------
         ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
         ! ---------------------------------------------------------------

#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,i1,i2,j1,j2,nrow,ncol,ncol_mod,L,LL,np) &
!$OMP PRIVATE(a_LAI_patches,a_area_patches,num_patches,isfirst)
#endif
         do j = 1, lat_points

            j1 = READ_row_UB(j)   ! read upper boundary of latitude 
            j2 = READ_row_LB(j)   ! read lower boundary of latitude

            do i = 1, lon_points

            i1 = READ_col_UB(i)   ! read upper boundary of longitude
            i2 = READ_col_LB(i)   ! read lower boundary of longitude

               num_patches(:) = 0

               do nrow = j1, j2
                  isfirst = 1
                  if(i1 > i2) i2 = i2 + nlon   ! for coarse grid crosses the dateline  
                  do ncol = i1, i2
                     ncol_mod = mod(ncol,nlon)
                     if(ncol_mod == 0) ncol_mod = nlon
                     
                     L = landtypes(ncol_mod,nrow)
                     num_patches (L) = num_patches(L) + 1
                     LL = num_patches(L)

                     a_LAI_patches (L,LL) = LAI(ncol_mod,nrow)
                     a_area_patches(L,LL) = find_min_area(lone_rad(i),lonw_rad(i),lone_rad_i(ncol_mod),&
                                    lonw_rad_i(ncol_mod),sinn(j),sins(j),sinn_i(nrow),sins_i(nrow),isfirst)
                     isfirst = 0 
                  enddo
               enddo
               
               do L = 0, N_land_classification
                  np = num_patches (L)
                  if (np == 0) then
                      LAI_patches(L,i,j) = 0.
                  else if(np == 1) then
                      LAI_patches(L,i,j) = a_LAI_patches(L,1)
                  else
                      LAI_patches(L,i,j) = sum(a_LAI_patches(L,1:np) * &
                      (a_area_patches(L,1:np)/sum(a_area_patches(L,1:np))))
                  end if
               enddo

            enddo
         enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

         ! ---------------------------------------------------
         ! write out the plant leaf area index of grid patches
         ! ---------------------------------------------------
         lndname = trim(dir_model_landdata)//'model_LAI_patches.'//trim(c)//'.bin'
         print*,lndname
         print*, maxval(LAI_patches),minval(LAI_patches)
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) LAI_patches
         close (iunit)

      ENDDO

      deallocate ( LAI )
      deallocate ( a_LAI_patches )
      deallocate ( LAI_patches )
      deallocate ( a_area_patches )
      deallocate ( num_patches )
      deallocate ( landtypes )

      go to 1000
100   print 101,nrow,lndname
101   format(' record =',i8,',  error occured on file: ',a50)
1000  continue

END SUBROUTINE aggregation_LAI
