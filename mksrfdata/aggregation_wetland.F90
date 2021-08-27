#include <define.h>
SUBROUTINE aggregation_wetland( dir_rawdata,dir_model_landdata, &
                                lon_points,lat_points, &
                                nrow_start,nrow_end,ncol_start,ncol_end, &
                                nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                                sinn,sins,lonw_rad,lone_rad,sinn_i,sins_i,lonw_rad_i,lone_rad_i,&
                                READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB)
! ----------------------------------------------------------------------
! 1. Global land cover types (updated with the specific dataset)
!
! 2. Global Lake and Wetlands Types
!     (http://www.wwfus.org/science/data.cfm)
!       1       Lake
!       2       Reservoir
!       3       River
!       4       Freshwater Marsh, Floodplain
!       5       Swamp Forest, Flooded Forest
!       6       Coastal Wetland (incl. Mangrove, Estuary, Delta, Lagoon)
!       7       Pan, Brackish/Saline Wetland
!       8       Bog, Fen, Mire (Peatland)
!       9       Intermittent Wetland/Lake
!       10      50-100% Wetland
!       11      25-50% Wetland
!       12      Wetland Complex (0-25% Wetland)
!
! Created by Yongjiu Dai, 02/2014
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
      integer,  intent(in) :: READ_row_UB(lat_points) ! north boundary index of a gridcell in fine grids
      integer,  intent(in) :: READ_col_UB(lon_points) ! west boundary index of a gridcell in fine grids
      integer,  intent(in) :: READ_row_LB(lat_points) ! south boundary index of a gridcell in fine grids
      integer,  intent(in) :: READ_col_LB(lon_points) ! east boundary index of a gridcell in fine grids

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
      integer LL, np, n, nn
      integer n_fine_gridcell

      integer, allocatable :: landtypes(:,:) ! GLCC USGS / MODIS IGBP land cover types 
      integer, allocatable :: lakewetland(:,:)  ! lake and wetland types
      integer, allocatable :: num_patches(:)
      integer, allocatable :: n_wetland_patches(:)
      real(r8), allocatable :: f_wetland(:)
      real(r8), allocatable :: area_wetland_patches(:)
      real(r8), allocatable :: fraction_wetland_patches(:,:,:)
      real(r8) area_for_sum

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
! ... (2) global lakes and wetland
! ................................................
      iunit = 100
      inquire(iolength=length) land_chr1
      lndname = trim(dir_rawdata)//'lake_wetland/glwd.bin'
      print*,lndname
      allocate ( lakewetland (nlon,nlat) )

      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old') 
      do nrow = nrow_start, nrow_end
         read(iunit,rec=nrow,err=100) land_chr1
         lakewetland(:,nrow) = ichar(land_chr1(:))
      enddo 
      close (iunit)
      print*,minval(lakewetland(:,nrow_start:nrow_end)), maxval(lakewetland(:,nrow_start:nrow_end))

!   ---------------------------------------------------------------
!   aggregate the wetland from the resolution of raw data to modelling resolution
!   ---------------------------------------------------------------
      n_fine_gridcell = nx_fine_gridcell * ny_fine_gridcell
      allocate ( num_patches(0:N_land_classification) )
      allocate ( fraction_wetland_patches(4:12,1:lon_points,1:lat_points) )
      allocate ( n_wetland_patches(n_fine_gridcell) )
      allocate ( f_wetland(4:12) )
      allocate ( area_wetland_patches(n_fine_gridcell) ) 

#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,i1,i2,j1,j2,nrow,ncol,ncol_mod,L,LL,num_patches,np,n,nn) &
!$OMP PRIVATE(n_wetland_patches,area_wetland_patches) &
!$OMP PRIVATE(f_wetland,area_for_sum,isfirst) 
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
                  
                  area_for_sum = find_min_area(lone_rad(i),lonw_rad(i),lone_rad_i(ncol_mod),&
                                 lonw_rad_i(ncol_mod),sinn(j),sins(j),sinn_i(nrow),sins_i(nrow),isfirst)
                  isfirst = 0               

                  L = landtypes(ncol_mod,nrow)
#if(defined USGS_CLASSIFICATION)
                  if(L==17)then    ! WETLAND (17,18), all herbaceous and wooded wetland types have been => 17
#endif
#if(defined IGBP_CLASSIFICATION)
                  if(L==11)then    ! WETLAND (11)
#endif
                     num_patches(L) = num_patches(L) + 1
                     LL = num_patches(L)
                     n_wetland_patches (LL) = lakewetland(ncol_mod,nrow)
                     area_wetland_patches(LL) = area_for_sum
                  endif
               enddo
            enddo

#if(defined USGS_CLASSIFICATION)
            np = num_patches(17)
#endif
#if(defined IGBP_CLASSIFICATION)
            np = num_patches(11)
#endif
          ! [Freshwater Marsh, Floodplain (4)] -> [Wetland Complex (12)]
            if(np == 0)then
               fraction_wetland_patches(:,i,j) = 0.0
            else if(np == 1) then
               fraction_wetland_patches(:,i,j) = 0.0
               nn = n_wetland_patches(1)
               if(nn>=4 .and. nn<=12) fraction_wetland_patches(nn,i,j) = 1.0
            else
               f_wetland(:) = 0.
               do n = 1, np
                  nn = n_wetland_patches(n)
                  if(nn>=4 .and. nn<=12)then
                     f_wetland(nn) = f_wetland(nn) + area_wetland_patches(n)
                  endif
               enddo
               fraction_wetland_patches(:,i,j) = f_wetland(:) / sum(area_wetland_patches(1:np))
            endif
         enddo
      enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

! ---------------------------------------------------
! write out the fraction of wetland patches
! ---------------------------------------------------
      lndname = trim(dir_model_landdata)//'model_wetland_types.bin'
      print*,lndname
      open(iunit,file=trim(lndname),form='unformatted',status='unknown')
      write(iunit,err=100) fraction_wetland_patches
      close (iunit)

      deallocate ( landtypes )
      deallocate ( lakewetland )
      deallocate ( num_patches )
      deallocate ( n_wetland_patches )
      deallocate ( f_wetland )
      deallocate ( area_wetland_patches ) 
      deallocate ( fraction_wetland_patches )

      go to 1000
100   print 101,nrow,lndname
101   format(' record =',i8,',  error occured on file: ',a50)
1000  continue

END SUBROUTINE aggregation_wetland
