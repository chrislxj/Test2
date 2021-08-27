
#include <define.h>

SUBROUTINE info_gridcell ( lon_points,lat_points,edgen,edgee,edges,edgew, &
                           nrow_start,nrow_end,ncol_start,ncol_end, &
                           nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                           latn,lats,lonw,lone,&
                           sinn,sins,lonw_rad,lone_rad,sinn_i,sins_i,lonw_rad_i,lone_rad_i,&
                           READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB )
! ----------------------------------------------------------------------
! Creates land gridcell infomation of original "raw" data files -
!     data with 30 x 30 arc seconds resolution
!     and model gridcell
!
! Created by Yongjiu Dai, 02/2014
!  
! ________________
! REVISION HISTORY:
!   /07/2014, Siguang Zhu & Xiangxiang Zhang: calculate the start/end
!               row/column number of fine resolution grid for user defined
!               model grid (low resolution).
!
!   Corrected by Nan Wei, 2020.02
! ----------------------------------------------------------------------
use precision

IMPLICIT NONE
! arguments:

      integer,  intent(in) :: lon_points        ! number of model longitude grid points
      integer,  intent(in) :: lat_points        ! number of model latitude grid points
      real(r8), intent(in) :: edgen             ! northern edge of the user defined domain (degrees)
      real(r8), intent(in) :: edgee             ! eastern edge of the user defined domain (degrees)
      real(r8), intent(in) :: edges             ! southern edge of the user defined domain (degrees)
      real(r8), intent(in) :: edgew             ! western edge of the user defined domain (degrees)
      real(r8), intent(in) :: latn(lat_points)  ! grid cell latitude, northern edge (deg)
      real(r8), intent(in) :: lats(lat_points)  ! grid cell latitude, southern edge (deg)
      real(r8), intent(in) :: lonw(lon_points)  ! grid cell longitude, western edge (deg)
      real(r8), intent(in) :: lone(lon_points)  ! grid cell longitude, eastern edge (deg)

      integer,  parameter :: nlat = 21600       ! 180*(60*2)
      integer,  parameter :: nlon = 43200       ! 360*(60*2)
      real(r8), parameter :: edgen_i = 90.      ! northern edge of fine grids (deg)
      real(r8), parameter :: edges_i = -90.     ! southern edge of fine grids (deg)
      real(r8), parameter :: edgew_i = -180.    ! western edge of fine grids (deg)
      real(r8), parameter :: edgee_i = 180.     ! eastern edge of fine grids (deg)

      integer,  intent(out) :: nrow_start
      integer,  intent(out) :: nrow_end
      integer,  intent(out) :: ncol_start
      integer,  intent(out) :: ncol_end
      integer,  intent(out) :: nx_fine_gridcell
      integer,  intent(out) :: ny_fine_gridcell
      real(r8), intent(out) :: sinn(lat_points)               ! grid cell latitude, northern edge(sin)
      real(r8), intent(out) :: sins(lat_points)               ! grid cell latitude, southern edge(sin)
      real(r8), intent(out) :: lonw_rad(lon_points)           ! grid cell longitude, western edge (radian)
      real(r8), intent(out) :: lone_rad(lon_points)           ! grid cell longitude, eastern edge (radian)
      real(r8), intent(out) :: sinn_i(nlat)                   ! fine grid cell latitude, northern edge(sin)
      real(r8), intent(out) :: sins_i(nlat)                   ! fine grid cell latitude, southern edge(sin)
      real(r8), intent(out) :: lonw_rad_i(nlon)               ! fine grid cell longitude, western edge (radian)
      real(r8), intent(out) :: lone_rad_i(nlon)               ! fine grid cell longitude, eastern edge (radian)
      real(r8), intent(out) :: area_fine_gridcell(nlon,nlat)  ! fine grid cell area (km**2)
      integer,  intent(out) :: READ_row_UB(lat_points)        ! north boundary index of a gridcell in fine grids
      integer,  intent(out) :: READ_col_UB(lon_points)        ! west boundary index of a gridcell in fine grids
      integer,  intent(out) :: READ_row_LB(lat_points)        ! south boundary index of a gridcell in fine grids
      integer,  intent(out) :: READ_col_LB(lon_points)        ! east boundary index of a gridcell in fine grids

      real(r8) :: latn_i(nlat)                  ! fine grid cell latitude, northern edge (deg)
      real(r8) :: lats_i(nlat)                  ! fine grid cell latitude, southern edge (deg)
      real(r8) :: lonw_i(nlon)                  ! fine grid cell longitude, western edge (deg)
      real(r8) :: lone_i(nlon)                  ! fine grid cell longitude, eastern edge (deg)
      real(r8) :: lat_i(nlat)                   ! fine grid cell latitude, center location (deg)
      real(r8) :: lon_i(nlon)                   ! fine grid cell longitude, center location (deg)

! ---------------------------------------------------------------
      real(r8) dx
      real(r8) dy
 
      real(r8) deg2rad                          ! pi/180
      real(r8) pi                               ! 3.14159265358979323846

      integer i, j
      integer tmp_lat(lat_points) 
      integer tmp_lon(lon_points) 

      integer, external :: nearest_boundary

! ---------------------------------------------------------------
! define the edges and area of the RAW DATA gridcells at 30 arc-seconds resolution
! ---------------------------------------------------------------
      pi = 4.*atan(1.)
      deg2rad = pi/180.
      dx = (edgee_i-edgew_i)/nlon   ! = 1./120.
      dy = (edgen_i-edges_i)/nlat   ! = 1./120.

#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(j)
#endif
      do j = 1, nlat
         lat_i(j)  = edgen_i - (2*j-1)*dy/2.
         latn_i(j) = edgen_i - (j-1)*dy
         lats_i(j) = edgen_i - j*dy
         sinn_i(j) = sin(latn_i(j)*deg2rad)
         sins_i(j) = sin(lats_i(j)*deg2rad)
      enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
      
#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i)
#endif
      do i = 1, nlon
         lon_i(i)  = edgew_i + (2*i-1)*dx/2.
         lonw_i(i) = edgew_i + (i-1)*dx
         lone_i(i) = edgew_i + i*dx
         lonw_rad_i(i) =  lonw_i(i)*deg2rad
         lone_rad_i(i) =  lone_i(i)*deg2rad
      enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      do j = 1,lat_points
         sinn(j) = sin(latn(j)*deg2rad)
         sins(j) = sin(lats(j)*deg2rad)
      enddo

      do i = 1,lon_points
         lonw_rad(i) =  lonw(i)*deg2rad
         lone_rad(i) =  lone(i)*deg2rad
      enddo


! calculate the area of the grid-cell of RAW DATA (fine-gridcell)
      call cellarea(nlat,nlon,latn_i,lats_i,lonw_i,lone_i,&
                    edgen_i,edgee_i,edges_i,edgew_i,&
                    area_fine_gridcell)

! --------------------------------------------------------
! define the starting and ending points and the numbers of 
! the RAW DATA fine gridcell in model grids
! --------------------------------------------------------

#if(defined USE_POINT_DATA)

      nrow_start = min(floor((90.-edgen)/dy) + 1, nlat)
      nrow_end   = nrow_start

      ncol_start = mod(floor((180.+edgew)/dx),nlon) + 1
      ncol_end   = ncol_start 

      nx_fine_gridcell = 1
      ny_fine_gridcell = 1

      READ_row_UB(1) = nrow_start
      READ_row_LB(1) = nrow_end
      READ_col_UB(1) = ncol_start
      READ_col_LB(1) = ncol_end

#else
        
      do i = 1,lat_points
         READ_row_UB(i) = nearest_boundary(latn(i),lat_i,nlat,'U')
      enddo
      
      do i = 1,lat_points
         READ_row_LB(i) = nearest_boundary(lats(i),lat_i,nlat,'L')
      enddo
    
      do i = 1,lon_points
         READ_col_UB(i) = nearest_boundary(lonw(i),lon_i,nlon,'U')
      enddo
      
      do i = 1,lon_points
         READ_col_LB(i) = nearest_boundary(lone(i),lon_i,nlon,'L')
         if (READ_col_LB(i) == READ_col_UB(i) .and. lone(i) < lonw(i)) then ! the model grid is too broad so that the starting
             READ_col_LB(i) = READ_col_LB(i) + nlon                         ! and ending points lie in the same fine grid
         end if
      enddo
      
      ! find the maximum of the number of fine grids in a model grid
      tmp_lat = READ_row_LB - READ_row_UB
      
      do i = 1,lon_points
         tmp_lon(i) = READ_col_LB(i) - READ_col_UB(i)
         if(READ_col_LB(i) < READ_col_UB(i)) then  ! gridcell across the dateline
            tmp_lon(i) = tmp_lon(i) + nlon
         endif
      end do

      nx_fine_gridcell = maxval(tmp_lon) + 1
      ny_fine_gridcell = maxval(tmp_lat) + 1

      nrow_start = minval(READ_row_UB)
      nrow_end   = maxval(READ_row_LB)
      ncol_start = 1
      ncol_end   = nlon
#endif

END SUBROUTINE info_gridcell


integer function nearest_boundary(degree,degree_i,dim_i,edge)

!=======================================================================
! Find the boundary index of a model gridcell in RAW data resolution for aggregation.
! The fine grids which are partially or fully overlapped with a coarse grid are all included.
! edited by zsg 20140823
! corrected by Nan Wei 2020.02
!======================================================================
   use precision
   implicit none
   
   integer  dim_i
   integer  tmp(1)
   real(r8) degree   
   real(r8) degree_i(dim_i)
   real(r8) diff_degree(dim_i)
   character edge

   ! find the first nearest fine grid as the model grid boundary's index
   diff_degree(:) = abs(degree - degree_i(:))

   if (edge == 'U')then
       diff_degree(1:dim_i) = diff_degree(dim_i:1:-1)
       tmp = minloc(diff_degree)
       if(tmp(1) > 1 .and. abs(diff_degree(tmp(1)) - diff_degree(tmp(1)-1)) < 1.0e-8) then
          tmp(1) = tmp(1) - 1     ! adjust the index to avoid the rounding error in calculating lon/lat info of fine grids
       end if
       nearest_boundary = dim_i - tmp(1) + 1
   else
       tmp = minloc(diff_degree)
       if(tmp(1) > 1 .and. abs(diff_degree(tmp(1)) - diff_degree(tmp(1)-1)) < 1.0e-8) then
          tmp(1) = tmp(1) - 1
       end if
       nearest_boundary = tmp(1) 
   end if

end function nearest_boundary 
!-----------------------------------------------------------------------
!EOP
