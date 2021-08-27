#include <define.h>

SUBROUTINE aggregation_soil_parameters ( dir_rawdata,dir_model_landdata, &
                                         lon_points,lat_points, &
                                         nrow_start,nrow_end,ncol_start,ncol_end, &
                                         nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                                         sinn,sins,lonw_rad,lone_rad,sinn_i,sins_i,lonw_rad_i,lone_rad_i,&
                                         READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB )
! ----------------------------------------------------------------------
! Creates land model surface dataset from original "raw" data files -
!     data with 30 arc seconds resolution to model resolution
!
! Yongjiu Dai, 02/2014, 05/2018
! Nan Wei, 07/2018, 02/2020
! options of upscaling methods for soil parameters: 
!            (1)median, (2)area-weighted arithmetic average, (3)AAVW (Wang et al., 2004), 
!            (4)area-weighted geometric average, (5)fitted values based on SW retention curves/Ke-Sr relationship 
! ----------------------------------------------------------------------
use precision

#if(defined SOILPAR_UPS_FIT)
use par_fitting
#endif

IMPLICIT NONE

! arguments:

#if(defined USGS_CLASSIFICATION)
      integer, parameter :: N_land_classification = 24 ! GLCC USGS number of land cover category
#endif
#if(defined IGBP_CLASSIFICATION)
      integer, parameter :: N_land_classification = 17 ! MODIS IGBP number of land cover category
#endif
      integer, parameter :: nlat = 21600               ! 180*(60*2)
      integer, parameter :: nlon = 43200               ! 360*(60*2)

      character(LEN=256), intent(in) :: dir_rawdata
      character(LEN=256), intent(in) :: dir_model_landdata

      integer,  intent(in) :: lon_points       ! number of model longitude grid cells
      integer,  intent(in) :: lat_points       ! number of model latitude grid cells
      integer,  intent(in) :: nrow_start       ! start latitudinal point in rawdata files for the model domain
      integer,  intent(in) :: nrow_end         ! end latitudinal point in rawdata files for the model domain
      integer,  intent(in) :: ncol_start       ! start longitudinal point in rawdata files for the model domain
      integer,  intent(in) :: ncol_end         ! end longitudinal point in rawdata files for the model domain
      integer,  intent(in) :: nx_fine_gridcell ! maximum of number of fine longitude cells for a model grid cell 
      integer,  intent(in) :: ny_fine_gridcell ! maximum of number of fine latitude cells for a model grid cell

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
! ---------------------------------------------------------------
      character(len=256) lndname
      character(len=1) land_chr1(nlon)
      character(len=2) land_chr2(nlon)
      integer(kind=1)  land_int1(nlon)
      integer(kind=2)  land_int2(nlon)
      real(r8) tmp(nlon)

!----------weinan------------
#if(defined Naqu  || defined Gravels0)
      real(r8) soil_grav_l(8)  ! gravel content   (% of volume)
      real(r8) soil_sand_l(8)  ! sand percentage  (% of weight)
      real(r8) soil_clay_l(8)  ! clay percentage  (% of weight)
      real(r8) soil_oc_l(8)    ! organic carbon percentage (% of weight)
      real(r8) soil_bd_l(8)    ! bulk density     (g/cm3)

      real(r8) vf_quartz_mineral_s ! volumetric fraction of quartz within mineral soil
      real(r8) wf_gravels_s ! weight fraction of gravel
      real(r8) wf_om_s      ! weight fraction of organic matter
      real(r8) wf_sand_s    ! weight fraction of sand
      real(r8) wf_clay_s    ! weight fraction of clay

      real(r8) vf_gravels_s ! volumetric fraction of gravels
      real(r8) vf_om_s      ! volumetric fraction of organic matter
      real(r8) vf_sand_s    ! volumetric fraction of sand
      real(r8) vf_clay_s    ! volumetric fraction of clay
      real(r8) vf_silt_s    ! volumetric fraction of silt
      real(r8) vf_pores_s   ! volumetric pore space of the soil

      real(r8) BD_mineral_s ! bulk density of mineral soil (g/cm^3)

      real(r8) theta_s      ! saturated water content (cm3/cm3)
      real(r8) psi_s        ! matric potential at saturation (cm)
      real(r8) lambda       ! pore size distribution index (dimensionless)
      real(r8) k_s          ! saturated hydraulic conductivity (cm/day)

      real(r8) VGM_theta_r  ! residual moisture content
      real(r8) VGM_alpha    ! a parameter corresponding approximately to the inverse of the air-entry value
      real(r8) VGM_n        ! a shape parameter
      real(r8) VGM_L        ! pore-connectivity parameter

      real(r8) csol         ! heat capacity of soil solids [J/(m3 K)]
      real(r8) k_solids     ! thermal conductivity of soil solids [W/m/K]
      real(r8) tksatu       ! thermal conductivity of unfrozen saturated soil [W/m-K]
      real(r8) tksatf       ! thermal conductivity of frozen saturated soil [W/m-K]
      real(r8) tkdry        ! thermal conductivity for dry soil  [W/(m-K)]
      real(r8) soildepth
      real(r8) a,SOM

      integer nl_soil
      real(r8), allocatable ::  zsoi(:)  ! soil layer depth [m]
      real(r8), allocatable ::  dzsoi(:) ! soil node thickness [m]
      real(r8), allocatable ::  zsoih(:) ! interface level below a zsoi level [m]

      integer(kind=1), allocatable :: int_soil_sand_l (:,:) ! Sand content (50-2000 micro meter) mass fraction in %
      integer(kind=1), allocatable :: int_soil_clay_l (:,:) ! Clay content (0-2 micro meter) mass fraction in %
      integer(kind=2), allocatable :: int_soil_oc_l   (:,:) ! Soil organic carbon content (fine earth fraction) in g/kg
      integer(kind=2), allocatable :: int_soil_bd_l   (:,:) ! Bulk density (fine earth) in kg/m3

! soil hydraulic parameters from Rosetta3-H3w
      real(r8), allocatable :: VGM_theta_r_Rose (:,:) ! residual moisture content
      real(r8), allocatable :: VGM_alpha_Rose   (:,:) ! a parameter corresponding approximately to the inverse of the air-entry value
      real(r8), allocatable :: VGM_n_Rose       (:,:) ! a shape parameter
      real(r8), allocatable :: k_s_Rose         (:,:) ! saturated hydraulic conductivity (cm/day)
#endif
!----------------------------

      character(len=256) c
      integer iunit
      integer length
      integer i, j, L
      integer i1, i2, j1, j2, isfirst
      integer nrow, ncol, ncol_mod
      integer LL, np
      integer n_fine_gridcell
      integer nsl, MODEL_SOIL_LAYER

      integer, allocatable :: landtypes (:,:)
      integer, allocatable :: num_patches(:) 

      real(r8), allocatable :: vf_quartz_mineral_s_l      (:,:)
      real(r8), allocatable :: vf_gravels_s_l             (:,:)
      real(r8), allocatable :: vf_om_s_l                  (:,:)
      real(r8), allocatable :: vf_sand_s_l                (:,:)
      real(r8), allocatable :: wf_gravels_s_l             (:,:)
      real(r8), allocatable :: wf_sand_s_l                (:,:)
      real(r8), allocatable :: theta_s_l                  (:,:)
      real(r8), allocatable :: psi_s_l                    (:,:)
      real(r8), allocatable :: lambda_l                   (:,:)
      real(r8), allocatable :: k_s_l                      (:,:)
      real(r8), allocatable :: csol_l                     (:,:)
      real(r8), allocatable :: k_solids_l                 (:,:)
      real(r8), allocatable :: tksatu_l                   (:,:)
      real(r8), allocatable :: tksatf_l                   (:,:)
      real(r8), allocatable :: tkdry_l                    (:,:)
      real(r8), allocatable :: VGM_theta_r_l              (:,:)
      real(r8), allocatable :: VGM_alpha_l                (:,:)
      real(r8), allocatable :: VGM_n_l                    (:,:)
      real(r8), allocatable :: VGM_L_l                    (:,:)

      real(r8), allocatable :: a_vf_quartz_mineral_s_l    (:,:)
      real(r8), allocatable :: a_vf_gravels_s_l           (:,:)
      real(r8), allocatable :: a_vf_om_s_l                (:,:)
      real(r8), allocatable :: a_vf_sand_s_l              (:,:)
      real(r8), allocatable :: a_wf_gravels_s_l           (:,:)
      real(r8), allocatable :: a_wf_sand_s_l              (:,:)     
      real(r8), allocatable :: a_theta_s_l                (:,:) 
      real(r8), allocatable :: a_psi_s_l                  (:,:) 
      real(r8), allocatable :: a_lambda_l                 (:,:) 
      real(r8), allocatable :: a_k_s_l                    (:,:) 
      real(r8), allocatable :: a_csol_l                   (:,:) 
      real(r8), allocatable :: a_k_solids_l               (:,:)
      real(r8), allocatable :: a_tksatu_l                 (:,:) 
      real(r8), allocatable :: a_tksatf_l                 (:,:)
      real(r8), allocatable :: a_tkdry_l                  (:,:)
      real(r8), allocatable :: a_VGM_theta_r_l            (:,:)
      real(r8), allocatable :: a_VGM_alpha_l              (:,:)
      real(r8), allocatable :: a_VGM_n_l                  (:,:)
      real(r8), allocatable :: a_VGM_L_l                  (:,:)
      real(r8), allocatable :: a_BA_alpha_l               (:,:)
      real(r8), allocatable :: a_BA_beta_l                (:,:)
      real(r8), allocatable :: a_area_fine_gridcell       (:,:)

      real(r8), allocatable :: soil_vf_quartz_mineral_s_l (:,:,:)
      real(r8), allocatable :: soil_vf_gravels_s_l        (:,:,:)
      real(r8), allocatable :: soil_vf_om_s_l             (:,:,:)
      real(r8), allocatable :: soil_vf_sand_s_l           (:,:,:)
      real(r8), allocatable :: soil_wf_gravels_s_l        (:,:,:)
      real(r8), allocatable :: soil_wf_sand_s_l           (:,:,:)
      real(r8), allocatable :: soil_theta_s_l             (:,:,:) 
      real(r8), allocatable :: soil_psi_s_l               (:,:,:) 
      real(r8), allocatable :: soil_lambda_l              (:,:,:) 
      real(r8), allocatable :: soil_k_s_l                 (:,:,:) 
      real(r8), allocatable :: soil_csol_l                (:,:,:) 
      real(r8), allocatable :: soil_k_solids_l            (:,:,:)
      real(r8), allocatable :: soil_tksatu_l              (:,:,:) 
      real(r8), allocatable :: soil_tksatf_l              (:,:,:)
      real(r8), allocatable :: soil_tkdry_l               (:,:,:)
      real(r8), allocatable :: soil_VGM_theta_r_l         (:,:,:)
      real(r8), allocatable :: soil_VGM_alpha_l           (:,:,:)
      real(r8), allocatable :: soil_VGM_n_l               (:,:,:)
      real(r8), allocatable :: soil_VGM_L_l               (:,:,:)
      real(r8), allocatable :: soil_BA_alpha_l            (:,:,:)
      real(r8), allocatable :: soil_BA_beta_l             (:,:,:) 

      real(r8), external :: median
      real(r8), external :: find_min_area

#if(defined SOILPAR_UPS_FIT)
! local variables for estimating the upscaled soil parameters using the Levenberg–Marquardt fitting method
! ---------------------------------------------------------------
      integer, parameter   :: npointw  = 24      
      integer, parameter   :: npointb  = 20
      real(r8),parameter   :: xdat(npointw) = (/1.,5.,10.,20.,30.,40.,50.,60.,70.,90.,110.,130.,150.,&
                                   170.,210.,300.,345.,690.,1020.,5100.,15300.,20000.,100000.,1000000./)
                              !  points of soil pressure heads used for fitting SW retention curves

      real(r8),parameter   :: xdatsr(npointb)=(/0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,&
                                                0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0/)
                              !  points of soil saturation levels (Sr) used for fitting Ke-Sr relationship

      real(r8),allocatable :: ydatc(:,:)     ! the Campbell SW retentions at fine grids 
      real(r8),allocatable :: ydatv(:,:)     ! the van Genuchten SW retentions at fine grids
      real(r8),allocatable :: ydatb(:,:)     ! the Balland and Arp (2005) Ke-Sr relationship at fine grids

      integer, parameter   :: nc = 2         ! the number of fitted parameters in Campbell SW retention curve (psi and lambda)
      integer, parameter   :: nv = 3         ! the number of fitted parameters in van Genuchten SW retention curve 
                                             ! (theta_r, alpha and n)
      integer, parameter   :: nb = 2         ! the number of fitted parameters in Ke-Sr relationship (alpha and beta)

! Variables needed for Levenberg–Marquardt algorithm in MINPACK library      
      real(r8),parameter   :: factor = 0.01
      real(r8),parameter   :: ftol = 1.0e-3
      real(r8),parameter   :: xtol = 1.0e-4
      real(r8),parameter   :: gtol = 0.0
      integer, parameter   :: mode = 1
      integer, parameter   :: nprint = 0
      integer              :: ldfjac,info,ipvtc(nc),ipvtv(nv),ipvtb(nb),maxfev,nfev,njev
      real(r8)             :: xc(nc),xv(nv),xb(nb),diagc(nc),diagv(nv),diagb(nb),qtfc(nc),qtfv(nv),qtfb(nb)
      real(r8),allocatable :: fjacc(:,:),fvecc(:),fjacv(:,:),fvecv(:),fjacb(:,:),fvecb(:)
      integer isiter                         ! flags to tell whether the iteration is completed, 1=Yes, 0=No

      external SW_CB_dist                    ! the objective function to be fitted for Campbell SW retention curve
      external SW_VG_dist                    ! the objective function to be fitted for van Genuchten SW retention curve
      external Ke_Sr_dist                    ! the objective function to be fitted for Balland and Arp (2005) Ke-Sr relationship
#endif

! ........................................
! ... [1] gloabl land cover types
! ........................................
      iunit = 100
      inquire(iolength=length) land_chr1
      allocate (landtypes(nlon,nlat))

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
         do ncol = 1, nlon
            landtypes(ncol,nrow) = ichar(land_chr1(ncol))
         enddo
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
         do ncol = 1, nlon
            landtypes(ncol,nrow) = ichar(land_chr1(ncol))
         enddo
      enddo
      close (iunit)
#endif

#endif

! ........................................
! ... [2] aggregate the soil parameters from the resolution of raw data to modelling resolution
! ........................................
      n_fine_gridcell = nx_fine_gridcell * ny_fine_gridcell

      allocate ( num_patches (0:N_land_classification) )

      allocate ( vf_quartz_mineral_s_l     (nlon,nlat) )
      allocate ( vf_gravels_s_l            (nlon,nlat) )
      allocate ( vf_om_s_l                 (nlon,nlat) )
      allocate ( vf_sand_s_l               (nlon,nlat) )
      allocate ( wf_gravels_s_l            (nlon,nlat) )
      allocate ( wf_sand_s_l               (nlon,nlat) )
      allocate ( theta_s_l                 (nlon,nlat) )
      allocate ( psi_s_l                   (nlon,nlat) )
      allocate ( lambda_l                  (nlon,nlat) )
      allocate ( k_s_l                     (nlon,nlat) )
      allocate ( csol_l                    (nlon,nlat) )
      allocate ( k_solids_l                (nlon,nlat) )
      allocate ( tksatu_l                  (nlon,nlat) )
      allocate ( tksatf_l                  (nlon,nlat) )
      allocate ( tkdry_l                   (nlon,nlat) )
      allocate ( VGM_theta_r_l             (nlon,nlat) )
      allocate ( VGM_alpha_l               (nlon,nlat) )
      allocate ( VGM_n_l                   (nlon,nlat) )
      allocate ( VGM_L_l                   (nlon,nlat) )

      allocate ( a_vf_quartz_mineral_s_l (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_vf_gravels_s_l        (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_vf_om_s_l             (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_vf_sand_s_l           (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_wf_gravels_s_l        (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_wf_sand_s_l           (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_theta_s_l             (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_psi_s_l               (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_lambda_l              (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_k_s_l                 (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_csol_l                (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_k_solids_l            (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_tksatu_l              (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_tksatf_l              (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_tkdry_l               (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_VGM_theta_r_l         (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_VGM_alpha_l           (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_VGM_n_l               (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_VGM_L_l               (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_BA_alpha_l            (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_BA_beta_l             (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_area_fine_gridcell    (0:N_land_classification,1:n_fine_gridcell) )
#if(defined SOILPAR_UPS_FIT)
      allocate ( ydatc                   (1:n_fine_gridcell,npointw) )
      allocate ( ydatv                   (1:n_fine_gridcell,npointw) )
      allocate ( ydatb                   (1:n_fine_gridcell,npointb) )
! the jacobian matrix required in Levenberg–Marquardt fitting method
      allocate ( fjacc                   (1:n_fine_gridcell,     nc) )           ! calculated in SW_CB_dist
      allocate ( fjacv                   (1:n_fine_gridcell,     nv) )           ! calculated in SW_VG_dist
      allocate ( fjacb                   (1:n_fine_gridcell,     nb) )           ! calculated in Ke_Sr_dist
! the values of objective functions to be fitted 
      allocate ( fvecc                   (1:n_fine_gridcell)         )           ! calculated in SW_CB_dist
      allocate ( fvecv                   (1:n_fine_gridcell)         )           ! calculated in SW_VG_dist
      allocate ( fvecb                   (1:n_fine_gridcell)         )           ! calculated in Ke_Sr_dist
#endif

      allocate ( soil_vf_quartz_mineral_s_l (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_vf_gravels_s_l        (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_vf_om_s_l             (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_vf_sand_s_l           (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_wf_gravels_s_l        (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_wf_sand_s_l           (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_theta_s_l             (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_psi_s_l               (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_lambda_l              (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_k_s_l                 (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_csol_l                (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_k_solids_l            (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_tksatu_l              (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_tksatf_l              (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_tkdry_l               (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_VGM_theta_r_l         (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_VGM_alpha_l           (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_VGM_n_l               (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_VGM_L_l               (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_BA_alpha_l            (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_BA_beta_l             (0:N_land_classification,1:lon_points,1:lat_points) )

!----------weinan------------
#if(defined Naqu || defined Gravels0)
      allocate ( int_soil_sand_l       (nlon,nlat) ,&
                 int_soil_clay_l       (nlon,nlat) ,&
                 int_soil_oc_l         (nlon,nlat) ,&
                 int_soil_bd_l         (nlon,nlat)  )

      allocate ( VGM_theta_r_Rose      (nlon,nlat) ,&
                 VGM_alpha_Rose        (nlon,nlat) ,&
                 VGM_n_Rose            (nlon,nlat) ,&
                 k_s_Rose              (nlon,nlat)  )

      nl_soil = 10
      allocate ( zsoi(1:nl_soil), dzsoi(1:nl_soil), zsoih(0:nl_soil) )

      ! ----------------------------------
      ! soil layer thickness, depths (m)
      ! ----------------------------------
      do nsl = 1, nl_soil
        zsoi(nsl) = 0.025*(exp(0.5*(nsl-0.5))-1.)  ! node depths
      end do

      dzsoi(1) = 0.5*(zsoi(1)+zsoi(2))         ! =zsoih(1)
      dzsoi(nl_soil) = zsoi(nl_soil)-zsoi(nl_soil-1)
      do nsl = 2, nl_soil-1
         dzsoi(nsl) = 0.5*(zsoi(nsl+1)-zsoi(nsl-1))  ! thickness b/n two interfaces
      end do

      zsoih(0) = 0.
      zsoih(nl_soil) = zsoi(nl_soil) + 0.5*dzsoi(nl_soil)
      do nsl = 1, nl_soil-1
         zsoih(nsl) = 0.5*(zsoi(nsl)+zsoi(nsl+1))    ! interface depths
      enddo
#endif 
!---------------------------- 

      iunit = 100
      DO nsl = 1, 8
         MODEL_SOIL_LAYER = nsl
         write(c,'(i1)') MODEL_SOIL_LAYER

! (1) Read in the volumetric fraction of quartz within mineral soil
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/vf_quartz_mineral_s_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) vf_quartz_mineral_s_l(:,nrow)
         enddo
         close(iunit)

! (2) Read in the volumetric fraction of gravels
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/vf_gravels_s_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) vf_gravels_s_l(:,nrow)
         enddo
         close(iunit)

! (3) Read in the volumetric fraction of organic matter
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/vf_om_s_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) vf_om_s_l(:,nrow)
         enddo
         close(iunit)

! (4) Read in volumetric fraction of sand
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/vf_sand_s_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) vf_sand_s_l(:,nrow)
         enddo
         close(iunit)

! (5) Read in the gravimetric fraction of gravels
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/wf_gravels_s_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) wf_gravels_s_l(:,nrow)
         enddo
         close(iunit)

! (6) Read in gravimetric fraction of sand
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/wf_sand_s_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) wf_sand_s_l(:,nrow) 
         enddo
         close(iunit)

! (7) Read in the saturated water content [cm3/cm3]
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/theta_s_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) theta_s_l(:,nrow)
         enddo
         close(iunit)

! (8) Read in the matric potential at saturation [cm]
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/psi_s_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) psi_s_l(:,nrow)
         enddo
         close(iunit)

! (9) Read in the pore size distribution index [dimensionless]
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/lambda_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) lambda_l(:,nrow)
         enddo
         close(iunit)

! (10) Read in the saturated hydraulic conductivity [cm/day]
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/k_s_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) k_s_l(:,nrow)
         enddo
         close(iunit)

! (11) Read in the heat capacity of soil solids [J/(m3 K)]
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/csol_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) csol_l(:,nrow)
         enddo
         close(iunit)

! (12) Read in the thermal conductivity of soil solids [W/m-K]
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/k_solids_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) k_solids_l(:,nrow)
         enddo
         close(iunit)

! (13) Read in the thermal conductivity of unfrozen saturated soil [W/m-K]
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/tksatu_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) tksatu_l(:,nrow)
         enddo
         close(iunit)

! (14) Read in the thermal conductivity of frozen saturated soil [W/m-K]
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/tksatf_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) tksatf_l(:,nrow)
         enddo
         close(iunit)

! (15) Read in the thermal conductivity for dry soil [W/(m-K)]
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/tkdry_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) tkdry_l(:,nrow)
         enddo
         close(iunit)

! (16) Read in the VGM's residual moisture content
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/VGM_theta_r_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) VGM_theta_r_l(:,nrow)
         enddo
         close(iunit)

! (17) Read in the VGM's parameter corresponding approximately to the inverse of the air-entry value
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/VGM_alpha_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) VGM_alpha_l(:,nrow)
         enddo
         close(iunit)

! (18) Read in the VGM's shape parameter
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/VGM_n_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) VGM_n_l(:,nrow)
         enddo
         close(iunit)

! (19) Read in the VGM's pore-connectivity parameter
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated/VGM_L_l'//trim(c)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) VGM_L_l(:,nrow)
         enddo
         close(iunit)

!-------------------weinan---------------------
#if(defined Gravels0)
! (20) Read in percentage of sand (% weight)
         inquire(iolength=length) land_int1
         lndname = trim(dir_rawdata)//'soil/SAND_L'//trim(c)
         print*,lndname

         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) int_soil_sand_l(:,nrow)
         enddo
         close (iunit)

! (21) Read in percentage of clay (% weight)
         inquire(iolength=length) land_int1
         lndname = trim(dir_rawdata)//'soil/CLAY_L'//trim(c)
         print*,lndname

         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) int_soil_clay_l(:,nrow)
         enddo
         close (iunit)

! (22) Read in percentage of organic carbon (% weight)
         inquire(iolength=length) land_int2
         lndname = trim(dir_rawdata)//'soil/OC_L'//trim(c)
         print*,lndname

         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) int_soil_oc_l(:,nrow)
         enddo
         close (iunit)

! (23) Read in bulk density (g/cm3) 
         inquire(iolength=length) land_int2
         lndname = trim(dir_rawdata)//'soil/BD_L'//trim(c)
         print*,lndname

         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) int_soil_bd_l(:,nrow)
         enddo
         close (iunit)

! (24) Rosetta parameters generated by yonggen Zhang
         inquire(iolength=length) VGM_theta_r_Rose(:,1)
         lndname = '/work/ygzhang/data/CLMrawdata_2021/Rosetta_VGM/SSCBD_L'//trim(c)//'_VGM_merged_thr_binary'
         print*,lndname

         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) VGM_theta_r_Rose(:,nrow)
         end do
         close(iunit)

         inquire(iolength=length) VGM_alpha_Rose(:,1)
         lndname = '/work/ygzhang/data/CLMrawdata_2021/Rosetta_VGM/SSCBD_L'//trim(c)//'_VGM_merged_alpha_binary'
         print*,lndname

         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) VGM_alpha_Rose(:,nrow)
         end do
         close(iunit)

         inquire(iolength=length) VGM_n_Rose(:,1)
         lndname = '/work/ygzhang/data/CLMrawdata_2021/Rosetta_VGM/SSCBD_L'//trim(c)//'_VGM_merged_n_binary'
         print*,lndname

         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) VGM_n_Rose(:,nrow)
         end do
         close(iunit)

         inquire(iolength=length) k_s_Rose(:,1)
         lndname = '/work/ygzhang/data/CLMrawdata_2021/Rosetta_VGM/SSCBD_L'//trim(c)//'_VGM_merged_Ks_binary'
         print*,lndname

         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) k_s_Rose(:,nrow)
         end do
         close(iunit)

#endif
!----------------------------------------------------------------

#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,i1,i2,j1,j2,nrow,ncol,ncol_mod,L,LL,num_patches,np) &
!$OMP PRIVATE(a_vf_quartz_mineral_s_l,a_vf_gravels_s_l,a_vf_om_s_l,a_vf_sand_s_l,a_wf_gravels_s_l,a_wf_sand_s_l) &
#if(defined Naqu  || defined Gravels0)
!$OMP PRIVATE(soil_grav_l,soil_sand_l,soil_clay_l,soil_oc_l,soil_bd_l) &
!$OMP PRIVATE(vf_quartz_mineral_s,wf_gravels_s,wf_om_s,wf_sand_s,wf_clay_s) &
!$OMP PRIVATE(vf_gravels_s,vf_om_s,vf_sand_s,vf_clay_s,vf_silt_s,vf_pores_s) &
!$OMP PRIVATE(BD_mineral_s,theta_s,psi_s,lambda,k_s,csol,k_solids,tksatu,tksatf,tkdry) &
!$OMP PRIVATE(VGM_theta_r,VGM_alpha,VGM_n,VGM_L) &
!$OMP PRIVATE(soildepth,a,SOM) &
#endif
#if(defined SOILPAR_UPS_FIT)
!$OMP PRIVATE(ydatc,ydatv,ydatb,ldfjac,info,ipvtc,ipvtv,ipvtb,maxfev,nfev,njev) &
!$OMP PRIVATE(xc,xv,xb,diagc,diagv,diagb,qtfc,qtfv,qtfb,fjacc,fvecc,fjacv,fvecv,fjacb,fvecb,isiter) &
#endif
!$OMP PRIVATE(a_theta_s_l,a_psi_s_l,a_lambda_l,a_k_s_l,a_csol_l,a_k_solids_l,a_tksatu_l,a_tksatf_l,a_tkdry_l) &
!$OMP PRIVATE(a_VGM_theta_r_l,a_VGM_alpha_l,a_VGM_n_l,a_VGM_L_l,a_BA_alpha_l,a_BA_beta_l,a_area_fine_gridcell,isfirst)
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
                  num_patches(L) = num_patches(L) + 1
                  LL = num_patches(L) 

                  a_vf_quartz_mineral_s_l (L,LL) = vf_quartz_mineral_s_l (ncol_mod,nrow)
                  a_vf_gravels_s_l        (L,LL) = vf_gravels_s_l        (ncol_mod,nrow)
                  a_vf_om_s_l             (L,LL) = vf_om_s_l             (ncol_mod,nrow)
                  a_vf_sand_s_l           (L,LL) = vf_sand_s_l           (ncol_mod,nrow)
                  a_wf_gravels_s_l        (L,LL) = wf_gravels_s_l        (ncol_mod,nrow)
                  a_wf_sand_s_l           (L,LL) = wf_sand_s_l           (ncol_mod,nrow)
                  
                  a_theta_s_l             (L,LL) = theta_s_l             (ncol_mod,nrow) 
                  a_psi_s_l               (L,LL) = psi_s_l               (ncol_mod,nrow)
                  a_lambda_l              (L,LL) = lambda_l              (ncol_mod,nrow)
                  a_k_s_l                 (L,LL) = k_s_l                 (ncol_mod,nrow)
                  a_csol_l                (L,LL) = csol_l                (ncol_mod,nrow)
                  a_k_solids_l            (L,LL) = k_solids_l            (ncol_mod,nrow)
                  a_tksatu_l              (L,LL) = tksatu_l              (ncol_mod,nrow)
                  a_tksatf_l              (L,LL) = tksatf_l              (ncol_mod,nrow)
                  a_tkdry_l               (L,LL) = tkdry_l               (ncol_mod,nrow)
                  a_VGM_theta_r_l         (L,LL) = VGM_theta_r_l         (ncol_mod,nrow)
                  a_VGM_alpha_l           (L,LL) = VGM_alpha_l           (ncol_mod,nrow)
                  a_VGM_n_l               (L,LL) = VGM_n_l               (ncol_mod,nrow)
                  a_VGM_L_l               (L,LL) = VGM_L_l               (ncol_mod,nrow)

                  ! the parameter values of Balland and Arp (2005) Ke-Sr relationship, 
                  ! modified by Barry-Macaulay et al.(2015), Evaluation of soil thermal conductivity models
                  if ((a_vf_gravels_s_l(L,LL) + a_vf_sand_s_l(L,LL)) > 0.4)  then
                     a_BA_alpha_l         (L,LL) = 0.38
                     a_BA_beta_l          (L,LL) = 35.0
                  else if ((a_vf_gravels_s_l(L,LL) + a_vf_sand_s_l(L,LL)) > 0.25)  then
                     a_BA_alpha_l         (L,LL) = 0.24
                     a_BA_beta_l          (L,LL) = 26.0
                  else
                     a_BA_alpha_l         (L,LL) = 0.2
                     a_BA_beta_l          (L,LL) = 10.0
                  end if
                  a_area_fine_gridcell    (L,LL) = find_min_area(lone_rad(i),lonw_rad(i),lone_rad_i(ncol_mod),&
                                       lonw_rad_i(ncol_mod),sinn(j),sins(j),sinn_i(nrow),sins_i(nrow),isfirst)
                  isfirst = 0

!----------weinan------------
                if(L/=0)then  ! NOT OCEAN(0)
#if(defined Gravels0)
                  soil_grav_l(nsl) = 0.0 ! Volumetric in %
                  soil_sand_l(nsl) = int_soil_sand_l(ncol_mod,nrow) ! Gravimetric fraction in %
                  soil_clay_l(nsl) = int_soil_clay_l(ncol_mod,nrow) ! Gravimetric fraction in %
                  soil_oc_l(nsl)   = int_soil_oc_l  (ncol_mod,nrow)*0.01 ! Gravimetric fraction (fine earth) in %
                  soil_bd_l(nsl)   = int_soil_bd_l  (ncol_mod,nrow)*0.01 ! (fine earth) in g/cm3

#if(defined USGS_CLASSIFICATION)
                  if(landtypes(ncol_mod,nrow)==16)then   !WATER BODIES(16)
#endif
#if(defined IGBP_CLASSIFICATION)
                  if(landtypes(ncol_mod,nrow)==17)then   !WATER BODIES(17)
#endif
                     soil_grav_l(nsl) = 0.
                     soil_sand_l(nsl) = 10.
                     soil_clay_l(nsl) = 45.
                     soil_oc_l(nsl)   = 3.0
                     soil_bd_l(nsl)   = 1.2
                  endif

#if(defined USGS_CLASSIFICATION)
                  if(landtypes(ncol_mod,nrow)==24)then   !GLACIER and ICESHEET(24)
#endif
#if(defined IGBP_CLASSIFICATION)
                  if(landtypes(ncol_mod,nrow)==15)then   !GLACIER and ICE SHEET(15)
#endif
                     soil_grav_l(nsl) = 0.
                     soil_sand_l(nsl) = 89.
                     soil_clay_l(nsl) = 1.
                     soil_oc_l(nsl)   = 0.
                     soil_bd_l(nsl)   = 2.0
                  endif

                  if( soil_sand_l(nsl) < 0.0 ) soil_sand_l(nsl) = 43.   ! missing value = -100
                  if( soil_clay_l(nsl) < 0.0 ) soil_clay_l(nsl) = 18.   ! missing value = -100
                  if( soil_oc_l(nsl)   < 0.0   ) soil_oc_l(nsl) = 1.0     ! missing value = -999
                  if( soil_bd_l(nsl)   < 0.0   ) soil_bd_l(nsl) = 1.2     ! missing value = -999

                  if( soil_sand_l(nsl) < 1.0   ) soil_sand_l(nsl) = 1.
                  if( soil_clay_l(nsl) < 1.0   ) soil_clay_l(nsl) = 1.

                  a = soil_sand_l(nsl) + soil_clay_l(nsl)
                  if( a >= 96. ) then
                      soil_sand_l(nsl) = soil_sand_l(nsl) * 96. / a
                      soil_clay_l(nsl) = soil_clay_l(nsl) * 96. / a
                  endif
                  if( soil_oc_l(nsl) < 0.01 ) soil_oc_l(nsl) = 0.01
                  if( soil_oc_l(nsl) > 58.0 ) soil_oc_l(nsl) = 58.0
                  if( soil_bd_l(nsl) < 0.1  ) soil_bd_l(nsl) = 0.1

                  if (soil_bd_l(nsl) < 0.111 .or. soil_bd_l(nsl) > 2.0 .or. soil_oc_l(nsl) > 10.0) then
                     SOM=1.724*soil_oc_l(nsl)
                     soil_bd_l(nsl) = 0.111*2.0/(2.0*SOM/100.+0.111*(100.-SOM)/100.)
                  end if
#endif
#if(defined Naqu)
                  soil_grav_l(1:8) = (/10.0,10.0,10.0,10.0,19.03,28.46,28.46,28.46/)
                  soil_sand_l(1:8) = (/63.68,63.68,63.68,43.5,71.94,67.08,64.75,64.75/)
                  soil_clay_l(1:8) = (/4.13,4.13,4.13,10.99,3.58,0.88,1.87,1.87/)
                  soil_oc_l(1:8)   = (/100.4,100.4,70.53,45.15,24.91,13.52,3.0,0.0/)  !soil organic matters(kg/m3)

                  soil_bd_l(nsl) = max(2.7*(1-a_theta_s_l (L,LL)),0.1)
                  soil_oc_l(nsl) = min(soil_oc_l(nsl)/1000.0/soil_bd_l(nsl)/1.724*100.0,58.0)

                  if( soil_sand_l(nsl) < 1.0   ) soil_sand_l(nsl) = 1.
                  if( soil_clay_l(nsl) < 1.0   ) soil_clay_l(nsl) = 1.
                  if( soil_oc_l(nsl)   < 0.01  ) soil_oc_l(nsl)   = 0.01

                  if (soil_bd_l(nsl) < 0.111 .or. soil_bd_l(nsl) > 2.0 .or. soil_oc_l(nsl) > 10.0) then
                     SOM=1.724*soil_oc_l(nsl)
                     soil_bd_l(nsl) = 0.111*2.0/(2.0*SOM/100.+0.111*(100.-SOM)/100.)
                  end if
#if(defined Gravels0)
                  soil_grav_l(nsl) = 0.0
#endif
#endif
#if(defined Naqu || defined Gravels0)
                  soildepth = zsoih(MODEL_SOIL_LAYER + 1)*100.0

                 ! --------------------------------------------------
                 ! The weight and volumetric fractions of soil solids
                 ! --------------------------------------------------
                  CALL soil_solids_fractions(&
                       soil_bd_l(nsl),soil_grav_l(nsl),soil_oc_l(nsl),soil_sand_l(nsl),soil_clay_l(nsl),&
                       wf_gravels_s,wf_om_s,wf_sand_s,wf_clay_s,&
                       vf_gravels_s,vf_om_s,vf_sand_s,vf_clay_s,vf_silt_s,vf_pores_s,&
                       vf_quartz_mineral_s,BD_mineral_s)

                       theta_s = vf_pores_s

                 ! ---------------------------------------------------------------------
                 ! The volumetric heat capacity and  thermal conductivity of soil solids
                 ! ---------------------------------------------------------------------
                  CALL soil_thermal_parameters(&
                       wf_gravels_s,wf_sand_s,wf_clay_s,&
                       vf_gravels_s,vf_om_s,vf_sand_s,vf_clay_s,vf_silt_s,vf_pores_s,&
                       vf_quartz_mineral_s,BD_mineral_s,k_solids,&
                       csol,tkdry,tksatu,tksatf)

                 ! -----------------------------
                 ! The soil hydraulic properties
                 ! -----------------------------
                  CALL soil_hydraulic_parameters(soil_bd_l(nsl),soil_sand_l(nsl),soil_clay_l(nsl),soil_oc_l(nsl),soildepth,&
                       vf_gravels_s,theta_s,psi_s,lambda,k_s,&
                       VGM_theta_r,VGM_alpha,VGM_n,VGM_L,&
                       VGM_theta_r_Rose(ncol_mod,nrow),VGM_alpha_Rose(ncol_mod,nrow),VGM_n_Rose(ncol_mod,nrow),&
                       k_s_Rose(ncol_mod,nrow))

                  vf_gravels_s = vf_gravels_s/(1 - vf_pores_s)
                  vf_om_s      = vf_om_s     /(1 - vf_pores_s)
                  vf_sand_s    = vf_sand_s   /(1 - vf_pores_s)
                  wf_sand_s    = soil_sand_l(nsl) / 100.0
    
                  a_vf_quartz_mineral_s_l (L,LL) = vf_quartz_mineral_s
                  a_vf_gravels_s_l        (L,LL) = vf_gravels_s
                  a_vf_om_s_l             (L,LL) = vf_om_s
                  a_vf_sand_s_l           (L,LL) = vf_sand_s
                  a_wf_gravels_s_l        (L,LL) = wf_gravels_s
                  a_wf_sand_s_l           (L,LL) = wf_sand_s

                  a_theta_s_l             (L,LL) = theta_s
                  a_psi_s_l               (L,LL) = psi_s
                  a_lambda_l              (L,LL) = lambda
                  a_k_s_l                 (L,LL) = k_s
                  a_csol_l                (L,LL) = csol
                  a_k_solids_l            (L,LL) = k_solids
                  a_tksatu_l              (L,LL) = tksatu
                  a_tksatf_l              (L,LL) = tksatf
                  a_tkdry_l               (L,LL) = tkdry
                  a_VGM_theta_r_l         (L,LL) = VGM_theta_r
                  a_VGM_alpha_l           (L,LL) = VGM_alpha
                  a_VGM_n_l               (L,LL) = VGM_n
                  a_VGM_L_l               (L,LL) = VGM_L
                  if ((a_vf_gravels_s_l(L,LL) + a_vf_sand_s_l(L,LL)) > 0.4)  then
                     a_BA_alpha_l         (L,LL) = 0.38
                     a_BA_beta_l          (L,LL) = 35.0
                  else if ((a_vf_gravels_s_l(L,LL) + a_vf_sand_s_l(L,LL)) > 0.25)  then
                     a_BA_alpha_l         (L,LL) = 0.24
                     a_BA_beta_l          (L,LL) = 26.0
                  else
                     a_BA_alpha_l         (L,LL) = 0.2
                     a_BA_beta_l          (L,LL) = 10.0
                  end if

#endif          
                 end if    
!----------------------------  
               enddo
            enddo
               
               do L = 0, N_land_classification 
!#if(defined USGS_CLASSIFICATION)
!                 if(L/=0 .and. L/=16 .and. L/=24)then  ! NOT OCEAN(0)/WATER BODIES(16)/GLACIER and ICESHEET(24)
!#endif
!#if(defined IGBP_CLASSIFICATION)
!                 if(L/=0 .and. L/=17 .and. L/=15)then  ! NOT OCEAN(0)/WATER BODIES(17)/GLACIER and ICE SHEET(15)
!#endif
                  if(L/=0)then  ! NOT OCEAN(0)
                     np = num_patches(L) 
                     if(np == 0)then
                        soil_vf_quartz_mineral_s_l (L,i,j) = -1.0e36
                        soil_vf_gravels_s_l        (L,i,j) = -1.0e36
                        soil_vf_om_s_l             (L,i,j) = -1.0e36
                        soil_vf_sand_s_l           (L,i,j) = -1.0e36
                        soil_wf_gravels_s_l        (L,i,j) = -1.0e36
                        soil_wf_sand_s_l           (L,i,j) = -1.0e36
                        soil_theta_s_l             (L,i,j) = -1.0e36
                        soil_psi_s_l               (L,i,j) = -1.0e36
                        soil_lambda_l              (L,i,j) = -1.0e36
                        soil_k_s_l                 (L,i,j) = -1.0e36
                        soil_csol_l                (L,i,j) = -1.0e36
                        soil_k_solids_l            (L,i,j) = -1.0e36
                        soil_tksatu_l              (L,i,j) = -1.0e36
                        soil_tksatf_l              (L,i,j) = -1.0e36
                        soil_tkdry_l               (L,i,j) = -1.0e36
                        soil_VGM_theta_r_l         (L,i,j) = -1.0e36
                        soil_VGM_alpha_l           (L,i,j) = -1.0e36
                        soil_VGM_n_l               (L,i,j) = -1.0e36
                        soil_VGM_L_l               (L,i,j) = -1.0e36
                        soil_BA_alpha_l            (L,i,j) = -1.0e36
                        soil_BA_beta_l             (L,i,j) = -1.0e36
                     else if(np == 1) then
                        soil_vf_quartz_mineral_s_l (L,i,j) = a_vf_quartz_mineral_s_l (L,1)
                        soil_vf_gravels_s_l        (L,i,j) = a_vf_gravels_s_l        (L,1)
                        soil_vf_om_s_l             (L,i,j) = a_vf_om_s_l             (L,1)
                        soil_vf_sand_s_l           (L,i,j) = a_vf_sand_s_l           (L,1)
                        soil_wf_gravels_s_l        (L,i,j) = a_wf_gravels_s_l        (L,1)
                        soil_wf_sand_s_l           (L,i,j) = a_wf_sand_s_l           (L,1)
                        soil_theta_s_l             (L,i,j) = a_theta_s_l             (L,1)
                        soil_psi_s_l               (L,i,j) = a_psi_s_l               (L,1)
                        soil_lambda_l              (L,i,j) = a_lambda_l              (L,1)
                        soil_k_s_l                 (L,i,j) = a_k_s_l                 (L,1)
                        soil_csol_l                (L,i,j) = a_csol_l                (L,1)
                        soil_k_solids_l            (L,i,j) = a_k_solids_l            (L,1)
                        soil_tksatu_l              (L,i,j) = a_tksatu_l              (L,1)
                        soil_tksatf_l              (L,i,j) = a_tksatf_l              (L,1)
                        soil_tkdry_l               (L,i,j) = a_tkdry_l               (L,1)
                        soil_VGM_theta_r_l         (L,i,j) = a_VGM_theta_r_l         (L,1)
                        soil_VGM_alpha_l           (L,i,j) = a_VGM_alpha_l           (L,1)
                        soil_VGM_n_l               (L,i,j) = a_VGM_n_l               (L,1)
                        soil_VGM_L_l               (L,i,j) = a_VGM_L_l               (L,1) 
                        soil_BA_alpha_l            (L,i,j) = a_BA_alpha_l            (L,1)
                        soil_BA_beta_l             (L,i,j) = a_BA_beta_l             (L,1)
                     else
#if(defined SOILPAR_UPS_MEDIAN)
            soil_vf_quartz_mineral_s_l (L,i,j) = median ( a_vf_quartz_mineral_s_l (L,1:np), np)
            soil_vf_gravels_s_l        (L,i,j) = median ( a_vf_gravels_s_l        (L,1:np), np)
            soil_vf_om_s_l             (L,i,j) = median ( a_vf_om_s_l             (L,1:np), np)
            soil_vf_sand_s_l           (L,i,j) = median ( a_vf_sand_s_l           (L,1:np), np)
            soil_wf_gravels_s_l        (L,i,j) = median ( a_wf_gravels_s_l        (L,1:np), np)
            soil_wf_sand_s_l           (L,i,j) = median ( a_wf_sand_s_l           (L,1:np), np)
            soil_theta_s_l             (L,i,j) = median ( a_theta_s_l             (L,1:np), np)
            soil_psi_s_l               (L,i,j) = median ( a_psi_s_l               (L,1:np), np)
            soil_lambda_l              (L,i,j) = median ( a_lambda_l              (L,1:np), np)
            soil_k_s_l                 (L,i,j) = median ( a_k_s_l                 (L,1:np), np)
            soil_csol_l                (L,i,j) = median ( a_csol_l                (L,1:np), np)
            soil_k_solids_l            (L,i,j) = median ( a_k_solids_l            (L,1:np), np)
            soil_tksatu_l              (L,i,j) = median ( a_tksatu_l              (L,1:np), np)
            soil_tksatf_l              (L,i,j) = median ( a_tksatf_l              (L,1:np), np)
            soil_tkdry_l               (L,i,j) = median ( a_tkdry_l               (L,1:np), np)
            soil_VGM_theta_r_l         (L,i,j) = median ( a_VGM_theta_r_l         (L,1:np), np)
            soil_VGM_alpha_l           (L,i,j) = median ( a_VGM_alpha_l           (L,1:np), np)
            soil_VGM_n_l               (L,i,j) = median ( a_VGM_n_l               (L,1:np), np)
            soil_VGM_L_l               (L,i,j) = median ( a_VGM_L_l               (L,1:np), np)
            soil_BA_alpha_l            (L,i,j) = median ( a_BA_alpha_l            (L,1:np), np)
            soil_BA_beta_l             (L,i,j) = median ( a_BA_beta_l             (L,1:np), np)
#elif(defined SOILPAR_UPS_ARITHMETIC || defined SOILPAR_UPS_AAVW)
            a_area_fine_gridcell(L,1:np)       = a_area_fine_gridcell(L,1:np)/sum(a_area_fine_gridcell(L,1:np))
            soil_vf_quartz_mineral_s_l (L,i,j) = sum(a_vf_quartz_mineral_s_l(L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_vf_gravels_s_l        (L,i,j) = sum(a_vf_gravels_s_l       (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_vf_om_s_l             (L,i,j) = sum(a_vf_om_s_l            (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_vf_sand_s_l           (L,i,j) = sum(a_vf_sand_s_l          (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_wf_gravels_s_l        (L,i,j) = sum(a_wf_gravels_s_l       (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_wf_sand_s_l           (L,i,j) = sum(a_wf_sand_s_l          (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_theta_s_l             (L,i,j) = sum(a_theta_s_l            (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_psi_s_l               (L,i,j) = sum(a_psi_s_l              (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_lambda_l              (L,i,j) = sum(a_lambda_l             (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_k_s_l                 (L,i,j) = 10**(sum(log10(a_k_s_l    (L,1:np))*a_area_fine_gridcell(L,1:np)))
            soil_csol_l                (L,i,j) = sum(a_csol_l               (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_k_solids_l            (L,i,j) = sum(a_k_solids_l           (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_tksatu_l              (L,i,j) = sum(a_tksatu_l             (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_tksatf_l              (L,i,j) = sum(a_tksatf_l             (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_tkdry_l               (L,i,j) = sum(a_tkdry_l              (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_VGM_theta_r_l         (L,i,j) = sum(a_VGM_theta_r_l        (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_VGM_alpha_l           (L,i,j) = 10**(sum(log10(a_VGM_alpha_l(L,1:np))*a_area_fine_gridcell(L,1:np)))
            soil_VGM_n_l               (L,i,j) = 10**(sum(log10(a_VGM_n_l  (L,1:np))*a_area_fine_gridcell(L,1:np)))
            soil_VGM_L_l               (L,i,j) = sum(a_VGM_L_l              (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_BA_alpha_l            (L,i,j) = sum(a_BA_alpha_l           (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_BA_beta_l             (L,i,j) = sum(a_BA_beta_l            (L,1:np)*a_area_fine_gridcell(L,1:np))
#if(defined SOILPAR_UPS_AAVW)
              a_area_fine_gridcell(L,1:np)     =     1./(a_vf_quartz_mineral_s_l(L,1:np)-soil_vf_quartz_mineral_s_l(L,i,j))**2/&
                                                 sum(1./(a_vf_quartz_mineral_s_l(L,1:np)-soil_vf_quartz_mineral_s_l(L,i,j))**2)
            soil_vf_quartz_mineral_s_l (L,i,j) = sum(a_vf_quartz_mineral_s_l(L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(a_vf_gravels_s_l(L,1:np)-soil_vf_gravels_s_l(L,i,j))**2/&
                                                 sum(1./(a_vf_gravels_s_l(L,1:np)-soil_vf_gravels_s_l(L,i,j))**2)
            soil_vf_gravels_s_l        (L,i,j) = sum(a_vf_gravels_s_l       (L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(a_vf_om_s_l(L,1:np)-soil_vf_om_s_l(L,i,j))**2/&
                                                 sum(1./(a_vf_om_s_l(L,1:np)-soil_vf_om_s_l(L,i,j))**2)
            soil_vf_om_s_l             (L,i,j) = sum(a_vf_om_s_l            (L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(a_vf_sand_s_l(L,1:np)-soil_vf_sand_s_l(L,i,j))**2/&
                                                 sum(1./(a_vf_sand_s_l(L,1:np)-soil_vf_sand_s_l(L,i,j))**2)
            soil_vf_sand_s_l           (L,i,j) = sum(a_vf_sand_s_l          (L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(a_wf_gravels_s_l(L,1:np)-soil_wf_gravels_s_l(L,i,j))**2/&
                                                 sum(1./(a_wf_gravels_s_l(L,1:np)-soil_wf_gravels_s_l(L,i,j))**2)
            soil_wf_gravels_s_l        (L,i,j) = sum(a_wf_gravels_s_l       (L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(a_wf_sand_s_l(L,1:np)-soil_wf_sand_s_l(L,i,j))**2/&
                                                 sum(1./(a_wf_sand_s_l(L,1:np)-soil_wf_sand_s_l(L,i,j))**2)
            soil_wf_sand_s_l           (L,i,j) = sum(a_wf_sand_s_l          (L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(a_theta_s_l(L,1:np)-soil_theta_s_l(L,i,j))**2/&
                                                 sum(1./(a_theta_s_l(L,1:np)-soil_theta_s_l(L,i,j))**2)
            soil_theta_s_l             (L,i,j) = sum(a_theta_s_l            (L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(a_psi_s_l(L,1:np)-soil_psi_s_l(L,i,j))**2/&
                                                 sum(1./(a_psi_s_l(L,1:np)-soil_psi_s_l(L,i,j))**2)
            soil_psi_s_l               (L,i,j) = sum(a_psi_s_l              (L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(a_lambda_l(L,1:np)-soil_lambda_l(L,i,j))**2/&
                                                 sum(1./(a_lambda_l(L,1:np)-soil_lambda_l(L,i,j))**2)
            soil_lambda_l              (L,i,j) = sum(a_lambda_l             (L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(log10(a_k_s_l(L,1:np))-log10(soil_k_s_l(L,i,j)))**2/&
                                                 sum(1./(log10(a_k_s_l(L,1:np))-log10(soil_k_s_l(L,i,j)))**2)
            soil_k_s_l                 (L,i,j) = 10**(sum(log10(a_k_s_l    (L,1:np))*a_area_fine_gridcell(L,1:np)))
              a_area_fine_gridcell(L,1:np)     =     1./(a_csol_l(L,1:np)-soil_csol_l(L,i,j))**2/&
                                                 sum(1./(a_csol_l(L,1:np)-soil_csol_l(L,i,j))**2)
            soil_csol_l                (L,i,j) = sum(a_csol_l               (L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(a_k_solids_l(L,1:np)-soil_k_solids_l(L,i,j))**2/&
                                                 sum(1./(a_k_solids_l(L,1:np)-soil_k_solids_l(L,i,j))**2)
            soil_k_solids_l            (L,i,j) = sum(a_k_solids_l           (L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(a_tksatu_l(L,1:np)-soil_tksatu_l(L,i,j))**2/&
                                                 sum(1./(a_tksatu_l(L,1:np)-soil_tksatu_l(L,i,j))**2)
            soil_tksatu_l              (L,i,j) = sum(a_tksatu_l             (L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(a_tksatf_l(L,1:np)-soil_tksatf_l(L,i,j))**2/&
                                                 sum(1./(a_tksatf_l(L,1:np)-soil_tksatf_l(L,i,j))**2)
            soil_tksatf_l              (L,i,j) = sum(a_tksatf_l             (L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(a_tkdry_l(L,1:np)-soil_tkdry_l(L,i,j))**2/&
                                                 sum(1./(a_tkdry_l(L,1:np)-soil_tkdry_l(L,i,j))**2)
            soil_tkdry_l               (L,i,j) = sum(a_tkdry_l              (L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(a_VGM_theta_r_l(L,1:np)-soil_VGM_theta_r_l(L,i,j))**2/&
                                                 sum(1./(a_VGM_theta_r_l(L,1:np)-soil_VGM_theta_r_l(L,i,j))**2)
            soil_VGM_theta_r_l         (L,i,j) = sum(a_VGM_theta_r_l        (L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(log10(a_VGM_alpha_l(L,1:np))-log10(soil_VGM_alpha_l(L,i,j)))**2/&
                                                 sum(1./(log10(a_VGM_alpha_l(L,1:np))-log10(soil_VGM_alpha_l(L,i,j)))**2)
            soil_VGM_alpha_l           (L,i,j) = 10**(sum(log10(a_VGM_alpha_l(L,1:np))*a_area_fine_gridcell(L,1:np)))
              a_area_fine_gridcell(L,1:np)     =     1./(log10(a_VGM_n_l(L,1:np))-log10(soil_VGM_n_l(L,i,j)))**2/&
                                                 sum(1./(log10(a_VGM_n_l(L,1:np))-log10(soil_VGM_n_l(L,i,j)))**2)
            soil_VGM_n_l               (L,i,j) = 10**(sum(log10(a_VGM_n_l  (L,1:np))*a_area_fine_gridcell(L,1:np)))
              a_area_fine_gridcell(L,1:np)     =     1./(a_VGM_L_l(L,1:np)-soil_VGM_L_l(L,i,j))**2/&
                                                 sum(1./(a_VGM_L_l(L,1:np)-soil_VGM_L_l(L,i,j))**2)
            soil_VGM_L_l               (L,i,j) = sum(a_VGM_L_l              (L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(a_BA_alpha_l(L,1:np)-soil_BA_alpha_l(L,i,j))**2/&
                                                 sum(1./(a_BA_alpha_l(L,1:np)-soil_BA_alpha_l(L,i,j))**2)
            soil_BA_alpha_l            (L,i,j) = sum(a_BA_alpha_l           (L,1:np)*a_area_fine_gridcell(L,1:np))
              a_area_fine_gridcell(L,1:np)     =     1./(a_BA_beta_l(L,1:np)-soil_BA_beta_l(L,i,j))**2/&
                                                 sum(1./(a_BA_beta_l(L,1:np)-soil_BA_beta_l(L,i,j))**2)
            soil_BA_beta_l             (L,i,j) = sum(a_BA_beta_l            (L,1:np)*a_area_fine_gridcell(L,1:np))
#endif
#elif(defined SOILPAR_UPS_GEOMETRIC)
            a_area_fine_gridcell(L,1:np)       = a_area_fine_gridcell(L,1:np)/sum(a_area_fine_gridcell(L,1:np))
            soil_vf_quartz_mineral_s_l (L,i,j) = product(a_vf_quartz_mineral_s_l(L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_vf_gravels_s_l        (L,i,j) = product(a_vf_gravels_s_l       (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_vf_om_s_l             (L,i,j) = product(a_vf_om_s_l            (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_vf_sand_s_l           (L,i,j) = product(a_vf_sand_s_l          (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_wf_gravels_s_l        (L,i,j) = product(a_wf_gravels_s_l       (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_wf_sand_s_l           (L,i,j) = product(a_wf_sand_s_l          (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_theta_s_l             (L,i,j) = product(a_theta_s_l            (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_psi_s_l               (L,i,j) = product(a_psi_s_l              (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_lambda_l              (L,i,j) = product(a_lambda_l             (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_k_s_l                 (L,i,j) = product(a_k_s_l                (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_csol_l                (L,i,j) = product(a_csol_l               (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_k_solids_l            (L,i,j) = product(a_k_solids_l           (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_tksatu_l              (L,i,j) = product(a_tksatu_l             (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_tksatf_l              (L,i,j) = product(a_tksatf_l             (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_tkdry_l               (L,i,j) = product(a_tkdry_l              (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_VGM_theta_r_l         (L,i,j) = product(a_VGM_theta_r_l        (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_VGM_alpha_l           (L,i,j) = product(a_VGM_alpha_l          (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_VGM_n_l               (L,i,j) = product(a_VGM_n_l              (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_VGM_L_l               (L,i,j) = product(a_VGM_L_l              (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_BA_alpha_l            (L,i,j) = product(a_BA_alpha_l           (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_BA_beta_l             (L,i,j) = product(a_BA_beta_l            (L,1:np)**a_area_fine_gridcell(L,1:np))
#elif(defined SOILPAR_UPS_FIT)
            a_area_fine_gridcell(L,1:np)       = a_area_fine_gridcell(L,1:np)/sum(a_area_fine_gridcell(L,1:np))

! Initial values for the fitting methods            
            soil_vf_quartz_mineral_s_l (L,i,j) = sum(a_vf_quartz_mineral_s_l(L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_vf_gravels_s_l        (L,i,j) = sum(a_vf_gravels_s_l       (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_vf_om_s_l             (L,i,j) = sum(a_vf_om_s_l            (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_vf_sand_s_l           (L,i,j) = sum(a_vf_sand_s_l          (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_wf_gravels_s_l        (L,i,j) = sum(a_wf_gravels_s_l       (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_wf_sand_s_l           (L,i,j) = sum(a_wf_sand_s_l          (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_theta_s_l             (L,i,j) = sum(a_theta_s_l            (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_csol_l                (L,i,j) = sum(a_csol_l               (L,1:np)*a_area_fine_gridcell(L,1:np))
            soil_k_solids_l            (L,i,j) = product(a_k_solids_l      (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_tksatu_l              (L,i,j) = product(a_tksatu_l        (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_tksatf_l              (L,i,j) = product(a_tksatf_l        (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_tkdry_l               (L,i,j) = product(a_tkdry_l         (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_k_s_l                 (L,i,j) = product(a_k_s_l           (L,1:np)**a_area_fine_gridcell(L,1:np))
            soil_psi_s_l               (L,i,j) = median (a_psi_s_l         (L,1:np), np)
            soil_lambda_l              (L,i,j) = median (a_lambda_l        (L,1:np), np)
            soil_VGM_theta_r_l         (L,i,j) = median (a_VGM_theta_r_l   (L,1:np), np)
            soil_VGM_alpha_l           (L,i,j) = median (a_VGM_alpha_l     (L,1:np), np)
            soil_VGM_n_l               (L,i,j) = median (a_VGM_n_l         (L,1:np), np)
            soil_VGM_L_l               (L,i,j) = median (a_VGM_L_l         (L,1:np), np)
            soil_BA_alpha_l            (L,i,j) = median (a_BA_alpha_l      (L,1:np), np)
            soil_BA_beta_l             (L,i,j) = median (a_BA_beta_l       (L,1:np), np)

! SW retentions and Ke-Sr relationship at fine grids (for a specific LCT)
            do LL = 1,np
               ydatc(LL,:) = (-1.0*xdat/a_psi_s_l(L,LL))**(-1.0*a_lambda_l(L,LL)) * a_theta_s_l(L,LL)
               ydatv(LL,:) = a_VGM_theta_r_l(L,LL)+(a_theta_s_l(L,LL) - a_VGM_theta_r_l(L,LL)) &
                             *(1+(a_VGM_alpha_l(L,LL)*xdat)**a_VGM_n_l(L,LL))**(1.0/a_VGM_n_l(L,LL)-1)
               ydatb(LL,:) = xdatsr**(0.5*(1.0+a_vf_om_s_l(L,LL)-a_BA_alpha_l(L,LL)*a_vf_sand_s_l(L,LL) &
                             -a_vf_gravels_s_l(L,LL))) * ((1.0/(1.0+exp(-a_BA_beta_l(L,LL)*xdatsr)))**3 &
                             -((1.0-xdatsr)/2.0)**3)**(1.0-a_vf_om_s_l(L,LL))
            end do

! Fitting the Campbell SW retention parameters             
            ldfjac = np
            xc(1)  = soil_psi_s_l (L,i,j)
            xc(2)  = soil_lambda_l(L,i,j)
            maxfev = 100 * ( nc + 1 )               ! maximum number of iteration
            isiter = 1

            ! the Levenberg–Marquardt fitting method
            call lmder ( SW_CB_dist, np, nc, xc, fvecc(1:np), fjacc(1:np,:), ldfjac, ftol, xtol, gtol, maxfev, &
                   diagc, mode, factor, nprint, info, nfev, njev, ipvtc, qtfc,&
                   xdat, npointw, ydatc(1:np,:), np, soil_theta_s_l(L,i,j), isiter)

            if( xc(1) >= -300. .and. xc(1) < 0.0 .and. xc(2) > 0.0 .and. xc(2) <= 1.0 .and. isiter == 1)then
                soil_psi_s_l (L,i,j) = xc(1)
                soil_lambda_l(L,i,j) = xc(2)
            end if

! Fitting the van Genuchten SW retention parameters              
!            xv(1) = soil_VGM_theta_r_l(L,i,j)
!            xv(2) = soil_VGM_alpha_l  (L,i,j)
!            xv(3) = soil_VGM_n_l      (L,i,j)
!            maxfev = 100 * ( nv + 1 )
!            isiter = 1

!            call lmder ( SW_VG_dist, np, nv, xv, fvecv(1:np), fjacv(1:np,:), ldfjac, ftol, xtol, gtol, maxfev, &
!                   diagv, mode, factor, nprint, info, nfev, njev, ipvtv, qtfv,&
!                   xdat, npointw, ydatv(1:np,:), np, soil_theta_s_l(L,i,j), isiter)

!            if ( xv(1) >= 0.0 .and. xv(1) <= soil_theta_s_l(L,i,j) .and. xv(2) >= 1.0e-5 .and. xv(2) <= 1.0 .and. &
!                 xv(3) >= 1.1 .and. xv(3) <= 10.0 .and. isiter == 1) then
!                 soil_VGM_theta_r_l(L,i,j) = xv(1)
!                 soil_VGM_alpha_l  (L,i,j) = xv(2)
!                 soil_VGM_n_l      (L,i,j) = xv(3)
!            end if

! Fitting the parameters in the Balland and Arp (2005) Ke-Sr relationship            
!            xb(1) = soil_BA_alpha_l(L,i,j)
!            xb(2) = soil_BA_beta_l (L,i,j) 
!            maxfev = 100 * ( nb + 1 )
!            isiter = 1

!            call lmder ( Ke_Sr_dist, np, nb, xb, fvecb(1:np), fjacb(1:np,:), ldfjac, ftol, xtol, gtol, maxfev, &
!                   diagb, mode, factor, nprint, info, nfev, njev, ipvtb, qtfb,&
!                   xdatsr,npointb,ydatb(1:np,:), np, soil_theta_s_l(L,i,j), isiter,&
!                   soil_vf_om_s_l(L,i,j),soil_vf_sand_s_l(L,i,j),soil_vf_gravels_s_l(L,i,j))

!            if ( 1+soil_vf_om_s_l(L,i,j)-xb(1)*soil_vf_sand_s_l(L,i,j)-soil_vf_gravels_s_l(L,i,j) > 0. .and. &
!                 xb(2) > 0.0 .and. isiter == 1  ) then
!                 soil_BA_alpha_l(L,i,j) = xb(1)
!                 soil_BA_beta_l (L,i,j) = xb(2)
!            end if
#endif
                     endif

                  else          ! OCEAN
                     soil_vf_quartz_mineral_s_l (L,i,j) = -1.0e36
                     soil_vf_gravels_s_l        (L,i,j) = -1.0e36
                     soil_vf_om_s_l             (L,i,j) = -1.0e36
                     soil_vf_sand_s_l           (L,i,j) = -1.0e36
                     soil_wf_gravels_s_l        (L,i,j) = -1.0e36
                     soil_wf_sand_s_l           (L,i,j) = -1.0e36
                     soil_theta_s_l             (L,i,j) = -1.0e36
                     soil_psi_s_l               (L,i,j) = -1.0e36
                     soil_lambda_l              (L,i,j) = -1.0e36
                     soil_k_s_l                 (L,i,j) = -1.0e36
                     soil_csol_l                (L,i,j) = -1.0e36
                     soil_k_solids_l            (L,i,j) = -1.0e36
                     soil_tksatu_l              (L,i,j) = -1.0e36
                     soil_tksatf_l              (L,i,j) = -1.0e36
                     soil_tkdry_l               (L,i,j) = -1.0e36
                     soil_VGM_theta_r_l         (L,i,j) = -1.0e36
                     soil_VGM_alpha_l           (L,i,j) = -1.0e36
                     soil_VGM_n_l               (L,i,j) = -1.0e36
                     soil_VGM_L_l               (L,i,j) = -1.0e36
                     soil_BA_alpha_l            (L,i,j) = -1.0e36
                     soil_BA_beta_l             (L,i,j) = -1.0e36
                  endif
               enddo

            enddo
         enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

! (1) Write-out the volumetric fraction of quartz within mineral soil
         lndname = trim(dir_model_landdata)//'model_vf_quartz_mineral_s_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_vf_quartz_mineral_s_l
         close(iunit)

! (2) Write-out the volumetric fraction of gravels
         lndname = trim(dir_model_landdata)//'model_vf_gravels_s_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_vf_gravels_s_l
         close(iunit)

! (3) Write-out the volumetric fraction of organic matter
         lndname = trim(dir_model_landdata)//'model_vf_om_s_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_vf_om_s_l
         close(iunit)

! (4) Write-out volumetric fraction of sand
         lndname = trim(dir_model_landdata)//'model_vf_sand_s_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_vf_sand_s_l
         close(iunit)

! (5) Write-out the gravimetric fraction of gravels
         lndname = trim(dir_model_landdata)//'model_wf_gravels_s_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_wf_gravels_s_l
         close(iunit)

! (6) Write-out gravimetric fraction of sand
         lndname = trim(dir_model_landdata)//'model_wf_sand_s_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_wf_sand_s_l
         close(iunit)

! (7) Write-out the saturated water content [cm3/cm3]
         lndname = trim(dir_model_landdata)//'model_theta_s_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_theta_s_l
         close(iunit)

! (8) Write-out the matric potential at saturation [cm]
         lndname = trim(dir_model_landdata)//'model_psi_s_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_psi_s_l
         close(iunit)

! (9) Write-out the pore size distribution index [dimensionless]
         lndname = trim(dir_model_landdata)//'model_lambda_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_lambda_l
         close(iunit)

! (10) Write-out the saturated hydraulic conductivity [cm/day]
         lndname = trim(dir_model_landdata)//'model_k_s_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_k_s_l
         close(iunit)

! (11) Write-out the heat capacity of soil solids [J/(m3 K)]
         lndname = trim(dir_model_landdata)//'model_csol_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_csol_l
         close(iunit)

! (12) Write-out the thermal conductivity of soil solids [W/m-K]
         lndname = trim(dir_model_landdata)//'model_k_solids_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_k_solids_l
         close(iunit)

! (13) Write-out the thermal conductivity of unfrozen saturated soil [W/m-K]
         lndname = trim(dir_model_landdata)//'model_tksatu_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_tksatu_l
         close(iunit)

! (14) Write-out the thermal conductivity of frozen saturated soil [W/m-K]
         lndname = trim(dir_model_landdata)//'model_tksatf_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_tksatf_l
         close(iunit)

! (15) Write-out the thermal conductivity for dry soil [W/(m-K)]
         lndname = trim(dir_model_landdata)//'model_tkdry_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_tkdry_l
         close(iunit)

! (16) Write-out the VGM's residual moisture content
         lndname = trim(dir_model_landdata)//'model_VGM_theta_r_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_VGM_theta_r_l
         close(iunit)

! (17) Write-out the VGM's parameter corresponding approximately to the inverse of the air-entry value
         lndname = trim(dir_model_landdata)//'model_VGM_alpha_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_VGM_alpha_l
         close(iunit) 

! (18) Write-out the VGM's shape parameter
         lndname = trim(dir_model_landdata)//'model_VGM_n_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_VGM_n_l
         close(iunit)

! (19) Write-out the VGM's pore-connectivity parameter
         lndname = trim(dir_model_landdata)//'model_VGM_L_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_VGM_L_l
         close(iunit)

! (20) Write-out the parameter alpha in the Balland V. and P. A. Arp (2005) model
         lndname = trim(dir_model_landdata)//'model_BA_alpha_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_BA_alpha_l
         close(iunit)

! (21) Write-out the parameter beta in the Balland V. and P. A. Arp (2005) model
         lndname = trim(dir_model_landdata)//'model_BA_beta_l'//trim(c)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_BA_beta_l
         close(iunit)

      ENDDO

! Deallocate the allocatable array
! --------------------------------
      deallocate ( num_patches )

      deallocate ( vf_quartz_mineral_s_l      )
      deallocate ( vf_gravels_s_l             )
      deallocate ( vf_om_s_l                  )
      deallocate ( vf_sand_s_l                )
      deallocate ( wf_gravels_s_l             )
      deallocate ( wf_sand_s_l                )
      deallocate ( theta_s_l                  )
      deallocate ( psi_s_l                    )
      deallocate ( lambda_l                   )
      deallocate ( k_s_l                      )
      deallocate ( csol_l                     )
      deallocate ( k_solids_l                 )
      deallocate ( tksatu_l                   )
      deallocate ( tksatf_l                   )
      deallocate ( tkdry_l                    )
      deallocate ( VGM_theta_r_l              )
      deallocate ( VGM_alpha_l                )
      deallocate ( VGM_n_l                    )
      deallocate ( VGM_L_l                    )

      deallocate ( a_vf_quartz_mineral_s_l    )
      deallocate ( a_vf_gravels_s_l           )
      deallocate ( a_vf_om_s_l                )
      deallocate ( a_vf_sand_s_l              )
      deallocate ( a_wf_gravels_s_l           )
      deallocate ( a_wf_sand_s_l              )
      deallocate ( a_theta_s_l                )
      deallocate ( a_psi_s_l                  )
      deallocate ( a_lambda_l                 )
      deallocate ( a_k_s_l                    )
      deallocate ( a_csol_l                   )
      deallocate ( a_k_solids_l               )
      deallocate ( a_tksatu_l                 )
      deallocate ( a_tksatf_l                 )
      deallocate ( a_tkdry_l                  )
      deallocate ( a_VGM_theta_r_l            )
      deallocate ( a_VGM_alpha_l              )
      deallocate ( a_VGM_n_l                  )
      deallocate ( a_VGM_L_l                  )
      deallocate ( a_BA_alpha_l               )
      deallocate ( a_BA_beta_l                )
      deallocate ( a_area_fine_gridcell       )

#if(defined Naqu || defined Gravels0)
      deallocate ( int_soil_sand_l    ,&
                   int_soil_clay_l    ,&
                   int_soil_oc_l      ,&
                   int_soil_bd_l        )
      deallocate ( VGM_theta_r_Rose   ,&
                   VGM_alpha_Rose     ,&
                   VGM_n_Rose         ,&
                   k_s_Rose             )
      deallocate ( zsoi, dzsoi, zsoih )
#endif 

#if(defined SOILPAR_UPS_FIT)
      deallocate ( ydatc                      )
      deallocate ( ydatv                      )
      deallocate ( ydatb                      )
      deallocate ( fjacc                      )
      deallocate ( fjacv                      )
      deallocate ( fjacb                      )
      deallocate ( fvecc                      )
      deallocate ( fvecv                      )
      deallocate ( fvecb                      )
#endif

      deallocate ( soil_vf_quartz_mineral_s_l )
      deallocate ( soil_vf_gravels_s_l        )
      deallocate ( soil_vf_om_s_l             )
      deallocate ( soil_vf_sand_s_l           )
      deallocate ( soil_wf_gravels_s_l        )
      deallocate ( soil_wf_sand_s_l           )
      deallocate ( soil_theta_s_l             )
      deallocate ( soil_psi_s_l               )
      deallocate ( soil_lambda_l              )
      deallocate ( soil_k_s_l                 )
      deallocate ( soil_csol_l                )
      deallocate ( soil_k_solids_l            )
      deallocate ( soil_tksatu_l              )
      deallocate ( soil_tksatf_l              )
      deallocate ( soil_tkdry_l               )
      deallocate ( soil_VGM_theta_r_l         )
      deallocate ( soil_VGM_alpha_l           )
      deallocate ( soil_VGM_n_l               )
      deallocate ( soil_VGM_L_l               )
      deallocate ( soil_BA_alpha_l            )
      deallocate ( soil_BA_beta_l             )

      deallocate ( landtypes )

      go to 1000
100   print 101,nrow,lndname
101   format(' record =',i8,',  error occured on file: ',a50)
1000  continue

END SUBROUTINE aggregation_soil_parameters
!-----------------------------------------------------------------------
!EOP
