! adapted to V6.2 L2C

module  l2_module_smap
implicit none
save

real(8), parameter      :: pi=3.14159265358979323846d0
real(8), parameter      :: rad=pi/180.d0


! geolocation constants 
real(8), parameter      :: re=6378.137d3  !earth equatorial radius (meters)
real(8), parameter      :: rp=6356.752d3  !earth      polar radius (meters)
real(8), parameter      :: ffac=(rp/re)*(rp/re)
      
integer(4), parameter   :: nlat= 720  !for resampling arrays
integer(4), parameter   :: nlon=1560  !360 deg plus 30 deg
integer(4), parameter   :: ndir=   2  !1=for 2=aft

	
integer(4)              :: klon_end
real(4)                 :: xlon_node1,xlon_node2
	

real(4), parameter      :: missing_val4 = -9999.0
real(8), parameter      :: missing_val8 = -9999.d0

real(4), parameter      :: tb_cos=3.0
! Frank has subtracted the isotropic tb_cos = 3K when deriving the reflected galactic tables
! The tb_cos is included in the calculation of tb_gal_dir

real(4), dimension(2), parameter    :: emiss_refl_adj=(/0.01012, 0.01012/)


real(4), parameter :: sss_val_min=1.5, sss_val_max=43.5

real(4), parameter :: sea_ice_adj_1 = 0.72 ! adjust sea-ice correction and gice in zones 1,2,3
real(4), parameter :: sea_ice_adj_2 = 0.7  ! reduce sea-ice error in zone 3
real(4), parameter :: sea_ice_adj_3 = 0.4  ! reduce sea-ice error in zone 4


! opt parameters to switch on(=1) or off(=0) correction
integer(4),dimension(10)  :: iopt
! 1: land intrusion
! 2: gal direct
! 3: gal refl
! 4: sun dir
! 5: sun ref
! 6: swh in roughness correction
! 7: WindSat wind speed in roughness correction


integer(4)                                  :: vernum ! processing version     

integer(4)                                  :: igal_wspd
! =1 use provided wind speed for galaxy roughness
! =2 add 3 m/s to provided wind speed for galaxy roughness according to Frank 

! L2A variables      
integer(1), dimension(  2,nlon,nlat  )      :: ifill
integer(4), dimension(  2,nlon,nlat  )      :: iscan
real(8),    dimension(  2,nlon,nlat  )      :: time
real(4),    dimension(  2,nlon,nlat  )      :: zang,cellat,cellon,alpha,eia,eaa,pra,sunglt,monglt,gallat,gallon,frdrot
real(4),    dimension(  2,nlon,nlat  )      :: teclat,teclon,fland,gland,dist,wt_sum

! added 09/24/2015: sun angles: beta = sun zenith angle   alpha = sun azimuth angle
real(4),    dimension(  2,nlon,nlat  )      :: sun_alpha, sun_beta

real(4),    dimension(4,2,nlon,nlat  )      :: ta_ant, ta_ant_filtered, ta_ant_calibrated, ta_ant_exp  
! this was called ta_erth in Frank's module. the filtered ta were added on 05/10/2015. 
! in VH basis
! ta_ant TA from L1A files. no RFI filter
! ta_ant_filtered TA from L1A files, RFI filtered.
! ta_ant_calibrated TA filtered after correcting emissive reflector and Ta bias correction

real(4),    dimension(3,2,nlon,nlat  )      :: ta_sun_dir,ta_sun_ref,ta_gal_dir,ta_gal_ref ! IQ basis
real(4),    dimension(3,2,nlon,nlat,5)      :: ta_gal_ref_tab 
! IQ basis. This is Frank's pre-computed table for 5 wind speeds

real(4),    dimension(2,2,nlon,nlat  )      :: ta_lnd


real(4),    dimension(2,2,nlon,nlat  )      :: loss_ant
real(4),    dimension(2,2,nlon,nlat  )      :: temp_ant, dtemp_ant
! antenna loss and physical temperature
! 1=V 2=H

! ancillary data (L2B variables)
real(4), dimension(2,nlon,nlat)             :: tec
real(4), dimension  (nlon,nlat)             :: sss_ref,surtep,winspd_anc,windir,tran,tbup,tbdw,lst,sm,solar_flux
real(4), dimension(nlon,nlat)               :: absp_oxy, absp_vap, absp_liq
real(4), dimension(nlon,nlat)               :: winspd, rain


! A MATRIX (APC) and inverse in IQ basis
real(8), dimension(4,4)                     :: amat_IQ, amat_inv_IQ


! new variables in L2C processing
real(4), dimension(4,2,nlon,nlat  )         :: ta_earth,   ta_earth_exp     ! Ta coming from Earth only 
                                                                            ! after removing galaxy, sun, moon, CS  (VH basis)

real(4), dimension(4,2,nlon,nlat  )         :: tb_toi,     tb_toi_exp       ! top of the ionosphere 
                                                                            ! after applying APC (VH basis)
real(4), dimension(4,2,nlon,nlat  )         ::             tb_toi_exp_2     ! top of the ionosphere exp using the TEC FRA + geometric PRA 
 

real(4), dimension(4,2,nlon,nlat  )         :: tb_toa,     tb_toa_exp       ! top of the atmosphere  
                                                                            ! after applying Farday rotation correction (VH basis)
                                                                            
real(4), dimension(4,2,nlon,nlat  )         :: tb_toa_lc                    ! top of the atmpshere land corrected  
                                                                            ! after applying land correction (VH basis)
                                                                            
real(4), dimension(4,2,nlon,nlat  )         :: tb_sur,     tb_sur_exp       ! at rough surface  
                                                                            ! after applying atmospheric correction (VH basis)
                                                                            
real(4), dimension(4,2,nlon,nlat  )         :: tb_sur0,    tb_sur0_exp      ! at flat surface  
                                                                            ! before applying surface roughness correction (VH basis)  

! added in V5
real(4), dimension(4,2,nlon,nlat  )         :: tb_sur0_sic                  ! at flat surface  
                                                                            ! after applying sea-ice correction (VH basis)  
                                                                            
real(4), dimension(2,nlon,nlat)             :: sss_smap, sss_smap_40km      ! retrieved SMAP salinity  
integer(4), dimension(2,nlon,nlat)          :: nsmooth                      ! number of observation in 70-km smoothing window. added in V5.
                                                                         
integer(4), dimension(2,nlon,nlat)          :: iflag_sss_conv               ! SSS convergence flag
real(4), dimension(2,nlon,nlat)             :: tb_consistency               ! TB consistency = sqrt(chisq) from SSS MLE

real(4), dimension(2,nlon,nlat)             :: pra_smap                     ! measured polarization rotation angle [deg]
real(4), dimension(2,nlon,nlat)             :: fra_exp, pratot_exp          ! expected fra and total pra (=fra + pra geo) [deg]

real(4), dimension(4,2,nlon,nlat)           :: dtb_rough                    ! surface roughness correction [K]

integer(4), dimension(2,nlon,nlat)          :: iqc_flag                     ! Q/C flag


real(4), dimension(4)                       :: dtb_bias_orbit
real(4), dimension(4)                       :: dtb_stitch_orbit
 
real(4), dimension(2)                       :: tf_ave_orbit

! changed for V3 processing
real(4), parameter                          :: offset_S3=+0.22,       offset_S4=-0.43 
real(4), parameter                          :: offset_S3_prime=+0.43, offset_S4_prime=-0.17 


!     added in V5	  
real(4), allocatable, dimension(:,:)         :: xobs_sat_ccmp, sst_err, gice_est
real(4), allocatable, dimension(:,:,:)       :: dtb_sea_ice_corr, dtb_sea_ice_error_est  ! V,H
integer(4), allocatable, dimension(:,:)      :: icezone
integer(4), allocatable, dimension(:,:,:)    :: iceflag_amsr2
!1= mask4 (climatology) 2=AMSR2 V8.2 flag   3=AMSR2 discriminant flag

! V5 perturbed fields
! perturbed retrievals

integer(4), parameter   ::  nuncertainties=9
! 1 = ancillary wind speed in roughness correction random
! 2 = NEDT v-pol
! 3 = NEDT h-pol
! 4 = ancillary SST
! 5 = wind direction 
! 6=  reflected galaxy. use perturbations 2 + 3
! 7=  land.             use perturbations 2 + 3
! 8=  sea-ice. .        use perturbations 2 + 3
! 9=  ancillary wind speed in roughness correction systematic

real(4), parameter, dimension(5)                :: dpert = (/0.1, 0.1, 0.1, 0.1, 5.0/)  

real(4), allocatable, dimension(:,:,:)          :: sss_P, sss_M, sss_0, tb_consistency_save
real(4), allocatable, dimension(:,:,:,:)        :: dsss_40km, dsss_70km       
real(4), allocatable, dimension(:,:,:)          :: dsss_tot_40km, dsss_tot_70km                                          

real(4), allocatable, dimension(:,:,:,:)        :: tb_sur0_save
real(4), allocatable, dimension(:,:)            :: surtep_save, windir_save, winspd_save

real(4), allocatable, dimension(:,:,:)          :: dsss_v, dsss_h

                                                                          
end module  l2_module_smap

