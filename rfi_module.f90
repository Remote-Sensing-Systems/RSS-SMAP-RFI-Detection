module rfi_module

use l2_module_smap

implicit none
save

byte, allocatable, dimension(:,:,:)         :: rfi_flag_chi2_cluster, rfi_flag_FA_cluster, rfi_flag_S4_cluster
byte, allocatable, dimension(:,:,:,:)       :: l_flag

byte, allocatable, dimension(:,:,:)         :: rfi_flag


! RFI parameters. The standard values are indicated as comment.
! The values need to be set in the calling routine.

! 28 JAN 2025
integer(4), parameter                   :: n_chi2_cases  = 3   ! cases considedered for chi2 clustering
real(4), dimension(n_chi2_cases)        :: thr_chi2            ! = (/1.0, 2.0, 3.0/)
integer(4), dimension(n_chi2_cases)     :: ncluster_chi2_1     ! = (/4,   2,   1/)
integer(4), dimension(n_chi2_cases)     :: ncluster_chi2_2     ! = (/6,   3,   2/)

integer(4), parameter                   :: n_FA_cases  = 5     ! cases considedered for FA clustering
real(4), dimension(n_FA_cases)          :: thr_FA              ! = (/0.8, 1.0, 1.5, 2.0, 2.5/)
integer(4), dimension(n_FA_cases)       :: ncluster_FA_1       ! = (/8,   7,   5,   3,   1/)
integer(4), dimension(n_FA_cases)       :: ncluster_FA_2       ! = (/9,   8,   6,   4,   2/)

integer(4), parameter                   :: n_S4_cases  = 5     ! cases considedered for S4 clustering
real(4), dimension(n_S4_cases)          :: thr_S4              ! = (/1.0, 1.3,  1.5, 2.0, 2.3/)
integer(4), dimension(n_S4_cases)       :: ncluster_S4_1       ! = (/7,   5,    4,   2,   1 /)
integer(4), dimension(n_S4_cases)       :: ncluster_S4_2       ! = (/8,   6,    5,   3,   2 /)


integer(4)                              :: iWW                 != 1 window size for counting clusters

integer(4)                              :: iNN_1               != 1 window for NN check low  RFI level
integer(4)                              :: iNN_2               != 3 window for NN check high RFI level

real(4)                                 :: xland_thr           ! = 0.01 land fraction threshold for FA and S4 checks


contains


subroutine allocate_rfi_arrays
use l2_module_smap
implicit none

allocate(rfi_flag(2,nlon,nlat))

allocate(rfi_flag_chi2_cluster(2,nlon,nlat))
allocate(rfi_flag_FA_cluster(2,nlon,nlat))
allocate(rfi_flag_S4_cluster(2,nlon,nlat))
allocate(l_flag(2,nlon,nlat,5))

return
end subroutine allocate_rfi_arrays


subroutine find_rfi_flag
use l2_module_smap
implicit none

byte, dimension(2,nlon,nlat)       :: rfi_flag_0
logical(1), dimension(2,nlon,nlat) :: rfi_flag_ext

integer(4)                         :: idir, ix, iy

call find_rfi_flag_chi2_cluster
call find_rfi_flag_FA_cluster
call find_rfi_flag_S4_cluster

do idir=1,2
do ix=1,nlon
do iy=1,nlat

rfi_flag_0(idir,ix,iy)=0

if (rfi_flag_chi2_cluster(idir,ix,iy)==2 .or. rfi_flag_FA_cluster(idir,ix,iy)==2 .or. rfi_flag_S4_cluster(idir,ix,iy)==2) then
    rfi_flag_0(idir,ix,iy)=2
else if (rfi_flag_chi2_cluster(idir,ix,iy)==1 .or. rfi_flag_FA_cluster(idir,ix,iy)==1 .or. rfi_flag_S4_cluster(idir,ix,iy)==1) then
    rfi_flag_0(idir,ix,iy)=1
endif

enddo
enddo
enddo

call find_rfi_flag_NN (rfi_flag_0, rfi_flag_ext)

! bit flag
do idir=1,2
do ix=1,nlon
do iy=1,nlat

rfi_flag(idir,ix,iy) = 0

if (rfi_flag_ext(idir,ix,iy)==.TRUE.)        rfi_flag(idir,ix,iy) = ibset(rfi_flag(idir,ix,iy),0)     ! Master Flag

if (rfi_flag_chi2_cluster(idir,ix,iy)==1)    rfi_flag(idir,ix,iy) = ibset(rfi_flag(idir,ix,iy),1)     ! CHI2 low RFI 
if (rfi_flag_chi2_cluster(idir,ix,iy)==2)    rfi_flag(idir,ix,iy) = ibset(rfi_flag(idir,ix,iy),2)     ! CHI2 high RFI

if (rfi_flag_FA_cluster(idir,ix,iy)==1)      rfi_flag(idir,ix,iy) = ibset(rfi_flag(idir,ix,iy),3)     ! FA low RFI
if (rfi_flag_FA_cluster(idir,ix,iy)==2)      rfi_flag(idir,ix,iy) = ibset(rfi_flag(idir,ix,iy),4)     ! FA high RFI

if (rfi_flag_S4_cluster(idir,ix,iy)==1)      rfi_flag(idir,ix,iy) = ibset(rfi_flag(idir,ix,iy),5)     ! S4 low RFI
if (rfi_flag_S4_cluster(idir,ix,iy)==2)      rfi_flag(idir,ix,iy) = ibset(rfi_flag(idir,ix,iy),6)     ! S4 high RFI

if (rfi_flag_ext(idir,ix,iy)==.TRUE.  .AND. rfi_flag_0(idir,ix,iy)==0) rfi_flag(idir,ix,iy) = ibset(rfi_flag(idir,ix,iy),7) ! NN only RFI 


enddo
enddo
enddo


return
end subroutine find_rfi_flag



subroutine find_rfi_flag_chi2_cluster
use l2_module_smap
implicit none

integer(4) :: idir, ix, iy, icase, kx, ky, jx, jy, iflag
real(4)    :: chi, chi2

integer(4) :: nn


l_flag = 0
rfi_flag_chi2_cluster = 0

do icase =1,n_chi2_cases

do idir=1,2
do ix=1,nlon
do iy=1,nlat

    if(btest(iqc_flag(idir,ix,iy), 0)) cycle ! no data
    if(btest(iqc_flag(idir,ix,iy), 1)) cycle ! OI wt 
    if(btest(iqc_flag(idir,ix,iy), 2)) cycle ! strong land
    if(btest(iqc_flag(idir,ix,iy), 3)) cycle ! strong sea ice
    if(btest(iqc_flag(idir,ix,iy), 5)) cycle ! sunglint
    if(btest(iqc_flag(idir,ix,iy), 6)) cycle ! moonglint
    
    chi = tb_consistency(idir,ix,iy)
    if ( abs(chi-missing_val4)<0.1) cycle
    chi2 = chi*chi
    
    if (chi2 >= thr_chi2(icase)) then ! count # of cluster
    nn=0
    do jx=ix-iww,ix+iww,1
    do jy=iy-iww,iy+iww,1 
    
        kx=jx
        ky=jy
        if (kx<1)  cycle
        if (kx>nlon) cycle
        if (ky<1)  cycle
        if (ky>nlat) cycle    
        
        chi = tb_consistency(idir,kx,ky)
        if ( abs(chi-missing_val4)<0.1) cycle
        chi2 = chi*chi
        
        if (chi2 >= thr_chi2(icase)) nn = nn + 1
    
    enddo !jy
    enddo !jx
    
    if (nn >= ncluster_chi2_1(icase)) l_flag(idir,ix,iy,icase) = 1 
    if (nn >= ncluster_chi2_2(icase)) l_flag(idir,ix,iy,icase) = 2 
   
    endif ! chi2 exceeds threshold
    
enddo !iy
enddo !ix
enddo !idir

enddo !icase

! aggregate the n_chi2_case

do idir=1,2
do ix=1,nlon
do iy=1,nlat

    iflag=0
    do icase =1,n_chi2_cases
        if (l_flag(idir,ix,iy,icase)==1) then
            iflag=1
            exit
        endif       
    enddo    ! icase loop
    
    if (iflag==1) then
        rfi_flag_chi2_cluster(idir,ix,iy) = 1
    endif

enddo !iy
enddo !ix
enddo !idir

do idir=1,2
do ix=1,nlon
do iy=1,nlat

    iflag=0
    do icase =1,n_chi2_cases
        if (l_flag(idir,ix,iy,icase)==2) then
            iflag=2
            exit
        endif       
    enddo    ! icase loop
    
    if (iflag==2) then
        rfi_flag_chi2_cluster(idir,ix,iy) = 2
    endif

enddo !iy
enddo !ix
enddo !idir

return
end subroutine find_rfi_flag_chi2_cluster



subroutine find_rfi_flag_FA_cluster
use l2_module_smap
implicit none

integer(4) :: idir, ix, iy, icase, kx, ky, jx, jy, iflag
integer(4) :: nn

real(4) :: TBV_fore, TBV_aft, TBH_fore, TBH_aft, DTBV, DTBH, DTBV_1, DTBH_1, DTBV_2, DTBH_2

l_flag = 0
rfi_flag_FA_cluster = 0

! the fore-aft check should only be done over open ocean
! near land or sea-ice fore/aft can differ

do icase =1,n_FA_cases

do ix =1,nlon
do iy =1,nlat

    ! fore
    if(btest(iqc_flag(1,ix,iy), 0)) cycle ! no data
    if(btest(iqc_flag(1,ix,iy), 1)) cycle ! OI wt 
    if(btest(iqc_flag(1,ix,iy), 2)) cycle ! strong land
    if(btest(iqc_flag(1,ix,iy), 3)) cycle ! strong sea ice
    if(btest(iqc_flag(1,ix,iy), 5)) cycle ! sunglint
    if(btest(iqc_flag(1,ix,iy), 6)) cycle ! moonglint    
    if(btest(iqc_flag(1,ix,iy), 8)) cycle ! moderate land
    if(btest(iqc_flag(1,ix,iy), 9)) cycle ! moderate sea-ice
    if(btest(iqc_flag(1,ix,iy),14)) cycle ! light sea-ice    
    if (gland(1,ix,iy) >= xland_thr) cycle

    ! aft
    if(btest(iqc_flag(2,ix,iy), 0)) cycle ! no data
    if(btest(iqc_flag(2,ix,iy), 1)) cycle ! OI wt 
    if(btest(iqc_flag(2,ix,iy), 2)) cycle ! strong land
    if(btest(iqc_flag(2,ix,iy), 3)) cycle ! strong sea ice
    if(btest(iqc_flag(2,ix,iy), 5)) cycle ! sunglint
    if(btest(iqc_flag(2,ix,iy), 6)) cycle ! moonglint
    if(btest(iqc_flag(2,ix,iy), 8)) cycle ! moderate land
    if(btest(iqc_flag(2,ix,iy), 9)) cycle ! moderate sea-ice
    if(btest(iqc_flag(2,ix,iy),14)) cycle ! light sea-ice  
    if (gland(2,ix,iy) >= xland_thr) cycle    
    
    TBV_fore =   tb_sur0(1,1,ix,iy)
    if (abs(TBV_fore - missing_val4) < 0.1) cycle
    TBV_aft  =   tb_sur0(1,2,ix,iy)
    if (abs(TBV_aft - missing_val4)  < 0.1) cycle

    TBH_fore =   tb_sur0(2,1,ix,iy)
    if (abs(TBH_fore - missing_val4) < 0.1) cycle
    TBH_aft  =   tb_sur0(2,2,ix,iy)
    if (abs(TBH_aft - missing_val4)  < 0.1) cycle
 
    DTBV = TBV_fore - TBV_aft
    DTBH = TBH_fore - TBH_aft   
    
    ! RFI in fore look
    if (DTBV >= thr_FA(icase) .or. DTBH >= thr_FA(icase)) then ! Fa diff exceeds threshold
    nn=0
    do jx=ix-iww,ix+iww,1
    do jy=iy-iww,iy+iww,1 
    
        kx=jx
        ky=jy
        if (kx<1)  cycle
        if (kx>nlon) cycle
        if (ky<1)  cycle
        if (ky>nlat) cycle    
 
        ! fore
        if(btest(iqc_flag(1,kx,ky), 0)) cycle ! no data
        if(btest(iqc_flag(1,kx,ky), 1)) cycle ! OI wt 
        if(btest(iqc_flag(1,kx,ky), 2)) cycle ! strong land
        if(btest(iqc_flag(1,kx,ky), 3)) cycle ! strong sea ice
        if(btest(iqc_flag(1,kx,ky), 5)) cycle ! sunglint
        if(btest(iqc_flag(1,kx,ky), 6)) cycle ! moonglint    
        if(btest(iqc_flag(1,kx,ky), 8)) cycle ! moderate land
        if(btest(iqc_flag(1,kx,ky), 9)) cycle ! moderate sea-ice
        if(btest(iqc_flag(1,kx,ky),14)) cycle ! light sea-ice    
        if (gland(1,kx,ky) >= xland_thr) cycle

        ! aft
        if(btest(iqc_flag(2,kx,ky), 0)) cycle ! no data
        if(btest(iqc_flag(2,kx,ky), 1)) cycle ! OI wt 
        if(btest(iqc_flag(2,kx,ky), 2)) cycle ! strong land
        if(btest(iqc_flag(2,kx,ky), 3)) cycle ! strong sea ice
        if(btest(iqc_flag(2,kx,ky), 5)) cycle ! sunglint
        if(btest(iqc_flag(2,kx,ky), 6)) cycle ! moonglint
        if(btest(iqc_flag(2,kx,ky), 8)) cycle ! moderate land
        if(btest(iqc_flag(2,kx,ky), 9)) cycle ! moderate sea-ice
        if(btest(iqc_flag(2,kx,ky),14)) cycle ! light sea-ice  
        if (gland(2,kx,ky) >= xland_thr) cycle
        
        TBV_fore =   tb_sur0(1,1,kx,ky)
        if (abs(TBV_fore - missing_val4) < 0.1) cycle
        TBV_aft  =   tb_sur0(1,2,kx,ky)
        if (abs(TBV_aft - missing_val4)  < 0.1) cycle

        TBH_fore =   tb_sur0(2,1,kx,ky)
        if (abs(TBH_fore - missing_val4) < 0.1) cycle
        TBH_aft  =   tb_sur0(2,2,kx,ky)
        if (abs(TBH_aft - missing_val4)  < 0.1) cycle
 
        DTBV_1 = TBV_fore - TBV_aft
        DTBH_1 = TBH_fore - TBH_aft  
        
        if (DTBV_1 >= thr_FA(icase) .or. DTBH_1 >= thr_FA(icase)) nn=nn+1    
    
    enddo !jy
    enddo !jx
    
    if (nn >= ncluster_FA_1(icase)) l_flag(1,ix,iy,icase) = 1 
    if (nn >= ncluster_FA_2(icase)) l_flag(1,ix,iy,icase) = 2
    
    endif ! fore case 
 
     ! RFI in aft look
    if (DTBV <= -thr_FA(icase) .or. DTBH <= -thr_FA(icase)) then ! Fa diff exceeds threshold
    nn=0
    do jx=ix-iww,ix+iww,1
    do jy=iy-iww,iy+iww,1 
    
        kx=jx
        ky=jy
        if (kx<1)  cycle
        if (kx>nlon) cycle
        if (ky<1)  cycle
        if (ky>nlat) cycle    
 
        ! fore
        if(btest(iqc_flag(1,kx,ky), 0)) cycle ! no data
        if(btest(iqc_flag(1,kx,ky), 1)) cycle ! OI wt 
        if(btest(iqc_flag(1,kx,ky), 2)) cycle ! strong land
        if(btest(iqc_flag(1,kx,ky), 3)) cycle ! strong sea ice
        if(btest(iqc_flag(1,kx,ky), 5)) cycle ! sunglint
        if(btest(iqc_flag(1,kx,ky), 6)) cycle ! moonglint    
        if(btest(iqc_flag(1,kx,ky), 8)) cycle ! moderate land
        if(btest(iqc_flag(1,kx,ky), 9)) cycle ! moderate sea-ice
        if(btest(iqc_flag(1,kx,ky),14)) cycle ! light sea-ice    
        if (gland(1,kx,ky) >= xland_thr) cycle

        ! aft
        if(btest(iqc_flag(2,kx,ky), 0)) cycle ! no data
        if(btest(iqc_flag(2,kx,ky), 1)) cycle ! OI wt 
        if(btest(iqc_flag(2,kx,ky), 2)) cycle ! strong land
        if(btest(iqc_flag(2,kx,ky), 3)) cycle ! strong sea ice
        if(btest(iqc_flag(2,kx,ky), 5)) cycle ! sunglint
        if(btest(iqc_flag(2,kx,ky), 6)) cycle ! moonglint
        if(btest(iqc_flag(2,kx,ky), 8)) cycle ! moderate land
        if(btest(iqc_flag(2,kx,ky), 9)) cycle ! moderate sea-ice
        if(btest(iqc_flag(2,kx,ky),14)) cycle ! light sea-ice  
        if (gland(2,kx,ky) >= xland_thr) cycle
        
        TBV_fore =   tb_sur0(1,1,kx,ky)
        if (abs(TBV_fore - missing_val4) < 0.1) cycle
        TBV_aft  =   tb_sur0(1,2,kx,ky)
        if (abs(TBV_aft - missing_val4)  < 0.1) cycle

        TBH_fore =   tb_sur0(2,1,kx,ky)
        if (abs(TBH_fore - missing_val4) < 0.1) cycle
        TBH_aft  =   tb_sur0(2,2,kx,ky)
        if (abs(TBH_aft - missing_val4)  < 0.1) cycle
 
        DTBV_2 = TBV_fore - TBV_aft
        DTBH_2 = TBH_fore - TBH_aft  
        
        if (DTBV_2 <= -thr_FA(icase) .or. DTBH_2 <= -thr_FA(icase)) nn=nn+1    
    
    enddo !jy
    enddo !jx
    
    if (nn >= ncluster_FA_1(icase)) l_flag(2,ix,iy,icase) = 1 
    if (nn >= ncluster_FA_2(icase)) l_flag(2,ix,iy,icase) = 2
    
    endif ! aft case 
 
enddo !iy
enddo !ix

enddo ! icase

! aggregate the n_FA_case

do idir=1,2
do ix=1,nlon
do iy=1,nlat

    iflag=0
    do icase =1,n_FA_cases
        if (l_flag(idir,ix,iy,icase)==1) then
            iflag=1
            exit
        endif       
    enddo    ! icase loop
    
    if (iflag==1) then
        rfi_flag_FA_cluster(idir,ix,iy) = 1
    endif

enddo !iy
enddo !ix
enddo !idir

do idir=1,2
do ix=1,nlon
do iy=1,nlat

    iflag=0
    do icase =1,n_FA_cases
        if (l_flag(idir,ix,iy,icase)==2) then
            iflag=2
            exit
        endif       
    enddo    ! icase loop
    
    if (iflag==2) then
        rfi_flag_FA_cluster(idir,ix,iy) = 2
    endif

enddo !iy
enddo !ix
enddo !idir

return
end subroutine find_rfi_flag_FA_cluster



subroutine find_rfi_flag_S4_cluster
use l2_module_smap

implicit none


integer(4) :: idir, ix, iy, icase, kx, ky, jx, jy, iflag
integer(4) :: nn
real(4) :: S4_meas,  S4

l_flag = 0
rfi_flag_S4_cluster = 0

do icase =1,n_S4_cases

do idir=1,2
do ix  =1,nlon
do iy  =1,nlat

    if(btest(iqc_flag(idir,ix,iy), 0)) cycle ! no data
    if(btest(iqc_flag(idir,ix,iy), 1)) cycle ! OI wt 
    if(btest(iqc_flag(idir,ix,iy), 2)) cycle ! strong land
    if(btest(iqc_flag(idir,ix,iy), 3)) cycle ! strong sea ice
    if(btest(iqc_flag(idir,ix,iy), 5)) cycle ! sunglint
    if(btest(iqc_flag(idir,ix,iy), 6)) cycle ! moonglint

    ! The S4 check should be done only over open ocean.
    ! Near land or sea-ice there is the squint effect, which has not been corrected.
    if(btest(iqc_flag(idir,ix,iy), 8)) cycle ! moderate land
    if(btest(iqc_flag(idir,ix,iy), 9)) cycle ! moderate sea-ice
    if(btest(iqc_flag(idir,ix,iy),14)) cycle ! light sea-ice 
    
    if (gland(idir,ix,iy) >= xland_thr) cycle  ! use same threshold as in FORE - AFT flagging       
    
    S4_meas = ta_ant_calibrated(4,idir,ix,iy) 
    if ( abs(S4_meas-missing_val4)<0.1) cycle

    S4 = S4_meas ! S4_exp = 0
        
    if (abs(S4) >= thr_S4(icase)) then ! S4 exceeds threshold
    nn=0
    do jx=ix-iww,ix+iww,1
    do jy=iy-iww,iy+iww,1 
    
        kx=jx
        ky=jy
        if (kx<1)  cycle
        if (kx>nlon) cycle
        if (ky<1)  cycle
        if (ky>nlat) cycle    
        
        if(btest(iqc_flag(idir,kx,ky), 0)) cycle ! no data
        if(btest(iqc_flag(idir,kx,ky), 1)) cycle ! OI wt 
        if(btest(iqc_flag(idir,kx,ky), 2)) cycle ! strong land
        if(btest(iqc_flag(idir,kx,ky), 3)) cycle ! strong sea ice
        if(btest(iqc_flag(idir,kx,ky), 5)) cycle ! sunglint
        if(btest(iqc_flag(idir,kx,ky), 6)) cycle ! moonglint

        ! The S4 check should be doen only over open ocean.
        ! Near land or sea-ice there is the squint effect, which has not been corrected.
        if(btest(iqc_flag(idir,kx,ky), 8)) cycle ! moderate land
        if(btest(iqc_flag(idir,kx,ky), 9)) cycle ! moderate sea-ice
        if(btest(iqc_flag(idir,kx,ky),14)) cycle ! light sea-ice     
        if (gland(idir,kx,ky) >= xland_thr) cycle  ! use same threshold as in FORE - AFT flagging
        
        S4_meas = ta_ant_calibrated(4,idir,kx,ky) 
        if ( abs(S4_meas-missing_val4)<0.1) cycle

        S4 = S4_meas ! S4_exp = 0        
        
        if (abs(S4) >= thr_S4(icase)) nn = nn + 1
    
    enddo !jy
    enddo !jx
    
    if (nn >= ncluster_S4_1(icase)) l_flag(idir,ix,iy,icase) = 1 
    if (nn >= ncluster_S4_2(icase)) l_flag(idir,ix,iy,icase) = 2    
    
    endif ! S4 exceeds threshold 

enddo ! iy
enddo ! ix
enddo ! idir

enddo ! icase

! aggregate the n_S4_case

do idir=1,2
do ix=1,nlon
do iy=1,nlat

    iflag=0
    do icase =1,n_S4_cases
        if (l_flag(idir,ix,iy,icase)==1) then
            iflag=1
            exit
        endif       
    enddo    ! icase loop
    
    if (iflag==1) then
        rfi_flag_S4_cluster(idir,ix,iy) = 1
    endif

enddo !iy
enddo !ix
enddo !idir

do idir=1,2
do ix=1,nlon
do iy=1,nlat

    iflag=0
    do icase =1,n_S4_cases
        if (l_flag(idir,ix,iy,icase)==2) then
            iflag=2
            exit
        endif       
    enddo    ! icase loop
    
    if (iflag==2) then
        rfi_flag_S4_cluster(idir,ix,iy) = 2
    endif

enddo !iy
enddo !ix
enddo !idir

return
end subroutine find_rfi_flag_S4_cluster


! NN RFI flag
! it flags the NN of a cell that has been already RFI flagged with the level rfi_flag_0 = 0,1,2
! level 1: NN_1
! level 2: NN_2

subroutine find_rfi_flag_NN (rfi_flag_0, rfi_flag_ext)
use l2_module_smap
implicit none

byte, dimension(2,nlon,nlat), intent(in)       :: rfi_flag_0
logical(1), dimension(2,nlon,nlat), intent(out):: rfi_flag_ext

integer(4) :: idir, ix, iy, jx, jy, kx, ky, iNN

rfi_flag_ext = .FALSE.  ! default

do idir=1,2
do ix  =1,nlon
do iy  =1,nlat

    if (rfi_flag_0(idir,ix,iy)/=0) then
    
    if (rfi_flag_0(idir,ix,iy)==1) iNN = iNN_1 
    if (rfi_flag_0(idir,ix,iy)==2) iNN = iNN_2 
   
    do jx=ix-iNN,ix+iNN,1
    do jy=iy-iNN,iy+iNN,1 
    
        kx=jx
        ky=jy
        if (kx<1)    kx=1
        if (kx>nlon) kx=nlon
        if (ky<1)    ky=1
        if (ky>nlat) ky=nlat  

        rfi_flag_ext(idir,kx,ky) = .TRUE.
        
    enddo !ky
    enddo !kx        
    
    endif ! rfi flag 0 is set

enddo ! iy
enddo ! ix
enddo ! idir

return
end subroutine find_rfi_flag_NN


end module rfi_module