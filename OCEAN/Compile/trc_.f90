










!======================================================================
! MARCO is a branch of ROMS/CROCO developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III, IIT, IISC, Indian Univ.
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! MARCO specific routines (AGRIF) are under CeCILL-C license.
!
! MARCO website : http://xxxx
!======================================================================
!







                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      

                      
                      
                      
                      
                      
                      

                      


                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                     
                      
!    Parallel reproducibility test
!    RVTK test (Restartability or Parallel reproducibility)








                      
                      
                      
                      
                      
                      
                      


!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!




























































































































!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!





















!









!-# define float dfloat
!-# define FLoaT dfloat
!-# define FLOAT dfloat
!-# define sqrt dsqrt
!-# define SQRT dsqrt
!-# define exp dexp
!-# define EXP dexp
!-# define dtanh dtanh
!-# define TANH dtanh









MODULE trc
   !!======================================================================
   !!                      ***  MODULE  trc  ***
   !! Passive tracers   :  module for tracers defined
   !!======================================================================
   !! History :   OPA  !  1996-01  (M. Levy)  Original code
   !!              -   !  2000-04  (O. Aumont, M.A. Foujols)  HAMOCC3 and P3ZD
   !!   NEMO      1.0  !  2004-03  (C. Ethe)  Free form and module
   !!----------------------------------------------------------------------

!   USE par_oce
!   USE par_trc
!   USE bdy_oce, only: jp_bdy, ln_bdy, nb_bdy, OBC_DATA
   USE par_pisces
   USE oce_trc

   IMPLICIT NONE
   PUBLIC

   PUBLIC   trc_alloc   ! called by nemogcm.F90


   !                                     !!- logical units of passive tracers
   INTEGER, PUBLIC ::   numont     = -1   !: reference passive tracer namelist output output.namelist.top
   INTEGER, PUBLIC ::   numonr     = -1   !: reference passive tracer namelist output output.namelist.top
   INTEGER, PUBLIC ::   numstr            !: tracer statistics
   INTEGER, PUBLIC ::   numrtr            !: tracer statistics
   INTEGER, PUBLIC ::   numrtw            !: tracer statistics
   CHARACTER(:), ALLOCATABLE, PUBLIC ::   numnat_ref   !: character buffer for reference passive tracer namelist_top_ref
   CHARACTER(:), ALLOCATABLE, PUBLIC ::   numnat_cfg   !: character buffer for configuration specific passive tracer namelist_top_cfg
   CHARACTER(:), ALLOCATABLE, PUBLIC ::   numtrc_ref   !: character buffer for reference passive tracer namelist_trc_ref
   CHARACTER(:), ALLOCATABLE, PUBLIC ::   numtrc_cfg   !: character buffer for configuration specific passive tracer namelist_trc_cfg

   !! passive tracers fields (before,now,after)
   !! --------------------------------------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)       ::  trai           !: initial total tracer
   REAL(wp), PUBLIC                                        ::  areatot        !: total volume 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:  ) ::  cvol           !: volume correction -degrad option- 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:,:) ::  tr           !: tracer concentration 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:  ) ::  sbc_trc_b      !: Before sbc fluxes for tracers
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:  ) ::  sbc_trc        !: Now sbc fluxes for tracers

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:  ) ::  trc_i          !: prescribed tracer concentration in sea ice for SBC
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:  ) ::  trc_o          !: prescribed tracer concentration in ocean for SBC
   INTEGER             , PUBLIC                            ::  nn_ice_tr      !: handling of sea ice tracers
   INTEGER             , PUBLIC                            ::  nn_ais_tr      !: handling of Antarctic Ice Sheet tracers

   !! interpolated gradient
   !!--------------------------------------------------  
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  gtru           !: hor. gradient at u-points at bottom ocean level
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  gtrv           !: hor. gradient at v-points at bottom ocean level
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  gtrui          !: hor. gradient at u-points at top    ocean level
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  gtrvi          !: hor. gradient at v-points at top    ocean level
   
   !! passive tracers  (input and output)
   !! ------------------------------------------  
   LOGICAL             , PUBLIC ::   ln_rsttr           !: boolean term for restart i/o for passive tracers (namelist)
   LOGICAL             , PUBLIC ::   lrst_trc           !: logical to control the trc restart write
   INTEGER             , PUBLIC ::   nn_writetrc        !: time step frequency for concentration outputs (namelist)
   INTEGER             , PUBLIC ::   nutwrs             !: output FILE for passive tracers restart
   INTEGER             , PUBLIC ::   nutrst             !: logical unit for restart FILE for passive tracers
   INTEGER             , PUBLIC ::   nn_rsttr           !: control of the time step ( 0 or 1 ) for pass. tr.
   CHARACTER(len = 80) , PUBLIC ::   cn_trcrst_in       !: suffix of pass. tracer restart name (input)
   CHARACTER(len = 256), PUBLIC ::   cn_trcrst_indir    !: restart input directory
   CHARACTER(len = 80) , PUBLIC ::   cn_trcrst_out      !: suffix of pass. tracer restart name (output)
   CHARACTER(len = 256), PUBLIC ::   cn_trcrst_outdir   !: restart output directory
   REAL(wp)            , PUBLIC ::   rDt_trc            !: = 2*rn_Dt except at nit000 (=rn_Dt) if l_1st_euler=.true.
   LOGICAL             , PUBLIC ::   ln_top_euler       !: boolean term for euler integration 
   LOGICAL             , PUBLIC ::   ln_trcdta          !: Read inputs data from files
   LOGICAL             , PUBLIC ::   ln_trcbc           !: Enable surface, lateral or open boundaries conditions
   LOGICAL             , PUBLIC ::   ln_trcais          !: Enable Antarctic Ice Sheet nutrient supply
   LOGICAL             , PUBLIC ::   ln_trcdmp          !: internal damping flag
   LOGICAL             , PUBLIC ::   ln_trcdmp_clo      !: internal damping flag on closed seas
   INTEGER             , PUBLIC ::   nittrc000          !: first time step of passive tracers model
!   LOGICAL             , PUBLIC ::   l_trcdm2dc         !: Diurnal cycle for TOP

   !! Information for the ice module for tracers
   !! ------------------------------------------
!   TYPE, PUBLIC ::   TRC_I_NML         !: Ice tracer namelist structure
!         REAL(wp)         :: trc_ratio    ! ice-ocean trc ratio
!         REAL(wp)         :: trc_prescr   ! prescribed ice trc cc
!         CHARACTER(len=2) :: ctrc_o       ! choice of ocean trc cc
!   END TYPE
!   !
!   REAL(wp)        , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   trc_ice_ratio    !: ice-ocean tracer ratio
!   REAL(wp)        , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   trc_ice_prescr   !: prescribed ice trc cc

!   CHARACTER(len=lca), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   cn_trc_o         !: choice of ocean tracer cc

   !! Information for the optics module
   !! ---------------------------------
   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  neln       !: number of T-levels + 1 in the euphotic layer
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  heup       !: euphotic layer depth
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  heup_01    !: Absolute euphotic layer depth
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  etot       !: par (photosynthetic available radiation)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  etot_ndcy  !: PAR over 24h in case of diurnal cycle


   !! information for outputs
   !! --------------------------------------------------
   TYPE, PUBLIC ::   PTRACER        !: Passive tracer type
      CHARACTER(len=20) ::   clsname   ! short name
      CHARACTER(len=80) ::   cllname   ! long name
      CHARACTER(len=20) ::   clunit    ! unit
!      LOGICAL           ::   llinit    ! read in a file or not
!      LOGICAL           ::   llsbc     ! read in a file or not
!      LOGICAL           ::   llcbc     ! read in a file or not
!      LOGICAL           ::   llobc     ! read in a file or not
!      LOGICAL           ::   llais     ! read in a file or not
   END TYPE PTRACER
   !
   CHARACTER(len=20), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ctrcnm   !: tracer name 
   CHARACTER(len=80), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ctrcnl   !: trccer field long name
   CHARACTER(len=20), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ctrcnu   !: tracer unit
   !
   TYPE, PUBLIC ::   DIAG         !: Passive trcacer ddditional diagnostic type
      CHARACTER(len=20) ::   sname   ! short name
      CHARACTER(len=80) ::   lname   ! long name
      CHARACTER(len=20) ::   units   ! unit
   END TYPE DIAG
   !
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   trc3d   !: 3D diagnostics for tracers
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   trc2d   !: 2D diagnostics for tracers

   !! information for inputs
   !! --------------------------------------------------
   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ln_trc_ini    !: Initialisation from data input file
   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ln_trc_obc    !: Use open boundary condition data
   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ln_trc_sbc    !: Use surface boundary condition data
   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ln_trc_cbc    !: Use coastal boundary condition data
   LOGICAL , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   ln_trc_ais    !: Use Antarctic Ice Sheet boundary condition data
   LOGICAL , PUBLIC                                  ::   ln_rnf_ctl    !: remove runoff dilution on tracers
   REAL(wp), PUBLIC                                  ::   rn_sbc_time   !: Time scaling factor for SBC data (seconds in a day)
   REAL(wp), PUBLIC                                  ::   rn_cbc_time   !: Time scaling factor for CBC data (seconds in a day)
   LOGICAL , PUBLIC                                  ::   lltrcbc       !: Applying one of the boundary conditions 
   !
!   CHARACTER(len=20), PUBLIC, DIMENSION(jp_bdy) :: cn_trc_dflt   ! Default OBC condition for all tracers
!   CHARACTER(len=20), PUBLIC, DIMENSION(jp_bdy) :: cn_trc        ! Choice of boundary condition for tracers
!   INTEGER,           PUBLIC, DIMENSION(jp_bdy) :: nn_trcdmp_bdy !: =T Tracer damping
!   LOGICAL,           PUBLIC, DIMENSION(jp_bdy) :: ln_zintobc    !: =T obc data requires a vertical interpolation
   !
   ! Vertical axis used in the sediment module
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   profsed

   CHARACTER(len = 80) , PUBLIC   ::  cn_pisrst_in   !: suffix of pass. tracer restart name (input)
   CHARACTER(len = 256), PUBLIC   ::  cn_pisrst_indir  !: restart input directory
   CHARACTER(len = 80) , PUBLIC   ::  cn_pisrst_out  !: suffix of pass. tracer restart name (output)
   CHARACTER(len = 256), PUBLIC   ::  cn_pisrst_outdir  !: restart output directory

!$AGRIF_DO_NOT_TREAT
   ! External data structure of BDY for TOP. Available elements: cn_obc, ll_trc, trcnow, dmp
!   TYPE(OBC_DATA), PUBLIC, ALLOCATABLE, DIMENSION(:,:), TARGET ::   trcdta_bdy   !: bdy external data (local process)
!$AGRIF_END_DO_NOT_TREAT
   !
   !! Substitutions





























   
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trc.F90 15420 2021-10-20 15:26:33Z lovato $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trc_alloc()
      !!-------------------------------------------------------------------
      !!                    *** ROUTINE trc_alloc ***
      !!-------------------------------------------------------------------
!      USE lib_mpp, ONLY: ctl_stop
      !!-------------------------------------------------------------------
      INTEGER :: ierr(4)
      !!-------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( neln(Istrp:Iendp,Jstrp:Jendp), heup(Istrp:Iendp,Jstrp:Jendp), heup_01(Istrp:Iendp,Jstrp:Jendp), &
          &    cvol(Istrp:Iendp,Jstrp:Jendp,N),  etot(Istrp:Iendp,Jstrp:Jendp,N), &
          &    etot_ndcy(Istrp:Iendp,Jstrp:Jendp,N), STAT = ierr(1)  )
      ! 
      trc_alloc = MAXVAL( ierr )
      IF( trc_alloc /= 0 )   CALL ctl_stop( 'STOP', 'trc_alloc: failed to allocate arrays' )
      !
   END FUNCTION trc_alloc


   !!======================================================================
END MODULE trc
