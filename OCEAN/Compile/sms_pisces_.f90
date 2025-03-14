










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









MODULE sms_pisces   
   !!----------------------------------------------------------------------
   !!                     ***  sms_pisces.F90  ***  
   !! TOP :    Source Minus Sink variables are declared and allocated
   !!----------------------------------------------------------------------
   !! History :   1.0  !  2000-02 (O. Aumont) original code
   !!             3.2  !  2009-04 (C. Ethe & NEMO team) style
   !!----------------------------------------------------------------------

!   USE par_oce
!   USE par_trc
    USE oce_trc

   IMPLICIT NONE
   PUBLIC

!   CHARACTER(50) :: numnatp_ref   !! Character buffer for reference namelist pisces
!   CHARACTER(50) :: numnatp_cfg   !! Character buffer for configuration namelist pisces
   INTEGER ::   numonp      = -1                !! Logical unit for namelist pisces output
   INTEGER ::   numnatp_ref = -1                !! Logical unit for namelist pisces output
   INTEGER ::   numnatp_cfg = -1                !! Logical unit for namelist pisces output

   !!* Model used
   LOGICAL  ::  ln_p2z            !: Flag to use   reduced model
   LOGICAL  ::  ln_p4z            !: Flag to use   model
   LOGICAL  ::  ln_p5z            !: Flag to use   quota model
   LOGICAL  ::  ln_ligand         !: Flag to enable organic ligands
   LOGICAL  ::  ln_sediment       !: Flag to enable sediment module

   !!*  Time variables
   INTEGER  ::   nrdttrc          !: ???
   REAL(wp) ::   rfact , rfactr   !: time step duration (in seconds)
   REAL(wp) ::   rfact2, rfact2r  !: time step duration (in seconds) when timesplitting is activated for 
   REAL(wp) ::   xstep            !: Time step duration for biology
!   REAL(wp) ::   ryyss            !: number of seconds per year 
   REAL(wp) ::   r1_ryyss         !: inverse number of seconds per year 

   !!*  Biological parameters 
   REAL(wp) ::   rno3             !: C/N stoichiometric ratio
   REAL(wp) ::   o2ut             !: O2/N stoichiometric ratio for ammonification
   REAL(wp) ::   po4r             !: C/P stoichiometric ratio
   REAL(wp) ::   rdenit           !: C/N ratio for denitrification
   REAL(wp) ::   rdenita          !: C/N ratio for denitrification
   REAL(wp) ::   o2nit            !: O2/N ratio for nitrification
   REAL(wp) ::   wsbio, wsbio2    !: Sinking speeds of particles
   REAL(wp) ::   wsbio2max        !: Maximum sinking speed of the largest particles
   REAL(wp) ::   wsbio2scale      !: Length scale for the variations of wsbio2
   REAL(wp) ::   oxymin           !:  half saturation constant for anoxia
   REAL(wp) ::   xkmort           !: Mortality half-saturation constant
   REAL(wp) ::   feratz           !: Fe/C in microzooplankton
   REAL(wp) ::   feratm           !: Fe/C in mesozooplankton
   REAL(wp) ::   ldocp            !: Ligand production ratio during PP
   REAL(wp) ::   ldocz            !: Ligand production ratio by grazing
   REAL(wp) ::   lthet            !: Uptake of ligand by phytoplankton
   REAL(wp) ::   no3rat3          !: C/N ratio of zooplankton
   REAL(wp) ::   po4rat3          !: C/P ratio of zooplankton

   !!*  diagnostic parameters 
   REAL(wp) ::  tpp               !: total primary production
   REAL(wp) ::  t_oce_co2_exp     !: total carbon export
   REAL(wp) ::  t_oce_co2_flx     !: Total ocean carbon flux
   REAL(wp) ::  t_oce_co2_flx_cum !: Cumulative Total ocean carbon flux
   REAL(wp) ::  t_atm_co2_flx     !: global mean of atmospheric pco2

   !!* restoring
   LOGICAL  ::  ln_pisdmp         !: restoring or not of nutrients to a mean value
   LOGICAL  ::  ln_pisdmp_alk     !: restoring or not of alkalinity to a mean value

   LOGICAL, PUBLIC ::   ln_ironice   !: boolean for Fe input from sea ice

   !!* Diurnal cycle in 
   LOGICAL  ::  ln_p4z_dcyc       !: Flag to activate diurnal cycle in 

   !!*  Biological fluxes for light : variables shared by pisces & lobster
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:  ) ::  strn  !: Day duration in hours
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  enano, ediat   !: PAR for phyto, nano and diat 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  enanom, ediatm !: mean PAR for phyto, nano and diat 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  epico          !: PAR for pico
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  epicom         !: mean PAR for pico
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  emoy, etotm    !: averaged PAR in the mixed layer
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::  xksi  !:  Half-saturation con,stant for diatoms

   !!*  Biological fluxes for primary production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)    ::   xksimax    !: Maximum half-saturation constant over the year (Si)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   biron      !: bioavailable fraction of iron
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   plig       !: proportion of iron organically complexed

   !!*  Sinking speed
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio3   !: POC sinking speed 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   wsbio4   !: GOC sinking speed

   !!*  SMS for the organic matter
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xfracal    !: Fraction of nanophytoplankton that are calcifying organisms
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   nitrfac    !: OMZ 
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   nitrfac2   !: N depleted indice
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   orem       !: oxic remineralisation
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xdiss      !: Shear rate
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prodcal    !: Calcite production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prodpoc    !: POC production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   conspoc    !: POC consumption
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   prodgoc    !: GOC production
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   consgoc    !: GOC consumption
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   consfe3    !: GOC consumption
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   blim       !: bacterial production factor
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sizen, logsizen    !: size of nanophyto
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sizep, logsizep    !: size of picophyto
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sized, logsized    !: size of diatoms
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sizena     !: size of nanophytoplankton, after
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sizepa     !: size of picophyto, after
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sizeda     !: size of diatomss, after
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   thetanano  !: size of diatomss, after
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xfecolagg  !: Refractory diagnostic concentration of ligands
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   xcoagfe    !: Coagulation rate of colloidal Fe/ligands

   !!* Variable for chemistry of the CO2 cycle
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak13       !: Carbonate chemistry constant
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak23       !: Carbonate chemistry constant
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aksp       !: Solubility product of CaCO3
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   hi         !: Proton concentration
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   excess     !: CO3 saturation
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aphscale   !: 

   !!* Temperature dependancy of SMS terms
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfunc    !: Temp. dependancy of various biological rates
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tgfunc2   !: Temp. dependancy of mesozooplankton rates

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   remintpoc  ! Mean remineralisation rate of POC
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   remintgoc  ! Mean remineralisation rate of GOC

   LOGICAL, SAVE :: lk_sed
   LOGICAL, SAVE :: l_diaadd

   LOGICAL, SAVE :: ln_bio, ln_lys, ln_sed, ln_flx, ln_sink
   LOGICAL, SAVE :: ln_fechem, ln_micro, ln_meso, ln_mort
   LOGICAL, SAVE :: ln_prod, ln_agg, ln_rem, ln_poc, ln_diaz

!! * Substitutions





























   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: sms_pisces.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sms_pisces_alloc()
      !!----------------------------------------------------------------------
      !!        *** ROUTINE sms_pisces_alloc ***
      !!----------------------------------------------------------------------
      INTEGER ::   ierr(19)        ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !*  Biological fluxes for light : shared variables for pisces & lobster
      ALLOCATE( strn(Istrp:Iendp,Jstrp:Jendp),  STAT=ierr(1) )

      !* Optics
      ALLOCATE(  enano(Istrp:Iendp,Jstrp:Jendp,N) , enanom(Istrp:Iendp,Jstrp:Jendp,N) ,   &
         &       emoy(Istrp:Iendp,Jstrp:Jendp,N)  , etotm(Istrp:Iendp,Jstrp:Jendp,N)  ,      STAT=ierr(2) )

      ! Biological SMS
      ALLOCATE( orem     (Istrp:Iendp,Jstrp:Jendp,N), xdiss   (Istrp:Iendp,Jstrp:Jendp,N),  &
         &      nitrfac  (Istrp:Iendp,Jstrp:Jendp,N), nitrfac2(Istrp:Iendp,Jstrp:Jendp,N),  &
         &      prodcal  (Istrp:Iendp,Jstrp:Jendp,N), prodpoc (Istrp:Iendp,Jstrp:Jendp,N),  &
         &      conspoc  (Istrp:Iendp,Jstrp:Jendp,N), xfracal (Istrp:Iendp,Jstrp:Jendp,N),   STAT=ierr(3) )

      !* Carbonate chemistry
      ALLOCATE( ak13(Istrp:Iendp,Jstrp:Jendp,N),                         &
         &      ak23(Istrp:Iendp,Jstrp:Jendp,N), aksp  (Istrp:Iendp,Jstrp:Jendp,N) ,    &
         &      hi  (Istrp:Iendp,Jstrp:Jendp,N), excess(Istrp:Iendp,Jstrp:Jendp,N) ,    &
         &      aphscale(Istrp:Iendp,Jstrp:Jendp,N),                          STAT=ierr(4) )
      !
      !* Temperature dependency of SMS terms
      ALLOCATE( tgfunc (Istrp:Iendp,Jstrp:Jendp,N) , tgfunc2(Istrp:Iendp,Jstrp:Jendp,N),     STAT=ierr(5) )
      !
      !* Sinking speed
      ALLOCATE( wsbio3 (Istrp:Iendp,Jstrp:Jendp,N) , wsbio4 (Istrp:Iendp,Jstrp:Jendp,N),     STAT=ierr(6) )

      !*  Size of phytoplankton cells
      ALLOCATE( sizen (Istrp:Iendp,Jstrp:Jendp,N), sizena(Istrp:Iendp,Jstrp:Jendp,N),       &
         &      logsizen(Istrp:Iendp,Jstrp:Jendp,N),                          STAT=ierr(7) )

      ALLOCATE( blim     (Istrp:Iendp,Jstrp:Jendp,N), consfe3 (Istrp:Iendp,Jstrp:Jendp,N),  &
         &      xfecolagg(Istrp:Iendp,Jstrp:Jendp,N), xcoagfe (Istrp:Iendp,Jstrp:Jendp,N),   STAT=ierr(8) )
      !
      ALLOCATE( plig(Istrp:Iendp,Jstrp:Jendp,N)  ,   biron(Istrp:Iendp,Jstrp:Jendp,N)    ,   STAT=ierr(9) )

      ALLOCATE( remintpoc(Istrp:Iendp,Jstrp:Jendp,N),                         STAT=ierr(10) )

      IF( ln_p2z )   &
         &   ALLOCATE( thetanano (Istrp:Iendp,Jstrp:Jendp,N),                 STAT=ierr(10) )

      IF( ln_p4z .OR. ln_p5z ) THEN
         !* Optics
         ALLOCATE(  ediat(Istrp:Iendp,Jstrp:Jendp,N) , ediatm(Istrp:Iendp,Jstrp:Jendp,N),    STAT=ierr(11) )

         !* Biological SMS
         ALLOCATE( xksimax(Istrp:Iendp,Jstrp:Jendp)  ,                          STAT=ierr(12) )

         ! Biological SMS
         ALLOCATE( prodgoc(Istrp:Iendp,Jstrp:Jendp,N), consgoc(Istrp:Iendp,Jstrp:Jendp,N),   STAT=ierr(13) )
         !
         !* Si 1/2 saturation constant 
         ALLOCATE( xksi (Istrp:Iendp,Jstrp:Jendp)  ,                            STAT=ierr(14) )

         !*  Size of phytoplankton cells
         ALLOCATE( sized (Istrp:Iendp,Jstrp:Jendp,N), sizeda(Istrp:Iendp,Jstrp:Jendp,N),     &
            &      logsized(Istrp:Iendp,Jstrp:Jendp,N),                       STAT=ierr(16) )

         ALLOCATE( remintgoc(Istrp:Iendp,Jstrp:Jendp,N),                      STAT=ierr(17) )
         ! 
      ENDIF
      !
      IF( ln_p5z ) THEN
         ! -QUOTA specific part      
         ALLOCATE( epico(Istrp:Iendp,Jstrp:Jendp,N)   , epicom(Istrp:Iendp,Jstrp:Jendp,N),   STAT=ierr(16) )

        !*  Size of phytoplankton cells
         ALLOCATE( sizep(Istrp:Iendp,Jstrp:Jendp,N), sizepa(Istrp:Iendp,Jstrp:Jendp,N),      &
            &      logsizep(Istrp:Iendp,Jstrp:Jendp,N),                       STAT=ierr(19) )
      ENDIF
      !
      sms_pisces_alloc = MAXVAL( ierr )
      !
      IF( sms_pisces_alloc /= 0 )   CALL ctl_stop( 'STOP', 'sms_pisces_alloc: failed to allocate arrays' ) 
      !
   END FUNCTION sms_pisces_alloc


   !!======================================================================   
END MODULE sms_pisces    
