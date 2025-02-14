










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









MODULE p4zmeso
   !!======================================================================
   !!                         ***  MODULE p4zmeso  ***
   !! TOP :    Compute the sources/sinks for mesozooplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!             3.*  !  2025-02  (S. Maishal, R. Person)  and optimization
   !!----------------------------------------------------------------------
   !!   p4z_meso        : Compute the sources/sinks for mesozooplankton
   !!   p4z_meso_init   : Initialization of the parameters for mesozooplankton
   !!   p4z_meso_alloc  : Allocate variables for mesozooplankton 
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE p4zprod         ! production
   USE prtctl          ! print control for debugging
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_meso              ! called in p4zbio.F90
   PUBLIC   p4z_meso_init         ! called in trcsms_pisces.F90
   PUBLIC   p4z_meso_alloc        ! called in trcini_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  part2        !: part of calcite not dissolved in mesozoo guts
   REAL(wp), PUBLIC ::  xpref2d      !: mesozoo preference for diatoms
   REAL(wp), PUBLIC ::  xpref2n      !: mesozoo preference for nanophyto
   REAL(wp), PUBLIC ::  xpref2z      !: mesozoo preference for microzooplankton
   REAL(wp), PUBLIC ::  xpref2c      !: mesozoo preference for POC 
   REAL(wp), PUBLIC ::  xpref2m      !: mesozoo preference for mesozoo
   REAL(wp), PUBLIC ::  xthresh2zoo  !: zoo feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2dia  !: diatoms feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2phy  !: nanophyto feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2poc  !: poc feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2mes  !: mesozoo feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2     !: feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  resrat2      !: exsudation rate of mesozooplankton
   REAL(wp), PUBLIC ::  mzrat2       !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::  grazrat2     !: maximal mesozoo grazing rate
   REAL(wp), PUBLIC ::  xkgraz2      !: non assimilated fraction of P by mesozoo 
   REAL(wp), PUBLIC ::  unass2       !: Efficicency of mesozoo growth 
   REAL(wp), PUBLIC ::  sigma2       !: Fraction of mesozoo excretion as DOM 
   REAL(wp), PUBLIC ::  epsher2      !: growth efficiency
   REAL(wp), PUBLIC ::  epsher2min   !: minimum growth efficiency at high food for grazing 2
   REAL(wp), PUBLIC ::  xsigma2      !: Width of the predation window
   REAL(wp), PUBLIC ::  xsigma2del   !: Maximum width of the predation window at low food density
   REAL(wp), PUBLIC ::  grazflux     !: mesozoo flux feeding rate
   REAL(wp), PUBLIC ::  xfracmig     !: Fractional biomass of meso that performs DVM
   LOGICAL , PUBLIC ::  ln_dvm_meso  !: Boolean to activate DVM of mesozooplankton
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: depmig  !: DVM of mesozooplankton : migration depth
   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:) :: kmig    !: Vertical indice of the the migration depth

   REAL(wp)         ::  xfracmigm1     !: Fractional biomass of meso that performs DVM
   LOGICAL          :: l_dia_graz, l_dia_lprodz

   !! * Substitutions












































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zmeso.F90 15482 2021-11-08 20:02:22Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_meso( kt, knt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_meso  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for mesozooplankton
      !!                This includes ingestion and assimilation, flux feeding
      !!                and mortality. We use a passive prey switching  
      !!                parameterization.
      !!                All living compartments smaller than mesozooplankton
      !!                are potential preys of mesozooplankton as well as small
      !!                sinking particles 
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step and ???
      INTEGER, INTENT(in)  ::  Kbb, kmm, Krhs ! time level indices
      !
      INTEGER  :: mi, mj, jk, jkt
      REAL(wp) :: zcompadi, zcompaph, zcompapoc, zcompaz, zcompam, zcompames
      REAL(wp) :: zgraze2, zdenom, zdenom2, zfact, zfood, zfoodlim, zproport, zbeta
      REAL(wp) :: zmortzgoc, zfrac, zfracfe, zratio, zratio2, zfracal, zgrazcal
      REAL(wp) :: zepsherf, zepshert, zepsherq, zepsherv, zgraztotc, zgraztotn, zgraztotf
      REAL(wp) :: zmigreltime, zprcaca, zmortz, zgrasratf, zgrasratn
      REAL(wp) :: zrespz, ztortz, zgrazdc, zgrazz, zgrazpof, zgraznc, zgrazpoc, zgraznf, zgrazdf
      REAL(wp) :: zgrazm, zgrazfffp, zgrazfffg, zgrazffep, zgrazffeg, zdep
      REAL(wp) :: zsigma, zsigma2, zsizedn, zdiffdn, ztmp1, ztmp2, ztmp3, ztmp4, ztmp5, ztmptot, zmigthick 
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N) :: zgrarem, zgraref, zgrapoc, zgrapof, zgrabsi
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   zgramigrem, zgramigref, zgramigpoc, zgramigpof
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   zgramigbsi
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zgrazing2, zzligprod, zw3d
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p4z_meso')
      !
      IF( kt == nittrc000 )  THEN
         l_dia_graz    = iom_use( "GRAZ2" ) .OR. iom_use( "FEZOO2" ) .OR. iom_use( "MesoZo2" )
         l_dia_lprodz  = ln_ligand .AND. iom_use( "LPRODZ2" )
         l_dia_graz = l_dia_graz .OR. l_diaadd
      ENDIF
      IF( l_dia_graz ) THEN
         ALLOCATE( zgrazing2(Istrp:Iendp,Jstrp:Jendp,N) )     ;    zgrazing2(Istrp:Iendp,Jstrp:Jendp,N) = 0.
      ENDIF
      !
      IF( l_dia_lprodz ) THEN
         ALLOCATE( zzligprod(Istrp:Iendp,Jstrp:Jendp,N) )
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk) &
         !$OMP SHARED(zzligprod, tr, jplgw, Krhs, N)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zzligprod(mi,mj,jk) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jplgw)
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
      ENDIF
      !
      zgrapoc(:,:,:) = 0._wp    ;  zgrarem(:,:,:) = 0._wp
      zgraref (:,:,:) = 0._wp   ;  zgrapof(:,:,:) = 0._wp
      zgrabsi (:,:,:) = 0._wp
      !
      !
      ! Diurnal vertical migration of mesozooplankton
      ! Computation of the migration depth
      ! ---------------------------------------------
      IF (ln_dvm_meso) CALL p4z_meso_depmig( Kbb, Kmm )
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zcompam, zfact, zrespz, ztortz, zcompadi, &
                    zcompaz, zcompapoc, zcompames, zcompaph, zfood, zfoodlim, &
                    zdenom, zgraze2, zdenom2, zsigma, zsigma2, zsizedn, zdiffdn, &
                    ztmp1, ztmp2, ztmp3, ztmp4, ztmp5, ztmptot, zgrazdc, zgraznc, &
                    zgrazpoc, zgrazz, zgrazm, zratio, zgraznf, zgrazdf, zgrazpof, &
                    zgrazffeg, zgrazfffg, zgrazffep, zgrazfffp, zgraztotc, zproport, &
                    zfrac, zfracfe, zgrazcal, zprcaca, zgraztotn, zgraztotf, zgrasratf, &
                    zgrasratn, zepshert, zbeta, zepsherf, zepsherq, zepsherv, zmortz, &
                    zmortzgoc, zgrarem, zgraref, zgrapoc, zgrapof) &
      !$OMP SHARED(tr, nitrfac, tgfunc2, xstep, resrat2, xkmort, mzrat2, xpref2d, &
                    xpref2z, xpref2n, xpref2c, xpref2m, xthresh2, grazrat2, xkgraz2, &
                    xsigma2, xsigma2del, 0.5*EPSILON(1.e0), feratm, epsher2, epsher2min, unass2, &
                    part2, prodcal, N)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zcompam   = MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpmes) - 1.e-9 ), 0.e0 )
         zfact     = xstep * tgfunc2(mi,mj,jk) * zcompam

         
         !  linear mortality of mesozooplankton
         !  A michaelis menten modulation term is used to avoid extinction of 
         !  mesozooplankton at very low food concentration. Mortality is

         !  enhanced in low O2 waters
         !  -----------------------------------------------------------------
         zrespz    = resrat2 * zfact * ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpmes) &
           &         / ( xkmort + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpmes) )  &
           &           + 3. * nitrfac(mi,mj,jk) )

         ! Zooplankton quadratic mortality. A square function has been selected with
         !  to mimic predation and disease (density dependent mortality). It also tends
         !  to stabilise the model
         !  -------------------------------------------------------------------------
         ztortz    = mzrat2 * 1.e6 * zfact * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpmes)  &
           &       * (1. - nitrfac(mi,mj,jk) )
         !
         !   Computation of the abundance of the preys
         !   A threshold can be specified in the namelist
         !   --------------------------------------------
         zcompadi  = MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) - xthresh2dia ), 0.e0 )
         zcompaz   = MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo) - xthresh2zoo ), 0.e0 )
         zcompapoc = MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) - xthresh2poc ), 0.e0 )
         zcompames = MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpmes) - xthresh2mes ), 0.e0 )
         zcompaph  = MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) - xthresh2phy ), 0.e0 )

         ! Mesozooplankton grazing
         ! The total amount of food is the sum of all preys accessible to mesozooplankton 
         ! multiplied by their food preference
         ! A threshold can be specified in the namelist (xthresh2). However, when food 
         ! concentration is close to this threshold, it is decreased to avoid the 
         ! accumulation of food in the mesozoopelagic domain
         ! -------------------------------------------------------------------------------
         zfood     = xpref2d * zcompadi + xpref2z * zcompaz &
           &         + xpref2n * zcompaph + xpref2c * zcompapoc &
           &         + xpref2m * zcompames
         zfoodlim  = MAX( 0., zfood - MIN( 0.5 * zfood, xthresh2 ) )
         zdenom    = zfoodlim / ( xkgraz2 + zfoodlim )
         zgraze2   = grazrat2 * xstep * tgfunc2(mi,mj,jk) &
           &       * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpmes) * (1. - nitrfac(mi,mj,jk))

         ! An active switching parameterization is used here.
         ! We don't use the KTW parameterization proposed by 
         ! Vallina et al. because it tends to produce too steady biomass
         ! composition and the variance of Chl is too low as it grazes
         ! too strongly on winning organisms. We use a generalized
         ! switching parameterization proposed by Morozov and Petrovskii (2013)
         ! ------------------------------------------------------------  
         ! The width of the selection window is increased when preys
         ! have low abundance, .i.e. zooplankton become less specific 
         ! to avoid starvation.
         ! ----------------------------------------------------------
         zdenom2 = zdenom * zdenom
         zsigma  = 1.0 - zdenom2/(0.05*0.05+zdenom2)
         zsigma  = xsigma2 + xsigma2del * zsigma
         zsigma2 = zsigma * zsigma
         ! Nanophytoplankton and diatoms are the only preys considered
         ! to be close enough to have potential interference
         ! -----------------------------------------------------------
         zsizedn = ABS(LOG(1.67 * sizen(mi,mj,jk) / (5.0 * sized(mi,mj,jk) + 0.5*EPSILON(1.e0) )) )
         zdiffdn = EXP( -zsizedn * zsizedn / zsigma2 )
         ztmp1 = xpref2n * zcompaph * ( zcompaph + zdiffdn * zcompadi )
         ztmp2 = xpref2m * zcompames * zcompames
         ztmp3 = xpref2c * zcompapoc * zcompapoc
         ztmp4 = xpref2d * zcompadi * ( zcompadi + zdiffdn * zcompaph )
         ztmp5 = xpref2z * zcompaz * zcompaz
         ztmptot = ztmp1 + ztmp2 + ztmp3 + ztmp4 + ztmp5 + 0.5*EPSILON(1.e0)
         ztmp1 = ztmp1 / ztmptot
         ztmp2 = ztmp2 / ztmptot
         ztmp3 = ztmp3 / ztmptot
         ztmp4 = ztmp4 / ztmptot
         ztmp5 = ztmp5 / ztmptot

         !   Mesozooplankton regular grazing on the different preys
         !   ------------------------------------------------------
         zgrazdc   = zgraze2 * ztmp4 * zdenom  ! diatoms
         zgraznc   = zgraze2 * ztmp1 * zdenom  ! nanophytoplankton
         zgrazpoc  = zgraze2 * ztmp3 * zdenom  ! small POC
         zgrazz    = zgraze2 * ztmp5 * zdenom  ! microzooplankton
         zgrazm    = zgraze2 * ztmp2 * zdenom  ! sm

         ! Ingestion rates of the Fe content of the different preys
         zratio = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnfe) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) + 0.5*EPSILON(1.e0))
         zgraznf   = zgraznc  * zratio 
         zratio = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdfe) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) + 0.5*EPSILON(1.e0))
         zgrazdf   = zgrazdc  * zratio
         zratio = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsfe) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) + 0.5*EPSILON(1.e0))
         zgrazpof  = zgrazpoc *  zratio

         !  Mesozooplankton flux feeding on GOC and POC. The feeding pressure
         ! is proportional to the flux
         !  ------------------------------------------------------------------
         zgrazffeg = grazflux  * xstep * wsbio4(mi,mj,jk) &
         &           * tgfunc2(mi,mj,jk) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc) &
         &           * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpmes) &
         &           * (1. - nitrfac(mi,mj,jk))
         zratio = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpbfe) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc) + 0.5*EPSILON(1.e0))
         zgrazfffg = zgrazffeg * zratio
         zgrazffep = grazflux  * xstep *  wsbio3(mi,mj,jk) &
         &           * tgfunc2(mi,mj,jk) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) &
         &           * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpmes) &
         &           * (1. - nitrfac(mi,mj,jk))
         zratio = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsfe) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) + 0.5*EPSILON(1.e0))
         zgrazfffp = zgrazffep * zratio
         !
         zgraztotc = zgrazdc + zgrazz + zgraznc &
             &     + zgrazm + zgrazpoc + zgrazffep + zgrazffeg

         ! Compute the proportion of filter feeders. It is assumed steady state.
         ! ---------------------------------------------------------------------
         zproport  = (zgrazffep + zgrazffeg)/(0.5*EPSILON(1.e0) + zgraztotc)
         zproport = zproport**2

         ! Compute fractionation of aggregates. It is assumed that 
         ! diatoms based aggregates are more prone to fractionation
         ! since they are more porous (marine snow instead of fecal pellets)
         ! -----------------------------------------------------------------

         ! Compute fractionation of aggregates. It is assumed that 
         ! diatoms based aggregates are more prone to fractionation
         ! since they are more porous (marine snow instead of fecal pellets)
         zratio    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgsi) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc) + 0.5*EPSILON(1.e0) )
         zratio2   = zratio * zratio
         zfrac     = zproport * grazflux  * xstep * wsbio4(mi,mj,jk) &
         &          * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpmes) &
         &          * ( 0.4 + 3.6 * zratio2 / ( 1.**2 + zratio2 ) )
         zratio    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpbfe) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc) + 0.5*EPSILON(1.e0) )
         zfracfe   = zfrac * zratio

         ! Flux feeding is multiplied by the fractional biomass of flux feeders
         zgrazffep = zproport * zgrazffep
         zgrazffeg = zproport * zgrazffeg
         zgrazfffp = zproport * zgrazfffp
         zgrazfffg = zproport * zgrazfffg
         zgrazdc   = (1.0 - zproport) * zgrazdc
         zgraznc   = (1.0 - zproport) * zgraznc
         zgrazz    = (1.0 - zproport) * zgrazz
         zgrazpoc  = (1.0 - zproport) * zgrazpoc
         zgrazm    = (1.0 - zproport) * zgrazm
         zgrazdf   = (1.0 - zproport) * zgrazdf
         zgraznf   = (1.0 - zproport) * zgraznf
         zgrazpof  = (1.0 - zproport) * zgrazpof

         ! Total ingestion rates in C, N, Fe
         zgraztotc = zgrazdc + zgrazz + zgraznc + zgrazpoc &
              &    + zgrazm + zgrazffep + zgrazffeg ! grazing by mesozooplankton
         IF( l_dia_graz ) zgrazing2(mi,mj,jk) = zgraztotc

         zgraztotn = zgrazdc * quotad(mi,mj,jk) + zgrazz + zgraznc * quotan(mi,mj,jk) &
         &          + zgrazm + zgrazpoc + zgrazffep + zgrazffeg
         zgraztotf = zgrazdf + zgraznf + zgrazz * feratz &
             &     + zgrazpof + zgrazfffp + zgrazfffg + zgrazm * feratm

         !   Stoichiometruc ratios of the food ingested by zooplanton 
         !   --------------------------------------------------------
         zgrasratf =  ( zgraztotf + 0.5*EPSILON(1.e0) )/ ( zgraztotc + 0.5*EPSILON(1.e0) )
         zgrasratn =  ( zgraztotn + 0.5*EPSILON(1.e0) )/ ( zgraztotc + 0.5*EPSILON(1.e0) )

         ! Mesozooplankton efficiency. 
         ! We adopt a formulation proposed by Mitra et al. (2007)
         ! The gross growth efficiency is controled by the most limiting nutrient.
         ! Growth is also further decreased when the food quality is poor. This is currently
         ! hard coded : it can be decreased by up to 50% (zepsherq)
         ! GGE can also be decreased when food quantity is high, zepsherf (Montagnes and Fulton, 2012)
         ! -----------------------------------------------------------------------------------
         zepshert  = MIN( 1., zgrasratn, zgrasratf / feratm)
         zbeta     = MAX(0., (epsher2 - epsher2min) )
         ! Food quantity deprivation of GGE
         zepsherf  = epsher2min + zbeta / ( 1.0 + 0.04E6 * 12. * zfood * zbeta )
         ! Food quality deprivation of GGE
         zepsherq  = 0.5 + (1.0 - 0.5) * zepshert * ( 1.0 + 1.0 ) / ( zepshert + 1.0 )
         ! Actual GGE
         zepsherv  = zepsherf * zepshert * zepsherq
         ! 
         ! Impact of grazing on the prognostic variables
         ! ---------------------------------------------
         zmortz = ztortz + zrespz
         ! Mortality induced by the upper trophic levels, ztortz, is allocated 
         ! according to a infinite chain of predators (Anderson et al., 2013)
         zmortzgoc = unass2 / ( 1. - epsher2 ) * ztortz + zrespz

         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpmes) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpmes) - zmortz &
                &                + zepsherv * zgraztotc - zgrazm 
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdia) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdia) - zgrazdc
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpzoo) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpzoo) - zgrazz
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpphy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpphy) - zgraznc
         zratio    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnch) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) + 0.5*EPSILON(1.e0) )
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnch) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnch) - zgraznc * zratio
         zratio    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdch) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) + 0.5*EPSILON(1.e0) )
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdch) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdch) - zgrazdc * zratio
         zratio    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdsi) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) + 0.5*EPSILON(1.e0) )
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdsi) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdsi) - zgrazdc * zratio
         zgrabsi(mi,mj,jk)       = zgrazdc * zratio
         !
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnfe) - zgraznf
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdfe) - zgrazdf
         !
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) - zgrazpoc - zgrazffep + zfrac
         prodpoc(mi,mj,jk) = prodpoc(mi,mj,jk) + zfrac
         conspoc(mi,mj,jk) = conspoc(mi,mj,jk) - zgrazpoc - zgrazffep
         !
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgoc) - zgrazffeg - zfrac
         consgoc(mi,mj,jk) = consgoc(mi,mj,jk) - zgrazffeg - zfrac
         !
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) - zgrazpof - zgrazfffp + zfracfe
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) - zgrazfffg - zfracfe

         ! Calcite remineralization due to zooplankton activity
         ! part2 of the ingested calcite is not dissolving in the 
         ! acidic gut
         ! ------------------------------------------------------
         zfracal = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpcal) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc) + 0.5*EPSILON(1.e0) )
         zgrazcal = zgrazffeg * (1. - part2) * zfracal
         ! calcite production by zooplankton activity
         zprcaca = xfracal(mi,mj,jk) * zgraznc
         prodcal(mi,mj,jk) = prodcal(mi,mj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
         !
         zprcaca = part2 * zprcaca
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) + zgrazcal - zprcaca
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) + 2. * ( zgrazcal - zprcaca )
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpcal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpcal) - zgrazcal + zprcaca

         ! Computation of total excretion and egestion by mesozoo.
         ! ---------------------------------------------------------
         zgrarem(mi,mj,jk) = zgraztotc * ( 1. - zepsherv - unass2 ) &
                 &         + ( 1. - epsher2 - unass2 ) / ( 1. - epsher2 ) * ztortz
         zgraref(mi,mj,jk) = zgraztotc * MAX( 0. , ( 1. - unass2 ) &
                 &                 * zgrasratf - feratm * zepsherv ) &
                 &         + feratm * ( ( 1. - epsher2 - unass2 ) &
                 &         /( 1. - epsher2 ) * ztortz )
         zgrapoc(mi,mj,jk) = zgraztotc * unass2 + zmortzgoc
         zgrapof(mi,mj,jk) = zgraztotf * unass2 + feratm * zmortzgoc
         !
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      ! Computation of the effect of DVM by mesozooplankton
      ! This part is only activated if ln_dvm_meso is set to true
      ! The parameterization has been published in Gorgues et al. (2019).
      ! -----------------------------------------------------------------
      IF (ln_dvm_meso) THEN
         ALLOCATE( zgramigrem(Istrp:Iendp,Jstrp:Jendp), zgramigref(Istrp:Iendp,Jstrp:Jendp), &
              &    zgramigpoc(Istrp:Iendp,Jstrp:Jendp), zgramigpof(Istrp:Iendp,Jstrp:Jendp) )
         ALLOCATE( zgramigbsi(Istrp:Iendp,Jstrp:Jendp) )
         zgramigrem(:,:) = 0.0    ;   zgramigref(:,:) = 0.0
         zgramigpoc(:,:) = 0.0    ;   zgramigpof(:,:) = 0.0
         zgramigbsi(:,:) = 0.0

        ! Compute the amount of materials that will go into vertical migration
        ! This fraction is sumed over the euphotic zone and is removed from 
        ! the fluxes driven by mesozooplankton in the euphotic zone.
        ! --------------------------------------------------------------------
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(mi, mj, jk, zmigreltime, zmigthick) &
        !$OMP SHARED(zgrarem, zgraref, zgrapoc, zgrapof, zgrabsi, zgramigrem, zgramigref, &
                     zgramigpoc, zgramigpof, zgramigbsi, strn, e3t, tmask, gdept, heup, &
                     xfracmig, xfracmigm1, N)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zmigreltime = (1. - strn(mi,mj) / 24.)
            zmigthick   = (1. - zmigreltime ) * Hz(mi,mj,N+1-jk) * tmask(mi,mj,jk)
            IF ( ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) <= heup(mi,mj) ) THEN
               zgramigrem(mi,mj) = zgramigrem(mi,mj) + xfracmig * zgrarem(mi,mj,jk) * zmigthick
               zgramigref(mi,mj) = zgramigref(mi,mj) + xfracmig * zgraref(mi,mj,jk) * zmigthick
               zgramigpoc(mi,mj) = zgramigpoc(mi,mj) + xfracmig * zgrapoc(mi,mj,jk) * zmigthick
               zgramigpof(mi,mj) = zgramigpof(mi,mj) + xfracmig * zgrapof(mi,mj,jk) * zmigthick
               zgramigbsi(mi,mj) = zgramigbsi(mi,mj) + xfracmig * zgrabsi(mi,mj,jk) * zmigthick

               zgrarem(mi,mj,jk) = zgrarem(mi,mj,jk) * ( xfracmigm1 + xfracmig * zmigreltime )
               zgraref(mi,mj,jk) = zgraref(mi,mj,jk) * ( xfracmigm1 + xfracmig * zmigreltime )
               zgrapoc(mi,mj,jk) = zgrapoc(mi,mj,jk) * ( xfracmigm1 + xfracmig * zmigreltime )
               zgrapof(mi,mj,jk) = zgrapof(mi,mj,jk) * ( xfracmigm1 + xfracmig * zmigreltime )
               zgrabsi(mi,mj,jk) = zgrabsi(mi,mj,jk) * ( xfracmigm1 + xfracmig * zmigreltime )
            ENDIF
         END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
         ! The inorganic and organic fluxes induced by migrating organisms are added at the 
         ! the migration depth (corresponding indice is set by kmig)
         ! --------------------------------------------------------------------------------
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jkt, zdep) &
         !$OMP SHARED(tmask, kmig, e3t, zgrarem, zgraref, zgrapoc, zgrapof, &
                      zgrabsi, zgramigrem, zgramigref, zgramigpoc, &
                      zgramigpof, zgramigbsi)
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            IF( tmask(mi,mj,1) == 1.) THEN
               jkt = kmig(mi,mj)
               zdep = 1. / Hz(mi,mj,N+1-jkt)
               zgrarem(mi,mj,jkt) = zgrarem(mi,mj,jkt) + zgramigrem(mi,mj) * zdep
               zgraref(mi,mj,jkt) = zgraref(mi,mj,jkt) + zgramigref(mi,mj) * zdep
               zgrapoc(mi,mj,jkt) = zgrapoc(mi,mj,jkt) + zgramigpoc(mi,mj) * zdep
               zgrapof(mi,mj,jkt) = zgrapof(mi,mj,jkt) + zgramigpof(mi,mj) * zdep
               zgrabsi(mi,mj,jkt) = zgrabsi(mi,mj,jkt) + zgramigbsi(mi,mj) * zdep
            ENDIF
         END DO   ;   END DO
         !$OMP END PARALLEL DO
         !
         ! Deallocate temporary variables
         ! ------------------------------
         DEALLOCATE( zgramigrem, zgramigref, zgramigpoc, zgramigpof, zgramigbsi )
      ! End of the ln_dvm_meso part
      ENDIF

      !   Update the arrays TRA which contain the biological sources and sinks
      !   This only concerns the variables which are affected by DVM (inorganic 
      !   nutrients, DOC agands, and particulate organic carbon).
      !$OMP PARALLEL DO &
      !$OMP DO PRIVATE(mi, mj, jk) &
      !$OMP SHARED(tr, zgrarem, zgraref, zgrapoc, zgrapof, zgrabsi, &
                   sigma2, ldocz, o2ut, rno3, ln_ligand)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) + zgrarem(mi,mj,jk) * sigma2
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnh4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnh4) + zgrarem(mi,mj,jk) * sigma2
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) &
               &                 + zgrarem(mi,mj,jk) * ( 1. - sigma2 )
         !
         IF( ln_ligand ) & 
          &  t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jplgw) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jplgw) &
          &                           + zgrarem(mi,mj,jk) * ( 1. - sigma2 ) * ldocz
         !
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) - o2ut * zgrarem(mi,mj,jk) * sigma2
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) + zgraref(mi,mj,jk)
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) + zgrarem(mi,mj,jk) * sigma2
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) + rno3 * zgrarem(mi,mj,jk) * sigma2
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgoc) + zgrapoc(mi,mj,jk)
         prodgoc(mi,mj,jk)   = prodgoc(mi,mj,jk)   + zgrapoc(mi,mj,jk)
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) + zgrapof(mi,mj,jk)
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgsi) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgsi) + zgrabsi(mi,mj,jk)
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      ! Write the output
      IF( .false. .AND. knt == nrdttrc ) THEN
        !
        IF( l_dia_graz ) THEN  !  
            ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(mi, mj, N+1-jk) &
            !$OMP SHARED(zw3d, zgrazing2, rfact2r, tmask)
            DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zw3d(mi,mj,N+1-jk) =  zgrazing2(mi,mj,jk) * 1.e+3 * rfact2r * tmask(mi,mj,jk)
            END DO   ;   END DO   ;   END DO
            !$OMP END PARALLEL DO
            CALL iom_put( "GRAZ2" , zw3d )  ! Total grazing of phyto by zooplankton
            CALL iom_put( "MesoZo2" , zw3d * ( 1. - epsher2 - unass2 ) * (-o2ut) * sigma2 ) ! o2 consumption by Microzoo
            !
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(mi, mj, N+1-jk) &
            !$OMP SHARED(zw3d, zgraref, rfact2r, tmask)
            DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zw3d(mi,mj,N+1-jk) =  zgraref(mi,mj,jk) * 1e9 * 1.e+3 * rfact2r * tmask(mi,mj,jk)
            END DO   ;   END DO   ;   END DO
            !$OMP END PARALLEL DO
            CALL iom_put( "FEZOO2", zw3d )  ! 
           DEALLOCATE( zw3d)
        ENDIF
        !
        IF( l_dia_lprodz ) THEN
            ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(mi,mj,N+1-jk) &
            !$OMP SHARED(zw3d,tr,zzligprod,rfact2r,tmask)
            DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zw3d(mi,mj,N+1-jk) = ( t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jplgw) - zzligprod(mi,mj,jk) ) &
                   &          * 1e9 * 1.e+3 * rfact2r * tmask(mi,mj,jk) ! conversion in nmol/m2/s
            END DO   ;   END DO   ;   END DO
            !$OMP END PARALLEL DO
           CALL iom_put( "LPRODZ2", zw3d )
           DEALLOCATE( zzligprod, zw3d)
        ENDIF
        !
      ENDIF
      !
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi,mj,N+1-jk) &
      !$OMP SHARED(trc3d,zgrazing2,epsher2,unass2,o2ut,sigma2,rfact2r,tmask)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         bioFlux(mi,mj,N+1-jk,Ngrapoc2) = zgrazing2(mi,mj,jk) * 1.e+3 &
                                     * rfact2r * tmask(mi,mj,jk) !  grazing of phyto by mesozoo
         bioFlux(mi,mj,N+1-jk,Nmeso2)   = zgrazing2(mi,mj,jk) * ( 1. - epsher2 - unass2 ) &
            &                        * (-o2ut) * sigma2 * 1.e+3 * rfact2r &
                                     * tmask(mi,mj,jk) ! o2 consumption by Mesozoo
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
     IF( l_dia_graz ) DEALLOCATE( zgrazing2 )
      
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('meso')")
        CALL prt_ctl_info( charout, cdcomp = 'top' )
  !      CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( .false. )   CALL timing_stop('p4z_meso')
      !
   END SUBROUTINE p4z_meso


   SUBROUTINE p4z_meso_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_meso_init  ***
      !!
      !! ** Purpose :   Initialization of mesozooplankton parameters
      !!
      !! ** Method  :   Read the namp4zmes namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampismes
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp4zmes/ part2, grazrat2, resrat2, mzrat2, xpref2n, xpref2d, xpref2z, &
         &                xpref2c, xpref2m, xthresh2dia, xthresh2phy, xthresh2zoo, xthresh2poc, xthresh2mes, &
         &                xthresh2, xkgraz2, epsher2, epsher2min, sigma2, unass2, grazflux, ln_dvm_meso, &
         &                xsigma2, xsigma2del, xfracmig
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*) 
         WRITE(stdout,*) 'p4z_meso_init : Initialization of mesozooplankton parameters'
         WRITE(stdout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND(numnatp_ref);READ(numnatp_ref,namp4zmes,IOSTAT=ios);CALL ctl_nam(ios,"namp4zmes (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,namp4zmes,IOSTAT=ios);CALL ctl_nam(ios,"namp4zmes (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE( numonp, namp4zmes )
      !
      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*) '   Namelist : namp4zmes'
         WRITE(stdout,*) '      part of calcite not dissolved in mesozoo guts  part2        =', part2
         WRITE(stdout,*) '      mesozoo preference for phyto                   xpref2n      =', xpref2n
         WRITE(stdout,*) '      mesozoo preference for diatoms                 xpref2d      =', xpref2d
         WRITE(stdout,*) '      mesozoo preference for zoo                     xpref2z      =', xpref2z
         WRITE(stdout,*) '      mesozoo preference for poc                     xpref2c      =', xpref2c
         WRITE(stdout,*) '      mesozoo preference for mesozoo                 xpref2m      = ', xpref2m
         WRITE(stdout,*) '      microzoo feeding threshold  for mesozoo        xthresh2zoo  =', xthresh2zoo
         WRITE(stdout,*) '      diatoms feeding threshold  for mesozoo         xthresh2dia  =', xthresh2dia
         WRITE(stdout,*) '      nanophyto feeding threshold for mesozoo        xthresh2phy  =', xthresh2phy
         WRITE(stdout,*) '      poc feeding threshold for mesozoo              xthresh2poc  =', xthresh2poc
         WRITE(stdout,*) '      mesozoo feeding threshold for mesozoo          xthresh2mes  = ', xthresh2mes
         WRITE(stdout,*) '      feeding threshold for mesozooplankton          xthresh2     =', xthresh2
         WRITE(stdout,*) '      exsudation rate of mesozooplankton             resrat2      =', resrat2
         WRITE(stdout,*) '      mesozooplankton mortality rate                 mzrat2       =', mzrat2
         WRITE(stdout,*) '      maximal mesozoo grazing rate                   grazrat2     =', grazrat2
         WRITE(stdout,*) '      mesozoo flux feeding rate                      grazflux     =', grazflux
         WRITE(stdout,*) '      non assimilated fraction of P by mesozoo       unass2       =', unass2
         WRITE(stdout,*) '      Efficiency of Mesozoo growth                   epsher2      =', epsher2
         WRITE(stdout,*) '      Minimum Efficiency of Mesozoo growth           epsher2min   =', epsher2min
         WRITE(stdout,*) '      Fraction of mesozoo excretion as DOM           sigma2       =', sigma2
         WRITE(stdout,*) '      half sturation constant for grazing 2          xkgraz2      =', xkgraz2
         WRITE(stdout,*) '      Width of the grazing window                     xsigma2     =', xsigma2
         WRITE(stdout,*) '      Maximum additional width of the grazing window  xsigma2del  =', xsigma2del
         WRITE(stdout,*) '      Diurnal vertical migration of mesozoo.         ln_dvm_meso  =', ln_dvm_meso
         WRITE(stdout,*) '      Fractional biomass of meso  that performs DVM  xfracmig     =', xfracmig
      ENDIF
      !
      xfracmigm1 = 1.0 - xfracmig
      !
   END SUBROUTINE p4z_meso_init

   SUBROUTINE p4z_meso_depmig( Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_meso_depmig  ***
      !!
      !! ** Purpose :   Computation the migration depth of mesozooplankton
      !!
      !! ** Method  :   Computes the DVM depth of mesozooplankton from oxygen
      !!      temperature and chlorophylle following the parameterization 
      !!      proposed by Bianchi et al. (2013)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)  ::  Kbb, kmm ! time level indices
      !
      INTEGER  :: mi, mj, jk, jkp1
      !
      REAL(wp) :: ztotchl, z1dep
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp) :: oxymoy, tempmoy, zdepmoy

      !!---------------------------------------------------------------------
      !
      IF( .false. )  CALL timing_start('p4z_meso_depmig')
      !
      oxymoy(:,:)  = 0.
      tempmoy(:,:) = 0.
      zdepmoy(:,:) = 0.
      depmig (:,:) = 5.
      kmig   (:,:) = 1
      !
      ! Compute the averaged values of oxygen, temperature over the domain 
      ! 150m to 500 m depth.
      ! ----------------------------------------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi,mj,jk) &
      !$OMP SHARED(oxymoy,tempmoy,zdepmoy,tr,ts,tmask,gdept,e3t)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF( tmask(mi,mj,jk) == 1.) THEN
            IF( ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) >= 150. .AND. ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) <= 500.) THEN
               oxymoy(mi,mj)  = oxymoy(mi,mj)  + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpoxy) * 1E6 * Hz(mi,mj,N+1-jk)
               tempmoy(mi,mj) = tempmoy(mi,mj) + t(mi,mj,N+1-jk,nnew,itemp) * Hz(mi,mj,N+1-jk)
               zdepmoy(mi,mj) = zdepmoy(mi,mj) + Hz(mi,mj,N+1-jk)
            ENDIF
         ENDIF
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      ! Compute the difference between surface values and the mean values in the mesopelagic
      ! domain
      ! ------------------------------------------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi,mj,z1dep) &
      !$OMP SHARED(oxymoy,tempmoy,tr,ts,zdepmoy,0.5*EPSILON(1.e0))
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         z1dep = 1. / ( zdepmoy(mi,mj) + 0.5*EPSILON(1.e0) )
         oxymoy(mi,mj)  = t(mi,mj,N+1-1,Kbb,itemp+ntrc_salt+jpoxy) * 1E6 - oxymoy(mi,mj)  * z1dep
         tempmoy(mi,mj) = t(mi,mj,N+1-1,nnew,itemp)      - tempmoy(mi,mj) * z1dep
      END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      ! Computation of the migration depth based on the parameterization of 
      ! Bianchi et al. (2013)
      ! -------------------------------------------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, ztotchl) &
      !$OMP SHARED(tmask, tr, oxymoy, depmig, hmld, tempmoy)
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF( tmask(mi,mj,1) == 1. ) THEN
            ztotchl = ( t(mi,mj,N+1-1,Kbb,itemp+ntrc_salt+jpnch) + t(mi,mj,N+1-1,Kbb,itemp+ntrc_salt+jpdch) ) * 1E6
            depmig(mi,mj) = 398. - 0.56 * oxymoy(mi,mj) -115. * log10(ztotchl) &
                          + 0.36 * hbl(mi,mj) -2.4 * tempmoy(mi,mj)
         ENDIF
      END DO   ;   END DO
      !$OMP END PARALLEL DO
      ! 
      ! Computation of the corresponding jk indice 
      ! -------------------------------------------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(depmig, gdepw, kmig)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF( depmig(mi,mj) >= gdepw(mi,mj,jk,Kmm) .AND. depmig(mi,mj) < gdepw(mi,mj,jk+1,Kmm) ) THEN
             kmig(mi,mj) = jk
          ENDIF
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      ! Correction of the migration depth and indice based on O2 levels
      ! If O2 is too low, imposing a migration depth at this low O2 levels
      ! would lead to negative O2 concentrations (respiration while O2 is close
      ! to 0. Thus, to avoid that problem, the migration depth is adjusted so
      ! that it falls above the OMZ
      ! ---------------------------------------------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, jkp1) &
      !$OMP SHARED(tr, kmig, depmig, gdept)
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF( t(mi,mj,N+1-kmig(mi,mj),Kbb,itemp+ntrc_salt+jpoxy) < 5E-6 ) THEN
            DO jk = kmig(mi,mj),1,-1
               jkp1 = jk+1
               IF( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpoxy) >= 5E-6 .AND. t(mi,mj,N+1-jkp1,Kbb,itemp+ntrc_salt+jpoxy) < 5E-6) THEN
                  kmig(mi,mj) = jk
                  depmig(mi,mj) = ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N)))
               ENDIF
            END DO
         ENDIF
      END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      IF( .false. )   CALL timing_stop('p4z_meso_depmig')
      !
   END SUBROUTINE p4z_meso_depmig

   INTEGER FUNCTION p4z_meso_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_meso_alloc  ***
      !!----------------------------------------------------------------------
      !
      ALLOCATE( depmig(Istrp:Iendp,Jstrp:Jendp), kmig(Istrp:Iendp,Jstrp:Jendp), STAT= p4z_meso_alloc  )
      !
      IF( p4z_meso_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p4z_meso_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p4z_meso_alloc


   !!======================================================================
END MODULE p4zmeso
