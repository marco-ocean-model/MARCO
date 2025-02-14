










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









MODULE p5zmicro
   !!======================================================================
   !!                         ***  MODULE p5zmicro  ***
   !! TOP :    Compute the sources/sinks for microzooplankton
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!             3.6  !  2015-05  (O. Aumont)  quota
   !!----------------------------------------------------------------------
   !!   p5z_micro       :   Compute the sources/sinks for microzooplankton
   !!   p5z_micro_init  :   Initialize and read the appropriate namelist
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !   Source Minus Sink variables
   USE iom             !  I/O manager
   USE prtctl          !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p5z_micro         ! called in p5zbio.F90
   PUBLIC   p5z_micro_init    ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  part        !: part of calcite not dissolved in microzoo guts
   REAL(wp), PUBLIC ::  xprefc     !: microzoo preference for POC 
   REAL(wp), PUBLIC ::  xprefn     !: microzoo preference for nanophyto
   REAL(wp), PUBLIC ::  xprefp     !: microzoo preference for picophyto
   REAL(wp), PUBLIC ::  xprefd     !: microzoo preference for diatoms
   REAL(wp), PUBLIC ::  xprefz     !: microzoo preference for microzoo
   REAL(wp), PUBLIC ::  xthreshdia  !: diatoms feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::  xthreshpic  !: picophyto feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::  xthreshphy  !: nanophyto threshold for microzooplankton 
   REAL(wp), PUBLIC ::  xthreshzoo  !: microzoo threshold for microzooplankton 
   REAL(wp), PUBLIC ::  xthreshpoc  !: poc threshold for microzooplankton 
   REAL(wp), PUBLIC ::  xthresh     !: feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::  resrat      !: exsudation rate of microzooplankton
   REAL(wp), PUBLIC ::  lmzrat      !: linear microzooplankton mortality rate 
   REAL(wp), PUBLIC ::  mzrat       !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::  grazrat     !: maximal microzoo grazing rate
   REAL(wp), PUBLIC ::  xkgraz      !: Half-saturation constant of assimilation
   REAL(wp), PUBLIC ::  unassc      !: Non-assimilated part of food
   REAL(wp), PUBLIC ::  unassn      !: Non-assimilated part of food
   REAL(wp), PUBLIC ::  unassp      !: Non-assimilated part of food
   REAL(wp), PUBLIC ::  epsher      !: Growth efficiency for microzoo
   REAL(wp), PUBLIC ::  epshermin   !: Minimum growth efficiency for microzoo
   REAL(wp), PUBLIC ::  srespir     !: half sturation constant for grazing 1 
   REAL(wp), PUBLIC ::  ssigma      !: Fraction excreted as semi-labile DOM
   REAL(wp), PUBLIC ::  xsigma      !: Width of the grazing window
   REAL(wp), PUBLIC ::  xsigmadel   !: Maximum additional width of the grazing window at low food density
   LOGICAL,  PUBLIC ::  bmetexc     !: Use of excess carbon for respiration

   REAL(wp)         :: rlogfactdn, rlogfactpn, rlogfactdp
   LOGICAL          :: l_dia_graz, l_dia_lprodz

   !! * Substitutions












































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p5zmicro.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p5z_micro( kt, knt, Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_micro  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for microzooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::  kt  ! ocean time step
      INTEGER, INTENT(in) ::  knt 
      INTEGER, INTENT(in) ::  Kbb, Krhs      ! time level indices
      !
      INTEGER  :: mi, mj, jk
      REAL(wp) :: zcompadi, zcompaz , zcompaph, zcompapoc, zcompapon, zcompapop
      REAL(wp) :: zcompapi, zgraze  , zdenom, zdenom2, zfact, zfood, zfoodlim
      REAL(wp) :: ztmp1, ztmp2, ztmp3, ztmp4, ztmp5, ztmptot
      REAL(wp) :: zepsherf, zepshert, zepsherq, zepsherv, zrespirc, zrespirn, zrespirp, zbasresb, zbasresi
      REAL(wp) :: zgraztotc, zgraztotn, zgraztotp, zgraztotf, zbasresn, zbasresp, zbasresf
      REAL(wp) :: zgradoc, zgradon, zgradop, zgraref, zgradoct, zgradont, zgradopt, zgrareft
      REAL(wp) :: zexcess, zgraren, zgrarep, zgrarem
      REAL(wp) :: zgrapoc, zgrapon, zgrapop, zgrapof, zprcaca, zmortz
      REAL(wp) :: zrespz, zltortz, ztortz, zgrasratf, zgrasratn, zgrasratp
      REAL(wp) :: zgraznc, zgraznn, zgraznp, zgrazpoc, zgrazpon, zgrazpop, zgrazpof
      REAL(wp) :: zgrazdc, zgrazdn, zgrazdp, zgrazdf, zgraznf, zgrazz
      REAL(wp) :: zgrazpc, zgrazpn, zgrazpp, zgrazpf, zbeta, zrfact2, zmetexcess
      REAL(wp) :: zsigma , zsigma2, zr_poc, zr_pic, zr_phy, zr_dia 
      REAL(wp) :: zsizepn, zsizedn, zsizedp, zdiffdn, zdiffpn, zdiffdp
      REAL(wp) :: zpexpod, zmaxsi, ztra
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N) :: zproportn, zproportd, zprodlig
      
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zgrazing, zfezoo, zzligprod, zw3d
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p5z_micro')
      !
      IF( kt == nittrc000 )  THEN
         l_dia_graz    = iom_use( "GRAZ1" ) .OR. iom_use( "FEZOO" ) .OR. iom_use( "MicroZo2" ) 
         l_dia_lprodz  = ln_ligand .AND. iom_use( "LPRODZ" )
         l_dia_graz = l_dia_graz .OR. l_diaadd
      ENDIF
      !
      IF( l_dia_graz ) THEN
         ALLOCATE( zgrazing(Istrp:Iendp,Jstrp:Jendp,N), zfezoo(Istrp:Iendp,Jstrp:Jendp,N) )
         zgrazing(Istrp:Iendp,Jstrp:Jendp,:) = 0.
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zfezoo(mi,mj,jk) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer)
         END DO   ;   END DO   ;   END DO
      ENDIF
      IF( l_dia_lprodz ) THEN
         ALLOCATE( zzligprod(Istrp:Iendp,Jstrp:Jendp,N) )
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zzligprod(mi,mj,jk) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jplgw)
         END DO   ;   END DO   ;   END DO
      ENDIF

      ! Use of excess carbon for metabolism
      zmetexcess = 0.0
      IF ( bmetexc ) zmetexcess = 1.0
      !
      ! Variables used to compute the preferences
      ! Proportion of nano and diatoms that are within the size range
      ! accessible to microzooplankton. 
      ! -------------------------------------------------------------
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF ( tmask(mi,mj,jk) == 1 ) THEN
            zproportd(mi,mj,jk) = EXP( 0.067 - 0.033 * sized(mi,mj,jk) * 6.0) / ( EXP( 0.067 -0.033 * 6.0 ) )
            zproportn(mi,mj,jk) = EXP( 0.131 - 0.047 * sizen(mi,mj,jk) * 4.0) / ( EXP( 0.131 -0.047 * 4.0 ) )
         ELSE
            zproportn(mi,mj,jk) = 1.0
            zproportd(mi,mj,jk) = 1.0
         ENDIF
      END DO   ;   END DO   ;   END DO
      !      
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zcompaz = MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo) - 1.e-9 ), 0.e0 )
         zfact     = xstep * tgfunc2(mi,mj,jk) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo)
         !  linear mortality of mesozooplankton
         !  A michaelis menten modulation term is used to avoid extinction of 
         !  microzooplankton at very low food concentrations. Mortality is 
         !  enhanced in low O2 waters
         !  -----------------------------------------------------------------

         !   Michaelis-Menten respiration rates of microzooplankton
         !   -----------------------------------------------------
         zrespz = resrat * zfact * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo) &
                 & / ( xkmort + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo) )

         !   Michaelis-Menten linear mortality rate of microzooplankton
         !   ----------------------------------------------------------
         zltortz = lmzrat * zfact * ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo) &
         &         / ( xkmort + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo) )  &
         &        + 3. * nitrfac(mi,mj,jk) )

         !  Zooplankton quadratic mortality. A square function has been selected with
         !  to mimic predation and disease (density dependent mortality). It also tends
         !  to stabilise the model
         !  -------------------------------------------------------------------------
         ztortz = mzrat * 1.e6 * zfact * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo) &
                 &  * (1. - nitrfac(mi,mj,jk))
         zmortz = ztortz + zltortz

         !   Computation of the abundance of the preys
         !   A threshold can be specified in the namelist
         !   Nanophyto and diatoms have a specific treatment with 
         !   teir preference decreasing with size.
         !   --------------------------------------------------------
         zcompadi  = zproportd(mi,mj,jk) * MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) - xthreshdia ), 0.e0 )
         zcompaph  = zproportn(mi,mj,jk) * MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) - xthreshphy ), 0.e0 )
         zcompaz   = MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo) - xthreshzoo ), 0.e0 )
         zcompapi  = MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppic) - xthreshpic ), 0.e0 )
         zcompapoc = MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) - xthreshpoc ), 0.e0 )
         zr_phy    = 1.0 / (t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) + 0.5*EPSILON(1.e0)) 
         zr_pic    = 1.0 / (t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppic) + 0.5*EPSILON(1.e0))
         zr_dia    = 1.0 / (t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) + 0.5*EPSILON(1.e0))
         zr_poc    = 1.0 / (t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) + 0.5*EPSILON(1.e0))
               
         ! Microzooplankton grazing
         ! The total amount of food is the sum of all preys accessible to mesozooplankton 
         ! multiplied by their food preference
         ! A threshold can be specified in the namelist (xthresh). However, when food 
         ! concentration is close to this threshold, it is decreased to avoid the 
         ! accumulation of food in the mesozoopelagic domain
         ! -------------------------------------------------------------------------------
         zfood     = xprefn * zcompaph + xprefc * zcompapoc + xprefd * zcompadi   &
         &           + xprefz * zcompaz + xprefp * zcompapi
         zfoodlim  = MAX( 0. , zfood - min(xthresh,0.5*zfood) )
         zdenom    = zfoodlim / ( xkgraz + zfoodlim )
         zgraze    = grazrat * zfact * (1. - nitrfac(mi,mj,jk)) 

         ! An active switching parameterization is used here.
         ! We don't use the KTW parameterization proposed by 
         ! Vallina et al. because it tends to produce too steady biomass
         ! composition and the variance of Chl is too low as it grazes
         ! too strongly on winning organisms. We use a generalized
         ! switching parameterization proposed by Morozov and 
         ! Petrovskii (2013)
         ! ------------------------------------------------------------  
         ! The width of the selection window is increased when preys
         ! have low abundance, .i.e. zooplankton become less specific 
         ! to avoid starvation.
         ! ----------------------------------------------------------
         zdenom2   = zdenom * zdenom
         zsigma    = 1.0 - zdenom2/( 0.05 * 0.05 + zdenom2 )
         zsigma    = xsigma + xsigmadel * zsigma
         zsigma2   = 2.0 * zsigma * zsigma
         !
         zsizepn   = rlogfactpn + ( logsizep(mi,mj,jk) - logsizen(mi,mj,jk) )
         zsizedn   = rlogfactdn + ( logsizen(mi,mj,jk) - logsized(mi,mj,jk) )
         zsizedp   = rlogfactdp + ( logsizep(mi,mj,jk) - logsized(mi,mj,jk) )
         zdiffpn   = EXP( -zsizepn * zsizepn / zsigma2 )
         zdiffdn   = EXP( -zsizedn * zsizedn / zsigma2 )
         zdiffdp   = EXP( -zsizedp * zsizedp / zsigma2 )
         !
         ztmp1     = xprefn * zcompaph * ( zcompaph + zdiffdn * zcompadi + zdiffpn * zcompapi )
         ztmp2     = xprefp * zcompapi * ( zcompapi + zdiffpn * zcompaph + zdiffdp * zcompadi )
         ztmp3     = xprefc * zcompapoc * zcompapoc
         ztmp4     = xprefd * zcompadi * ( zdiffdp * zcompapi + zdiffdn * zcompaph + zcompadi )
         ztmp5     = xprefz * zcompaz * zcompaz
         ztmptot   = ztmp1 + ztmp2 + ztmp3 + ztmp4 + ztmp5 + 0.5*EPSILON(1.e0)
         ztmp1     = ztmp1 / ztmptot
         ztmp2     = ztmp2 / ztmptot
         ztmp3     = ztmp3 / ztmptot
         ztmp4     = ztmp4 / ztmptot
         ztmp5     = ztmp5 / ztmptot

         !   Microzooplankton regular grazing on the different preys
         !   -------------------------------------------------------
               !   Nanophytoplankton
         zgraznc   = zgraze  * ztmp1  * zdenom
         zgraznn   = zgraznc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnph) * zr_phy
         zgraznp   = zgraznc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppph) * zr_phy
         zgraznf   = zgraznc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnfe) * zr_phy

               ! Picophytoplankton
         zgrazpc   = zgraze  * ztmp2  * zdenom
         zgrazpn   = zgrazpc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnpi) * zr_pic
         zgrazpp   = zgrazpc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpppi) * zr_pic
         zgrazpf   = zgrazpc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppfe) * zr_pic

               ! Microzooplankton
         zgrazz    = zgraze  * ztmp5   * zdenom

               ! small POC
         zgrazpoc  = zgraze   * ztmp3   * zdenom
         zgrazpon  = zgrazpoc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppon) * zr_poc
         zgrazpop  = zgrazpoc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppop) * zr_poc
         zgrazpof  = zgrazpoc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsfe) * zr_poc

               ! Diatoms
         zgrazdc   = zgraze  * ztmp4  * zdenom
         zgrazdn   = zgrazdc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpndi) * zr_dia
         zgrazdp   = zgrazdc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppdi) * zr_dia
         zgrazdf   = zgrazdc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdfe) * zr_dia
         !
               ! Total ingestion rates in C, P, Fe, N
         zgraztotc = zgraznc + zgrazpoc + zgrazdc + zgrazz + zgrazpc ! Grazing by microzooplankton
         IF( l_dia_graz ) zgrazing(mi,mj,jk) = zgraztotc

         zgraztotn = zgraznn + zgrazpn + zgrazpon + zgrazdn + zgrazz * no3rat3
         zgraztotp = zgraznp + zgrazpp + zgrazpop + zgrazdp + zgrazz * po4rat3
         zgraztotf = zgraznf + zgrazpf + zgrazpof + zgrazdf + zgrazz * feratz
         !
         !   Stoichiometruc ratios of the food ingested by zooplanton 
         !   --------------------------------------------------------
         zgrasratf = (zgraztotf + 0.5*EPSILON(1.e0)) / ( zgraztotc + 0.5*EPSILON(1.e0) )
         zgrasratn = (zgraztotn + 0.5*EPSILON(1.e0)) / ( zgraztotc + 0.5*EPSILON(1.e0) )
         zgrasratp = (zgraztotp + 0.5*EPSILON(1.e0)) / ( zgraztotc + 0.5*EPSILON(1.e0) )

         ! Mesozooplankton efficiency. 
         ! We adopt a formulation proposed by Mitra et al. (2007)
         ! The gross growth efficiency is controled by the most limiting nutrient.
         ! Growth is also further decreased when the food quality is poor. This is currently
         ! hard coded : it can be decreased by up to 50% (zepsherq)
         ! GGE can also be decreased when food quantity is high, zepsherf (Montagnes and 
         ! Fulton, 2012)
         ! -----------------------------------------------------------------------------------
         zepshert  = MIN( 1., zgrasratn/ no3rat3, zgrasratp/ po4rat3, zgrasratf / feratz)
         zbeta     = MAX( 0., (epsher - epshermin) )
         ! Food density deprivation of GGE
         zepsherf  = epshermin + zbeta / ( 1.0 + 0.04E6 * 12. * zfood * zbeta )
         ! Food quality deprivation of GGE
         zepsherq  = 0.5 + (1.0 - 0.5) * zepshert * ( 1.0 + 1.0 ) / ( zepshert + 1.0 )
         ! Actual GGE
         zepsherv  = zepsherf * zepshert * zepsherq

         ! Respiration of microzooplankton
         ! Excess carbon in the food is used preferentially
         ! when activated by zmetexcess
         ! ------------------------------------------------
         zexcess   = zgraztotc * zepsherq * zepsherf * (1.0 - zepshert) * zmetexcess
         zbasresb  = MAX(0., zrespz - zexcess)
         zbasresi  = zexcess + MIN(0., zrespz - zexcess)  
         zrespirc  = srespir * zepsherv * zgraztotc + zbasresb
         
         ! When excess carbon is used, the other elements in excess
         ! are also used proportionally to their abundance
         ! --------------------------------------------------------
         zexcess   = ( zgrasratn/ no3rat3 - zepshert ) / ( 1.0 - zepshert + 0.5*EPSILON(1.e0))
         zbasresn  = zbasresi * zexcess * zgrasratn 
         zexcess   = ( zgrasratp/ po4rat3 - zepshert ) / ( 1.0 - zepshert + 0.5*EPSILON(1.e0))
         zbasresp  = zbasresi * zexcess * zgrasratp
         zexcess   = ( zgrasratf/ feratz  - zepshert ) / ( 1.0 - zepshert + 0.5*EPSILON(1.e0))
         zbasresf  = zbasresi * zexcess * zgrasratf

         ! Voiding of the excessive elements as DOM
         ! ----------------------------------------
         zgradoct  = (1. - unassc - zepsherv) * zgraztotc - zbasresi  
         zgradont  = (1. - unassn) * zgraztotn - zepsherv * no3rat3 * zgraztotc - zbasresn
         zgradopt  = (1. - unassp) * zgraztotp - zepsherv * po4rat3 * zgraztotc - zbasresp
         zgrareft  = (1. - unassc) * zgraztotf - zepsherv * feratz  * zgraztotc - zbasresf

         ! Since only semilabile DOM is represented in 
         ! part of DOM is in fact labile and is then released
         ! as dissolved inorganic compounds (ssigma)
         ! --------------------------------------------------
         zgradoc   =  (1.0 - ssigma) * zgradoct
         zgradon   =  (1.0 - ssigma) * zgradont
         zgradop   =  (1.0 - ssigma) * zgradopt
         zgrarem   = ssigma * zgradoct
         zgraren   = ssigma * zgradont
         zgrarep   = ssigma * zgradopt
         zgraref   = zgrareft

         ! Defecation as a result of non assimilated products
         ! --------------------------------------------------
         zgrapoc   = zgraztotc * unassc
         zgrapon   = zgraztotn * unassn
         zgrapop   = zgraztotp * unassp
         zgrapof   = zgraztotf * unassc

         ! Addition of respiration to the release of inorganic nutrients
         ! -------------------------------------------------------------
         zgrarem   = zgrarem + zbasresi + zrespirc
         zgraren   = zgraren + zbasresn + zrespirc * no3rat3
         zgrarep   = zgrarep + zbasresp + zrespirc * po4rat3
         zgraref   = zgraref + zbasresf + zrespirc * feratz


         !   Update of the TRA arrays
         !   ------------------------
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) + zgrarep
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnh4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnh4) + zgraren
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) + zgradoc
         !
         IF( ln_ligand ) THEN 
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jplgw) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jplgw) + zgradoc * ldocz
         ENDIF
         !
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdon) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdon) + zgradon
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdop) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdop) + zgradop
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) - o2ut * zgrarem 
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) + zgraref
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpzoo) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpzoo) &
                 &                + zepsherv * zgraztotc - zrespirc - zmortz - zgrazz
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpphy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpphy) - zgraznc
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnph) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnph) - zgraznn
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppph) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppph) - zgraznp
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppic) - zgrazpc
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnpi) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnpi) - zgrazpn
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpppi) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpppi) - zgrazpp
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdia) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdia) - zgrazdc
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpndi) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpndi) - zgrazdn
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppdi) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppdi) - zgrazdp
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnch) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnch) &
                 &                - zgraznc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnch) * zr_phy
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppch) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppch) &
                 &                - zgrazpc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppch) * zr_pic
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdch) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdch) &
                 &                - zgrazdc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdch) * zr_dia
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdsi) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdsi) &
                 &                - zgrazdc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdsi) * zr_dia
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgsi) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgsi) &
                 &                + zgrazdc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdsi) * zr_dia
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnfe) - zgraznf
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppfe) - zgrazpf
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdfe) - zgrazdf
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) + zmortz + zgrapoc - zgrazpoc 
         prodpoc(mi,mj,jk) = prodpoc(mi,mj,jk) + zmortz + zgrapoc
         conspoc(mi,mj,jk) = conspoc(mi,mj,jk) - zgrazpoc
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppon) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppon) &
                 &                + no3rat3 * zmortz + zgrapon - zgrazpon
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppop) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppop) &
                 &                + po4rat3 * zmortz + zgrapop - zgrazpop
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) &
                 &                + feratz * zmortz  + zgrapof - zgrazpof
         !
         ! calcite production
         zprcaca = xfracal(mi,mj,jk) * zgraznc
         prodcal(mi,mj,jk) = prodcal(mi,mj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
         !
         zprcaca = part * zprcaca
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) + zgrarem - zprcaca
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) - 2. * zprcaca + rno3 * zgraren
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpcal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpcal) + zprcaca
      END DO   ;   END DO   ;   END DO
      !
      IF( .false. .AND. knt == nrdttrc ) THEN
        !
        IF( l_dia_graz ) THEN  !   Total grazing of phyto by zooplankton
            ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
            DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zw3d(mi,mj,N+1-jk) =  zgrazing(mi,mj,jk) * 1.e+3 * rfact2r * tmask(mi,mj,jk)
            END DO   ;   END DO   ;   END DO
            CALL iom_put( "GRAZ1" , zw3d )  ! conversion in mol/m2/s
            CALL iom_put( "MicroZo2" , zw3d * ( 1. - epsher - unassc ) * (-o2ut) * ssigma ) ! o2 consumption by Microzoo
            !
            DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               ztra = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer)
               zw3d(mi,mj,N+1-jk) = ( ztra - zfezoo(mi,mj,jk) ) &
                  &             * 1e9 * 1.e+3 * rfact2r * tmask(mi,mj,jk) ! conversion in nmol/m2/s
            END DO   ;   END DO   ;   END DO
           CALL iom_put( "FEZOO", zw3d )
            DEALLOCATE( zw3d, zfezoo )
        ENDIF
        !
        IF( l_dia_lprodz ) THEN
            ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
            DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               ztra = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jplgw)
               zw3d(mi,mj,N+1-jk) = ( ztra - zzligprod(mi,mj,jk) ) &
                   &             * 1e9 * 1.e+3 * rfact2r * tmask(mi,mj,jk) ! conversion in nmol/m2/s
            END DO   ;   END DO   ;   END DO
           CALL iom_put( "LPRODZ", zw3d )
           DEALLOCATE( zzligprod, zw3d )
        ENDIF
        !
      ENDIF
      !
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         bioFlux(mi,mj,N+1-jk,Ngrapoc) = zgrazing(mi,mj,jk) * 1.e+3 * rfact2r * tmask(mi,mj,jk) !  grazing of phyto by microzoo
         bioFlux(mi,mj,N+1-jk,Nmico2) = zgrazing(mi,mj,jk) * ( 1.- epsher - unassc ) &
          &                  * (-o2ut) * ssigma * 1.e+3 * rfact2r * tmask(mi,mj,jk)   ! o2 consumption by Microzoo
      END DO   ;   END DO   ;   END DO
      IF( l_dia_graz ) DEALLOCATE( zgrazing )
      !
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('micro')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
  !       CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( .false. )   CALL timing_stop('p5z_micro')
      !
   END SUBROUTINE p5z_micro


   SUBROUTINE p5z_micro_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p5z_micro_init  ***
      !!
      !! ** Purpose :   Initialization of microzooplankton parameters
      !!
      !! ** Method  :   Read the namp5zzoo namelist and check the parameters
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namp5zzoo
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !!
      NAMELIST/namp5zzoo/ part, grazrat, bmetexc, resrat, lmzrat, mzrat, xprefc, xprefn, &
         &                xprefp, xprefd, xprefz, xthreshdia, xthreshphy, &
         &                xthreshpic, xthreshpoc, xthreshzoo, xthresh, xkgraz, &
         &                epsher, epshermin, ssigma, srespir, unassc, unassn, unassp,   &
         &                xsigma, xsigmadel   
      !!----------------------------------------------------------------------
      !
      REWIND(numnatp_ref);READ(numnatp_ref,namp5zzoo,IOSTAT=ios);CALL ctl_nam(ios,"namp5zzoo (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,namp5zzoo,IOSTAT=ios);CALL ctl_nam(ios,"namp5zzoo (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE ( numonp, namp5zzoo )
      !
      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*) ' '
         WRITE(stdout,*) ' Namelist parameters for microzooplankton, nampiszooq'
         WRITE(stdout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(stdout,*) '    part of calcite not dissolved in microzoo guts  part        =', part
         WRITE(stdout,*) '    microzoo preference for POC                     xprefc     =', xprefc
         WRITE(stdout,*) '    microzoo preference for nano                    xprefn     =', xprefn
         WRITE(stdout,*) '    microzoo preference for pico                    xprefp     =', xprefp
         WRITE(stdout,*) '    microzoo preference for diatoms                 xprefd     =', xprefd
         WRITE(stdout,*) '    microzoo preference for microzoo                xprefz     =', xprefz
         WRITE(stdout,*) '    diatoms feeding threshold  for microzoo         xthreshdia  =', xthreshdia
         WRITE(stdout,*) '    nanophyto feeding threshold for microzoo        xthreshphy  =', xthreshphy
         WRITE(stdout,*) '    picophyto feeding threshold for microzoo        xthreshpic  =', xthreshpic
         WRITE(stdout,*) '    poc feeding threshold for microzoo              xthreshpoc  =', xthreshpoc
         WRITE(stdout,*) '    microzoo feeding threshold for microzoo         xthreshzoo  =', xthreshzoo
         WRITE(stdout,*) '    feeding threshold for microzooplankton          xthresh     =', xthresh
         WRITE(stdout,*) '    exsudation rate of microzooplankton             resrat      =', resrat
         WRITE(stdout,*) '    linear microzooplankton mortality rate          lmzrat      =', lmzrat
         WRITE(stdout,*) '    microzooplankton mortality rate                 mzrat       =', mzrat
         WRITE(stdout,*) '    maximal microzoo grazing rate                   grazrat     =', grazrat
         WRITE(stdout,*) '    C egested fraction of fodd by microzoo          unassc      =', unassc
         WRITE(stdout,*) '    N egested fraction of fodd by microzoo          unassn      =', unassn
         WRITE(stdout,*) '    P egested fraction of fodd by microzoo          unassp      =', unassp
         WRITE(stdout,*) '    Efficicency of microzoo growth                  epsher      =', epsher
         WRITE(stdout,*) '    Minimum Efficiency of Microzoo growth           epshermin   =', epshermin
         WRITE(stdout,*) '    Fraction excreted as semi-labile DOM            ssigma      =', ssigma
         WRITE(stdout,*) '    Active respiration                              srespir     =', srespir
         WRITE(stdout,*) '    half sturation constant for grazing 1           xkgraz      =', xkgraz
         WRITE(stdout,*) '    Use of excess carbon for respiration            bmetexc     =', bmetexc
         WRITE(stdout,*) '      Width of the grazing window                     xsigma      =', xsigma
         WRITE(stdout,*) '      Maximum additional width of the grazing window  xsigmadel   =', xsigmadel
      ENDIF
      !
      rlogfactpn = LOG(0.8 / 4.0)
      rlogfactdn = LOG(4.0 / 6.0)
      rlogfactdp = LOG(0.8 / 6.0)
      !
   END SUBROUTINE p5z_micro_init


   !!======================================================================
END MODULE p5zmicro
