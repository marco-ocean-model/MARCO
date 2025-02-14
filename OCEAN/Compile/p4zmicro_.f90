










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









MODULE p4zmicro
   !!======================================================================
   !!                         ***  MODULE p4zmicro  ***
   !! TOP :    Compute the sources/sinks for microzooplankton
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!             3.*  !  2025-02  (S. Maishal, R. Person)  and optimization
   !!----------------------------------------------------------------------
   !!   p4z_micro      : Compute the sources/sinks for microzooplankton
   !!   p4z_micro_init : Initialize and read the appropriate namelist
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE p4zprod         ! production
   USE iom             ! I/O manager
   USE prtctl          ! print control for debugging

   IMPLICIT NONE
   PRIVATE

   !! * Shared module variables
   PUBLIC   p4z_micro         ! called in p4zbio.F90
   PUBLIC   p4z_micro_init    ! called in trcsms_pisces.F90

   REAL(wp), PUBLIC ::   part        !: part of calcite not dissolved in microzoo guts
   REAL(wp), PUBLIC ::   xprefc      !: microzoo preference for POC 
   REAL(wp), PUBLIC ::   xprefn      !: microzoo preference for nanophyto
   REAL(wp), PUBLIC ::   xprefd      !: microzoo preference for diatoms
   REAL(wp), PUBLIC ::   xprefz      !: microzoo preference for microzooplankton
   REAL(wp), PUBLIC ::   xthreshdia  !: diatoms feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::   xthreshphy  !: nanophyto threshold for microzooplankton 
   REAL(wp), PUBLIC ::   xthreshpoc  !: poc threshold for microzooplankton 
   REAL(wp), PUBLIC ::   xthreshzoo  !: microzoo threshold for microzooplankton 
   REAL(wp), PUBLIC ::   xthresh     !: feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::   resrat      !: exsudation rate of microzooplankton
   REAL(wp), PUBLIC ::   mzrat       !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::   grazrat     !: maximal microzoo grazing rate
   REAL(wp), PUBLIC ::   xkgraz      !: Half-saturation constant of assimilation
   REAL(wp), PUBLIC ::   unass       !: Non-assimilated part of food
   REAL(wp), PUBLIC ::   sigma1      !: Fraction of microzoo excretion as DOM 
   REAL(wp), PUBLIC ::   epsher      !: growth efficiency for grazing 1 
   REAL(wp), PUBLIC ::   epshermin   !: minimum growth efficiency for grazing 1
   REAL(wp), PUBLIC ::   xsigma      !: Width of the grazing window
   REAL(wp), PUBLIC ::   xsigmadel   !: Maximum additional width of the grazing window at low food density 

   LOGICAL          :: l_dia_graz, l_dia_lprodz

   !! * Substitutions












































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zmicro.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_micro( kt, knt, Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_micro  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for microzooplankton
      !!                This includes ingestion and assimilation, flux feeding
      !!                and mortality. We use a passive prey switching  
      !!                parameterization.
      !!                All living compartments smaller than microzooplankton
      !!                are potential preys of microzooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time step
      INTEGER, INTENT(in) ::   knt   ! ??? 
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      !
      INTEGER  :: mi, mj, jk
      REAL(wp) :: zcompadi, zcompaz , zcompaph, zcompapoc, zratio
      REAL(wp) :: zgraze, zdenom, zdenom2, zfact, zfood, zfoodlim, zbeta
      REAL(wp) :: zepsherf, zepshert, zepsherq, zepsherv, zgrarsig, zgraztotc, zgraztotn, zgraztotf
      REAL(wp) :: zgrarem, zgrafer, zgrapoc, zprcaca, zmortz
      REAL(wp) :: zrespz, ztortz, zgrasratf, zgrasratn
      REAL(wp) :: zgraznc, zgrazz, zgrazpoc, zgrazdc, zgrazpof, zgrazdf, zgraznf
      REAL(wp) :: zsigma, zsigma2, zsizedn, zdiffdn, ztmp1, ztmp2, ztmp3, ztmp4, ztmptot, zproport
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zgrazing, zfezoo, zzligprod, zw3d
      CHARACTER (len=25) :: charout

      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p4z_micro')
      !
      IF( kt == nittrc000 )  THEN 
         l_dia_graz    = iom_use( "GRAZ1" ) .OR. iom_use( "FEZOO" ) .OR. iom_use( "MicroZo2" ) 
         l_dia_lprodz  = ln_ligand .AND. iom_use( "LPRODZ" ) 
         l_dia_graz = l_dia_graz .OR. l_diaadd
      ENDIF
      !
      IF( l_dia_graz ) THEN
         ALLOCATE( zgrazing(Istrp:Iendp,Jstrp:Jendp,N), zfezoo(Istrp:Iendp,Jstrp:Jendp,N) )
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk) &
         !$OMP SHARED(zfezoo, tr)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zfezoo(mi,mj,jk) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer)
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
      ENDIF
      !
      IF( l_dia_lprodz ) THEN
         ALLOCATE( zzligprod(Istrp:Iendp,Jstrp:Jendp,N) )
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk) &
         !$OMP SHARED(zzligprod, tr)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zzligprod(mi,mj,jk) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jplgw)
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
      ENDIF
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(tr, xstep, tgfunc2, sized, nitrfac, 0.5*EPSILON(1.e0), xprefn, &
                   xprefc, xprefd, xprefz, xthresh, sigma1, grazrat, &
                   epsher, epshermin, feratz, unass, xfracal, part)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zcompaz = MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo) - 1.e-9 ), 0.e0 )
         zfact   = xstep * tgfunc2(mi,mj,jk) * zcompaz

         ! Proportion of diatoms that are within the size range
         ! accessible to microzooplankton. 
         zproport  = sized(mi,mj,jk)**(-0.48) * ( 1.0 - ( sized(mi,mj,jk)**1.6 - 1.0 ) / 14.0 )

         !  linear mortality of mesozooplankton
         !  A michaelis menten modulation term is used to avoid extinction of 
         !  microzooplankton at very low food concentrations. Mortality is 
         !  enhanced in low O2 waters
         !  -----------------------------------------------------------------
         zrespz = resrat * zfact * ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo) &
            &    / ( xkmort + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo) )  &
            &   + 3. * nitrfac(mi,mj,jk) )

         !  Zooplankton quadratic mortality. A square function has been selected with
         !  to mimic predation and disease (density dependent mortality). It also tends
         !  to stabilise the model
         !  -------------------------------------------------------------------------
         ztortz = mzrat * 1.e6 * zfact * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo) &
             &        * (1. - nitrfac(mi,mj,jk))
         zmortz = ztortz + zrespz

         !   Computation of the abundance of the preys
         !   A threshold can be specified in the namelist
         !   Diatoms have a specific treatment. WHen concentrations 
         !   exceed a certain value, diatoms are suppposed to be too 
         !   big for microzooplankton.
         !   --------------------------------------------------------
         zcompadi  = zproport * MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) - xthreshdia ), 0.e0 )
         zcompaph  = MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) - xthreshphy ), 0.e0 )
         zcompapoc = MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) - xthreshpoc ), 0.e0 )
         zcompaz   = MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo) - xthreshzoo ), 0.e0 )
 
         ! Microzooplankton grazing
         ! The total amount of food is the sum of all preys accessible to mesozooplankton 
         ! multiplied by their food preference
         ! A threshold can be specified in the namelist (xthresh). However, when food 
         ! concentration is close to this threshold, it is decreased to avoid the 
         ! accumulation of food in the mesozoopelagic domain
         ! -------------------------------------------------------------------------------
         zfood     = xprefn * zcompaph + xprefc * zcompapoc + xprefd * zcompadi + xprefz * zcompaz
         zfoodlim  = MAX( 0. , zfood - MIN(xthresh,0.5*zfood) )
         zdenom    = zfoodlim / ( xkgraz + zfoodlim )
         zgraze    = grazrat * xstep * tgfunc2(mi,mj,jk) &
             &      * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo) * (1. - nitrfac(mi,mj,jk))

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
         zdenom2 = zdenom * zdenom
         zsigma = 1.0 - zdenom2/(0.05 * 0.05 + zdenom2)
         zsigma = xsigma + xsigmadel * zsigma
         zsigma2 = zsigma * zsigma
         !
         zsizedn = ABS(LOG(1.67 * sizen(mi,mj,jk) / (5.0 * sized(mi,mj,jk) + 0.5*EPSILON(1.e0) )) )
         zdiffdn = EXP( -zsizedn * zsizedn / zsigma2 )
         ztmp1 = xprefn * zcompaph * ( zcompaph + zdiffdn * zcompadi ) 
         ztmp2 = xprefd * zcompadi * ( zdiffdn * zcompaph + zcompadi )
         ztmp3 = xprefc * zcompapoc * zcompapoc 
         ztmp4 = xprefz * zcompaz * zcompaz
         ztmptot = ztmp1 + ztmp2 + ztmp3 + ztmp4 + 0.5*EPSILON(1.e0)
         ztmp1 = ztmp1 / ztmptot
         ztmp2 = ztmp2 / ztmptot
         ztmp3 = ztmp3 / ztmptot
         ztmp4 = ztmp4 / ztmptot

         ! Ingestion terms on the different preys of microzooplankton
         zgraznc   = zgraze   * ztmp1 * zdenom  ! Nanophytoplankton
         zgrazdc   = zgraze   * ztmp2 * zdenom  ! Diatoms
         zgrazpoc  = zgraze   * ztmp3 * zdenom  ! POC
         zgrazz    = zgraze   * ztmp4 * zdenom  ! Microzoo

         ! Ingestion terms on the iron content of the different preys
         zgraznf   = zgraznc  * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnfe) &
            &        / (t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) + 0.5*EPSILON(1.e0))
         zgrazpof  = zgrazpoc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsfe) &
            &        / (t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) + 0.5*EPSILON(1.e0))
         zgrazdf   = zgrazdc  * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdfe) &
            &        / (t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) + 0.5*EPSILON(1.e0))
         !
         ! Total ingestion rate in C, Fe, N units
         zgraztotc = zgraznc + zgrazpoc + zgrazdc + zgrazz
         IF( l_dia_graz )   zgrazing(mi,mj,jk) = zgraztotc

         zgraztotf = zgraznf + zgrazdf  + zgrazpof + zgrazz * feratz
         zgraztotn = zgraznc * quotan(mi,mj,jk) + zgrazpoc &
             &      + zgrazdc * quotad(mi,mj,jk) + zgrazz

         !   Stoichiometruc ratios of the food ingested by zooplankton 
         !   --------------------------------------------------------
         zgrasratf = ( zgraztotf + 0.5*EPSILON(1.e0) ) / ( zgraztotc + 0.5*EPSILON(1.e0) )
         zgrasratn = ( zgraztotn + 0.5*EPSILON(1.e0) ) / ( zgraztotc + 0.5*EPSILON(1.e0) )

         ! Microzooplankton efficiency. 
         ! We adopt a formulation proposed by Mitra et al. (2007)
         ! The gross growth efficiency is controled by the most limiting nutrient.
         ! Growth is also further decreased when the food quality is poor. This is currently
         ! hard coded : it can be decreased by up to 50% (zepsherq)
         ! GGE can also be decreased when food quantity is high, zepsherf (Montagnes and 
         ! Fulton, 2012)
         ! -----------------------------------------------------------------------------
         zepshert  =  MIN( 1., zgrasratn, zgrasratf / feratz)
         zbeta     =  MAX(0., (epsher - epshermin) )
         ! Food quantity deprivation of the GGE
         zepsherf  = epshermin + zbeta / ( 1.0 + 0.04E6 * 12. * zfood * zbeta )
         ! Food quality deprivation of the GGE
         zepsherq  = 0.5 + (1.0 - 0.5) * zepshert * ( 1.0 + 1.0 ) / ( zepshert + 1.0 )
         ! Actual GGE of microzooplankton
         zepsherv  = zepsherf * zepshert * zepsherq
         ! Excretion of Fe
         zgrafer   = zgraztotc * MAX( 0. , ( 1. - unass ) * zgrasratf - feratz * zepsherv ) 
         ! Excretion of C, N, P
         zgrarem   = zgraztotc * ( 1. - zepsherv - unass )
         ! Egestion of C, N, P
         zgrapoc   = zgraztotc * unass

         !  Update of the TRA arrays
         !  ------------------------
         ! Fraction of excretion as inorganic nutrients and DIC
         zgrarsig  = zgrarem * sigma1
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) + zgrarsig
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnh4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnh4) + zgrarsig
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) + zgrarem - zgrarsig
         !
         IF( ln_ligand ) THEN
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jplgw) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jplgw) + (zgrarem - zgrarsig) * ldocz
         ENDIF
         !
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) - o2ut * zgrarsig
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) + zgrafer
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) + zgrapoc
         prodpoc(mi,mj,jk)   = prodpoc(mi,mj,jk) + zgrapoc
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) + zgraztotf * unass
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) + zgrarsig
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) + rno3 * zgrarsig
         !   Update the arrays TRA which contain the biological sources and sinks
         !   --------------------------------------------------------------------
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpzoo) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpzoo) - zmortz &
                 &                + zepsherv * zgraztotc - zgrazz 
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpphy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpphy) - zgraznc
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdia) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdia) - zgrazdc
         !
         zratio = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnch)/(t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy)+0.5*EPSILON(1.e0))
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnch) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnch) - zgraznc * zratio
         zratio = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdch)/(t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia)+0.5*EPSILON(1.e0))
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdch) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdch) - zgrazdc * zratio
         zratio = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdsi)/(t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia)+0.5*EPSILON(1.e0))
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdsi) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdsi) - zgrazdc * zratio
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgsi) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgsi) + zgrazdc * zratio
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnfe) - zgraznf
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdfe) - zgrazdf
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) + zmortz - zgrazpoc
         prodpoc(mi,mj,jk) = prodpoc(mi,mj,jk) + zmortz
         conspoc(mi,mj,jk) = conspoc(mi,mj,jk) - zgrazpoc
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) + feratz * zmortz - zgrazpof
         !
         ! Calcite remineralization due to zooplankton activity
         ! part of the ingested calcite is not dissolving in the acidic gut
         ! ----------------------------------------------------------------
         zprcaca = xfracal(mi,mj,jk) * zgraznc
         prodcal(mi,mj,jk) = prodcal(mi,mj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)

         !
         zprcaca = part * zprcaca
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) - zprcaca
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) - 2. * zprcaca
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpcal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpcal) + zprcaca
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      IF( .false. .AND. knt == nrdttrc ) THEN
        !
        IF( l_dia_graz ) THEN  !   Total grazing of phyto by zooplankton
            ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) ) ; zw3d(:,:,:) = 0._wp
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(mi, mj, N+1-jk) &
            !$OMP SHARED(zw3d, zgrazing, rfact2r, tmask)
            DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zw3d(mi,mj,N+1-jk) = zgrazing(mi,mj,jk) * 1.e+3 * rfact2r * tmask(mi,mj,jk) ! conversion in nmol/m2/s
            END DO   ;   END DO   ;   END DO
            !$OMP END PARALLEL DO
            CALL iom_put( "GRAZ1" , zw3d )  ! conversion in mol/m2/s
            CALL iom_put( "MicroZo2" , zw3d * ( 1. - epsher - unass ) * (-o2ut) * sigma1 ) ! o2 consumption by Microzoo
            !
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(mi, mj, N+1-jk) &
            !$OMP SHARED(zw3d, tr, jpfer, Krhs, zfezoo, rfact2r, tmask)
            DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zw3d(mi,mj,N+1-jk) = ( t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) - zfezoo(mi,mj,jk) ) &
                  &              * 1e9 * 1.e+3 * rfact2r * tmask(mi,mj,jk) ! conversion in nmol/m2/s
            END DO   ;   END DO   ;   END DO
            !$OMP END PARALLEL DO
           CALL iom_put( "FEZOO", zw3d )
            DEALLOCATE( zw3d, zfezoo )
        ENDIF
        !
        IF( l_dia_lprodz ) THEN
            ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(mi, mj, N+1-jk) &
            !$OMP SHARED(zw3d, tr, jplgw, Krhs, zzligprod, rfact2r, tmask)
            DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zw3d(mi,mj,N+1-jk) = ( t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jplgw) - zzligprod(mi,mj,jk) ) &
                   &                * 1e9 * 1.e+3 * rfact2r * tmask(mi,mj,jk) ! conversion in nmol/m2/s
            END DO   ;   END DO   ;   END DO
            !$OMP END PARALLEL DO
           CALL iom_put( "LPRODZ", zw3d )
           DEALLOCATE( zzligprod, zw3d )
        ENDIF
        !
      ENDIF
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, N+1-jk) &
      !$OMP SHARED(trc3d, zgrazing, rfact2r, tmask, epsher, unass, o2ut, sigma1)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         bioFlux(mi,mj,N+1-jk,Ngrapoc) = zgrazing(mi,mj,jk) * 1.e+3 * rfact2r &
                                    * tmask(mi,mj,jk) !  grazing of phyto by microzoo
         bioFlux(mi,mj,N+1-jk,Nmico2)  = zgrazing(mi,mj,jk) * ( 1.- epsher - unass ) &
          &                         * (-o2ut) * sigma1 * 1.e+3 * rfact2r &
                                    * tmask(mi,mj,jk) ! o2 consumption by Microzoo
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      IF( l_dia_graz ) DEALLOCATE( zgrazing )
      !
      IF(sn_cfctl%l_prttrc) THEN ! print mean trends (used for debugging)
         WRITE(charout, FMT="('micro')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
    !     CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( .false. )   CALL timing_stop('p4z_micro')
      !
   END SUBROUTINE p4z_micro


   SUBROUTINE p4z_micro_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_micro_init  ***
      !!
      !! ** Purpose :   Initialization of microzooplankton parameters
      !!
      !! ** Method  :   Read the namp4zzoo namelist and check the parameters
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namp4zzoo
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios ! Local integer
      !
      NAMELIST/namp4zzoo/ part, grazrat, resrat, mzrat, xprefn, xprefc, &
         &                xprefd, xprefz, xthreshdia, xthreshphy, xthreshpoc, &
         &                xthreshzoo, xthresh, xkgraz, epsher, epshermin, &
         &                sigma1, unass, xsigma, xsigmadel
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*) 
         WRITE(stdout,*) 'p4z_micro_init : Initialization of microzooplankton parameters'
         WRITE(stdout,*) '~~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND(numnatp_ref);READ(numnatp_ref,namp4zzoo,IOSTAT=ios);CALL ctl_nam(ios,"namp4zzoo (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,namp4zzoo,IOSTAT=ios);CALL ctl_nam(ios,"namp4zzoo (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE( numonp, namp4zzoo )
      !
      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*) '   Namelist : namp4zzoo'
         WRITE(stdout,*) '      part of calcite not dissolved in microzoo guts  part        =', part
         WRITE(stdout,*) '      microzoo preference for POC                     xprefc      =', xprefc
         WRITE(stdout,*) '      microzoo preference for nano                    xprefn      =', xprefn
         WRITE(stdout,*) '      microzoo preference for diatoms                 xprefd      =', xprefd
         WRITE(stdout,*) '      microzoo preference for microzooplankton        xprefz      =', xprefz
         WRITE(stdout,*) '      diatoms feeding threshold  for microzoo         xthreshdia  =', xthreshdia
         WRITE(stdout,*) '      nanophyto feeding threshold for microzoo        xthreshphy  =', xthreshphy
         WRITE(stdout,*) '      poc feeding threshold for microzoo              xthreshpoc  =', xthreshpoc
         WRITE(stdout,*) '      microzoo feeding threshold for microzoo         xthreshzoo  =', xthreshzoo
         WRITE(stdout,*) '      feeding threshold for microzooplankton          xthresh     =', xthresh
         WRITE(stdout,*) '      exsudation rate of microzooplankton             resrat      =', resrat
         WRITE(stdout,*) '      microzooplankton mortality rate                 mzrat       =', mzrat
         WRITE(stdout,*) '      maximal microzoo grazing rate                   grazrat     =', grazrat
         WRITE(stdout,*) '      non assimilated fraction of P by microzoo       unass       =', unass
         WRITE(stdout,*) '      Efficicency of microzoo growth                  epsher      =', epsher
         WRITE(stdout,*) '      Minimum efficicency of microzoo growth          epshermin   =', epshermin
         WRITE(stdout,*) '      Fraction of microzoo excretion as DOM           sigma1      =', sigma1
         WRITE(stdout,*) '      half saturation constant for grazing 1          xkgraz      =', xkgraz
         WRITE(stdout,*) '      Width of the grazing window                     xsigma      =', xsigma
         WRITE(stdout,*) '      Maximum additional width of the grazing window  xsigmadel   =', xsigmadel

      ENDIF
      !
   END SUBROUTINE p4z_micro_init


   !!======================================================================
END MODULE p4zmicro
