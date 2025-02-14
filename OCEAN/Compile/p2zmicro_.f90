










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









MODULE p2zmicro
   !!======================================================================
   !!                         ***  MODULE p2zmicro  ***
   !! TOP :   REDUCED  Compute the sources/sinks for microzooplankton
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!             3.*  !  2025-01  (S. Maishal, R. Person) Change to High Performance
   !!----------------------------------------------------------------------
   !!   p2z_micro      : Compute the sources/sinks for microzooplankton
   !!   p2z_micro_init : Initialize and read the appropriate namelist
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE p2zprod         ! production
   USE p4zsink         ! sedimentation of particles
   USE iom             ! I/O manager
   USE prtctl          ! print control for debugging

   IMPLICIT NONE
   PRIVATE

   !! * Shared module variables
   PUBLIC   p2z_micro         ! called in p2zbio.F90
   PUBLIC   p2z_micro_init    ! called in trcsms_pisces.F90

   REAL(wp), PUBLIC ::   part        !: part of calcite not dissolved in microzoo guts
   REAL(wp), PUBLIC ::   xprefc      !: microzoo preference for POC 
   REAL(wp), PUBLIC ::   xprefn      !: microzoo preference for nanophyto
   REAL(wp), PUBLIC ::   xprefz      !: microzoo preference for microzooplankton
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

   LOGICAL          ::   l_dia_graz, l_dia_lprodz

   !! * Substitutions












































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p2zmicro.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p2z_micro( kt, knt, Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_micro  ***
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
      REAL(wp) :: zcompaz , zcompaph, zcompapoc
      REAL(wp) :: zgraze, zdenom, zfact, zfood, zfoodlim, zbeta
      REAL(wp) :: zepsherf, zepshert, zepsherq, zepsherv, zgrarsig, zgraztotc, zgraztotn
      REAL(wp) :: zgrarem, zgrapoc, zprcaca, zmortz
      REAL(wp) :: zrespz, ztortz, zgrasratn
      REAL(wp) :: zgraznc, zgrazz, zgrazpoc
      REAL(wp) :: ztmp1, ztmp2, ztmp3, ztmptot, zproport, zproport2
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zgrazing, zw3d
      CHARACTER (len=25) :: charout

      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p2z_micro')
      !
      IF( kt == nittrc000 )  THEN
         l_dia_graz  = iom_use( "GRAZ1" ) .OR. iom_use( "MicroZo2" )
         l_dia_graz = l_dia_graz .OR. l_diaadd
      ENDIF

      IF( l_dia_graz )  ALLOCATE( zgrazing(Istrp:Iendp,Jstrp:Jendp,N) ) 
      !
      !$OMP PARALLEL DO & 
      !$OMP PRIVATE(mi, mj, jk, zcompaz, zfact, zproport, zproport2, zrespz, &
        ztortz, zmortz, zcompaph, zcompapoc, zfood, zfoodlim, zdenom, zgraze, &
        ztmp1, ztmp2, ztmp3, ztmptot, zgraznc, zgrazpoc, zgrazz, zgraztotc, &
        zgraztotn, zgrasratn, zepshert, zbeta, zepsherf, zepsherq, zepsherv, &
        zgrarem, zgrapoc, zgrarsig, zprcaca) &
      !$OMP SHARED(tr, tgfunc2, nitrfac, prodpoc, conspoc, prodcal, xfracal, &
        sizen, wsbio3, quotan, grazrat, xstep, 0.5*EPSILON(1.e0), unass, part, &
        l_dia_graz, sigma1, feratz, o2ut, o2nit, rno3)

      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ! Computation of auxiliary variables
         zcompaz = MAX((t(mi, mj,N+1- jk, Kbb,itemp+ntrc_salt+ jpzoo) - 1.e-9), 0.e0)
         zfact = xstep * tgfunc2(mi, mj, jk) * zcompaz

         zproport = MAX(sizen(mi, mj, jk) / 3.0, 1.0)**(-0.48) * &
                   (1.0 - (sizen(mi, mj, jk)**2.0 - 1.0) / 160.0)
         zproport2 = MIN(1.0, (wsbio2 - wsbio3(mi, mj, jk)) / &
                            (wsbio2 - wsbio))

         ! Linear mortality and grazing computation
         zrespz = resrat * zfact * (t(mi, mj,N+1- jk, Kbb,itemp+ntrc_salt+ jpzoo) / &
                 (xkmort + t(mi, mj,N+1- jk, Kbb,itemp+ntrc_salt+ jpzoo)) + &
                  3.0 * nitrfac(mi, mj, jk))

         ztortz = mzrat * 1.e6 * zfact * t(mi, mj,N+1- jk, Kbb,itemp+ntrc_salt+ jpzoo) * &
                 (1.0 - nitrfac(mi, mj, jk))
         zmortz = ztortz + zrespz

         ! Compute total food availability
         zcompaph = zproport * MAX((t(mi, mj,N+1- jk, Kbb,itemp+ntrc_salt+ jpphy) &
                    - xthreshphy), 0.e0)
         zcompapoc = zproport2 * MAX((t(mi, mj,N+1- jk, Kbb,itemp+ntrc_salt+ jppoc) &
                     - xthreshpoc), 0.e0)
         zcompaz = MAX((t(mi, mj,N+1- jk, Kbb,itemp+ntrc_salt+ jpzoo) &
                   - xthreshzoo), 0.e0)
         zfood = xprefn * zcompaph + xprefc * &
                 zcompapoc + xprefz * zcompaz
         zfoodlim = MAX(0.0, zfood - MIN(xthresh, 0.5 * zfood))
         zdenom = zfoodlim / (xkgraz + zfoodlim)
         zgraze = grazrat * xstep * tgfunc2(mi, mj, jk) * &
                  t(mi, mj,N+1- jk, Kbb,itemp+ntrc_salt+ jpzoo) * (1.0 - nitrfac(mi, mj, jk))
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
         ztmp1 = xprefn * zcompaph**2
         ztmp2 = xprefc * zcompapoc**2
         ztmp3 = xprefz * zcompaz**2
         ztmptot = ztmp1 + ztmp2 + ztmp3 + 0.5*EPSILON(1.e0)
         ztmp1 = ztmp1 / ztmptot
         ztmp2 = ztmp2 / ztmptot
         ztmp3 = ztmp3 / ztmptot
          ! Ingestion terms on the different preys of microzooplankton
         zgraznc   = zgraze   * ztmp1 * zdenom  ! Nanophytoplankton
         zgrazpoc  = zgraze   * ztmp2 * zdenom  ! POC
         zgrazz    = zgraze   * ztmp3 * zdenom  ! Microzoo

         ! Ingestion terms on the iron content of the different preys
         ! Total ingestion rate in C, Fe, N units
         zgraztotc = zgraznc + zgrazpoc + zgrazz
         IF( l_dia_graz )   zgrazing(mi,mj,jk) = zgraztotc
         zgraztotn = zgraznc * quotan(mi,mj,jk) + zgrazpoc + zgrazz

         !   Stoichiometruc ratios of the food ingested by zooplanton 
         !   --------------------------------------------------------
         zgrasratn = ( zgraztotn + 0.5*EPSILON(1.e0) ) / ( zgraztotc + 0.5*EPSILON(1.e0) )

         ! Microzooplankton efficiency. 
         ! We adopt a formulation proposed by Mitra et al. (2007)
         ! The gross growth efficiency is controled by the most limiting nutrient.
         ! Growth is also further decreased when the food quality is poor. This is currently
         ! hard coded : it can be decreased by up to 50% (zepsherq)
         ! GGE can also be decreased when food quantity is high, zepsherf (Montagnes and 
         ! Fulton, 2012)
         ! -----------------------------------------------------------------------------
         zepshert  =  MIN( 1., zgrasratn )
         zbeta     =  MAX(0., (epsher - epshermin) )
         ! Food quantity deprivation of the GGE
         zepsherf  = epshermin + zbeta / ( 1.0 + 0.04E6 * 12. * zfood * zbeta )
         ! Food quality deprivation of the GGE
         zepsherq  = 0.5 + (1.0 - 0.5) * zepshert * ( 1.0 + 1.0 ) / ( zepshert + 1.0 )
         ! Actual GGE of microzooplankton
         zepsherv  = zepsherf * zepshert * zepsherq
         ! Excretion of C, N, P
         zgrarem   = zgraztotc * ( 1. - zepsherv - unass ) &
                 &  + ( 1. - epsher - unass ) / ( 1. - epsher ) * ztortz
         ! Egestion of C, N, P
         zgrapoc   = zgraztotc * unass + unass / ( 1. - epsher ) * ztortz + zrespz
         !  Update of the TRA arrays
         !  ------------------------
         ! Fraction of excretion as inorganic nutrients and DIC
         zgrarsig  = zgrarem * sigma1
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpno3) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpno3) + zgrarsig
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) + zgrarem - zgrarsig
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) + zgrarem * feratz
         !
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) - (o2ut + o2nit) &
                                   * zgrarsig
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) + zgrarsig
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) - rno3 * zgrarsig
         !   Update the arrays TRA which contain the biological sources and sinks
         !   --------------------------------------------------------------------
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpzoo) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpzoo) &
                                  - zmortz + zepsherv * zgraztotc - zgrazz 
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpphy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpphy) - zgraznc
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) + zgrapoc - zgrazpoc
         prodpoc(mi,mj,jk) = prodpoc(mi,mj,jk) + zgrapoc
         conspoc(mi,mj,jk) = conspoc(mi,mj,jk) - zgrazpoc
         !
         ! Calcite remineralization due to zooplankton activity
         ! part of the ingested calcite is not dissolving in the acidic gut
         ! ----------------------------------------------------------------
         zprcaca = xfracal(mi,mj,jk) * zgraznc
         ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
         prodcal(mi,mj,jk) = prodcal(mi,mj,jk) + zprcaca * part
         !
         zprcaca = part * zprcaca
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) - zprcaca
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) - 2. * zprcaca
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO

!                      sm update    22-01-2025 end 

      !
      IF( .false. .AND. knt == nrdttrc ) THEN
        !
        IF( l_dia_graz ) THEN  !   Total grazing of phyto by zooplankton
            ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
            DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zw3d(mi,mj,N+1-jk) =  zgrazing(mi,mj,jk) * 1.e+3 * rfact2r * tmask(mi,mj,jk)
            END DO   ;   END DO   ;   END DO
            CALL iom_put( "GRAZ1" , zw3d )  ! conversion in mol/m2/s
            CALL iom_put( "MicroZo2" , zw3d * ( 1. - epsher - unass ) * (-o2ut) * sigma1 ) ! o2 consumption by Microzoo
            DEALLOCATE( zw3d )
        ENDIF
        !
      ENDIF
      !
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         bioFlux(mi,mj,N+1-jk,Ngrapoc) = zgrazing(mi,mj,jk) * 1.e+3 * rfact2r * tmask(mi,mj,jk) !  grazing of phyto by microzoo
         bioFlux(mi,mj,N+1-jk,Nmico2)  = zgrazing(mi,mj,jk) * ( 1. -  epsher - unass ) &
           &                      * (-o2ut) * sigma1 * 1.e+3 * rfact2r * tmask(mi,mj,jk)   ! o2 consumption by Microzoo
      END DO   ;   END DO   ;   END DO
      IF( l_dia_graz )   DEALLOCATE( zgrazing )
      !
      IF(sn_cfctl%l_prttrc) THEN      ! print mean trends (used for debugging)
         WRITE(charout, FMT="('micro')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
  !       CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( .false. )   CALL timing_stop('p2z_micro')
      !
   END SUBROUTINE p2z_micro


   SUBROUTINE p2z_micro_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p2z_micro_init  ***
      !!
      !! ** Purpose :   Initialization of microzooplankton parameters
      !!
      !! ** Method  :   Read the namp2zzoo namelist and check the parameters
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namp2zzoo
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp2zzoo/ part, grazrat, resrat, mzrat, xprefn, xprefc, &
         &                xprefz, xthreshphy, xthreshpoc, xthreshzoo, &
         &                xthresh, xkgraz, epsher, epshermin, sigma1, unass
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*) 
         WRITE(stdout,*) 'p2z_micro_init : Initialization of microzooplankton parameters'
         WRITE(stdout,*) '~~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND(numnatp_ref);READ(numnatp_ref,namp2zzoo,IOSTAT=ios);CALL ctl_nam(ios,"namp2zzoo (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,namp2zzoo,IOSTAT=ios);CALL ctl_nam(ios,"namp2zzoo (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE( numonp, namp2zzoo )
      !
      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*) '   Namelist : namp2zzoo'
         WRITE(stdout,*) '      part of calcite not dissolved in microzoo guts  part        =', part
         WRITE(stdout,*) '      microzoo preference for POC                     xprefc      =', xprefc
         WRITE(stdout,*) '      microzoo preference for nano                    xprefn      =', xprefn
         WRITE(stdout,*) '      microzoo preference for microzooplankton        xprefz      =', xprefz
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
      ENDIF
      !
   END SUBROUTINE p2z_micro_init


   !!======================================================================
END MODULE p2zmicro
