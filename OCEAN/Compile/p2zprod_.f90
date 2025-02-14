










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









MODULE p2zprod
   !!======================================================================
   !!                         ***  MODULE p2zprod  ***
   !! TOP :  Growth Rate of phytoplankton 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-05  (O. Aumont, C. Ethe) New parameterization of light limitation
   !!             3.*  !  2025-01  (S. Maishal, R. Person) Change to High Performance
   !!----------------------------------------------------------------------
   !!   p2z_prod       : Compute the growth Rate of the two phytoplanktons groups
   !!   p2z_prod_init  : Initialization of the parameters for growth
   !!   p2z_prod_alloc : Allocate variables for growth
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE p2zlim          ! Co-limitations of differents nutrients
   USE prtctl          ! print control for debugging
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p2z_prod         ! called in p4zbio.F90
   PUBLIC   p2z_prod_init    ! called in trcsms_pisces.F90
   PUBLIC   p2z_prod_alloc   ! called in trcini_pisces.F90

   REAL(wp), PUBLIC ::   pislopen     !:  P-I slope of nanophytoplankton
   REAL(wp), PUBLIC ::   xadap        !:  Adaptation factor to low light 
   REAL(wp), PUBLIC ::   excretn      !:  Excretion ratio of nanophyto
   REAL(wp), PUBLIC ::   bresp        !:  Basal respiration rate
   REAL(wp), PUBLIC ::   chlcnm       !:  Maximum Chl/C ratio of nano
   REAL(wp), PUBLIC ::   chlcmin      !:  Minimum Chl/C ratio of phytoplankton

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotan   !: proxy of N quota in Nanophyto
   
   REAL(wp) ::   r1_rday    ! 1 / rday
   REAL(wp) ::   texcretn   ! 1 - excretn 

   LOGICAL  :: l_dia_pp, l_dia_mu, l_dia_light, l_dia_lprod

   !! * Substitutions












































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p2zprod.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p2z_prod( kt , knt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_prod  ***
      !!
      !! ** Purpose :   Computes phytoplankton production depending on
      !!                light, temperature and nutrient availability
      !!                Computes also the chlorophyll content 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   !
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      !
      INTEGER  ::   mi, mj, jk
      REAL(wp) ::   znanotot, zpislopen, zfact
      REAL(wp) ::   zlimfac, zlimfac3, zsizetmp, zprodfer, zprprod
      REAL(wp) ::   zprod, zval, zmxl_fac, zmxl_chl
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N) :: zprmax, zmxl
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N) :: zprbio, zprchln
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N) :: zprorcan
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p2z_prod')
      !
      IF( kt == nittrc000 ) THEN
         l_dia_pp    = iom_use( "PPPHYN" ) .OR. iom_use( "TPP"  ) .OR. iom_use( "PPNEWo2" )  &
                       .OR.  iom_use( "THETANANO" ) .OR. iom_use( "CHL" )
         l_dia_mu    = iom_use( "Mumax"  ) .OR. iom_use( "MuN"    )
         l_dia_light = iom_use( "LNlight")
      ENDIF

      ! Initialize the local arrays
      zprorcan(:,:,:) = 0._wp
      zprbio  (:,:,:) = 0._wp
      zmxl    (:,:,:) = 0._wp

      ! Computation of the maximimum production. Based on a Q10 description
      ! of the thermal dependency. Parameters are taken from Bissinger et al. (2008)
      zprmax(:,:,:) = 0.65_wp * r1_rday * tgfunc(:,:,:)

      ! Intermittency is supposed to have a similar effect on production as 
      ! day length (Shatwell et al., 2012). The correcting factor is zmxl_fac. 
      ! zmxl_chl is the fractional day length and is used to compute the mean
      ! PAR during daytime. The effect of mixing is computed using the 
      ! absolute light level definition of the euphotic zone
      ! ------------------------------------------------------------------------- 
      IF ( ln_p4z_dcyc ) THEN    ! Diurnal cycle in 
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk, zval) &
         !$OMP SHARED(etot_ndcy, gdepw, hmld, heup_01, zmxl, 0.5*EPSILON(1.e0))
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            IF( etot_ndcy(mi,mj,jk) > 1.E-3 ) THEN
               zval = 24.0
               IF( gdepw(mi,mj,jk+1,Kmm) <= hbl(mi,mj) ) THEN
                  zval = zval * MIN(1., heup_01(mi,mj) / ( hbl(mi,mj) + 0.5*EPSILON(1.e0) ))
               ENDIF
               zmxl(mi,mj,jk) = zval
            ENDIF
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
      ELSE ! No diurnal cycle in 
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk, zval) &
         !$OMP SHARED(etot_ndcy, gdepw, hmld, heup_01, strn, zmxl, 0.5*EPSILON(1.e0))
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            IF( etot_ndcy(mi,mj,jk) > 1.E-3 ) THEN
               zval = MAX( 1., strn(mi,mj) )
               IF( gdepw(mi,mj,jk+1,Kmm) <= hbl(mi,mj) ) THEN
                  zval = zval * MIN(1., heup_01(mi,mj) / ( hbl(mi,mj) + 0.5*EPSILON(1.e0) ))
               ENDIF
               zmxl(mi,mj,jk) = zval
            ENDIF
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
      ENDIF

      ! The formulation proposed by Geider et al. (1997) has been modified 
      ! to exclude the effect of nutrient limitation and temperature in the PI
      ! curve following Vichi et al. (2007)
      ! -----------------------------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zmxl_fac, zmxl_chl, zprbio, zpislopen) &
      !$OMP SHARED(etot_ndcy, zmxl, zprmax, pislopen, thetanano, &
                   xlimphy, day2sec, 0.5*EPSILON(1.e0), enano, enanom, zprchln)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF( etot_ndcy(mi,mj,jk) > 1.E-3 ) THEN
            zmxl_fac = 1.0 - EXP( -0.26 * zmxl(mi,mj,jk) )
            zmxl_chl = zmxl(mi,mj,jk) / 24.
            !
            zprbio(mi,mj,jk) = zprmax(mi,mj,jk) * zmxl_fac            
            !
            ! The initial slope of the PI curve can be increased for nano
            ! to account for photadaptation, for instance in the DCM
            ! This parameterization is adhoc and should be either 
            ! improved or removed in future versions of the model
            ! Nanophytoplankton
            ! Computation of production function for Carbon
            ! Actual light levels are used here 
            ! ----------------------------------------------
            zpislopen = pislopen * thetanano(mi,mj,jk) / &
              &    ( zprmax(mi,mj,jk) * xlimphy(mi,mj,jk) * zmxl_fac * day2sec + 0.5*EPSILON(1.e0) )
            zprbio(mi,mj,jk)  = zprbio(mi,mj,jk) * ( 1.- EXP( -zpislopen * enano(mi,mj,jk) )  )
            !
            zpislopen = zpislopen * zmxl_fac / ( zmxl_chl + 0.5*EPSILON(1.e0) )
            zprchln(mi,mj,jk) = ( 1.- EXP( -zpislopen * enanom(mi,mj,jk) ) )
         ENDIF
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !  Computation of a proxy of the N/C quota from nutrient limitation 
      !  and light limitation. Steady state is assumed to allow the computation
      !  ----------------------------------------------------------------------
      !thetanano(:,:,:) = chlcnm
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zval, zmxl_fac, thetanano) &
      !$OMP SHARED(xnanono3, zprmax, zprbio, 0.5*EPSILON(1.e0), quotan, &
                   zmxl, chlcnm, pislopen, enanom, xlimphy, &
                   day2sec, chlcmin)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
          zval = xnanono3(mi,mj,jk) * zprmax(mi,mj,jk) / ( zprbio(mi,mj,jk) + 0.5*EPSILON(1.e0) )
          quotan(mi,mj,jk) = MIN( 1., 0.3 + 0.7 * zval )

          ! Diagnostic Chl/C ratio according to Geider et al. (1997)
          ! --------------------------------------------------------
          zmxl_fac = 1.0 - EXP( -0.26 * zmxl(mi,mj,jk) )
          thetanano(mi,mj,jk) = chlcnm / ( 1.0 + pislopen * chlcnm * enanom(mi,mj,jk)   &
             &                  / ( 2.0 * zprmax(mi,mj,jk) * zmxl_fac * xlimphy(mi,mj,jk) * day2sec + 0.5*EPSILON(1.e0) ) )
          thetanano(mi,mj,jk) = MAX( chlcmin, thetanano(mi,mj,jk) )
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      ! Sea-ice effect on production
      ! No production is assumed below sea ice
      ! -------------------------------------- 
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(zprbio, fr_i)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zprbio(mi,mj,jk) = zprbio(mi,mj,jk) * ( 1. - fr_i(mi,mj) )
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      ! Computation of the various production  and nutrient uptake terms
      ! -----------------------------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zlimfac, zlimfac3, zsizetmp) &
      !$OMP SHARED(etot_ndcy, zprorcan, zprbio, xlimphy, tr, &
                   jpphy, Kbb, rfact2, zprchln, zprmax, 0.5*EPSILON(1.e0), &
                   sizena, xsizern)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF( etot_ndcy(mi,mj,jk) > 1.E-3 ) THEN
            !  production terms for nanophyto. (C)
            zprorcan(mi,mj,jk) = zprbio(mi,mj,jk)  * xlimphy(mi,mj,jk) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) * rfact2
            !
            ! Size computation
            ! Size is made a function of the limitation of of phytoplankton growth
            ! Strongly limited cells are supposed to be smaller. sizena is the 
            ! size at time step t+1 and is thus updated at the end of the 
            ! current time step
            ! --------------------------------------------------------------------
         !   zlimfac = xlimphy(mi,mj,jk) * zprchln(mi,mj,jk) / ( zprmax(mi,mj,jk) + 0.5*EPSILON(1.e0) )
            zlimfac = xlimphy(mi,mj,jk) * zprchln(mi,mj,jk) 
            zlimfac3 = zlimfac * zlimfac * zlimfac
            zsizetmp = 1.0 + 1.3 * ( xsizern - 1.0 ) * zlimfac3/(0.3 + zlimfac3)
            sizena(mi,mj,jk) = min(xsizern, max( sizena(mi,mj,jk), zsizetmp ) )
         ENDIF
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !   Update the arrays TRA which contain the biological sources and sinks
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zprodfer) &
      !$OMP SHARED(etot_ndcy, zprorcan, feratz, texcretn, &
                   tr, jpno3, jpfer, consfe3, 0.5*EPSILON(1.e0), plig, &
                   jpdoc, jpoxy, excretn, o2ut, o2nit, rno3, &
                   jpfer, jpphy, jpdic, jptal, Kbb, rfact2)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF( etot_ndcy(mi,mj,jk) > 1.E-3 ) THEN
            zprodfer = zprorcan(mi,mj,jk) * feratz * texcretn
            !
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpno3) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpno3) - zprorcan(mi,mj,jk)
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) &
                    &               - zprorcan(mi,mj,jk) * feratz * texcretn
            consfe3(mi,mj,jk)   = zprodfer * 75.0 / ( 0.5*EPSILON(1.e0) + ( plig(mi,mj,jk) + 75.0 * (1.0 - plig(mi,mj,jk) ) )   &
            &                   * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer) ) / rfact2
            !
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpphy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpphy) + zprorcan(mi,mj,jk) * texcretn
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) + excretn * zprorcan(mi,mj,jk)
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) &
                    &                + ( o2ut + o2nit ) * zprorcan(mi,mj,jk)
            !
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) - zprorcan(mi,mj,jk)
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) + rno3 * zprorcan(mi,mj,jk)
         ENDIF
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
    ! Total primary production per year
    IF( l_dia_pp )  THEN
       ALLOCATE( zw3d(Istrp:Iendp,Jstrp:Jendp,N) )  ;  zw3d(Istrp:Iendp,Jstrp:Jendp,N) = 0._wp
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(zprorcan, cvol, zw3d)
       DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
          zw3d(mi,mj,jk) = zprorcan(mi,mj,jk) * cvol(mi,mj,jk)
       END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
       tpp = glob_sum( 'p2zprod', zw3d )
       DEALLOCATE ( zw3d )
    ENDIF

    IF( .false. .AND.  knt == nrdttrc ) THEN
       !
       zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
       IF( l_dia_pp ) THEN
          ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, N+1-jk) &
         !$OMP SHARED(thetanano, zfact, tmask, zw3d)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = thetanano(mi,mj,jk) * zfact * tmask(mi,mj,jk) 
          END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
          CALL iom_put( "THETANANO", zw3d ) ! Diagnostic Chl:C ratio
          !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, N+1-jk) &
         !$OMP SHARED(thetanano, tr, tmask, zw3d)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = thetanano(mi,mj,jk) * 12. &
               &         * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) * 1.0e+6 * tmask(mi,mj,jk) 
          END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
          CALL iom_put( "CHL", zw3d ) ! total Chloropyll
          !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, N+1-jk) &
         !$OMP SHARED(zprorcan, zfact, tmask, zw3d)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = zprorcan(mi,mj,jk) * zfact * tmask(mi,mj,jk) 
          END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
          CALL iom_put( "PPPHYN", zw3d )  ! primary production by nanophyto
          CALL iom_put( "TPP", zw3d ) ! total primary production
          CALL iom_put( "PPNEWo2", ( o2ut + o2nit ) * zw3d ) ! Oxygen production by the New Produc
          CALL iom_put( "tintpp"  , tpp * zfact )  !  global total integrated primary production molC/s
          DEALLOCATE ( zw3d )
       ENDIF
       !
       IF( l_dia_mu ) THEN
          ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, N+1-jk) &
         !$OMP SHARED(zprmax, tmask, zw3d)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = zprmax(mi,mj,jk) * tmask(mi,mj,jk) 
          END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
          CALL iom_put( "Mumax", zw3d )
          ! Realized growth rate for nanophyto
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, N+1-jk) &
         !$OMP SHARED(zprbio, xlimphy, tmask, zw3d)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = zprbio(mi,mj,jk) * xlimphy(mi,mj,jk) &
                 &          * tmask(mi,mj,jk) 
          END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
          CALL iom_put( "MuN", zw3d )
          DEALLOCATE ( zw3d )
       ENDIF
       !
       IF( l_dia_light ) THEN
          ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
          ! light limitation term for nano
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, N+1-jk) &
         !$OMP SHARED(zprbio, zprmax, 0.5*EPSILON(1.e0), tmask, zw3d)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = zprbio(mi,mj,jk) &
                 &          / ( zprmax(mi,mj,jk) + 0.5*EPSILON(1.e0) ) &
                 &          * tmask(mi,mj,jk) 
          END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
          CALL iom_put( "LNlight", zw3d )
          DEALLOCATE ( zw3d )
       ENDIF
       !
     ENDIF      

      !   Supplementary diagnostics
     zfact = 1.e3 * rfact2r
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(mi, mj, N+1-jk) &
    !$OMP SHARED(zprorcan, zfact, tmask, trc3d, &
                 thetanano, tr, o2ut, o2nit)
     DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
        bioFlux(mi,mj,N+1-jk,Nprorca  ) = zprorcan(mi,mj,jk) * zfact * tmask(mi,mj,jk)  ! primary production by nanophyto
        bioFlux(mi,mj,N+1-jk,Npronew  ) = thetanano(mi,mj,jk) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) * 1.0e+6 * tmask(mi,mj,jk) ! Total chloro.
        bioFlux(mi,mj,N+1-jk,Npronewo2) = ( o2ut + o2nit ) * zprorcan(mi,mj,jk) * zfact * tmask(mi,mj,jk) ! Oxygen production by the New Produc.
     END DO   ;   END DO   ;   END DO
    !$OMP END PARALLEL DO
     IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('prod')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
  !       CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
     ENDIF
      !
      IF( .false. )  CALL timing_stop('p2z_prod')
      !
   END SUBROUTINE p2z_prod


   SUBROUTINE p2z_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p2z_prod_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton production parameters
      !!
      !! ** Method  :   Read the namp2zprod namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namp2zprod
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      ! Namelist block
      NAMELIST/namp2zprod/ pislopen, bresp, excretn,  &
         &                 chlcnm, chlcmin
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*)
         WRITE(stdout,*) 'p2z_prod_init : phytoplankton growth'
         WRITE(stdout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND(numnatp_ref);READ(numnatp_ref,namp2zprod,IOSTAT=ios);CALL ctl_nam(ios,"namp2zprod (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,namp2zprod,IOSTAT=ios);CALL ctl_nam(ios,"namp2zprod (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE( numonp, namp2zprod )

      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*) '   Namelist : namp2zprod'
         WRITE(stdout,*) '      P-I slope                                 pislopen     =', pislopen
         WRITE(stdout,*) '      excretion ratio of nanophytoplankton      excretn      =', excretn
         WRITE(stdout,*) '      basal respiration in phytoplankton        bresp        =', bresp
         WRITE(stdout,*) '      Maximum Chl/C in phytoplankton            chlcmin      =', chlcmin
         WRITE(stdout,*) '      Minimum Chl/C in nanophytoplankton        chlcnm       =', chlcnm
      ENDIF
      !
      r1_rday   = 1._wp / day2sec 
      texcretn  = 1._wp - excretn
      tpp       = 0._wp
      !
   END SUBROUTINE p2z_prod_init

   INTEGER FUNCTION p2z_prod_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_prod_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( quotan(Istrp:Iendp,Jstrp:Jendp,N), STAT = p2z_prod_alloc )
      !
      IF( p2z_prod_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p2z_prod_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p2z_prod_alloc



   !!======================================================================
END MODULE p2zprod
