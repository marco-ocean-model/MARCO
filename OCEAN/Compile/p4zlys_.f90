










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









MODULE p4zlys
   !!======================================================================
   !!                         ***  MODULE p4zlys  ***
   !! TOP :    
   !!======================================================================
   !! History :    -   !  1988-07  (E. MAIER-REIMER) Original code
   !!              -   !  1998     (O. Aumont) additions
   !!              -   !  1999     (C. Le Quere) modifications
   !!             1.0  !  2004     (O. Aumont) modifications
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J. Orr)  Calcon salinity dependence
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Improvment of calcite dissolution
   !!             3.6  !  2015-05  (O. Aumont)  quota
   !!             4.2  !  2020     (J. ORR )  rhop is replaced by "in situ density" rhd
   !!             4.*  !  2025-02  (S. Maishal, R. Person)  and optimization
   !!----------------------------------------------------------------------
   !!   p4z_lys        :   Compute the CaCO3 dissolution 
   !!   p4z_lys_init   :   Read the namelist parameters
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !   Source Minus Sink variables
   USE p4zche          !  Chemical model
   USE p4zsink         !  sinking of particles
   USE prtctl          !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p2z_lys         ! called in trcsms_pisces.F90
   PUBLIC   p4z_lys         ! called in trcsms_pisces.F90
   PUBLIC   p4z_lys_init    ! called in trcsms_pisces.F90

   REAL(wp), PUBLIC ::   kdca   !: diss. rate constant calcite
   REAL(wp), PUBLIC ::   nca    !: order of reaction for calcite dissolution

   INTEGER  :: rmtss            ! number of seconds per month 
   LOGICAL  :: l_dia

   !! * Module variables
   REAL(wp) :: calcon = 1.03E-2           !: mean calcite concentration [Ca2+]  in sea water [mole/kg solution]
   ! J. ORR: Made consistent with mocsy's choice based on literature review from Munhoven
   !   REAL(wp) :: calcon = 1.0287E-2     !: mean calcite concentration [Ca2+] in sea water [mole/kg solution]

   !! * Substitutions
   !! * Substitutions












































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zlys.F90 15532 2021-11-24 11:47:32Z techene $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p2z_lys( kt, knt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_lys  ***
      !!
      !! ** Purpose :   CALCULATES DEGREE OF CACO3 SATURATION IN THE WATER
      !!                COLUMN, DISSOLUTION/PRECIPITATION OF CACO3 AND LOSS
      !!                OF CACO3 TO THE CACO3 SEDIMENT POOL.
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt        ! ocean time step and ???
      INTEGER, INTENT(in)  ::  Kbb, Kmm, Krhs ! time level indices
      !
      INTEGER  ::   mi, mj, jk, jn
      REAL(wp) ::   zdispot, zfact, zcalcon, zdepexp, zdissol
      REAL(wp) ::   zomegaca, zexcess, zexcess0, zkd, zwsbio
      CHARACTER (len=25) ::   charout
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N) :: zhinit, zhi, zco3, zcaco3, ztra
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)  :: zw3d, zcaldiss
      INTEGER :: ik100
      !!---------------------------------------------------------------------
      !
      IF( .false. )  CALL timing_start('p2z_lys')
      !
      IF( kt == nittrc000 )  &
           & l_dia = iom_use( "PH" ) .OR. iom_use( "CO3" ) .OR. iom_use( "CO3sat" ) &
           &    .OR. iom_use( "DCAL" ) .OR. iom_use( "PCAL" )  &
           &    .OR. iom_use( "EPC100" ) .OR. iom_use( "EXPCAL" )
      !
      IF( l_dia )   THEN                  !* Save ta and sa trends
         ALLOCATE( zcaldiss(Istrp:Iendp,Jstrp:Jendp,N) )  
         zcaldiss(Istrp:Iendp,Jstrp:Jendp,:) = 0._wp   
         zco3(Istrp:Iendp,Jstrp:Jendp,N) = 0._wp
      ENDIF
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(zhinit, hi, rhop, 0.5*EPSILON(1.e0))
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zhinit(mi,mj,jk) = hi(mi,mj,jk) * 1000._wp / ( (rho0+rho1(mi,mj,N+1-jk)) + 0.5*EPSILON(1.e0) )
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      !     -------------------------------------------
      !     COMPUTE [CO3--] and [H+] CONCENTRATIONS
      !     -------------------------------------------

      CALL solve_at_general( zhinit, zhi, Kbb )
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zco3, zhi) &
      !$OMP SHARED(tr, jpdic, Kbb, ak13, ak23, hi, rhop, 0.5*EPSILON(1.e0))
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zco3(mi,mj,jk) = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdic) &
            &             * ak13(mi,mj,jk) * ak23(mi,mj,jk) / (zhi(mi,mj,jk)**2 &
            &             + ak13(mi,mj,jk) * zhi(mi,mj,jk) + ak13(mi,mj,jk) &
            &             * ak23(mi,mj,jk) + 0.5*EPSILON(1.e0) )
         hi (mi,mj,jk) = zhi(mi,mj,jk) * (rho0+rho1(mi,mj,N+1-jk)) / 1000._wp
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !     ---------------------------------------------------------
      !        CALCULATE DEGREE OF CACO3 SATURATION AND CORRESPONDING
      !        DISSOLOUTION AND PRECIPITATION OF CACO3 (BE AWARE OF
      !        MGCO3)
      !     ---------------------------------------------------------

      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zcalcon, zfact, zomegaca, excess, &
                    zexcess0, zexcess, zdispot, zkd, ztra) &
      !$OMP SHARED(salinprac, rhop, calcon, aksp, 0.5*EPSILON(1.e0), nca, kdca, rmtss)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp

         ! DEVIATION OF [CO3--] FROM SATURATION VALUE
         ! Salinity dependance in zomegaca and divide by rhop to have good units
         zcalcon  = calcon * ( salinprac(mi,mj,jk) / 35._wp )
         zfact    = (rho0+rho1(mi,mj,N+1-jk)) / 1000._wp
         zomegaca = ( zcalcon * zco3(mi,mj,jk) ) / ( aksp(mi,mj,jk) * zfact + 0.5*EPSILON(1.e0) )

         ! SET DEGREE OF UNDER-/SUPERSATURATION
         excess(mi,mj,jk) = 1._wp - zomegaca
         zexcess0 = MAX( 0., excess(mi,mj,jk) )

         IF( zomegaca < 0.8 ) THEN
            zexcess = zexcess0**nca
            ! AMOUNT OF CACO3 THAT RE-ENTERS SOLUTION
            zdispot = kdca * zexcess 
         ELSE
            zkd = kdca * 0.2**(nca - 0.2)
            zexcess = zexcess0**0.2
            zdispot = zkd * zexcess
        ENDIF

        !  CHANGE OF [CO3--] , [ALK], PARTICULATE [CACO3],
        !  AND [SUM(CO2)] DUE TO CACO3 DISSOLUTION/PRECIPITATION
        ztra(mi,mj,jk)  = zdispot / rmtss ! calcite dissolution
        !
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      zcaco3(:,:,:) = 0._wp
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj) &
      !$OMP SHARED(zcaco3, prodcal, rfact2r, wsbio4, e3t, day2sec, ztra, Kmm)
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zcaco3(mi,mj,1) = prodcal(mi,mj,1) * rfact2r &
             &        / ( wsbio4(mi,mj,1) / Hz(mi,mj,N+1-1) / day2sec + ztra(mi,mj,1) )
      END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zdissol, zwsbio, zdepexp) &
      !$OMP SHARED(zcaco3, prodcal, rfact2r, wsbio4, day2sec, tmask, &
                   ztra, e3t, Kmm, tr, Krhs, l_dia, zcaldiss)
      DO jk= 2, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zdissol = 0.0
         zwsbio = wsbio4(mi,mj,jk) / day2sec
         IF ( tmask(mi,mj,1) == 1. ) THEN
            IF ( ztra(mi,mj,jk) == 0.0 ) THEN
               zcaco3(mi,mj,jk) = zcaco3(mi,mj,jk-1) &
                 &      + prodcal(mi,mj,jk) * rfact2r / zwsbio * Hz(mi,mj,N+1-jk)
            ELSE
               zdepexp = exp( - ztra(mi,mj,jk) * Hz(mi,mj,N+1-jk) / zwsbio )
               zcaco3(mi,mj,jk) = prodcal(mi,mj,jk) * rfact2r / ztra(mi,mj,jk) &
                  & * (1.0 - zdepexp ) + zcaco3(mi,mj,jk-1) * zdepexp
               zdissol = prodcal(mi,mj,jk) * Hz(mi,mj,N+1-jk) + prodcal(mi,mj,jk) &
                  &      * zwsbio / ztra(mi,mj,jk) * ( zdepexp - 1.0 ) &
                  &      + zwsbio * zcaco3(mi,mj,jk-1) * ( 1.0 - zdepexp ) * rfact2
               zdissol = zdissol / Hz(mi,mj,N+1-jk)
            ENDIF
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) + zdissol
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) + 2.0 * zdissol
            !
            IF( l_dia ) zcaldiss(mi,mj,jk) = zdissol
            !
         ENDIF
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj) &
      !$OMP SHARED(sinkcalb, wsbio4, mbkt, zcaco3, rfact2, day2sec)
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         sinkcalb(mi,mj) = wsbio4(mi,mj,N) &
            &            * zcaco3(mi,mj,N) * rfact2 / day2sec
      END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('lys ')")
        CALL prt_ctl_info( charout, cdcomp = 'top' )
 !       CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( l_dia .AND. knt == nrdttrc ) THEN
         ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,1:N) )   ;   zw3d(:,:,:) = 0.
         zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, N+1-jk, zfact) &
         !$OMP SHARED(zw3d, prodcal, tmask)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zw3d(mi,mj,N+1-jk) = prodcal(mi,mj,jk) * zfact * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
         CALL iom_put( "PCAL", zw3d )   ! calcite production
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, N+1-jk, zfact) &
         !$OMP SHARED(zw3d, zcaldiss, tmask)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zw3d(mi,mj,N+1-jk) = zcaldiss(mi,mj,jk) * zfact * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
         CALL iom_put( "DCAL", zw3d )  ! calcite dissolution
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, N+1-jk) &
         !$OMP SHARED(zw3d, wsbio4, zcaco3, day2sec, tmask)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zw3d(mi,mj,N+1-jk) = wsbio4(mi,mj,jk) * zcaco3(mi,mj,jk) &
              &       * 1.e+3 / day2sec  * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
         !
         CALL iom_put( "EPCAL100", zw3d(:,:,ik100) ) ! Export of calcite at 100m
         CALL iom_put( "EXPCAL" , zw3d ) ! Export of calcite in the water column
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, N+1-jk) &
         !$OMP SHARED(zw3d, hi, 0.5*EPSILON(1.e0), tmask)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zw3d(mi,mj,N+1-jk) = -1. * LOG10( MAX( hi(mi,mj,jk),0.5*EPSILON(1.e0) ) ) & 
                &                  * tmask(mi,mj,jk)         
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
         CALL iom_put( "PH" , zw3d )  ! PH
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, N+1-jk) &
         !$OMP SHARED(zw3d, zco3, tmask)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zw3d(mi,mj,N+1-jk) = zco3(mi,mj,jk) * 1.e+3 * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
         CALL iom_put( "CO3", zw3d ) ! ion carbonate
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, N+1-jk, zcalcon, zfact) &
         !$OMP SHARED(zw3d, aksp, rhop, calcon, salinprac, tmask, 0.5*EPSILON(1.e0))
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zcalcon         = calcon * ( salinprac(mi,mj,jk) / 35._wp )
            zfact           = (rho0+rho1(mi,mj,N+1-jk)) / 1000._wp
            zw3d(mi,mj,N+1-jk) = aksp(mi,mj,jk) * zfact &
              &             / ( zcalcon + 0.5*EPSILON(1.e0) ) * 1.e+3 * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
         !
         CALL iom_put( "CO3sat", zw3d )  ! calcite saturation
         !
         DEALLOCATE( zcaldiss, zw3d )
      ENDIF
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, N+1-jk, zcalcon, zfact) &
      !$OMP SHARED(trc3d, hi, tmask, zco3, aksp, rhop, calcon, salinprac, 0.5*EPSILON(1.e0))
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         bioFlux(mi,mj,N+1-jk,Nhi    ) = -1. * LOG10( MAX( hi(mi,mj,jk), + 1.0e-12, 0.5*EPSILON(1.e0)) ) * tmask(mi,mj,jk) ! PH
         bioFlux(mi,mj,N+1-jk,NCo3   ) = zco3(mi,mj,jk)  * 1.e+3 * tmask(mi,mj,jk) ! Ion carbonate
         bioFlux(mi,mj,N+1-jk,Naksp) = aksp(mi,mj,jk) * (rho0+rho1(mi,mj,N+1-jk)) / 1000._wp &
              &   / calcon * ( salinprac(mi,mj,jk) / 35._wp ) * tmask(mi,mj,jk) ! calcite saturation
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      IF( .false. )   CALL timing_stop('p2z_lys')
      !
   END SUBROUTINE p2z_lys

   SUBROUTINE p4z_lys( kt, knt, Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_lys  ***
      !!
      !! ** Purpose :   CALCULATES DEGREE OF CACO3 SATURATION IN THE WATER
      !!                COLUMN, DISSOLUTION/PRECIPITATION OF CACO3 AND LOSS
      !!                OF CACO3 TO THE CACO3 SEDIMENT POOL.
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step and ???
      INTEGER, INTENT(in)  ::  Kbb, Krhs ! time level indices
      !
      INTEGER  ::   mi, mj, jk, jn
      REAL(wp) ::   zdispot, zfact, zcalcon, ztra, zdissol
      REAL(wp) ::   zomegaca, zexcess, zexcess0, zkd
      CHARACTER (len=25) ::   charout
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N) :: zhinit, zhi, zco3
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)  :: zw3d, zcaldiss
      !!---------------------------------------------------------------------
      !
      IF( .false. )  CALL timing_start('p4z_lys')
      !
     IF( kt == nittrc000 ) &
           & l_dia = iom_use( "PH" ) .OR. iom_use( "CO3" ) .OR. iom_use( "CO3sat" ) &
           &                  .OR. iom_use( "DCAL" ) .OR. iom_use( "PCAL" )

      IF( l_dia )   THEN
         ALLOCATE( zcaldiss(Istrp:Iendp,Jstrp:Jendp,N) )
         zcaldiss(Istrp:Iendp,Jstrp:Jendp,:) = 0._wp
         zco3(Istrp:Iendp,Jstrp:Jendp,N) = 0._wp
      ENDIF
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(zhinit, hi, rhop, 0.5*EPSILON(1.e0))
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zhinit(mi,mj,jk) = hi(mi,mj,jk) * 1000._wp / ( (rho0+rho1(mi,mj,N+1-jk)) + 0.5*EPSILON(1.e0) )
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      !     -------------------------------------------
      !     COMPUTE [CO3--] and [H+] CONCENTRATIONS
      !     -------------------------------------------

      CALL solve_at_general( zhinit, zhi, Kbb )
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(zco3, tr, jpdic, Kbb, ak13, ak23, zhi, 0.5*EPSILON(1.e0), hi, rhop)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zco3(mi,mj,jk) = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdic) * ak13(mi,mj,jk) * ak23(mi,mj,jk) &
                          / (zhi(mi,mj,jk)**2 + ak13(mi,mj,jk) * zhi(mi,mj,jk) &
                          + ak13(mi,mj,jk) * ak23(mi,mj,jk) + 0.5*EPSILON(1.e0) )
         hi  (mi,mj,jk) = zhi(mi,mj,jk) * (rho0+rho1(mi,mj,N+1-jk)) / 1000._wp
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !     ---------------------------------------------------------
      !        CALCULATE DEGREE OF CACO3 SATURATION AND CORRESPONDING
      !        DISSOLOUTION AND PRECIPITATION OF CACO3 (BE AWARE OF
      !        MGCO3)
      !     ---------------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zcalcon, zfact, zomegaca, excess, &
                           zexcess0, zexcess, zdispot, zkd, ztra) &
      !$OMP SHARED(tr, jptal, jpcal, jpdic, Krhs, zco3, aksp, 0.5*EPSILON(1.e0), &
                           rhop, calcon, salinprac, nca, kdca, rfact2, &
                           rmtss, l_dia, zcaldiss)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp

         ! DEVIATION OF [CO3--] FROM SATURATION VALUE
         ! Salinity dependance in zomegaca and divide by rhd to have good units
         zcalcon  = calcon * ( salinprac(mi,mj,jk) / 35._wp )
         zfact    = (rho0+rho1(mi,mj,N+1-jk)) / 1000._wp
         zomegaca = ( zcalcon * zco3(mi,mj,jk) ) / ( aksp(mi,mj,jk) * zfact + 0.5*EPSILON(1.e0) )

         ! SET DEGREE OF UNDER-/SUPERSATURATION
         excess(mi,mj,jk) = 1._wp - zomegaca
         zexcess0 = MAX( 0., excess(mi,mj,jk) )

         IF( zomegaca < 0.8 ) THEN
            zexcess = zexcess0**nca
            ! AMOUNT CACO3 THAT RE-ENTERS SOLUTION
            zdispot = kdca * zexcess * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpcal)
         ELSE
            zkd = kdca * 0.2**(nca - 0.2)
            zexcess = zexcess0**0.2
            zdispot = zkd * zexcess * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpcal)
        ENDIF

        !  CHANGE OF [CO3--] , [ALK], PARTICULATE [CACO3],
        !       AND [SUM(CO2)] DUE TO CACO3 DISSOLUTION/PRECIPITATION
        ztra  = zdispot * rfact2 / rmtss ! calcite dissolution
        !
        t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) + 2. * ztra
        t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpcal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpcal) -      ztra
        t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) +      ztra
        !
        IF( l_dia ) zcaldiss(mi,mj,jk) = zdissol
        !
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      IF( l_dia .AND. knt == nrdttrc ) THEN
         ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,1:N) )  ;  zw3d(:,:,:) = 0._wp
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk, zfact) &
         !$OMP SHARED(zw3d, prodcal, tmask, N)
         zfact = 1.e+3 * rfact2r  ! conversion from mol/l/kt to  mol/m3/s
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zw3d(mi,mj,N+1-jk) = prodcal(mi,mj,jk) * zfact * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
         CALL iom_put( "PCAL", zw3d ) ! calcite production
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk, zfact) &
         !$OMP SHARED(zw3d, zcaldiss, tmask, N)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zw3d(mi,mj,N+1-jk) = zcaldiss(mi,mj,jk) * zfact * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
         CALL iom_put( "DCAL", zw3d ) ! calcite dissolution
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk) &
         !$OMP SHARED(zw3d, hi, 0.5*EPSILON(1.e0), tmask, N)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zw3d(mi,mj,N+1-jk) = -1. * LOG10( MAX( hi(mi,mj,jk),0.5*EPSILON(1.e0) ) ) &
                &                 * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
         CALL iom_put( "PH" , zw3d ) ! PH
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk) &
         !$OMP SHARED(zw3d, zco3, tmask, N)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zw3d(mi,mj,N+1-jk) = zco3(mi,mj,jk) * 1.e+3 * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
         CALL iom_put( "CO3", zw3d ) ! ion carbonate
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk, zcalcon, zfact) &
         !$OMP SHARED(zw3d, aksp, rhop, salinprac, calcon, 0.5*EPSILON(1.e0), tmask, N)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zcalcon        = calcon * ( salinprac(mi,mj,jk) / 35._wp )
             zfact          = (rho0+rho1(mi,mj,N+1-jk)) / 1000._wp
             zw3d(mi,mj,N+1-jk) = aksp(mi,mj,jk) * zfact &
                 &        / ( zcalcon + 0.5*EPSILON(1.e0) )  * 1.e+3 * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
         CALL iom_put( "CO3sat", zw3d )  ! calcite saturation
         !
         DEALLOCATE( zcaldiss, zw3d )
      ENDIF
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(trc3d, hi, zco3, aksp, rhop, &
                   calcon, salinprac, 0.5*EPSILON(1.e0), tmask, N)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ! PH
         bioFlux(mi,mj,N+1-jk,Nhi) = -1. * LOG10( MAX( hi(mi,mj,jk), + 1.0e-12, 0.5*EPSILON(1.e0)) ) &
                                      * tmask(mi,mj,jk)
         ! Ion carbonate
         bioFlux(mi,mj,N+1-jk,NCo3) = zco3(mi,mj,jk) * 1.e+3 &
                                      * tmask(mi,mj,jk)
         ! calcite saturation
         bioFlux(mi,mj,N+1-jk,Naksp) = aksp(mi,mj,jk) &
                                      * (rho0+rho1(mi,mj,N+1-jk)) / 1000._wp / calcon &
                                      * ( salinprac(mi,mj,jk) / 35._wp ) &
                                      * tmask(mi,mj,jk)
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('lys ')")
        CALL prt_ctl_info( charout, cdcomp = 'top' )
!        CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( .false. )   CALL timing_stop('p4z_lys')
      !
   END SUBROUTINE p4z_lys


   SUBROUTINE p4z_lys_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_lys_init  ***
      !!
      !! ** Purpose :   Initialization of CaCO3 dissolution parameters
      !!
      !! ** Method  :   Read the nampiscal namelist and check the parameters
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampiscal
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios  ! Local integer
      !
      NAMELIST/nampiscal/ kdca, nca
      !!----------------------------------------------------------------------
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*)
         WRITE(stdout,*) 'p4z_lys_init : initialization of CaCO3 dissolution'
         WRITE(stdout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND(numnatp_ref);READ(numnatp_ref,nampiscal,IOSTAT=ios);CALL ctl_nam(ios,"nampiscal (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,nampiscal,IOSTAT=ios);CALL ctl_nam(ios,"nampiscal (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE( numonp, nampiscal )
      !
      IF(mynode .eq. 0) THEN          ! control print
         WRITE(stdout,*) '   Namelist : nampiscal'
         WRITE(stdout,*) '   diss. rate constant calcite (per month)        kdca =', kdca
         WRITE(stdout,*) '   order of reaction for calcite dissolution      nca  =', nca
      ENDIF
      !
      ! Number of seconds per month 
      rmtss =  year2day * day2sec / 12.
      !
      ! CE not really needed ; temporary, shoud be removed when quotan( Istrp:Iendp,Jstrp:Jendp,N )
      excess(:,:,:) = 0._wp
      !
   END SUBROUTINE p4z_lys_init


   !!======================================================================
END MODULE p4zlys
