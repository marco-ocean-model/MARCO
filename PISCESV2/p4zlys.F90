#include "cppdefs.h"

MODULE p4zlys
   !!======================================================================
   !!                         ***  MODULE p4zlys  ***
   !! TOP :   PISCES 
   !!======================================================================
   !! History :    -   !  1988-07  (E. MAIER-REIMER) Original code
   !!              -   !  1998     (O. Aumont) additions
   !!              -   !  1999     (C. Le Quere) modifications
   !!             1.0  !  2004     (O. Aumont) modifications
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J. Orr)  Calcon salinity dependence
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Improvment of calcite dissolution
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!             4.2  !  2020     (J. ORR )  rhop is replaced by "in situ density" rhd
   !!             4.*  !  2025-02  (S. Maishal, R. Person) MPI and optimization
#if defined key_pisces
   !!----------------------------------------------------------------------
   !!   p4z_lys        :   Compute the CaCO3 dissolution 
   !!   p4z_lys_init   :   Read the namelist parameters
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
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
#  include "ocean2pisces.h90"
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"
#  include "read_nml_substitute.h90"

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
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zdispot, zfact, zcalcon, zdepexp, zdissol
      REAL(wp) ::   zomegaca, zexcess, zexcess0, zkd, zwsbio
      CHARACTER (len=25) ::   charout
      REAL(wp), DIMENSION(A2D(0),jpk) :: zhinit, zhi, zco3, zcaco3, ztra
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)  :: zw3d, zcaldiss
      INTEGER :: ik100
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('p2z_lys')
      !
      IF( kt == nittrc000 )  &
           & l_dia = iom_use( "PH" ) .OR. iom_use( "CO3" ) .OR. iom_use( "CO3sat" ) &
           &    .OR. iom_use( "DCAL" ) .OR. iom_use( "PCAL" )  &
           &    .OR. iom_use( "EPC100" ) .OR. iom_use( "EXPCAL" )
      !
      IF( l_dia )   THEN                  !* Save ta and sa trends
         ALLOCATE( zcaldiss(A2D(0),jpk) )  
         zcaldiss(A2D(0),:) = 0._wp   
         zco3(A2D(0),jpk) = 0._wp
      ENDIF
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj, jk) &
      !$OMP SHARED(zhinit, hi, rhop, rtrn)
      DO_3D( 0, 0, 0, 0, 1, jpkm1)
         zhinit(ji,jj,jk) = hi(ji,jj,jk) * 1000._wp / ( rhop(ji,jj,jk) + rtrn )
      END_3D
      !$OMP END PARALLEL DO
      !
      !     -------------------------------------------
      !     COMPUTE [CO3--] and [H+] CONCENTRATIONS
      !     -------------------------------------------

      CALL solve_at_general( zhinit, zhi, Kbb )
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj, jk, zco3, zhi) &
      !$OMP SHARED(tr, jpdic, Kbb, ak13, ak23, hi, rhop, rtrn)
      DO_3D( 0, 0, 0, 0, 1, jpkm1)
         zco3(ji,jj,jk) = tr(ji,jj,jk,jpdic,Kbb) &
            &             * ak13(ji,jj,jk) * ak23(ji,jj,jk) / (zhi(ji,jj,jk)**2 &
            &             + ak13(ji,jj,jk) * zhi(ji,jj,jk) + ak13(ji,jj,jk) &
            &             * ak23(ji,jj,jk) + rtrn )
         hi (ji,jj,jk) = zhi(ji,jj,jk) * rhop(ji,jj,jk) / 1000._wp
      END_3D
      !$OMP END PARALLEL DO
      !     ---------------------------------------------------------
      !        CALCULATE DEGREE OF CACO3 SATURATION AND CORRESPONDING
      !        DISSOLOUTION AND PRECIPITATION OF CACO3 (BE AWARE OF
      !        MGCO3)
      !     ---------------------------------------------------------

      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj, jk, zcalcon, zfact, zomegaca, excess, &
                    zexcess0, zexcess, zdispot, zkd, ztra) &
      !$OMP SHARED(salinprac, rhop, calcon, aksp, rtrn, nca, kdca, rmtss)
      DO_3D( 0, 0, 0, 0, 1, jpkm1)

         ! DEVIATION OF [CO3--] FROM SATURATION VALUE
         ! Salinity dependance in zomegaca and divide by rhop to have good units
         zcalcon  = calcon * ( salinprac(ji,jj,jk) / 35._wp )
         zfact    = rhop(ji,jj,jk) / 1000._wp
         zomegaca = ( zcalcon * zco3(ji,jj,jk) ) / ( aksp(ji,jj,jk) * zfact + rtrn )

         ! SET DEGREE OF UNDER-/SUPERSATURATION
         excess(ji,jj,jk) = 1._wp - zomegaca
         zexcess0 = MAX( 0., excess(ji,jj,jk) )

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
        ztra(ji,jj,jk)  = zdispot / rmtss ! calcite dissolution
        !
      END_3D
      !$OMP END PARALLEL DO
      !
      zcaco3(:,:,:) = 0._wp
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj) &
      !$OMP SHARED(zcaco3, prodcal, rfact2r, wsbio4, e3t, rday, ztra, Kmm)
      DO_2D( 0, 0, 0, 0 )
         zcaco3(ji,jj,1) = prodcal(ji,jj,1) * rfact2r &
             &        / ( wsbio4(ji,jj,1) / e3t(ji,jj,1,Kmm) / rday + ztra(ji,jj,1) )
      END_2D
      !$OMP END PARALLEL DO
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj, jk, zdissol, zwsbio, zdepexp) &
      !$OMP SHARED(zcaco3, prodcal, rfact2r, wsbio4, rday, tmask, &
                   ztra, e3t, Kmm, tr, Krhs, l_dia, zcaldiss)
      DO_3D( 0, 0, 0, 0, 2, jpkm1)
         zdissol = 0.0
         zwsbio = wsbio4(ji,jj,jk) / rday
         IF ( tmask(ji,jj,1) == 1. ) THEN
            IF ( ztra(ji,jj,jk) == 0.0 ) THEN
               zcaco3(ji,jj,jk) = zcaco3(ji,jj,jk-1) &
                 &      + prodcal(ji,jj,jk) * rfact2r / zwsbio * e3t(ji,jj,jk,Kmm)
            ELSE
               zdepexp = exp( - ztra(ji,jj,jk) * e3t(ji,jj,jk,Kmm) / zwsbio )
               zcaco3(ji,jj,jk) = prodcal(ji,jj,jk) * rfact2r / ztra(ji,jj,jk) &
                  & * (1.0 - zdepexp ) + zcaco3(ji,jj,jk-1) * zdepexp
               zdissol = prodcal(ji,jj,jk) * e3t(ji,jj,jk,Kmm) + prodcal(ji,jj,jk) &
                  &      * zwsbio / ztra(ji,jj,jk) * ( zdepexp - 1.0 ) &
                  &      + zwsbio * zcaco3(ji,jj,jk-1) * ( 1.0 - zdepexp ) * rfact2
               zdissol = zdissol / e3t(ji,jj,jk,Kmm)
            ENDIF
            tr(ji,jj,jk,jpdic,Krhs) = tr(ji,jj,jk,jpdic,Krhs) + zdissol
            tr(ji,jj,jk,jptal,Krhs) = tr(ji,jj,jk,jptal,Krhs) + 2.0 * zdissol
            !
            IF( l_dia ) zcaldiss(ji,jj,jk) = zdissol
            !
         ENDIF
      END_3D
      !$OMP END PARALLEL DO
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj) &
      !$OMP SHARED(sinkcalb, wsbio4, mbkt, zcaco3, rfact2, rday)
      DO_2D( 0, 0, 0, 0 )
         sinkcalb(ji,jj) = wsbio4(ji,jj,mbkt(ji,jj)) &
            &            * zcaco3(ji,jj,mbkt(ji,jj)) * rfact2 / rday
      END_2D
      !$OMP END PARALLEL DO
      !
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('lys ')")
        CALL prt_ctl_info( charout, cdcomp = 'top' )
 !       CALL prt_ctl(tab4d_1=tr(:,:,:,:,Krhs), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( l_dia .AND. knt == nrdttrc ) THEN
         ALLOCATE( zw3d(GLOBAL_2D_ARRAY,1:jpk) )   ;   zw3d(:,:,:) = 0.
         zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jkR, zfact) &
         !$OMP SHARED(zw3d, prodcal, tmask)
         DO_3D( 0, 0, 0, 0, 1, jpk)
            zw3d(ji,jj,jkR) = prodcal(ji,jj,jk) * zfact * tmask(ji,jj,jk)
         END_3D
         !$OMP END PARALLEL DO
         CALL iom_put( "PCAL", zw3d )   ! calcite production
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jkR, zfact) &
         !$OMP SHARED(zw3d, zcaldiss, tmask)
         DO_3D( 0, 0, 0, 0, 1, jpk)
            zw3d(ji,jj,jkR) = zcaldiss(ji,jj,jk) * zfact * tmask(ji,jj,jk)
         END_3D
         !$OMP END PARALLEL DO
         CALL iom_put( "DCAL", zw3d )  ! calcite dissolution
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jkR) &
         !$OMP SHARED(zw3d, wsbio4, zcaco3, rday, tmask)
         DO_3D( 0, 0, 0, 0, 1, jpk)
            zw3d(ji,jj,jkR) = wsbio4(ji,jj,jk) * zcaco3(ji,jj,jk) &
              &       * 1.e+3 / rday  * tmask(ji,jj,jk)
         END_3D
         !$OMP END PARALLEL DO
         !
         CALL iom_put( "EPCAL100", zw3d(:,:,ik100) ) ! Export of calcite at 100m
         CALL iom_put( "EXPCAL" , zw3d ) ! Export of calcite in the water column
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jkR) &
         !$OMP SHARED(zw3d, hi, rtrn, tmask)
         DO_3D( 0, 0, 0, 0, 1, jpk)
            zw3d(ji,jj,jkR) = -1. * LOG10( MAX( hi(ji,jj,jk),rtrn ) ) & 
                &                  * tmask(ji,jj,jk)         
         END_3D
         !$OMP END PARALLEL DO
         CALL iom_put( "PH" , zw3d )  ! PH
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jkR) &
         !$OMP SHARED(zw3d, zco3, tmask)
         DO_3D( 0, 0, 0, 0, 1, jpk)
            zw3d(ji,jj,jkR) = zco3(ji,jj,jk) * 1.e+3 * tmask(ji,jj,jk)
         END_3D
         !$OMP END PARALLEL DO
         CALL iom_put( "CO3", zw3d ) ! ion carbonate
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jkR, zcalcon, zfact) &
         !$OMP SHARED(zw3d, aksp, rhop, calcon, salinprac, tmask, rtrn)
         DO_3D( 0, 0, 0, 0, 1, jpk)
            zcalcon         = calcon * ( salinprac(ji,jj,jk) / 35._wp )
            zfact           = rhop(ji,jj,jk) / 1000._wp
            zw3d(ji,jj,jkR) = aksp(ji,jj,jk) * zfact &
              &             / ( zcalcon + rtrn ) * 1.e+3 * tmask(ji,jj,jk)
         END_3D
         !$OMP END PARALLEL DO
         !
         CALL iom_put( "CO3sat", zw3d )  ! calcite saturation
         !
         DEALLOCATE( zcaldiss, zw3d )
      ENDIF
# if defined key_trc_diaadd
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj, jkR, zcalcon, zfact) &
      !$OMP SHARED(trc3d, hi, tmask, zco3, aksp, rhop, calcon, salinprac, rtrn)
      DO_3D( 0, 0, 0, 0, 1, jpk)
         trc3d(ji,jj,jkR,jp_hi    ) = -1. * LOG10( MAX( hi(ji,jj,jk), + 1.0e-12, rtrn) ) * tmask(ji,jj,jk) ! PH
         trc3d(ji,jj,jkR,jp_co3   ) = zco3(ji,jj,jk)  * 1.e+3 * tmask(ji,jj,jk) ! Ion carbonate
         trc3d(ji,jj,jkR,jp_co3sat) = aksp(ji,jj,jk) * rhop(ji,jj,jk) / 1000._wp &
              &   / calcon * ( salinprac(ji,jj,jk) / 35._wp ) * tmask(ji,jj,jk) ! calcite saturation
      END_3D
      !$OMP END PARALLEL DO
# endif
      !
      IF( ln_timing )   CALL timing_stop('p2z_lys')
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
      INTEGER  ::   ji, jj, jk, jn
      REAL(wp) ::   zdispot, zfact, zcalcon, ztra, zdissol
      REAL(wp) ::   zomegaca, zexcess, zexcess0, zkd
      CHARACTER (len=25) ::   charout
      REAL(wp), DIMENSION(A2D(0),jpk) :: zhinit, zhi, zco3
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:)  :: zw3d, zcaldiss
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('p4z_lys')
      !
     IF( kt == nittrc000 ) &
           & l_dia = iom_use( "PH" ) .OR. iom_use( "CO3" ) .OR. iom_use( "CO3sat" ) &
           &                  .OR. iom_use( "DCAL" ) .OR. iom_use( "PCAL" )

      IF( l_dia )   THEN
         ALLOCATE( zcaldiss(A2D(0),jpk) )
         zcaldiss(A2D(0),:) = 0._wp
         zco3(A2D(0),jpk) = 0._wp
      ENDIF
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj, jk) &
      !$OMP SHARED(zhinit, hi, rhop, rtrn)
      DO_3D( 0, 0, 0, 0, 1, jpkm1)
         zhinit(ji,jj,jk) = hi(ji,jj,jk) * 1000._wp / ( rhop(ji,jj,jk) + rtrn )
      END_3D
      !$OMP END PARALLEL DO
      !
      !     -------------------------------------------
      !     COMPUTE [CO3--] and [H+] CONCENTRATIONS
      !     -------------------------------------------

      CALL solve_at_general( zhinit, zhi, Kbb )
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj, jk) &
      !$OMP SHARED(zco3, tr, jpdic, Kbb, ak13, ak23, zhi, rtrn, hi, rhop)
      DO_3D( 0, 0, 0, 0, 1, jpkm1)
         zco3(ji,jj,jk) = tr(ji,jj,jk,jpdic,Kbb) * ak13(ji,jj,jk) * ak23(ji,jj,jk) &
                          / (zhi(ji,jj,jk)**2 + ak13(ji,jj,jk) * zhi(ji,jj,jk) &
                          + ak13(ji,jj,jk) * ak23(ji,jj,jk) + rtrn )
         hi  (ji,jj,jk) = zhi(ji,jj,jk) * rhop(ji,jj,jk) / 1000._wp
      END_3D
      !$OMP END PARALLEL DO
      !     ---------------------------------------------------------
      !        CALCULATE DEGREE OF CACO3 SATURATION AND CORRESPONDING
      !        DISSOLOUTION AND PRECIPITATION OF CACO3 (BE AWARE OF
      !        MGCO3)
      !     ---------------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj, jk, zcalcon, zfact, zomegaca, excess, &
                           zexcess0, zexcess, zdispot, zkd, ztra) &
      !$OMP SHARED(tr, jptal, jpcal, jpdic, Krhs, zco3, aksp, rtrn, &
                           rhop, calcon, salinprac, nca, kdca, rfact2, &
                           rmtss, l_dia, zcaldiss)
      DO_3D( 0, 0, 0, 0, 1, jpkm1)

         ! DEVIATION OF [CO3--] FROM SATURATION VALUE
         ! Salinity dependance in zomegaca and divide by rhd to have good units
         zcalcon  = calcon * ( salinprac(ji,jj,jk) / 35._wp )
         zfact    = rhop(ji,jj,jk) / 1000._wp
         zomegaca = ( zcalcon * zco3(ji,jj,jk) ) / ( aksp(ji,jj,jk) * zfact + rtrn )

         ! SET DEGREE OF UNDER-/SUPERSATURATION
         excess(ji,jj,jk) = 1._wp - zomegaca
         zexcess0 = MAX( 0., excess(ji,jj,jk) )

         IF( zomegaca < 0.8 ) THEN
            zexcess = zexcess0**nca
            ! AMOUNT CACO3 THAT RE-ENTERS SOLUTION
            zdispot = kdca * zexcess * tr(ji,jj,jk,jpcal,Kbb)
         ELSE
            zkd = kdca * 0.2**(nca - 0.2)
            zexcess = zexcess0**0.2
            zdispot = zkd * zexcess * tr(ji,jj,jk,jpcal,Kbb)
        ENDIF

        !  CHANGE OF [CO3--] , [ALK], PARTICULATE [CACO3],
        !       AND [SUM(CO2)] DUE TO CACO3 DISSOLUTION/PRECIPITATION
        ztra  = zdispot * rfact2 / rmtss ! calcite dissolution
        !
        tr(ji,jj,jk,jptal,Krhs) = tr(ji,jj,jk,jptal,Krhs) + 2. * ztra
        tr(ji,jj,jk,jpcal,Krhs) = tr(ji,jj,jk,jpcal,Krhs) -      ztra
        tr(ji,jj,jk,jpdic,Krhs) = tr(ji,jj,jk,jpdic,Krhs) +      ztra
        !
        IF( l_dia ) zcaldiss(ji,jj,jk) = zdissol
        !
      END_3D
      !$OMP END PARALLEL DO
      !
      IF( l_dia .AND. knt == nrdttrc ) THEN
         ALLOCATE( zw3d(GLOBAL_2D_ARRAY,1:jpk) )  ;  zw3d(:,:,:) = 0._wp
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jk, zfact) &
         !$OMP SHARED(zw3d, prodcal, tmask, jpk)
         zfact = 1.e+3 * rfact2r  ! conversion from mol/l/kt to  mol/m3/s
         DO_3D( 0, 0, 0, 0, 1, jpk)
            zw3d(ji,jj,jkR) = prodcal(ji,jj,jk) * zfact * tmask(ji,jj,jk)
         END_3D
         !$OMP END PARALLEL DO
         CALL iom_put( "PCAL", zw3d ) ! calcite production
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jk, zfact) &
         !$OMP SHARED(zw3d, zcaldiss, tmask, jpk)
         DO_3D( 0, 0, 0, 0, 1, jpk)
            zw3d(ji,jj,jkR) = zcaldiss(ji,jj,jk) * zfact * tmask(ji,jj,jk)
         END_3D
         !$OMP END PARALLEL DO
         CALL iom_put( "DCAL", zw3d ) ! calcite dissolution
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jk) &
         !$OMP SHARED(zw3d, hi, rtrn, tmask, jpk)
         DO_3D( 0, 0, 0, 0, 1, jpk)
            zw3d(ji,jj,jkR) = -1. * LOG10( MAX( hi(ji,jj,jk),rtrn ) ) &
                &                 * tmask(ji,jj,jk)
         END_3D
         !$OMP END PARALLEL DO
         CALL iom_put( "PH" , zw3d ) ! PH
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jk) &
         !$OMP SHARED(zw3d, zco3, tmask, jpk)
         DO_3D( 0, 0, 0, 0, 1, jpk)
            zw3d(ji,jj,jkR) = zco3(ji,jj,jk) * 1.e+3 * tmask(ji,jj,jk)
         END_3D
         !$OMP END PARALLEL DO
         CALL iom_put( "CO3", zw3d ) ! ion carbonate
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jk, zcalcon, zfact) &
         !$OMP SHARED(zw3d, aksp, rhop, salinprac, calcon, rtrn, tmask, jpk)
         DO_3D( 0, 0, 0, 0, 1, jpk)
             zcalcon        = calcon * ( salinprac(ji,jj,jk) / 35._wp )
             zfact          = rhop(ji,jj,jk) / 1000._wp
             zw3d(ji,jj,jkR) = aksp(ji,jj,jk) * zfact &
                 &        / ( zcalcon + rtrn )  * 1.e+3 * tmask(ji,jj,jk)
         END_3D
         !$OMP END PARALLEL DO
         CALL iom_put( "CO3sat", zw3d )  ! calcite saturation
         !
         DEALLOCATE( zcaldiss, zw3d )
      ENDIF
      !
# if defined key_trc_diaadd
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj, jk) &
      !$OMP SHARED(trc3d, hi, zco3, aksp, rhop, &
                   calcon, salinprac, rtrn, tmask, jpk)
      DO_3D( 0, 0, 0, 0, 1, jpk)
         ! PH
         trc3d(ji,jj,jkR,jp_hi) = -1. * LOG10( MAX( hi(ji,jj,jk), + 1.0e-12, rtrn) ) &
                                      * tmask(ji,jj,jk)
         ! Ion carbonate
         trc3d(ji,jj,jkR,jp_co3) = zco3(ji,jj,jk) * 1.e+3 &
                                      * tmask(ji,jj,jk)
         ! calcite saturation
         trc3d(ji,jj,jkR,jp_co3sat) = aksp(ji,jj,jk) &
                                      * rhop(ji,jj,jk) / 1000._wp / calcon &
                                      * ( salinprac(ji,jj,jk) / 35._wp ) &
                                      * tmask(ji,jj,jk)
      END_3D
      !$OMP END PARALLEL DO
# endif
      
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('lys ')")
        CALL prt_ctl_info( charout, cdcomp = 'top' )
!        CALL prt_ctl(tab4d_1=tr(:,:,:,:,Krhs), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_lys')
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
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_lys_init : initialization of CaCO3 dissolution'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      READ_NML_REF(numnatp,nampiscal)
      READ_NML_CFG(numnatp,nampiscal)
      IF(lwm) WRITE( numonp, nampiscal )
      !
      IF(lwp) THEN          ! control print
         WRITE(numout,*) '   Namelist : nampiscal'
         WRITE(numout,*) '   diss. rate constant calcite (per month)        kdca =', kdca
         WRITE(numout,*) '   order of reaction for calcite dissolution      nca  =', nca
      ENDIF
      !
      ! Number of seconds per month 
      rmtss =  nyear_len(1) * rday / raamo
      !
      ! CE not really needed ; temporary, shoud be removed when quotan( A2D(0),jpk )
      excess(:,:,:) = 0._wp
      !
   END SUBROUTINE p4z_lys_init

#else
   !!======================================================================
   !!  Dummy module :                                   No PISCES bio-model
   !!======================================================================
CONTAINS
   SUBROUTINE p4z_lys( kt )                   ! Empty routine
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'p4z_lys: You should not have seen this print! error?', kt
   END SUBROUTINE p4z_lys
#endif

   !!======================================================================
END MODULE p4zlys
