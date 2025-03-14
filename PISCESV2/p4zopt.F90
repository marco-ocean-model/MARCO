#include "cppdefs.h"

MODULE p4zopt
   !!======================================================================
   !!                         ***  MODULE p4zopt  ***
   !! TOP - PISCES : Compute the light availability in the water column
   !!======================================================================
   !! History :  1.0  !  2004     (O. Aumont) Original code
   !!            2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!            3.2  !  2009-04  (C. Ethe, G. Madec)  optimisation
   !!            3.4  !  2011-06  (O. Aumont, C. Ethe) Improve light availability of nano & diat
   !!            3.*  !  2025-02  (S. Maishal, R. Person) MPI and optimization
#if defined  key_pisces   
   !!----------------------------------------------------------------------
   !!   p4z_opt       : light availability in the water column
   !!----------------------------------------------------------------------
   USE trc            ! tracer variables
   USE oce_trc        ! tracer-ocean share variables
   USE sms_pisces     ! Source Minus Sink of PISCES
   USE iom            ! I/O manager
   USE prtctl         ! print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_opt        ! called in p4zbio.F90 module
   PUBLIC   p4z_opt_init   ! called in trcsms_pisces.F90 module
   PUBLIC   p4z_opt_alloc

   !! * Shared module variables

   LOGICAL  ::   ln_varpar   ! boolean for variable PAR fraction
   REAL(wp) ::   parlux      ! Fraction of shortwave as PAR
   REAL(wp) ::   xparsw      ! parlux/3
   REAL(wp) ::   xsi0r       ! 1. /rn_si0

!# ifdef NEMO
!   TYPE(FLD), ALLOCATABLE, DIMENSION(:) :: sf_par ! structure of input par
!# endif
   INTEGER , PARAMETER :: nbtimes = 366  !: maximum number of times record in a file
   INTEGER  :: ntimes_par                ! number of time steps in a file
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   par_varsw      ! PAR fraction of shortwave
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ekb, ekg, ekr  ! wavelength (Red-Green-Blue)
 
   LOGICAL  :: l_dia_heup, l_dia_par 

      !! * Substitutions
#  include "ocean2pisces.h90"   
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
#  include "domzgr_substitute.h90"

      !!----------------------------------------------------------------------
      !! NEMO/TOP 4.0 , NEMO Consortium (2018)
      !! $Id: p4zopt.F90 15459 2021-10-29 08:19:18Z cetlod $ 
      !! Software governed by the CeCILL license (see ./LICENSE)
      !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_opt( kt, knt, Kbb, Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_opt  ***
      !!
      !! ** Purpose :   Compute the light availability in the water column
      !!              depending on the depth and the chlorophyll concentration
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm  ! time level indices
      !
      INTEGER  ::   ji, jj, jk
      INTEGER  ::   irgb
      REAL(wp) ::   zchl
      REAL(wp) ::   zc0 , zc1 , zc2, zc3, z1_dep
      REAL(wp), DIMENSION(A2D(0)    ) :: zdepmoy, zetmp1
      REAL(wp), DIMENSION(A2D(0)    ) :: zqsr100, zqsr_corr
      REAL(wp), DIMENSION(A2D(0),jpk) :: ze0, ze1, ze2, ze3
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zpar
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zw2d
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_opt')

      IF( kt == nittrc000 )  THEN
         l_dia_heup = iom_use( "Heup") 
         l_dia_par  = iom_use( "PAR" ) 
      ENDIF

      IF( knt == 1 .AND. ln_varpar )   CALL p4z_opt_sbc( kt )

      !
      ! Attenuation coef. function of Chlorophyll and wavelength (Red-Green-Blue)
      ! Thus the light penetration scheme is based on a decomposition of PAR
      ! into three wave length domains. This was first officially published
      ! in Lengaigne et al. (2007).
      ! --------------------------------------------------------
      !
      ! Computation of the light attenuation parameters based on a 
      ! look-up table
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj, jk, irgb, zchl) &
      !$OMP SHARED(tr, thetanano, rtrn, ln_p2z, ln_p5z, rkrgb, e3t)
      DO_3D( 0, 0, 0, 0, 1, jpkm1)
         IF( ln_p2z ) THEN
            zchl = ( tr(ji,jj,jk,jpphy,Kbb) &
              &     * 12.0 * thetanano(ji,jj,jk) + rtrn ) * 1.e6
         ELSE
            zchl =  ( tr(ji,jj,jk,jpnch,Kbb) &
              &     + tr(ji,jj,jk,jpdch,Kbb) + rtrn ) * 1.e6
            IF( ln_p5z ) zchl = zchl + tr(ji,jj,jk,jppch,Kbb) * 1.e6
         ENDIF
         zchl = MIN( 10. , MAX( 0.05, zchl ) )
         irgb = NINT( 41 + 20.* LOG10( zchl ) + rtrn )
         !
         ekb(ji,jj,jk) = rkrgb(1,irgb) * e3t(ji,jj,jk,Kmm)
         ekg(ji,jj,jk) = rkrgb(2,irgb) * e3t(ji,jj,jk,Kmm)
         ekr(ji,jj,jk) = rkrgb(3,irgb) * e3t(ji,jj,jk,Kmm)
      END_3D
      !$OMP END PARALLEL DO
      ! Photosynthetically Available Radiation (PAR)
      ! Two cases are considered in the following : 
      ! (1) An explicit diunal cycle is activated. In that case, mean 
      ! QSR is used as PISCES in its current state has not been parameterized
      ! for an explicit diurnal cycle
      ! (2) no diurnal cycle of SW is active and in that case, QSR is used.
      ! --------------------------------------------
      IF( ln_trcdc2dm ) THEN                     !  diurnal cycle
         IF ( ln_p4z_dcyc ) THEN   ! Diurnal cycle in PISCES
            !
            !
            ! SW over the ice free zone of the grid cell. This assumes that
            ! SW is zero below sea ice which is a very crude assumption that is 
            ! not fully correct with LIM3 and SI3 but no information is 
            ! currently available to do a better job. SHould be improved in the 
            ! (near) future.
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(ji, jj) &
            !$OMP SHARED(qsr_mean, fr_i, rtrn, zqsr_corr)
            DO_2D( 0, 0, 0, 0 )
               zqsr_corr(ji,jj) = qsr_mean(ji,jj) / ( 1.-fr_i(ji,jj) + rtrn )
            END_2D
            !$OMP END PARALLEL DO
            !
            CALL p4z_opt_par( kt, Kmm, zqsr_corr, ze1, ze2, ze3, pqsr100 = zqsr100 )
            !
            ! Used PAR is computed for each phytoplankton species
            ! etot_ndcy is PAR at level jk averaged over 24h.
            ! Due to their size, they have different light absorption characteristics
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(ji, jj, jk) &
            !$OMP SHARED(etot_ndcy, ze1, ze2, ze3)
            DO_3D( 0, 0, 0, 0, 1, nksr )
               etot_ndcy(ji,jj,jk) = ze1(ji,jj,jk) + ze2(ji,jj,jk) + ze3(ji,jj,jk)
            END_3D
            !$OMP END PARALLEL DO
            !
            ! SW over the ice free zone of the grid cell. This assumes that
            ! SW is zero below sea ice which is a very crude assumption that is 
            ! not fully correct with LIM3 and SI3 but no information is 
            ! currently available to do a better job. SHould be improved in the 
            ! (near) future.
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(ji, jj) &
            !$OMP SHARED(qsr, fr_i, rtrn, zqsr_corr)
            DO_2D( 0, 0, 0, 0 )
               zqsr_corr(ji,jj) = qsr(ji,jj) / ( 1.-fr_i(ji,jj) + rtrn )
            END_2D
            !$OMP END PARALLEL DO
            !
            CALL p4z_opt_par( kt, Kmm, zqsr_corr, ze1, ze2, ze3 )
            !
            ! Total PAR computation at level jk that includes the diurnal cycle
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(ji, jj, jk) &
            !$OMP SHARED(etot, enano, ze1, ze2, ze3)
            DO_3D( 0, 0, 0, 0, 1, nksr )
               etot (ji,jj,jk) = ze1(ji,jj,jk) + ze2(ji,jj,jk) + ze3(ji,jj,jk)
               enano(ji,jj,jk) = 1.85 * ze1(ji,jj,jk) + 0.69 * ze2(ji,jj,jk) &
                               + 0.46 * ze3(ji,jj,jk)
            END_3D
            !$OMP END PARALLEL DO
            IF( .NOT. ln_p2z ) THEN
               !$OMP PARALLEL DO &
               !$OMP PRIVATE(ji, jj, jk) &
               !$OMP SHARED(ediat, ze1, ze2, ze3)
               DO_3D( 0, 0, 0, 0, 1, nksr )
                  ediat(ji,jj,jk) = 1.62 * ze1(ji,jj,jk) + 0.74 * ze2(ji,jj,jk) &
                                  + 0.63 * ze3(ji,jj,jk)
               END_3D
               !$OMP END PARALLEL DO
               IF( ln_p5z ) THEN
                  !$OMP PARALLEL DO &
                  !$OMP PRIVATE(ji, jj, jk) &
                  !$OMP SHARED(epico, ze1, ze2, ze3)
                  DO_3D( 0, 0, 0, 0, 1, nksr )
                     epico(ji,jj,jk) = 1.94 * ze1(ji,jj,jk) + 0.66 * ze2(ji,jj,jk) &
                                     + 0.4 * ze3(ji,jj,jk)
                  END_3D
                  !$OMP END PARALLEL DO
               ENDIF
            ENDIF

         ELSE ! No diurnal cycle in PISCES

            !
            !
            ! SW over the ice free zone of the grid cell. This assumes that
            ! SW is zero below sea ice which is a very crude assumption that is 
            ! not fully correct with LIM3 and SI3 but no information is 
            ! currently available to do a better job. SHould be improved in the 
            ! (near) future.
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(ji, jj) &
            !$OMP SHARED(qsr_mean, fr_i, rtrn, zqsr_corr)
            DO_2D( 0, 0, 0, 0 )
               zqsr_corr(ji,jj) = qsr_mean(ji,jj) / ( 1.-fr_i(ji,jj) + rtrn )
            END_2D
            !$OMP END PARALLEL DO
            !
            CALL p4z_opt_par( kt, Kmm, zqsr_corr, ze1, ze2, ze3, pqsr100 = zqsr100 )
            !
            ! Used PAR is computed for each phytoplankton species
            ! etot_ndcy is PAR at level jk averaged over 24h.
            ! Due to their size, they have different light absorption characteristics
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(ji, jj, jk) &
            !$OMP SHARED(etot_ndcy, enano, ze1, ze2, ze3)
            DO_3D( 0, 0, 0, 0, 1, nksr )
               etot_ndcy(ji,jj,jk) = ze1(ji,jj,jk) + ze2(ji,jj,jk) + ze3(ji,jj,jk)
               enano    (ji,jj,jk) = 1.85 * ze1(ji,jj,jk) + 0.69 * ze2(ji,jj,jk) &
                                   + 0.46 * ze3(ji,jj,jk)
            END_3D
            !$OMP END PARALLEL DO
            IF( .NOT. ln_p2z ) THEN
               !$OMP PARALLEL DO &
               !$OMP PRIVATE(ji, jj, jk) &
               !$OMP SHARED(ediat, ze1, ze2, ze3)
               DO_3D( 0, 0, 0, 0, 1, nksr )
                  ediat(ji,jj,jk) = 1.62 * ze1(ji,jj,jk) + 0.74 * ze2(ji,jj,jk) &
                                  + 0.63 * ze3(ji,jj,jk)
               END_3D
               !$OMP END PARALLEL DO
               IF( ln_p5z ) THEN
                  !$OMP PARALLEL DO &
                  !$OMP PRIVATE(ji, jj, jk) &
                  !$OMP SHARED(epico, ze1, ze2, ze3)
                  DO_3D( 0, 0, 0, 0, 1, nksr )
                     epico(ji,jj,jk) = 1.94 * ze1(ji,jj,jk) + 0.66 * ze2(ji,jj,jk) &
                                     + 0.4 * ze3(ji,jj,jk)
                  END_3D
                  !$OMP END PARALLEL DO
               ENDIF
            ENDIF
            !
            ! SW over the ice free zone of the grid cell. This assumes that
            ! SW is zero below sea ice which is a very crude assumption that is -*-
            ! not fully correct with LIM3 and SI3 but no information is 
            ! currently available to do a better job. Should be improved in the 
            ! (near) future.-*-
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(ji, jj) &
            !$OMP SHARED(qsr, fr_i, rtrn, zqsr_corr)
            DO_2D( 0, 0, 0, 0 )
               zqsr_corr(ji,jj) = qsr(ji,jj) / ( 1.-fr_i(ji,jj) + rtrn )
            END_2D
            !$OMP END PARALLEL DO
            !
            CALL p4z_opt_par( kt, Kmm, zqsr_corr, ze1, ze2, ze3 ) 
            !
            ! Total PAR computation at level jk that includes the diurnal cycle
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(ji, jj, jk) &
            !$OMP SHARED(etot, ze1, ze2, ze3)
            DO_3D( 0, 0, 0, 0, 1, nksr )
               etot(ji,jj,jk) = ze1(ji,jj,jk) + ze2(ji,jj,jk) + ze3(ji,jj,jk)
            END_3D
            !$OMP END PARALLEL DO
         ENDIF
         !
      ELSE   ! no diurnal cycle
         !
         !
         ! SW over the ice free zone of the grid cell. This assumes that
         ! SW is zero below sea ice which is a very crude assumption that is
         ! not fully correct with LIM3 and SI3 but no information is 
         ! currently available to do a better job. Should be improved in the -*-
         ! (near) future._*_
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj) &
         !$OMP SHARED(qsr, fr_i, rtrn, zqsr_corr)
         DO_2D( 0, 0, 0, 0 )
            zqsr_corr(ji,jj) = qsr(ji,jj) / ( 1.-fr_i(ji,jj) + rtrn )
         END_2D
         !$OMP END PARALLEL DO
         !
         CALL p4z_opt_par( kt, Kmm, zqsr_corr, ze1, ze2, ze3, pqsr100 = zqsr100 ) 
         !

         ! Used PAR is computed for each phytoplankton species
         ! Due to their size, they have different light absorption characteristics
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jk) &
         !$OMP SHARED(etot, enano, ze1, ze2, ze3)
         DO_3D( 0, 0, 0, 0, 1, nksr )
            etot (ji,jj,jk) = ze1(ji,jj,jk) + ze2(ji,jj,jk) + ze3(ji,jj,jk)
            enano(ji,jj,jk) = 1.85 * ze1(ji,jj,jk) + 0.69 * ze2(ji,jj,jk) &
                            + 0.46 * ze3(ji,jj,jk)
         END_3D
         !$OMP END PARALLEL DO
         IF( .NOT. ln_p2z ) THEN
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(ji, jj, jk) &
            !$OMP SHARED(ediat, ze1, ze2, ze3)
            DO_3D( 0, 0, 0, 0, 1, nksr )
               ediat(ji,jj,jk) = 1.62 * ze1(ji,jj,jk) + 0.74 &
                               * ze2(ji,jj,jk) + 0.63 * ze3(ji,jj,jk)
            END_3D
            !$OMP END PARALLEL DO
            IF( ln_p5z ) THEN
               !$OMP PARALLEL DO &
               !$OMP PRIVATE(ji, jj, jk) &
               !$OMP SHARED(epico, ze1, ze2, ze3)
               DO_3D( 0, 0, 0, 0, 1, nksr )
                  epico(ji,jj,jk) = 1.94 * ze1(ji,jj,jk) + 0.66 &
                                  * ze2(ji,jj,jk) + 0.4 * ze3(ji,jj,jk)
               END_3D
               !$OMP END PARALLEL DO
            ENDIF
         ENDIF
         etot_ndcy(:,:,:) =  etot(:,:,:) 
      ENDIF


      ! Biophysical feedback part (computation of vertical penetration of SW)
      IF( ln_qsr_bio ) THEN!* heat flux accros w-level (used in the dynamics)
         !------------------------------------------ ------------------------
         DO_2D( 0, 0, 0, 0 )
            zqsr_corr(ji,jj) = qsr(ji,jj) 
         END_2D
         !
         CALL p4z_opt_par( kt, Kmm, zqsr_corr, ze1, ze2, ze3, pe0=ze0 )
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj) &
         !$OMP SHARED(etot3, qsr, tmask)
         DO_2D( 0, 0, 0, 0 )
            etot3(ji,jj,1) = qsr(ji,jj) * tmask(ji,jj,1)
         END_2D
         !$OMP END PARALLEL DO
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jk) &
         !$OMP SHARED(etot3, ze0, ze1, ze2, ze3, tmask)
         DO_3D( 0, 0, 0, 0, 2, nksr+1 )
            etot3(ji,jj,jk) = ( ze0(ji,jj,jk) + ze1(ji,jj,jk) + ze2(ji,jj,jk) &
                            + ze3(ji,jj,jk) ) * tmask(ji,jj,jk)
         END_3D
         !$OMP END PARALLEL DO
         !-------------------------------------------------------------------
      ENDIF
      
      ! Euphotic depth and level
      ! Two definitions of the euphotic zone are used here. 
      ! (1) The classical definition based on the relative threshold value
      ! (2) An alternative definition based on a absolute threshold value.
      ! ---------------------------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj) &
      !$OMP SHARED(neln, heup, heup_01, gdepw)
      DO_2D( 0, 0, 0, 0 )
         neln   (ji,jj) = 1
         heup   (ji,jj) = gdepw(ji,jj,2,Kmm)
         heup_01(ji,jj) = gdepw(ji,jj,2,Kmm)
      END_2D
      !$OMP END PARALLEL DO
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj, jk) &
      !$OMP SHARED(etot_ndcy, tmask, zqsr100, neln, heup, heup_01, gdepw, Kmm)
      DO_3D( 0, 0, 0, 0, 2, nksr)
        IF( etot_ndcy(ji,jj,jk) * tmask(ji,jj,jk) >=  zqsr100(ji,jj) )  THEN
           neln(ji,jj) = jk+1                  ! Euphotic level : 1rst T-level strictly below Euphotic layer
           !                                   ! nb: ensure the compatibility with nmld_trc definition in trd_mld_trc_zint
           heup(ji,jj) = gdepw(ji,jj,jk+1,Kmm) ! Euphotic layer depth
        ENDIF
        IF( etot_ndcy(ji,jj,jk) * tmask(ji,jj,jk) >= 0.10 ) THEN
           heup_01(ji,jj) = gdepw(ji,jj,jk+1,Kmm) ! Euphotic layer depth (light level definition)
        ENDIF
      END_3D
      !$OMP END PARALLEL DO
      !
      ! The euphotic depth can not exceed 300 meters.
      heup   (:,:) = MIN( 300., heup(:,:) )
      heup_01(:,:) = MIN( 300., heup_01(:,:) )
      ! Mean PAR over the mixed layer
      ! -----------------------------
      zdepmoy(:,:) = 0.e0
      zetmp1 (:,:) = 0.e0
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj, jk) &
      !$OMP SHARED(gdepw, hmld, zetmp1, zdepmoy, etot, e3t, Kmm)
      DO_3D( 0, 0, 0, 0, 1, nksr)
         IF( gdepw(ji,jj,jk+1,Kmm) <= hmld(ji,jj) ) THEN
            zetmp1 (ji,jj) = zetmp1 (ji,jj) + etot(ji,jj,jk) &
                           * e3t(ji,jj,jk,Kmm) ! Actual PAR for remineralisation
            zdepmoy(ji,jj) = zdepmoy(ji,jj) + e3t(ji,jj,jk,Kmm)
         ENDIF
      END_3D
      !$OMP END PARALLEL DO
      !
      emoy(:,:,:) = etot(:,:,:) ! remineralisation
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj, jk, z1_dep) &
      !$OMP SHARED(gdepw, hmld, zetmp1, zdepmoy, rtrn, emoy, Kmm)
      DO_3D( 0, 0, 0, 0, 1, nksr)
         IF( gdepw(ji,jj,jk+1,Kmm) <= hmld(ji,jj) ) THEN
            z1_dep = 1. / ( zdepmoy(ji,jj) + rtrn )
            emoy (ji,jj,jk) = zetmp1(ji,jj) * z1_dep
         ENDIF
      END_3D
      !$OMP END PARALLEL DO
      ! Computation of the mean usable light for the different phytoplankton
      ! groups based on their absorption characteristics.
      zetmp1(:,:)   = 0.e0
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj, jk) &
      !$OMP SHARED(gdepw, hmld, heup_01, zetmp1, enano, e3t, Kmm)
      DO_3D( 0, 0, 0, 0, 1, nksr)
         IF( gdepw(ji,jj,jk+1,Kmm) <= MIN(hmld(ji,jj), heup_01(ji,jj)) ) THEN
            zetmp1(ji,jj) = zetmp1(ji,jj) + enano(ji,jj,jk) &
                          * e3t(ji,jj,jk,Kmm) ! Nanophytoplankton
         ENDIF
      END_3D
      !$OMP END PARALLEL DO
      enanom(:,:,:) = enano(:,:,:)
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(ji, jj, jk, z1_dep) &
      !$OMP SHARED(gdepw, hmld, zetmp1, zdepmoy, rtrn, enanom, Kmm)
      DO_3D( 0, 0, 0, 0, 1, nksr)
         IF( gdepw(ji,jj,jk+1,Kmm) <= hmld(ji,jj) ) THEN
            z1_dep = 1. / ( zdepmoy(ji,jj) + rtrn )
            enanom(ji,jj,jk) = zetmp1(ji,jj) * z1_dep
         ENDIF
      END_3D
      !$OMP END PARALLEL DO
      !
      IF( .NOT. ln_p2z ) THEN
         ! Diatoms when using PISCES-operational or PISCES-QUOTA
         zetmp1(:,:) = 0.e0
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jk) &
         !$OMP SHARED(gdepw, hmld, heup_01, zetmp1, ediat, e3t, Kmm)
         DO_3D( 0, 0, 0, 0, 1, nksr)
            IF( gdepw(ji,jj,jk+1,Kmm) <= MIN(hmld(ji,jj), heup_01(ji,jj)) ) THEN
               zetmp1(ji,jj) = zetmp1(ji,jj) + ediat(ji,jj,jk) * e3t(ji,jj,jk,Kmm) ! Diatoms
            ENDIF
         END_3D
         !$OMP END PARALLEL DO
         !
         ediatm(:,:,:) = ediat(:,:,:)
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jk, z1_dep) &
         !$OMP SHARED(gdepw, hmld, zetmp1, zdepmoy, rtrn, ediatm, Kmm)
         DO_3D( 0, 0, 0, 0, 1, nksr)
            IF( gdepw(ji,jj,jk+1,Kmm) <= hmld(ji,jj) ) THEN
               z1_dep = 1. / ( zdepmoy(ji,jj) + rtrn )
               ediatm(ji,jj,jk) = zetmp1(ji,jj) * z1_dep
            ENDIF
         END_3D
         !$OMP END PARALLEL DO
      ENDIF
      IF( ln_p5z ) THEN
         ! Picophytoplankton when using PISCES-QUOTA
         zetmp1(:,:) = 0.e0
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jk) &
         !$OMP SHARED(gdepw, hmld, heup_01, zetmp1, epico, e3t, Kmm)
         DO_3D( 0, 0, 0, 0, 1, nksr)
            IF( gdepw(ji,jj,jk+1,Kmm) <= MIN(hmld(ji,jj), heup_01(ji,jj)) ) THEN
               zetmp1(ji,jj) = zetmp1(ji,jj) + epico(ji,jj,jk) &
                             * e3t(ji,jj,jk,Kmm)
            ENDIF
         END_3D
         !$OMP END PARALLEL DO
         !
         epicom(:,:,:) = epico(:,:,:)
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(ji, jj, jk, z1_dep) &
         !$OMP SHARED(gdepw, hmld, zetmp1, zdepmoy, rtrn, epicom, Kmm)
         DO_3D( 0, 0, 0, 0, 1, nksr)
            IF( gdepw(ji,jj,jk+1,Kmm) <= hmld(ji,jj) ) THEN
               z1_dep = 1. / ( zdepmoy(ji,jj) + rtrn )
               epicom(ji,jj,jk) = zetmp1(ji,jj) * z1_dep
            ENDIF
         END_3D
         !$OMP END PARALLEL DO
      ENDIF
      !
      IF( l_dia_par .OR. l_diaadd ) THEN
        zetmp1(:,:)   = 0.e0
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(ji, jj, jk) &
        !$OMP SHARED(gdepw, hmld, zetmp1, etot_ndcy, e3t, Kmm)
        DO_3D( 0, 0, 0, 0, 1, nksr)
           IF( gdepw(ji,jj,jk+1,Kmm) <= hmld(ji,jj) ) THEN
              zetmp1(ji,jj) = zetmp1(ji,jj) &
                &       + etot_ndcy(ji,jj,jk) * e3t(ji,jj,jk,Kmm) ! Par averaged over 24h for production
           ENDIF
        END_3D
        !$OMP END PARALLEL DO
      ENDIF
      IF( lk_iomput .AND.  knt == nrdttrc ) THEN
         IF( l_dia_heup ) THEN
           ALLOCATE( zw2d(GLOBAL_2D_ARRAY) )  ;  zw2d(:,:) = 0._wp
           zw2d(A2D(0)) = heup(A2D(0)) * tmask(A2D(0),1)
           CALL iom_put( "Heup", zw2d )  ! Euphotic layer depth
           DEALLOCATE( zw2d ) 
        ENDIF
        IF( l_dia_par ) THEN   ! diagnostic : PAR with no diurnal cycle
           ALLOCATE( zw3d(GLOBAL_2D_ARRAY,jpk) )  ;  zw3d(:,:,:) = 0._wp
           DO_3D( 0, 0, 0, 0, 1, jpk)
              zw3d(ji,jj,jkR) = etot_ndcy(ji,jj,jk) * tmask(ji,jj,jk)
           END_3D
           DO_3D( 0, 0, 0, 0, 1, nksr)
              IF( gdepw(ji,jj,jk+1,Kmm) <= hmld(ji,jj) ) THEN
                 z1_dep = 1. / ( zdepmoy(ji,jj) + rtrn )
                 zw3d(ji,jj,jkR) = zetmp1(ji,jj) * z1_dep
              ENDIF
           END_3D
           DO_2D( 0, 0, 0, 0 )
              zw3d(ji,jj,ktop) = zw3d(ji,jj,ktop) * ( 1._wp - fr_i(ji,jj) )
           END_2D
           CALL iom_put( "PAR", zw3d ) 
           DEALLOCATE( zw3d ) 
        ENDIF
      ENDIF
      !
#if defined key_trc_diaadd
      DO_3D( 0, 0, 0, 0, 1, jpk )
         trc3d(ji,jj,jkR,jp_etot) = etot_ndcy(ji,jj,jk) * tmask(ji,jj,jk) ! PAR
      END_3D
      DO_3D( 0, 0, 0, 0, 1, nksr)
         IF( gdepw(ji,jj,jk+1,Kmm) <= hmld(ji,jj) ) THEN
            z1_dep = 1. / ( zdepmoy(ji,jj) + rtrn )
            trc3d(ji,jj,jkR,jp_etot) = zetmp1(ji,jj) * z1_dep
          ENDIF
      END_3D
      DO_2D( 0, 0, 0, 0 )
         trc3d(ji,jj,ktop,jp_etot) = trc3d(ji,jj,ktop,jp_etot) * ( 1._wp - fr_i(ji,jj) )
      END_2D
      !
      DO_2D( 0, 0, 0, 0 )
         trc2d(ji,jj,jp_heup) = heup(ji,jj) * tmask(ji,jj,1) ! euphotic layer
      END_2D
#endif      
      IF( ln_timing )   CALL timing_stop('p4z_opt')
      !
   END SUBROUTINE p4z_opt


   SUBROUTINE p4z_opt_par( kt, Kmm, pqsr, pe1, pe2, pe3, pe0, pqsr100 )
      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_opt_par  ***
      !!
      !! ** purpose :   compute PAR of each wavelength (Red-Green-Blue)
      !!                for a given shortwave radiation
      !!
      !!----------------------------------------------------------------------
      INTEGER                        , INTENT(in)              ::   kt                ! ocean time-step
      INTEGER                        , INTENT(in)              ::   Kmm               ! ocean time-index
      REAL(wp), DIMENSION(A2D(0))    , INTENT(in   )           ::   pqsr              ! shortwave
      REAL(wp), DIMENSION(A2D(0),jpk), INTENT(inout)           ::   pe1 , pe2 , pe3   ! PAR ( R-G-B)
      REAL(wp), DIMENSION(A2D(0),jpk), INTENT(inout), OPTIONAL ::   pe0               !
      REAL(wp), DIMENSION(A2D(0))    , INTENT(  out), OPTIONAL ::   pqsr100           !
      !
      INTEGER    ::   ji, jj, jk, jkm1     ! dummy loop indices _*_
      REAL(wp), DIMENSION(A2D(0)) ::  zqsr ! shortwave
      !!----------------------------------------------------------------------

      !  Real shortwave
!# ifdef NEMO
!      IF( ln_varpar ) THEN  ;  zqsr(:,:) = par_varsw(:,:) * pqsr(:,:)
!      ELSE                  ;  zqsr(:,:) = xparsw         * pqsr(:,:)
!      ENDIF
!# endif
      zqsr(:,:) = xparsw         * pqsr(:,:) 

      !  Light at the euphotic depth 
      IF( PRESENT( pqsr100 ) )   pqsr100(:,:) = 0.01 * 3. * zqsr(:,:)

      IF( PRESENT( pe0 ) ) THEN     !  W-level
         !
         pe0(:,:,1) = pqsr(:,:) - 3. * zqsr(:,:)    !   ( 1 - 3 * alpha ) * q
         pe1(:,:,1) = zqsr(:,:)         
         pe2(:,:,1) = zqsr(:,:)
         pe3(:,:,1) = zqsr(:,:)
         !
         DO_3D( 0, 0, 0, 0, 2, nksr + 1)
            jkm1 = jk-1
            pe0(ji,jj,jk) = pe0(ji,jj,jk-1) * EXP( -e3t(ji,jj,jkm1,Kmm) * xsi0r )
            pe1(ji,jj,jk) = pe1(ji,jj,jk-1) * EXP( -ekb  (ji,jj,jk-1 )        )
            pe2(ji,jj,jk) = pe2(ji,jj,jk-1) * EXP( -ekg  (ji,jj,jk-1 )        )
            pe3(ji,jj,jk) = pe3(ji,jj,jk-1) * EXP( -ekr  (ji,jj,jk-1 )        )
        END_3D
        !
      ELSE   ! T- level
        !
        pe1(:,:,1) = zqsr(:,:) * EXP( -0.5 * ekb(:,:,1) )
        pe2(:,:,1) = zqsr(:,:) * EXP( -0.5 * ekg(:,:,1) )
        pe3(:,:,1) = zqsr(:,:) * EXP( -0.5 * ekr(:,:,1) )
        !
        DO_3D( 0, 0, 0, 0, 2, nksr)
           pe1(ji,jj,jk) = pe1(ji,jj,jk-1) * EXP( -0.5 * ( ekb(ji,jj,jk-1) + ekb(ji,jj,jk) ) )
           pe2(ji,jj,jk) = pe2(ji,jj,jk-1) * EXP( -0.5 * ( ekg(ji,jj,jk-1) + ekg(ji,jj,jk) ) )
           pe3(ji,jj,jk) = pe3(ji,jj,jk-1) * EXP( -0.5 * ( ekr(ji,jj,jk-1) + ekr(ji,jj,jk) ) )
        END_3D
        !
      ENDIF
      ! 
   END SUBROUTINE p4z_opt_par


   SUBROUTINE p4z_opt_sbc( kt )
      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_opt_sbc  ***
      !!
      !! ** purpose :   read and interpolate the variable PAR fraction
      !!                of shortwave radiation
      !!
      !! ** method  :   read the files and interpolate the appropriate variables
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !
      INTEGER  :: ji,jj
      REAL(wp) :: zcoef
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('p4z_optsbc')
      !
# ifdef NEMO
      ! Compute par_varsw at nit000 or only if there is more than 1 time record in par coefficient file
      IF( ln_varpar ) THEN
         IF( kt == nit000 .OR. ( kt /= nit000 .AND. ntimes_par > 1 ) ) THEN
            CALL fld_read( kt, 1, sf_par )
            DO_2D( 0, 0, 0, 0 )
               par_varsw(ji,jj) = ( sf_par(1)%fnow(ji,jj,1) ) / 3.0
            END_2D
         ENDIF
      ENDIF
# endif 
      !
      IF( ln_timing )  CALL timing_stop('p4z_optsbc')
      !
   END SUBROUTINE p4z_opt_sbc


   SUBROUTINE p4z_opt_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_opt_init  ***
      !!
      !! ** Purpose :   Initialization of tabulated attenuation coef
      !!                and of the percentage of PAR in Shortwave
      !!
      !! ** Input   :   external ascii and netcdf files
      !!----------------------------------------------------------------------
      INTEGER :: numpar, ierr, ios   ! Local integer 
      !
      CHARACTER(len=100) ::  cn_dir          ! Root directory for location of ssr files
# ifdef NEMO
      TYPE(FLD_N) ::   sn_par                ! informations about the fields to be read
      !
      NAMELIST/nampisopt/cn_dir, sn_par, ln_varpar, parlux, ln_p4z_dcyc
# endif
      !!----------------------------------------------------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_opt_init : '
         WRITE(numout,*) '~~~~~~~~~~~~ '
      ENDIF
!# ifdef NEMO
!      READ_NML_REF(numnatp,nampisopt)
!      READ_NML_CFG(numnatp,nampisopt)
!      IF(lwm) WRITE ( numonp, nampisopt )
!# endif
      ln_p4z_dcyc =  .false. ! Diurnal cycle in PISCES
      ln_varpar   = .FALSE.  !
      parlux      =  0.43    ! Fraction of shortwave as PAR
      IF(lwp) THEN
         WRITE(numout,*) '      Namelist : nampisopt '
         WRITE(numout,*) '      PAR as a variable fraction of SW       ln_varpar   = ', ln_varpar
         WRITE(numout,*) '      Default value for the PAR fraction     parlux      = ', parlux
         WRITE(numout,*) '      Activate the diurnal cycle in PISCES   ln_p4z_dcyc = ', ln_p4z_dcyc
      ENDIF
      !
      !
      xparsw = parlux / 3.0
      xsi0r  = 1.e0 / rn_si0

      ! Warning : activate the diurnal cycle with no diurnal cycle in the forcing fields makes no sense
      ! That does not produce a bug because the model does not use the flag but a warning is necessary
      ! ----------------------------------------------------------------------------------------------
      IF ( ln_p4z_dcyc .AND. l_trcdm2dc ) THEN
         IF (lwp) WRITE(numout,*) 'No diurnal cycle in the PAR forcing field '
         IF (lwp) WRITE(numout,*) 'Activating the diurnal cycle in PISCES has no effect'
      ENDIF
      !
      ! Variable PAR at the surface of the ocean
      ! ---- ------------------------------------------------------------------------------------------
!# ifdef NEMO
!      IF( ln_varpar ) THEN
!         IF(lwp) WRITE(numout,*)
!         IF(lwp) WRITE(numout,*) '   ==>>>   initialize variable par fraction (ln_varpar=T)'
!         !
!         ALLOCATE( par_varsw(A2D(0)) )
!         !
!         ALLOCATE( sf_par(1), STAT=ierr )           !* allocate and fill sf_sst (forcing structure) with sn_sst
!         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'p4z_opt_init: unable to allocate sf_par structure' )
!         !
!         CALL fld_fill( sf_par, (/ sn_par /), cn_dir, 'p4z_opt_init', 'Variable PAR fraction ', 'nampisopt' )
!                              ALLOCATE( sf_par(1)%fnow(jpi,jpj,1) )
!         IF( sn_par%ln_tint ) ALLOCATE( sf_par(1)%fdta(jpi,jpj,1,2) )
!
!         CALL iom_open (  TRIM( sn_par%clname ) , numpar )
!         ntimes_par = iom_getszuld( numpar )   ! get number of record in file
!      ENDIF
!# endif
      !
                         ekr       (:,:,:) = 0._wp
                         ekb       (:,:,:) = 0._wp
                         ekg       (:,:,:) = 0._wp
                         etot      (:,:,:) = 0._wp
                         etot_ndcy (:,:,:) = 0._wp
                         enano     (:,:,:) = 0._wp
      IF( .NOT. ln_p2z)  ediat     (:,:,:) = 0._wp
      IF( ln_p5z      )  epico     (:,:,:) = 0._wp
      IF( ln_qsr_bio  )  etot3     (:,:,:) = 0._wp
      ! 
   END SUBROUTINE p4z_opt_init


   INTEGER FUNCTION p4z_opt_alloc()
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE p4z_opt_alloc  ***
      !!----------------------------------------------------------------------
      !
      ALLOCATE( ekb(A2D(0),jpk), ekr(A2D(0),jpk), &
          &     ekg(A2D(0),jpk), STAT= p4z_opt_alloc  ) 
      !
      IF( p4z_opt_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p4z_opt_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p4z_opt_alloc

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                   No PISCES bio-model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE p4z_opt                                  ! Empty routine
   END SUBROUTINE p4z_opt
#endif

   !!======================================================================
END MODULE p4zopt