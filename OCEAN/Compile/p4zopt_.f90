










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









MODULE p4zopt
   !!======================================================================
   !!                         ***  MODULE p4zopt  ***
   !! TOP -  : Compute the light availability in the water column
   !!======================================================================
   !! History :  1.0  !  2004     (O. Aumont) Original code
   !!            2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!            3.2  !  2009-04  (C. Ethe, G. Madec)  optimisation
   !!            3.4  !  2011-06  (O. Aumont, C. Ethe) Improve light availability of nano & diat
   !!            3.*  !  2025-02  (S. Maishal, R. Person)  and optimization
   !!----------------------------------------------------------------------
   !!   p4z_opt       : light availability in the water column
   !!----------------------------------------------------------------------
   USE trc            ! tracer variables
   USE oce_trc        ! tracer-ocean share variables
   USE sms_pisces     ! Source Minus Sink of 
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
      INTEGER  ::   mi, mj, jk
      INTEGER  ::   irgb
      REAL(wp) ::   zchl
      REAL(wp) ::   zc0 , zc1 , zc2, zc3, z1_dep
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp    ) :: zdepmoy, zetmp1
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp    ) :: zqsr100, zqsr_corr
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N) :: ze0, ze1, ze2, ze3
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zpar
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zw2d
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p4z_opt')

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
      !$OMP PRIVATE(mi, mj, jk, irgb, zchl) &
      !$OMP SHARED(tr, thetanano, 0.5*EPSILON(1.e0), ln_p2z, ln_p5z, rkrgb, e3t)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF( ln_p2z ) THEN
            zchl = ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) &
              &     * 12.0 * thetanano(mi,mj,jk) + 0.5*EPSILON(1.e0) ) * 1.e6
         ELSE
            zchl =  ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnch) &
              &     + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdch) + 0.5*EPSILON(1.e0) ) * 1.e6
            IF( ln_p5z ) zchl = zchl + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppch) * 1.e6
         ENDIF
         zchl = MIN( 10. , MAX( 0.05, zchl ) )
         irgb = NINT( 41 + 20.* LOG10( zchl ) + 0.5*EPSILON(1.e0) )
         !
         ekb(mi,mj,jk) = rkrgb(1,irgb) * Hz(mi,mj,N+1-jk)
         ekg(mi,mj,jk) = rkrgb(2,irgb) * Hz(mi,mj,N+1-jk)
         ekr(mi,mj,jk) = rkrgb(3,irgb) * Hz(mi,mj,N+1-jk)
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      ! Photosynthetically Available Radiation (PAR)
      ! Two cases are considered in the following : 
      ! (1) An explicit diunal cycle is activated. In that case, mean 
      ! QSR is used as  in its current state has not been parameterized
      ! for an explicit diurnal cycle
      ! (2) no diurnal cycle of SW is active and in that case, QSR is used.
      ! --------------------------------------------
      IF( ln_trcdc2dm ) THEN                     !  diurnal cycle
         IF ( ln_p4z_dcyc ) THEN   ! Diurnal cycle in 
            !
            !
            ! SW over the ice free zone of the grid cell. This assumes that
            ! SW is zero below sea ice which is a very crude assumption that is 
            ! not fully correct with LIM3 and SI3 but no information is 
            ! currently available to do a better job. SHould be improved in the 
            ! (near) future.
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(mi, mj) &
            !$OMP SHARED(qsr_mean, fr_i, 0.5*EPSILON(1.e0), zqsr_corr)
            DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zqsr_corr(mi,mj) = max(1.e-10,rho0*Cp*srflx(mi,mj)) / ( 1.-fr_i(mi,mj) + 0.5*EPSILON(1.e0) )
            END DO   ;   END DO
            !$OMP END PARALLEL DO
            !
            CALL p4z_opt_par( kt, Kmm, zqsr_corr, ze1, ze2, ze3, pqsr100 = zqsr100 )
            !
            ! Used PAR is computed for each phytoplankton species
            ! etot_ndcy is PAR at level jk averaged over 24h.
            ! Due to their size, they have different light absorption characteristics
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(mi, mj, jk) &
            !$OMP SHARED(etot_ndcy, ze1, ze2, ze3)
            DO jk= 1, nksr  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               etot_ndcy(mi,mj,jk) = ze1(mi,mj,jk) + ze2(mi,mj,jk) + ze3(mi,mj,jk)
            END DO   ;   END DO   ;   END DO
            !$OMP END PARALLEL DO
            !
            ! SW over the ice free zone of the grid cell. This assumes that
            ! SW is zero below sea ice which is a very crude assumption that is 
            ! not fully correct with LIM3 and SI3 but no information is 
            ! currently available to do a better job. SHould be improved in the 
            ! (near) future.
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(mi, mj) &
            !$OMP SHARED(qsr, fr_i, 0.5*EPSILON(1.e0), zqsr_corr)
            DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zqsr_corr(mi,mj) = max(1.e-10,rho0*Cp*srflx(mi,mj)) / ( 1.-fr_i(mi,mj) + 0.5*EPSILON(1.e0) )
            END DO   ;   END DO
            !$OMP END PARALLEL DO
            !
            CALL p4z_opt_par( kt, Kmm, zqsr_corr, ze1, ze2, ze3 )
            !
            ! Total PAR computation at level jk that includes the diurnal cycle
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(mi, mj, jk) &
            !$OMP SHARED(etot, enano, ze1, ze2, ze3)
            DO jk= 1, nksr  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               etot (mi,mj,jk) = ze1(mi,mj,jk) + ze2(mi,mj,jk) + ze3(mi,mj,jk)
               enano(mi,mj,jk) = 1.85 * ze1(mi,mj,jk) + 0.69 * ze2(mi,mj,jk) &
                               + 0.46 * ze3(mi,mj,jk)
            END DO   ;   END DO   ;   END DO
            !$OMP END PARALLEL DO
            IF( .NOT. ln_p2z ) THEN
               !$OMP PARALLEL DO &
               !$OMP PRIVATE(mi, mj, jk) &
               !$OMP SHARED(ediat, ze1, ze2, ze3)
               DO jk= 1, nksr  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                  ediat(mi,mj,jk) = 1.62 * ze1(mi,mj,jk) + 0.74 * ze2(mi,mj,jk) &
                                  + 0.63 * ze3(mi,mj,jk)
               END DO   ;   END DO   ;   END DO
               !$OMP END PARALLEL DO
               IF( ln_p5z ) THEN
                  !$OMP PARALLEL DO &
                  !$OMP PRIVATE(mi, mj, jk) &
                  !$OMP SHARED(epico, ze1, ze2, ze3)
                  DO jk= 1, nksr  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                     epico(mi,mj,jk) = 1.94 * ze1(mi,mj,jk) + 0.66 * ze2(mi,mj,jk) &
                                     + 0.4 * ze3(mi,mj,jk)
                  END DO   ;   END DO   ;   END DO
                  !$OMP END PARALLEL DO
               ENDIF
            ENDIF

         ELSE ! No diurnal cycle in 

            !
            !
            ! SW over the ice free zone of the grid cell. This assumes that
            ! SW is zero below sea ice which is a very crude assumption that is 
            ! not fully correct with LIM3 and SI3 but no information is 
            ! currently available to do a better job. SHould be improved in the 
            ! (near) future.
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(mi, mj) &
            !$OMP SHARED(qsr_mean, fr_i, 0.5*EPSILON(1.e0), zqsr_corr)
            DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zqsr_corr(mi,mj) = max(1.e-10,rho0*Cp*srflx(mi,mj)) / ( 1.-fr_i(mi,mj) + 0.5*EPSILON(1.e0) )
            END DO   ;   END DO
            !$OMP END PARALLEL DO
            !
            CALL p4z_opt_par( kt, Kmm, zqsr_corr, ze1, ze2, ze3, pqsr100 = zqsr100 )
            !
            ! Used PAR is computed for each phytoplankton species
            ! etot_ndcy is PAR at level jk averaged over 24h.
            ! Due to their size, they have different light absorption characteristics
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(mi, mj, jk) &
            !$OMP SHARED(etot_ndcy, enano, ze1, ze2, ze3)
            DO jk= 1, nksr  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               etot_ndcy(mi,mj,jk) = ze1(mi,mj,jk) + ze2(mi,mj,jk) + ze3(mi,mj,jk)
               enano    (mi,mj,jk) = 1.85 * ze1(mi,mj,jk) + 0.69 * ze2(mi,mj,jk) &
                                   + 0.46 * ze3(mi,mj,jk)
            END DO   ;   END DO   ;   END DO
            !$OMP END PARALLEL DO
            IF( .NOT. ln_p2z ) THEN
               !$OMP PARALLEL DO &
               !$OMP PRIVATE(mi, mj, jk) &
               !$OMP SHARED(ediat, ze1, ze2, ze3)
               DO jk= 1, nksr  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                  ediat(mi,mj,jk) = 1.62 * ze1(mi,mj,jk) + 0.74 * ze2(mi,mj,jk) &
                                  + 0.63 * ze3(mi,mj,jk)
               END DO   ;   END DO   ;   END DO
               !$OMP END PARALLEL DO
               IF( ln_p5z ) THEN
                  !$OMP PARALLEL DO &
                  !$OMP PRIVATE(mi, mj, jk) &
                  !$OMP SHARED(epico, ze1, ze2, ze3)
                  DO jk= 1, nksr  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                     epico(mi,mj,jk) = 1.94 * ze1(mi,mj,jk) + 0.66 * ze2(mi,mj,jk) &
                                     + 0.4 * ze3(mi,mj,jk)
                  END DO   ;   END DO   ;   END DO
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
            !$OMP PRIVATE(mi, mj) &
            !$OMP SHARED(qsr, fr_i, 0.5*EPSILON(1.e0), zqsr_corr)
            DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zqsr_corr(mi,mj) = max(1.e-10,rho0*Cp*srflx(mi,mj)) / ( 1.-fr_i(mi,mj) + 0.5*EPSILON(1.e0) )
            END DO   ;   END DO
            !$OMP END PARALLEL DO
            !
            CALL p4z_opt_par( kt, Kmm, zqsr_corr, ze1, ze2, ze3 ) 
            !
            ! Total PAR computation at level jk that includes the diurnal cycle
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(mi, mj, jk) &
            !$OMP SHARED(etot, ze1, ze2, ze3)
            DO jk= 1, nksr  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               etot(mi,mj,jk) = ze1(mi,mj,jk) + ze2(mi,mj,jk) + ze3(mi,mj,jk)
            END DO   ;   END DO   ;   END DO
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
         !$OMP PRIVATE(mi, mj) &
         !$OMP SHARED(qsr, fr_i, 0.5*EPSILON(1.e0), zqsr_corr)
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zqsr_corr(mi,mj) = max(1.e-10,rho0*Cp*srflx(mi,mj)) / ( 1.-fr_i(mi,mj) + 0.5*EPSILON(1.e0) )
         END DO   ;   END DO
         !$OMP END PARALLEL DO
         !
         CALL p4z_opt_par( kt, Kmm, zqsr_corr, ze1, ze2, ze3, pqsr100 = zqsr100 ) 
         !

         ! Used PAR is computed for each phytoplankton species
         ! Due to their size, they have different light absorption characteristics
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk) &
         !$OMP SHARED(etot, enano, ze1, ze2, ze3)
         DO jk= 1, nksr  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            etot (mi,mj,jk) = ze1(mi,mj,jk) + ze2(mi,mj,jk) + ze3(mi,mj,jk)
            enano(mi,mj,jk) = 1.85 * ze1(mi,mj,jk) + 0.69 * ze2(mi,mj,jk) &
                            + 0.46 * ze3(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
         IF( .NOT. ln_p2z ) THEN
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(mi, mj, jk) &
            !$OMP SHARED(ediat, ze1, ze2, ze3)
            DO jk= 1, nksr  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               ediat(mi,mj,jk) = 1.62 * ze1(mi,mj,jk) + 0.74 &
                               * ze2(mi,mj,jk) + 0.63 * ze3(mi,mj,jk)
            END DO   ;   END DO   ;   END DO
            !$OMP END PARALLEL DO
            IF( ln_p5z ) THEN
               !$OMP PARALLEL DO &
               !$OMP PRIVATE(mi, mj, jk) &
               !$OMP SHARED(epico, ze1, ze2, ze3)
               DO jk= 1, nksr  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                  epico(mi,mj,jk) = 1.94 * ze1(mi,mj,jk) + 0.66 &
                                  * ze2(mi,mj,jk) + 0.4 * ze3(mi,mj,jk)
               END DO   ;   END DO   ;   END DO
               !$OMP END PARALLEL DO
            ENDIF
         ENDIF
         etot_ndcy(:,:,:) =  etot(:,:,:) 
      ENDIF


      ! Biophysical feedback part (computation of vertical penetration of SW)
      IF( ln_qsr_bio ) THEN!* heat flux accros w-level (used in the dynamics)
         !------------------------------------------ ------------------------
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zqsr_corr(mi,mj) = max(1.e-10,rho0*Cp*srflx(mi,mj)) 
         END DO   ;   END DO
         !
         CALL p4z_opt_par( kt, Kmm, zqsr_corr, ze1, ze2, ze3, pe0=ze0 )
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj) &
         !$OMP SHARED(etot3, qsr, tmask)
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            etot3(mi,mj,1) = max(1.e-10,rho0*Cp*srflx(mi,mj)) * tmask(mi,mj,1)
         END DO   ;   END DO
         !$OMP END PARALLEL DO
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk) &
         !$OMP SHARED(etot3, ze0, ze1, ze2, ze3, tmask)
         DO jk= 2, nksr+1  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            etot3(mi,mj,jk) = ( ze0(mi,mj,jk) + ze1(mi,mj,jk) + ze2(mi,mj,jk) &
                            + ze3(mi,mj,jk) ) * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
         !-------------------------------------------------------------------
      ENDIF
      
      ! Euphotic depth and level
      ! Two definitions of the euphotic zone are used here. 
      ! (1) The classical definition based on the relative threshold value
      ! (2) An alternative definition based on a absolute threshold value.
      ! ---------------------------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj) &
      !$OMP SHARED(neln, heup, heup_01, gdepw)
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         neln   (mi,mj) = 1
         heup   (mi,mj) = gdepw(mi,mj,2,Kmm)
         heup_01(mi,mj) = gdepw(mi,mj,2,Kmm)
      END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(etot_ndcy, tmask, zqsr100, neln, heup, heup_01, gdepw, Kmm)
      DO jk= 2, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
        IF( etot_ndcy(mi,mj,jk) * tmask(mi,mj,jk) >=  zqsr100(mi,mj) )  THEN
           neln(mi,mj) = jk+1                  ! Euphotic level : 1rst T-level strictly below Euphotic layer
           !                                   ! nb: ensure the compatibility with nmld_trc definition in trd_mld_trc_zint
           heup(mi,mj) = gdepw(mi,mj,jk+1,Kmm) ! Euphotic layer depth
        ENDIF
        IF( etot_ndcy(mi,mj,jk) * tmask(mi,mj,jk) >= 0.10 ) THEN
           heup_01(mi,mj) = gdepw(mi,mj,jk+1,Kmm) ! Euphotic layer depth (light level definition)
        ENDIF
      END DO   ;   END DO   ;   END DO
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
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(gdepw, hmld, zetmp1, zdepmoy, etot, e3t, Kmm)
      DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF( gdepw(mi,mj,jk+1,Kmm) <= hbl(mi,mj) ) THEN
            zetmp1 (mi,mj) = zetmp1 (mi,mj) + etot(mi,mj,jk) &
                           * Hz(mi,mj,N+1-jk) ! Actual PAR for remineralisation
            zdepmoy(mi,mj) = zdepmoy(mi,mj) + Hz(mi,mj,N+1-jk)
         ENDIF
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      emoy(:,:,:) = etot(:,:,:) ! remineralisation
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, z1_dep) &
      !$OMP SHARED(gdepw, hmld, zetmp1, zdepmoy, 0.5*EPSILON(1.e0), emoy, Kmm)
      DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF( gdepw(mi,mj,jk+1,Kmm) <= hbl(mi,mj) ) THEN
            z1_dep = 1. / ( zdepmoy(mi,mj) + 0.5*EPSILON(1.e0) )
            emoy (mi,mj,jk) = zetmp1(mi,mj) * z1_dep
         ENDIF
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      ! Computation of the mean usable light for the different phytoplankton
      ! groups based on their absorption characteristics.
      zetmp1(:,:)   = 0.e0
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(gdepw, hmld, heup_01, zetmp1, enano, e3t, Kmm)
      DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF( gdepw(mi,mj,jk+1,Kmm) <= MIN(hbl(mi,mj), heup_01(mi,mj)) ) THEN
            zetmp1(mi,mj) = zetmp1(mi,mj) + enano(mi,mj,jk) &
                          * Hz(mi,mj,N+1-jk) ! Nanophytoplankton
         ENDIF
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      enanom(:,:,:) = enano(:,:,:)
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, z1_dep) &
      !$OMP SHARED(gdepw, hmld, zetmp1, zdepmoy, 0.5*EPSILON(1.e0), enanom, Kmm)
      DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF( gdepw(mi,mj,jk+1,Kmm) <= hbl(mi,mj) ) THEN
            z1_dep = 1. / ( zdepmoy(mi,mj) + 0.5*EPSILON(1.e0) )
            enanom(mi,mj,jk) = zetmp1(mi,mj) * z1_dep
         ENDIF
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      IF( .NOT. ln_p2z ) THEN
         ! Diatoms when using -operational or -QUOTA
         zetmp1(:,:) = 0.e0
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk) &
         !$OMP SHARED(gdepw, hmld, heup_01, zetmp1, ediat, e3t, Kmm)
         DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            IF( gdepw(mi,mj,jk+1,Kmm) <= MIN(hbl(mi,mj), heup_01(mi,mj)) ) THEN
               zetmp1(mi,mj) = zetmp1(mi,mj) + ediat(mi,mj,jk) * Hz(mi,mj,N+1-jk) ! Diatoms
            ENDIF
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
         !
         ediatm(:,:,:) = ediat(:,:,:)
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk, z1_dep) &
         !$OMP SHARED(gdepw, hmld, zetmp1, zdepmoy, 0.5*EPSILON(1.e0), ediatm, Kmm)
         DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            IF( gdepw(mi,mj,jk+1,Kmm) <= hbl(mi,mj) ) THEN
               z1_dep = 1. / ( zdepmoy(mi,mj) + 0.5*EPSILON(1.e0) )
               ediatm(mi,mj,jk) = zetmp1(mi,mj) * z1_dep
            ENDIF
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
      ENDIF
      IF( ln_p5z ) THEN
         ! Picophytoplankton when using -QUOTA
         zetmp1(:,:) = 0.e0
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk) &
         !$OMP SHARED(gdepw, hmld, heup_01, zetmp1, epico, e3t, Kmm)
         DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            IF( gdepw(mi,mj,jk+1,Kmm) <= MIN(hbl(mi,mj), heup_01(mi,mj)) ) THEN
               zetmp1(mi,mj) = zetmp1(mi,mj) + epico(mi,mj,jk) &
                             * Hz(mi,mj,N+1-jk)
            ENDIF
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
         !
         epicom(:,:,:) = epico(:,:,:)
         !
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk, z1_dep) &
         !$OMP SHARED(gdepw, hmld, zetmp1, zdepmoy, 0.5*EPSILON(1.e0), epicom, Kmm)
         DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            IF( gdepw(mi,mj,jk+1,Kmm) <= hbl(mi,mj) ) THEN
               z1_dep = 1. / ( zdepmoy(mi,mj) + 0.5*EPSILON(1.e0) )
               epicom(mi,mj,jk) = zetmp1(mi,mj) * z1_dep
            ENDIF
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
      ENDIF
      !
      IF( l_dia_par .OR. l_diaadd ) THEN
        zetmp1(:,:)   = 0.e0
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(mi, mj, jk) &
        !$OMP SHARED(gdepw, hmld, zetmp1, etot_ndcy, e3t, Kmm)
        DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           IF( gdepw(mi,mj,jk+1,Kmm) <= hbl(mi,mj) ) THEN
              zetmp1(mi,mj) = zetmp1(mi,mj) &
                &       + etot_ndcy(mi,mj,jk) * Hz(mi,mj,N+1-jk) ! Par averaged over 24h for production
           ENDIF
        END DO   ;   END DO   ;   END DO
        !$OMP END PARALLEL DO
      ENDIF
      IF( .false. .AND.  knt == nrdttrc ) THEN
         IF( l_dia_heup ) THEN
           ALLOCATE( zw2d(-2:Lm+3+padd_X,-2:Mm+3+padd_E) )  ;  zw2d(:,:) = 0._wp
           zw2d(Istrp:Iendp,Jstrp:Jendp) = heup(Istrp:Iendp,Jstrp:Jendp) * tmask(Istrp:Iendp,Jstrp:Jendp,1)
           CALL iom_put( "Heup", zw2d )  ! Euphotic layer depth
           DEALLOCATE( zw2d ) 
        ENDIF
        IF( l_dia_par ) THEN   ! diagnostic : PAR with no diurnal cycle
           ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
           DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
              zw3d(mi,mj,N+1-jk) = etot_ndcy(mi,mj,jk) * tmask(mi,mj,jk)
           END DO   ;   END DO   ;   END DO
           DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
              IF( gdepw(mi,mj,jk+1,Kmm) <= hbl(mi,mj) ) THEN
                 z1_dep = 1. / ( zdepmoy(mi,mj) + 0.5*EPSILON(1.e0) )
                 zw3d(mi,mj,N+1-jk) = zetmp1(mi,mj) * z1_dep
              ENDIF
           END DO   ;   END DO   ;   END DO
           DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
              zw3d(mi,mj,N) = zw3d(mi,mj,N) * ( 1._wp - fr_i(mi,mj) )
           END DO   ;   END DO
           CALL iom_put( "PAR", zw3d ) 
           DEALLOCATE( zw3d ) 
        ENDIF
      ENDIF
      !
      DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         bioFlux(mi,mj,N+1-jk,Netot) = etot_ndcy(mi,mj,jk) * tmask(mi,mj,jk) ! PAR
      END DO   ;   END DO   ;   END DO
      DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF( gdepw(mi,mj,jk+1,Kmm) <= hbl(mi,mj) ) THEN
            z1_dep = 1. / ( zdepmoy(mi,mj) + 0.5*EPSILON(1.e0) )
            bioFlux(mi,mj,N+1-jk,Netot) = zetmp1(mi,mj) * z1_dep
          ENDIF
      END DO   ;   END DO   ;   END DO
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         bioFlux(mi,mj,N,Netot) = bioFlux(mi,mj,N,Netot) * ( 1._wp - fr_i(mi,mj) )
      END DO   ;   END DO
      !
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         bioVSink(mi,mj,Nheup) = heup(mi,mj) * tmask(mi,mj,1) ! euphotic layer
      END DO   ;   END DO
      IF( .false. )   CALL timing_stop('p4z_opt')
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
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp)    , INTENT(in   )           ::   pqsr              ! shortwave
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N), INTENT(inout)           ::   pe1 , pe2 , pe3   ! PAR ( R-G-B)
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N), INTENT(inout), OPTIONAL ::   pe0               !
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp)    , INTENT(  out), OPTIONAL ::   pqsr100           !
      !
      INTEGER    ::   mi, mj, jk, jkm1     ! dummy loop indices _*_
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp) ::  zqsr ! shortwave
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
         DO jk= 2, nksr + 1 ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            jkm1 = jk-1
            pe0(mi,mj,jk) = pe0(mi,mj,jk-1) * EXP( -Hz(mi,mj,N+1-jkm1) * xsi0r )
            pe1(mi,mj,jk) = pe1(mi,mj,jk-1) * EXP( -ekb  (mi,mj,jk-1 )        )
            pe2(mi,mj,jk) = pe2(mi,mj,jk-1) * EXP( -ekg  (mi,mj,jk-1 )        )
            pe3(mi,mj,jk) = pe3(mi,mj,jk-1) * EXP( -ekr  (mi,mj,jk-1 )        )
        END DO   ;   END DO   ;   END DO
        !
      ELSE   ! T- level
        !
        pe1(:,:,1) = zqsr(:,:) * EXP( -0.5 * ekb(:,:,1) )
        pe2(:,:,1) = zqsr(:,:) * EXP( -0.5 * ekg(:,:,1) )
        pe3(:,:,1) = zqsr(:,:) * EXP( -0.5 * ekr(:,:,1) )
        !
        DO jk= 2, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           pe1(mi,mj,jk) = pe1(mi,mj,jk-1) * EXP( -0.5 * ( ekb(mi,mj,jk-1) + ekb(mi,mj,jk) ) )
           pe2(mi,mj,jk) = pe2(mi,mj,jk-1) * EXP( -0.5 * ( ekg(mi,mj,jk-1) + ekg(mi,mj,jk) ) )
           pe3(mi,mj,jk) = pe3(mi,mj,jk-1) * EXP( -0.5 * ( ekr(mi,mj,jk-1) + ekr(mi,mj,jk) ) )
        END DO   ;   END DO   ;   END DO
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
      INTEGER  :: mi,mj
      REAL(wp) :: zcoef
      !!---------------------------------------------------------------------
      !
      IF( .false. )  CALL timing_start('p4z_optsbc')
      !
      !
      IF( .false. )  CALL timing_stop('p4z_optsbc')
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
      !!----------------------------------------------------------------------
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*)
         WRITE(stdout,*) 'p4z_opt_init : '
         WRITE(stdout,*) '~~~~~~~~~~~~ '
      ENDIF
!# ifdef NEMO
!      REWIND(numnatp_ref);READ(numnatp_ref,nampisopt,IOSTAT=ios);CALL ctl_nam(ios,"nampisopt (ref)",.TRUE.)
!      REWIND(numnatp_cfg);READ(numnatp_cfg,nampisopt,IOSTAT=ios);CALL ctl_nam(ios,"nampisopt (cfg)",.FALSE.)
!      IF(mynode .eq. 0) WRITE ( numonp, nampisopt )
!# endif
      ln_p4z_dcyc =  .false. ! Diurnal cycle in 
      ln_varpar   = .FALSE.  !
      parlux      =  0.43    ! Fraction of shortwave as PAR
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*) '      Namelist : nampisopt '
         WRITE(stdout,*) '      PAR as a variable fraction of SW       ln_varpar   = ', ln_varpar
         WRITE(stdout,*) '      Default value for the PAR fraction     parlux      = ', parlux
         WRITE(stdout,*) '      Activate the diurnal cycle in PISCES   ln_p4z_dcyc = ', ln_p4z_dcyc
      ENDIF
      !
      !
      xparsw = parlux / 3.0
      xsi0r  = 1.e0 / rn_si0

      ! Warning : activate the diurnal cycle with no diurnal cycle in the forcing fields makes no sense
      ! That does not produce a bug because the model does not use the flag but a warning is necessary
      ! ----------------------------------------------------------------------------------------------
      IF ( ln_p4z_dcyc .AND. .false. ) THEN
         IF (mynode .eq. 0) WRITE(stdout,*) 'No diurnal cycle in the PAR forcing field '
         IF (mynode .eq. 0) WRITE(stdout,*) 'Activating the diurnal cycle in PISCES has no effect'
      ENDIF
      !
      ! Variable PAR at the surface of the ocean
      ! ---- ------------------------------------------------------------------------------------------
!# ifdef NEMO
!      IF( ln_varpar ) THEN
!         IF(mynode .eq. 0) WRITE(stdout,*)
!         IF(mynode .eq. 0) WRITE(stdout,*) '   ==>>>   initialize variable par fraction (ln_varpar=T)'
!         !
!         ALLOCATE( par_varsw(Istrp:Iendp,Jstrp:Jendp) )
!         !
!         ALLOCATE( sf_par(1), STAT=ierr )           !* allocate and fill sf_sst (forcing structure) with sn_sst
!         IF( ierr > 0 )   CALL ctl_stop( 'STOP', 'p4z_opt_init: unable to allocate sf_par structure' )
!         !
!         CALL fld_fill( sf_par, (/ sn_par /), cn_dir, 'p4z_opt_init', 'Variable PAR fraction ', 'nampisopt' )
!                              ALLOCATE( sf_par(1)%fnow(Lm,Mm,1) )
!         IF( sn_par%ln_tint ) ALLOCATE( sf_par(1)%fdta(Lm,Mm,1,2) )
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
      ALLOCATE( ekb(Istrp:Iendp,Jstrp:Jendp,N), ekr(Istrp:Iendp,Jstrp:Jendp,N), &
          &     ekg(Istrp:Iendp,Jstrp:Jendp,N), STAT= p4z_opt_alloc  ) 
      !
      IF( p4z_opt_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p4z_opt_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p4z_opt_alloc


   !!======================================================================
END MODULE p4zopt
