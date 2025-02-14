










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









MODULE trcsink
   !!======================================================================
   !!                         ***  MODULE trcsink  ***
   !! TOP :  vertical flux of particulate matter due to gravitational sinking
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!             3.5  !  2012-07  (O. Aumont) Introduce potential time-splitting
   !!             4.0  !  2018-12  (O. Aumont) Generalize the  code to make it usable by any model
   !!             5.0  !  2023-10  (C. Ethe ) Introduce semi-lagragian sinking scheme
   !!----------------------------------------------------------------------
   !!   trc_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE lib_mpp
   USE sms_pisces

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_sink
   PUBLIC trc_sink_ini

   LOGICAL, PUBLIC :: ln_sink_mus    !: MUSCL sinkin scheme
   LOGICAL, PUBLIC :: ln_sink_slg    !: Semi-Lagrangian sinkin scheme
   INTEGER, PUBLIC :: nitermax       !: Maximum number of iterations for sinking ( ln_sink_mus )
   INTEGER, PUBLIC :: nn_sink_lbc    !: Type of boundary conditons for sinking ( ln_sink_slg )

   INTEGER, PARAMETER ::   np_MUS = 1   ! MUSCL sinking scheme 
   INTEGER, PARAMETER ::   np_SLG = 2   ! Semi-Lagrangian sinking scheme

   INTEGER  ::   nsnk     ! user choice of the type of sinking scheme

   !! * Substitutions











































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcsink.F90 10069 2018-08-28 14:12:24Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   !!----------------------------------------------------------------------
   !!   'standard sinking parameterisation'                  ???
   !!----------------------------------------------------------------------

   SUBROUTINE trc_sink ( kt, Kbb, Kmm, pwsink, psinkflx, jp_tra, rsfact )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER , INTENT(in)  :: kt
      INTEGER , INTENT(in)  :: Kbb, Kmm
      INTEGER , INTENT(in)  :: jp_tra    ! tracer index index      
      REAL(wp), INTENT(in)  :: rsfact    ! time step duration
      REAL(wp), INTENT(in)   , DIMENSION(Istrp:Iendp,Jstrp:Jendp,N) :: pwsink
      REAL(wp), INTENT(inout), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N) :: psinkflx
      !
      INTEGER  ::   mi, mj, jk
      INTEGER , ALLOCATABLE, DIMENSION(:,:)   :: iiter
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zwsink
      REAL(wp) ::   zfact, zwsmax, zmax
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('trc_sink')
      !
      !

      !                       !----------------------------!
      SELECT CASE( nsnk )     !  Sinking type              !
      !                       !----------------------------!
      CASE( np_MUS )                                     !==  MUSCL sinking scheme  ==!
         !
         ALLOCATE( iiter( Istrp:Iendp,Jstrp:Jendp ) )    ;     ALLOCATE( zwsink(Istrp:Iendp,Jstrp:Jendp,N+1) )
         ! OA This is (I hope) a temporary solution for the problem that may 
         ! OA arise in specific situation where the CFL criterion is broken 
         ! OA for vertical sedimentation of particles. To avoid this, a time
         ! OA splitting algorithm has been coded. A specific maximum
         ! OA iteration number is provided and may be specified in the namelist 
         ! OA This is to avoid very large iteration number when explicit free
         ! OA surface is used (for instance). When niter?max is set to 1, 
         ! OA this computation is skipped. The crude old threshold method is 
         ! OA then applied. This also happens when niter exceeds nitermax.
         IF( nitermax == 1 ) THEN
            iiter(:,:) = 1
         ELSE
            DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               iiter(mi,mj) = 1
               DO jk = 1, N-1
                  IF( tmask(mi,mj,jk) == 1.0 ) THEN
                      zwsmax =  0.5 * Hz(mi,mj,N+1-jk) * day2sec / rsfact
                      iiter(mi,mj) =  MAX( iiter(mi,mj), INT( pwsink(mi,mj,jk) / zwsmax ) + 1 )
                  ENDIF
               END DO
            END DO   ;   END DO
            iiter(:,:) = MIN( iiter(:,:), nitermax )
         ENDIF

         DO jk= 1, N-1  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zwsmax = 0.5 * Hz(mi,mj,N+1-jk) * day2sec / rsfact
            zwsink(mi,mj,jk+1) = -MIN( pwsink(mi,mj,jk), zwsmax * REAL( iiter(mi,mj), wp ) ) / day2sec
         END DO   ;   END DO   ;   END DO
         zwsink(:,:,1)     = 0._wp
         zwsink(:,:,N+1) = 0._wp

         !  Initializa to zero all the sinking arrays 
         !  -----------------------------------------
         psinkflx(:,:,:) = 0.e0

         !   Compute the sedimentation term using trc_sink2_mus for the considered sinking particle
         !   -----------------------------------------------------
         CALL trc_sink2_mus( Kbb, Kmm, zwsink, psinkflx, jp_tra, iiter, rsfact )
         !
         DEALLOCATE( iiter )   ;     DEALLOCATE( zwsink ) 
         !
      CASE( np_SLG )                                     !==  Semi-Lagrangian sinking scheme ==!
         !
         !  Initializa to zero all the sinking arrays 
         !  -----------------------------------------
         psinkflx(:,:,:) = 0.e0

         !   Compute the sedimentation term using trc_sink2_slg for the considered sinking particle
         !   -----------------------------------------------------
         CALL trc_sink2_slg( Kbb, Kmm, pwsink, psinkflx, jp_tra, rsfact )
         !
      END SELECT
      !
      IF( .false. )   CALL timing_stop('trc_sink')
      !
   END SUBROUTINE trc_sink

   SUBROUTINE trc_sink2_mus( Kbb, Kmm, pwsink, psinkflx, jp_tra, kiter, rsfact )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sink2_mus  ***
      !!
      !! ** Purpose :   Compute the sedimentation terms for the various sinking
      !!     particles. The scheme used to compute the trends is based
      !!     on MUSCL.
      !!
      !! ** Method  : - this ROUTINE compute not exactly the advection but the
      !!      transport term, i.e.  div(u*tra).
      !!---------------------------------------------------------------------
      INTEGER,  INTENT(in   )                          ::   Kbb, Kmm  ! time level indices
      INTEGER,  INTENT(in   )                          ::   jp_tra    ! tracer index index      
      REAL(wp), INTENT(in   )                          ::   rsfact    ! duration of time step
      INTEGER,  INTENT(in   ), DIMENSION(Istrp:Iendp,Jstrp:Jendp)       ::   kiter     ! number of iterations for time-splitting 
      REAL(wp), INTENT(in   ), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N+1) ::   pwsink    ! sinking speed
      REAL(wp), INTENT(inout), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N+1) ::   psinkflx  ! sinking fluxe
      !
      INTEGER  ::  mi, mj, jk, jn, jt
      REAL(wp) ::  zigma,z0w,zign, zflx, zstep, zzwx, zzwy, zalpha
      REAL(wp) ::  ztraz, ztraz_km1
      REAL(wp), DIMENSION(N)   :: ztrb 
      REAL(wp), DIMENSION(N+1) :: zakz, zsinking 
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('trc_sink2_mus')
      !
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ! Vertical advective flux
         zstep = rsfact / REAL( kiter(mi,mj), wp ) / 2.
         DO jt = 1, kiter(mi,mj)
            zakz (:) = 0.e0
            DO jk = 1, N
               ztrb(jk) = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jp_tra)
            ENDDO
            DO jn = 1, 2
               !              
               ztraz_km1 = ( ztrb(1) - ztrb(2) ) * tmask(mi,mj,2)
               DO jk = 2, N-1
                  ztraz     = ( ztrb(jk) - ztrb(jk+1) ) * tmask(mi,mj,jk+1)
                  zign      = 0.25 + SIGN( 0.25_wp, ztraz_km1 * ztraz )
                  zakz(jk)  = ( ztraz_km1 + ztraz ) * zign
                  zakz(jk)  = SIGN( 1.0_wp, zakz(jk) ) *        &
                     &        MIN( ABS( zakz(jk) ), 2. * ABS(ztraz), 2. * ABS(ztraz_km1) )
                  ztraz_km1 = ztraz
               END DO
      
               ! vertical advective flux
               zsinking(1)     = 0.e0
               zsinking(N+1) = 0.e0
               DO jk = 1, N-1
                  z0w      = SIGN( 0.5_wp, pwsink(mi,mj,jk+1) )
                  zalpha   = 0.5 + z0w 
                  zigma    = z0w - 0.5 * pwsink(mi,mj,jk+1) * zstep / e3w(mi,mj,jk+1,Kmm)
                  zzwx     = ztrb(jk+1) + zigma * zakz(jk+1)
                  zzwy     = ztrb(jk) + zigma * zakz(jk)
                  zsinking(jk+1) = -pwsink(mi,mj,jk+1) * ( zalpha * zzwx + (1.0 - zalpha) * zzwy ) * zstep
                  zflx     = ( zsinking(jk) - zsinking(jk+1) ) / Hz(mi,mj,N+1-jk)
                  ztrb(jk) = ztrb(jk) + zflx * tmask(mi,mj,jk)
               END DO
            END DO
            DO jk = 1, N-1
               zflx = ( zsinking(jk) - zsinking(jk+1) ) / Hz(mi,mj,N+1-jk)
               t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jp_tra) = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jp_tra) + 2. * zflx * tmask(mi,mj,jk)
               psinkflx(mi,mj,jk)      = psinkflx(mi,mj,jk) + 2. * zsinking(jk)
            END DO
            psinkflx(mi,mj,N) = psinkflx(mi,mj,N) + 2. * zsinking(N)
         END DO
      END DO   ;   END DO
      !
      IF( .false. )  CALL timing_stop('trc_sink2_mus')
      !
   END SUBROUTINE trc_sink2_mus

   SUBROUTINE trc_sink2_slg( Kbb, Kmm, pwsink, psinkflx, jp_tra, rsfact )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_sink2_slg  ***
      !!
      !! ** Purpose :   Compute the sedimentation terms for the various sinking particles.
      !!                The scheme used to compute the trends is based on 
      !!                a semi-Lagrangian advective flux algorithm
      !!
      !! ** Method  : - uses a parabolic,  vertical reconstructuion of the suspended particle 
      !!                in  the water column with PPT/WENO constraints to avoid oscillation
      !!                                                                     
      !!  References:                                                         
      !!                                                                      
      !!  Colella, P. and P. Woodward, 1984: The piecewise parabolic method   
      !!    (PPM) for gas-dynamical simulations, J. Comp. Phys., 54, 174-201. 
      !!                                                                      
      !!  Liu, X.D., S. Osher, and T. Chan, 1994: Weighted essentially        
      !!    nonoscillatory shemes, J. Comp. Phys., 115, 200-212.              
      !!                                                                      
      !!  Warner, J.C., C.R. Sherwood, R.P. Signell, C.K. Harris, and H.G.   
      !!    Arango, 2008:  Development of a three-dimensional,  regional,    
      !!    coupled wave, current, and sediment-transport model, Computers   
      !!    & Geosciences, 34, 1284-1306.                                     
      !!---------------------------------------------------------------------
      INTEGER,  INTENT(in   )                          ::   Kbb, Kmm  ! time level indices
      INTEGER,  INTENT(in   )                          ::   jp_tra    ! tracer index index      
      REAL(wp), INTENT(in   ), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N+1) ::   pwsink    ! sinking speed
      REAL(wp), INTENT(inout), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N+1) ::   psinkflx  ! sinking fluxe
      REAL(wp), INTENT(in   )                          ::   rsfact    ! duration of time step
      !
      INTEGER  :: mi, mj, jk, ik, jkm1, jkp1
      REAL(wp) :: zcff, zcu, zcffL, zcffR, zdltL, zdltR, zflx
      REAL(wp) :: zHz_inv, zHz_inv2, zHz_inv3
      !
      INTEGER , DIMENSION(Istrp:Iendp,N) :: ksource
      REAL(wp), DIMENSION(Istrp:Iendp,N+1) :: zFC
      REAL(wp), DIMENSION(Istrp:Iendp,N) :: zqR
      REAL(wp), DIMENSION(Istrp:Iendp,N) :: zqL
      REAL(wp), DIMENSION(Istrp:Iendp,N) :: zWR
      REAL(wp), DIMENSION(Istrp:Iendp,N) :: zWL
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('trc_sink2_slg')
      !
      zqR(:,:) = 0._wp
      zqL(:,:) = 0._wp
      zWR(:,:) = 0._wp
      zWL(:,:) = 0._wp
      zFC(:,:) = 0._wp

      !-----------------------------------------------------------------------
      !  Vertical sinking of particle concentration
      !-----------------------------------------------------------------------
      DO mj=Jstrp,Jendp                                  !  i-k slices loop  !
         !  Compute semi-Lagrangian flux due to sinking.
         DO jk= 2, N, 1  ; DO mi=Istrp,Iendp
            jkm1 = jk-1
            zHz_inv2   = 1._wp / ( Hz(mi,mj,N+1-jk) + Hz(mi,mj,N+1-jkm1)  )
            zFC(mi,jk) = ( t(mi,mj,N+1-jkm1,Kbb,itemp+ntrc_salt+jp_tra) - t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jp_tra) ) * zHz_inv2  
         END DO   ;   END DO
         !
         DO jk= 2, N-1, 1  ; DO mi=Istrp,Iendp
            !
            jkp1 = jk+1
            zdltR = Hz(mi,mj,N+1-jk) * zFC(mi,jk)
            zdltL = Hz(mi,mj,N+1-jk) * zFC(mi,jk+1)
            zcff  = Hz(mi,mj,N+1-jkp1) + 2. * Hz(mi,mj,N+1-jk) + Hz(mi,mj,N+1-jkm1)
            zcffR = zcff * zFC(mi,jk)
            zcffL = zcff * zFC(mi,jk+1)
            !
            !  Apply PPM monotonicity constraint to prevent oscillations within the grid box.
            IF( zdltR * zdltL <= 0._wp ) THEN
                zdltR = 0._wp
                zdltL = 0._wp
            ELSE IF( ABS( zdltR ) >= zcffL ) THEN
                zdltR = zcffL
            ELSE IF( ABS( zdltL ) > ABS( zcffR ) ) THEN
                zdltL = zcffR
            ENDIF
            !
            !  Compute right and left side values (qR,qL) of parabolic segments
            !  within grid box Hz(k); (WR,WL) are measures of quadratic variations.
            !
            !  NOTE: Although each parabolic segment is monotonic within its grid
            !        box, monotonicity of the whole profile is not guaranteed,
            !        because qL(k+1)-qR(k) may still have different sign than
            !        trb(k+1)-trb(k).  This possibility is excluded, after qL and qR
            !        are reconciled using WENO procedure.
            !
            zHz_inv3   = 1._wp / ( Hz(mi,mj,N+1-jk) + Hz(mi,mj,N+1-jkm1) + Hz(mi,mj,N+1-jkp1) )
            zcff       = ( zdltR - zdltL ) * zHz_inv3
            zdltR      = zdltR - zcff * Hz(mi,mj,N+1-jkm1)
            zdltL      = zdltL + zcff * Hz(mi,mj,N+1-jkp1)
            zqR(mi,jk) = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jp_tra) + zdltR
            zqL(mi,jk) = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jp_tra) - zdltL
            zWR(mi,jk) = ( 2._wp * zdltR - zdltL )**2
            zWL(mi,jk) = ( zdltR - 2._wp * zdltL )**2
         END DO   ;   END DO
         !
         zcff = 1.e-14
         DO jk= 2, N-2, 1  ; DO mi=Istrp,Iendp
            zdltL        = MAX( zcff, zWL(mi,jk  ))
            zdltR        = MAX( zcff, zWR(mi,jk-1))
            zqR(mi,jk)   = ( zdltR * zqR(mi,jk) + zdltL * zqL(mi,jk-1) ) / ( zdltR + zdltL )
            zqL(mi,jk-1) = zqR(mi,jk)
         END DO   ;   END DO
         !
         SELECT CASE( nn_sink_lbc )     
      
         CASE( 1 )         !  linear continuation
            DO mi=Istrp,Iendp                                  
               zFC(mi,1)   = 0.              ! no-flux boundary condition
               !
               zqL(mi,1)   = zqR(mi,2)
               zqR(mi,1)   = 2._wp * t(mi,mj,N+1-1,Kbb,itemp+ntrc_salt+jp_tra) - zqL(mi,1)
               !
               zqR(mi,N) = zqL(mi,N-1)
               zqL(mi,N) = 2._wp * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jp_tra) - zqR(mi,N)
            END DO
         CASE( 2 )         !  Neumann conditions
            DO mi=Istrp,Iendp                                  
               zFC(mi,1)   = 0.              ! no-flux boundary condition
               !
               zqL(mi,1) = zqR(mi,2)
               zqR(mi,1) = 1.5_wp * t(mi,mj,N+1-1,Kbb,itemp+ntrc_salt+jp_tra) - 0.5_wp * zqL(mi,1)
               !
               zqR(mi,N) = zqL(mi,N-1)
               zqL(mi,N) = 1.5_wp * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jp_tra) - 0.5_wp * zqR(mi,N)
            END DO
         CASE DEFAULT    !  default strictly monotonic conditions
            DO mi=Istrp,Iendp                                  
               zFC(mi,1)=0.              ! no-flux boundary condition
               !
               zqR(mi,1) = t(mi,mj,N+1-1,Kbb,itemp+ntrc_salt+jp_tra)         ! default strictly monotonic
               zqL(mi,1) = t(mi,mj,N+1-1,Kbb,itemp+ntrc_salt+jp_tra)         ! conditions
               zqR(mi,2) = t(mi,mj,N+1-1,Kbb,itemp+ntrc_salt+jp_tra)
               !
               ik =  N
               IF( ik > 1 ) THEN
                  zqR(mi,ik  ) = t(mi,mj,N+1-ik,Kbb,itemp+ntrc_salt+jp_tra)               
                  zqL(mi,ik-1) = t(mi,mj,N+1-ik,Kbb,itemp+ntrc_salt+jp_tra)  ! bottom grid boxes are re-assumed to be piecewise constant.
                  zqL(mi,ik  ) = t(mi,mj,N+1-ik,Kbb,itemp+ntrc_salt+jp_tra)           
               ENDIF
            END DO
         END SELECT
!
         !  Apply monotonicity constraint again, since the reconciled interfacial
         !  values may cause a non-monotonic behavior of the parabolic segments
         !  inside the grid box.
         DO jk= 1, N, 1  ; DO mi=Istrp,Iendp
            zdltR = zqR(mi,jk) - t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jp_tra)
            zdltL = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jp_tra) - zqL(mi,jk)
            zcffR = 2._wp * zdltR
            zcffL = 2._wp * zdltL
            IF( zdltR * zdltL < 0._wp ) THEN
               zdltR = 0._wp
               zdltL = 0._wp
            ELSE IF( ABS( zdltR ) > ABS( zcffL ) ) THEN
               zdltR = zcffL
            ELSE IF( ABS( zdltL ) > ABS( zcffR ) ) THEN
               zdltL = zcffR
            ENDIF
            zqR(mi,jk) = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jp_tra) + zdltR
            zqL(mi,jk) = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jp_tra) - zdltL
         END DO   ;   END DO

         !  After this moment reconstruction is considered complete. The next
         !  stage is to compute vertical advective fluxes, FC. It is expected
         !  that sinking may occurs relatively fast, the algorithm is designed
         !  to be free of CFL criterion, which is achieved by allowing
         !  integration bounds for semi-Lagrangian advective flux to use as
         !  many grid boxes in upstream direction as necessary.

         !  In the two code segments below, WL is the z-coordinate of the
         !  departure point for grid box interface z_w with the same indices;
         !  FC is the finite volume flux; ksource(:,k) is index of vertical
         !  grid box which contains the departure point (restricted by N(ng)).
         !  During the search: also add in content of whole grid boxes
         !  participating in FC.

         DO jk= 1, N-1, 1  ; DO mi=Istrp,Iendp
            zcff           = rsfact * ABS( pwsink(mi,mj,jk) ) / day2sec * tmask(mi,mj,jk)
            zFC(mi,jk+1)   = 0._wp
            zWL(mi,jk)     = -gdepw(mi,mj,jk+1,Kmm) + zcff 
            zWR(mi,jk)     = Hz(mi,mj,N+1-jk) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jp_tra)
            ksource(mi,jk) = jk
         END DO   ;   END DO

         DO jk = 1, N
            DO ik = 2, jk
               DO mi=Istrp,Iendp                                  !  i- loop
                  IF( zWL(mi,jk) > -gdepw(mi,mj,ik,Kmm) ) THEN
                     ksource(mi,jk) = ik - 1
                     zFC(mi,jk+1)   = zFC(mi,jk+1) + zWR(mi,ik)
                  ENDIF
               END DO
            ENDDO
         END DO
         !
         !  Finalize computation of flux: add fractional part.
         !
         DO jk= 1, N-1, 1  ; DO mi=Istrp,Iendp
            ik           = ksource(mi,jk)
            zHz_inv      = 1._wp / Hz(mi,mj,N+1-jk)
            zcu          = MIN( 1._wp, ( zWL(mi,jk) + gdepw(mi,mj,ik+1,Kmm) ) * zHz_inv )
            zFC(mi,jk+1) = zFC(mi,jk+1)                                  & 
               &         + Hz(mi,mj,N+1-ik) * zcu                      &
               &         * ( zqL(mi,ik) + zcu                           &
               &               * ( 0.5_wp * ( zqR(mi,ik) - zqL(mi,ik) ) &
               &                    - ( 1.5_wp - zcu ) * ( zqR(mi,ik) + zqL(mi,ik) &
               &                            - 2._wp * t(mi,mj,N+1-ik,Kbb,itemp+ntrc_salt+jp_tra) ) ) ) 
         END DO   ;   END DO
         !
         DO jk= 1, N-1, 1  ; DO mi=Istrp,Iendp
            zHz_inv = 1._wp / Hz(mi,mj,N+1-jk)
            zflx    = ( zFC(mi,jk) - zFC(mi,jk+1) ) * zHz_inv
            t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jp_tra) = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jp_tra) + zflx
            psinkflx(mi,mj,jk)      = zFC(mi,jk)
         END DO   ;   END DO
         !
      END DO
      !
      IF( .false. )  CALL timing_stop('trc_sink2_slg')
      !
   END SUBROUTINE trc_sink2_slg


  SUBROUTINE trc_sink_ini
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trc_sink_ini ***
      !!
      !! ** Purpose :   read  namelist options 
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ios   ! Local integer output status for namelist read
      !!
      NAMELIST/nampis_snk/ ln_sink_slg, ln_sink_mus, nitermax, nn_sink_lbc
      !!----------------------------------------------------------------------
      !
      REWIND(numnatp_ref);READ(numnatp_ref,nampis_snk,IOSTAT=ios);CALL ctl_nam(ios,"nampis_snk (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,nampis_snk,IOSTAT=ios);CALL ctl_nam(ios,"nampis_snk (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE( numonp, nampis_snk )

      IF(mynode .eq. 0) THEN                     !   ! Control print
         WRITE(stdout,*)
         WRITE(stdout,*) 'trc_sink : Sedimentation of particles '
         WRITE(stdout,*) '~~~~~~~ '
         WRITE(stdout,*) '   Namelist namtrc_snk : sedimentation of particles'
         WRITE(stdout,*) '   Use MUSCL sinking scheme              ln_sink_mus     = ', ln_sink_mus
         WRITE(stdout,*) '       Maximum number of iterations          nitermax    = ', nitermax
         WRITE(stdout,*) '   Use Semi-Lagrangian sinking scheme    ln_sink_slg     = ', ln_sink_slg
         WRITE(stdout,*) '       Type of boundary conditions           nn_sink_lbc = ', nn_sink_lbc
         WRITE(stdout,*)
         IF( ln_sink_slg ) THEN
            SELECT CASE( nn_sink_lbc )             ! Type of boundary conditions
            CASE( 1 )    ;  WRITE(stdout,*) '   ==>>>   Dirichlet condition : linear continuation' 
            CASE( 2 )    ;  WRITE(stdout,*) '   ==>>>   Neumann condition '
            CASE DEFAULT ;  WRITE(stdout,*) '   ==>>>   Strictly monotonic conditions'
            END SELECT
         ENDIF
      ENDIF

      ioptio = 0              !**  Parameter control  **!
      IF( ln_sink_mus  )   ioptio = ioptio + 1
      IF( ln_sink_slg  )   ioptio = ioptio + 1
      !
      IF( ioptio /= 1 )   CALL ctl_stop( 'STOP','Choose ONE type of sinking scheme namelist namtrc_snk' )
      !
      IF( ln_sink_mus )  nsnk = np_MUS
      IF( ln_sink_slg )  nsnk = np_SLG
      !
   END SUBROUTINE trc_sink_ini


   !!======================================================================
END MODULE trcsink
