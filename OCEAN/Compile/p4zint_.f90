










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









MODULE p4zint
   !!=========================================================================
   !!                         ***  MODULE p4zint  ***
   !! TOP :    interpolation and computation of various accessory fields
   !!=========================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.*  !  2025-02  (S. Maishal, R. Person)  and optimization
   !!----------------------------------------------------------------------
   !!   p4z_int        :  interpolation and computation of various accessory fields
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !   Source Minus Sink variables

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_int  
   REAL(wp) ::   xksilim = 16.5e-6_wp   ! Half-saturation constant for the Si half-saturation constant computation

!! * Substitutions






























   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zint.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_int( kt, Kbb, Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_int  ***
      !!
      !! ** Purpose :   interpolation and computation of various accessory fields
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
      INTEGER, INTENT( in ) ::   Kbb, Kmm ! time level indices
      !
      INTEGER  :: mi, mj, jk              ! dummy loop indices
      REAL(wp) :: zrum, zcodel, zargu, zvar
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p4z_int')
      !
      ! Computation of phyto and zoo metabolic rate
      ! -------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(ts, tgfunc, tgfunc2)
      DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ! Generic temperature dependence (Eppley, 1972)
         tgfunc (mi,mj,jk) = EXP( 0.0631 * t(mi,mj,N+1-jk,nnew,itemp) )
         ! Temperature dependence of mesozooplankton (Buitenhuis et al. (2005))
         tgfunc2(mi,mj,jk) = EXP( 0.0761 * t(mi,mj,N+1-jk,nnew,itemp) )
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO

      IF( ln_p4z .OR. ln_p5z ) THEN
         ! Computation of the silicon dependant half saturation  constant for silica uptake
         ! This is based on an old study by Pondaven et al. (1998)
         ! --------------------------------------------------------------------------------
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, zvar) &
         !$OMP SHARED(tr, xksimax, xksilim)
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zvar = t(mi,mj,N+1-1,Kbb,itemp+ntrc_salt+jpsil) * t(mi,mj,N+1-1,Kbb,itemp+ntrc_salt+jpsil)
            xksimax(mi,mj) = MAX( xksimax(mi,mj), ( 1.+ 7.* zvar / ( xksilim * xksilim + zvar ) ) * 1e-6 )
         END DO   ;   END DO
         !$OMP END PARALLEL DO
         !
         ! At the end of each year, the half saturation constant for silica is 
         ! updated as this is based on the highest concentration reached over 
         ! the year
         ! -------------------------------------------------------------------
         IF( (int(tdays)+1) == year2day ) THEN
            xksi   (:,:) = xksimax(:,:)
            xksimax(:,:) = 0._wp
         ENDIF
      ENDIF
         !
         ! compute the day length depending on latitude and the day
         ! Astronomical parameterization taken from HAMOCC3
      zrum = REAL( (int(tdays)+1) - 80, wp ) / REAL( year2day, wp )
      zcodel = ASIN(  SIN( zrum * pi * 2._wp ) * SIN( pi * 23.5_wp )  )

      ! day length in hours
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, zargu) &
      !$OMP SHARED(zcodel, gphit, pi, strn)
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zargu = TAN( zcodel ) * TAN( latr(mi,mj) * pi )
         zargu = MAX( -1., MIN(  1., zargu ) )
         strn(mi,mj) = MAX( 0.0, 24. - 2. * ACOS( zargu ) / pi / 15. )
      END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(nitrfac, nitrfac2, tr, oxymin, jpoxy, jpno3, Kbb)
      DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
        ! denitrification factor computed from O2 levels
         ! This factor diagnoses below which level of O2 denitrification
         ! is active
         nitrfac(mi,mj,jk) = MAX(  0.e0, 0.4 * ( 6.e-6  - t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpoxy) )  &
            &                                / ( oxymin + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpoxy) )  )
         nitrfac(mi,mj,jk) = MIN( 1., nitrfac(mi,mj,jk) )
         !
         ! redox factor computed from NO3 levels
         ! This factor diagnoses below which level of NO3 additional redox
         ! reactions are taking place.
         nitrfac2(mi,mj,jk) = MAX( 0.e0,       ( 1.E-6 - t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) )  &
            &                                / ( 1.E-6 + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) ) )
         nitrfac2(mi,mj,jk) = MIN( 1., nitrfac2(mi,mj,jk) )
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      IF( .false. )   CALL timing_stop('p4z_int')
      !
   END SUBROUTINE p4z_int


   !!======================================================================
END MODULE p4zint
