










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









MODULE p4zagg
   !!======================================================================
   !!                         ***  MODULE p4zagg  ***
   !! TOP :    aggregation of particles (DOC, POC, GOC)
   !!        This module is the same for both  and -QUOTA
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!             3.5  !  2012-07  (O. Aumont) Introduce potential time-splitting
   !!             3.6  !  2015-05  (O. Aumont)  quota
   !!             3.*  !  2025-02  (S. Maishal, R. Person)  and optimization
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   p4z_agg       :  Compute aggregation of particles
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !   Source Minus Sink variables
   USE prtctl          !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_agg         ! called in p4zbio.F90

   !! * Substitutions





























   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zagg.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_agg ( kt, knt, Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_agg  ***
      !!
      !! ** Purpose :   Compute aggregation of particle. Aggregation by 
      !!                brownian motion, differential settling and shear
      !!                are considered.
      !!
      !! ** Method  : - Aggregation rates are computed assuming a fixed and 
      !!                constant size spectrum in the different particulate 
      !!                pools. The coagulation rates have been computed 
      !!                externally using dedicated programs (O. Aumont). They 
      !!                are hard-coded because they can't be changed 
      !!                independently of each other. 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   !
      INTEGER, INTENT(in) ::   Kbb, Krhs ! time level indices
      !
      INTEGER  ::   mi, mj, jk
      REAL(wp) ::   zagg, zagg1, zagg2, zagg3, zagg4
      REAL(wp) ::   zaggpoc1, zaggpoc2, zaggpoc3, zaggpoc4
      REAL(wp) ::   zaggpoc , zaggfe, zaggdoc, zaggdoc2, zaggdoc3
      REAL(wp) ::   zaggpon , zaggdon, zaggdon2, zaggdon3
      REAL(wp) ::   zaggpop, zaggdop, zaggdop2, zaggdop3
      REAL(wp) ::   zaggtmp, zfact
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p4z_agg')
      !
      !  Exchange between organic matter compartments due to 
      !  coagulation/disaggregation
      !  ---------------------------------------------------

      !  part
      IF( ln_p4z ) THEN
         !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zfact, zagg1, zagg2, zagg3, zagg4, &
                        zagg, zaggfe, zaggdoc, zaggdoc2, zaggdoc3) &
      !$OMP SHARED(xstep, xdiss, tr, tmask, 0.5*EPSILON(1.e0), conspoc, prodgoc)
         DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            !
            zfact = xstep * xdiss(mi,mj,jk)
            ! Part I : Coagulation dependent on turbulence
            !  The stickiness has been assumed to be 0.1
            !  Part I : Coagulation dependent on turbulence
            zagg1 = 12.5  * zfact * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc)
            zagg2 = 169.7 * zfact * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc)

            ! Part II : Differential settling
            ! Aggregation of small into large particles
            ! The stickiness has been assumed to be 0.1
            zagg3 =  8.63  * xstep * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc)
            zagg4 =  132.8 * xstep * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc)

            zagg   = zagg1 + zagg2 + zagg3 + zagg4
            zaggfe = zagg * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsfe) &
                    &   / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) + 0.5*EPSILON(1.e0) )

            ! Aggregation of DOC to POC : 
            ! 1st term is shear aggregation of DOC-DOC
            ! 2nd term is shear aggregation of DOC-POC
            ! 3rd term is differential settling of DOC-POC
            ! 1/3 of DOC is supposed to experience aggregation (HMW)
            zaggdoc  = ( ( 12.0 * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) &
            &            + 9.05 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) ) * zfact       &
            &            + 2.49 * xstep * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) ) &
            &             * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc)
            ! transfer of DOC to GOC : 
            ! 1st term is shear aggregation
            ! 1/3 of DOC is supposed to experience aggregation (HMW)
            zaggdoc2 = ( 1.94 * zfact + 1.37 * xstep ) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc) &
                    &   * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc)
            ! tranfer of DOC to POC due to brownian motion
            ! The temperature dependency has been omitted.
            zaggdoc3 =  ( 127.8 * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) &
                    &  + 725.7 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) ) &
                    &   * xstep * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc)

            !  Update the trends
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) - zagg + zaggdoc + zaggdoc3
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgoc) + zagg + zaggdoc2
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) - zaggfe
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) + zaggfe
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) - zaggdoc - zaggdoc2 - zaggdoc3
            !
            conspoc(mi,mj,jk) = conspoc(mi,mj,jk) - zagg + zaggdoc + zaggdoc3
            prodgoc(mi,mj,jk) = prodgoc(mi,mj,jk) + zagg + zaggdoc2
            !
         END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      ELSE    ! ln_p5z
        ! -QUOTA part
        !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zaggtmp, zfact, zaggpoc1, zaggpoc2, &
                    zaggpoc3, zaggpoc4, zaggpoc, zaggpon, zaggpop, zaggfe, &
                    zaggdoc, zaggdoc2, zaggdoc3, zaggdon, zaggdon2, zaggdon3, &
                    zaggdop, zaggdop2, zaggdop3) &
      !$OMP SHARED(tr, xdiss, xstep, tmask, 0.5*EPSILON(1.e0), conspoc, prodgoc)
         DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            !
            zfact = xstep * xdiss(mi,mj,jk)
            !  Part I : Coagulation dependent on turbulence
            ! The stickiness has been assumed to be 0.1
            zaggtmp = 12.5  * zfact * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc)
            zaggpoc1 = zaggtmp * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc)
            zaggtmp = 169.7 * zfact * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc)
            zaggpoc2 = zaggtmp * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc)
                  
            ! Part II : Differential settling
            ! The stickiness has been assumed to be 0.1
   
            !  Aggregation of small into large particles
            zaggtmp =  8.63  * xstep * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc)
            zaggpoc3 = zaggtmp * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc)
            zaggtmp =  132.8 * xstep * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc)
            zaggpoc4 = zaggtmp * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc)

            zaggpoc = zaggpoc1 + zaggpoc2 + zaggpoc3 + zaggpoc4
            zaggpon = zaggpoc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppon) &
                    &  / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) + 0.5*EPSILON(1.e0))
            zaggpop = zaggpoc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppop) &
                    & / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) + 0.5*EPSILON(1.e0))
            zaggfe  = zaggpoc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsfe) &
                    & / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc)  + 0.5*EPSILON(1.e0) )

            ! Aggregation of DOC to POC : 
            ! 1st term is shear aggregation of DOC-DOC
            ! 2nd term is shear aggregation of DOC-POC
            ! 3rd term is differential settling of DOC-POC
            ! 1/3 of DOC is supposed to experience aggregation (HMW)
            zaggtmp = ( ( 12.0 * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) &
                    &  + 9.05 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) ) * zfact       &
            &            + 2.49 * xstep * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) )
            zaggdoc  = zaggtmp * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc)
            zaggdon  = zaggtmp * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdon)
            zaggdop  = zaggtmp * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdop)

            ! transfer of DOC to GOC : 
            ! 1st term is shear aggregation
            ! 2nd term is differential settling 
            ! 1/3 of DOC is supposed to experience aggregation (HMW)
            zaggtmp = ( 1.94 * zfact + 1.37 * xstep ) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc)
            zaggdoc2 = zaggtmp * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc)
            zaggdon2 = zaggtmp * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdon)
            zaggdop2 = zaggtmp * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdop)

            ! tranfer of DOC to POC due to brownian motion
            ! 1/3 of DOC is supposed to experience aggregation (HMW)
            zaggtmp = ( 127.8 * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) &
                    & +  725.7 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) ) * xstep
            zaggdoc3 =  zaggtmp * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc)
            zaggdon3 =  zaggtmp * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdon)
            zaggdop3 =  zaggtmp * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdop)

            !  Update the trends
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) - zaggpoc + zaggdoc + zaggdoc3
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppon) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppon) - zaggpon + zaggdon + zaggdon3
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppop) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppop) - zaggpop + zaggdop + zaggdop3
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgoc) + zaggpoc + zaggdoc2
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgon) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgon) + zaggpon + zaggdon2
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgop) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgop) + zaggpop + zaggdop2
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) - zaggfe
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) + zaggfe
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) - zaggdoc - zaggdoc2 - zaggdoc3
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdon) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdon) - zaggdon - zaggdon2 - zaggdon3
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdop) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdop) - zaggdop - zaggdop2 - zaggdop3
            !
            conspoc(mi,mj,jk) = conspoc(mi,mj,jk) - zaggpoc + zaggdoc + zaggdoc3
            prodgoc(mi,mj,jk) = prodgoc(mi,mj,jk) + zaggpoc + zaggdoc2
            !
         END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO

         !
      ENDIF
      !
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('agg')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
 !        CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( .false. )   CALL timing_stop('p4z_agg')
      !
   END SUBROUTINE p4z_agg


   !!======================================================================
END MODULE p4zagg
