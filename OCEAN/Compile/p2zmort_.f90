










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









MODULE p2zmort
   !!======================================================================
   !!                         ***  MODULE p2zmort  ***
   !! TOP :    Compute the mortality terms for phytoplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont)  Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.*  !  2025-01  (S. Maishal, R. Person) Change to High Performance
   !!----------------------------------------------------------------------
   !!   p4z_mort       : Compute the mortality terms for phytoplankton
   !!   p4z_mort_init  : Initialize the mortality params for phytoplankton
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE p2zlim          ! Phytoplankton limitation terms
   USE prtctl          ! print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p2z_mort           ! Called from p4zbio.F90 
   PUBLIC   p2z_mort_init      ! Called from trcini_pisces.F90 

   REAL(wp), PUBLIC ::   wchln    !: Quadratic mortality rate of nanophytoplankton
   REAL(wp), PUBLIC ::   mpratn   !: Linear mortality rate of nanophytoplankton

   !! * Substitutions











































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zmort.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p2z_mort( kt, Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_mort_nano  ***
      !!
      !! ** Purpose :   Compute the mortality terms for nanophytoplankton
      !!
      !! ** Method  :   Both quadratic and simili linear mortality terms
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      !!---------------------------------------------------------------------
      !
      INTEGER  ::   mi, mj, jk
      REAL(wp) ::   zcompaph, zprcaca
      REAL(wp) ::   ztortp , zrespp , zmortp, zlim1, zlim2 
      CHARACTER (len=25) ::   charout
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p2z_mort')
      !
      prodcal(:,:,:) = 0._wp   ! calcite production variable set to zero

      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zcompaph, zlim2, zlim1, &
                       zrespp, ztortp, zmortp, zprcaca) &
      !$OMP SHARED(tr, prodcal, xfracal, xlimphy, xdiss, &
                       xkmort, mpratn, wchln, prodpoc)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zcompaph = MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) - 1e-9 ), 0.e0 )

         ! Quadratic mortality of nano due to aggregation during
         ! blooms (Doney et al. 1996)
         ! -----------------------------------------------------
         zlim2   = xlimphy(mi,mj,jk) * xlimphy(mi,mj,jk)
         zlim1   = 0.25 * ( 1. - zlim2 ) / ( 0.25 + zlim2 ) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy)
         zrespp  = wchln * 1.e6 * xstep * zlim1 * xdiss(mi,mj,jk) * zcompaph

         ! Phytoplankton linear mortality
         ! A michaelis-menten like term is introduced to avoid 
         ! extinction of nanophyto in highly limited areas
         ! ----------------------------------------------------
         ztortp = mpratn * xstep * zcompaph / ( xkmort + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) ) &
                 &   * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy)

         zmortp = zrespp + ztortp
         
         !   Update the arrays TRA which contains the biological sources and sinks
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpphy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpphy) - zmortp

         ! Production PIC particles due to mortality
         zprcaca = xfracal(mi,mj,jk) * zmortp
         prodcal(mi,mj,jk) = prodcal(mi,mj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)

         ! POC associated with the shell is supposed to be routed to 
         ! big particles because of the ballasting effect
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) - zprcaca
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) - 2. * zprcaca
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) + zmortp
         prodpoc(mi,mj,jk) = prodpoc(mi,mj,jk) + zmortp
         !
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
       IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nano')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
 !        CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( .false. )   CALL timing_stop('p2z_mort')
      !
   END SUBROUTINE p2z_mort

   SUBROUTINE p2z_mort_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p2z_mort_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton parameters
      !!
      !! ** Method  :   Read the namp4zmort namelist and check the parameters
      !!              called at the first timestep
      !!
      !! ** input   :   Namelist namp2zmort
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp2zmort/ wchln, mpratn
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*) 
         WRITE(stdout,*) 'p2z_mort_init : Initialization of phytoplankton mortality parameters'
         WRITE(stdout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND(numnatp_ref);READ(numnatp_ref,namp2zmort,IOSTAT=ios);CALL ctl_nam(ios,"namp2zmort (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,namp2zmort,IOSTAT=ios);CALL ctl_nam(ios,"namp2zmort (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE( numonp, namp2zmort )
      !
      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*) '   Namelist : namp4zmort'
         WRITE(stdout,*) '      quadratic mortality of phytoplankton        wchln  =', wchln
         WRITE(stdout,*) '      phytoplankton mortality rate                mpratn =', mpratn
      ENDIF
      !
   END SUBROUTINE p2z_mort_init


   !!======================================================================
END MODULE p2zmort
