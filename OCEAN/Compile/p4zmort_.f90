










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









MODULE p4zmort
   !!======================================================================
   !!                         ***  MODULE p4zmort  ***
   !! TOP :    Compute the mortality terms for phytoplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont)  Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec) F90
   !!             3.*  !  2025-02  (S. Maishal, R. Person)  and optimization
   !!----------------------------------------------------------------------
   !!   p4z_mort       : Compute the mortality terms for phytoplankton
   !!   p4z_mort_init  : Initialize the mortality params for phytoplankton
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE p4zprod         ! Primary productivity 
   USE p2zlim          ! Phytoplankton limitation terms
   USE p4zlim          ! Phytoplankton limitation terms
   USE prtctl          ! print control for debugging

   !! * Substitutions












































































   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_mort           ! Called from p4zbio.F90 
   PUBLIC   p4z_mort_init      ! Called from trcini_pisces.F90 

   REAL(wp), PUBLIC ::   wchln    !: Quadratic mortality rate of nanophytoplankton
   REAL(wp), PUBLIC ::   wchld    !: Quadratic mortality rate of diatoms
   REAL(wp), PUBLIC ::   mpratn   !: Linear mortality rate of nanophytoplankton
   REAL(wp), PUBLIC ::   mpratd   !: Linear mortality rate of diatoms

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zmort.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_mort( kt, Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_mort  ***
      !!
      !! ** Purpose :   Calls the different subroutine to compute
      !!                the different phytoplankton mortality terms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      !!---------------------------------------------------------------------
      !
      CALL p4z_mort_nano( Kbb, Krhs )            ! nanophytoplankton
      CALL p4z_mort_diat( Kbb, Krhs )            ! diatoms
      !
   END SUBROUTINE p4z_mort


   SUBROUTINE p4z_mort_nano( Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_mort_nano  ***
      !!
      !! ** Purpose :   Compute the mortality terms for nanophytoplankton
      !!
      !! ** Method  :   Both quadratic and simili linear mortality terms
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Krhs  ! time level indices
      INTEGER  ::   mi, mj, jk
      REAL(wp) ::   zcompaph
      REAL(wp) ::   zfactfe, zfactch, zprcaca, zfracal
      REAL(wp) ::   ztortp , zrespp , zmortp, zlim1, zlim2 
      CHARACTER (len=25) ::   charout
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p4z_mort_nano')
      !
      prodcal(:,:,:) = 0._wp ! calcite production variable set to zero
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(tr, xlimphy, wchln, xstep, xdiss, mpratn, xkmort, &
                   0.5*EPSILON(1.e0), xfracal, prodcal, prodpoc, zp, zfactfe, zfactch)
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
             &    * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy)

         zmortp = zrespp + ztortp
         
         ! Update the arrays TRA which contains the biological sources and sinks
         zfactfe = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnfe)/(t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy)+0.5*EPSILON(1.e0))
         zfactch = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnch)/(t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy)+0.5*EPSILON(1.e0))
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpphy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpphy) - zmortp
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnch) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnch) - zmortp * zfactch
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnfe) - zmortp * zfactfe

         ! Production PIC particles due to mortality
         zprcaca = xfracal(mi,mj,jk) * zmortp
         ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
         prodcal(mi,mj,jk) = prodcal(mi,mj,jk) + zprcaca

         ! POC associated with the shell is supposed to be routed to 
         ! big particles because of the ballasting effect
         zfracal = 0.5 * xfracal(mi,mj,jk)
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) - zprcaca
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) - 2. * zprcaca
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpcal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpcal) + zprcaca
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) + zmortp
         prodpoc(mi,mj,jk) = prodpoc(mi,mj,jk) + zmortp

         ! Update the arrays TRA which contains the biological sources and sinks
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) + zmortp * zfactfe
         !
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
       IF(sn_cfctl%l_prttrc) THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nano')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
      !  CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( .false. ) CALL timing_stop('p4z_mort_nano')
      !
   END SUBROUTINE p4z_mort_nano


   SUBROUTINE p4z_mort_diat( Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_mort_diat  ***
      !!
      !! ** Purpose :   Compute the mortality terms for diatoms
      !!
      !! ** Method  : - Both quadratic and simili linear mortality terms
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: Kbb, Krhs ! time level indices
      INTEGER  :: mi, mj, jk
      REAL(wp) :: zfactfe,zfactsi,zfactch, zcompadi
      REAL(wp) :: zrespp2, ztortp2, zmortp2
      REAL(wp) :: zlim2, zlim1
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p4z_mort_diat')
      !
      ! Aggregation term for diatoms is increased in case of nutrient
      ! stress as observed in reality. The stressed cells become more
      ! sticky and coagulate to sink quickly out of the euphotic zone
      ! This is due to the production of EPS by stressed cells
      ! -------------------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(tr, xlimdia, wchld, xstep, xdiss, mpratd, 0.5*EPSILON(1.e0), &
                   prodpoc, prodgoc, zfactch, zfactfe, zfactsi)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp

         zcompadi = MAX( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) - 1e-9), 0. )

         ! Aggregation term for diatoms is increased in case of nutrient
         ! stress as observed in reality. The stressed cells become more
         ! sticky and coagulate to sink quickly out of the euphotic zone
         ! ------------------------------------------------------------
         zlim2   = xlimdia(mi,mj,jk) * xlimdia(mi,mj,jk)
         zlim1   = 0.25 * ( 1. - zlim2 ) / ( 0.25 + zlim2 ) 
         zrespp2 = 1.e6 * xstep * wchld * zlim1 * xdiss(mi,mj,jk) &
                 * zcompadi * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia)

         ! Phytoplankton linear mortality
         ! A michaelis-menten like term is introduced to avoid 
         ! extinction of diatoms in highly limited areas
         !  ---------------------------------------------------
         ztortp2 = mpratd * xstep * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) &
            &     / ( xkmort + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) ) * zcompadi

         zmortp2 = zrespp2 + ztortp2

         !   Update the arrays trends which contains the biological sources and sinks
         !   ---------------------------------------------------------------------
         zfactch = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdch) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) + 0.5*EPSILON(1.e0) )
         zfactfe = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdfe) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) + 0.5*EPSILON(1.e0) )
         zfactsi = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdsi) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) + 0.5*EPSILON(1.e0) )
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdia) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdia) - zmortp2
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdch) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdch) - zmortp2 * zfactch
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdfe) - zmortp2 * zfactfe
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdsi) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdsi) - zmortp2 * zfactsi
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgsi) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgsi) + zmortp2 * zfactsi

         ! Half of the linear mortality term is routed to big particles
         ! becaue of the ballasting effect
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgoc) + zrespp2
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) + ztortp2
         prodpoc(mi,mj,jk) = prodpoc(mi,mj,jk) + ztortp2
         prodgoc(mi,mj,jk) = prodgoc(mi,mj,jk) + zrespp2
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) + ztortp2 * zfactfe
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) + zrespp2 * zfactfe
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      IF(sn_cfctl%l_prttrc) THEN      ! print mean trends (used for debugging)
         WRITE(charout, FMT="('diat')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
   !      CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( .false. )   CALL timing_stop('p4z_mort_diat')
      !
   END SUBROUTINE p4z_mort_diat


   SUBROUTINE p4z_mort_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_mort_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton parameters
      !!
      !! ** Method  :   Read the namp4zmort namelist and check the parameters
      !!              called at the first timestep
      !!
      !! ** input   :   Namelist namp4zmort
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ios ! Local integer
      !
      NAMELIST/namp4zmort/ wchln, wchld, mpratn, mpratd
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*) 
         WRITE(stdout,*) 'p4z_mort_init : Initialization of phytoplankton mortality parameters'
         WRITE(stdout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND(numnatp_ref);READ(numnatp_ref,namp4zmort,IOSTAT=ios);CALL ctl_nam(ios,"namp4zmort (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,namp4zmort,IOSTAT=ios);CALL ctl_nam(ios,"namp4zmort (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE( numonp, namp4zmort )
      !
      IF(mynode .eq. 0) THEN              ! control print
         WRITE(stdout,*) '      Namelist : namp4zmort'
         WRITE(stdout,*) '      quadratic mortality of phytoplankton        wchln  =', wchln
         WRITE(stdout,*) '      maximum quadratic mortality of diatoms      wchld  =', wchld
         WRITE(stdout,*) '      phytoplankton mortality rate                mpratn =', mpratn
         WRITE(stdout,*) '      Diatoms mortality rate                      mpratd =', mpratd
      ENDIF
      !
   END SUBROUTINE p4z_mort_init


   !!======================================================================
END MODULE p4zmort
