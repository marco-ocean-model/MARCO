










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









MODULE p4zpoc
   !!======================================================================
   !!                         ***  MODULE p4zpoc  ***
   !! TOP :    Compute remineralization of organic particles
   !!         Same module for both  and -QUOTA
   !!=========================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!             3.6  !  2016-03  (O. Aumont) Quota model and diverse
   !!             4.0  !  2018     (O. Aumont) Variable lability parameterization
   !!----------------------------------------------------------------------
   !!   p4z_poc       :  Compute remineralization/dissolution of organic compounds
   !!   p4z_poc_init  :  Initialisation of parameters for remineralisation
   !!   alngam and gamain : computation of the incomplete gamma function
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !   Source Minus Sink variables
   USE prtctl          !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_poc         ! called in p4zbio.F90
   PUBLIC   p4z_poc_init    ! called in trcini_pisces.F90

   REAL(wp), PUBLIC ::   xremipc    !: remineralisation rate of DOC
   REAL(wp), PUBLIC ::   xremipn    !: remineralisation rate of DON
   REAL(wp), PUBLIC ::   xremipp    !: remineralisation rate of DOP
   INTEGER , PUBLIC ::   jcpoc      !: number of lability classes
   REAL(wp), PUBLIC ::   rshape     !: shape factor of the gamma distribution

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)       ::  alphan, reminp   !: variable lability of POC and initial distribution
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::  alphap, alphag    !: lability distribution of small particles
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::  orem3    !: lability distribution of small particles
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::  pdep

   REAL(wp ) :: solgoc, rbound
   INTEGER   :: ndayflx 
   LOGICAL   :: l_dia_remin

   LOGICAL, PUBLIC  :: ll_poc_lab

   !! * Substitutions











































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zpoc.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_poc( kt, knt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_poc  ***
      !!
      !! ** Purpose :   Compute remineralization of organic particles
      !!                A reactivity-continuum parameterization is chosen
      !!                to describe the lability of the organic particles
      !!                As a consequence, the remineralisation rates of the 
      !!                the different pools change with time as a function of 
      !!                the lability distribution
      !!
      !! ** Method  : - Computation of the remineralisation rates is performed
      !!                according to reactivity continuum formalism described
      !!                in Aumont et al. (2017). 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt         ! ocean time step and ???
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      !
      INTEGER  ::   mi, mj, jk, jn
      REAL(wp) ::   zremip, zremig, zorem, zorem2, zofer, zfact
      REAL(wp) ::   zopon, zopop, zopon2, zopop2
      REAL(wp) ::   zofer2, zofer3, zreminp1, zreminp2
      CHARACTER (len=25) :: charout
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zremipoc, zremigoc, zfolimi, zw3d
      !!---------------------------------------------------------------------
      !
      IF( .false. )  CALL timing_start('p4z_poc')
      !
      IF( kt == nittrc000 .AND. knt == 1 ) THEN
         l_dia_remin = iom_use( "REMINP" ) .OR. iom_use( "REMING" ) .OR. iom_use( "REMINF" )
         ALLOCATE( pdep(Istrp:Iendp,Jstrp:Jendp) )
      ENDIF
      !
      IF( l_dia_remin ) THEN
         ALLOCATE( zfolimi (Istrp:Iendp,Jstrp:Jendp,N) )  ;  zfolimi (Istrp:Iendp,Jstrp:Jendp,N) = 0._wp
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zfolimi (mi,mj,jk) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer)
         END DO   ;   END DO   ;   END DO
      ENDIF
      ! Initialisation of arrays
      orem(:,:,N) = 0._wp
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         pdep(mi,mj) = MAX( hbl(mi,mj), ((-1)*(z_r(mi,mj,N+1-1)-z_w(mi,mj,N))) )
      END DO   ;   END DO
      !
      IF( ndayflx /= (int(tdays)+1) ) THEN   ! New day
         ll_poc_lab  = .TRUE.     
         ndayflx = (int(tdays)+1)
      ELSE
         ll_poc_lab = .FALSE.
      ENDIF    
      !
      IF( kt == nittrc000 .AND. .NOT. ln_rsttr )   ll_poc_lab = .TRUE.
      !
      ll_poc_lab = ll_poc_lab .AND. knt == 1

      IF( .NOT. ln_p2z ) THEN
         !
         IF( ll_poc_lab ) THEN
            ! Initialisation of the lability distributions that are set to the distribution of newly produced organic particles
            IF(mynode .eq. 0) write(stdout,*)
            IF(mynode .eq. 0) write(stdout,*) ' Compute variable lability for GOC at kt =  ',  kt, '  day = ', (int(tdays)+1)
            IF(mynode .eq. 0) write(stdout,*) '~~~~~~'
            !
            CALL p4z_goc_lab( kt, Kbb, Kmm )
            !
         ENDIF
         ! Lability parameterization. This is the big particles part (GOC)
         ! ----------------------------------------------------------------- 
         IF( l_dia_remin ) THEN
           ALLOCATE( zremigoc(Istrp:Iendp,Jstrp:Jendp,N) )  ;  zremigoc(Istrp:Iendp,Jstrp:Jendp,N) = 0._wp
           DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
              zremigoc(mi,mj,jk) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc)
           END DO   ;   END DO   ;   END DO
         ENDIF
         !
         ! The standard  part
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            ! POC degradation by bacterial activity. It is a function
            ! of the mean lability and of temperature. This also includes
            ! shrinking of particles due to the bacterial activity
            ! -----------------------------------------------------------
            zremig  = remintgoc(mi,mj,jk) * xstep * tgfunc(mi,mj,jk)
            zorem2  = zremig * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc)
            orem(mi,mj,jk)  = zorem2
            orem3(mi,mj,jk) = zremig * solgoc * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc)
            zofer2 = zremig * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpbfe)

            ! update of the TRA arrays
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) + orem3(mi,mj,jk)
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgoc) - zorem2 - orem3(mi,mj,jk)
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) + solgoc * zofer2
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) - (1. + solgoc) * zofer2
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) + zorem2
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) + zofer2
         END DO   ;   END DO   ;   END DO
         !
         IF ( ln_p5z ) THEN
            DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               ! POC degradation by bacterial activity. It is a function
               ! of the mean lability and of temperature. This also includes
               ! shrinking of particles due to the bacterial activity
               ! --------------------------------------------------------
               zremig = remintgoc(mi,mj,jk) * xstep * tgfunc(mi,mj,jk)
               zopon2 = xremipn / xremipc * zremig * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgon)
               zopop2 = xremipp / xremipc * zremig * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgop)

               ! update of the TRA arrays
               t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppon) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppon) + solgoc * zopon2
               t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppop) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppop) + solgoc * zopop2
               t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdon) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdon) + zopon2
               t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdop) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdop) + zopop2
               t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgon) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgon) - zopon2 * (1. + solgoc)
               t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgop) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgop) - zopop2 * (1. + solgoc)
            END DO   ;   END DO   ;   END DO
         ENDIF
         ! Remineralisation rate of large particles diag.
         IF( l_dia_remin ) THEN
           DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
              zremigoc(mi,mj,jk) = ( t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) - zremigoc(mi,mj,jk) )  &
              &       / ( xstep * tgfunc(mi,mj,jk) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc) + 0.5*EPSILON(1.e0) ) &
              &             *  tmask(mi,mj,jk) ! =zremipart
           END DO   ;   END DO   ;   END DO
         ENDIF
         !
         IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
            WRITE(charout, FMT="('poc1')")
            CALL prt_ctl_info( charout, cdcomp = 'top' )
!            CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
         ENDIF
         !
      ENDIF
      !
      ! Lability parameterization for the small OM particles.
      ! ---------------------------------------------------------
      IF( ll_poc_lab ) THEN
         ! Initialisation of the lability distributions that are set to the distribution of newly produced organic particles
         IF(mynode .eq. 0) write(stdout,*)
         IF(mynode .eq. 0) write(stdout,*) ' Compute variable lability for POC at kt =  ',  kt, '  day = ', (int(tdays)+1)
         IF(mynode .eq. 0) write(stdout,*) '~~~~~~'
         !
         CALL p4z_poc_lab( kt, Kbb, Kmm )
         !
      ENDIF
      !
      IF( l_dia_remin ) THEN
         ALLOCATE( zremipoc(Istrp:Iendp,Jstrp:Jendp,N) )  ;  zremipoc(Istrp:Iendp,Jstrp:Jendp,N) = 0._wp
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zremipoc(mi,mj,jk) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc)
         END DO   ;   END DO   ;   END DO
      ENDIF
      !
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
        ! POC disaggregation by turbulence and bacterial activity.It is a function
        ! of the mean lability and of temperature  
        ! --------------------------------------------------------
        zremip = remintpoc(mi,mj,jk) * xstep * tgfunc(mi,mj,jk)
        zorem  = zremip * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc)
              
        ! Update of the TRA arrays
        t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) + zorem
        t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) - zorem
     END DO   ;   END DO   ;   END DO
     IF( ln_p2z ) THEN
        DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           ! POC disaggregation by turbulence and bacterial activity.It is a function
           ! of the mean lability and of temperature  
           ! --------------------------------------------------------
           zremip = remintpoc(mi,mj,jk) * xstep * tgfunc(mi,mj,jk)
           zorem  = zremip * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc)
           orem(mi,mj,jk)  =  zorem
           ! Update of the TRA arrays
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) + zorem * feratz
        END DO   ;   END DO   ;   END DO
     ELSE
        DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           ! POC disaggregation by turbulence and bacterial activity.It is a function
           ! of the mean lability and of temperature  
           ! --------------------------------------------------------
           zremip = remintpoc(mi,mj,jk) * xstep * tgfunc(mi,mj,jk)
           zorem  = zremip * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc)
           orem(mi,mj,jk)  = orem(mi,mj,jk) + zorem
           zofer  = zremip * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsfe)

           ! Update of the TRA arrays
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) + zofer
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) - zofer
        END DO   ;   END DO   ;   END DO
     ENDIF
     IF ( ln_p5z ) THEN
        DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           ! POC disaggregation by turbulence and bacterial activity.It is a function
           ! of the mean lability and of temperature  
           !--------------------------------------------------------
           zremip = remintpoc(mi,mj,jk) * xstep * tgfunc(mi,mj,jk)
           zopon  = xremipn / xremipc * zremip * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppon)
           zopop  = xremipp / xremipc * zremip * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppop)
           
           ! Update of the TRA arrays
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppon) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppon) - zopon
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppop) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppop) - zopop
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdon) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdon) + zopon 
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdop) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdop) + zopop 
        END DO   ;   END DO   ;   END DO
     ENDIF
     !
     IF( l_dia_remin ) THEN
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zremipoc(mi,mj,jk) = ( t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) - zremipoc(mi,mj,jk) ) &
               &     / ( xstep * tgfunc(mi,mj,jk) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) + 0.5*EPSILON(1.e0) ) &
               &          * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zfolimi (mi,mj,jk) = ( t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) - zfolimi (mi,mj,jk) ) * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
     ENDIF

     IF( .false. .AND. l_dia_remin .AND. knt == nrdttrc ) THEN
        ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
        IF( .NOT. ln_p2z ) THEN
           DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
              zw3d(mi,mj,N+1-jk) = zremigoc(mi,mj,jk)
           END DO   ;   END DO   ;   END DO
           CALL iom_put( "REMING", zw3d ) ! Remineralisation rate of large particles
           DEALLOCATE ( zremigoc )
        ENDIF
        DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           zw3d(mi,mj,N+1-jk) = zremipoc(mi,mj,jk)
        END DO   ;   END DO   ;   END DO
        CALL iom_put( "REMINP", zw3d )  ! Remineralisation rate of small particles
        !
        DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           zw3d(mi,mj,N+1-jk) = zfolimi(mi,mj,jk) * 1.e+9 * 1.e3 * rfact2r
        END DO   ;   END DO   ;   END DO
        CALL iom_put( "REMINF", zw3d ) ! Remineralisation of biogenic particulate iron
        DEALLOCATE ( zremipoc, zfolimi )
        DEALLOCATE ( zw3d )
     ENDIF

      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('poc2')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
  !       CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      !
      IF( .false. )   CALL timing_stop('p4z_poc')
      !
   END SUBROUTINE p4z_poc

   SUBROUTINE p4z_goc_lab( kt, Kbb, Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_poc  ***
      !!
      !! ** Purpose :   Compute remineralization of organic particles
      !!                A reactivity-continuum parameterization is chosen
      !!                to describe the lability of the organic particles
      !!                As a consequence, the remineralisation rates of the 
      !!                the different pools change with time as a function of 
      !!                the lability distribution
      !!
      !! ** Method  : - Computation of the remineralisation rates is performed
      !!                according to reactivity continuum formalism described
      !!                in Aumont et al. (2017). 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt         ! ocean time step and ???
      INTEGER, INTENT(in) ::   Kbb, Kmm  ! time level indices
      !
      INTEGER  ::   mi, mj, jk, jn
      REAL(wp) ::   zsizek, alphat, zremint, alphatm1
      REAL(wp) ::   zpoc1, zpoc2, zfact, zconcpoc
      REAL(wp) ::   zreminp1, zreminp2
      REAL(wp) ::   ztemp, ztemp1, ztemp2
      REAL(wp), DIMENSION(jcpoc) :: alpham1

      DO jn = 1, jcpoc
         alphag(:,:,:,jn) = alphan(jn)
      END DO

     ! Lability parameterization. This is the big particles part (GOC)
     ! This lability parameterization is always active. However, if only one
     ! lability class is specified in the namelist, this is equivalent to 
     ! a standard parameterisation with a constant lability
     ! -----------------------------------------------------------------------
     remintgoc(:,:,:) = xremipc
     DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
        alpham1(:) = alphan(:)
        DO jk = 2, N-1
           IF( tmask(mi,mj,jk) == 1. .AND. ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) > pdep(mi,mj) ) THEN
              !
              ! In the case of GOC, lability is constant in the mixed layer 
              ! It is computed only below the mixed layer depth
              ! ------------------------------------------------------------
              zsizek = Hz(mi,mj,N+1-jk) / 2. / (wsbio4(mi,mj,jk) + 0.5*EPSILON(1.e0))
              !
              ! standard algorithm in the rest of the water column
              ! See the comments in the previous block.
              ! ---------------------------------------------------
              zfact  = day2sec / rfact2 / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc) + 0.5*EPSILON(1.e0) )
              zpoc1  = MIN(rbound, MAX(-rbound, consgoc(mi,mj,jk) * zfact ) )
              zpoc2  = prodgoc(mi,mj,jk) * day2sec / rfact2
              !
              zconcpoc = ( Hz(mi,mj,N+1-jk-1) * t(mi,mj,N+1-jk-1,Kbb,itemp+ntrc_salt+jpgoc) &
                      &  + Hz(mi,mj,N+1-jk) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc) )   &
                      &        / ( Hz(mi,mj,N+1-jk-1) + Hz(mi,mj,N+1-jk) )
              !
              DO jn = 1, jcpoc
                 zreminp1 = reminp(jn) * tgfunc(mi,mj,jk) - zpoc1
                 ztemp    = MIN(rbound, MAX(-rbound,  zreminp1 ) )
                 ztemp1   = EXP( -ztemp * zsizek)
                 ztemp2   = zpoc2 * ( 1. - ztemp1 ) / ztemp * alphan(jn)
                 alphag(mi,mj,jk,jn) = alpham1(jn) * ztemp1 * zconcpoc + ztemp2
                 alpham1(jn) = alphag(mi,mj,jk,jn) * ztemp1 + ztemp2
              END DO

              alphatm1 = SUM( alpham1(:) ) + 0.5*EPSILON(1.e0)
              alphat   = SUM( alphag(mi,mj,jk,:) ) + 0.5*EPSILON(1.e0)
              alphag(mi,mj,jk,:) = alphag(mi,mj,jk,:) / alphat
              alpham1(:)         = alpham1(:) / alphatm1
              !
              ! The contribution of each lability class at the current level is computed
              zremint = SUM( alphag(mi,mj,jk,:) * reminp(:) )
              ! Computation of the mean remineralisation rate
              remintgoc(mi,mj,jk) = MIN( xremipc, zremint )
              !
           ENDIF
        END DO
     END DO   ;   END DO
     !
   END SUBROUTINE p4z_goc_lab

   SUBROUTINE p4z_poc_lab( kt, Kbb, Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_poc  ***
      !!
      !! ** Purpose :   Compute remineralization of organic particles
      !!                A reactivity-continuum parameterization is chosen
      !!                to describe the lability of the organic particles
      !!                As a consequence, the remineralisation rates of the 
      !!                the different pools change with time as a function of 
      !!                the lability distribution
      !!
      !! ** Method  : - Computation of the remineralisation rates is performed
      !!                according to reactivity continuum formalism described
      !!                in Aumont et al. (2017). 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt         ! ocean time step and ???
      INTEGER, INTENT(in) ::   Kbb, Kmm  ! time level indices
      !
      INTEGER  ::   mi, mj, jk, jn
      REAL(wp) ::   zsizek, alphat, zremint, alphatm1
      REAL(wp) ::   zpoc1, zpoc2, zpoc3, zfact, zconcpoc
      REAL(wp) ::   zreminp1, zreminp2
      REAL(wp) ::   ztemp, ztemp1, ztemp2, ztemp3
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp  )   :: ztotprod, ztotthick, ztotcons
      REAL(wp), DIMENSION(jcpoc) :: alpham1

      DO jn = 1, jcpoc
         alphap(:,:,:,jn) = alphan(jn)
      END DO

     ! Lability parameterization for the small OM particles. This param
     ! is based on the same theoretical background as the big particles.
     ! However, because of its low sinking speed, lability is not supposed
     ! to be equal to its initial value (the value of the freshly produced
     ! organic matte) in the MLD. It is however uniform in the mixed layer.
     ! ---------------------------------------------------------------------
     ztotprod (:,:)    = 0.
     ztotthick(:,:)    = 0.
     ztotcons (:,:)    = 0.
     remintpoc(:,:,:) = xremipc

     ! intregrated production and consumption of POC in the mixed layer
     ! ----------------------------------------------------------------
     DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
        IF (tmask(mi,mj,jk) == 1. .AND. ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) <= pdep(mi,mj) ) THEN
          zfact = Hz(mi,mj,N+1-jk) * day2sec / rfact2
          ztotprod(mi,mj)  = ztotprod(mi,mj) + prodpoc(mi,mj,jk) * zfact
          ! The temperature effect is included here
          ztotthick(mi,mj) = ztotthick(mi,mj) + Hz(mi,mj,N+1-jk) * tgfunc(mi,mj,jk)
          ztotcons(mi,mj)  = ztotcons(mi,mj) - conspoc(mi,mj,jk) &
                  &    * zfact / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) + 0.5*EPSILON(1.e0) )
        ENDIF
     END DO   ;   END DO   ;   END DO

     ! Computation of the lability spectrum in the mixed layer. In the mixed
     ! layer, this spectrum is supposed to be uniform as a result of intense
     ! mixing.
     ! ---------------------------------------------------------------------
     DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
        IF (tmask(mi,mj,jk) == 1. .AND. ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) <= pdep(mi,mj) ) THEN
           DO jn = 1, jcpoc
             ! For each lability class, the system is supposed to be
             ! at equilibrium: Prod - Sink - w alphap = 0.
             alphap(mi,mj,jk,jn) = ztotprod(mi,mj) * alphan(jn) / ( reminp(jn)    &
             &                     * ztotthick(mi,mj) + ztotcons(mi,mj) + wsbio3(mi,mj,jk) + 0.5*EPSILON(1.e0) )
             alphap(mi,mj,jk,jn) = MAX(0., alphap(mi,mj,jk,jn) )
          END DO
          alphat = SUM( alphap(mi,mj,jk,:) ) + 0.5*EPSILON(1.e0)
          alphap(mi,mj,jk,:)  = alphap(mi,mj,jk,:) / alphat
          zremint = SUM( alphap(mi,mj,jk,:) * reminp(:) )
          ! Mean remineralization rate in the mixed layer
          remintpoc(mi,mj,jk) =  MIN( xremipc, zremint )
        ENDIF
     END DO   ;   END DO   ;   END DO
     !
     !
     ! The lability parameterization is used here. The code is here
     ! almost identical to what is done for big particles. The only difference
     ! is that an additional source from GOC to POC is included. This means
     ! that since we need the lability spectrum of GOC, GOC spectrum
     ! should be determined before.
     ! -----------------------------------------------------------------------
     DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
        alpham1(:) = alphap(mi,mj,1,:)
        DO jk = 2, N-1
           IF (tmask(mi,mj,jk) == 1. .AND. ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) > pdep(mi,mj) ) THEN
              ! the scale factors are corrected with temperature
              zsizek  = Hz(mi,mj,N+1-jk) / 2. / (wsbio3(mi,mj,jk) + 0.5*EPSILON(1.e0))
              !
              ! Special treatment of the level just below the MXL
              ! See the comments in the GOC section
              zfact  = day2sec / rfact2 / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) + 0.5*EPSILON(1.e0) )
              zpoc1  = MIN(rbound, MAX(-rbound, conspoc(mi,mj,jk) * zfact ) )
              zpoc2  = prodpoc(mi,mj,jk) * day2sec / rfact2
              zpoc3  = orem3 (mi,mj,jk) * day2sec / rfact2
              !
              zconcpoc = ( Hz(mi,mj,N+1-jk-1) * t(mi,mj,N+1-jk-1,Kbb,itemp+ntrc_salt+jppoc) &
                &        + Hz(mi,mj,N+1-jk) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) )   &
                &        / ( Hz(mi,mj,N+1-jk-1) + Hz(mi,mj,N+1-jk) )
              !
              DO jn = 1, jcpoc
                 zreminp1 = reminp(jn) * tgfunc(mi,mj,jk) - zpoc1
                 ztemp    = MIN(rbound, MAX(-rbound,  zreminp1 ) )
                 ztemp1   = EXP( MIN(rbound,-ztemp * zsizek ) )
                 ztemp2   = zpoc2 * ( 1. - ztemp1 ) / ztemp * alphan(jn)
                 ztemp3   = zpoc3 * ( 1. - ztemp1 ) / ztemp * alphag(mi,mj,jk,jn)
                 alphap(mi,mj,jk,jn) = alpham1(jn) * ztemp1 * zconcpoc + ztemp2 + ztemp3
                 alpham1(jn) = alphap(mi,mj,jk,jn) * ztemp1 + ztemp2 + ztemp3
              END DO
              !
              alphat = SUM( alphap(mi,mj,jk,:) ) + 0.5*EPSILON(1.e0)
              alphatm1 = SUM( alpham1(:) ) + 0.5*EPSILON(1.e0)
              alphap(mi,mj,jk,:) = alphap(mi,mj,jk,:) / alphat
              alpham1(:) = alpham1(:) / alphatm1
              ! The contribution of each lability class at the current level is computed
              zremint = SUM( alphap(mi,mj,jk,:) * reminp(:) )
              ! Computation of the mean remineralisation rate
              remintpoc(mi,mj,jk) =  MIN( xremipc, zremint )
           ENDIF
        END DO
     END DO   ;   END DO
     !
   END SUBROUTINE p4z_poc_lab


   SUBROUTINE p4z_poc_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_poc_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampispoc namelist and check the parameters
      !!              called at the first timestep
      !!
      !! ** input   :   Namelist nampispoc
      !!----------------------------------------------------------------------
      INTEGER ::   jn            ! dummy loop index
      INTEGER ::   ios           ! Local integer
      REAL(wp)::   zremindelta, zreminup, zremindown
      REAL(wp)::   zup, zup1, zdown, zdown1
      !!
      NAMELIST/nampispoc/ jcpoc  , rshape,  &
         &                xremipc, xremipn, xremipp
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*)
         WRITE(stdout,*) 'p4z_poc_init : Initialization of remineralization parameters'
         WRITE(stdout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND(numnatp_ref);READ(numnatp_ref,nampispoc,IOSTAT=ios);CALL ctl_nam(ios,"nampispoc (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,nampispoc,IOSTAT=ios);CALL ctl_nam(ios,"nampispoc (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE( numonp, nampispoc )

      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*) '   Namelist : nampispoc'
         WRITE(stdout,*) '      remineralisation rate of POC              xremipc   =', xremipc
         IF( ln_p5z ) THEN 
            WRITE(stdout,*) '      remineralisation rate of PON              xremipn   =', xremipn
            WRITE(stdout,*) '      remineralisation rate of POP              xremipp   =', xremipp
         ENDIF
         WRITE(stdout,*) '      Number of lability classes for POC        jcpoc     =', jcpoc
         WRITE(stdout,*) '      Shape factor of the gamma distribution    rshape    =', rshape
      ENDIF
      !
      ! Discretization along the lability space
      ! ---------------------------------------
      !
     !
      ALLOCATE( alphan(jcpoc) , reminp(jcpoc) )
      ALLOCATE( alphap(Istrp:Iendp,Jstrp:Jendp,N,jcpoc), alphag(Istrp:Iendp,Jstrp:Jendp,N,jcpoc) )
      ALLOCATE( orem3(Istrp:Iendp,Jstrp:Jendp,N) )
      !
      IF (jcpoc > 1) THEN  ! Case when more than one lability class is used
         !
         zremindelta = LOG(4. * 1000. ) / REAL(jcpoc-1, wp)
         zreminup = 1./ 400. * EXP(zremindelta)
         !
         ! Discretization based on incomplete gamma functions
         ! As incomplete gamma functions are not available in standard 
         ! fortran 95, they have been coded as functions in this module (gamain)
         ! ---------------------------------------------------------------------
         !
         CALL gamain(zreminup, rshape, zup )
         CALL gamain(zreminup, rshape+1.0, zup1 )
         alphan(1) = zup
         reminp(1) =  zup1 * xremipc / alphan(1)
         DO jn = 2, jcpoc-1
            zreminup = 1./ 400. * EXP( REAL(jn, wp) * zremindelta)
            zremindown = 1. / 400. * EXP( REAL(jn-1, wp) * zremindelta)
            CALL gamain(zreminup, rshape, zup )
            CALL gamain(zremindown, rshape, zdown )
            alphan(jn) = zup - zdown
            CALL gamain(zreminup, rshape+1.0, zup1 )
            CALL gamain(zremindown, rshape+1.0, zdown1 )
            reminp(jn) = zup1 -zdown1
            reminp(jn) = reminp(jn) * xremipc / alphan(jn) 
         END DO
         zremindown = 1. / 400. * EXP( REAL(jcpoc-1, wp) * zremindelta)
         CALL gamain(zremindown, rshape, zdown )
         CALL gamain(zremindown, rshape+1.0, zdown1 )
         alphan(jcpoc) = 1.0 - zdown
         reminp(jcpoc) = 1.0 - zdown1
         reminp(jcpoc) = reminp(jcpoc) * xremipc / alphan(jcpoc)

      ELSE  ! Only one lability class is used
         alphan(jcpoc) = 1.
         reminp(jcpoc) = xremipc
      ENDIF

      DO jn = 1, jcpoc
         alphap(:,:,:,jn) = alphan(jn)
         alphag(:,:,:,jn) = alphan(jn)
      END DO

      ! Here we compute the GOC -> POC rate due to the shrinking
      ! of the fecal pellets/aggregates as a result of bacterial
      ! solubilization
      ! This is based on a fractal dimension of 2.56 and a spectral
      ! slope of -3.6 (identical to what is used in p4zsink to compute
      ! aggregation
      solgoc = 0.04/ 2.56 * 1./ ( 1.-50**(-0.04) )
      !
      rbound = 1.e+01_wp
      !
      orem3(:,:,:) = 0.
      !
      ndayflx = (int(tdays)+1)  ! Initialize a counter of the current day
      !

   END SUBROUTINE p4z_poc_init


   SUBROUTINE alngam( xvalue, xres )
      !*****************************************************************************80
      !
      !! ALNGAM computes the logarithm of the gamma function.
      !
      !  Modified:    13 January 2008
      !
      !  Author  :    Allan Macleod
      !               FORTRAN90 version by John Burkardt
      !
      !  Reference:
      !    Allan Macleod, Algorithm AS 245,
      !    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
      !    Applied Statistics,
      !    Volume 38, Number 2, 1989, pages 397-402.
      !
      !  Parameters:
      !
      !    Input, real ( kind = 8 ) XVALUE, the argument of the Gamma function.
      !
      !    0, no error occurred.
      !    1, XVALUE is less than or equal to 0.
      !    2, XVALUE is too big.
      !
      !    Output, real ( kind = 8 ) ALNGAM, the logarithm of the gamma function of X.
      !*****************************************************************************80
  real(wp), parameter :: alr2pi = 0.918938533204673E+00
  real(wp), parameter, dimension ( 9 ) :: r1 = (/ &
    -2.66685511495E+00, &
    -24.4387534237E+00, &
    -21.9698958928E+00, &
     11.1667541262E+00, &
     3.13060547623E+00, &
     0.607771387771E+00, &
     11.9400905721E+00, &
     31.4690115749E+00, &
     15.2346874070E+00 /)
  real(wp), parameter, dimension ( 9 ) :: r2 = (/ &
    -78.3359299449E+00, &
    -142.046296688E+00, &
     137.519416416E+00, &
     78.6994924154E+00, &
     4.16438922228E+00, &
     47.0668766060E+00, &
     313.399215894E+00, &
     263.505074721E+00, &
     43.3400022514E+00 /)
  real(wp), parameter, dimension ( 9 ) :: r3 = (/ &
    -2.12159572323E+05, &
     2.30661510616E+05, &
     2.74647644705E+04, &
    -4.02621119975E+04, &
    -2.29660729780E+03, &
    -1.16328495004E+05, &
    -1.46025937511E+05, &
    -2.42357409629E+04, &
    -5.70691009324E+02 /)
  real(wp), parameter, dimension ( 5 ) :: r4 = (/ &
     0.279195317918525E+00, &
     0.4917317610505968E+00, &
     0.0692910599291889E+00, &
     3.350343815022304E+00, &
     6.012459259764103E+00 /)
  real (wp) :: x
  real (wp) :: x1
  real (wp) :: x2
  real (wp), parameter :: xlge = 5.10E+05
  real (wp), parameter :: xlgst = 1.0E+30
  REAL (wp), INTENT(in) :: xvalue
  REAL (wp), INTENT(inout) :: xres
  real (wp) :: y

  x = xvalue
  xres = 0.0E+00
!
!  Check the input.
!
  if ( xlgst <= x ) then
    return
  end if
  if ( x <= 0.0E+00 ) then
    return
  end if

!
!  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
!
  if ( x < 1.5E+00 ) then

    if ( x < 0.5E+00 ) then
      xres = - log ( x )
      y = x + 1.0E+00
!
!  Test whether X < machine epsilon.
!
      if ( y == 1.0E+00 ) then
        return
      end if

    else

      xres = 0.0E+00
      y = x
      x = ( x - 0.5E+00 ) - 0.5E+00

    end if

    xres = xres + x * (((( &
        r1(5)   * y &
      + r1(4) ) * y &
      + r1(3) ) * y &
      + r1(2) ) * y &
      + r1(1) ) / (((( &
                  y &
      + r1(9) ) * y &
      + r1(8) ) * y &
      + r1(7) ) * y &
      + r1(6) )

    return

  end if
!
!  Calculation for 1.5 <= X < 4.0.
!
  if ( x < 4.0E+00 ) then

    y = ( x - 1.0E+00 ) - 1.0E+00

    xres = y * (((( &
        r2(5)   * x &
      + r2(4) ) * x &
      + r2(3) ) * x &
      + r2(2) ) * x &
      + r2(1) ) / (((( &
                  x &
      + r2(9) ) * x &
      + r2(8) ) * x &
      + r2(7) ) * x &
      + r2(6) )
!
!  Calculation for 4.0 <= X < 12.0.
!
  else if ( x < 12.0E+00 ) then

    xres = (((( &
        r3(5)   * x &
      + r3(4) ) * x &
      + r3(3) ) * x &
      + r3(2) ) * x &
      + r3(1) ) / (((( &
                  x &
      + r3(9) ) * x &
      + r3(8) ) * x &
      + r3(7) ) * x &
      + r3(6) )
!
!  Calculation for 12.0 <= X.
!
  else

    y = log ( x )
    xres = x * ( y - 1.0E+00 ) - 0.5E+00 * y + alr2pi

    if ( x <= xlge ) then

      x1 = 1.0E+00 / x
      x2 = x1 * x1

      xres = xres + x1 * ( ( &
             r4(3)   * &
        x2 + r4(2) ) * &
        x2 + r4(1) ) / ( ( &
        x2 + r4(5) ) * &
        x2 + r4(4) )

    end if

   end if

   END SUBROUTINE alngam


   SUBROUTINE gamain( x, p, xres )
!*****************************************************************************80
!
!! GAMAIN computes the incomplete gamma ratio.
!
!  Discussion:
!
!    A series expansion is used if P > X or X <= 1.  Otherwise, a
!    continued fraction approximation is used.
!
!  Modified:
!
!    17 January 2008
!
!  Author:
!
!    G Bhattacharjee
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    G Bhattacharjee,
!    Algorithm AS 32:
!    The Incomplete Gamma Integral,
!    Applied Statistics,
!    Volume 19, Number 3, 1970, pages 285-287.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, P, the parameters of the incomplete 
!    gamma ratio.  0 <= X, and 0 < P.
!
!    0, no errors.
!    1, P <= 0.
!    2, X < 0.
!    3, underflow.
!    4, error return from the Log Gamma routine.
!
!    Output, real ( kind = 8 ) GAMAIN, the value of the incomplete
!    gamma ratio.
!

  real (wp), intent(out) :: xres
  real (wp) a
  real (wp), parameter :: acu = 1.0E-08
  real (wp) an
  real (wp) arg
  real (wp) b
  real (wp) dif
  real (wp) factor
  real (wp) g
  real (wp) gin
  integer i
  real (wp), parameter :: oflo = 1.0E+37
  REAL (wp), INTENT(in) :: p
  real (wp) pn(6)
  real (wp) rn
  real (wp) term
  real (wp), parameter :: uflo = 1.0E-37
  REAL (wp), intent(in) :: x
!
!  Check the input.
!
  if ( p <= 0.0E+00 ) then
    xres = 0.0E+00
    return
  end if

  if ( x < 0.0E+00 ) then
    xres = 0.0E+00
    return
  end if

  if ( x == 0.0E+00 ) then
    xres = 0.0E+00
    return
  end if

  CALL alngam ( p, g )

  arg = p * log ( x ) - x - g

  if ( arg < log ( uflo ) ) then
    xres = 0.0E+00
    return
  end if

  factor = exp ( arg )
!
!  Calculation by series expansion.
!
  if ( x <= 1.0E+00 .or. x < p ) then

    gin = 1.0E+00
    term = 1.0E+00
    rn = p

    do

      rn = rn + 1.0E+00
      term = term * x / rn
      gin = gin + term

      if ( term <= acu ) then
        exit
      end if

    end do

    xres = gin * factor / p
    return

  end if
!
!  Calculation by continued fraction.
!
  a = 1.0E+00 - p
  b = a + x + 1.0E+00
  term = 0.0E+00

  pn(1) = 1.0E+00
  pn(2) = x
  pn(3) = x + 1.0E+00
  pn(4) = x * b

  gin = pn(3) / pn(4)

  do

    a = a + 1.0E+00
    b = b + 2.0E+00
    term = term + 1.0E+00
    an = a * term
    do i = 1, 2
      pn(i+4) = b * pn(i+2) - an * pn(i)
    end do

    if ( pn(6) /= 0.0E+00 ) then

      rn = pn(5) / pn(6)
      dif = abs ( gin - rn )
!
!  Absolute error tolerance satisfied?
!
      if ( dif <= acu ) then
!
!  Relative error tolerance satisfied?
!
        if ( dif <= acu * rn ) then
          xres = 1.0E+00 - factor * gin
          exit
        end if

      end if

      gin = rn

    end if

    do i = 1, 4
      pn(i) = pn(i+2)
    end do
    if ( oflo <= abs ( pn(5) ) ) then

      do i = 1, 4
        pn(i) = pn(i) / oflo
      end do

    end if

  end do

  END SUBROUTINE gamain


   !!======================================================================
END MODULE p4zpoc
