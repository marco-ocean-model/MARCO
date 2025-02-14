










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









MODULE p4zsink
   !!======================================================================
   !!                         ***  MODULE p4zsink  ***
   !! TOP :    vertical flux of particulate matter due to 
   !!        gravitational sinking
   !!        This module is the same for both  and -QUOTA
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!             3.5  !  2012-07  (O. Aumont) Introduce potential time-splitting
   !!             4.0  !  2019     (O. Aumont) an external subroutine is called
   !!                                          to compute the impact of sinking
   !!----------------------------------------------------------------------
   !!   p4z_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   p4z_sink_init  :  Unitialisation of sinking speed parameters
   !!   p4z_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !   Source Minus Sink variables
   USE trcsink         !  General routine to compute sedimentation
   USE prtctl          !  print control for debugging
   USE iom             !  I/O manager
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sink         ! called in p4zbio.F90
   PUBLIC   p4z_sink_init    ! called in trcini_pisces.F90
   PUBLIC   p4z_sink_alloc   ! called in trcini_pisces.F90

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  sinkpocb  !: POC sinking fluxes at bottom 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  sinkcalb  !: CaCO3 sinking fluxes at bottom
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  sinksilb  !: BSi sinking fluxes at bottom
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  sinkponb  !: POC sinking fluxes at bottom 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  sinkpopb  !: POC sinking fluxes at bottom 
   INTEGER,  PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: nlev100

   REAL(wp) :: xfact
   LOGICAL  :: l_dia_sink, l_diag

   !! * Substitutions






























   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zsink.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_sink ( kt, knt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking. 
      !!
      !! ** Method  : - An external advection subroutine is called to compute
      !!                the impact of sinking on the particles. The tracers
      !!                concentrations are updated in this subroutine which
      !!                is mandatory to deal with negative concentrations
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      INTEGER, INTENT(in) :: Kbb, Kmm, Krhs  ! time level indices
      INTEGER  ::   mi, mj, jk, ikb, ik1
      CHARACTER (len=25) :: charout
      REAL(wp) :: zmax, zfact
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N+1) :: zsinking, zsinking2
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zw2d
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p4z_sink')

      IF( kt == nittrc000 )  THEN
         l_dia_sink = iom_use( "EPC100" ) .OR. iom_use( "EPFE100" ) .OR. iom_use( "EPCAL100" ) .OR. iom_use( "EPSI100" )  &
           &    .OR.  iom_use( "EXPC"   ) .OR. iom_use( "EXPFE"   ) .OR. iom_use( "EXPCAL"   ) .OR. iom_use( "EXPSI"   )  &
           &    .OR.  iom_use( "tcexp"  ) 
         !
         xfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
      ENDIF
      l_diag = l_dia_sink .AND. knt == nrdttrc

      !
      nlev100(:,:) = N        !  last level where depth less than 100 m
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         DO jk = N, 1, -1
            IF( ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) > 100. )  nlev100(mi,mj) = jk - 1
         END DO
      END DO   ;   END DO

      ! Initialize the local arrays
      zsinking(:,:,:) = 0._wp ; zsinking2(:,:,:) = 0._wp

      ! Sinking speeds of big detritus is increased with depth as shown
      ! by data and from the coagulation theory. This is controled by
      ! wsbio2max and wsbio2scale. If wsbio2max is set to wsbio2, then
      ! sinking speed is constant with depth.
      ! CaCO3 and bSi are supposed to sink at the big particles speed 
      ! due to their high density
      ! ---------------------------------------------------------------
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zmax  = MAX( heup_01(mi,mj), hbl(mi,mj) )
         zfact = MAX( 0., gdepw(mi,mj,jk+1,Kmm) - zmax ) / wsbio2scale
         wsbio4(mi,mj,jk) = wsbio2 + MAX(0., ( wsbio2max - wsbio2 )) * zfact
      END DO   ;   END DO   ;   END DO
      wsbio4(:,:,N) = wsbio2

      ! Sinking speed of the small particles is always constant
      ! except in  reduced
      IF( ln_p2z ) THEN
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zmax  = MAX( heup_01(mi,mj), hbl(mi,mj) )
            zfact = MAX( 0., gdepw(mi,mj,jk+1,Kmm) - zmax ) &
              &      * tgfunc(mi,mj,jk) * ( 0.03 / wsbio2 - 0.015/ wsbio)
            zfact = MIN(3.0, EXP(-zfact) )
            wsbio3(mi,mj,jk) = ( wsbio + wsbio2 * ( sizen(mi,mj,1) - 1.0 ) * 0.05 * zfact )    &
               &               / ( 1.0 + ( sizen(mi,mj,1) - 1.0 ) * 0.05 * zfact )
         END DO   ;   END DO   ;   END DO
         wsbio3(:,:,N) = wsbio
      ELSE
         wsbio3(:,:,:) = wsbio
      ENDIF

      ! Compute the sedimentation term using trc_sink for all the sinking particles
      ! ---------------------------------------------------------------------------
      CALL trc_sink( kt, Kbb, Kmm, wsbio3, zsinking , jppoc, rfact2 )
      IF( .NOT. ln_p2z ) THEN
         CALL trc_sink( kt, Kbb, Kmm, wsbio4, zsinking2, jpgoc, rfact2 )
      ENDIF
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ikb = N+1
         sinkpocb(mi,mj) = zsinking(mi,mj,ikb) + zsinking2(mi,mj,ikb)
      END DO   ;   END DO
      !
      IF( l_diag ) THEN
         ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp     
         ALLOCATE( zw2d(-2:Lm+3+padd_X,-2:Mm+3+padd_E)     )  ;  zw2d(:,:) = 0._wp
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zw3d(mi,mj,N+1-jk) = ( zsinking(mi,mj,jk) + zsinking2(mi,mj,jk) ) &
                   &        * xfact * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            ik1 = nlev100(mi,mj)
            zw2d(mi,mj) = ( zsinking(mi,mj,ik1) + zsinking2(mi,mj,ik1) ) &
                   &        * xfact * tmask(mi,mj,ik1)
         END DO   ;   END DO
         CALL iom_put( "EPC100",  zw2d )  ! Export of carbon at 100m
         CALL iom_put( "EXPC"  ,  zw3d )             ! Export of carbon in the water column
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zw2d(mi,mj) = zw2d(mi,mj) * e1e2t(mi,mj)
         END DO   ;   END DO 
         t_oce_co2_exp  = glob_sum( 'p4zsink',  zw2d )
         CALL iom_put( "tcexp", t_oce_co2_exp )      ! Total cabon exort
      ENDIF
      IF( .NOT. ln_p2z ) THEN
         ! Compute the sedimentation term using trc_sink for all the sinking particles
         ! ---------------------------------------------------------------------------
         CALL trc_sink( kt, Kbb, Kmm, wsbio4, zsinking, jpcal, rfact2 )
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            ikb = N+1
            sinkcalb(mi,mj) = zsinking(mi,mj,ikb)
         END DO   ;   END DO
         IF( l_diag ) THEN
            DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zw3d(mi,mj,N+1-jk) =  zsinking(mi,mj,jk) * xfact * tmask(mi,mj,jk)
            END DO   ;   END DO   ;   END DO
            DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               ik1 = nlev100(mi,mj)
               zw2d(mi,mj) = zsinking(mi,mj,ik1) * xfact * tmask(mi,mj,ik1)
            END DO   ;   END DO
            CALL iom_put( "EPCAL100",  zw2d )  ! Export of calcite at 100m
            CALL iom_put( "EXPCAL"  ,  zw3d )  ! Export of calcite in the water column
         ENDIF
         !
         CALL trc_sink( kt, Kbb, Kmm, wsbio4, zsinking, jpgsi, rfact2 )
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            ikb = N+1
            sinksilb(mi,mj) = zsinking(mi,mj,ikb)
         END DO   ;   END DO
         IF( l_diag ) THEN
            DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zw3d(mi,mj,N+1-jk) =  zsinking(mi,mj,jk) * xfact * tmask(mi,mj,jk)
            END DO   ;   END DO   ;   END DO
            DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               ik1 = nlev100(mi,mj)
               zw2d(mi,mj) = zsinking(mi,mj,ik1) * xfact * tmask(mi,mj,ik1)
            END DO   ;   END DO
            CALL iom_put( "EPSI100",  zw2d )  ! Export of Silicate at 100m
            CALL iom_put( "EXPSI"  ,  zw3d )             ! Export of Silicate in the water column
         ENDIF
         !
         CALL trc_sink( kt, Kbb, Kmm, wsbio3, zsinking, jpsfe, rfact2 )
         CALL trc_sink( kt, Kbb, Kmm, wsbio4, zsinking2, jpbfe, rfact2 )
         IF( l_diag ) THEN
            DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
              zw3d(mi,mj,N+1-jk) = ( zsinking(mi,mj,jk) + zsinking2(mi,mj,jk) ) &
               &        * xfact * tmask(mi,mj,jk)
            END DO   ;   END DO   ;   END DO
            DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               ik1 = nlev100(mi,mj)
               zw2d(mi,mj) = ( zsinking(mi,mj,ik1) + zsinking2(mi,mj,ik1) ) &
                      &        * xfact * tmask(mi,mj,ik1)
            END DO   ;   END DO
            CALL iom_put( "EPFE100",  zw2d )  ! Export of iron at 100m
            CALL iom_put( "EXPFE"  ,  zw3d )             ! Export of iron in the water column
         ENDIF
      ENDIF
      !
      ! -QUOTA part
      IF( ln_p5z ) THEN
         !
         CALL trc_sink( kt, Kbb, Kmm, wsbio3, zsinking, jppon, rfact2 )
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            ikb = N+1
            sinkponb(mi,mj) = zsinking(mi,mj,ikb)
         END DO   ;   END DO
         !
         CALL trc_sink( kt, Kbb, Kmm, wsbio4, zsinking, jpgon, rfact2 )
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            ikb = N+1
            sinkponb(mi,mj) = sinkponb(mi,mj) + zsinking(mi,mj,ikb)
         END DO   ;   END DO
         !
         CALL trc_sink( kt, Kbb, Kmm, wsbio3, zsinking, jppop, rfact2 )
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            ikb = N+1
            sinkpopb(mi,mj) = zsinking(mi,mj,ikb)
         END DO   ;   END DO
         !
         CALL trc_sink( kt, Kbb, Kmm, wsbio4, zsinking, jpgop, rfact2 )
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            ikb = N+1
            sinkpopb(mi,mj) = sinkpopb(mi,mj) + zsinking(mi,mj,ikb)
         END DO   ;   END DO
      ENDIF
      !
      IF( l_diag )  DEALLOCATE( zw3d, zw2d )
      !
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
      !   CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Kbb,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( .false. )   CALL timing_stop('p4z_sink')
      !
   END SUBROUTINE p4z_sink


   SUBROUTINE p4z_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sink_init  ***
      !!
      !! ** Purpose :   Initialization of sinking parameters
      !!
      !! ** Method  :   
      !!
      !! ** input   :   
      !!----------------------------------------------------------------------
      !
      wsbio3(:,:,:) = wsbio
      wsbio4(:,:,:) = wsbio2
      sinkpocb(:,:) = 0.
      sinkcalb(:,:) = 0.
      IF( .NOT. ln_p2z )  sinksilb(:,:) = 0.
      IF( ln_p5z ) THEN
         sinkponb(:,:) = 0.
         sinkpopb(:,:) = 0.
      ENDIF

   END SUBROUTINE p4z_sink_init

   INTEGER FUNCTION p4z_sink_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(4)
      !!----------------------------------------------------------------------
      !
      ierr(:) = 0
      !
      ALLOCATE( sinkpocb(Istrp:Iendp,Jstrp:Jendp), sinkcalb(Istrp:Iendp,Jstrp:Jendp), STAT=ierr(1) )
      ALLOCATE( nlev100(Istrp:Iendp,Jstrp:Jendp), STAT=ierr(2) )
      !
      IF( .NOT. ln_p2z ) THEN
         ALLOCATE( sinksilb(Istrp:Iendp,Jstrp:Jendp), STAT=ierr(3) )
         !
         IF( ln_p5z ) &
             &   ALLOCATE( sinkponb(Istrp:Iendp,Jstrp:Jendp), sinkpopb(Istrp:Iendp,Jstrp:Jendp), STAT=ierr(4) )
      ENDIF
      !
      p4z_sink_alloc = MAXVAL( ierr )
      IF( p4z_sink_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p4z_sink_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p4z_sink_alloc

   
   !!======================================================================
END MODULE p4zsink
