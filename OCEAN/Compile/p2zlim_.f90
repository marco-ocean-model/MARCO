










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









MODULE p2zlim
   !!======================================================================
   !!                         ***  MODULE p2zlim  ***
   !! TOP :   Computes the nutrient limitation terms of phytoplankton
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-04  (O. Aumont, C. Ethe) Limitation for iron modelled in quota 
   !!             3.*  !  2025-01  (S. Maishal, R. Person) Change to High Performance
   !!----------------------------------------------------------------------
   !!   p2z_lim        :   Compute the nutrients limitation terms 
   !!   p2z_lim_init   :   Read the namelist 
   !!----------------------------------------------------------------------

   USE oce_trc         ! Shared ocean-passive tracers variables
   USE trc             ! Tracers defined
   USE sms_pisces      !  variables
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC p2z_lim           ! called in p4zbio.F90 
   PUBLIC p2z_lim_init      ! called in trcsms_pisces.F90 
   PUBLIC p2z_lim_alloc     ! called in trcini_pisces.F90


   !! * Shared module variables
   REAL(wp), PUBLIC ::  concnno3    !:  NO3, PO4 half saturation   
   REAL(wp), PUBLIC ::  concbno3    !:  NO3 half saturation  for bacteria 
   REAL(wp), PUBLIC ::  concnfer    !:  Iron half saturation for nanophyto 
   REAL(wp), PUBLIC ::  xsizephy    !:  Minimum size criteria for nanophyto
   REAL(wp), PUBLIC ::  xsizern     !:  Size ratio for nanophytoplankton
   REAL(wp), PUBLIC ::  xkdoc       !:  2nd half-sat. of DOC remineralization  
   REAL(wp), PUBLIC ::  concbfe     !:  Fe half saturation for bacteria 
   REAL(wp), PUBLIC ::  caco3r      !:  mean rainratio 

   !!* Phytoplankton limitation terms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanono3   !: Nanophyto limitation by NO3
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimphy    !: Nutrient limitation term of nanophytoplankton
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimbac    !: Bacterial limitation term
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimbacl   !: Bacterial limitation term
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimnfe    !: Nanophyto limitation by Iron

   LOGICAL  :: l_dia

   !! * Substitutions












































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p2zlim.F90 10069 2018-08-28 14:12:24Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p2z_lim( kt, knt, Kbb, Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_lim  ***
      !!
      !! ** Purpose :   Compute the co-limitations by the various nutrients
      !!                for the unique phytoplankton species 
      !!
      !! ** Method  : - Limitation is computed according to Monod formalism
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)  :: kt, knt
      INTEGER, INTENT(in)  :: Kbb, Kmm      ! time level indices
      !
      INTEGER  ::   mi, mj, jk
      REAL(wp) ::   zcoef, zconc0n, zconcnf, zlim1, zlim2, zlim3
      REAL(wp) ::   zbiron, ztem1, ztem2, zetot1, zetot2, zsize
      REAL(wp) ::   zferlim, zno3
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p2z_lim')
      !
      IF( kt == nittrc000 )  &
         &     l_dia = iom_use( "LNnut" ) .OR. iom_use( "SIZEN" ) .OR. iom_use( "xfracal" )

      !
      sizena(:,:,:) = 1.0
      !
      !! Parallelize the 3D loop for nutrient limitations using OpenMP
      !$OMP PARALLEL DO & 
   PRIVATE(mi, mj, jk, zcoef, zconc0n, zconcnf, zlim1, zlim2, zlim3, &
           zbiron, ztem1, ztem2, zetot1, zetot2, zsize) &
   SHARED(tr, sizen, plig, biron, xlimbacl, xlimbac, xnanono3, &
          xlimphy, xlimnfe, concnno3, concnfer, concbno3, concbfe, xkdoc)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp

         ! Tuning of the iron concentration to a minimum level that is set to the detection limit
         !-------------------------------------
         zno3    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) / 40.e-6
         zferlim = MAX( 5e-11 * zno3 * zno3, 2e-11 )
         zferlim = MIN( zferlim, 5e-11 )
         t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer) = MAX( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer), zferlim )
         
         ! Computation of a variable Ks for NO3 on phyto taking into account
         ! that increasing biomass is made of generally bigger cells
         ! The allometric relationship is classical.
         !------------------------------------------------
         zsize    = sizen(mi,mj,jk)**0.81
         zconc0n  = concnno3 * zsize
         zconcnf  = concnfer * zsize

         ! Nanophytoplankton
         zbiron = ( 75.0 * ( 1.0 - plig(mi,mj,jk) ) + plig(mi,mj,jk) ) * biron(mi,mj,jk)

         ! Michaelis-Menten Limitation term by nutrients of
         ! heterotrophic bacteria
         ! -------------------------------------------------
         zlim1  = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) / ( concbno3 + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) )
         zlim2  = zbiron / ( concbfe + zbiron )
         zlim3  = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) / ( xkdoc   + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) )

         ! Xlimbac is used for DOC solubilization whereas xlimbacl
         ! is used for all the other bacterial-dependent terms
         ! -------------------------------------------------------
         xlimbacl(mi,mj,jk) = MIN( zlim1, zlim2)
         xlimbac (mi,mj,jk) = xlimbacl(mi,mj,jk) * zlim3

         ! Michaelis-Menten Limitation term by nutrients: Nanophyto
         ! Optimal parameterization by Smith and Pahlow series of 
         ! papers is used. Optimal allocation is supposed independant
         ! for all nutrients. 
         ! --------------------------------------------------------

         ! Limitation of nanophytoplankton growth
         xnanono3(mi,mj,jk) = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) / ( zconc0n + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) )
         xlimnfe (mi,mj,jk) = zbiron / ( zbiron + zconcnf )
         xlimphy (mi,mj,jk) = MIN( xlimnfe(mi,mj,jk), xnanono3(mi,mj,jk) )
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO

      ! Size estimation of phytoplankton based on total biomass
      ! Assumes that larger biomass implies addition of larger cells
      ! ------------------------------------------------------------
      !$OMP PARALLEL DO & 
   PRIVATE(mi, mj, jk, zcoef) &
   SHARED(tr, sizen, sizena, xsizephy, xsizern)

      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zcoef = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) - MIN(xsizephy, t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) )
         sizena(mi,mj,jk) = 1. + ( xsizern -1.0 ) * zcoef / ( xsizephy + zcoef )
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO

      ! Compute the fraction of nanophytoplankton that is made of calcifiers
      ! This is a purely adhoc formulation described in Aumont et al. (2015)
      ! This fraction depends on nutrient limitation, light, temperature
      ! --------------------------------------------------------------------
      !$OMP PARALLEL DO & 
   PRIVATE(mi, mj, jk, ztem1, ztem2, zetot1, zetot2) &
   SHARED(xlimphy, caco3r, ts, etot_ndcy, xfracal, tmask)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ztem1  = MAX( 0., t(mi,mj,N+1-jk,nnew,itemp) + 1.8)
         ztem2  = t(mi,mj,N+1-jk,nnew,itemp) - 10.
         zetot1 = MAX( 0., etot_ndcy(mi,mj,jk) - 1.) / ( 4. + etot_ndcy(mi,mj,jk) ) 
         zetot2 = 30. / ( 30.0 + etot_ndcy(mi,mj,jk) )

         xfracal(mi,mj,jk) = caco3r * xlimphy(mi,mj,jk)                              &
            &                       * ztem1 / ( 0.1 + ztem1 )                        &
            &                       * MAX( 1., t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) / xsizephy )  &
            &                       * zetot1 * zetot2                                &
            &                       * ( 1. + EXP(-ztem2 * ztem2 / 25. ) )            &
            &                       * MIN( 1., 50. / ( hbl(mi,mj) + 0.5*EPSILON(1.e0) ) )
         xfracal(mi,mj,jk) = MIN( 0.8 , xfracal(mi,mj,jk) )
         xfracal(mi,mj,jk) = MAX( 0.02, xfracal(mi,mj,jk) )
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      IF( l_dia .AND. knt == nrdttrc ) THEN        ! save output diagnostics
        !
        ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
        ! fraction of calcifiers
        DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           zw3d(mi,mj,N+1-jk) = xfracal(mi,mj,jk) * tmask(mi,mj,jk)
        END DO   ;   END DO   ;   END DO
        CALL iom_put( "xfracal",  zw3d)
        ! Nutrient limitation term
        DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           zw3d(mi,mj,N+1-jk) = xlimphy(mi,mj,jk) * tmask(mi,mj,jk)
        END DO   ;   END DO   ;   END DO
        CALL iom_put( "LNnut",  zw3d)
        ! Size limitation term
        DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           zw3d(mi,mj,N+1-jk) = sizen(mi,mj,jk) * tmask(mi,mj,jk)
        END DO   ;   END DO   ;   END DO
        CALL iom_put( "SIZEN",  zw3d)
        !
        DEALLOCATE( zw3d )
        !
      ENDIF
      !
      IF( .false. )   CALL timing_stop('p2z_lim')
      !
   END SUBROUTINE p2z_lim


   SUBROUTINE p2z_lim_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p2z_lim_init  ***
      !!
      !! ** Purpose :   Initialization of the nutrient limitation parameters
      !!
      !! ** Method  :   Read the namp2zlim namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namp2zlim
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer

      ! Namelist block
      NAMELIST/namp2zlim/ concnno3, concbno3, concnfer, xsizephy, xsizern,  &
         &                concbfe, xkdoc, caco3r, oxymin
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*)
         WRITE(stdout,*) 'p2z_lim_init : initialization of nutrient limitations'
         WRITE(stdout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND(numnatp_ref);READ(numnatp_ref,namp2zlim,IOSTAT=ios);CALL ctl_nam(ios,"namp2zlim (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,namp2zlim,IOSTAT=ios);CALL ctl_nam(ios,"namp2zlim (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE( numonp, namp2zlim )

      !
      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*) '   Namelist : namp2zlim'
         WRITE(stdout,*) '      mean rainratio                           caco3r    = ', caco3r
         WRITE(stdout,*) '      NO3 half saturation of phyto             concnno3  = ', concnno3
         WRITE(stdout,*) '      Iron half saturation for nanophyto       concnfer  = ', concnfer
         WRITE(stdout,*) '      Fe half saturation for bacteria          concbfe   = ', concbfe
         WRITE(stdout,*) '      half-sat. of DOC remineralization        xkdoc     = ', xkdoc
         WRITE(stdout,*) '      size ratio for phytoplankton             xsizern   = ', xsizern
         WRITE(stdout,*) '      NO3 half saturation of bacteria          concbno3  = ', concbno3
         WRITE(stdout,*) '      Minimum size criteria for phyto          xsizephy  = ', xsizephy
         WRITE(stdout,*) '      halk saturation constant for anoxia      oxymin    =' , oxymin
      ENDIF
      !
      xnanono3(:,:,N) = 0._wp
      xlimphy (:,:,N) = 0._wp
      xlimnfe (:,:,N) = 0._wp
      xlimbac (:,:,N) = 0._wp
      xlimbacl(:,:,N) = 0._wp
      !
   END SUBROUTINE p2z_lim_init


   INTEGER FUNCTION p2z_lim_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_lim_alloc  ***
      !! 
      !            Allocation of the arrays used in this module
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_stop
      !!----------------------------------------------------------------------

      !*  Biological arrays for phytoplankton growth
      ALLOCATE( xnanono3(Istrp:Iendp,Jstrp:Jendp,N), xlimphy (Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xlimnfe (Istrp:Iendp,Jstrp:Jendp,N), xlimbac (Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xlimbacl(Istrp:Iendp,Jstrp:Jendp,N),                       STAT=p2z_lim_alloc )
 
      !
      IF( p2z_lim_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p2z_lim_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p2z_lim_alloc


   !!======================================================================
END MODULE p2zlim
