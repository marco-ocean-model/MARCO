










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









MODULE p4zligand
   !!======================================================================
   !!                         ***  MODULE p4zligand  ***
   !! TOP :    Compute remineralization/dissolution of organic ligands
   !!=========================================================================
   !! History :   3.6  !  2016-03  (O. Aumont, A. Tagliabue) Quota model and reorganization
   !!             3.*  !  2025-02  (S. Maishal, R. Person)  and optimization
   !!----------------------------------------------------------------------
   !!   p4z_ligand     :  Compute remineralization/dissolution of organic ligands
   !!   p4z_ligand_init:  Initialisation of parameters for remineralisation
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE prtctl          ! print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_ligand         ! called in p4zbio.F90
   PUBLIC   p4z_ligand_init    ! called in trcsms_pisces.F90

   REAL(wp), PUBLIC ::  rlgw     !: lifetime (years) of weak ligands
   REAL(wp), PUBLIC ::  rlgs     !: lifetime (years) of strong ligands
   REAL(wp), PUBLIC ::  rlig     !: Remin ligand production
   REAL(wp), PUBLIC ::  prlgw    !: Photochemical of weak ligand
   REAL(wp), PUBLIC ::  xklig    !: 1/2 saturation constant of photolysis

   REAL(wp) ::  xklig2    !: xklig * xklig
   LOGICAL  ::  l_dia_ligand

   !! * Substitutions











































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zligand.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_ligand( kt, knt, Kbb, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_ligand  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic ligands
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step
      INTEGER, INTENT(in)  ::  Kbb, Krhs ! time level indices
      !
      INTEGER  :: mi, mj, jk
      REAL(wp) :: zlig, zlig2, zlig3
      REAL(wp) :: zlgwp, zlgwpr, zlgwr, zlablgw 
      REAL(wp) :: zfecoag, zaggliga, zligco
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      CHARACTER (len=25) ::   charout
      !!---------------------------------------------------------------------
      !
      IF( kt == nittrc000 )  &
          &  l_dia_ligand = iom_use( "LIGREM" ) .OR. iom_use( "LIGPR" ) &
          &            .OR. iom_use( "LPRODR" ) .OR. iom_use( "LGWCOLL" )

      IF( .false. )   CALL timing_start('p4z_ligand')
      !
      ! production from remineralisation of organic matter
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zlig  = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jplgw)
         zlig2 = zlig  * zlig
         zlig3 = zlig2 * zlig
         zfecoag = xfecolagg(mi,mj,jk) * 1.0E-9
         zligco  = MAX(0., zlig -  zfecoag)
         !
         ! production from remineralisation of organic matter
         zlgwp = orem(mi,mj,jk) * rlig
         ! Decay of weak ligand
         ! This is based on the idea that as LGW is lower
         ! there is a larger fraction of refractory OM
         !
         zlgwr = ( 1.0 / rlgs * zligco + 1.0 / rlgw * zfecoag ) / ( 0.5*EPSILON(1.e0) + zlig )
         zlgwr = zlgwr * tgfunc(mi,mj,jk) * ( xstep / year2day ) * blim(mi,mj,jk) * zlig 
         !
         ! photochem loss of weak ligand
         zlgwpr = prlgw * xstep * etot(mi,mj,jk) * zlig3 * (1. - fr_i(mi,mj) ) / ( zlig2 + xklig2 )
         !
         ! Coagulation of ligands due to various processes (Brownian, shear, diff. sedimentation
         ! xcoagfe is computed in p4zfechem
         ! 50% of the ligands are supposed to be in the colloidal size fraction as for FeL
         zaggliga = xcoagfe(mi,mj,jk) * xstep * 0.5 * zligco
         !
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jplgw) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jplgw) &
                 &                 + zlgwp  - zlgwr - zlgwpr - zaggliga
         !
      END DO   ;   END DO   ;   END DO
      
      IF( l_dia_ligand .AND. ( .false. .AND. knt == nrdttrc ) ) THEN
         ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp     
         !
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zw3d(mi,mj,N+1-jk) =  orem(mi,mj,jk) * rlig * 1e9 * 1.e+3 * rfact2r * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO   
         CALL iom_put( "LPRODR", zw3d )
         !
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zlig = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jplgw)
            zw3d(mi,mj,N+1-jk)  = ( 1.0 / rlgs * MAX(0., zlig - xfecolagg(mi,mj,jk) * 1.0E-9 )    &
              &                  + 1.0 / rlgw * xfecolagg(mi,mj,jk) * 1.0E-9 ) / ( 0.5*EPSILON(1.e0) + zlig ) &
              &                 * tgfunc(mi,mj,jk) * ( xstep / year2day ) * blim(mi,mj,jk) * zlig 
         END DO   ;   END DO   ;   END DO   
         CALL iom_put( "LIGREM", zw3d )
         !
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zlig = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jplgw)
            zw3d(mi,mj,N+1-jk)  = prlgw * xstep * etot(mi,mj,jk) * zlig*zlig*zlig &
              &            * fr_i(mi,mj) / ( zlig*zlig + xklig2 )
         END DO   ;   END DO   ;   END DO   
         CALL iom_put( "LIGPR", zw3d )
         !
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zlig = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jplgw)
            zw3d(mi,mj,N+1-jk)  =  xcoagfe(mi,mj,jk) * xstep &
              &        * 0.5 * MAX(0., zlig - xfecolagg(mi,mj,jk) * 1.0E-9 )
         END DO   ;   END DO   ;   END DO   
         CALL iom_put( "LGWCOLL", zw3d )
         !
         DEALLOCATE( zw3d ) 
      ENDIF

      !
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('ligand1')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
 !        CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( .false. )   CALL timing_stop('p4z_ligand')
      !
   END SUBROUTINE p4z_ligand


   SUBROUTINE p4z_ligand_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_ligand_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampislig namelist and check the parameters
      !!
      !! ** input   :   Namelist nampislig
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer 
      !
      NAMELIST/nampislig/ rlgw, prlgw, rlgs, rlig, xklig
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*)
         WRITE(stdout,*) 'p4z_ligand_init : remineralization/scavenging of organic ligands'
         WRITE(stdout,*) '~~~~~~~~~~~~~~~'
      ENDIF
      REWIND(numnatp_ref);READ(numnatp_ref,nampislig,IOSTAT=ios);CALL ctl_nam(ios,"nampislig (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,nampislig,IOSTAT=ios);CALL ctl_nam(ios,"nampislig (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE ( numonp, nampislig )
      !
      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*) '   Namelist : nampislig'
         WRITE(stdout,*) '      Lifetime (years) of weak ligands             rlgw  =', rlgw
         WRITE(stdout,*) '      Remin ligand production per unit C           rlig  =', rlig
         WRITE(stdout,*) '      Photolysis of weak ligand                    prlgw =', prlgw
         WRITE(stdout,*) '      Lifetime (years) of strong ligands           rlgs  =', rlgs
         WRITE(stdout,*) '      1/2 saturation for photolysis                xklig =', xklig
      ENDIF
      !
      xklig2 = xklig * xklig
      !
   END SUBROUTINE p4z_ligand_init


   !!======================================================================
END MODULE p4zligand
