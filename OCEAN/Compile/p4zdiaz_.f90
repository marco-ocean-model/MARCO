










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









MODULE p4zdiaz
   !!======================================================================
   !!                         ***  MODULE p4zdiaz  ***
   !! TOP :    Compute Nitrogen Fixation ( Diazotrophy )
   !!         This module is common to both  and -QUOTA
   !!=========================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             5.0  !  2023-12  (O. Aumont, C. Ethe)
   !!             3.*  !  2025-02  (S. Maishal, R. Person)  and optimization
   !!----------------------------------------------------------------------
   !!   p4z_rem       :  Compute remineralization/dissolution of organic compounds
   !!   p4z_rem_init  :  Initialisation of parameters for remineralisation
   !!   p4z_rem_alloc :  Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !   Source Minus Sink variables
   USE p4zche          !  chemical model
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p2zlim          !  Co-limitations of differents nutrients
   USE p4zlim          !  Co-limitations of differents nutrients
   USE prtctl          !  print control for debugging
   USE iom             !  I/O manager


   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_diaz         ! called in p4zbio.F90
   PUBLIC   p4z_diaz_init    ! called in trcini_pisces.F90
   PUBLIC   p4z_diaz_alloc   ! called in trcini_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::   nitrfix      !: Nitrogen fixation rate
   REAL(wp), PUBLIC ::   diazolight   !: Nitrogen fixation sensitivty to light
   REAL(wp), PUBLIC ::   concfediaz   !: Fe half-saturation Cste for diazotrophs

   REAL(wp), SAVE :: r1_rday

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: nitrpot    !: Nitrogen fixation

   LOGICAL         :: l_dia_nfix

   !! * Substitutions












































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zrem.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS


   SUBROUTINE p4z_diaz( kt, knt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_diaz  ***
      !!
      !! ** Purpose : - Compute Nitrogen Fixation 
      !!                Small source iron from particulate inorganic iron
      !!
      !! ** Method  : - Potential nitrogen fixation dependant on temperature and iron
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt         ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      !
      INTEGER  ::   mi, mj, jk
      !
      REAL(wp) ::  ztrfer, ztrpo4s, ztrdp, zwdust, zmudia, ztemp
      REAL(wp) ::  zsoufer, zlight, ztrpo4, ztrdop, zratpo4
      REAL(wp) ::  zfact, zlim, zdiano3, zdianh4
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zw3d
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p4z_diaz')
      !
      IF( kt == nittrc000 )  l_dia_nfix   = iom_use( "Nfix" ) .OR. iom_use( "Nfixo2" )

      ! Nitrogen fixation process
      !$OMP PARALLEL DO & 
      !$OMP PRIVATE(mi, mj, jk, zlight, ztemp, zmudia, zdiano3, zdianh4, &
                                    zlim, zfact, ztrfer, ztrpo4, ztrdop) &
      !$OMP SHARED(etot_ndcy, diazolight, fr_i, ts, rno3, tr, jpno3, jpnh4, &
                   jppo4, jpdop, biron, concnno3, concnnh4, concfediaz, &
                   nitrpot, rfact2, r1_rday, 0.5*EPSILON(1.e0), ln_p2z, ln_p5z)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         !                      ! Potential nitrogen fixation dependant on temperature and iron
         zlight  =  ( 1.- EXP( -etot_ndcy(mi,mj,jk) / diazolight ) ) &
                 &  * ( 1. - fr_i(mi,mj) )
         !
         ztemp = t(mi,mj,N+1-jk,nnew,itemp)
         zmudia = MAX( 0.,-0.001096*ztemp*ztemp + 0.057*ztemp -0.637 ) / rno3
         !       Potential nitrogen fixation dependant on temperature and iron
         IF( ln_p2z ) THEN
            zdiano3 = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) &
                &    / ( concnno3 + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) )
            zlim    = 1.- zdiano3
            zfact   = zlim * rfact2
            ztrfer  = biron(mi,mj,jk) / ( concfediaz + biron(mi,mj,jk) )
            nitrpot(mi,mj,jk) =  zmudia * r1_rday * zfact * ztrfer * zlight
         ELSE
            zdianh4 = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) &
                &   / ( concnnh4 + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) )
            zdiano3 = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) &
                &   / ( concnno3 + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) ) * (1. - zdianh4)
            zlim    = ( 1.- zdiano3 - zdianh4 )
            zfact   = zlim * rfact2
            ztrfer  = biron(mi,mj,jk) / ( concfediaz + biron(mi,mj,jk) )
            ztrpo4  = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) &
               &    / ( 1E-6 + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) )
            IF (ln_p5z) THEN
               ztrdop  = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdop) &
                  &   / ( 1E-6 + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdop) ) * (1. - ztrpo4)
               ztrpo4  =  ztrpo4 + ztrdop + 0.5*EPSILON(1.e0)
            ENDIF
            nitrpot(mi,mj,jk) =  zmudia * r1_rday * zfact * MIN( ztrfer, ztrpo4 ) * zlight
         ENDIF
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      ! Nitrogen change due to nitrogen fixation
      ! ----------------------------------------
      IF( ln_p2z ) THEN
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk, zfact, zlight, zsoufer) &
         !$OMP SHARED(nitrpot, nitrfix, etot_ndcy, diazolight, fr_i, biron, &
                      tr, jpno3, jptal, jpdic, jpdoc, jppoc, jpfer, jpoxy, &
                      Krhs, rno3, o2ut, o2nit, feratz, rfact2, day2sec)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zfact = nitrpot(mi,mj,jk) * nitrfix
            zlight  =  ( 1.- EXP( -etot_ndcy(mi,mj,jk) / diazolight ) ) &
                  &  * ( 1. - fr_i(mi,mj) )
            zsoufer = zlight * 2E-11 / ( 2E-11 + biron(mi,mj,jk) )
            !
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpno3) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpno3) + zfact / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) - rno3 * zfact / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) - zfact * 2.0 / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) + zfact * 1.0 / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) + zfact * 1.0 / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) - zfact * 1.0 / 3.0 * feratz
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) + ( o2ut + o2nit ) * zfact
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) &
                    &                + 0.005 * 4E-10 * zsoufer * rfact2 / day2sec
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
      ELSE
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk, zfact, zlight, zsoufer) &
         !$OMP SHARED(nitrpot, nitrfix, etot_ndcy, diazolight, &
                      fr_i, biron, tr, jpnh4, jptal, jpdic, jpdoc, &
                      jppoc, jpgoc, jpoxy, jpfer, jpsfe, jpbfe, Krhs, &
                      rno3, o2ut, o2nit, rfact2, day2sec)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zfact = nitrpot(mi,mj,jk) * nitrfix
            zlight  =  ( 1.- EXP( -etot_ndcy(mi,mj,jk) / diazolight ) ) &
                  &  * ( 1. - fr_i(mi,mj) )
            zsoufer = zlight * 2E-11 / ( 2E-11 + biron(mi,mj,jk) )
            !
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnh4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnh4) + zfact / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) + rno3 * zfact / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) - zfact * 2.0 / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) + zfact * 1.0 / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppoc) &
                    &                 + zfact * 1.0 / 3.0 * 2.0 / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgoc) &
                    &                 + zfact * 1.0 / 3.0 * 1.0 / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) &
                    &           + ( o2ut + o2nit ) * zfact * 2.0 / 3.0 + o2nit * zfact / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) &
                    &                - 30E-6 * zfact * 1.0 / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) &
                    &                 + 30E-6 * zfact * 1.0 / 3.0 * 2.0 / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) &
                    &                 + 30E-6 * zfact * 1.0 / 3.0 * 1.0 / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) &
                    &                + 0.003 * 4E-10 * zsoufer * rfact2 / day2sec
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
      ENDIF
      !
      IF( ln_p4z ) THEN
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk, zfact) &
         !$OMP SHARED(nitrpot, nitrfix, tr, jppo4, jpdoc, &
                      Kbb, Krhs, concdnh4, xstep)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zfact = nitrpot(mi,mj,jk) * nitrfix
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) - zfact * 2.0 / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) + concdnh4 &
                 &                  / ( concdnh4 + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) ) &
                 &                    * 0.001 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) * xstep
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
      ENDIF
      !
      IF( ln_p5z ) THEN
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk, ztrpo4, ztrdop, zratpo4, zfact) &
         !$OMP SHARED(nitrpot, nitrfix, tr, jppo4, jpdon, jpdop, &
                      jppon, jppop, jpgon, jpgop, Krhs, Kbb, 0.5*EPSILON(1.e0))
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            ztrpo4  = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) &
              &      / ( 1E-6 + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) )
            ztrdop  = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdop) &
               &     / ( 1E-6 + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdop) ) * (1. - ztrpo4)
            zratpo4 = ztrpo4 / (ztrpo4 + ztrdop + 0.5*EPSILON(1.e0))
            !
            zfact = nitrpot(mi,mj,jk) * nitrfix
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) &
            &                     - 16.0 / 46.0 * zfact * 2.0 / 3.0 * zratpo4
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdon) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdon) + zfact * 1.0 / 3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdop) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdop) &
            &                         + 16.0 / 46.0 * zfact / 3.0  &
            &                        - 16.0 / 46.0 * zfact * 2.0 / 3.0 * (1.0 - zratpo4)
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppon) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppon) &
            &                         + zfact * 1.0 / 3.0 * 2.0 /3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppop) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppop) &
            &                      + 16.0 / 46.0 * zfact * 1.0 / 3.0 * 2.0 /3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgon) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgon) &
            &                       + zfact * 1.0 / 3.0 * 1.0 /3.0
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgop) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgop) &
            &                       + 16.0 / 46.0 * zfact * 1.0 / 3.0 * 1.0 /3.0
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
      ENDIF
         
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('diaz')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
  !       CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF

      IF( l_dia_nfix .AND. .false. .AND. knt == nrdttrc ) THEN 
         ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
         zfact = rno3 * 1.e+3 * rfact2r !  conversion from molC/l/kt  to molN/m3/s
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, N+1-jk, zfact) &
         !$OMP SHARED(nitrpot, zw3d, tmask)
         DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zw3d(mi,mj,N+1-jk) =  nitrpot(mi,mj,jk) * zfact * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
         CALL iom_put( "Nfix", zw3d ) ! diazotrophy
         CALL iom_put( "Nfixo2", zw3d * o2nit) ! O2 production by diazotrophy
         DEALLOCATE( zw3d ) 
      ENDIF

      zfact = 1.e+3 * rfact2r
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
        bioVSink(mi,mj,Nnitrpot) = 0.
      END DO   ;   END DO
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zfact) &
      !$OMP SHARED(trc2d, nitrpot, nitrfix, rno3, e3t, tmask, Nnitrpot)
      DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         bioVSink(mi,mj,Nnitrpot) = bioVSink(mi,mj,Nnitrpot ) &
            &                 +  nitrpot(mi,mj,jk) * nitrfix * rno3    &
            &                 * zfact * Hz(mi,mj,N+1-jk) * tmask(mi,mj,jk) ! nitrogen fixation at surface
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zfact) &
      !$OMP SHARED(trc3d, nitrpot, nitrfix, rno3, o2nit, tmask, Nfixo2)
      DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         bioFlux(mi,mj,N+1-jk,Nfixo2 ) = nitrpot(mi,mj,jk) * nitrfix * rno3 &
           &                      * zfact * o2nit * tmask(mi,mj,jk)  ! O2 production by Nfix
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      IF( .false. )   CALL timing_stop('p4z_diaz')
      !
   END SUBROUTINE p4z_diaz


   SUBROUTINE p4z_diaz_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_diaz_init  ***
      !!
      !! ** Purpose :   Initialization of diazotrophy parameters
      !!
      !! ** Method  :   Read the nampisdiaz namelist and check the parameters
      !!
      !! ** input   :   Namelist nampisdiaz
      !!
      !!----------------------------------------------------------------------
      NAMELIST/nampisdiaz/nitrfix, diazolight, concfediaz
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*)
         WRITE(stdout,*) 'p4z_diaz_init : Initialization of diazotrophy parameters'
         WRITE(stdout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND(numnatp_ref);READ(numnatp_ref,nampisdiaz,IOSTAT=ios);CALL ctl_nam(ios,"nampisdiaz (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,nampisdiaz,IOSTAT=ios);CALL ctl_nam(ios,"nampisdiaz (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE( numonp, nampisdiaz )

      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*) '   Namelist parameters for diazotrophy, nampisdiaz'
         WRITE(stdout,*) '      nitrogen fixation rate                       nitrfix = ', nitrfix
         WRITE(stdout,*) '      nitrogen fixation sensitivty to light    diazolight  = ', diazolight
         IF( .NOT. ln_p2z ) &
           &  WRITE(stdout,*) '      Fe half-saturation cste for diazotrophs  concfediaz  = ', concfediaz
      ENDIF
      !
      nitrpot(:,:,:) = 0.
      r1_rday  = 1. / day2sec
      !
   END SUBROUTINE p4z_diaz_init


   INTEGER FUNCTION p4z_diaz_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_diaz_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( nitrpot(Istrp:Iendp,Jstrp:Jendp,N), STAT=p4z_diaz_alloc )
      !
      IF( p4z_diaz_alloc /= 0 )   CALL ctl_stop( 'STOP', 'p4z_diaz_alloc: failed to allocate arrays' )
      !
   END FUNCTION p4z_diaz_alloc


   !!======================================================================
END MODULE p4zdiaz
