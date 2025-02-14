










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









MODULE p4zrem
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP :    Compute remineralization/dissolution of organic compounds
   !!         except for POC which is treated in p4zpoc.F90
   !!         This module is common to both  and -QUOTA
   !!=========================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
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
   USE p2zlim          !  Nutrient limitation terms
   USE p4zlim          !  Nutrient limitation terms
   USE prtctl          !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_rem         ! called in p4zbio.F90
   PUBLIC   p2z_rem
   PUBLIC   p4z_rem_init    ! called in trcini_pisces.F90
   PUBLIC   p4z_rem_alloc   ! called in trcini_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::   xremikc    !: remineralisation rate of DOC (p5z) 
   REAL(wp), PUBLIC ::   xremikn    !: remineralisation rate of DON (p5z) 
   REAL(wp), PUBLIC ::   xremikp    !: remineralisation rate of DOP (p5z) 
   REAL(wp), PUBLIC ::   nitrif     !: NH4 nitrification rate 
   REAL(wp), PUBLIC ::   xsirem     !: remineralisation rate of biogenic silica
   REAL(wp), PUBLIC ::   xsiremlab  !: fast remineralisation rate of BSi
   REAL(wp), PUBLIC ::   xsilab     !: fraction of labile biogenic silica 
   REAL(wp), PUBLIC ::   feratb     !: Fe/C quota in bacteria
   REAL(wp), PUBLIC ::   xkferb     !: Half-saturation constant for bacterial Fe/C

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr   !: denitrification array

   LOGICAL         :: l_dia_remin, l_dia_bact, l_dia_denit

   !! * Substitutions












































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zrem.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p2z_rem( kt, knt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p2z_rem  ***
      !!
      !! ** Purpose :   Compute remineralization/dissolution of organic compounds
      !!                Computes also nitrification of ammonium 
      !!                The solubilization/remineralization of POC is treated 
      !!                in p4zpoc.F90. The dissolution of calcite is processed
      !!                in p4zlys.F90. 
      !!
      !! ** Method  : - Bacterial biomass is computed implicitely based on a 
      !!                parameterization developed from an explicit modeling
      !!                of  in an alternative version 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt         ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      !
      INTEGER  ::   mi, mj, jk
      REAL(wp) ::   zremik, zremikc, zfact
      REAL(wp) ::   zdep, zdepmin, zfactdep
      REAL(wp) ::   zammonic, zoxyremc, zolimic
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N) :: zdepbac
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp    ) :: ztempbac
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zw3d, zolimi
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p2z_rem')
      !
      IF( kt == nittrc000 )  THEN
         l_dia_remin  = iom_use( "REMIN" )
         l_dia_denit  = iom_use( "DENIT" )
         l_dia_bact   = iom_use( "BACT" )
      ENDIF
      IF( l_dia_remin ) THEN
         ALLOCATE( zolimi(Istrp:Iendp,Jstrp:Jendp,N) )    ;   zolimi(Istrp:Iendp,Jstrp:Jendp,N) = 0._wp
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zolimi(mi,mj,jk) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy)
         END DO   ;   END DO   ;   END DO
      ENDIF

      ! Computation of the mean bacterial concentration
      ! this parameterization has been deduced from a model version
      ! that was modeling explicitely bacteria. This is a very old parame
      ! that will be very soon updated based on results from a much more
      ! recent version of  with bacteria.
      ! ----------------------------------------------------------------
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zdep = MAX( hbl(mi,mj), heup_01(mi,mj), ((-1)*(z_r(mi,mj,N+1-1)-z_w(mi,mj,N))) )
         IF ( ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) <= zdep ) THEN
            zdepbac(mi,mj,jk) = 0.6 * ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo) * 1.0E6 )**0.6 * 1.E-6
            ztempbac(mi,mj)   = zdepbac(mi,mj,jk)
         ELSE
            zdepmin = zdep / ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N)))
            zdepbac(mi,mj,jk) = zdepmin**0.683 * ztempbac(mi,mj)
         ENDIF
      END DO   ;   END DO   ;   END DO

      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ! DOC ammonification. Depends on depth, phytoplankton biomass
         ! and a limitation term which is supposed to be a parameterization of the bacterial activity. 
         ! --------------------------------------------------------------------------
         zremik = xstep / 1.e-6 * xlimbac(mi,mj,jk) * zdepbac(mi,mj,jk)
         zremik = MAX( zremik, 2.74e-4 * xstep / xremikc )
         zremikc = xremikc * zremik
         ! Ammonification in oxic waters with oxygen consumption
         ! -----------------------------------------------------
         zolimic = zremikc * ( 1.- nitrfac(mi,mj,jk) ) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc)
         zolimic = MAX(0., MIN( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpoxy) - 0.5*EPSILON(1.e0) ) / o2ut, zolimic ) )

         ! Ammonification in suboxic waters with denitrification
         ! -----------------------------------------------------
         zammonic = zremikc * nitrfac(mi,mj,jk) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc)
         denitr(mi,mj,jk)  = zammonic * ( 1. - nitrfac2(mi,mj,jk) )
         denitr(mi,mj,jk)  = MAX(0., MIN(  ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) - 0.5*EPSILON(1.e0) ) &
              &            / rdenit, denitr(mi,mj,jk) ) )

         ! Ammonification in waters depleted in O2 and NO3 based on 
         ! other redox processes
         ! --------------------------------------------------------
         zoxyremc = MAX(0., zammonic - denitr(mi,mj,jk) )

         ! Update of the the trends arrays
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpno3) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpno3) - denitr (mi,mj,jk) * rdenit
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) &
            &                    - ( zolimic + denitr(mi,mj,jk) + zoxyremc )
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) &
            &                    - zolimic * (o2ut + o2nit)
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) &
            &                    + zolimic + denitr(mi,mj,jk) + zoxyremc
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpno3) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpno3) &
            &                    + zolimic + denitr(mi,mj,jk) + zoxyremc
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) - &
            &                      rno3 * ( zolimic + zoxyremc - &
            &                      ( rdenit - 1.) * denitr(mi,mj,jk) )
      END DO   ;   END DO   ;   END DO

      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem1')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
      !   CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF

      IF( .false. .AND. knt == nrdttrc ) THEN
          !
          IF( l_dia_remin ) THEN    ! Remineralisation rate
             ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp      
             !
             DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                zw3d(mi,mj,N+1-jk) = ( zolimi(mi,mj,jk) - t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) ) / o2ut &
                   &               * rfact2r * tmask(mi,mj,jk) ! 
             END DO   ;   END DO   ;   END DO
             CALL iom_put( "REMIN", zw3d )
             DEALLOCATE( zolimi, zw3d )
          ENDIF
          !
          IF( l_dia_bact ) THEN   ! Bacterial biomass
             ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp      
             DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                zw3d(mi,mj,N+1-jk) = zdepbac(mi,mj,jk) * 1.e+6 * tmask(mi,mj,jk)
             END DO   ;   END DO   ;   END DO
             CALL iom_put( "BACT", zw3d )
             DEALLOCATE( zw3d )
          ENDIF
          !
          IF( l_dia_denit )  THEN ! Denitrification
             ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp      
             DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                zw3d(mi,mj,N+1-jk) = denitr(mi,mj,jk) * 1.e+3 &
                      &          * rfact2r * rno3 * tmask(mi,mj,jk)
             END DO   ;   END DO   ;   END DO
             CALL iom_put( "DENIT", zw3d )
             DEALLOCATE( zw3d )
          ENDIF
          !
      ENDIF
      !
      IF( .false. )   CALL timing_stop('p2z_rem')
      !
   END SUBROUTINE p2z_rem

   SUBROUTINE p4z_rem( kt, knt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem  ***
      !!
      !! ** Purpose :   Compute remineralization/dissolution of organic compounds
      !!                Computes also nitrification of ammonium 
      !!                The solubilization/remineralization of POC is treated 
      !!                in p4zpoc.F90. The dissolution of calcite is processed
      !!                in p4zlys.F90. 
      !!
      !! ** Method  : - Bacterial biomass is computed implicitely based on a 
      !!                parameterization developed from an explicit modeling
      !!                of  in an alternative version 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt         ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      !
      INTEGER  ::   mi, mj, jk
      REAL(wp) ::   zremik, zremikc, zremikn, zremikp, zsiremin, zfact 
      REAL(wp) ::   zsatur, zsatur2, znusil, znusil2, zdep, zdepmin, zfactdep
      REAL(wp) ::   zbactfer, zonitr, zootot, zratio, zremtrd
      REAL(wp) ::   zammonic, zoxyremc, zosil, ztem, zdenitnh4, zolimic
      REAL(wp) ::   zfacsi, zdepeff
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N) :: zdepbac, zfacsib
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp    ) :: ztempbac
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zw3d, zolimi, zfebact
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p4z_rem')
      !
      IF( kt == nittrc000 )  THEN
         l_dia_remin  = iom_use( "REMIN" )  .OR. iom_use( "Remino2" )
         l_dia_bact   = iom_use( "FEBACT" ) .OR. iom_use( "BACT" )
         l_dia_denit  = iom_use( "DENIT" )
      ENDIF
      IF( l_dia_remin ) THEN
         ALLOCATE( zolimi(Istrp:Iendp,Jstrp:Jendp,N) )    ;   zolimi(Istrp:Iendp,Jstrp:Jendp,N) = 0._wp
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zolimi(mi,mj,jk) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy)
         END DO   ;   END DO   ;   END DO
      ENDIF
      IF( l_dia_bact ) THEN
         ALLOCATE( zfebact(Istrp:Iendp,Jstrp:Jendp,N) )   ;   zfebact(Istrp:Iendp,Jstrp:Jendp,N) = 0._wp
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zfebact(mi,mj,jk) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer)
         END DO   ;   END DO   ;   END DO
      ENDIF
      ! Initialisation of arrays
      zfacsib(:,:,:)  = xsilab / ( 1.0 - xsilab )

      ! Computation of the mean bacterial concentration
      ! this parameterization has been deduced from a model version
      ! that was modeling explicitely bacteria. This is a very old parame
      ! that will be very soon updated based on results from a much more
      ! recent version of  with bacteria.
      ! ----------------------------------------------------------------
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zdep = MAX( hbl(mi,mj), heup_01(mi,mj), ((-1)*(z_r(mi,mj,N+1-1)-z_w(mi,mj,N))) )
         IF ( ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) <= zdep ) THEN
            zootot = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpzoo) + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpmes)
            zdepbac(mi,mj,jk) = 0.6 * ( MAX(0.0, zootot ) * 1.0E6 )**0.6 * 1.E-6
            ztempbac(mi,mj)   = zdepbac(mi,mj,jk)
         ELSE
            zdepmin = MIN( 1., zdep / ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) )
            zdepbac(mi,mj,jk) = zdepmin**0.683 * ztempbac(mi,mj)
         ENDIF
      END DO   ;   END DO   ;   END DO

      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ! DOC ammonification. Depends on depth, phytoplankton biomass
         ! and a limitation term which is supposed to be a parameterization of the bacterial activity. 
         ! --------------------------------------------------------------------------
         zremik = xstep / 1.e-6 * xlimbac(mi,mj,jk) * zdepbac(mi,mj,jk) 
         zremik = MAX( zremik, 2.74e-4 * xstep / xremikc )
         zremikc = xremikc * zremik
         ! Ammonification in oxic waters with oxygen consumption
         ! -----------------------------------------------------
         zolimic = zremikc * ( 1.- nitrfac(mi,mj,jk) ) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) 
         zolimic = MAX(0., MIN( ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpoxy) - 0.5*EPSILON(1.e0) ) / o2ut, zolimic ) ) 

         ! Ammonification in suboxic waters with denitrification
         ! -----------------------------------------------------
         zammonic = zremikc * nitrfac(mi,mj,jk) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc)
         denitr(mi,mj,jk)  = zammonic * ( 1. - nitrfac2(mi,mj,jk) )
         denitr(mi,mj,jk)  = MAX(0., MIN(  ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) - 0.5*EPSILON(1.e0) ) &
             &             / rdenit, denitr(mi,mj,jk) ) )

         ! Ammonification in waters depleted in O2 and NO3 based on 
         ! other redox processes
         ! --------------------------------------------------------
         zoxyremc   = MAX(0., zammonic - denitr(mi,mj,jk) )
         zremtrd   =  zolimic + denitr(mi,mj,jk) + zoxyremc
         ! Update of the the trends arrays
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpno3) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpno3) - denitr (mi,mj,jk) * rdenit
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) - zremtrd
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) - zolimic * o2ut
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) + zremtrd
         IF( ln_p4z ) THEN ! -std
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) + zremtrd
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnh4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnh4) + zremtrd
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) &
                    &               + rno3 * ( zolimic + zoxyremc &
                    &                  + ( rdenit + 1.) * denitr(mi,mj,jk) )
         ELSE  ! -QUOTA (p5z)
            zratio = t(mi,mj,N+1-jk,kbb,itemp+ntrc_salt+jpdon) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) + 0.5*EPSILON(1.e0) )
            zremikn = xremikn / xremikc * zratio
            zratio = t(mi,mj,N+1-jk,kbb,itemp+ntrc_salt+jpdop) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) + 0.5*EPSILON(1.e0) )
            zremikp = xremikp / xremikc * zratio
            
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) + zremikp * zremtrd
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnh4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnh4) + zremikn * zremtrd
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdon) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdon) - zremikn * zremtrd
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdop) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdop) - zremikp * zremtrd
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) + rno3 * zremikn &
                 &             * ( zolimic + zoxyremc + ( rdenit + 1.) * denitr(mi,mj,jk) )
         ENDIF
      END DO   ;   END DO   ;   END DO

      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ! NH4 nitrification to NO3. Ceased for oxygen concentrations
         ! below 2 umol/L. Inhibited at strong light 
         ! ----------------------------------------------------------
         zonitr  = nitrif * xstep * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) * ( 1.- nitrfac(mi,mj,jk) )  &
         &         / ( 1.+ emoy(mi,mj,jk) ) * ( 1. + fr_i(mi,mj) * emoy(mi,mj,jk) ) 
         zdenitnh4 = nitrif * xstep * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) * nitrfac(mi,mj,jk)
         zdenitnh4 = MAX(0., MIN(  ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) - 0.5*EPSILON(1.e0) ) / rdenita, zdenitnh4 ) )
         ! Update of the tracers trends
         ! ----------------------------
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnh4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnh4) - zonitr - zdenitnh4
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpno3) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpno3) + zonitr - rdenita * zdenitnh4
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) - o2nit * zonitr
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) - 2 * rno3 * zonitr &
               &                  + rno3 * ( rdenita - 1. ) * zdenitnh4
      END DO   ;   END DO   ;   END DO

      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem1')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
      !   CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF

      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zdep = MAX( hbl(mi,mj), heup_01(mi,mj) )
         IF( ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) >= zdep ) THEN
            zdepmin = MIN( 1., zdep / ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) )
            zdepeff = 0.3_wp * zdepmin**0.6
         ELSE
            zdepeff = 0.3_wp
         ENDIF

         ! Bacterial uptake of iron. No iron is available in DOC. So
         ! Bacteries are obliged to take up iron from the water. Some
         ! studies (especially at Papa) have shown this uptake to be significant
         ! ----------------------------------------------------------
         zbactfer = feratb * 0.6_wp * xstep * tgfunc(mi,mj,jk) * xlimbacl(mi,mj,jk) * biron(mi,mj,jk)    &
           &       / ( xkferb + biron(mi,mj,jk) ) * zdepeff * zdepbac(mi,mj,jk)
         
         ! Only the transfer of iron from its dissolved form to particles
         ! is treated here. The GGE of bacteria supposed to be equal to 
         ! 0.33. This is hard-coded. 
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) - zbactfer*0.1
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) + zbactfer*0.08
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) + zbactfer*0.02
         blim(mi,mj,jk)          = xlimbacl(mi,mj,jk)  * zdepbac(mi,mj,jk) / 1.e-6
      END DO   ;   END DO   ;   END DO

       IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem2')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
     !    CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
       ENDIF

      ! Initialization of the array which contains the labile fraction
      ! of bSi. Set to a constant in the upper ocean
      ! ---------------------------------------------------------------
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ! Remineralization rate of BSi dependent on T and saturation
         ! The parameterization is taken from Ridgwell et al. (2002) 
         ! ---------------------------------------------------------
         zdep     = MAX( hbl(mi,mj), heup_01(mi,mj) )
         zsatur   = MAX( 0.5*EPSILON(1.e0), ( sio3eq(mi,mj,jk) - t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsil) ) &
             &     / ( sio3eq(mi,mj,jk) + 0.5*EPSILON(1.e0) ) )
         zsatur2  = ( 1. + t(mi,mj,N+1-jk,nnew,itemp) / 400.)**37
         znusil   = 0.225  * ( 1. + t(mi,mj,N+1-jk,nnew,itemp) / 15.) &
             &     * zsatur + 0.775 * zsatur2 * zsatur**9.25
 
         ! Two fractions of bSi are considered : a labile one and a more
         ! refractory one based on the commonly observed two step 
         ! dissolution of bSi (initial rapid dissolution followed by 
         ! more slowly dissolution).
         ! Computation of the vertical evolution of the labile fraction
         ! of bSi. This is computed assuming steady state.
         ! --------------------------------------------------------------
         IF ( ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) > zdep ) THEN
            zfactdep = EXP( -0.5 * ( xsiremlab - xsirem ) &
                &     * znusil * Hz(mi,mj,N+1-jk) / wsbio4(mi,mj,jk) )
            zfacsib(mi,mj,jk) = zfacsib(mi,mj,jk-1) * zfactdep
            zfacsi            = zfacsib(mi,mj,jk) / ( 1.0 + zfacsib(mi,mj,jk) )
            zfacsib(mi,mj,jk) = zfacsib(mi,mj,jk) * zfactdep
         ELSE
            zfacsi  = xsilab
         ENDIF
         zsiremin = ( xsiremlab * zfacsi + xsirem * ( 1. - zfacsi ) ) * xstep * znusil
         zosil    = zsiremin * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgsi)
         !
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgsi) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpgsi) - zosil
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsil) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsil) + zosil
      END DO   ;   END DO   ;   END DO

      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem3')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
       !  CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF

      IF( .false. .AND. knt == nrdttrc ) THEN
          !
          IF( l_dia_remin ) THEN    ! Remineralisation rate
             ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
             !
             DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                zw3d(mi,mj,N+1-jk) = ( zolimi(mi,mj,jk) - t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) )  &
                   &              / o2ut  * rfact2r * tmask(mi,mj,jk) !
             END DO   ;   END DO   ;   END DO
             CALL iom_put( "REMIN", zw3d )
             DEALLOCATE( zolimi, zw3d )
          ENDIF
          !
          IF( l_dia_bact ) THEN   ! Bacterial biomass
             ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
             DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                zw3d(mi,mj,N+1-jk) = zdepbac(mi,mj,jk) * 1.e+6 * tmask(mi,mj,jk)
             END DO   ;   END DO   ;   END DO
             CALL iom_put( "BACT", zw3d )
             !
             DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                zw3d(mi,mj,N+1-jk) = ( zfebact(mi,mj,jk) - t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) ) &
                   &              * 1e9 * rfact2r * tmask(mi,mj,jk) ! conversion in nmol/m2/s
             END DO   ;   END DO   ;   END DO
             CALL iom_put( "FEBACT", zw3d )
             DEALLOCATE( zfebact, zw3d )
          ENDIF
          !
          IF( l_dia_denit )  THEN ! Denitrification
             ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
             DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                zw3d(mi,mj,N+1-jk) = denitr(mi,mj,jk) * 1.e+3 &
                      &          * rfact2r * rno3 * tmask(mi,mj,jk)
             END DO   ;   END DO   ;   END DO
             CALL iom_put( "DENIT", zw3d )
             DEALLOCATE( zw3d )
          ENDIF
          !
      ENDIF          
      !
      IF( .false. )   CALL timing_stop('p4z_rem')
      !
   END SUBROUTINE p4z_rem


   SUBROUTINE p4z_rem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_rem_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampisrem namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisrem
      !!
      !!----------------------------------------------------------------------
      NAMELIST/nampisrem/ nitrif, xsirem, xsiremlab, xsilab, feratb, xkferb, & 
         &                xremikc, xremikn, xremikp
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*)
         WRITE(stdout,*) 'p4z_rem_init : Initialization of remineralization parameters'
         WRITE(stdout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND(numnatp_ref);READ(numnatp_ref,nampisrem,IOSTAT=ios);CALL ctl_nam(ios,"nampisrem (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,nampisrem,IOSTAT=ios);CALL ctl_nam(ios,"nampisrem (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE( numonp, nampisrem )

      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*) '   Namelist parameters for remineralization, nampisrem'
         WRITE(stdout,*) '      remineralization rate of DOC              xremikc   =', xremikc
         IF( ln_p5z ) THEN 
            WRITE(stdout,*) '      remineralization rate of DOC              xremikc   =', xremikc
            WRITE(stdout,*) '      remineralization rate of DON              xremikn   =', xremikn
            WRITE(stdout,*) '      remineralization rate of DOP              xremikp   =', xremikp
         ENDIF
         IF( ln_p5z .OR. ln_p4z ) THEN
            WRITE(stdout,*) '      remineralization rate of Si               xsirem    =', xsirem
            WRITE(stdout,*) '      fast remineralization rate of Si          xsiremlab =', xsiremlab
            WRITE(stdout,*) '      fraction of labile biogenic silica        xsilab    =', xsilab
            WRITE(stdout,*) '      NH4 nitrification rate                    nitrif    =', nitrif
            WRITE(stdout,*) '      Bacterial Fe/C ratio                      feratb    =', feratb
            WRITE(stdout,*) '      Half-saturation constant for bact. Fe/C   xkferb    =', xkferb
         ENDIF
      ENDIF
      !
      denitr(:,:,N) = 0._wp
      blim  (:,:,N) = 0._wp
      !
   END SUBROUTINE p4z_rem_init


   INTEGER FUNCTION p4z_rem_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( denitr(Istrp:Iendp,Jstrp:Jendp,N), STAT=p4z_rem_alloc )
      !
      IF( p4z_rem_alloc /= 0 )   CALL ctl_stop( 'STOP', 'p4z_rem_alloc: failed to allocate arrays' )
      !
   END FUNCTION p4z_rem_alloc


   !!======================================================================
END MODULE p4zrem
