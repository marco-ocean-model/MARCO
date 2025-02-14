










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









MODULE p4zfechem
   !!======================================================================
   !!                         ***  MODULE p4zfechem  ***
   !! TOP :    Compute iron chemistry and scavenging
   !!======================================================================
   !! History :   3.5  !  2012-07 (O. Aumont, A. Tagliabue, C. Ethe) Original code
   !!             3.6  !  2015-05  (O. Aumont)  quota
   !!             3.*  !  2025-02  (S. Maishal, R. Person)  and optimization
   !!----------------------------------------------------------------------
   !!   p4z_fechem       : Compute remineralization/scavenging of iron
   !!   p4z_fechem_init  : Initialisation of parameters for remineralisation
   !!   p4z_fechem_alloc : Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE p4zche          ! chemical model
   USE p4zbc           ! Boundary conditions from sediments
   USE prtctl          ! print control for debugging
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_fechem        ! called in p4zbio.F90
   PUBLIC   p4z_fechem_init   ! called in trcsms_pisces.F90

   LOGICAL          ::   ln_ligvar    !: boolean for variable ligand concentration following Tagliabue and voelker
   REAL(wp), PUBLIC ::   xlam1        !: scavenging rate of Iron 
   REAL(wp), PUBLIC ::   xlamdust     !: scavenging rate of Iron by dust 
   REAL(wp), PUBLIC ::   ligand       !: ligand concentration in the ocean 
   REAL(wp), PUBLIC ::   kfep         !: rate constant for nanoparticle formation
   REAL(wp), PUBLIC ::   scaveff      !: Fraction of scavenged iron that is considered as being subject to solubilization

   LOGICAL  :: l_dia_fechem
   REAL(wp) :: xpow = 5.0118e-7   ! 10**(-6.3)

   !! * Substitutions











































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zfechem.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_fechem( kt, knt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_fechem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of iron
      !!
      !! ** Method  :   A simple chemistry model of iron from Aumont and Bopp (2006)
      !!                based on one ligand and one inorganic form
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      !
      INTEGER  ::   mi, mj, jk, jic, jn
      REAL(wp) ::   zlam1a, zlam1b
      REAL(wp) ::   zkeq, zfesatur, fe3sol, zligco
      REAL(wp) ::   zscave, zaggdfea, zaggdfeb, ztrc, zdust, zklight
      REAL(wp) ::   ztfe, zhplus, zxlam, zaggliga, zaggligb
      REAL(wp) ::   zprecip, zprecipno3,  zconsfe, za1, ztl1, zfel1
      REAL(wp) ::   zrfact2
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N) ::   zFe3, ztotlig,  zfecoll
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zcoll3d, zscav3d, zfeprecip, zw3d
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p4z_fechem')
      !
      IF( kt == nittrc000 )  &
         l_dia_fechem  = iom_use( "Fe3" ) .OR. iom_use( "FeL1" ) .OR. iom_use( "TL1" ) .OR.  &
            &            iom_use( "Totlig" ) .OR. iom_use( "Biron" ) .OR. iom_use( "FESCAV" ) .OR.  &
            &            iom_use( "FECOLL" ) .OR. iom_use( "FEPREC" ) 

      IF( l_dia_fechem )  &
        & ALLOCATE( zcoll3d(Istrp:Iendp,Jstrp:Jendp,N), zscav3d(Istrp:Iendp,Jstrp:Jendp,N), zfeprecip(Istrp:Iendp,Jstrp:Jendp,N) ) 
      !
      ! Total ligand concentration : Ligands can be chosen to be constant or variable
      ! Parameterization from Pham and Ito (2018)
      ! -------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, xfecolagg) &
      !$OMP SHARED(ligand, chemo2, tr, jpoxy, Kbb)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         xfecolagg(mi,mj,jk) = ligand * 1E9 + 0.01 &
                 &  * MAX(0., (chemo2(mi,mj,jk) - t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpoxy) ) * 1E6 )**0.8
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      IF( ln_ligvar ) THEN
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            ztotlig(mi,mj,jk) =  0.07 * 0.667 * (t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) * 1E6 )**0.8  &
                    &   + xfecolagg(mi,mj,jk)
            ztotlig(mi,mj,jk) =  MIN( ztotlig(mi,mj,jk), 10. )
         END DO   ;   END DO   ;   END DO
      ELSE
        IF( ln_ligand ) THEN
           !$OMP PARALLEL DO &
           !$OMP PRIVATE(mi, mj, jk, ztotlig) &
           !$OMP SHARED(tr, jpdoc, Kbb, xfecolagg)
           DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
              ztotlig(mi,mj,jk) = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jplgw) * 1E9
           END DO   ;   END DO   ;   END DO
           !$OMP END PARALLEL DO
        ELSE
             ztotlig(:,:,:) = ligand * 1E9 
        ENDIF
      ENDIF

      ! ------------------------------------------------------------
      ! From Aumont and Bopp (2006)
      ! This model is based on one ligand, Fe2+ and Fe3+ 
      ! Chemistry is supposed to be fast enough to be at equilibrium
      ! ------------------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, ztl1, zkeq, zklight, zconsfe, zfesatur, ztfe, za1) &
      !$OMP SHARED(fekeq, etot, xpow, consfe3, tr, jpfer, Kbb, 0.5*EPSILON(1.e0), zFe3)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
          ztl1            = ztotlig(mi,mj,jk)
          zkeq            = fekeq(mi,mj,jk)
          zklight         = 4.77E-7 * etot(mi,mj,jk) * 0.5 / xpow
          zconsfe         = consfe3(mi,mj,jk) / xpow
          zfesatur        = ztl1 * 1E-9
          ztfe            = (1.0 + zklight) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer) 
          ! Fe' is the root of a 2nd order polynomial
          za1 =  1. + zfesatur * zkeq + zklight +  zconsfe - zkeq * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer)
          zFe3 (mi,mj,jk) = ( -1 * za1 + SQRT( za1**2 + 4. * ztfe * zkeq) ) / ( 2. * zkeq + 0.5*EPSILON(1.e0) )
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zfel1) &
      !$OMP SHARED(tr, jpfer, Kbb, 0.5*EPSILON(1.e0), zFe3, plig)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zfel1 = MAX( 0., t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer) - zFe3(mi,mj,jk) )
         plig(mi,mj,jk) =  MAX( 0., ( zfel1 / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer) + 0.5*EPSILON(1.e0) ) ) )
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      zdust = 0.         ! if no dust input available

      ! Computation of the colloidal fraction that is subjecto to coagulation
      ! The assumption is that 50% of complexed iron is colloidal. Furthermore
      ! The refractory part is supposed to be non sticky. The refractory
      ! fraction is supposed to equal to the background concentration + 
      ! the fraction that accumulates in the deep ocean. AOU is taken as a 
      ! proxy of that accumulation following numerous studies showing 
      ! some relationship between weak ligands and AOU.
      ! An issue with that parameterization is that when ligands are not
      ! prognostic or non variable, all the colloidal fraction is supposed
      ! to coagulate
      ! ----------------------------------------------------------------------
      IF (ln_ligand) THEN
         !$OMP PARALLEL DO &
         !$OMP PRIVATE(mi, mj, jk, zfel1) &
         !$OMP SHARED(tr, jpfer, Kbb, zFe3, zfecoll, ztotlig, xfecolagg, 0.5*EPSILON(1.e0))
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zfel1 = MAX( 0., t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer) - zFe3(mi,mj,jk) )
            zfecoll(mi,mj,jk) = 0.5 * zfel1 * MAX(0., ztotlig(mi,mj,jk) - xfecolagg(mi,mj,jk) ) &
                  &              / ( ztotlig(mi,mj,jk) + 0.5*EPSILON(1.e0) ) 
         END DO   ;   END DO   ;   END DO
         !$OMP END PARALLEL DO
      ELSE
         IF (ln_ligvar) THEN
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(mi, mj, jk, zfel1) &
            !$OMP SHARED(tr, jpfer, Kbb, zFe3, zfecoll, ztotlig, xfecolagg, 0.5*EPSILON(1.e0))
            DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zfel1 = MAX( 0., t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer) - zFe3(mi,mj,jk) )
               zfecoll(mi,mj,jk) = 0.5 * zfel1 * MAX(0., ztotlig(mi,mj,jk) - xfecolagg(mi,mj,jk) ) &
                  &              / ( ztotlig(mi,mj,jk) + 0.5*EPSILON(1.e0) ) 
            END DO   ;   END DO   ;   END DO
            !$OMP END PARALLEL DO
         ELSE
            !$OMP PARALLEL DO &
            !$OMP PRIVATE(mi, mj, jk, zfel1) &
            !$OMP SHARED(tr, jpfer, Kbb, zFe3, zfecoll)
            DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zfel1 = MAX( 0., t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer) - zFe3(mi,mj,jk) )
               zfecoll(mi,mj,jk) = 0.5 * zfel1
            END DO   ;   END DO   ;   END DO
            !$OMP END PARALLEL DO
         ENDIF
      ENDIF
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zhplus, fe3sol, zprecip, zprecipno3, zlam1a, &
                    zaggdfea, zdust, zxlam, ztrc, zlam1b, zscave, zaggdfeb, &
                    xcoagfe, zscav3d, zcoll3d, zfeprecip) &
      !$OMP SHARED(tr, jpfer, Kbb, zFe3, fesol, 0.5*EPSILON(1.e0), hi, kfep, nitrfac, &
                   xstep, xpow, xlamdust, tmask, ln_p2z, day2sec, scaveff, &
                   l_dia_fechem, xlam1, xdiss, tr)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ! Scavenging rate of iron. This scavenging rate depends on the load of particles of sea water. 
         ! This parameterization assumes a simple second order kinetics (k[Particles][Fe]).
         ! Scavenging onto dust is also included as evidenced from the DUNE experiments.
         ! --------------------------------------------------------------------------------------
         zhplus  = max( 0.5*EPSILON(1.e0), hi(mi,mj,jk) )
         fe3sol  = fesol(mi,mj,jk,1) * ( zhplus**3 + fesol(mi,mj,jk,2) * zhplus**2  &
         &         + fesol(mi,mj,jk,3) * zhplus + fesol(mi,mj,jk,4)     &
         &         + fesol(mi,mj,jk,5) / zhplus )
         !
         ! precipitation of Fe3+, creation of nanoparticles
         zprecip = MAX( 0., ( zFe3(mi,mj,jk) - fe3sol ) ) * kfep * xstep * ( 1.0 - nitrfac(mi,mj,jk) ) 
         ! Precipitation of Fe2+ due to oxidation by NO3 (Croot et al., 2019)
         ! This occurs in anoxic waters only
         zprecipno3 = 2.0 * 130.0 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) * nitrfac(mi,mj,jk) * xstep * zFe3(mi,mj,jk)
         !
         !  Compute the coagulation of colloidal iron. This parameterization 
         !  could be thought as an equivalent of colloidal pumping.
         !  It requires certainly some more work as it is very poorly constrained.
         !  ----------------------------------------------------------------
         zlam1a   = ( 12.0  * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) &
             &        + 9.05  * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) ) * xdiss(mi,mj,jk)    &
             &    + ( 2.49  * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) )     &
             &    + ( 127.8 * 0.3 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) &
             &         + 725.7 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) )
         zaggdfea = zlam1a * xstep * zfecoll(mi,mj,jk)
         !
         IF( ll_dust )  zdust  = dust(mi,mj) / ( wdust / day2sec ) * tmask(mi,mj,jk)
         zxlam  = MAX( 1.E-3, (1. - EXP(-2 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpoxy) / 100.E-6 ) ))

         IF( ln_p2z ) THEN
            ztrc = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) * 1e6
         ELSE
            ztrc = ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc) &
               &  + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpcal) + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgsi) ) * 1.e6
         ENDIF
         ztrc = MAX( 0.5*EPSILON(1.e0), ztrc )
         zlam1b = 3.e-5 + ( xlamdust * zdust + xlam1 * ztrc ) * zxlam
         zscave = zFe3(mi,mj,jk) * zlam1b * xstep

         !
         IF( ln_p2z ) THEN
            zaggdfeb = 0._wp
            xcoagfe(mi,mj,jk) = zlam1a
         ELSE
            zlam1b   = ( 1.94 * xdiss(mi,mj,jk) + 1.37 ) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpgoc)
            zaggdfeb = zlam1b * xstep * zfecoll(mi,mj,jk)
            xcoagfe(mi,mj,jk) =  zlam1a + zlam1b
            !
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) &
                    &     + zscave * scaveff * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) / ztrc
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) &
                    &     + zscave * scaveff * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppoc) / ztrc
            !
            ! Precipitated iron is supposed to be permanently lost.
            ! Scavenged iron is supposed to be released back to seawater
            ! when POM is solubilized. This is highly uncertain as probably
            ! a significant part of it may be rescavenged back onto 
            ! the particles. An efficiency factor is applied that is read
            ! in the namelist. 
            ! See for instance Tagliabue et al. (2019).
            ! Aggregated FeL is considered as biogenic Fe as it 
            ! probably remains  complexed when the particle is solubilized.
            ! -------------------------------------------------------------
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsfe) + zaggdfea
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpbfe) + zaggdfeb
            !
         ENDIF
         t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) - zscave - zaggdfea - zaggdfeb &
            &                    - ( zprecip + zprecipno3 )

         IF( l_dia_fechem ) THEN
            zscav3d(mi,mj,jk)   = zscave 
            zcoll3d(mi,mj,jk)   = zaggdfea + zaggdfeb
            zfeprecip(mi,mj,jk) = zprecip + zprecipno3
         ENDIF
         !
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      !  Define the bioavailable fraction of iron
      !  ----------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(tr, jpfer, Kbb, biron)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         biron(mi,mj,jk) = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer) 
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      !  Output of some diagnostics variables
      !     ---------------------------------
      IF( l_dia_fechem .AND. .false. .AND. knt == nrdttrc ) THEN
        !
        zrfact2 = 1.e3 * rfact2r  ! conversion from mol/L/timestep into mol/m3/s
        ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
        ! Fe3+
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(mi, mj, N+1-jk) &
        !$OMP SHARED(zFe3, tmask, zw3d)
        DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           zw3d(mi,mj,N+1-jk) = zFe3(mi,mj,jk) * tmask(mi,mj,jk)
        END DO   ;   END DO   ;   END DO
        !$OMP END PARALLEL DO
        CALL iom_put( "Fe3", zw3d )
        !  FeL1
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(mi, mj, N+1-jk) &
        !$OMP SHARED(tr, zFe3, tmask, zw3d)
        DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
          zw3d(mi,mj,N+1-jk) = MAX( 0., t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer) - zFe3(mi,mj,jk) ) &
                  &        * tmask(mi,mj,jk)
        END DO   ;   END DO   ;   END DO
        !$OMP END PARALLEL DO
        CALL iom_put( "FeL1", zw3d )
        ! TL1 = Totlig
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(mi, mj, N+1-jk) &
        !$OMP SHARED(ztotlig, tmask, zw3d)
        DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           zw3d(mi,mj,N+1-jk) = ztotlig(mi,mj,jk) * tmask(mi,mj,jk)
        END DO   ;   END DO   ;   END DO
        !$OMP END PARALLEL DO
        CALL iom_put( "TL1", zw3d )
        ! Totlig
        CALL iom_put( "Totlig", zw3d )
        ! biron
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(mi, mj, N+1-jk) &
        !$OMP SHARED(biron, tmask, zw3d)
        DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           zw3d(mi,mj,N+1-jk) = biron(mi,mj,jk) * tmask(mi,mj,jk)
        END DO   ;   END DO   ;   END DO
        !$OMP END PARALLEL DO
        CALL iom_put( "Biron", zw3d )
        ! FESCAV
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(mi, mj, N+1-jk) &
        !$OMP SHARED(zscav3d, tmask, zw3d, zrfact2)
        DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           zw3d(mi,mj,N+1-jk) = zscav3d(mi,mj,jk) * tmask(mi,mj,jk) * zrfact2
        END DO   ;   END DO   ;   END DO
        !$OMP END PARALLEL DO
        CALL iom_put( "FESCAV", zw3d )
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(mi, mj, N+1-jk) &
        !$OMP SHARED(zcoll3d, tmask, zw3d, zrfact2)
        DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           zw3d(mi,mj,N+1-jk) = zcoll3d(mi,mj,jk) * tmask(mi,mj,jk) * zrfact2
        END DO   ;   END DO   ;   END DO
        !$OMP END PARALLEL DO
        ! FECOLL
        CALL iom_put( "FECOLL", zw3d )
        ! FEPREC
        !$OMP PARALLEL DO &
        !$OMP PRIVATE(mi, mj, N+1-jk) &
        !$OMP SHARED(zfeprecip, tmask, zw3d, zrfact2)
        DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           zw3d(mi,mj,N+1-jk) = zfeprecip(mi,mj,jk) * tmask(mi,mj,jk) * zrfact2
        END DO   ;   END DO   ;   END DO
        !$OMP END PARALLEL DO
        CALL iom_put( "FEPREC", zw3d )
        !
        DEALLOCATE( zcoll3d, zscav3d, zfeprecip, zw3d )
        !
      ENDIF
      !
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('fechem')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
 !        CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( .false. )   CALL timing_stop('p4z_fechem')
      !
   END SUBROUTINE p4z_fechem


   SUBROUTINE p4z_fechem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_fechem_init  ***
      !!
      !! ** Purpose :   Initialization of iron chemistry parameters
      !!
      !! ** Method  :   Read the nampisfer namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisfer
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer 
      !!
      NAMELIST/nampisfer/ ln_ligvar, xlam1, xlamdust, ligand, kfep, scaveff 
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*)
         WRITE(stdout,*) 'p4z_rem_init : Initialization of iron chemistry parameters'
         WRITE(stdout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND(numnatp_ref);READ(numnatp_ref,nampisfer,IOSTAT=ios);CALL ctl_nam(ios,"nampisfer (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,nampisfer,IOSTAT=ios);CALL ctl_nam(ios,"nampisfer (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE( numonp, nampisfer )

      IF(mynode .eq. 0) THEN                     ! control print
         WRITE(stdout,*) '   Namelist : nampisfer'
         WRITE(stdout,*) '      variable concentration of ligand          ln_ligvar    =', ln_ligvar
         WRITE(stdout,*) '      scavenging rate of Iron                   xlam1        =', xlam1
         WRITE(stdout,*) '      scavenging rate of Iron by dust           xlamdust     =', xlamdust
         WRITE(stdout,*) '      ligand concentration in the ocean         ligand       =', ligand
         WRITE(stdout,*) '      rate constant for nanoparticle formation  kfep         =', kfep
         WRITE(stdout,*) '      Scavenged iron that is added to POFe      scaveff      =', scaveff
      ENDIF
      !
   END SUBROUTINE p4z_fechem_init

   
   !!======================================================================
END MODULE p4zfechem
