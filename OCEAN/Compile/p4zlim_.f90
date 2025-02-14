










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









MODULE p4zlim
   !!======================================================================
   !!                         ***  MODULE p4zlim  ***
   !! TOP :   Computes the nutrient limitation terms of phytoplankton
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-04  (O. Aumont, C. Ethe) Limitation for iron modelled in quota
   !!             3.*  !  2025-02  (S. Maishal, R. Person)  and optimization
   !!----------------------------------------------------------------------
   !!   p4z_lim        :   Compute the nutrients limitation terms 
   !!   p4z_lim_init   :   Read the namelist 
   !!----------------------------------------------------------------------

   USE oce_trc          ! Shared ocean-passive tracers variables
   USE trc              ! Tracers defined
   USE sms_pisces       !  variables
   USE p2zlim           ! Reduced  nutrient limitation
   USE iom              ! I/O manager
   IMPLICIT NONE
   PRIVATE
   PUBLIC p4z_lim       ! called in p4zbio.F90 
   PUBLIC p4z_lim_init  ! called in trcsms_pisces.F90 
   PUBLIC p4z_lim_alloc ! called in trcini_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  concdno3    !:  Phosphate half saturation for diatoms  
   REAL(wp), PUBLIC ::  concnnh4    !:  NH4 half saturation for nanophyto  
   REAL(wp), PUBLIC ::  concdnh4    !:  NH4 half saturation for diatoms
   REAL(wp), PUBLIC ::  concdfer    !:  Iron half saturation for diatoms  
   REAL(wp), PUBLIC ::  concbnh4    !:  NH4 half saturation for bacteria
   REAL(wp), PUBLIC ::  xsizedia    !:  Minimum size criteria for diatoms
   REAL(wp), PUBLIC ::  xsizerd     !:  Size ratio for diatoms
   REAL(wp), PUBLIC ::  xksi1       !:  half saturation constant for Si uptake 
   REAL(wp), PUBLIC ::  xksi2       !:  half saturation constant for Si/C 
   REAL(wp), PUBLIC ::  qnfelim     !:  optimal Fe quota for nanophyto
   REAL(wp), PUBLIC ::  qdfelim     !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  ratchl      !:  C associated with Chlorophyll

   !!* Phytoplankton limitation terms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatno3   !: Diatoms limitation by NO3
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanonh4   !: Nanophyto limitation by NH4
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatnh4   !: Diatoms limitation by NH4
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanopo4   !: Nanophyto limitation by PO4
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatpo4   !: Diatoms limitation by PO4
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdia    !: Nutrient limitation term of diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdfe    !: Diatoms limitation by iron
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimsi     !: Diatoms limitation by Si
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concdfe    !: Limitation of diatoms uptake of Fe
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   concnfe    !: Limitation of Nano uptake of Fe
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanofer   !: Limitation of Fe uptake by nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatfer   !: Limitation of Fe uptake by diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xqfuncfecd, xqfuncfecn

   ! Coefficient for iron limitation following Flynn and Hipkin (1999)
   REAL(wp) ::  xcoef1   = 0.0016  / 55.85  
   REAL(wp) ::  xcoef2   = 1.21E-5 * 14. / 55.85 / 7.3125 * 0.5 * 1.5
   REAL(wp) ::  xcoef3   = 1.15E-4 * 14. / 55.85 / 7.3125 * 0.5 

   LOGICAL  :: l_dia_nut_lim, l_dia_iron_lim, l_dia_size_lim, l_dia_fracal

   !! * Substitutions












































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zlim.F90 10069 2018-08-28 14:12:24Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_lim( kt, knt, Kbb, Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_lim  ***
      !!
      !! ** Purpose :   Compute the co-limitations by the various nutrients
      !!                for the various phytoplankton species
      !!
      !! ** Method  : - Limitation follows the Liebieg law of the minimum
      !!              - Monod approach for N, P and Si. Quota approach 
      !!                for Iron
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in)  :: kt, knt
      INTEGER, INTENT(in)  :: Kbb, Kmm      ! time level indices
      !
      INTEGER  ::   mi, mj, jk
      REAL(wp) ::   zlim1, zlim2, zlim3, zlim4, zcoef
      REAL(wp) ::   z1_trbdia, z1_trbphy, ztem1, ztem2, zetot1, zetot2
      REAL(wp) ::   zdenom, zratio, zironmin, zbactno3, zbactnh4
      REAL(wp) ::   zconc1d, zconc1dnh4, zconc0n, zconc0nnh4   
      REAL(wp) ::   fananof, fadiatf, znutlim, zfalim
      REAL(wp) ::   znutlimtot, zlimno3, zlimnh4, zbiron
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p4z_lim')
      !
      IF( kt == nittrc000 )  THEN
         l_dia_nut_lim  = iom_use( "LNnut"   ) .OR. iom_use( "LDnut" )  
         l_dia_iron_lim = iom_use( "LNFe"    ) .OR. iom_use( "LDFe"  )
         l_dia_size_lim = iom_use( "SIZEN"   ) .OR. iom_use( "SIZED" )
         l_dia_fracal   = iom_use( "xfracal" )
      ENDIF
      !
      sizena(:,:,:) = 1.0  ;  sizeda(:,:,:) = 1.0
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, z1_trbphy, z1_trbdia, concnfe, zconc0n, zconc0nnh4, &
                    concdfe, zconc1d, zconc1dnh4, zbiron, znutlim, fananof, fadiatf, &
                    zlimnh4, zlimno3, znutlimtot, zbactnh4, zbactno3, zlim1, zlim2, &
                    zlim3, zlim4, zratio, zironmin, zfalim) &
      !$OMP SHARED(tr, jpphy, jpdia, Kbb, concnfer, concnno3, concnnh4, concdfer, &
                   concdno3, concdnh4, sizen, sized, plig, biron, concbno3, concbnh4, &
                   concbfe, xkdoc, xksi, xcoef1, xcoef2, xcoef3, qnfelim, qdfelim, &
                   xnanonh4, xnanono3, xnanopo4, xlimnfe, xlimphy, xlimbac, xlimbacl, &
                   xdiatnh4, xdiatno3, xdiatpo4, xlimdfe, xlimdia, xlimsi, &
                   xnanofer, xdiatfer)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
          ! Computation of a variable Ks for iron on diatoms taking into account
          ! that increasing biomass is made of generally bigger cells
          ! The allometric relationship is classical.
          !------------------------------------------------
          z1_trbphy   = 1. / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) + 0.5*EPSILON(1.e0) )
          z1_trbdia   = 1. / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) + 0.5*EPSILON(1.e0) )

          concnfe(mi,mj,jk) = concnfer * sizen(mi,mj,jk)**0.81
          zconc0n           = concnno3 * sizen(mi,mj,jk)**0.81
          zconc0nnh4        = concnnh4 * sizen(mi,mj,jk)**0.81

          concdfe(mi,mj,jk) = concdfer * sized(mi,mj,jk)**0.81 
          zconc1d           = concdno3 * sized(mi,mj,jk)**0.81 
          zconc1dnh4        = concdnh4 * sized(mi,mj,jk)**0.81  

          ! Computation of the optimal allocation parameters
          ! Based on the different papers by Pahlow et al., and 
          ! Smith et al.
          ! ---------------------------------------------------

          ! Nanophytoplankton
          zbiron = ( 75.0 * ( 1.0 - plig(mi,mj,jk) ) + plig(mi,mj,jk) ) * biron(mi,mj,jk)
          znutlim = zbiron / concnfe(mi,mj,jk)
          fananof = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )

          ! Diatoms
          znutlim = zbiron / concdfe(mi,mj,jk)
          fadiatf = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )

          ! Michaelis-Menten Limitation term by nutrients of
          ! heterotrophic bacteria
          ! -------------------------------------------------
          zlimnh4 = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) / ( concbno3 + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) )
          zlimno3 = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) / ( concbno3 + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) )
          znutlimtot = ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) ) &
             &       / ( concbno3 + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) )
          zbactnh4 = znutlimtot * 5.0 * zlimnh4 / ( zlimno3 + 5.0 * zlimnh4 + 0.5*EPSILON(1.e0) )
          zbactno3 = znutlimtot * zlimno3 / ( zlimno3 + 5.0 * zlimnh4 + 0.5*EPSILON(1.e0) )
          !
          zlim1    = zbactno3 + zbactnh4
          zlim2    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) + concbnh4 )
          zlim3    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer) / ( concbfe + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer) )
          zlim4    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) / ( xkdoc   + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) )
          ! Xlimbac is used for DOC solubilization whereas xlimbacl
          ! is used for all the other bacterial-dependent terms
          ! -------------------------------------------------------
          xlimbacl(mi,mj,jk) = MIN( zlim1, zlim2, zlim3 )
          xlimbac (mi,mj,jk) = MIN( zlim1, zlim2, zlim3 ) * zlim4

          ! Michaelis-Menten Limitation term by nutrients: Nanophyto
          ! Optimal parameterization by Smith and Pahlow series of 
          ! papers is used. Optimal allocation is supposed independant
          ! for all nutrients. 
          ! --------------------------------------------------------

          ! Limitation of Fe uptake (Quota formalism)
          zfalim = (1.-fananof) / fananof
          xnanofer(mi,mj,jk) = (1. - fananof) * zbiron / ( zbiron + zfalim * concnfe(mi,mj,jk) )

          ! Limitation of nanophytoplankton growth
          zlimnh4 = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) / ( zconc0n + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) )
          zlimno3 = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) / ( zconc0n + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) )
          znutlimtot = ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) ) &
              &      / ( zconc0n + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) )
          xnanonh4(mi,mj,jk) = znutlimtot * 5.0 * zlimnh4 / ( zlimno3 + 5.0 * zlimnh4 + 0.5*EPSILON(1.e0) )
          xnanono3(mi,mj,jk) = znutlimtot * zlimno3 / ( zlimno3 + 5.0 * zlimnh4 + 0.5*EPSILON(1.e0) )
          !
          zlim1    = xnanono3(mi,mj,jk) + xnanonh4(mi,mj,jk)
          zlim2    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) + zconc0nnh4 )
          zratio   = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnfe) * z1_trbphy 

          ! The minimum iron quota depends on the size of PSU, respiration
          ! and the reduction of nitrate following the parameterization 
          ! proposed by Flynn and Hipkin (1999)
          zironmin = xcoef1 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnch) * z1_trbphy + xcoef2 * zlim1 + xcoef3 * xnanono3(mi,mj,jk)
          xqfuncfecn(mi,mj,jk) = zironmin + qnfelim
          zlim3    = MAX( 0.,( zratio - zironmin ) / qnfelim )
          xnanopo4(mi,mj,jk) = zlim2
          xlimnfe (mi,mj,jk) = MIN( 1., zlim3 )
          xlimphy (mi,mj,jk) = MIN( zlim1, zlim2, zlim3 )
               
          !   Michaelis-Menten Limitation term by nutrients : Diatoms
          !   -------------------------------------------------------
          ! Limitation of Fe uptake (Quota formalism)
          zfalim = (1.-fadiatf) / fadiatf
          xdiatfer(mi,mj,jk) = (1. - fadiatf) * zbiron / ( zbiron + zfalim * concdfe(mi,mj,jk) )

          ! Limitation of diatoms growth
          zlimnh4 = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) / ( zconc1d + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) )
          zlimno3 = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) / ( zconc1d + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) )
          znutlimtot = ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) ) &
              &     / ( zconc1d + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) )
          xdiatnh4(mi,mj,jk) = znutlimtot * 5.0 * zlimnh4 / ( zlimno3 + 5.0 * zlimnh4 + 0.5*EPSILON(1.e0) ) 
          xdiatno3(mi,mj,jk) = znutlimtot * zlimno3 / ( zlimno3 + 5.0 * zlimnh4 + 0.5*EPSILON(1.e0) )
          !
          zlim1    = xdiatno3(mi,mj,jk) + xdiatnh4(mi,mj,jk)
          zlim2    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) + zconc1dnh4  )
          zlim3    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsil) &
               &   / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsil) + xksi(mi,mj) + 0.5*EPSILON(1.e0) )
          zratio   = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdfe) * z1_trbdia

          ! The minimum iron quota depends on the size of PSU, respiration
          ! and the reduction of nitrate following the parameterization 
          ! proposed by Flynn and Hipkin (1999)
          zironmin = xcoef1 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdch) * z1_trbdia + xcoef2 * zlim1 + xcoef3 * xdiatno3(mi,mj,jk)
          xqfuncfecd(mi,mj,jk) = zironmin + qdfelim
          zlim4    = MAX( 0., ( zratio - zironmin ) / qdfelim )
          xdiatpo4(mi,mj,jk) = zlim2
          xlimdfe (mi,mj,jk) = MIN( 1., zlim4 )
          xlimdia (mi,mj,jk) = MIN( zlim1, zlim2, zlim3, zlim4 )
          xlimsi  (mi,mj,jk) = MIN( zlim1, zlim2, zlim4 )
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO

      ! Size estimation of phytoplankton based on total biomass
      ! Assumes that larger biomass implies addition of larger cells
      ! ------------------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zcoef) &
      !$OMP SHARED(tr, jpphy, jpdia, Kbb, xsizern, xsizephy, &
                   xsizerd, xsizedia, sizena, sizeda)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zcoef = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) - MIN(xsizephy, t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) )
         sizena(mi,mj,jk) = 1. + ( xsizern -1.0 ) * zcoef / ( xsizephy + zcoef )
         zcoef = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) - MIN(xsizedia, t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) )
         sizeda(mi,mj,jk) = 1. + ( xsizerd - 1.0 ) * zcoef / ( xsizedia + zcoef )
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO

      ! Compute the fraction of nanophytoplankton that is made of calcifiers
      ! This is a purely adhoc formulation described in Aumont et al. (2015)
      ! This fraction depends on nutrient limitation, light, temperature
      ! --------------------------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zlim1, zlim2, zlim3, ztem1, ztem2, zetot1, zetot2) &
      !$OMP SHARED(tr, jppo4, jpfer, itemp, Kbb, Kmm, xnanonh4, xnanono3, &
                   concnnh4, etot_ndcy, caco3r, xsizephy, hmld, 0.5*EPSILON(1.e0), xfracal)
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zlim1  = xnanonh4(mi,mj,jk) + xnanono3(mi,mj,jk) 
         zlim2  = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) + concnnh4 )
         zlim3  = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer) +  6.E-11   )
         ztem1  = MAX( 0., t(mi,mj,N+1-jk,nnew,itemp) + 1.8)
         ztem2  = t(mi,mj,N+1-jk,nnew,itemp) - 10.
         zetot1 = MAX( 0., etot_ndcy(mi,mj,jk) - 1.) / ( 4. + etot_ndcy(mi,mj,jk) ) 
         zetot2 = 30. / ( 30.0 + etot_ndcy(mi,mj,jk) )

         xfracal(mi,mj,jk) = caco3r * MIN( zlim1, zlim2, zlim3 )                  &
            &                       * ztem1 / ( 0.1 + ztem1 )                     &
            &                       * MAX( 1., t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) / xsizephy )  &
            &                       * zetot1 * zetot2               &
            &                       * ( 1. + EXP(-ztem2 * ztem2 / 25. ) )         &
            &                       * MIN( 1., 50. / ( hbl(mi,mj) + 0.5*EPSILON(1.e0) ) )
         xfracal(mi,mj,jk) = MIN( 0.8 , xfracal(mi,mj,jk) )
         xfracal(mi,mj,jk) = MAX( 0.02, xfracal(mi,mj,jk) )
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      IF( .false. .AND. knt == nrdttrc ) THEN        ! save output diagnostics
        !
        IF( l_dia_fracal ) THEN   ! fraction of calcifiers
          ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
          !$OMP PARALLEL DO &
          !$OMP PRIVATE(mi, mj, jk) &
          !$OMP SHARED(zw3d, xfracal, tmask)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = xfracal(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          !$OMP END PARALLEL DO
          CALL iom_put( "xfracal",  zw3d)
          DEALLOCATE( zw3d )
        ENDIF
        !
        IF( l_dia_nut_lim ) THEN   ! Nutrient limitation term
          ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
          !$OMP PARALLEL DO &
          !$OMP PRIVATE(mi, mj, jk) &
          !$OMP SHARED(zw3d, xlimphy, tmask)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = xlimphy(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          !$OMP END PARALLEL DO
          CALL iom_put( "LNnut",  zw3d)
          !$OMP PARALLEL DO &
          !$OMP PRIVATE(mi, mj, jk) &
          !$OMP SHARED(zw3d, xlimdia, tmask)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = xlimdia(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          !$OMP END PARALLEL DO
          CALL iom_put( "LDnut",  zw3d)
          DEALLOCATE( zw3d )
        ENDIF
        !
        IF( l_dia_iron_lim ) THEN   ! Iron limitation term
          ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
          !$OMP PARALLEL DO &
          !$OMP PRIVATE(mi, mj, jk) &
          !$OMP SHARED(zw3d, xlimnfe, tmask)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = xlimnfe(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          !$OMP END PARALLEL DO
          CALL iom_put( "LNFe",  zw3d)
          !$OMP PARALLEL DO &
          !$OMP PRIVATE(mi, mj, jk) &
          !$OMP SHARED(zw3d, xlimdfe, tmask)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = xlimdfe(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          !$OMP END PARALLEL DO
          CALL iom_put( "LDFe",  zw3d)
          DEALLOCATE( zw3d )
        ENDIF
        !
        IF( l_dia_size_lim ) THEN   ! Size limitation term
          ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
          !$OMP PARALLEL DO &
          !$OMP PRIVATE(mi, mj, jk) &
          !$OMP SHARED(zw3d, sizen, tmask)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = sizen(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          !$OMP END PARALLEL DO
          CALL iom_put( "SIZEN",  zw3d)
          !$OMP PARALLEL DO &
          !$OMP PRIVATE(mi, mj, jk) &
          !$OMP SHARED(zw3d, sized, tmask)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = sized(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          !$OMP END PARALLEL DO
          CALL iom_put( "SIZED",  zw3d)
          DEALLOCATE( zw3d )
        ENDIF
        !
      ENDIF
      !
      IF( .false. )   CALL timing_stop('p4z_lim')
      !
   END SUBROUTINE p4z_lim


   SUBROUTINE p4z_lim_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_lim_init  ***
      !!
      !! ** Purpose :   Initialization of the nutrient limitation parameters
      !!
      !! ** Method  :   Read the namp4zlim namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namp4zlim
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer

      ! Namelist block
      NAMELIST/namp4zlim/ concnno3, concdno3, concnnh4, concdnh4, concnfer, concdfer, concbfe, &
         &                concbno3, concbnh4, xsizedia, xsizephy, xsizern, xsizerd, &
         &                xksi1, xksi2, xkdoc, qnfelim, qdfelim, caco3r, oxymin, ratchl
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*)
         WRITE(stdout,*) 'p4z_lim_init : initialization of nutrient limitations'
         WRITE(stdout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND(numnatp_ref);READ(numnatp_ref,namp4zlim,IOSTAT=ios);CALL ctl_nam(ios,"namp4zlim (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,namp4zlim,IOSTAT=ios);CALL ctl_nam(ios,"namp4zlim (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE( numonp, namp4zlim )

      !
      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*) '      Namelist : namp4zlim'
         WRITE(stdout,*) '      mean rainratio                           caco3r    = ', caco3r
         WRITE(stdout,*) '      C associated with Chlorophyll            ratchl    = ', ratchl
         WRITE(stdout,*) '      NO3 half saturation of nanophyto         concnno3  = ', concnno3
         WRITE(stdout,*) '      NO3 half saturation of diatoms           concdno3  = ', concdno3
         WRITE(stdout,*) '      NH4 half saturation for phyto            concnnh4  = ', concnnh4
         WRITE(stdout,*) '      NH4 half saturation for diatoms          concdnh4  = ', concdnh4
         WRITE(stdout,*) '      half saturation constant for Si uptake   xksi1     = ', xksi1
         WRITE(stdout,*) '      half saturation constant for Si/C        xksi2     = ', xksi2
         WRITE(stdout,*) '      half-sat. of DOC remineralization        xkdoc     = ', xkdoc
         WRITE(stdout,*) '      Iron half saturation for nanophyto       concnfer  = ', concnfer
         WRITE(stdout,*) '      Iron half saturation for diatoms         concdfer  = ', concdfer
         WRITE(stdout,*) '      size ratio for nanophytoplankton         xsizern   = ', xsizern
         WRITE(stdout,*) '      size ratio for diatoms                   xsizerd   = ', xsizerd
         WRITE(stdout,*) '      NO3 half saturation of bacteria          concbno3  = ', concbno3
         WRITE(stdout,*) '      NH4 half saturation for bacteria         concbnh4  = ', concbnh4
         WRITE(stdout,*) '      Minimum size criteria for diatoms        xsizedia  = ', xsizedia
         WRITE(stdout,*) '      Minimum size criteria for nanophyto      xsizephy  = ', xsizephy
         WRITE(stdout,*) '      Fe half saturation for bacteria          concbfe   = ', concbfe
         WRITE(stdout,*) '      halk saturation constant for anoxia      oxymin    =' , oxymin
         WRITE(stdout,*) '      optimal Fe quota for nano.               qnfelim   = ', qnfelim
         WRITE(stdout,*) '      Optimal Fe quota for diatoms             qdfelim   = ', qdfelim
      ENDIF
      !
      xfracal (:,:,N) = 0._wp
      xlimphy (:,:,N) = 0._wp    ;      xlimdia (:,:,N) = 0._wp
      xlimnfe (:,:,N) = 0._wp    ;      xlimdfe (:,:,N) = 0._wp
      xnanono3(:,:,N) = 0._wp    ;      xdiatno3(:,:,N) = 0._wp
      xnanofer(:,:,N) = 0._wp    ;      xdiatfer(:,:,N) = 0._wp
      xnanonh4(:,:,N) = 0._wp    ;      xdiatnh4(:,:,N) = 0._wp
      xnanopo4(:,:,N) = 0._wp    ;      xdiatpo4(:,:,N) = 0._wp
      xdiatpo4(:,:,N) = 0._wp    ;      xdiatpo4(:,:,N) = 0._wp
      xlimdia (:,:,N) = 0._wp    ;      xlimdfe (:,:,N) = 0._wp
      concnfe (:,:,N) = 0._wp    ;      concdfe (:,:,N) = 0._wp
      xqfuncfecn(:,:,N) = 0._wp    ;    xqfuncfecd(:,:,N) = 0._wp
      xlimsi  (:,:,N) = 0._wp
      xlimbac (:,:,N) = 0._wp    ;      xlimbacl(:,:,N) = 0._wp
      !
   END SUBROUTINE p4z_lim_init


   INTEGER FUNCTION p4z_lim_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_lim_alloc  ***
      !! 
      !            Allocation of the arrays used in this module
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_stop
      !!----------------------------------------------------------------------

      !*  Biological arrays for phytoplankton growth
      ALLOCATE( xdiatno3(Istrp:Iendp,Jstrp:Jendp,N),                             &
         &      xnanonh4(Istrp:Iendp,Jstrp:Jendp,N), xdiatnh4(Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xnanopo4(Istrp:Iendp,Jstrp:Jendp,N), xdiatpo4(Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xnanofer(Istrp:Iendp,Jstrp:Jendp,N), xdiatfer(Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xlimdia (Istrp:Iendp,Jstrp:Jendp,N), xlimdfe (Istrp:Iendp,Jstrp:Jendp,N),       &
         &      concnfe (Istrp:Iendp,Jstrp:Jendp,N), concdfe (Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xqfuncfecn(Istrp:Iendp,Jstrp:Jendp,N), xqfuncfecd(Istrp:Iendp,Jstrp:Jendp,N),   &
         &      xlimsi (Istrp:Iendp,Jstrp:Jendp,N), STAT=p4z_lim_alloc )
      !
      IF( p4z_lim_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p4z_lim_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p4z_lim_alloc


   !!======================================================================
END MODULE p4zlim
