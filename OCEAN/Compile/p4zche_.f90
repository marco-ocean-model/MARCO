










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









MODULE p4zche
   !!======================================================================
   !!                         ***  MODULE p4zche  ***
   !! TOP :    Sea water chemistry computed following OCMIP protocol
   !!======================================================================
   !! History :   OPA  !  1988     (E. Maier-Reimer)  Original code
   !!              -   !  1998     (O. Aumont)  addition
   !!              -   !  1999     (C. Le Quere)  modification
   !!   NEMO      1.0  !  2004     (O. Aumont)  modification
   !!              -   !  2006     (R. Gangsto)  modification
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                  !  2011-02  (J. Simeon, J.Orr ) update O2 solubility constants
   !!             3.6  !  2016-03  (O. Aumont) Change chemistry to MOCSY standards
   !!             4.2  !  2020     (J. ORR )  rhop is replaced by "in situ  density" rhd
   !!             3.*  !  2025-02  (S. Maishal, R. Person)  and optimization
   !!----------------------------------------------------------------------
   !!   p4z_che      :  Sea water chemistry computed following OCMIP protocol
   !!----------------------------------------------------------------------
   USE oce_trc       !  shared variables between ocean and passive tracers
   USE trc           !  passive tracers common variables
   USE sms_pisces    !   Source Minus Sink variables
   USE lib_mpp       !  MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_che          !
   PUBLIC   p4z_che_alloc    !
   PUBLIC   ahini_for_at     !
   PUBLIC   solve_at_general !

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: sio3eq   ! chemistry of Si
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: fekeq    ! chemistry of Fe
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: chemc    ! Solubilities of O2 and CO2
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   :: chemo2    ! Solubilities of O2 and CO2
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) :: fesol    ! solubility of Fe
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   salinprac  ! Practical salinity
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tempis   ! In situ temperature

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akb3       !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akw3       !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   akf3       !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aks3       !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak1p3      !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak2p3      !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   ak3p3      !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   aksi3      !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   borat      !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   fluorid    !: ???
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sulfat     !: ???

   !!* Variable for chemistry of the CO2 cycle

   REAL(wp), PUBLIC ::   atcox  = 0.20946         ! units atm

   REAL(wp) ::   o2atm  = 1. / ( 1000. * 0.20946 )  

   REAL(wp) ::   rgas   = 83.14472      ! universal gas constants
   REAL(wp) ::   oxyco  = 1. / 22.4144  ! converts from liters of an ideal gas to moles
   !                                    ! coeff. for seawater pressure correction : millero 95
   !                                    ! AGRIF doesn't like the DATA instruction
   REAL(wp) :: devk10  = -25.5
   REAL(wp) :: devk11  = -15.82
   REAL(wp) :: devk12  = -29.48
   REAL(wp) :: devk13  = -20.02
   REAL(wp) :: devk14  = -18.03
   REAL(wp) :: devk15  = -9.78
   REAL(wp) :: devk16  = -48.76
   REAL(wp) :: devk17  = -14.51
   REAL(wp) :: devk18  = -23.12
   REAL(wp) :: devk19  = -26.57
   REAL(wp) :: devk110  = -29.48
   !
   REAL(wp) :: devk20  = 0.1271
   REAL(wp) :: devk21  = -0.0219
   REAL(wp) :: devk22  = 0.1622
   REAL(wp) :: devk23  = 0.1119
   REAL(wp) :: devk24  = 0.0466
   REAL(wp) :: devk25  = -0.0090
   REAL(wp) :: devk26  = 0.5304
   REAL(wp) :: devk27  = 0.1211
   REAL(wp) :: devk28  = 0.1758
   REAL(wp) :: devk29  = 0.2020
   REAL(wp) :: devk210  = 0.1622
   !
   REAL(wp) :: devk30  = 0.
   REAL(wp) :: devk31  = 0.
   REAL(wp) :: devk32  = 2.608E-3
   REAL(wp) :: devk33  = -1.409e-3
   REAL(wp) :: devk34  = 0.316e-3
   REAL(wp) :: devk35  = -0.942e-3
   REAL(wp) :: devk36  = 0.
   REAL(wp) :: devk37  = -0.321e-3
   REAL(wp) :: devk38  = -2.647e-3
   REAL(wp) :: devk39  = -3.042e-3
   REAL(wp) :: devk310  = -2.6080e-3
   !
   REAL(wp) :: devk40  = -3.08E-3
   REAL(wp) :: devk41  = 1.13E-3
   REAL(wp) :: devk42  = -2.84E-3
   REAL(wp) :: devk43  = -5.13E-3
   REAL(wp) :: devk44  = -4.53e-3
   REAL(wp) :: devk45  = -3.91e-3
   REAL(wp) :: devk46  = -11.76e-3
   REAL(wp) :: devk47  = -2.67e-3
   REAL(wp) :: devk48  = -5.15e-3
   REAL(wp) :: devk49  = -4.08e-3
   REAL(wp) :: devk410  = -2.84e-3
   !
   REAL(wp) :: devk50  = 0.0877E-3
   REAL(wp) :: devk51  = -0.1475E-3     
   REAL(wp) :: devk52  = 0.
   REAL(wp) :: devk53  = 0.0794E-3      
   REAL(wp) :: devk54  = 0.09e-3
   REAL(wp) :: devk55  = 0.054e-3
   REAL(wp) :: devk56  = 0.3692E-3
   REAL(wp) :: devk57  = 0.0427e-3
   REAL(wp) :: devk58  = 0.09e-3
   REAL(wp) :: devk59  = 0.0714e-3
   REAL(wp) :: devk510  = 0.0
   !
   ! General parameters
   REAL(wp), PARAMETER :: pp_rdel_ah_target = 1.E-4_wp
   REAL(wp), PARAMETER :: pp_ln10 = 2.302585092994045684018_wp

   ! Maximum number of iterations for each method
   INTEGER, PARAMETER :: jp_maxniter_atgen    = 20

   ! Bookkeeping variables for each method
   ! - SOLVE_AT_GENERAL
   INTEGER :: niter_atgen    = jp_maxniter_atgen

   !! * Substitutions






























   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zche.F90 15532 2021-11-24 11:47:32Z techene $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_che( Kbb, Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_che  ***
      !!
      !! ** Purpose :   Sea water chemistry computed following OCMIP protocol
      !!
      !! ** Method  : - ...
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   Kbb, Kmm  ! time level indices
      INTEGER  ::   mi, mj, jk
      REAL(wp) ::   ztkel, ztkel1, zt , zsal  , zsal2 , zbuf1 , zbuf2
      REAL(wp) ::   ztgg , ztgg2, ztgg3 , ztgg4 , ztgg5
      REAL(wp) ::   zpres, ztc  , zcl   , zcpexp, zoxy  , zcpexp2
      REAL(wp) ::   zsqrt, ztr  , zlogt , zcek1, zc1, zplat
      REAL(wp) ::   zis  , zis2 , zsal15, zisqrt, za1, za2
      REAL(wp) ::   zckb , zck1 , zck2  , zckw  , zak1 , zak2  , zakb , zaksp0, zakw
      REAL(wp) ::   zck1p, zck2p, zck3p, zcksi, zak1p, zak2p, zak3p, zaksi
      REAL(wp) ::   zst  , zft  , zcks  , zckf  , zaksp1
      REAL(wp) ::   total2free, free2SWS, total2SWS, SWS2total
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p4z_che')
      !
      ! Computation of chemical constants require practical salinity
      ! Thus, when TEOS08 is used, absolute salinity is converted to 
      ! practical salinity
      ! -------------------------------------------------------------
      IF (0 == -1) THEN
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(jk, mi, mj) &
      !$OMP SHARED(ts, salinprac, itemp+ntrc_salt, Kmm, N)
         DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            salinprac(mi,mj,jk) = t(mi,mj,N+1-jk,nnew,itemp+ntrc_salt) * 35.0 / 35.16504
         END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      ELSE
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(jk, mi, mj) &
      !$OMP SHARED(ts, salinprac, itemp+ntrc_salt, Kmm, N)
         DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            salinprac(mi,mj,jk) = t(mi,mj,N+1-jk,nnew,itemp+ntrc_salt)
         END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      ENDIF
      !
      ! Computations of chemical constants require in situ temperature
      ! Here a quite simple formulation is used to convert 
      ! potential temperature to in situ temperature. The errors is less than 
      ! 0.04°C relative to an exact computation
      ! ---------------------------------------------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(jk, mi, mj, zpres, za1, za2) &
      !$OMP SHARED(ts, salinprac, tempis, itemp, Kmm, N)
      DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zpres = ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) / 1000.
         za1 = 0.04 * ( 1.0 + 0.185 * t(mi,mj,N+1-jk,nnew,itemp) + 0.035 &
              * (salinprac(mi,mj,jk) - 35.0) )
         za2 = 0.0075 * ( 1.0 - t(mi,mj,N+1-jk,nnew,itemp) / 30.0 )
         tempis(mi,mj,jk) = t(mi,mj,N+1-jk,nnew,itemp) - za1 &
                            * zpres + za2 * zpres**2
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      ! CHEMICAL CONSTANTS - SURFACE LAYER
      ! ----------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, ztkel, zt, zsal, zcek1) &
      !$OMP SHARED(tempis, salinprac, tmask, chemc)
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         !                             ! SET ABSOLUTE 
         ztkel = tempis(mi,mj,1) + 273.15
         zt    = ztkel * 0.01
         zsal  = salinprac(mi,mj,1) + ( 1.- tmask(mi,mj,1) ) * 35.
         !                             ! LN(K0) OF SOLUBILITY OF CO2 (EQ. 12, WEISS, 1980)
         !                             !     AND FOR THE ATMOSPHERE FOR NON IDEAL GAS
         ! J. ORR: The previous code has been modified. It computed CO2 solubility in mol/(kg*atm), then converted that to mol/(L*atm).
         ! But Weiss (1974) provides sets of coefficients for each of those 2 units.
         ! Thus I have changed the code to use the coefficients for mol*L/atm.
         ! Hence I've eliminated using the conversion (which used the variable rhop)
         ! OLD - Coefficients for CO2 soulbility in mol/(kg*atm) (Weiss,1974, Table 1, column 2)
         !zcek1 = 9345.17/ztkel - 60.2409 + 23.3585 * LOG(zt) + zsal*(0.023517 - 0.00023656*ztkel    &
         !&       + 0.0047036e-4*ztkel**2)
         ! NEW - Coefficients for CO2 soulbility in mol/(L*atm) (Weiss, 1974, Table 1, column 1)
         zcek1 = 9050.69/ztkel - 58.0931 + 22.2940 * LOG(zt) + zsal*(0.027766 - 0.00025888*ztkel    &
                 &       + 0.0050578e-4*ztkel**2)
         !
         ! OLD:  chemc(mi,mj,1) = EXP( zcek1 ) * 1E-6 * (rho0+rho1(mi,mj,N+1-1)) / 1000. ! mol/(L atm)
         ! The units indicated in the above line are wrong. They are actually "mol/(L*uatm)"
         ! NEW:
         chemc(mi,mj,1) = EXP( zcek1 ) * 1E-6 ! mol/(L * uatm)
         chemc(mi,mj,2) = -1636.75 + 12.0408*ztkel - 0.0327957*ztkel**2 + 0.0000316528*ztkel**3
         chemc(mi,mj,3) = 57.7 - 0.118*ztkel
      END DO   ;   END DO
      !$OMP END PARALLEL DO
      ! OXYGEN SOLUBILITY - DEEP OCEAN
      ! -------------------------------
      !$OMP PARALLEL DO PRIVATE(mi, mj, jk, ztkel, zsal, &
                              zsal2, ztgg, ztgg2, ztgg3, &
                                     ztgg4, ztgg5, zoxy) &
      !$OMP SHARED(tempis, salinprac, tmask, chemo2, &
                           o2atm, oxyco, atcox, N)
      DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ztkel = tempis(mi,mj,jk) + 273.15
         zsal  = salinprac(mi,mj,jk) + ( 1.- tmask(mi,mj,jk) ) * 35.
         zsal2 = zsal * zsal
         ztgg  = LOG( ( 298.15 - tempis(mi,mj,jk) ) / ztkel )  ! Set the GORDON & GARCIA scaled temperature
         ztgg2 = ztgg  * ztgg
         ztgg3 = ztgg2 * ztgg
         ztgg4 = ztgg3 * ztgg
         ztgg5 = ztgg4 * ztgg
         zoxy  = 2.00856 + 3.22400 * ztgg + 3.99063 * ztgg2 + 4.80299 * ztgg3    &
         &       + 9.78188e-1 * ztgg4 + 1.71069 * ztgg5 + zsal * ( -6.24097e-3   &
         &       - 6.93498e-3 * ztgg - 6.90358e-3 * ztgg2 - 4.29155e-3 * ztgg3 )   &
         &       - 3.11680e-7 * zsal2
         chemo2(mi,mj,jk) = ( EXP( zoxy ) * o2atm ) * oxyco * atcox     ! mol/(L atm)
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      ! CHEMICAL CONSTANTS - DEEP OCEAN
      ! -------------------------------
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, ztkel, zsal, zsal2, ztgg, ztgg2, &
                               ztgg3, ztgg4, ztgg5, zoxy, zplat, zc1, zpres, &
                               zsqrt, zsal15, zlogt, ztr, zis, zis2, zisqrt, &
                                 ztc, zcl, zst, zft, zcks, zckf, zckb, zck1, &
                             zck2, zckw, zck1p, zck2p, zck3p, zcksi, zaksp0, &
                           total2free, free2SWS, total2SWS, SWS2total, zak1, &
                              zak2, zakb, zakw, zaksp1, zak1p, zak2p, zak3p, &
                                       zaksi, zcpexp, zcpexp2, zbuf1, zbuf2) &
      !$OMP SHARED(tempis, salinprac, tmask, chemo2, o2atm, oxyco, atcox, &
                           gdept, devk10, devk20, devk30, devk40, devk50, &
                  devk11, devk21, devk31, devk41, devk51, devk12, devk22, &
                  devk32, devk42, devk52, devk13, devk23, devk33, devk43, &
                  devk53, devk14, devk24, devk34, devk44, devk54, devk15, &
                  devk25, devk35, devk45, devk55, devk17, devk27, devk37, &
                  devk47, devk57, devk18, devk28, devk38, devk48, devk58, &
                        devk19, devk29, devk39, devk49, devk59)
      DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
          ! SET PRESSION ACCORDING TO SAUNDER (1980)
          zplat   = SIN ( ABS(latr(mi,mj)*3.141592654/180.) )
          zc1 = 5.92E-3 + zplat**2 * 5.25E-3
          zpres = ((1-zc1)-SQRT(((1-zc1)**2)-(8.84E-6*((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N)))))) / 4.42E-6
          zpres = zpres / 10.0

          ! SET ABSOLUTE 
          ztkel   = tempis(mi,mj,jk) + 273.15
          zsal    = salinprac(mi,mj,jk) + ( 1.-tmask(mi,mj,jk) ) * 35.
          zsqrt  = SQRT( zsal )
          zsal15  = zsqrt * zsal
          zlogt  = LOG( ztkel )
          ztr    = 1. / ztkel
          zis    = 19.924 * zsal / ( 1000.- 1.005 * zsal )
          zis2   = zis * zis
          zisqrt = SQRT( zis )
          ztc     = tempis(mi,mj,jk) + ( 1.- tmask(mi,mj,jk) ) * 20.

          ! CHLORINITY (WOOSTER ET AL., 1969)
          zcl     = zsal / 1.80655

          ! TOTAL SULFATE CONCENTR. [MOLES/kg soln]
          zst     = 0.14 * zcl /96.062

          ! TOTAL FLUORIDE CONCENTR. [MOLES/kg soln]
          zft     = 0.000067 * zcl /18.9984

          ! DISSOCIATION CONSTANT FOR SULFATES on free H scale (Dickson 1990)
          zcks    = EXP(-4276.1 * ztr + 141.328 - 23.093 * zlogt         &
          &         + (-13856. * ztr + 324.57 - 47.986 * zlogt) * zisqrt &
          &         + (35474. * ztr - 771.54 + 114.723 * zlogt) * zis    &
          &         - 2698. * ztr * zis**1.5 + 1776.* ztr * zis2         &
          &         + LOG(1.0 - 0.001005 * zsal))

          ! DISSOCIATION CONSTANT FOR FLUORIDES on free H scale (Dickson and Riley 79)
          zckf    = EXP( 1590.2*ztr - 12.641 + 1.525*zisqrt   &
          &         + LOG(1.0d0 - 0.001005d0*zsal)            &
          &         + LOG(1.0d0 + zst/zcks))

          ! DISSOCIATION CONSTANT FOR CARBONATE AND BORATE
          zckb=  (-8966.90 - 2890.53*zsqrt - 77.942*zsal        &
          &      + 1.728*zsal15 - 0.0996*zsal*zsal)*ztr         &
          &      + (148.0248 + 137.1942*zsqrt + 1.62142*zsal)   &
          &      + (-24.4344 - 25.085*zsqrt - 0.2474*zsal)      & 
          &      * zlogt + 0.053105*zsqrt*ztkel

          ! DISSOCIATION COEFFICIENT FOR CARBONATE ACCORDING TO 
          ! MEHRBACH (1973) REFIT BY MILLERO (1995), seawater scale
          zck1    = -1.0*(3633.86*ztr - 61.2172 + 9.6777*zlogt  &
             - 0.011555*zsal + 0.0001152*zsal*zsal)
          zck2    = -1.0*(471.78*ztr + 25.9290 - 3.16967*zlogt      &
             - 0.01781*zsal + 0.0001122*zsal*zsal)

          ! PKW (H2O) (MILLERO, 1995) from composite data
          zckw    = -13847.26 * ztr + 148.9652 - 23.6521 * zlogt + ( 118.67 * ztr    &
                    - 5.977 + 1.0495 * zlogt ) * zsqrt - 0.01615 * zsal

          ! CONSTANTS FOR PHOSPHATE (MILLERO, 1995)
         zck1p    = -4576.752*ztr + 115.540 - 18.453*zlogt   &
         &          + (-106.736*ztr + 0.69171) * zsqrt       &
         &          + (-0.65643*ztr - 0.01844) * zsal

         zck2p    = -8814.715*ztr + 172.1033 - 27.927*zlogt  &
         &          + (-160.340*ztr + 1.3566)*zsqrt          &
         &          + (0.37335*ztr - 0.05778)*zsal

         zck3p    = -3070.75*ztr - 18.126                    &
         &          + (17.27039*ztr + 2.81197) * zsqrt       &
         &          + (-44.99486*ztr - 0.09984) * zsal 

         ! CONSTANT FOR SILICATE, MILLERO (1995)
         zcksi    = -8904.2*ztr  + 117.400 - 19.334*zlogt   &
         &          + (-458.79*ztr + 3.5913) * zisqrt       &
         &          + (188.74*ztr - 1.5998) * zis           &
         &          + (-12.1652*ztr + 0.07871) * zis2       &
         &          + LOG(1.0 - 0.001005*zsal)

          ! APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE IN SEAWATER
          !       (S=27-43, T=2-25 DEG C) at pres =0 (atmos. pressure) (MUCCI 1983)
          zaksp0  = -171.9065 -0.077993*ztkel + 2839.319*ztr + 71.595*LOG10( ztkel )   &
             &      + (-0.77712 + 0.00284263*ztkel + 178.34*ztr) * zsqrt  &
             &      - 0.07711*zsal + 0.0041249*zsal15

          ! CONVERT FROM DIFFERENT PH SCALES
          total2free  = 1.0/(1.0 + zst/zcks)
          free2SWS    = 1. + zst/zcks + zft/(zckf*total2free)
          total2SWS   = total2free * free2SWS
          SWS2total   = 1.0 / total2SWS

          ! K1, K2 OF CARBONIC ACID, KB OF BORIC ACID, KW (H2O) (LIT.?)
          zak1    = 10**(zck1) * total2SWS
          zak2    = 10**(zck2) * total2SWS
          zakb    = EXP( zckb ) * total2SWS
          zakw    = EXP( zckw )
          zaksp1  = 10**(zaksp0)
          zak1p   = exp( zck1p )
          zak2p   = exp( zck2p )
          zak3p   = exp( zck3p )
          zaksi   = exp( zcksi )
          zckf    = zckf * total2SWS

          ! FORMULA FOR CPEXP AFTER EDMOND & GIESKES (1970)
          !        (REFERENCE TO CULBERSON & PYTKOQICZ (1968) AS MADE
          !        IN BROECKER ET AL. (1982) IS INCORRECT; HERE RGAS IS
          !        TAKEN TENFOLD TO CORRECT FOR THE NOTATION OF pres  IN
          !        DBAR INSTEAD OF BAR AND THE EXPRESSION FOR CPEXP IS
          !        MULTIPLIED BY LN(10.) TO ALLOW USE OF EXP-FUNCTION
          !        WITH BASIS E IN THE FORMULA FOR AKSPP (CF. EDMOND
          !        & GIESKES (1970), P. 1285-1286 (THE SMALL
          !        FORMULA ON P. 1286 IS RIGHT AND CONSISTENT WITH THE
          !        SIGN IN PARTIAL MOLAR VOLUME CHANGE AS SHOWN ON P. 1285))
          zcpexp  = zpres / (rgas*ztkel)
          zcpexp2 = zpres * zcpexp

          ! KB OF BORIC ACID, K1,K2 OF CARBONIC ACID PRESSURE
          !        CORRECTION AFTER CULBERSON AND PYTKOWICZ (1968)
          !        (CF. BROECKER ET AL., 1982)

          zbuf1  = -     ( devk10 + devk20 * ztc + devk30 * ztc * ztc )
          zbuf2  = 0.5 * ( devk40 + devk50 * ztc )
          ak13(mi,mj,jk) = zak1 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

          zbuf1  =     - ( devk11 + devk21 * ztc + devk31 * ztc * ztc )
          zbuf2  = 0.5 * ( devk41 + devk51 * ztc )
          ak23(mi,mj,jk) = zak2 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

          zbuf1  =     - ( devk12 + devk22 * ztc + devk32 * ztc * ztc )
          zbuf2  = 0.5 * ( devk42 + devk52 * ztc )
          akb3(mi,mj,jk) = zakb * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

          zbuf1  =     - ( devk13 + devk23 * ztc + devk33 * ztc * ztc )
          zbuf2  = 0.5 * ( devk43 + devk53 * ztc )
          akw3(mi,mj,jk) = zakw * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

          zbuf1  =     - ( devk14 + devk24 * ztc + devk34 * ztc * ztc )
          zbuf2  = 0.5 * ( devk44 + devk54 * ztc )
          aks3(mi,mj,jk) = zcks * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

          zbuf1  =     - ( devk15 + devk25 * ztc + devk35 * ztc * ztc )
          zbuf2  = 0.5 * ( devk45 + devk55 * ztc )
          akf3(mi,mj,jk) = zckf * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

          zbuf1  =     - ( devk17 + devk27 * ztc + devk37 * ztc * ztc )
          zbuf2  = 0.5 * ( devk47 + devk57 * ztc )
          ak1p3(mi,mj,jk) = zak1p * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

          zbuf1  =     - ( devk18 + devk28 * ztc + devk38 * ztc * ztc )
          zbuf2  = 0.5 * ( devk48 + devk58 * ztc )
          ak2p3(mi,mj,jk) = zak2p * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

          zbuf1  =     - ( devk19 + devk29 * ztc + devk39 * ztc * ztc )
          zbuf2  = 0.5 * ( devk49 + devk59 * ztc )
          ak3p3(mi,mj,jk) = zak3p * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

          zbuf1  =     - ( devk110 + devk210 * ztc + devk310 * ztc * ztc )
          zbuf2  = 0.5 * ( devk410 + devk510 * ztc )
          aksi3(mi,mj,jk) = zaksi * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

          ! CONVERT FROM DIFFERENT PH SCALES
          total2free  = 1.0/(1.0 + zst/aks3(mi,mj,jk))
          free2SWS    = 1. + zst/aks3(mi,mj,jk) + zft/akf3(mi,mj,jk)
          total2SWS   = total2free * free2SWS
          SWS2total   = 1.0 / total2SWS

          ! Convert to total scale
          ak13(mi,mj,jk)  = ak13(mi,mj,jk)  * SWS2total
          ak23(mi,mj,jk)  = ak23(mi,mj,jk)  * SWS2total
          akb3(mi,mj,jk)  = akb3(mi,mj,jk)  * SWS2total
          akw3(mi,mj,jk)  = akw3(mi,mj,jk)  * SWS2total
          ak1p3(mi,mj,jk) = ak1p3(mi,mj,jk) * SWS2total
          ak2p3(mi,mj,jk) = ak2p3(mi,mj,jk) * SWS2total
          ak3p3(mi,mj,jk) = ak3p3(mi,mj,jk) * SWS2total
          aksi3(mi,mj,jk) = aksi3(mi,mj,jk) * SWS2total
          akf3(mi,mj,jk)  = akf3(mi,mj,jk)  / total2free

          ! APPARENT SOLUBILITY PRODUCT K'SP OF CALCITE 
          !        AS FUNCTION OF PRESSURE FOLLOWING MILLERO
          !        (P. 1285) AND BERNER (1976)
          zbuf1  =     - ( devk16 + devk26 * ztc + devk36 * ztc * ztc )
          zbuf2  = 0.5 * ( devk46 + devk56 * ztc )
          aksp(mi,mj,jk) = zaksp1 * EXP( zbuf1 * zcpexp + zbuf2 * zcpexp2 )

          ! TOTAL F, S, and BORATE CONCENTR. [MOLES/L]
          borat(mi,mj,jk) = 0.0002414 * zcl / 10.811
          sulfat(mi,mj,jk) = zst
          fluorid(mi,mj,jk) = zft 

          fekeq (mi,mj,jk) = 10**( 17.27 - 1565.7 / ztkel ) 

          ! Liu and Millero (1999) only valid 5 - 50 degC
          ztkel1 = MAX( 5. , tempis(mi,mj,jk) ) + 273.16
          fesol(mi,mj,jk,1) = 10**(-13.486 - 0.1856* zis**0.5 + 0.3073*zis + 5254.0/ztkel1)
          fesol(mi,mj,jk,2) = 10**(2.517 - 0.8885*zis**0.5 + 0.2139 * zis - 1320.0/ztkel1 )
          fesol(mi,mj,jk,3) = 10**(0.4511 - 0.3305*zis**0.5 - 1996.0/ztkel1 )
          fesol(mi,mj,jk,4) = 10**(-0.2965 - 0.7881*zis**0.5 - 4086.0/ztkel1 )
          fesol(mi,mj,jk,5) = 10**(4.4466 - 0.8505*zis**0.5 - 7980.0/ztkel1 )
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO

      ! Iron and SIO3 saturation concentration from ...
      IF( .NOT. ln_p2z) THEN
      !$OMP PARALLEL DO & 
      !$OMP PRIVATE(mi, mj, jk, ztkel) & 
      !$OMP SHARED(tempis, sio3eq, N)
         DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             ! SET ABSOLUTE 
             ztkel   = tempis(mi,mj,jk) + 273.15
             sio3eq(mi,mj,jk) = EXP(  LOG( 10.) * ( 6.44 - 968. / ztkel )  ) * 1.e-6
            !
         END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      ENDIF
      !
      IF( .false. )  CALL timing_stop('p4z_che')
      !
   END SUBROUTINE p4z_che

   SUBROUTINE ahini_for_at(p_hini, Kbb )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE ahini_for_at  ***
      !!
      !! Subroutine returns the root for the 2nd order approximation of the
      !! DIC -- B_T -- A_CB equation for [H+] (reformulated as a cubic 
      !! polynomial) around the local minimum, if it exists.
      !! Returns * 1E-03_wp if p_alkcb <= 0
      !!         * 1E-10_wp if p_alkcb >= 2*p_dictot + p_bortot
      !!         * 1E-07_wp if 0 < p_alkcb < 2*p_dictot + p_bortot
      !!                    and the 2nd order approximation does not have 
      !!                    a solution
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N), INTENT(OUT)  ::  p_hini
      INTEGER,                          INTENT(in)   ::  Kbb      ! time level indices
      INTEGER  ::   mi, mj, jk
      REAL(wp)  ::  zca1, zba1
      REAL(wp)  ::  zd, zsqrtd, zhmin
      REAL(wp)  ::  za2, za1, za0, zrhd
      REAL(wp)  ::  p_dictot, p_bortot, p_alkcb 
      !!---------------------------------------------------------------------

      IF( .false. )  CALL timing_start('ahini_for_at')
      !
      !$OMP PARALLEL DO & 
      !$OMP PRIVATE(mi, mj, jk, zrhd, p_alkcb, p_dictot, &
                    p_bortot, zca1, zba1, za2, za1, za0, &
                    zd, zsqrtd, zhmin) & 
      !$OMP SHARED(rhd, tr, borat, p_hini, jptal, Kbb, &
                            jpdic, 0.5*EPSILON(1.e0), ak13, akb3, ak23)
      DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
      zrhd = 1._wp / ( rho(mi,mj,N+1-jk) + 1. )
      p_alkcb  = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jptal) * zrhd
      p_dictot = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdic) * zrhd
      p_bortot = borat(mi,mj,jk)
      IF (p_alkcb <= 0.) THEN
          p_hini(mi,mj,jk) = 1.e-3
      ELSEIF (p_alkcb >= (2.*p_dictot + p_bortot)) THEN
          p_hini(mi,mj,jk) = 1.e-10_wp
      ELSE
          zca1 = p_dictot/( p_alkcb + 0.5*EPSILON(1.e0) )
          zba1 = p_bortot/ (p_alkcb + 0.5*EPSILON(1.e0) )
     ! Coefficients of the cubic polynomial
          za2 = aKb3(mi,mj,jk)*(1. - zba1) + ak13(mi,mj,jk)*(1.-zca1)
          za1 = ak13(mi,mj,jk)*akb3(mi,mj,jk)*(1. - zba1 - zca1)    &
          &     + ak13(mi,mj,jk)*ak23(mi,mj,jk)*(1. - (zca1+zca1))
          za0 = ak13(mi,mj,jk)*ak23(mi,mj,jk)*akb3(mi,mj,jk)*(1. - zba1 - (zca1+zca1))
                                  ! Taylor expansion around the minimum
          zd = za2*za2 - 3.*za1   ! Discriminant of the quadratic equation
                                  ! for the minimum close to the root

          IF(zd > 0.) THEN        ! If the discriminant is positive
            zsqrtd = SQRT(zd)
            IF(za2 < 0) THEN
              zhmin = (-za2 + zsqrtd)/3.
            ELSE
              zhmin = -za1/(za2 + zsqrtd)
            ENDIF
            p_hini(mi,mj,jk) = zhmin + SQRT(-(za0 + zhmin*(za1 + zhmin*(za2 + zhmin)))/zsqrtd)
          ELSE
            p_hini(mi,mj,jk) = 1.e-7
          ENDIF
       !
       ENDIF
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      IF( .false. )  CALL timing_stop('ahini_for_at')
      !
   END SUBROUTINE ahini_for_at

   !===============================================================================

   SUBROUTINE anw_infsup( p_alknw_inf, p_alknw_sup, Kbb )

   ! Subroutine returns the lower and upper bounds of "non-water-selfionization"
   ! contributions to total alkalinity (the infimum and the supremum), i.e
   ! inf(TA - [OH-] + [H+]) and sup(TA - [OH-] + [H+])

   ! Argument variables
   REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N), INTENT(OUT) :: p_alknw_inf
   REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N), INTENT(OUT) :: p_alknw_sup
   INTEGER,                          INTENT(in)  ::  Kbb      ! time level indices
   INTEGER  ::   mi, mj, jk
   REAL(wp)  ::  zrhd

   IF( ln_p2z ) THEN
      !$OMP PARALLEL DO & 
      !$OMP PRIVATE(mi, mj, jk, zrhd) & 
      !$OMP SHARED(rhd, p_alknw_inf, p_alknw_sup, sulfat, &
                            fluorid, tr, jpdic, Kbb, borat)
      DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zrhd = 1._wp / ( rho(mi,mj,N+1-jk) + 1. )
         p_alknw_inf(mi,mj,jk) =  -2.174E-6 * zrhd - sulfat(mi,mj,jk) &
         &              - fluorid(mi,mj,jk)
         p_alknw_sup(mi,mj,jk) =   (2. * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdic) + 2. * 2.174E-6    &
         &               + 90.33E-6 ) * zrhd + borat(mi,mj,jk)
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
   ELSE
      !$OMP PARALLEL DO & 
      !$OMP PRIVATE(mi, mj, jk, zrhd) & 
      !$OMP SHARED(rhd, p_alknw_inf, p_alknw_sup, tr, jppo4, &
              jpdic, jpsil, Kbb, po4r, sulfat, fluorid, borat)
      DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zrhd = 1._wp / ( rho(mi,mj,N+1-jk) + 1. )
         p_alknw_inf(mi,mj,jk) =  -t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) * zrhd &
                                &      * po4r - sulfat(mi,mj,jk) &
                                &            - fluorid(mi,mj,jk)
         p_alknw_sup(mi,mj,jk) =   (2. * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdic) &
         &                + 2. * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) * po4r &
         &               + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsil) ) * zrhd + borat(mi,mj,jk)
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
   ENDIF

   END SUBROUTINE anw_infsup


   SUBROUTINE solve_at_general( p_hini, zhi, Kbb )

   ! Universal pH solver that converges from any given initial value,
   ! determines upper an lower bounds for the solution if required

   ! Argument variables
   !--------------------
   REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N), INTENT(IN)   :: p_hini
   REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N), INTENT(OUT)  :: zhi
   INTEGER,                          INTENT(in)   :: Kbb  ! time level indices

   ! Local variables
   !-----------------
   INTEGER   ::  mi, mj, jk, jn
   REAL(wp)  ::  zh_ini, zh, zh_prev, zh_lnfactor
   REAL(wp)  ::  zdelta, zh_delta
   REAL(wp)  ::  zeqn, zdeqndh, zalka
   REAL(wp)  ::  aphscale
   REAL(wp)  ::  znumer_dic, zdnumer_dic, zdenom_dic, zalk_dic, zdalk_dic
   REAL(wp)  ::  znumer_bor, zdnumer_bor, zdenom_bor, zalk_bor, zdalk_bor
   REAL(wp)  ::  znumer_po4, zdnumer_po4, zdenom_po4, zalk_po4, zdalk_po4
   REAL(wp)  ::  znumer_sil, zdnumer_sil, zdenom_sil, zalk_sil, zdalk_sil
   REAL(wp)  ::  znumer_so4, zdnumer_so4, zdenom_so4, zalk_so4, zdalk_so4
   REAL(wp)  ::  znumer_flu, zdnumer_flu, zdenom_flu, zalk_flu, zdalk_flu
   REAL(wp)  ::  zalk_wat, zdalk_wat
   REAL(wp)  ::  zrhd, p_alktot, zdic, zbot, zpt, zst, zft, zsit
   LOGICAL   ::  l_exitnow
   REAL(wp), PARAMETER :: pz_exp_threshold = 1.0
   REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N) :: zalknw_inf, zalknw_sup, rmask, zh_min, zh_max, zeqn_absmin

   IF( .false. )  CALL timing_start('solve_at_general')

   CALL anw_infsup( zalknw_inf, zalknw_sup, Kbb )

   rmask(Istrp:Iendp,Jstrp:Jendp,1:N) = tmask(Istrp:Iendp,Jstrp:Jendp,1:N)
   zhi(:,:,:)   = 0.

   ! TOTAL H+ scale: conversion factor for Htot = aphscale * Hfree
      !$OMP PARALLEL DO & 
      !$OMP PRIVATE(mi, mj, jk, zrhd, p_alktot, aphscale, zh_ini, zdelta) & 
      !$OMP SHARED(rmask, rhd, p_alktot, tr, jptal, Kbb, sulfat, aks3, &
                   p_hini, zalknw_inf, zalknw_sup, akw3, zh_min, zh_max, zhi)
      DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF (rmask(mi,mj,jk) == 1.) THEN
            zrhd = 1._wp / ( rho(mi,mj,N+1-jk) + 1. )
            p_alktot = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jptal) * zrhd
            aphscale = 1. + sulfat(mi,mj,jk)/aks3(mi,mj,jk)
            zh_ini = p_hini(mi,mj,jk)

            zdelta = (p_alktot-zalknw_inf(mi,mj,jk))**2 + 4.*akw3(mi,mj,jk)/aphscale

            IF(p_alktot >= zalknw_inf(mi,mj,jk)) THEN
              zh_min(mi,mj,jk) = 2.*akw3(mi,mj,jk) /( p_alktot-zalknw_inf(mi,mj,jk) + SQRT(zdelta) )
            ELSE
              zh_min(mi,mj,jk) = aphscale*(-(p_alktot-zalknw_inf(mi,mj,jk)) + SQRT(zdelta) ) / 2.
            ENDIF

            zdelta = (p_alktot-zalknw_sup(mi,mj,jk))**2 + 4.*akw3(mi,mj,jk)/aphscale

            IF(p_alktot <= zalknw_sup(mi,mj,jk)) THEN
              zh_max(mi,mj,jk) = aphscale*(-(p_alktot-zalknw_sup(mi,mj,jk)) + SQRT(zdelta) ) / 2.
            ELSE
              zh_max(mi,mj,jk) = 2.*akw3(mi,mj,jk) /( p_alktot-zalknw_sup(mi,mj,jk) + SQRT(zdelta) )
            ENDIF

            zhi(mi,mj,jk) = MAX(MIN(zh_max(mi,mj,jk), zh_ini), zh_min(mi,mj,jk))
         ENDIF
      END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO

   zeqn_absmin(:,:,:) = HUGE(1._wp)

   DO jn = 1, jp_maxniter_atgen
   !$OMP PARALLEL DO & 
   DO PRIVATE(mi, mj, jk, zrhd, p_alktot, zdic, zbot, zst, &
           zft, aphscale, zh, zh_prev, zsit, zpt, znumer_dic, &
           zdenom_dic, zalk_dic, zdnumer_dic, zdalk_dic, &
           znumer_bor, zdenom_bor, zalk_bor, zdnumer_bor, zdalk_bor, &
           znumer_po4, zdenom_po4, zalk_po4, zdnumer_po4, zdalk_po4, &
           znumer_sil, zdenom_sil, zalk_sil, zdnumer_sil, zdalk_sil, &
           znumer_so4, zdenom_so4, zalk_so4, zdnumer_so4, zdalk_so4, &
           znumer_flu, zdenom_flu, zalk_flu, zdnumer_flu, zdalk_flu, &
           zalk_wat, zdalk_wat, zeqn, zalka, zdeqndh, zh_min, zh_max, &
           zh_lnfactor, zh_delta, zeqn_absmin, pz_exp_threshold) &
   !$OMP SHARED(rmask, tr, Kbb, borat, sulfat, fluorid, aks3, zhi, ln_p2z, &
                po4r, ak13, ak23, akb3, ak1p3, ak2p3, ak3p3, aksi3, aks3, akw3)
   DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
      IF (rmask(mi,mj,jk) == 1.) THEN
         zrhd = 1._wp / ( rho(mi,mj,N+1-jk) + 1. )
         p_alktot = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jptal) * zrhd
         zdic  = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdic) * zrhd
         zbot  = borat(mi,mj,jk)
         zst = sulfat (mi,mj,jk)
         zft = fluorid(mi,mj,jk)
         aphscale = 1. + sulfat(mi,mj,jk)/aks3(mi,mj,jk)
         zh = zhi(mi,mj,jk)
         zh_prev = zh
         IF( ln_p2z ) THEN
            zsit = 90.33E-6 * zrhd
            zpt  = 2.174E-6 * zrhd
         ELSE
            zpt  = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) * zrhd * po4r
            zsit = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsil) * zrhd
         ENDIF

         ! H2CO3 - HCO3 - CO3 : n=2, m=0
         znumer_dic = 2.*ak13(mi,mj,jk)*ak23(mi,mj,jk) + zh*ak13(mi,mj,jk)
         zdenom_dic = ak13(mi,mj,jk)*ak23(mi,mj,jk) + zh*(ak13(mi,mj,jk) + zh)
         zalk_dic   = zdic * (znumer_dic/zdenom_dic)
         zdnumer_dic = ak13(mi,mj,jk)*ak13(mi,mj,jk)*ak23(mi,mj,jk) + zh     &
                       *(4.*ak13(mi,mj,jk)*ak23(mi,mj,jk) + zh*ak13(mi,mj,jk))
         zdalk_dic   = -zdic*(zdnumer_dic/zdenom_dic**2)


         ! B(OH)3 - B(OH)4 : n=1, m=0
         znumer_bor = akb3(mi,mj,jk)
         zdenom_bor = akb3(mi,mj,jk) + zh
         zalk_bor   = zbot * (znumer_bor/zdenom_bor)
         zdnumer_bor = akb3(mi,mj,jk)
         zdalk_bor   = -zbot*(zdnumer_bor/zdenom_bor**2)


         ! H3PO4 - H2PO4 - HPO4 - PO4 : n=3, m=1
         znumer_po4 = 3.*ak1p3(mi,mj,jk)*ak2p3(mi,mj,jk)*ak3p3(mi,mj,jk)  &
         &            + zh*(2.*ak1p3(mi,mj,jk)*ak2p3(mi,mj,jk) + zh* ak1p3(mi,mj,jk))
         zdenom_po4 = ak1p3(mi,mj,jk)*ak2p3(mi,mj,jk)*ak3p3(mi,mj,jk)     &
         &            + zh*( ak1p3(mi,mj,jk)*ak2p3(mi,mj,jk) + zh*(ak1p3(mi,mj,jk) + zh))
         zalk_po4   = zpt * (znumer_po4/zdenom_po4 - 1.) ! Zero level of H3PO4 = 1
         zdnumer_po4 = ak1p3(mi,mj,jk)*ak2p3(mi,mj,jk)*ak1p3(mi,mj,jk)*ak2p3(mi,mj,jk)*ak3p3(mi,mj,jk)  &
         &             + zh*(4.*ak1p3(mi,mj,jk)*ak1p3(mi,mj,jk)*ak2p3(mi,mj,jk)*ak3p3(mi,mj,jk)         &
         &             + zh*(9.*ak1p3(mi,mj,jk)*ak2p3(mi,mj,jk)*ak3p3(mi,mj,jk)                         &
         &             + ak1p3(mi,mj,jk)*ak1p3(mi,mj,jk)*ak2p3(mi,mj,jk)                                &
         &             + zh*(4.*ak1p3(mi,mj,jk)*ak2p3(mi,mj,jk) + zh * ak1p3(mi,mj,jk) ) ) )
         zdalk_po4   = -zpt * (zdnumer_po4/zdenom_po4**2)

         ! H4SiO4 - H3SiO4 : n=1, m=0
         znumer_sil = aksi3(mi,mj,jk)
         zdenom_sil = aksi3(mi,mj,jk) + zh
         zalk_sil   = zsit * (znumer_sil/zdenom_sil)
         zdnumer_sil = aksi3(mi,mj,jk)
         zdalk_sil   = -zsit * (zdnumer_sil/zdenom_sil**2)

         ! HSO4 - SO4 : n=1, m=1
         aphscale = 1.0 + zst/aks3(mi,mj,jk)
         znumer_so4 = aks3(mi,mj,jk) * aphscale
         zdenom_so4 = aks3(mi,mj,jk) * aphscale + zh
         zalk_so4   = zst * (znumer_so4/zdenom_so4 - 1.)
         zdnumer_so4 = aks3(mi,mj,jk)
         zdalk_so4   = -zst * (zdnumer_so4/zdenom_so4**2)

         ! HF - F : n=1, m=1
         znumer_flu =  akf3(mi,mj,jk)
         zdenom_flu =  akf3(mi,mj,jk) + zh
         zalk_flu   =  zft * (znumer_flu/zdenom_flu - 1.)
         zdnumer_flu = akf3(mi,mj,jk)
         zdalk_flu   = -zft * (zdnumer_flu/zdenom_flu**2)

         ! H2O - OH
         aphscale = 1.0 + zst/aks3(mi,mj,jk)
         zalk_wat   = akw3(mi,mj,jk)/zh - zh/aphscale
         zdalk_wat  = -akw3(mi,mj,jk)/zh**2 - 1./aphscale

         ! CALCULATE [ALK]([CO3--], [HCO3-])
         zeqn = zalk_dic + zalk_bor + zalk_po4 + zalk_sil   &
         &      + zalk_so4 + zalk_flu                       &
         &      + zalk_wat - p_alktot

         zalka = p_alktot - (zalk_bor + zalk_po4 + zalk_sil   &
         &       + zalk_so4 + zalk_flu + zalk_wat)

         zdeqndh = zdalk_dic + zdalk_bor + zdalk_po4 + zdalk_sil &
         &         + zdalk_so4 + zdalk_flu + zdalk_wat

         ! Adapt bracketing interval
         IF(zeqn > 0._wp) THEN
           zh_min(mi,mj,jk) = zh_prev
         ELSEIF(zeqn < 0._wp) THEN
           zh_max(mi,mj,jk) = zh_prev
         ENDIF

         IF(ABS(zeqn) >= 0.5_wp*zeqn_absmin(mi,mj,jk)) THEN
         ! if the function evaluation at the current point is
         ! not decreasing faster than with a bisection step (at least linearly)
         ! in absolute value take one bisection step on [ph_min, ph_max]
         ! ph_new = (ph_min + ph_max)/2d0
         !
         ! In terms of [H]_new:
         ! [H]_new = 10**(-ph_new)
         !         = 10**(-(ph_min + ph_max)/2d0)
         !         = SQRT(10**(-(ph_min + phmax)))
         !         = SQRT(zh_max * zh_min)
            zh = SQRT(zh_max(mi,mj,jk) * zh_min(mi,mj,jk))
            zh_lnfactor = (zh - zh_prev)/zh_prev ! Required to test convergence below
         ELSE
         ! dzeqn/dpH = dzeqn/d[H] * d[H]/dpH
         !           = -zdeqndh * LOG(10) * [H]
         ! \Delta pH = -zeqn/(zdeqndh*d[H]/dpH) = zeqn/(zdeqndh*[H]*LOG(10))
         !
         ! pH_new = pH_old + \deltapH
         !
         ! [H]_new = 10**(-pH_new)
         !         = 10**(-pH_old - \Delta pH)
         !         = [H]_old * 10**(-zeqn/(zdeqndh*[H]_old*LOG(10)))
         !         = [H]_old * EXP(-LOG(10)*zeqn/(zdeqndh*[H]_old*LOG(10)))
         !         = [H]_old * EXP(-zeqn/(zdeqndh*[H]_old))

            zh_lnfactor = -zeqn/(zdeqndh*zh_prev)

            IF(ABS(zh_lnfactor) > pz_exp_threshold) THEN
               zh          = zh_prev*EXP(zh_lnfactor)
            ELSE
               zh_delta    = zh_lnfactor*zh_prev
               zh          = zh_prev + zh_delta
            ENDIF

            IF( zh < zh_min(mi,mj,jk) ) THEN
               ! if [H]_new < [H]_min
               ! i.e., if ph_new > ph_max then
               ! take one bisection step on [ph_prev, ph_max]
               ! ph_new = (ph_prev + ph_max)/2d0
               ! In terms of [H]_new:
               ! [H]_new = 10**(-ph_new)
               !         = 10**(-(ph_prev + ph_max)/2d0)
               !         = SQRT(10**(-(ph_prev + phmax)))
               !         = SQRT([H]_old*10**(-ph_max))
               !         = SQRT([H]_old * zh_min)
               zh                = SQRT(zh_prev * zh_min(mi,mj,jk))
               zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below
            ENDIF

            IF( zh > zh_max(mi,mj,jk) ) THEN
               ! if [H]_new > [H]_max
               ! i.e., if ph_new < ph_min, then
               ! take one bisection step on [ph_min, ph_prev]
               ! ph_new = (ph_prev + ph_min)/2d0
               ! In terms of [H]_new:
               ! [H]_new = 10**(-ph_new)
               !         = 10**(-(ph_prev + ph_min)/2d0)
               !         = SQRT(10**(-(ph_prev + ph_min)))
               !         = SQRT([H]_old*10**(-ph_min))
               !         = SQRT([H]_old * zhmax)
               zh                = SQRT(zh_prev * zh_max(mi,mj,jk))
               zh_lnfactor       = (zh - zh_prev)/zh_prev ! Required to test convergence below
            ENDIF
         ENDIF

         zeqn_absmin(mi,mj,jk) = MIN( ABS(zeqn), zeqn_absmin(mi,mj,jk))

         ! Stop iterations once |\delta{[H]}/[H]| < rdel
         ! <=> |(zh - zh_prev)/zh_prev| = |EXP(-zeqn/(zdeqndh*zh_prev)) -1| < rdel
         ! |EXP(-zeqn/(zdeqndh*zh_prev)) -1| ~ |zeqn/(zdeqndh*zh_prev)|

         ! Alternatively:
         ! |\Delta pH| = |zeqn/(zdeqndh*zh_prev*LOG(10))|
         !             ~ 1/LOG(10) * |\Delta [H]|/[H]
         !             < 1/LOG(10) * rdel

         ! Hence |zeqn/(zdeqndh*zh)| < rdel

         ! rdel <-- pp_rdel_ah_target
         l_exitnow = (ABS(zh_lnfactor) < pp_rdel_ah_target)

         IF(l_exitnow) THEN 
            rmask(mi,mj,jk) = 0.
         ENDIF

         zhi(mi,mj,jk) =  zh

         IF(jn >= jp_maxniter_atgen) THEN
            zhi(mi,mj,jk) = -1._wp
         ENDIF

      ENDIF
   END DO   ;   END DO   ;   END DO
   !$OMP END PARALLEL DO
   END DO
   !

      IF( .false. )   CALL timing_stop('solve_at_general')
      !
   END SUBROUTINE solve_at_general


   INTEGER FUNCTION p4z_che_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_che_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER ::   ierr(4)        ! Local variables
      !!----------------------------------------------------------------------

      ierr(:) = 0

      ALLOCATE( fekeq(Istrp:Iendp,Jstrp:Jendp,N), chemc(Istrp:Iendp,Jstrp:Jendp,3), &
         &      chemo2(Istrp:Iendp,Jstrp:Jendp,N), STAT=ierr(1) )

      ALLOCATE( akb3(Istrp:Iendp,Jstrp:Jendp,N)     , tempis(Istrp:Iendp,Jstrp:Jendp,N),       &
         &      akw3(Istrp:Iendp,Jstrp:Jendp,N)     , borat (Istrp:Iendp,Jstrp:Jendp,N)  ,       &
         &      aks3(Istrp:Iendp,Jstrp:Jendp,N)     , akf3(Istrp:Iendp,Jstrp:Jendp,N)    ,       &
         &      ak1p3(Istrp:Iendp,Jstrp:Jendp,N)    , ak2p3(Istrp:Iendp,Jstrp:Jendp,N)   ,       &
         &      ak3p3(Istrp:Iendp,Jstrp:Jendp,N)    , aksi3(Istrp:Iendp,Jstrp:Jendp,N)   ,       &
         &      fluorid(Istrp:Iendp,Jstrp:Jendp,N)  , sulfat(Istrp:Iendp,Jstrp:Jendp,N)  ,       &
         &      salinprac(Istrp:Iendp,Jstrp:Jendp,N),                 STAT=ierr(2) )

      ALLOCATE( fesol(Istrp:Iendp,Jstrp:Jendp,N,5), STAT=ierr(3) )
      IF( .NOT. ln_p2z ) ALLOCATE( sio3eq(Istrp:Iendp,Jstrp:Jendp,N), STAT=ierr(4) )

      !* Variable for chemistry of the CO2 cycle
      p4z_che_alloc = MAXVAL( ierr )
      !
      IF( p4z_che_alloc /= 0 )   CALL ctl_stop( 'STOP', 'p4z_che_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p4z_che_alloc



   !!======================================================================
END MODULE p4zche
