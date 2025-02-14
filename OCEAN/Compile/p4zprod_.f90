










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









MODULE p4zprod
   !!======================================================================
   !!                         ***  MODULE p4zprod  ***
   !! TOP :  Growth Rate of the two phytoplankton groups of  
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-05  (O. Aumont, C. Ethe) New parameterization of light limitation
   !!----------------------------------------------------------------------
   !!   p4z_prod       : Compute the growth Rate of the two phytoplanktons groups
   !!   p4z_prod_init  : Initialization of the parameters for growth
   !!   p4z_prod_alloc : Allocate variables for growth
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      !  Source Minus Sink variables
   USE p2zlim          ! Co-limitations of different nutrients
   USE p4zlim          ! Co-limitations of differents nutrients
   USE prtctl          ! print control for debugging
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_prod         ! called in p4zbio.F90
   PUBLIC   p4z_prod_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_prod_alloc   ! called in trcini_pisces.F90

   REAL(wp), PUBLIC ::   pislopen     !:  P-I slope of nanophytoplankton
   REAL(wp), PUBLIC ::   pisloped     !:  P-I slope of diatoms
   REAL(wp), PUBLIC ::   excretn      !:  Excretion ratio of nanophyto
   REAL(wp), PUBLIC ::   excretd      !:  Excretion ratio of diatoms
   REAL(wp), PUBLIC ::   bresp        !:  Basal respiration rate
   REAL(wp), PUBLIC ::   chlcmin      !:  Minimum Chl/C ratio of phytoplankton
   REAL(wp), PUBLIC ::   fecnm        !:  Maximum Fe/C ratio of nano
   REAL(wp), PUBLIC ::   fecdm        !:  Maximum Fe/C ratio of diatoms
   REAL(wp), PUBLIC ::   grosip       !:  Mean Si/C ratio of diatoms

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotan   !: proxy of N quota in Nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotad   !: proxy of N quota in diatoms
   
   REAL(wp) ::   r1_rday    ! 1 / rday
   REAL(wp) ::   texcretn   ! 1 - excretn 
   REAL(wp) ::   texcretd   ! 1 - excretd        

   LOGICAL  :: l_dia_pp, l_dia_mu, l_dia_light, l_dia_lprod

   !! * Substitutions












































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zprod.F90 15459 2021-10-29 08:19:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_prod( kt, knt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod  ***
      !!
      !! ** Purpose :   Computes phytoplankton production depending on
      !!                light, temperature and nutrient availability
      !!                Computes also the uptake of Iron and Si as well 
      !!                as the chlorophyll content of the cells
      !!                 relies on a mixed Monod-Quota formalism 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   !
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      !
      INTEGER  ::   mi, mj, jk
      REAL(wp) ::   zsilfac, znanotot, zdiattot
      REAL(wp) ::   zratio, zratio2, zmax, zsilim, zlim, zsiborn
      REAL(wp) ::   zpptot, zpnewtot, zpregtot, zprochln, zprochld
      REAL(wp) ::   zproddoc, zprodsil, zprodfer, zprodlig, zprod1
      REAL(wp) ::   zpislopen, zpisloped, zfact, xksi2_3
      REAL(wp) ::   zratiosi, zratiosi_4, zmaxsi, zlimfac, zlimfac3, zsizetmp, zfecnm, zfecdm
      REAL(wp) ::   zprod, zval, zmxl_fac_nano, zmxl_fac_diat, zmxl_chl, zpronewn, zpronewd
      REAL(wp) ::   zpislopeadn, zpislopeadd
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp, N) :: zprmax, zmxl, zysopt
      !!!       
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp, N) :: zprdia, zprbio, zprchld, zprchln, ratchln, ratchld
      !!!   
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp, N) :: zprorcan, zprorcad
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp, N) :: zprofed, zprofen
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d

      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p4z_prod')
      !
      IF( kt == nittrc000 ) THEN
         l_dia_pp =      iom_use( "PPPHYN" ) .OR. iom_use( "PPPHYD" ) .OR. iom_use( "TPP" )   &
             &      .OR. iom_use( "PPNEWN" ) .OR. iom_use( "PPNEWD" ) .OR. iom_use( "TPNEW")  &
             &      .OR. iom_use( "PFeN"   ) .OR. iom_use( "PFeD"   ) .OR. iom_use( "TPBFE")  &
             &      .OR. iom_use( "PBSi"   ) .OR. iom_use( "PPNEWo2") .OR. iom_use( "PPRego2" )
         l_dia_mu    =   iom_use( "Mumax"  ) .OR. iom_use( "MuN"    ) .OR. iom_use( "MuD")
         l_dia_light =   iom_use( "LNlight") .OR. iom_use( "LDlight")
         l_dia_lprod =   ln_ligand .AND. ( iom_use( "LPRODP") .OR. iom_use( "LDETP") )
      ENDIF

      ! Initialize the local arrays
      zprorcan(:,:,:) = 0._wp
      zprorcad(:,:,:) = 0._wp
      zprofen (:,:,:) = 0._wp
      zprofed (:,:,:) = 0._wp
      zprchld (:,:,:) = 0._wp
      zprchln (:,:,:) = 0._wp 
      zprbio  (:,:,:) = 0._wp
      zprdia  (:,:,:) = 0._wp 
      zmxl    (:,:,:) = 0._wp
      zysopt  (:,:,:) = 0._wp

      ! Computation of the maximimum production. Based on a Q10 description
      ! of the thermal dependency. Parameters are taken from Bissinger et al. (2008)
      zprmax(:,:,:) = 0.65_wp * r1_rday * tgfunc(:,:,:)

      ! Intermittency is supposed to have a similar effect on production as 
      ! day length (Shatwell et al., 2012). The correcting factor is zmxl_fac. 
      ! zmxl_chl is the fractional day length and is used to compute the mean
      ! PAR during daytime. The effect of mixing is computed using the 
      ! absolute light level definition of the euphotic zone
      ! ------------------------------------------------------------------------- 
      IF ( ln_p4z_dcyc ) THEN    ! Diurnal cycle in 
         DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            IF( etot_ndcy(mi,mj,jk) > 1.E-3 ) THEN
               zval = 24.0
               IF( gdepw(mi,mj,jk+1,Kmm) <= hbl(mi,mj) ) THEN
                  zval = zval * MIN(1., heup_01(mi,mj) / ( hbl(mi,mj) + 0.5*EPSILON(1.e0) ))
               ENDIF
               zmxl(mi,mj,jk) = zval
            ENDIF
         END DO   ;   END DO   ;   END DO
      ELSE ! No diurnal cycle in 
         DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            IF( etot_ndcy(mi,mj,jk) > 1.E-3 ) THEN
               zval = MAX( 1., strn(mi,mj) )
               IF( gdepw(mi,mj,jk+1,Kmm) <= hbl(mi,mj) ) THEN
                  zval = zval * MIN(1., heup_01(mi,mj) / ( hbl(mi,mj) + 0.5*EPSILON(1.e0) ))
               ENDIF
               zmxl(mi,mj,jk) = zval
            ENDIF
         END DO   ;   END DO   ;   END DO
      ENDIF

      ! The formulation proposed by Geider et al. (1997) has been modified 
      ! to exclude the effect of nutrient limitation and temperature in the PI
      ! curve following Vichi et al. (2007)
      ! -----------------------------------------------------------------------
      DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF( etot_ndcy(mi,mj,jk) > 1.E-3 ) THEN
            zmxl_fac_nano = 1.0 - EXP( -0.25 * MAX(0., zmxl(mi,mj,jk) - 2. ) )
            zmxl_fac_diat = 1.0 - EXP( -0.25 * MAX(0., zmxl(mi,mj,jk) - 1. ) )
            zmxl_chl = zmxl(mi,mj,jk) / 24.
            zprbio(mi,mj,jk) = zprmax(mi,mj,jk) * zmxl_fac_nano
            zprdia(mi,mj,jk) = zprmax(mi,mj,jk) * zmxl_fac_diat
            !
            ! The initial slope of the PI curve could be increased for nano
            ! to account for photadaptation, for instance in the DCM

            ! Nanophytoplankton
            zpislopeadn = pislopen * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnch)   &
              &                    /( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) * 12. + 0.5*EPSILON(1.e0))

            ! Diatoms
            zpislopeadd = pisloped * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdch)   &
              &                    /( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) * 12. + 0.5*EPSILON(1.e0))
            !
            ! Computation of production function for Carbon
            ! Actual light levels are used here 
            ! ----------------------------------------------
            zpislopen = zpislopeadn / ( zprmax(mi,mj,jk) * xlimphy(mi,mj,jk) &
              &           * zmxl_fac_nano * day2sec + 0.5*EPSILON(1.e0))
            zpisloped = zpislopeadd / ( zprmax(mi,mj,jk) * xlimdia(mi,mj,jk) &
              &           * zmxl_fac_diat * day2sec + 0.5*EPSILON(1.e0))
            zprbio(mi,mj,jk) = zprbio(mi,mj,jk) * ( 1.- EXP( -zpislopen * enano(mi,mj,jk) )  )
            zprdia(mi,mj,jk) = zprdia(mi,mj,jk) * ( 1.- EXP( -zpisloped * ediat(mi,mj,jk) )  )

            !  Computation of production function for Chlorophyll
            !  Mean light level in the mixed layer (when appropriate)
            !  is used here (acclimation is in general slower than 
            !  the characteristic time scales of vertical mixing)
            !  ------------------------------------------------------
            zpislopen = zpislopen * zmxl_fac_nano / ( zmxl_chl + 0.5*EPSILON(1.e0) )
            zpisloped = zpisloped * zmxl_fac_diat / ( zmxl_chl + 0.5*EPSILON(1.e0) )
            zprchln(mi,mj,jk) = ( 1.- EXP( -zpislopen * enanom(mi,mj,jk) ) )
            zprchld(mi,mj,jk) = ( 1.- EXP( -zpisloped * ediatm(mi,mj,jk) ) )
         ENDIF
      END DO   ;   END DO   ;   END DO

      !  Computation of a proxy of the N/C quota from nutrient limitation 
      !  and light limitation. Steady state is assumed to allow the computation
      !  ----------------------------------------------------------------------
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
          zval = MIN( xnanopo4(mi,mj,jk), ( xnanonh4(mi,mj,jk) + xnanono3(mi,mj,jk) ) )   &
            &      * zprmax(mi,mj,jk) / ( zprbio(mi,mj,jk) + 0.5*EPSILON(1.e0) )
          quotan(mi,mj,jk) = MIN( 1., 0.3 + 0.7 * zval )
          zval = MIN( xdiatpo4(mi,mj,jk), ( xdiatnh4(mi,mj,jk) + xdiatno3(mi,mj,jk) ) )   &
            &      * zprmax(mi,mj,jk) / ( zprdia(mi,mj,jk) + 0.5*EPSILON(1.e0) )
          quotad(mi,mj,jk) = MIN( 1., 0.3 + 0.7 * zval )
      END DO   ;   END DO   ;   END DO


      DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
          IF( etot_ndcy(mi,mj,jk) > 1.E-3 ) THEN
             ! Si/C of diatoms
             ! ------------------------
             ! Si/C increases with iron stress and silicate availability
             ! Si/C is arbitrariliy increased for very high Si concentrations
             ! to mimic the very high ratios observed in the Southern Ocean (zsilfac)
             ! A parameterization derived from Flynn (2003) is used for the control
             ! when Si is not limiting which is similar to the parameterisation
             ! proposed by Gurney and Davidson (1999).
             ! -----------------------------------------------------------------------
            zlim  = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsil) &
               &      / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsil) + xksi1 )
            zsilim = xlimdia(mi,mj,jk) * zprdia(mi,mj,jk) / ( zprmax(mi,mj,jk) + 0.5*EPSILON(1.e0) )
            zsiborn =  t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsil) &
               &     * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsil) &
               &     * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsil)
            IF (latr(mi,mj) < -30 ) THEN
              zsilfac = 1. + 2. * zsiborn / ( zsiborn + xksi2_3 )
            ELSE
              zsilfac = 1. +      zsiborn / ( zsiborn + xksi2_3 )
            ENDIF
            zratiosi = 1.0 - t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdsi) &
               &     / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) + 0.5*EPSILON(1.e0) ) &
               &     / ( zsilfac * grosip * 3.0 + 0.5*EPSILON(1.e0) )
            zratiosi = MAX(0., MIN(1.0, zratiosi) )
            zratiosi_4 = zratiosi * zratiosi * zratiosi * zratiosi 
            zmaxsi     = ( 1.0 + 1.E-4 ) * zratiosi_4 / ( zratiosi_4 + 1.E-4 )
            IF( xlimsi(mi,mj,jk) /= xlimdia(mi,mj,jk) ) THEN
               zysopt(mi,mj,jk) = zlim * zsilfac * grosip * 1.0 * zmaxsi
            ELSE
               zysopt(mi,mj,jk) = zlim * zsilfac * grosip * 1.0 * zsilim**0.7 * zmaxsi
            ENDIF
        ENDIF
      END DO   ;   END DO   ;   END DO

      ! Sea-ice effect on production
      ! No production is assumed below sea ice
      ! -------------------------------------- 
      DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zprbio(mi,mj,jk) = zprbio(mi,mj,jk) * ( 1. - fr_i(mi,mj) )
         zprdia(mi,mj,jk) = zprdia(mi,mj,jk) * ( 1. - fr_i(mi,mj) )
      END DO   ;   END DO   ;   END DO

      !
      ! Size computation
      ! Size is made a function of the limitation of of phytoplankton growth
      ! Strongly limited cells are supposed to be smaller. sizena is the 
      ! size at time step t+1 and is thus updated at the end of the 
      ! current time step
      ! --------------------------------------------------------------------

      DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         zlimfac  = MIN(1.0, xlimphy(mi,mj,jk), zprbio(mi,mj,jk) &
              &    / ( zprmax(mi,mj,jk) + 0.5*EPSILON(1.e0) ) )
         zsizetmp = MIN( xsizern, 1.0 + xsizern * (1.0 - COS(3.141/2. * MAX(0., zlimfac - 0.2) / 0.8)) )
         sizena(mi,mj,jk) = sizen(mi,mj,jk) + MAX(zprbio(mi,mj,jk), 1E-2 / day2sec ) &
                   &        * xlimphy(mi,mj,jk) * ( 1.0 - xlimphy(mi,mj,jk) )   &
                   &        * rfact2 * ( zsizetmp - sizen(mi,mj,jk) )
         sizena(mi,mj,jk) = MIN( xsizern, sizena(mi,mj,jk) )
         !
         zlimfac  = MIN( 1.0, xlimdia(mi,mj,jk), zprdia(mi,mj,jk) / ( zprmax(mi,mj,jk) + 0.5*EPSILON(1.e0) ) )
         zsizetmp = MIN( xsizerd, 1.0 + xsizerd * (1.0 - COS(3.141/2. * MAX(0., zlimfac - 0.2) / 0.8)) )
         sizeda(mi,mj,jk) = sized(mi,mj,jk) + MAX(zprdia(mi,mj,jk), 1E-2 / day2sec ) &
                   &        * xlimdia(mi,mj,jk) * ( 1.0 - xlimdia(mi,mj,jk) )   &
                   &        * rfact2 * ( zsizetmp - sized(mi,mj,jk) )
         sizeda(mi,mj,jk) = MIN( xsizerd, sizeda(mi,mj,jk) )
      END DO   ;   END DO   ;   END DO


      ! Computation of the various production  and nutrient uptake terms
      ! ---------------------------------------------------------------
      DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF( etot_ndcy(mi,mj,jk) > 1.E-3 ) THEN
              !  production terms for nanophyto. (C)
            zprorcan(mi,mj,jk) = zprbio(mi,mj,jk) * xlimphy(mi,mj,jk) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) * rfact2

            ! Iron uptake rates of nanophytoplankton. Upregulation is  
            ! not parameterized at low iron concentrations as observations
            ! do not suggest it for accimated cells. Uptake is
            ! downregulated when the quota is close to the maximum quota
            zfecnm   = xqfuncfecn(mi,mj,jk) + ( fecnm - xqfuncfecn(mi,mj,jk) ) * ( xnanono3(mi,mj,jk) + xnanonh4(mi,mj,jk) )
            zratio   = 1.0 - MIN(1.0,t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnfe) &
                &      / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) * zfecnm + 0.5*EPSILON(1.e0) ) )
            zratio2  = zratio * zratio
            zmax     = zratio2 / (0.0025+zratio2)
            zprofen(mi,mj,jk) = zfecnm * zprmax(mi,mj,jk) * ( 1.0 - fr_i(mi,mj) )        &
              &          * (1. + 0.8 * xnanono3(mi,mj,jk) / ( 0.5*EPSILON(1.e0) + xnanono3(mi,mj,jk)  &
              &          + xnanonh4(mi,mj,jk) ) * (1. - xnanofer(mi,mj,jk) ) )           &
              &          * xnanofer(mi,mj,jk) * zmax * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) * rfact2
            ! production terms of diatoms (C)
            zprorcad(mi,mj,jk) = zprdia(mi,mj,jk) * xlimdia(mi,mj,jk) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) * rfact2

            ! Iron uptake rates of diatoms. Upregulation is  
            ! not parameterized at low iron concentrations as observations
            ! do not suggest it for accimated cells. Uptake is
            ! downregulated when the quota is close to the maximum quota
            zfecdm   = xqfuncfecd(mi,mj,jk) + ( fecdm - xqfuncfecd(mi,mj,jk) ) &
                 &    * ( xdiatno3(mi,mj,jk) + xdiatnh4(mi,mj,jk) )
            zratio   = 1.0 - MIN(1.0, t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdfe) &
                 &   / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) * zfecdm + 0.5*EPSILON(1.e0) ) )
            zratio2  = zratio * zratio 
            zmax     = zratio2/ (0.0025+zratio2)
            zprofed(mi,mj,jk) = zfecdm * zprmax(mi,mj,jk) * (1.0 - fr_i(mi,mj) )         &
              &          * (1. + 0.8 * xdiatno3(mi,mj,jk) / ( 0.5*EPSILON(1.e0) + xdiatno3(mi,mj,jk)  &
              &          + xdiatnh4(mi,mj,jk) ) * (1. - xdiatfer(mi,mj,jk) ) )           &
              &          * xdiatfer(mi,mj,jk) * zmax * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) * rfact2
               
         ENDIF
      END DO   ;   END DO   ;   END DO

      ! Computation of the chlorophyll production terms
      ! The parameterization is taken from Geider et al. (1997)
      ! -------------------------------------------------------
      DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         IF( etot_ndcy(mi,mj,jk) > 1.E-3 ) THEN
            zmxl_chl = zmxl(mi,mj,jk)  / 24.
            !  production terms for nanophyto. ( chlorophyll )
            znanotot = enanom(mi,mj,jk) / ( zmxl_chl + 0.5*EPSILON(1.e0) )
            zprod1   = zprorcan(mi,mj,jk) * texcretn / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) + 0.5*EPSILON(1.e0) )
            zprod    = zprod1 / ratchln(mi,mj,jk) &
              &        * ( pislopen * znanotot / ( zprmax(mi,mj,jk) * day2sec ) &
              &        * ( 1.0 - zprchln(mi,mj,jk) ) &
              &        * MAX(0.0, (1.0 - ratchln(mi,mj,jk) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnch)    &
              &        / ( 12. * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) * xlimphy(mi,mj,jk) + 0.5*EPSILON(1.e0) ) ) ) &
              &        - ratchln(mi,mj,jk) * zprchln(mi,mj,jk) ) + zprod1
            zprochln = MAX(zprod * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnch) , chlcmin * 12 * zprorcan(mi,mj,jk) )

            !  production terms for diatoms ( chlorophyll )
            zdiattot = ediatm(mi,mj,jk) / ( zmxl_chl + 0.5*EPSILON(1.e0) )
            zprod1   = zprorcad(mi,mj,jk) * texcretd / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) + 0.5*EPSILON(1.e0) )
            zprod    = zprod1 / ratchld(mi,mj,jk) &
              &        * ( pisloped * zdiattot / ( zprmax(mi,mj,jk) * day2sec ) &
              &        * ( 1.0 - zprchld(mi,mj,jk) ) &
              &        * MAX(0.0, (1.0 - ratchld(mi,mj,jk) * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdch)    &
              &        / ( 12. * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) * xlimdia(mi,mj,jk) + 0.5*EPSILON(1.e0) ) ) ) &
              &        - ratchld(mi,mj,jk) * zprchld(mi,mj,jk) ) + zprod1
            zprochld = MAX(zprod * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdch) , chlcmin * 12 * zprorcad(mi,mj,jk) )

            !   Update the arrays TRA which contain the Chla sources and sinks
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnch) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnch) + zprochln
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdch) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdch) + zprochld
         ENDIF
      END DO   ;   END DO   ;   END DO

      !   Update the arrays TRA which contain the biological sources and sinks
      DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
        IF( etot_ndcy(mi,mj,jk) > 1.E-3 ) THEN
           ! New production (uptake of NO3)
           zpronewn = zprorcan(mi,mj,jk) * xnanono3(mi,mj,jk) &
                &   / ( xnanono3(mi,mj,jk) + xnanonh4(mi,mj,jk) + 0.5*EPSILON(1.e0) )
           zpronewd = zprorcad(mi,mj,jk) * xdiatno3(mi,mj,jk) &
                &   / ( xdiatno3(mi,mj,jk) + xdiatnh4(mi,mj,jk) + 0.5*EPSILON(1.e0) )
           !
           zpptot   = zprorcan(mi,mj,jk) + zprorcad(mi,mj,jk)
           zpnewtot = zpronewn + zpronewd
           zpregtot = zpptot - zpnewtot
           zprodsil  = zprmax(mi,mj,jk) * zysopt(mi,mj,jk) * rfact2 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia)
           zproddoc  = excretd * zprorcad(mi,mj,jk) + excretn * zprorcan(mi,mj,jk)
           zprodfer  = texcretn * zprofen(mi,mj,jk) + texcretd * zprofed(mi,mj,jk)
           !
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jppo4) - zpptot
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpno3) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpno3) - zpnewtot
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnh4) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnh4) - zpregtot
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpphy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpphy) + zprorcan(mi,mj,jk) * texcretn
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpnfe) + zprofen(mi,mj,jk)  * texcretn
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdia) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdia) + zprorcad(mi,mj,jk) * texcretd
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdfe) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdfe) + zprofed(mi,mj,jk)  * texcretd
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdsi) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdsi) + zprodsil
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsil) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpsil) - zprodsil
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdoc) + zproddoc
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpoxy) &
                 &   + o2ut * zpregtot + ( o2ut + o2nit ) * zpnewtot
           !
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) - zprodfer
           consfe3(mi,mj,jk)   = zprodfer * 75.0 &
           &       / ( 0.5*EPSILON(1.e0) + ( plig(mi,mj,jk) + 75.0 * (1.0 - plig(mi,mj,jk) ) )   &
           &                   * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpfer) ) / rfact2
           !
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpdic) - zpptot
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jptal) + rno3 * ( zpnewtot - zpregtot )
        ENDIF
      END DO   ;   END DO   ;   END DO

     ! Production and uptake of ligands by phytoplankton. This part is activated 
     ! when ln_ligand is set to .true. in the namelist. Ligand uptake is small 
     ! and based on the FeL model by Morel et al. (2008) and on the study of
     ! Shaked et al. (2020)
     ! -------------------------------------------------------------------------
     IF( ln_ligand ) THEN
         DO jk= 1, nksr ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           IF( etot_ndcy(mi,mj,jk) > 1.E-3 ) THEN
              zproddoc = excretd * zprorcad(mi,mj,jk) + excretn * zprorcan(mi,mj,jk)
              zprodfer = texcretn * zprofen(mi,mj,jk) + texcretd * zprofed(mi,mj,jk)
              zprodlig = plig(mi,mj,jk) &
                 &      / ( 0.5*EPSILON(1.e0) + plig(mi,mj,jk) + 75.0 * (1.0 - plig(mi,mj,jk) ) ) * lthet
              !
              t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jplgw) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jplgw) &
                 &              + zproddoc * ldocp - zprodfer * zprodlig
           ENDIF
         END DO   ;   END DO   ;   END DO
     ENDIF

    ! Total primary production per year
    IF( l_dia_pp )  THEN
       ALLOCATE( zw3d(Istrp:Iendp,Jstrp:Jendp,N) )  ;  zw3d(Istrp:Iendp,Jstrp:Jendp,N) = 0._wp
       DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
          zw3d(mi,mj,jk) = ( zprorcan(mi,mj,jk) + zprorcad(mi,mj,jk) ) * cvol(mi,mj,jk)
       END DO   ;   END DO   ;   END DO
       tpp = glob_sum( 'p4zprod', zw3d )
       DEALLOCATE ( zw3d )
    ENDIF
    
    IF( .false. .AND.  knt == nrdttrc ) THEN
       !
       IF( l_dia_pp ) THEN
          ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          ! primary production by nanophyto
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = zprorcan(mi,mj,jk) * zfact * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "PPPHYN", zw3d )  
          ! primary production by diatomes
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = zprorcad(mi,mj,jk) * zfact * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "PPPHYD", zw3d )  
          ! total primary production
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = ( zprorcan(mi,mj,jk) + zprorcad(mi,mj,jk) ) &
                  &           * zfact * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "TPP", zw3d )  
          CALL iom_put( "PPNEWo2", zw3d * ( o2ut + o2nit ) ) ! Oxygen production by the New Produc
          CALL iom_put( "tintpp"  , tpp * zfact )  !  global total integrated primary production molC/s
          ! new primary production by nano
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = ( zprorcan(mi,mj,jk) + xnanono3(mi,mj,jk)  &
                &              / ( xnanono3(mi,mj,jk) + xnanonh4(mi,mj,jk) + 0.5*EPSILON(1.e0) ) ) &
                &           * zfact * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "PPNEWN", zw3d )  
          ! new primary production by diatomes
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = ( zprorcad(mi,mj,jk) + xdiatno3(mi,mj,jk)  &
                &              / ( xdiatno3(mi,mj,jk) + xdiatnh4(mi,mj,jk) + 0.5*EPSILON(1.e0) ) ) &
                &           * zfact * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "PPNEWD", zw3d )  
          ! total new production 
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = ( ( zprorcan(mi,mj,jk) + xnanono3(mi,mj,jk)  &
                &              / ( xnanono3(mi,mj,jk) + xnanonh4(mi,mj,jk) + 0.5*EPSILON(1.e0) ) ) &
                &              +  ( zprorcad(mi,mj,jk) + xdiatno3(mi,mj,jk)  &
                &              / ( xdiatno3(mi,mj,jk) + xdiatnh4(mi,mj,jk) + 0.5*EPSILON(1.e0) ) ) ) &
                &           * zfact * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "TPNEW", zw3d )  
          ! Regenerated production 
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = ( ( zprorcan(mi,mj,jk) + xnanonh4(mi,mj,jk)  &
                &              / ( xnanono3(mi,mj,jk) + xnanonh4(mi,mj,jk) + 0.5*EPSILON(1.e0) ) ) &
                &              +  ( zprorcad(mi,mj,jk) + xdiatnh4(mi,mj,jk)  &
                &              / ( xdiatno3(mi,mj,jk) + xdiatnh4(mi,mj,jk) + 0.5*EPSILON(1.e0) ) ) ) &
                &           * o2ut * zfact * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "PPRego2", zw3d )  
          !  biogenic silica production
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = zprorcad(mi,mj,jk) * zysopt(mi,mj,jk) &
                  &         * zfact * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "PBSi", zw3d )  
          ! biogenic iron production by nanophyto
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = zprofen(mi,mj,jk) * zfact * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "PFeN", zw3d )  
          ! biogenic iron production by diatomes
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = zprofed(mi,mj,jk) * zfact * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "PFeD", zw3d )  
          ! total biogenic iron production
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = ( zprofen(mi,mj,jk) + zprofed(mi,mj,jk) ) &
                  &           * zfact * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "TPBFE", zw3d )  
          !
          DEALLOCATE ( zw3d ) 
       ENDIF
       !
       IF( l_dia_mu ) THEN
          ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = zprmax(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "Mumax", zw3d )  
          ! Realized growth rate for nanophyto
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = zprbio(mi,mj,jk) * xlimphy(mi,mj,jk) &
                &             * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "MuN", zw3d )  
          ! Realized growth rate for diatoms
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = zprdia(mi,mj,jk) * xlimdia(mi,mj,jk) &
                &             * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "MuD", zw3d )  
          DEALLOCATE ( zw3d ) 
       ENDIF
       !
       IF( l_dia_light ) THEN
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
          ! light limitation term for nano
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = zprbio(mi,mj,jk) /( zprmax(mi,mj,jk) + 0.5*EPSILON(1.e0) ) &
                &             * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "LNlight", zw3d )  
          ! light limitation term for diatomes
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = zprdia(mi,mj,jk) /( zprmax(mi,mj,jk) + 0.5*EPSILON(1.e0) ) &
                &             * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "LDlight", zw3d )  
          DEALLOCATE ( zw3d ) 
       ENDIF
       !
       IF( l_dia_lprod ) THEN
          ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = ( excretd * zprorcad(mi,mj,jk) + excretn * zprorcan(mi,mj,jk) ) &
                &             * zfact * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "LPRODP"  , zw3d * ldocp * 1e9 )
          !
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = ( texcretd * zprofed(mi,mj,jk) + texcretn * zprofed(mi,mj,jk) ) &
            &                 * plig(mi,mj,jk) / ( 0.5*EPSILON(1.e0) + plig(mi,mj,jk) &
            &                                     + 75.0 * (1.0 - plig(mi,mj,jk) ) )  &
                &             * zfact * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "LDETP"   , zw3d * lthet * 1e9 )
          DEALLOCATE ( zw3d ) 
       ENDIF
       !
     ENDIF
      !   Supplementary diagnostics
     zfact = 1.e3 * rfact2r
     DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
        zfact = zfact * tmask(mi,mj,jk)
        bioFlux(mi,mj,N+1-jk,Nprorca  )  = zprorcan(mi,mj,jk) * zfact  ! primary production by nanophyto
        bioFlux(mi,mj,N+1-jk,Nprorcad )  = zprorcad(mi,mj,jk) * zfact  ! primary production by diatom
        bioFlux(mi,mj,N+1-jk,Npronew  )  = ( ( zprorcan(mi,mj,jk) * xnanono3(mi,mj,jk) ) &
           &                      / ( xnanono3(mi,mj,jk) + xnanonh4(mi,mj,jk) + 0.5*EPSILON(1.e0) ) ) * zfact ! new primary production by nanophyto
        bioFlux(mi,mj,N+1-jk,Npronewd )  = ( ( zprorcad(mi,mj,jk) * xdiatno3(mi,mj,jk) ) &
           &                      / ( xdiatno3(mi,mj,jk) + xdiatnh4(mi,mj,jk) + 0.5*EPSILON(1.e0) ) ) * zfact ! new primary production by nanophyto
        bioFlux(mi,mj,N+1-jk,Nprobsi  )  = zprorcad(mi,mj,jk) * zysopt(mi,mj,jk) * zfact ! biogenic silica production
        bioFlux(mi,mj,N+1-jk,Nprofed  )  = zprofed (mi,mj,jk) * zfact  ! biogenic iron production by diatom
        bioFlux(mi,mj,N+1-jk,Nprofen  )  = zprofen (mi,mj,jk) * zfact !  biogenic iron production by nanophyto

        bioFlux(mi,mj,N+1-jk,Npronewo2)  = ( o2ut + o2nit ) &  ! Oxygen production by the New Produc.
           &                        * ( zprorcan(mi,mj,jk) + zprorcad(mi,mj,jk) ) * zfact

        bioFlux(mi,mj,N+1-jk,Nprorego2  ) =  ( zprorcan(mi,mj,jk) * xnanono3(mi,mj,jk) ) &
           &                      / ( xnanono3(mi,mj,jk) + xnanonh4(mi,mj,jk) + 0.5*EPSILON(1.e0) )  &
           &                      + (  zprorcad(mi,mj,jk) * xdiatno3(mi,mj,jk) ) &
           &                      / ( xdiatno3(mi,mj,jk) + xdiatnh4(mi,mj,jk) + 0.5*EPSILON(1.e0) )  &
           &                      * o2ut * zfact ! Oxygen production by the Regen Produc
     END DO   ;   END DO   ;   END DO

     IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('prod')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
    !     CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
     ENDIF
      !
      IF( .false. )  CALL timing_stop('p4z_prod')
      !
   END SUBROUTINE p4z_prod


   SUBROUTINE p4z_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_prod_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton production parameters
      !!
      !! ** Method  :   Read the namp4zprod namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namp4zprod
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      ! Namelist block
      NAMELIST/namp4zprod/ pislopen, pisloped, bresp, excretn, excretd,  &
         &                 chlcmin, fecnm, fecdm, grosip
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*)
         WRITE(stdout,*) 'p4z_prod_init : phytoplankton growth'
         WRITE(stdout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND(numnatp_ref);READ(numnatp_ref,namp4zprod,IOSTAT=ios);CALL ctl_nam(ios,"namp4zprod (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,namp4zprod,IOSTAT=ios);CALL ctl_nam(ios,"namp4zprod (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE( numonp, namp4zprod )

      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*) '   Namelist : namp4zprod'
         WRITE(stdout,*) '      mean Si/C ratio                           grosip       =', grosip
         WRITE(stdout,*) '      P-I slope                                 pislopen     =', pislopen
         WRITE(stdout,*) '      excretion ratio of nanophytoplankton      excretn      =', excretn
         WRITE(stdout,*) '      excretion ratio of diatoms                excretd      =', excretd
         WRITE(stdout,*) '      basal respiration in phytoplankton        bresp        =', bresp
         WRITE(stdout,*) '      Maximum Chl/C in phytoplankton            chlcmin      =', chlcmin
         WRITE(stdout,*) '      P-I slope  for diatoms                    pisloped     =', pisloped
         WRITE(stdout,*) '      Maximum Fe/C in nanophytoplankton         fecnm        =', fecnm
         WRITE(stdout,*) '      Minimum Fe/C in diatoms                   fecdm        =', fecdm
      ENDIF
      !
      r1_rday   = 1._wp / day2sec 
      texcretn  = 1._wp - excretn
      texcretd  = 1._wp - excretd
      tpp       = 0._wp
      !
   END SUBROUTINE p4z_prod_init

   INTEGER FUNCTION p4z_prod_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( quotan(Istrp:Iendp,Jstrp:Jendp,N), quotad(Istrp:Iendp,Jstrp:Jendp,N), STAT = p4z_prod_alloc )
      !
      IF( p4z_prod_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p4z_prod_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p4z_prod_alloc


   !!======================================================================
END MODULE p4zprod
