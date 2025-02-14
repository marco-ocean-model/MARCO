










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









MODULE p4zflx
   !!======================================================================
   !!                         ***  MODULE p4zflx  ***
   !! TOP :    CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
   !!======================================================================
   !! History :   -   !  1988-07  (E. MAIER-REIMER) Original code
   !!             -   !  1998     (O. Aumont) additions
   !!             -   !  1999     (C. Le Quere) modifications
   !!            1.0  !  2004     (O. Aumont) modifications
   !!            2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!                 !  2011-02  (J. Simeon, J. Orr) Include total atm P correction 
   !!            4.2  !  2020     (J. ORR )  rhop is replaced by "in situ density" rhd
   !!            3.*  !  2025-02  (S. Maishal, R. Person)  and optimization
   !!----------------------------------------------------------------------
   !!   p4z_flx       :   CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
   !!   p4z_flx_init  :   Read the namelist
   !!   p4z_patm      :   Read sfc atm pressure [atm] for each grid cell
   !!----------------------------------------------------------------------
   USE oce_trc        !  shared variables between ocean and passive tracers 
   USE trc            !  passive tracers common variables
   USE sms_pisces     !   Source Minus Sink variables
   USE p4zche         !  Chemical model
   USE prtctl         !  print control for debugging
   USE iom            !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_flx  
   PUBLIC   p4z_flx_init  
   PUBLIC   p4z_flx_alloc  

   !                                  !! ** Namelist  nampisext  **
   REAL(wp)          ::   atcco2      !: pre-industrial atmospheric [co2] (ppm) 
   LOGICAL           ::   ln_co2int   !: flag to read in a file and interpolate atmospheric pco2 or not
   CHARACTER(len=34) ::   clname      !: filename of pco2 values
   INTEGER           ::   nn_offset   !: Offset model-data start year (default = 0) 

   !!  Variables related to reading atmospheric CO2 time history    
   INTEGER                                   ::   nmaxrec, numco2   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:) ::   atcco2h, years    !

   !                                   !!* nampisatm namelist (Atmospheric PRessure) *
   LOGICAL, PUBLIC ::   ln_presatm     !: ref. pressure: global mean Patm (F) or a constant (F)
   LOGICAL, PUBLIC ::   ln_presatmco2  !: accounting for spatial atm CO2 in the compuation of carbon flux (T) or not (F)

   REAL(wp) , ALLOCATABLE, SAVE, DIMENSION(:,:) ::   patm      ! atmospheric pressure at kt                 [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  satmco2   !: atmospheric pco2 

   REAL(wp) ::   xconv  = 0.01_wp / 3600._wp   !: coefficients for conversion 

   LOGICAL  :: l_dia_cflx, l_dia_tcflx
   LOGICAL  :: l_dia_oflx, l_dia_kg

   !! * Substitutions












































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zflx.F90 15532 2021-11-24 11:47:32Z techene $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_flx ( kt, knt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_flx  ***
      !!
      !! ** Purpose :   CALCULATES GAS EXCHANGE AND CHEMISTRY AT SEA SURFACE
      !!
      !! ** Method  : 
      !!              - Include total atm P correction via Esbensen & Kushnir (1981) 
      !!              - Remove Wanninkhof chemical enhancement;
      !!              - Add option for time-interpolation of atcco2.txt  
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   !
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs      ! time level indices
      !
      INTEGER  ::   mi, mj, jm, iind, iindm1
      REAL(wp) ::   ztc, ztc2, ztc3, ztc4, zws, zkgwan
      REAL(wp) ::   zfld, zflu, zfld16, zflu16, zrhd
      REAL(wp) ::   zvapsw, zsal, zfco2, zxc2, xCO2approx, ztkel, zfugcoeff
      REAL(wp) ::   zph, zdic, zsch_o2, zsch_co2
      REAL(wp) ::   zyr_dec, zdco2dt
      CHARACTER (len=25) ::   charout
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp) ::   zkgco2, zkgo2, zh2co3, zoflx,  zpco2atm, zpco2oce  
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zw2d
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p4z_flx')
      !
      IF( kt == nittrc000 )  THEN
         l_dia_cflx  = iom_use( "Cflx" ) .OR. iom_use( "Dpco2" ) &
            &     .OR. iom_use( "pCO2sea" ) .OR. iom_use( "AtmCo2" )
         l_dia_oflx  = iom_use( "Oflx" ) .OR. iom_use( "Dpo2" )  
         l_dia_tcflx = iom_use( "tcflx" ) .OR. iom_use( "tcflxcum" )
         l_dia_kg    = iom_use( "Kg" ) 
      ENDIF
      
      ! SURFACE CHEMISTRY (PCO2 AND [H+] IN
      !     SURFACE LAYER); THE RESULT OF THIS CALCULATION
      !     IS USED TO COMPUTE AIR-SEA FLUX OF CO2

      IF( kt /= ntstart .AND. .NOT.l_co2cpl .AND. knt == 1 )   CALL p4z_patm( kt )   ! Get sea-level pressure (E&K [1981] climatology) for use in flux calcs

      IF( ln_co2int .AND. .NOT.ln_presatmco2 .AND. .NOT.l_co2cpl ) THEN 
         ! Linear temporal interpolation  of atmospheric pco2.  atcco2.txt has annual values.
         ! Caveats: First column of .txt must be in years, decimal  years preferably. 
         ! For nn_offset, if your model year is iyy, nn_offset=(years(1)-iyy) 
         ! then the first atmospheric CO2 record read is at years(1)
         zyr_dec = REAL( (int(tdays*day2year*year2day)) + nn_offset, wp ) + REAL( (int(tdays)+1), wp ) / REAL( year2day, wp )
         jm = 1
         DO WHILE( jm <= nmaxrec .AND. years(jm) < zyr_dec ) ;  jm = jm + 1 ;  END DO
         iind = jm  ;   iindm1 = jm - 1
         zdco2dt = ( atcco2h(iind) - atcco2h(iindm1) ) / ( years(iind) - years(iindm1) + 0.5*EPSILON(1.e0) )
         atcco2  = zdco2dt * ( zyr_dec - years(iindm1) ) + atcco2h(iindm1)
         satmco2(:,:) = atcco2 
      ENDIF

      IF( l_co2cpl ) THEN
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            satmco2(mi,mj) = atm_co2(mi,mj)
         END DO   ;   END DO
      ENDIF

      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, zrhd, zdic, zph) &
      !$OMP SHARED(rhd, tr, hi, 0.5*EPSILON(1.e0), ak13, ak23, zh2co3)
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ! (DUMMY) VARIABLES FOR DIC, H+, AND BORATE
         zrhd = rho(mi,mj,N+1-1) + 1._wp
         zdic  = t(mi,mj,N+1-1,Kbb,itemp+ntrc_salt+jpdic)
         zph   = MAX( hi(mi,mj,1), 1.e-10 ) / ( zrhd + 0.5*EPSILON(1.e0) )
         ! CALCULATE [H2CO3]
         zh2co3(mi,mj) = zdic/(1. + ak13(mi,mj,1)/zph + ak13(mi,mj,1)*ak23(mi,mj,1)/zph**2)
      END DO   ;   END DO
      !$OMP END PARALLEL DO

      ! --------------
      ! COMPUTE FLUXES
      ! --------------

      ! FIRST COMPUTE GAS EXCHANGE COEFFICIENTS
      ! ---------------------------------------

      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, ztc, ztc2, ztc3, ztc4, zsch_co2, zsch_o2, zws, zkgwan) &
      !$OMP SHARED(ts, itemp, Kmm, wndm, fr_i, tmask, xconv, zkgco2, zkgo2)
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ztc  = MIN( 35., t(mi,mj,N+1-1,nnew,itemp) )
         ztc2 = ztc * ztc
         ztc3 = ztc * ztc2 
         ztc4 = ztc2 * ztc2 
         ! Compute the schmidt Number both O2 and CO2
         zsch_co2 = 2116.8 - 136.25 * ztc + 4.7353 * ztc2 - 0.092307 * ztc3 + 0.0007555 * ztc4
         zsch_o2  = 1920.4 - 135.6  * ztc + 5.2122 * ztc2 - 0.109390 * ztc3 + 0.0009377 * ztc4
         !  wind speed 
         zws =   sqrt(sqrt((sustr(mi,mj)*rho0)**2+(svstr(mi,mj)*rho0)**2)/1.25e-3) &
             & * sqrt(sqrt((sustr(mi,mj)*rho0)**2+(svstr(mi,mj)*rho0)**2)/1.25e-3)
         ! Compute the piston velocity for O2 and CO2
         zkgwan = 0.251 * zws
         zkgwan = zkgwan * xconv * ( 1.- fr_i(mi,mj) ) * tmask(mi,mj,1)
         ! compute gas exchange for CO2 and O2
         zkgco2(mi,mj) = zkgwan * SQRT( 660./ zsch_co2 )
         zkgo2 (mi,mj) = zkgwan * SQRT( 660./ zsch_o2 )
      END DO   ;   END DO
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, ztkel, zsal, zvapsw, zpco2atm, zxc2, zfugcoeff, &
                      zfco2, zfld, zflu, zpco2oce, oce_co2, zfld16, zflu16) &
      !$OMP SHARED(tempis, salinprac, tmask, satmco2, patm, chemc, zkgco2, &
                   zh2co3, 0.5*EPSILON(1.e0), rfact2, e3t, tr, jpoxy, Krhs, jpfer)
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ztkel = tempis(mi,mj,1) + 273.15
         zsal  = salinprac(mi,mj,1) + ( 1.- tmask(mi,mj,1) ) * 35.
         zvapsw    = EXP(24.4543 - 67.4509*(100.0/ztkel) - 4.8489*LOG(ztkel/100) - 0.000544*zsal)
         zpco2atm(mi,mj) = satmco2(mi,mj) * ( patm(mi,mj) - zvapsw )
         zxc2      = ( 1.0 - zpco2atm(mi,mj) * 1E-6 )**2
         zfugcoeff = EXP( patm(mi,mj) * (chemc(mi,mj,2) + 2.0 * zxc2 * chemc(mi,mj,3) ) &
         &           / ( 82.05736 * ztkel ))
         zfco2 = zpco2atm(mi,mj) * zfugcoeff

         ! Compute CO2 flux for the air-sea
         zfld = zfco2 * chemc(mi,mj,1) * zkgco2(mi,mj)  ! (mol/L) * (m/s)
         zflu = zh2co3(mi,mj) * zkgco2(mi,mj)                                   ! (mol/L) (m/s)
         zpco2oce(mi,mj) = zh2co3(mi,mj) / ( chemc(mi,mj,1) * zfugcoeff + 0.5*EPSILON(1.e0) )
         oce_co2(mi,mj)  = ( zfld - zflu ) * tmask(mi,mj,1) 
         ! compute the trend
         t(mi,mj,N+1-1,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-1,Krhs,itemp+ntrc_salt+jpdic) &
                 &         + oce_co2(mi,mj) * rfact2 / Hz(mi,mj,N+1-1)

         ! Compute O2 flux 
         zfld16 = patm(mi,mj) * chemo2(mi,mj,1) * zkgo2(mi,mj)          ! (mol/L)  (m/s)
         zflu16 = t(mi,mj,N+1-1,Kbb,itemp+ntrc_salt+jpoxy) * zkgo2(mi,mj)
         zoflx(mi,mj) = ( zfld16 - zflu16 ) * tmask(mi,mj,1)
         t(mi,mj,N+1-1,Krhs,itemp+ntrc_salt+jpoxy) = t(mi,mj,N+1-1,Krhs,itemp+ntrc_salt+jpoxy) &
                 &        + zoflx(mi,mj) * rfact2 / Hz(mi,mj,N+1-1)
      END DO   ;   END DO
      !$OMP END PARALLEL DO

      IF( l_dia_tcflx .OR. kt == nrst )  THEN
         ALLOCATE( zw2d(Istrp:Iendp,Jstrp:Jendp) )
         zw2d(Istrp:Iendp,Jstrp:Jendp) = oce_co2(Istrp:Iendp,Jstrp:Jendp) * e1e2t(Istrp:Iendp,Jstrp:Jendp) * 1000._wp
         t_oce_co2_flx  = glob_sum( 'p4zflx',  zw2d(:,:) )           !  Total Flux of Carbon
         t_oce_co2_flx_cum = t_oce_co2_flx_cum + t_oce_co2_flx       !  Cumulative Total Flux of Carbon
!        t_atm_co2_flx     = glob_sum( 'p4zflx', satmco2(:,:) * e1e2t(:,:) )  ! Total atmospheric pCO2
         t_atm_co2_flx     =  atcco2                                          ! Total atmospheric pCO2
         DEALLOCATE( zw2d )
      ENDIF
     
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('flx ')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
!         CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF

      IF( .false. .AND. knt == nrdttrc ) THEN
        !
        IF( l_dia_cflx ) THEN
           ALLOCATE( zw2d(-2:Lm+3+padd_X,-2:Mm+3+padd_E) )  ;  zw2d(:,:) = 0._wp
           ! Atmospheric CO2 concentration
           zw2d(Istrp:Iendp,Jstrp:Jendp) = satmco2(Istrp:Iendp,Jstrp:Jendp) * tmask(Istrp:Iendp,Jstrp:Jendp,1)
           CALL iom_put( "AtmCo2", zw2d )
           ! Carbon flux
           zw2d(Istrp:Iendp,Jstrp:Jendp) = oce_co2(Istrp:Iendp,Jstrp:Jendp) * 1000._wp
           CALL iom_put( "Cflx", zw2d )
           ! atmospheric Dpco2 
           zw2d(Istrp:Iendp,Jstrp:Jendp) =  ( zpco2atm(Istrp:Iendp,Jstrp:Jendp) - zpco2oce(Istrp:Iendp,Jstrp:Jendp) ) &
                &           * tmask(Istrp:Iendp,Jstrp:Jendp,1)
           CALL iom_put( "Dpco2", zw2d )
           ! oceanic Dpco2 
           zw2d(Istrp:Iendp,Jstrp:Jendp) =  zpco2oce(Istrp:Iendp,Jstrp:Jendp) * tmask(Istrp:Iendp,Jstrp:Jendp,1)
           CALL iom_put( "pCO2sea", zw2d )
           !
           DEALLOCATE( zw2d )
        ENDIF
        !
        IF( l_dia_oflx ) THEN
           ALLOCATE( zw2d(-2:Lm+3+padd_X,-2:Mm+3+padd_E) )  ;  zw2d(:,:) = 0._wp
           !  oxygen flux 
           zw2d(Istrp:Iendp,Jstrp:Jendp) = zoflx(Istrp:Iendp,Jstrp:Jendp) * 1000._wp
           CALL iom_put( "Oflx", zw2d )
           !  Dpo2 
           !$OMP PARALLEL DO &
           !$OMP PRIVATE(mi, mj) &
           !$OMP SHARED(atcox, patm, tr, jpoxy, Kbb, chemo2, 0.5*EPSILON(1.e0), tmask, zw2d)
           DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
              zw2d(mi,mj) =  ( atcox * patm(mi,mj) - atcox * t(mi,mj,N+1-1,Kbb,itemp+ntrc_salt+jpoxy) &
                            / ( chemo2(mi,mj,1) + 0.5*EPSILON(1.e0) ) ) * tmask(mi,mj,1) 
           END DO   ;   END DO
           !$OMP END PARALLEL DO
           CALL iom_put( "Dpo2", zw2d )
           DEALLOCATE( zw2d )
        ENDIF
        !
        IF( l_dia_kg ) THEN
           ALLOCATE( zw2d(-2:Lm+3+padd_X,-2:Mm+3+padd_E) )  ;  zw2d(:,:) = 0._wp
           zw2d(Istrp:Iendp,Jstrp:Jendp) = zkgco2(Istrp:Iendp,Jstrp:Jendp) * tmask(Istrp:Iendp,Jstrp:Jendp,1)
           CALL iom_put( "Kg", zw2d )
           DEALLOCATE( zw2d )
        ENDIF
        IF( l_dia_tcflx ) THEN
          CALL iom_put( "tcflx"   , t_oce_co2_flx )    ! global flux of carbon
          CALL iom_put( "tcflxcum", t_oce_co2_flx_cum )   !  Cumulative flux of carbon
        ENDIF
        !
      ENDIF
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj) &
      !$OMP SHARED(oce_co2, zoflx, zkgco2, tmask, zpco2atm, zpco2oce, trc2d)
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         bioVSink(mi,mj,Nfld) = oce_co2(mi,mj) * 1000.  !  carbon flux
         bioVSink(mi,mj,Nflu16 ) = zoflx(mi,mj)  * 1000.   !  O2 flux
         bioVSink(mi,mj,Nkgco2 ) = zkgco2(mi,mj) * tmask(mi,mj,1)              !  gas exchange for CO2
         bioVSink(mi,mj,Natcco2 ) = ( zpco2atm(mi,mj) - zpco2oce(mi,mj) ) * tmask(mi,mj,1) ! delta pco2
      END DO   ;   END DO
      !$OMP END PARALLEL DO
      IF( .false. )   CALL timing_stop('p4z_flx')
      !
   END SUBROUTINE p4z_flx


   SUBROUTINE p4z_flx_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_flx_init  ***
      !!
      !! ** Purpose :   Initialization of atmospheric conditions
      !!
      !! ** Method  :   Read the nampisext namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampisext
      !!----------------------------------------------------------------------
      INTEGER ::   jm, ios   ! Local integer 
      !!
      NAMELIST/nampisext/ln_co2int, atcco2, clname, nn_offset
      !!----------------------------------------------------------------------
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*)
         WRITE(stdout,*) ' p4z_flx_init : atmospheric conditions for air-sea flux calculation'
         WRITE(stdout,*) ' ~~~~~~~~~~~~'
      ENDIF
      !
      REWIND(numnatp_ref);READ(numnatp_ref,nampisext,IOSTAT=ios);CALL ctl_nam(ios,"nampisext (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,nampisext,IOSTAT=ios);CALL ctl_nam(ios,"nampisext (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE ( numonp, nampisext )
      !
      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*) '   Namelist : nampisext --- parameters for air-sea exchange'
         WRITE(stdout,*) '      reading in the atm pCO2 file or constant value   ln_co2int =', ln_co2int
      ENDIF
      !
      CALL p4z_patm( ntstart )
      !
      IF( .NOT.ln_co2int .AND. .NOT.ln_presatmco2 ) THEN
         IF(mynode .eq. 0) THEN                         ! control print
            WRITE(stdout,*) '         Constant Atmospheric pCO2 value               atcco2    =', atcco2
         ENDIF
         satmco2(:,:)  = atcco2      ! Initialisation of atmospheric pco2
      ELSEIF( ln_co2int .AND. .NOT.ln_presatmco2 ) THEN
         IF(mynode .eq. 0)  THEN
            WRITE(stdout,*) '         Constant Atmospheric pCO2 value               atcco2    =', atcco2
            WRITE(stdout,*) '         Atmospheric pCO2 value  from file             clname    =', TRIM( clname )
            WRITE(stdout,*) '         Offset model-data start year                  nn_offset =', nn_offset
         ENDIF
         CALL ctl_opn( numco2, TRIM( clname) , 'OLD', 'FORMATTED', 'SEQUENTIAL', -1 , stdout, mynode .eq. 0 )
         jm = 0                      ! Count the number of record in co2 file
         DO
           READ(numco2,*,END=100) 
           jm = jm + 1
         END DO
 100     nmaxrec = jm - 1 
         ALLOCATE( years  (nmaxrec) )   ;   years  (:) = 0._wp
         ALLOCATE( atcco2h(nmaxrec) )   ;   atcco2h(:) = 0._wp
         !
         REWIND(numco2)
         DO jm = 1, nmaxrec          ! get  xCO2 data
            READ(numco2, *)  years(jm), atcco2h(jm)
            IF(mynode .eq. 0) WRITE(stdout, '(f6.0,f7.2)')  years(jm), atcco2h(jm)
         END DO
         CLOSE(numco2)
      ELSEIF( .NOT.ln_co2int .AND. ln_presatmco2 ) THEN
         IF(mynode .eq. 0) WRITE(stdout,*) '    Spatialized Atmospheric pCO2 from an external file'
      ELSE
         IF(mynode .eq. 0) WRITE(stdout,*) '    Spatialized Atmospheric pCO2 from an external file'
      ENDIF
      !
!      oce_co2(:,:)  = 0._wp                ! Initialization of Flux of Carbon
      t_oce_co2_flx = 0._wp
      t_atm_co2_flx = 0._wp
      !
   END SUBROUTINE p4z_flx_init

   SUBROUTINE p4z_patm( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_atm  ***
      !!
      !! ** Purpose :   Read and interpolate the external atmospheric sea-level pressure
      !! ** Method  :   Read the files and interpolate the appropriate variables
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !
      INTEGER            ::   ierr, ios   ! Local integer
      CHARACTER(len=100) ::   cn_dir      ! Root directory for location of ssr files
! Need to check from RP
      INTEGER  :: mi, mj
      !!
      !!----------------------------------------------------------------------
      !
 ! Need to check from RP
         satmco2(:,:) = atcco2    ! Initialize atmco2 if no reading from a file
         patm(:,:) = 1._wp    ! Initialize patm if no reading from a file
      !
   END SUBROUTINE p4z_patm

   INTEGER FUNCTION p4z_flx_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_flx_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( satmco2(Istrp:Iendp,Jstrp:Jendp), patm(Istrp:Iendp,Jstrp:Jendp), STAT=p4z_flx_alloc )
      !
      IF( p4z_flx_alloc /= 0 )   CALL ctl_stop( 'STOP', 'p4z_flx_alloc : failed to allocate arrays' )
      !
   END FUNCTION p4z_flx_alloc

   !!======================================================================
END MODULE p4zflx
