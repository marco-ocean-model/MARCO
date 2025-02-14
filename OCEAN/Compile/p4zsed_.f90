










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









MODULE p4zsed
   !!======================================================================
   !!                         ***  MODULE p4sed  ***
   !! TOP :    Compute loss of organic matter in the sediments
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12 (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06 (C. Ethe) USE of fldread
   !!             3.5  !  2012-07 (O. Aumont) improvment of river input of nutrients 
   !!----------------------------------------------------------------------
   !!   p4z_sed        :  Compute loss of organic matter in the sediments
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !   Source Minus Sink variables
   USE p4zsink         !  Sinking fluxes
   USE sed             !  Sediment module
   USE iom             !  I/O manager
   USE prtctl          !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sed  
   PUBLIC   p4z_sed_init
   PUBLIC   p4z_sed_alloc

   REAL(wp) ::   bureffmin    !: Minimum burial efficiency
   REAL(wp) ::   bureffvar    !: Variable coef. for burial efficiency

   REAL(wp) :: sedsilfrac     !: percentage of silica loss in the sediments
   REAL(wp) :: sedcalfrac     !: percentage of calcite loss in the sediments
   REAL(wp) :: sedfactcalmin  !: Minimum value for dissolving calcite at the bottom
   REAL(wp) :: sedfactcalvar  !: Variable  value for dissolving calcite at the bottom
   
 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:  ) :: sdenit     !: Nitrate reduction in the sediments
   !
   LOGICAL  :: l_dia_sdenit, l_dia_sed

   !! * Substitutions












































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zsed.F90 15287 2021-09-24 11:11:02Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_sed( kt, knt, Kbb, Kmm, Krhs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed  ***
      !!
      !! ** Purpose :   Compute loss of organic matter in the sediments. This
      !!              is by no way a sediment model. The loss is simply 
      !!              computed to balance the inout from rivers and dust
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level indices
      INTEGER  ::  mi, mj, jk, ikt
      REAL(wp) ::  zbureff, zflx, zflx1
      REAL(wp) ::  zfact, zfactcal
      REAL(wp) ::  zo2, zno3, zpdenit, z1pdenit, zolimit
      REAL(wp) ::  zsiloss, zsiloss2, zcaloss, zdep
      REAL(wp) ::  zwstpoc, zwstpon, zwstpop
      !
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp) :: zdenit2d, zrivno3, zrivalk, zrivsil
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zw2d
      !!---------------------------------------------------------------------
      !
      IF( .false. )  CALL timing_start('p4z_sed')
      !
      IF( kt == nittrc000 )  THEN
         l_dia_sdenit = iom_use( "Sdenit" )
         l_dia_sed    = iom_use( "SedC" ) .OR. iom_use( "SedCal" ) .OR. iom_use( "SedSi" ) 
      ENDIF
      !
      zdenit2d(:,:) = 0.e0
      zrivno3 (:,:) = 0.e0
      zrivalk (:,:) = 0.e0
      zrivsil (:,:) = 0.e0

      ! Computation of the sediment denitrification proportion: The metamodel from midlleburg (2006) is being used
      ! Computation of the fraction of organic matter that is permanently buried from Dunne's model
      ! -------------------------------------------------------
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
        IF( tmask(mi,mj,1) == 1 ) THEN
           ikt = N
           zflx = sinkpocb(mi,mj) / xstep  * 1E6 
           zdep  = gdepw(mi,mj,ikt+1,Kmm) 
           !
           zflx1  = zflx * 1E3 / 1E4
           zflx1  = LOG10( MAX( 1E-3, zflx1 ) )
           zo2    = LOG10( MAX( 10. , t(mi,mj,N+1-ikt,Kbb,itemp+ntrc_salt+jpoxy) * 1E6 ) )
           zno3   = LOG10( MAX( 1.  , t(mi,mj,N+1-ikt,Kbb,itemp+ntrc_salt+jpno3) * 1E6 * rno3 ) )
           zpdenit = -2.2567 - 1.185 * zflx1 - 0.221 * zflx1 * zflx1 &
             &       - 0.3995 * zno3 * zo2 + 1.25 * zno3    &
             &       + 0.4721 * zo2 - 0.0996 * LOG10(zdep) + 0.4256 * zflx1 * zo2
           zdenit2d(mi,mj) = 10.0**zpdenit
           !
           zflx1 = ( 7.0 + zflx )
           zbureff = bureffmin + bureffvar * zflx * zflx &
              &     / ( zflx1 * zflx1 ) * MIN( zdep / 1000.00, 1.0 )
           zrivno3(mi,mj) = 1. - zbureff
        ENDIF
      END DO   ;   END DO

      ! This loss is scaled at each bottom grid cell for equilibrating the total budget of silica in the ocean.
      ! Thus, the amount of silica lost in the sediments equal the supply at the surface (dust+rivers)
      ! ------------------------------------------------------
      !
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            ikt  = N
            zdep = 1._wp / Hz(mi,mj,N+1-ikt)
            zcaloss = sinkcalb(mi,mj) * zdep
            !
            zfactcal = MAX(-0.1, MIN( excess(mi,mj,ikt), 0.2 ) )
            zfactcal = sedfactcalmin + sedfactcalvar &
                    &            * MIN( 1., (0.1 + zfactcal) / ( 0.5 - zfactcal ) )
            zrivalk(mi,mj) = sedcalfrac * zfactcal
            t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jptal) =  t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jptal) &
                    &                  + zcaloss * zrivalk(mi,mj) * 2.0
            t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpdic) =  t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpdic) &
                    &                  + zcaloss * zrivalk(mi,mj)
         END DO   ;   END DO

         IF( .NOT. ln_p2z ) THEN
            DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               ikt  = N
               zdep = 1._wp / Hz(mi,mj,N+1-ikt)
               zsiloss = sinksilb(mi,mj) * zdep
               zsiloss2 = sinksilb(mi,mj) / xstep * 365.0 * 1E3 * 1E-4 * 1E6
               zrivsil(mi,mj) = 1.0 - sedsilfrac * zsiloss2 / ( 15.0 + zsiloss2 )
               t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpsil) = t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpsil) &
                       &                 + zsiloss * zrivsil(mi,mj)
            END DO   ;   END DO
         ENDIF
      !
      ! The 0.5 factor in zpdenit is to avoid negative NO3 concentration after
      ! denitrification in the sediments. Not very clever, but simpliest option.
         IF( ln_p2z ) THEN
            DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               ikt  = N
               zwstpoc = sinkpocb(mi,mj) / Hz(mi,mj,N+1-ikt)
               zpdenit  = MIN( 0.5 * ( t(mi,mj,N+1-ikt,Kbb,itemp+ntrc_salt+jpno3) - 0.5*EPSILON(1.e0) ) &
                   &    / rdenit, zdenit2d(mi,mj) * zwstpoc * zrivno3(mi,mj) )
               z1pdenit = zwstpoc * zrivno3(mi,mj) - zpdenit
               zolimit = MIN( ( t(mi,mj,N+1-ikt,Kbb,itemp+ntrc_salt+jpoxy) - 0.5*EPSILON(1.e0) ) &
                   &     / (o2ut + o2nit), z1pdenit * ( 1.- nitrfac(mi,mj,ikt) ) )
               t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpdoc) = t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpdoc) &
                       &                  + z1pdenit - zolimit
               t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpno3) = t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpno3) &
                       &                  + zpdenit + zolimit - rdenit * zpdenit
               t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpoxy) = t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpoxy) &
                       &                 - zolimit * (o2ut + o2nit)
               t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jptal) &
                       &              - rno3 * (zolimit + (1.-rdenit) * zpdenit )
               t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpdic) &
                       &                 + zpdenit + zolimit
               sdenit(mi,mj) = rdenit * zpdenit * Hz(mi,mj,N+1-ikt)
            END DO   ;   END DO
         ELSE
            DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               ikt  = N
               zwstpoc = sinkpocb(mi,mj) / Hz(mi,mj,N+1-ikt)
               zpdenit  = MIN( 0.5 * ( t(mi,mj,N+1-ikt,Kbb,itemp+ntrc_salt+jpno3) - 0.5*EPSILON(1.e0) ) &
                  &     / rdenit, zdenit2d(mi,mj) * zwstpoc * zrivno3(mi,mj) )
               z1pdenit = zwstpoc * zrivno3(mi,mj) - zpdenit
               zolimit = MIN( ( t(mi,mj,N+1-ikt,Kbb,itemp+ntrc_salt+jpoxy) - 0.5*EPSILON(1.e0) ) &
                    &  / o2ut, z1pdenit * ( 1.- nitrfac(mi,mj,ikt) ) )
               t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpdoc) = t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpdoc) &
                       &                + z1pdenit - zolimit
               t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jppo4) = t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jppo4) &
                       &                + zpdenit + zolimit
               t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpnh4) = t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpnh4) &
                       &                + zpdenit + zolimit
               t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpno3) = t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpno3) &
                       &                 - rdenit * zpdenit
               t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpoxy) = t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpoxy) &
                       &                 - zolimit * o2ut
               t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jptal) &
                       &                  + rno3 * (zolimit + (1.+rdenit) * zpdenit )
               t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpdic) = t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpdic) &
                       &                  + zpdenit + zolimit 
               sdenit(mi,mj) = rdenit * zpdenit * Hz(mi,mj,N+1-ikt)
            END DO   ;   END DO
         ENDIF
         IF( ln_p5z ) THEN
            DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               ikt  = N
               zdep = 1._wp / Hz(mi,mj,N+1-ikt)
               zwstpoc = sinkpocb(mi,mj) * zdep
               zwstpop = sinkpopb(mi,mj) * zdep
               zwstpon = sinkponb(mi,mj) * zdep
               zpdenit  = MIN( 0.5 * ( t(mi,mj,N+1-ikt,Kbb,itemp+ntrc_salt+jpno3) - 0.5*EPSILON(1.e0) ) &
                  &       / rdenit, zdenit2d(mi,mj) * zwstpoc * zrivno3(mi,mj) )
               z1pdenit = zwstpoc * zrivno3(mi,mj) - zpdenit
               zolimit = MIN( ( t(mi,mj,N+1-ikt,Kbb,itemp+ntrc_salt+jpoxy) - 0.5*EPSILON(1.e0) ) &
                  &     / o2ut, z1pdenit * ( 1.- nitrfac(mi,mj,ikt) ) )
               t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpdon) = t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpdon) &
                       &              + ( z1pdenit - zolimit ) * zwstpon / (zwstpoc + 0.5*EPSILON(1.e0))
               t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpdop) = t(mi,mj,N+1-ikt,Krhs,itemp+ntrc_salt+jpdop) &
                       &              + ( z1pdenit - zolimit ) * zwstpop / (zwstpoc + 0.5*EPSILON(1.e0))
            END DO   ;   END DO
         ENDIF
      !
      IF( .false. .AND. knt == nrdttrc ) THEN
          zfact = 1.e+3 * rfact2r !  conversion from molC/l/kt  to molC/m3/s
          IF( l_dia_sdenit ) THEN
             ALLOCATE( zw2d(-2:Lm+3+padd_X,-2:Mm+3+padd_E) )  ;  zw2d(:,:) = 0._wp
             DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                zw2d(mi,mj) =  sdenit(mi,mj) * rno3 * zfact
             END DO   ;   END DO
             CALL iom_put( "Sdenit", zw2d )
             DEALLOCATE( zw2d )
          ENDIF        
          IF( l_dia_sed ) THEN
             ALLOCATE( zw2d(-2:Lm+3+padd_X,-2:Mm+3+padd_E) )  ;  zw2d(:,:) = 0._wp
             DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                zw2d(mi,mj) =  ( 1.0 - zrivalk(mi,mj) ) * sinkcalb(mi,mj) * zfact
             END DO   ;   END DO
             CALL iom_put( "SedCal", zw2d )
             !
             DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                zw2d(mi,mj) =  ( 1.0 - zrivno3(mi,mj) ) * sinkpocb(mi,mj) * zfact
             END DO   ;   END DO
             CALL iom_put( "SedC", zw2d )
             !
             IF( .NOT. ln_p2z ) THEN
                DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                   zw2d(mi,mj) =  ( 1.0 - zrivsil(mi,mj) ) * sinksilb(mi,mj) * zfact
               END DO   ;   END DO
               CALL iom_put( "SedSi", zw2d )
             ENDIF
             DEALLOCATE( zw2d )
          ENDIF
          !
      ENDIF

      !
      IF(sn_cfctl%l_prttrc) THEN  ! print mean trneds (USEd for debugging)
         WRITE(charout, fmt="('sed ')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
  !       CALL prt_ctl(tab4d_1=t(:,:,N+1-:,Krhs,itemp+ntrc_salt+:), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( .false. )  CALL timing_stop('p4z_sed')
      !
   END SUBROUTINE p4z_sed

   SUBROUTINE p4z_sed_init
      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_sed_init  ***
      !!
      !! ** purpose :   initialization of some parameters
      !!
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      INTEGER  :: ios                 ! Local integer output status for namelist read
      !
      !!
     NAMELIST/nampissed/bureffmin, bureffvar, &
         &               sedsilfrac, sedcalfrac, sedfactcalmin, sedfactcalvar
      
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*)
         WRITE(stdout,*) 'p4z_sed_init : initialization of sediment mobilisation '
         WRITE(stdout,*) '~~~~~~~~~~~~ '
      ENDIF
      !                            !* set file information
      REWIND(numnatp_ref);READ(numnatp_ref,nampissed,IOSTAT=ios);CALL ctl_nam(ios,"nampissed (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,nampissed,IOSTAT=ios);CALL ctl_nam(ios,"nampissed (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE ( numonp, nampissed )


      IF(mynode .eq. 0) THEN
         WRITE(stdout,*) '      Minimum burial efficiency                           bureffmin      = ', bureffmin
         WRITE(stdout,*) '      Variable coef. for burial efficiency                bureffvar      = ', bureffvar
         WRITE(stdout,*) '      percentage of silica loss in the sediments          sedsilfrac     = ', sedsilfrac
         WRITE(stdout,*) '      percentage of calcite loss in the sediments         sedcalfrac     = ', sedcalfrac
         WRITE(stdout,*) '      Minimum value for dissolving calcite at the bottom  sedfactcalmin  = ', sedfactcalmin
         WRITE(stdout,*) '      variable value for dissolving calcite at the bottom sedfactcalvar  = ', sedfactcalvar
      ENDIF
      !
      lk_sed = ln_sediment .AND. ln_sed_2way 
      !
   END SUBROUTINE p4z_sed_init


   INTEGER FUNCTION p4z_sed_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( sdenit(Istrp:Iendp,Jstrp:Jendp), STAT=p4z_sed_alloc )
      !
      IF( p4z_sed_alloc /= 0 )   CALL ctl_stop( 'STOP', 'p4z_sed_alloc: failed to allocate arrays' )
      !
   END FUNCTION p4z_sed_alloc

   
   !!======================================================================
END MODULE p4zsed
