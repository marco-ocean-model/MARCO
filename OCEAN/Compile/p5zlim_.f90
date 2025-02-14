










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









MODULE p5zlim
   !!======================================================================
   !!                         ***  MODULE p5zlim  ***
   !! TOP :   -QUOTA : Computes the various nutrient limitation terms
   !!                        of phytoplankton
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-04  (O. Aumont, C. Ethe) Limitation for iron modelled in quota 
   !!             3.6  !  2015-05  (O. Aumont)  quota
   !!----------------------------------------------------------------------
   !!   p5z_lim        :   Compute the nutrients limitation terms 
   !!   p5z_lim_init   :   Read the namelist 
   !!----------------------------------------------------------------------
   USE oce_trc         ! Shared ocean-passive tracers variables
   USE trc             ! Tracers defined
   USE p2zlim          ! Nutrient limitation
   USE p4zlim          ! Nutrient limitation 
   USE sms_pisces      !  variables
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC p5z_lim           ! called in p4zbio.F90  
   PUBLIC p5z_lim_init      ! called in trcsms_pisces.F90 
   PUBLIC p5z_lim_alloc     ! called in trcini_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  concpno3    !:  NO3 half saturation for picophyto  
   REAL(wp), PUBLIC ::  concpnh4    !:  NH4 half saturation for picophyto
   REAL(wp), PUBLIC ::  concnpo4    !:  PO4 half saturation for nanophyto
   REAL(wp), PUBLIC ::  concppo4    !:  PO4 half saturation for picophyto
   REAL(wp), PUBLIC ::  concdpo4    !:  PO4 half saturation for diatoms
   REAL(wp), PUBLIC ::  concpfer    !:  Iron half saturation for picophyto
   REAL(wp), PUBLIC ::  concbpo4    !:  PO4 half saturation for bacteria
   REAL(wp), PUBLIC ::  xsizepic    !:  Minimum size criteria for picophyto
   REAL(wp), PUBLIC ::  xsizerp     !:  Size ratio for picophytoplankton
   REAL(wp), PUBLIC ::  qfnopt      !:  optimal Fe quota for nanophyto
   REAL(wp), PUBLIC ::  qfpopt      !:  optimal Fe quota for picophyto
   REAL(wp), PUBLIC ::  qfdopt      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qnnmin      !:  minimum N  quota for nanophyto
   REAL(wp), PUBLIC ::  qnnmax      !:  maximum N quota for nanophyto
   REAL(wp), PUBLIC ::  qpnmin      !:  minimum P quota for nanophyto
   REAL(wp), PUBLIC ::  qpnmax      !:  maximum P quota for nanophyto
   REAL(wp), PUBLIC ::  qnpmin      !:  minimum N quota for nanophyto
   REAL(wp), PUBLIC ::  qnpmax      !:  maximum N quota for nanophyto
   REAL(wp), PUBLIC ::  qppmin      !:  minimum P quota for nanophyto
   REAL(wp), PUBLIC ::  qppmax      !:  maximum P quota for nanophyto
   REAL(wp), PUBLIC ::  qndmin      !:  minimum N quota for diatoms
   REAL(wp), PUBLIC ::  qndmax      !:  maximum N quota for diatoms
   REAL(wp), PUBLIC ::  qpdmin      !:  minimum P quota for diatoms
   REAL(wp), PUBLIC ::  qpdmax      !:  maximum P quota for diatoms
   REAL(wp), PUBLIC ::  qfnmax      !:  maximum Fe quota for nanophyto
   REAL(wp), PUBLIC ::  qfpmax      !:  maximum Fe quota for picophyto
   REAL(wp), PUBLIC ::  qfdmax      !:  maximum Fe quota for diatoms
   REAL(wp), PUBLIC ::  xpsinh4     !:  respiration cost of NH4 assimilation
   REAL(wp), PUBLIC ::  xpsino3     !:  respiration cost of NO3 assimilation
   REAL(wp), PUBLIC ::  xpsiuptk    !:  Mean respiration cost

   !!*  Allometric variations of the quotas
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqnnmin    !: Minimum N quota of nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqnnmax    !: Maximum N quota of nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpnmin    !: Minimum P quota of nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpnmax    !: Maximum P quota of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqnpmin    !: Minimum N quota of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqnpmax    !: Maximum N quota of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqppmin    !: Minimum P quota of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqppmax    !: Maximum P quota of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqndmin    !: Minimum N quota of diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqndmax    !: Maximum N quota of diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpdmin    !: Minimum P quota of diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpdmax    !: Maximum P quota of diatoms

   !!* Phytoplankton nutrient limitation terms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpicono3   !: Limitation of NO3 uptake by picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpiconh4   !: Limitation of NH4 uptake by picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpicopo4   !: Limitation of PO4 uptake by picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanodop   !: Limitation of DOP uptake by nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpicodop   !: Limitation of DOP uptake by picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatdop   !: Limitation of DOP uptake by diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpicofer   !: Limitation of Fe uptake by picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimpic    !: Limitation of picophyto PP by nutrients
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimpics   !: Limitation of picophyto PP by nutrients
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimphys   !: Limitation of nanophyto PP by nutrients
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdias   !: Limitation of diatoms PP by nutrients
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimpfe    !: Limitation of picophyto PP by Fe
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   fvnuptk    !: Maximum potential uptake rate of nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   fvpuptk    !: Maximum potential uptake rate of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   fvduptk    !: Maximum potential uptake rate of diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xqfuncfecp !: 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimnpn, xlimnpp, xlimnpd
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   ratchlp

   ! Coefficient for iron limitation following Flynn and Hipkin (1999)
   REAL(wp) ::  xcoef1   = 0.00167  / 55.85
   REAL(wp) ::  xcoef2   = 1.21E-5 * 14. / 55.85 / 7.625 * 0.5 * 1.5
   REAL(wp) ::  xcoef3   = 1.15E-4 * 14. / 55.85 / 7.625 * 0.5 
   REAL(wp) ::  rlogfactdp, rlogfactnp   

   LOGICAL  :: l_dia_nut_lim, l_dia_iron_lim, l_dia_fracal
   LOGICAL  :: l_dia_size_lim, l_dia_size_pro

   !! * Substitutions












































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p5zlim.F90 10070 2018-08-28 14:30:54Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p5z_lim( kt, knt, Kbb, Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_lim  ***
      !!
      !! ** Purpose :   Compute the co-limitations by the various nutrients
      !!                for the various phytoplankton species. Quota based
      !!                approach. The quota model is derived from theoretical
      !!                models proposed by Pahlow and Oschlies (2009) and 
      !!                Flynn (2001). Various adaptations from several 
      !!                publications by these authors have been also adopted. 
      !!
      !! ** Method  : Quota based approach. The quota model is derived from 
      !!              theoretical models by Pahlow and Oschlies (2009) and 
      !!              Flynn (2001). Various adaptations from several publications
      !!              by these authors have been also adopted.
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in)  :: kt, knt
      INTEGER, INTENT(in)  :: Kbb, Kmm  ! time level indices
      !
      INTEGER  ::   mi, mj, jk
      REAL(wp) ::   zlim1, zlim2, zlim3, zlim4, zno3, zferlim
      REAL(wp) ::   z1_trndia, z1_trnpic, z1_trnphy, ztem1, ztem2, zetot1
      REAL(wp) ::   zratio, zration, zratiof, znutlim, zfalim
      REAL(wp) ::   zconc1d, zconc1dnh4, zconc0n, zconc0nnh4, zconc0npo4, zconc0dpo4
      REAL(wp) ::   zconc0p, zconc0pnh4, zconc0ppo4, zconcpfe, zconcnfe, zconcdfe
      REAL(wp) ::   fanano, fananop, fananof, fadiat, fadiatp, fadiatf
      REAL(wp) ::   fapico, fapicop, fapicof, zlimpo4, zlimdop
      REAL(wp) ::   zrpho, zrass, zfuptk, ztrn, ztrp
      REAL(wp) ::   zproporteuk, zrassint
      REAL(wp) ::   zfvn, zfvp, zfvf, zsizen, zsizep, zsized, znanochl, zpicochl, zdiatchl
      REAL(wp) ::   zqfemn, zqfemp, zqfemd
      REAL(wp) ::   znutlimtot, zlimno3, zlimnh4, zlim1f, ztemp, zrphomin
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      !!---------------------------------------------------------------------
      !
      IF( .false. )   CALL timing_start('p5z_lim')

      IF( kt == nittrc000 )  THEN
         l_dia_nut_lim  = iom_use( "LNnut"   ) .OR. iom_use( "LDnut" ) .OR. iom_use( "LPnut" )
         l_dia_iron_lim = iom_use( "LNFe"    ) .OR. iom_use( "LDFe"  ) .OR. iom_use( "LPFe"  )
         l_dia_size_lim = iom_use( "SIZEN"   ) .OR. iom_use( "SIZED" ) .OR. iom_use( "SIZEP" )
         l_dia_size_pro = iom_use( "RASSN"   ) .OR. iom_use( "RASSP" ) .OR. iom_use( "RASSP" )
         l_dia_fracal   = iom_use( "xfracal" )
      ENDIF
      !
      sizena(:,:,:) = 0.0  ;  sizepa(:,:,:) = 0.0  ;  sizeda(:,:,:) = 0.0
      !
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         logsizen(mi,mj,jk) = LOG( sizen(mi,mj,jk) )
         logsizep(mi,mj,jk) = LOG( sizep(mi,mj,jk) )
         logsized(mi,mj,jk) = LOG( sized(mi,mj,jk) )
      END DO   ;   END DO   ;   END DO
      !      
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ! Computation of the Chl/C ratio of each phytoplankton group
         ! -------------------------------------------------------
         z1_trnphy   = 1. / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) + 0.5*EPSILON(1.e0) )
         z1_trnpic   = 1. / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppic) + 0.5*EPSILON(1.e0) )
         z1_trndia   = 1. / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) + 0.5*EPSILON(1.e0) )
         znanochl = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnch) * z1_trnphy
         zpicochl = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppch) * z1_trnpic
         zdiatchl = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdch) * z1_trndia

         ! Computation of a variable Ks for the different phytoplankton
         ! group as a function of their relative size. Allometry
         ! from Edwards et al. (2012)
         !------------------------------------------------

         ! diatoms
         zsized            = EXP(logsized(mi,mj,jk)*0.81)
         zconcdfe          = concdfer * zsized
         zconc1d           = concdno3 * zsized
         zconc1dnh4        = concdnh4 * zsized
         zconc0dpo4        = concdpo4 * zsized

         ! picophytoplankton
         zsizep            = EXP(logsizep(mi,mj,jk)*0.81)
         zconcpfe          = concpfer * zsizep
         zconc0p           = concpno3 * zsizep
         zconc0pnh4        = concpnh4 * zsizep
         zconc0ppo4        = concppo4 * zsizep

         ! nanophytoplankton
         zsizen            = EXP(logsizen(mi,mj,jk)*0.81)
         zconcnfe          = concnfer * zsizen
         zconc0n           = concnno3 * zsizen
         zconc0nnh4        = concnnh4 * zsizen
         zconc0npo4        = concnpo4 * zsizen

         ! Allometric scaling of the Chl/C ratio
         ratchln(mi,mj,jk) = ratchl * EXP( -0.078 * ( rlogfactnp + logsizen(mi,mj,jk) ) )
         ratchlp(mi,mj,jk) = ratchl * EXP( -0.078 * logsizep(mi,mj,jk) )
         ratchld(mi,mj,jk) = ratchl * EXP( -0.078 * ( rlogfactdp + logsized(mi,mj,jk) ) )

         ! Allometric variations of the minimum and maximum quotas
         ! From Talmy et al. (2014) and Maranon et al. (2013)
         ! -------------------------------------------------------
         xqnnmin(mi,mj,jk) = qnnmin * EXP(-0.18 * logsizen(mi,mj,jk) )
         xqnnmax(mi,mj,jk) = qnnmax
         xqndmin(mi,mj,jk) = qndmin * EXP(-0.18 * logsized(mi,mj,jk) )
         xqndmax(mi,mj,jk) = qndmax
         xqnpmin(mi,mj,jk) = qnpmin * EXP(-0.18 * logsizep(mi,mj,jk) )
         xqnpmax(mi,mj,jk) = qnpmax
         
         !
         ! Michaelis-Menten Limitation term for nutrients Small flagellates
         ! -----------------------------------------------
         ztrn    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3)
         ztrp    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdop) / 200.0

         ! Computation of the optimal allocation parameters
         ! Based on the different papers by Pahlow et al., and Smith et al.
         ! -----------------------------------------------------------------
         ! Nanophytoplankton
         znutlim = ztrn / zconc0n
         fanano  = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
         znutlim = ztrp / zconc0npo4
         fananop = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
         znutlim = biron(mi,mj,jk) / zconcnfe
         fananof = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )

         ! Picophytoplankton
         znutlim = ztrn / zconc0p
         fapico  = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
         znutlim = ztrp / zconc0npo4
         fapicop = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
         znutlim = biron(mi,mj,jk) / zconcpfe
         fapicof = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )

         ! Diatoms
         znutlim = ztrn / zconc1d
         fadiat  = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
         znutlim = ztrp / zconc0dpo4
         fadiatp = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
         znutlim = biron(mi,mj,jk) / zconcdfe
         fadiatf = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )

         !
         ! Michaelis-Menten Limitation term by nutrients of
         !  heterotrophic bacteria
         ! -------------------------------------------------------------
         zlim1    = ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) )  &
           &        / ( concbno3 + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) )
         !
         zlim2    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) + concbpo4)
         zlim3    = biron(mi,mj,jk) / ( concbfe + biron(mi,mj,jk) )
         zlim4    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) / ( xkdoc   + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdoc) )

         ! Xlimbac is used for DOC solubilization whereas xlimbacl
         ! is used for all the other bacterial-dependent term
         ! -------------------------------------------------------
         xlimbacl(mi,mj,jk) = MIN( zlim1, zlim2, zlim3 )
         xlimbac (mi,mj,jk) = xlimbacl(mi,mj,jk) * zlim4
         !
         ! Limitation of N based nutrients uptake (NO3 and NH4)
         zfalim     = (1.-fanano) / fanano
         zlimnh4    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) / ( zconc0n + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) )
         zlimno3    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) / ( zconc0n + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) )
         znutlimtot = (1. - fanano) * ztrn  / ( zfalim * zconc0n + ztrn )
         ztemp      = znutlimtot / ( zlimno3 + 5.0 * zlimnh4 + 0.5*EPSILON(1.e0) )
         xnanonh4(mi,mj,jk) = 5.0 * zlimnh4 * ztemp
         xnanono3(mi,mj,jk) = zlimno3 * ztemp
         !
         ! Limitation of P based nutrients (PO4 and DOP)
         zfalim     = (1.-fananop) / fananop
         zlimpo4    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) + zconc0npo4 )
         zlimdop    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdop) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdop) + zconc0npo4 )
         znutlimtot = (1. - fananop) * ztrp / ( zfalim * zconc0npo4 + ztrp )
         ztemp      = znutlimtot / ( zlimdop + 100.0 * zlimpo4 + 0.5*EPSILON(1.e0) )
         xnanopo4(mi,mj,jk) = 100.0 * zlimpo4 * ztemp
         xnanodop(mi,mj,jk) = zlimdop * ztemp
         !
         ! Limitation of Fe uptake
         zfalim     = (1.-fananof) / fananof
         xnanofer(mi,mj,jk) = (1. - fananof) * biron(mi,mj,jk) / ( biron(mi,mj,jk) + zfalim * zconcnfe )
         !
         ! The minimum iron quota depends on the size of PSU, respiration
         ! and the reduction of nitrate following the parameterization 
         ! proposed by Flynn and Hipkin (1999)
         zratiof   = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnfe) * z1_trnphy
         zqfemn    = xcoef1 * znanochl + xcoef2 + xcoef3 * xnanono3(mi,mj,jk)
         xqfuncfecn(mi,mj,jk) = zqfemn + qfnopt
         !
         zration   = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnph) * z1_trnphy
         zration   = MIN(xqnnmax(mi,mj,jk), MAX( xqnnmin(mi,mj,jk), zration ))
         fvnuptk(mi,mj,jk) = xpsiuptk / xpsinh4 * xqnnmin(mi,mj,jk) / (zration + 0.5*EPSILON(1.e0))  &
         &                   * MAX(0., (1. - ratchln(mi,mj,jk) * znanochl / 12. ) )
         !
         zlim1     = (zration - xqnnmin(mi,mj,jk) ) / (xqnnmax(mi,mj,jk) - xqnnmin(mi,mj,jk) )

         ! The value of the optimal quota in the formulation below
         ! has been found by solving a non linear equation
         zlim1f    = ( 1.13 - xqnnmin(mi,mj,jk) ) / (xqnnmax(mi,mj,jk) - xqnnmin(mi,mj,jk) )
         zlim3     = MAX( 0., ( zratiof - zqfemn ) / qfnopt )
         ! computation of the various limitation terms of nanophyto
         ! growth and PP
         xlimnfe (mi,mj,jk) = MIN( 1., zlim3 )
         xlimphy (mi,mj,jk) = MIN( 1., zlim1, zlim3 )
         xlimphys(mi,mj,jk) = MIN( 1., zlim1 / ( zlim1f + 0.5*EPSILON(1.e0) ), zlim3 )
         xlimnpn (mi,mj,jk) = MIN( 1., zlim1 )

         !
         ! Michaelis-Menten Limitation term for nutrients picophytoplankton
         ! ----------------------------------------------------------------
         ! Limitation of N based nutrients uptake (NO3 and NH4) 
         zfalim     = (1.-fapico) / fapico 
         zlimnh4    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) / ( zconc0p + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) )
         zlimno3    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) / ( zconc0p + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) )
         znutlimtot = (1. - fapico) * ztrn / ( zfalim * zconc0p + ztrn )
         ztemp      = znutlimtot / ( zlimno3 + 5.0 * zlimnh4 + 0.5*EPSILON(1.e0) )
         xpiconh4(mi,mj,jk) = 5.0 * zlimnh4 * ztemp
         xpicono3(mi,mj,jk) = zlimno3 * ztemp
         !
         ! Limitation of P based nutrients uptake (PO4 and DOP)
         zfalim     = (1.-fapicop) / fapicop 
         zlimpo4    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) + zconc0ppo4 )
         zlimdop    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdop) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdop) + zconc0ppo4 )
         znutlimtot = (1. - fapicop) * ztrp / ( zfalim * zconc0ppo4 + ztrp)
         ztemp      = znutlimtot / ( zlimdop + 100.0 * zlimpo4 + 0.5*EPSILON(1.e0) )
         xpicopo4(mi,mj,jk) = 100.0 * zlimpo4 * ztemp
         xpicodop(mi,mj,jk) = zlimdop * ztemp
         !
         zfalim     = (1.-fapicof) / fapicof
         xpicofer(mi,mj,jk) = (1. - fapicof) * biron(mi,mj,jk) / ( biron(mi,mj,jk) + zfalim * zconcpfe )
         !
         ! The minimum iron quota depends on the size of PSU, respiration
         ! and the reduction of nitrate following the parameterization 
         ! proposed by Flynn and Hipkin (1999)
         zratiof    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppfe) * z1_trnpic
         zqfemp     = xcoef1 * zpicochl + xcoef2 + xcoef3 * xpicono3(mi,mj,jk)
         xqfuncfecp(mi,mj,jk) = zqfemp + qfpopt
         !
         zration    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnpi) * z1_trnpic
         zration    = MIN(xqnpmax(mi,mj,jk), MAX( xqnpmin(mi,mj,jk), zration ))
         fvpuptk(mi,mj,jk) = xpsiuptk / xpsinh4 * xqnpmin(mi,mj,jk) / (zration + 0.5*EPSILON(1.e0))  &
           &                 * MAX(0., (1. - ratchlp(mi,mj,jk) * zpicochl / 12. ) ) 
         !
         zlim1      = (zration - xqnpmin(mi,mj,jk) ) / (xqnpmax(mi,mj,jk) - xqnpmin(mi,mj,jk) )

         ! The value of the optimal quota in the formulation below
         ! has been found by solving a non linear equation
         zlim1f     = (1.13 - xqnpmin(mi,mj,jk) ) / (xqnpmax(mi,mj,jk) - xqnpmin(mi,mj,jk) )
         zlim3      = MAX( 0.,( zratiof - zqfemp ) / qfpopt )

         ! computation of the various limitation terms of picophyto
         ! growth and PP
         xlimpfe (mi,mj,jk) = MIN( 1., zlim3 )
         xlimpic (mi,mj,jk) = MIN( 1., zlim1, zlim3 )
         xlimnpp (mi,mj,jk) = MIN( 1., zlim1 )
         xlimpics(mi,mj,jk) = MIN( 1., zlim1 / ( zlim1f + 0.5*EPSILON(1.e0) ), zlim3 )
         !
         !   Michaelis-Menten Limitation term for nutrients Diatoms
         !   ------------------------------------------------------
         !
         ! Limitation of N based nutrients uptake (NO3 and NH4)
         zfalim     = (1.-fadiat) / fadiat 
         zlimnh4    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) / ( zconc1d + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnh4) )
         zlimno3    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) / ( zconc1d + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpno3) )
         znutlimtot = (1.0 - fadiat) * ztrn / ( zfalim * zconc1d + ztrn )
         ztemp      = znutlimtot / ( zlimno3 + 5.0 * zlimnh4 + 0.5*EPSILON(1.e0) )
         xdiatnh4(mi,mj,jk) = 5.0 * zlimnh4 * ztemp
         xdiatno3(mi,mj,jk) = zlimno3 * ztemp
         !
         ! Limitation of P based nutrients uptake (PO4 and DOP)
         zfalim     = (1.-fadiatp) / fadiatp
         zlimpo4    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) + zconc0dpo4 )
         zlimdop    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdop) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdop) + zconc0dpo4 )
         znutlimtot = (1. - fadiatp) * ztrp / ( zfalim * zconc0dpo4 + ztrp )
         ztemp      = znutlimtot / ( zlimdop + 100.0 * zlimpo4 + 0.5*EPSILON(1.e0) )
         xdiatpo4(mi,mj,jk) = 100.0 * zlimpo4 * ztemp
         xdiatdop(mi,mj,jk) = zlimdop * ztemp
         !
         ! Limitation of Fe uptake
         zfalim     = (1.-fadiatf) / fadiatf
         xdiatfer(mi,mj,jk) = (1. - fadiatf) * biron(mi,mj,jk) / ( biron(mi,mj,jk) + zfalim * zconcdfe )
         !
         ! The minimum iron quota depends on the size of PSU, respiration
         ! and the reduction of nitrate following the parameterization 
         ! proposed by Flynn and Hipkin (1999)
         zratiof    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdfe) * z1_trndia
         zqfemd     = xcoef1 * zdiatchl + xcoef2 + xcoef3 * xdiatno3(mi,mj,jk)
         xqfuncfecd(mi,mj,jk) = zqfemd + qfdopt
         !
         zration    = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpndi) * z1_trndia
         zration    = MIN(xqndmax(mi,mj,jk), MAX( xqndmin(mi,mj,jk), zration ))
         fvduptk(mi,mj,jk) = xpsiuptk / xpsinh4 * xqndmin(mi,mj,jk) / (zration + 0.5*EPSILON(1.e0))   &
         &                   * MAX(0., (1. - ratchld(mi,mj,jk) * zdiatchl / 12. ) ) 
         !
         zlim1      = (zration - xqndmin(mi,mj,jk) ) / (xqndmax(mi,mj,jk) - xqndmin(mi,mj,jk) )
         ! The value of the optimal quota in the formulation below
         ! has been found by solving a non linear equation
         zlim1f     = (1.13 - xqndmin(mi,mj,jk) ) / (xqndmax(mi,mj,jk) - xqndmin(mi,mj,jk) )
         zlim3      = t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsil) / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpsil) + xksi(mi,mj) )
         zlim4      = MAX( 0., ( zratiof - zqfemd ) / qfdopt )
         ! computation of the various limitation terms of diatoms
         ! growth and PP
         xlimdfe(mi,mj,jk)  = MIN( 1., zlim4 )
         xlimdia(mi,mj,jk)  = MIN( 1., zlim1, zlim3, zlim4 )
         xlimdias(mi,mj,jk) = MIN( 1., zlim1 / (zlim1f + 0.5*EPSILON(1.e0) ), zlim3, zlim4 )
         xlimsi(mi,mj,jk)   = MIN( 1., zlim1, zlim4 )
         xlimnpd(mi,mj,jk)  = MIN( 1., zlim1 )
         !
      END DO   ;   END DO   ;   END DO

      !
      ! Compute the phosphorus quota values. It is based on Litchmann et al., 2004 and Daines et al, 2013.
      ! The relative contribution of three fonctional pools are computed: light harvesting apparatus, 
      ! nutrient uptake pool and assembly machinery. DNA is assumed to represent 1% of the dry mass of 
      ! phytoplankton (see Daines et al., 2013). 
      ! --------------------------------------------------------------------------------------------------
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ztrp     = MAX(1.E-6,t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppo4) + t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdop) / 200.0)
         ! N/P ratio of nanophytoplankton
         ! ------------------------------
         zfuptk   = 0.2 + 0.12 / ( 4.0 * sizen(mi,mj,jk) + 0.5*EPSILON(1.e0) )
         ! Computed from Inomura et al. (2020) using Pavlova Lutheri
         zrpho    = 11.55 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnch) &
                 &    / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) * 12. + 0.5*EPSILON(1.e0) )
         zrphomin = 11.55 * 0.0025
         zrass    = 0.62 * (0.25 + 0.75 * ( 1. - zrpho - zfuptk ) * xlimnpn(mi,mj,jk) )
         xqpnmin(mi,mj,jk) = ( 0.0078 + 0.62 * 0.25 * 0.0783 + zrphomin * 0.0089 ) * 16.
         xqpnmax(mi,mj,jk) = ( zrpho * 0.0089 + zrass * 0.0783 + 0.0078 + 0.022 ) * 16.
         xqpnmax(mi,mj,jk) = MIN( qpnmax, xqpnmax(mi,mj,jk) + 4000 * ztrp )

         ! N/P ratio of picophytoplankton
         ! ------------------------------
         zproporteuk = 0.15 + ( sizep(mi,mj,jk) - 1.0 ) * 0.96
         zfuptk   = 0.2 + 0.12 / ( 0.8 * sizep(mi,mj,jk) + 0.5*EPSILON(1.e0) )
         ! Computed from Inomura et al. (2020) using a synechococcus
         zrpho    = 13.4 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppch) &
                &    / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppic) * 12. + 0.5*EPSILON(1.e0) )
         zrphomin = 13.4 * 0.0025
         zrassint = ( 0.4 * 0.0517 * (1.0 - zproporteuk) + 0.62 * 0.0783 * zproporteuk )
         xqppmin(mi,mj,jk) = ( 0.0078 + zrassint * 0.25 + zrphomin * 0.0078 ) * 16.
         zrass    = zrassint * ( 0.25 + 0.75 * ( 1. - zrpho - zfuptk ) * xlimnpp(mi,mj,jk) )
         xqppmax(mi,mj,jk) = ( zrpho * 0.0076 + zrass + 0.0078 + 0.022 ) * 16.
         xqppmax(mi,mj,jk) = MIN( qppmax, xqppmax(mi,mj,jk) + ( 2000 * (1.0 - zproporteuk) + 4000.0 * zproporteuk) * ztrp )

         ! N/P ratio of diatoms
         ! --------------------
         zfuptk   = 0.2 + 0.12 / ( 6.0 * sized(mi,mj,jk) + 0.5*EPSILON(1.e0) )
         ! Computed from Inomura et al. (2020) using a synechococcus
         zrpho    = 8.08 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdch) &
                 & / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdia) * 12. + 0.5*EPSILON(1.e0) )
         zrphomin = 8.08 * 0.0025
         zrass    = 0.66 * ( 0.25 + 0.75 * ( 1. - zrpho - zfuptk ) * xlimnpd(mi,mj,jk) )
         xqpdmin(mi,mj,jk) = ( 0.0078 + 0.66 * 0.25 * 0.0783 + zrphomin * 0.0135 ) * 16.
         xqpdmax(mi,mj,jk) = ( zrpho * 0.0135 + zrass * 0.0783 + 0.0078 + 0.022 ) * 16.
         xqpdmax(mi,mj,jk) = MIN( qpdmax, xqpdmax(mi,mj,jk) + 5500 * ztrp )
      END DO   ;   END DO   ;   END DO

      ! Compute the fraction of nanophytoplankton that is made of calcifiers
      ! This is a purely adhoc formulation described in Aumont et al. (2015)
      ! This fraction depends on nutrient limitation, light, temperature
      ! --------------------------------------------------------------------
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ztem1  = MAX( 0., t(mi,mj,N+1-jk,nnew,itemp) + 1.8 )
         ztem2  = t(mi,mj,N+1-jk,nnew,itemp) - 10.
         zetot1 = MAX( 0., etot_ndcy(mi,mj,jk) - 1.) &
            &           / ( 4. + etot_ndcy(mi,mj,jk) ) * 30. / ( 30. + etot_ndcy(mi,mj,jk) ) 

         xfracal(mi,mj,jk) = caco3r * ztem1 / ( 0.1 + ztem1 )     &
            &                * MAX( 1., t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) / xsizephy )           &
            &                * ( 1. + EXP(-ztem2 * ztem2 / 25. ) )                    &
            &                * zetot1 * MIN( 1., 50. / ( hbl(mi,mj) + 0.5*EPSILON(1.e0) ) )
         xfracal(mi,mj,jk) = MAX( 0.02, MIN( 0.8 , xfracal(mi,mj,jk) ) )
      END DO   ;   END DO   ;   END DO
      !
      IF( .false. .AND. knt == nrdttrc ) THEN        ! save output diagnostics
        !
        IF( l_dia_fracal ) THEN   ! fraction of calcifiers
          ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp      
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = xfracal(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "xfracal",  zw3d)
          DEALLOCATE( zw3d )
        ENDIF
        !
        IF( l_dia_nut_lim ) THEN   ! Nutrient limitation term
          ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp      
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = xlimphy(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "LNnut",  zw3d)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = xlimdia(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "LDnut",  zw3d)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = xlimpic(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "LPnut",  zw3d)
          DEALLOCATE( zw3d )
        ENDIF
        !
        IF( l_dia_iron_lim ) THEN   ! Iron limitation term
          ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp      
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = xlimnfe(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "LNFe",  zw3d)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = xlimdfe(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "LDFe",  zw3d)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = xlimpfe(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "LPFe",  zw3d)
          DEALLOCATE( zw3d )
        ENDIF
        !
        IF( l_dia_size_lim ) THEN   ! Size limitation term
          ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp      
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = sizen(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "SIZEN",  zw3d)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = sized(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "SIZED",  zw3d)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zw3d(mi,mj,N+1-jk) = sizep(mi,mj,jk) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "SIZEP",  zw3d)
          DEALLOCATE( zw3d )
        ENDIF
        !
        IF( l_dia_size_pro ) THEN   ! Size of the protein machinery
          ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )  ;  zw3d(:,:,:) = 0._wp      
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zfuptk = 0.2 + 0.12 / ( 3.0 * sizen(mi,mj,jk) + 0.5*EPSILON(1.e0) )
             zrpho  = 11.55 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpnch) &
                     &  / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpphy) * 12. + 0.5*EPSILON(1.e0) )
             zw3d(mi,mj,N+1-jk) = MAX(0.62/4., ( 1. - zrpho - zfuptk ) &
                     &       * xlimnpn(mi,mj,jk) ) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "RASSN",  zw3d)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zfuptk = 0.2 + 0.12 / ( 3.0 * sizep(mi,mj,jk) + 0.5*EPSILON(1.e0) )
            zrpho  = 11.55 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppch) &
                &   / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jppic) * 12. + 0.5*EPSILON(1.e0) )
            zw3d(mi,mj,N+1-jk) = MAX(0.62/4., ( 1. - zrpho - zfuptk ) &
                &            * xlimnpp(mi,mj,jk) ) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "RASSP",  zw3d)
          DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
             zfuptk = 0.2 + 0.12 / ( 3.0 * sized(mi,mj,jk) + 0.5*EPSILON(1.e0) )
             zrpho  = 11.55 * t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpdch) &
                  &      / ( t(mi,mj,N+1-jk,Kbb,itemp+ntrc_salt+jpndi) * 12. + 0.5*EPSILON(1.e0) )
             zw3d(mi,mj,N+1-jk) = MAX(0.62/4., ( 1. - zrpho - zfuptk ) &
                  &  * xlimnpd(mi,mj,jk) ) * tmask(mi,mj,jk)
          END DO   ;   END DO   ;   END DO
          CALL iom_put( "RASSD",  zw3d)
          DEALLOCATE( zw3d )
        ENDIF
        !
      ENDIF
      !
      IF( .false. )  CALL timing_stop('p5z_lim')
      !
   END SUBROUTINE p5z_lim


   SUBROUTINE p5z_lim_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p5z_lim_init  ***
      !!
      !! ** Purpose :   Initialization of nutrient limitation parameters
      !!
      !! ** Method  :   Read the namp5zlim and nampisquota namelists and check
      !!      the parameters called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namp5zlim
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!
      NAMELIST/namp5zlim/ concnno3, concpno3, concdno3, concnnh4, concpnh4, concdnh4,  &
         &                concnfer, concpfer, concdfer, concbfe, concnpo4, concppo4,   &
         &                concdpo4, concbno3, concbnh4, concbpo4, xsizedia, xsizepic,  &
         &                xsizephy, xsizern, xsizerp, xsizerd, xksi1, xksi2, xkdoc,    &
         &                caco3r, oxymin, ratchl
         !
      NAMELIST/namp5zquota/ qnnmin, qnnmax, qpnmin, qpnmax, qnpmin, qnpmax, qppmin,      &
         &                  qppmax, qndmin, qndmax, qpdmin, qpdmax, qfnmax, qfpmax, qfdmax,  &
         &                  qfnopt, qfpopt, qfdopt
      !!----------------------------------------------------------------------
      !
      REWIND(numnatp_ref);READ(numnatp_ref,namp5zlim,IOSTAT=ios);CALL ctl_nam(ios,"namp5zlim (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,namp5zlim,IOSTAT=ios);CALL ctl_nam(ios,"namp5zlim (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE ( numonp, namp5zlim )
      !
      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*) ' '
         WRITE(stdout,*) ' Namelist parameters for nutrient limitations, namp5zlim'
         WRITE(stdout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(stdout,*) '    mean rainratio                           caco3r    = ', caco3r
         WRITE(stdout,*) '    C associated with Chlorophyll            ratchl    = ', ratchl
         WRITE(stdout,*) '    NO3 half saturation of nanophyto         concnno3  = ', concnno3
         WRITE(stdout,*) '    NO3 half saturation of picophyto         concpno3  = ', concpno3
         WRITE(stdout,*) '    NO3 half saturation of diatoms           concdno3  = ', concdno3
         WRITE(stdout,*) '    NH4 half saturation for phyto            concnnh4  = ', concnnh4
         WRITE(stdout,*) '    NH4 half saturation for pico             concpnh4  = ', concpnh4
         WRITE(stdout,*) '    NH4 half saturation for diatoms          concdnh4  = ', concdnh4
         WRITE(stdout,*) '    PO4 half saturation for phyto            concnpo4  = ', concnpo4
         WRITE(stdout,*) '    PO4 half saturation for pico             concppo4  = ', concppo4
         WRITE(stdout,*) '    PO4 half saturation for diatoms          concdpo4  = ', concdpo4
         WRITE(stdout,*) '    half saturation constant for Si uptake   xksi1     = ', xksi1
         WRITE(stdout,*) '    half saturation constant for Si/C        xksi2     = ', xksi2
         WRITE(stdout,*) '    half-sat. of DOC remineralization        xkdoc     = ', xkdoc
         WRITE(stdout,*) '    Iron half saturation for nanophyto       concnfer  = ', concnfer
         WRITE(stdout,*) '    Iron half saturation for picophyto       concpfer  = ', concpfer
         WRITE(stdout,*) '    Iron half saturation for diatoms         concdfer  = ', concdfer
         WRITE(stdout,*) '    size ratio for nanophytoplankton         xsizern   = ', xsizern
         WRITE(stdout,*) '    size ratio for picophytoplankton         xsizerp   = ', xsizerp
         WRITE(stdout,*) '    size ratio for diatoms                   xsizerd   = ', xsizerd
         WRITE(stdout,*) '    NO3 half saturation of bacteria          concbno3  = ', concbno3
         WRITE(stdout,*) '    NH4 half saturation for bacteria         concbnh4  = ', concbnh4
         WRITE(stdout,*) '    Minimum size criteria for diatoms        xsizedia  = ', xsizedia
         WRITE(stdout,*) '    Minimum size criteria for picophyto      xsizepic  = ', xsizepic
         WRITE(stdout,*) '    Minimum size criteria for nanophyto      xsizephy  = ', xsizephy
         WRITE(stdout,*) '    Fe half saturation for bacteria          concbfe   = ', concbfe
         WRITE(stdout,*) '    halk saturation constant for anoxia       oxymin   =' , oxymin
      ENDIF
      !
      REWIND(numnatp_ref);READ(numnatp_ref,namp5zquota,IOSTAT=ios);CALL ctl_nam(ios,"namp5zquota (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,namp5zquota,IOSTAT=ios);CALL ctl_nam(ios,"namp5zquota (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE ( numonp, namp5zquota )
      !
      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*) ' '
         WRITE(stdout,*) ' Namelist parameters for nutrient limitations, namp5zquota'
         WRITE(stdout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(stdout,*) '    optimal Fe quota for nano.               qfnopt    = ', qfnopt
         WRITE(stdout,*) '    optimal Fe quota for pico.               qfpopt    = ', qfpopt
         WRITE(stdout,*) '    Optimal Fe quota for diatoms             qfdopt    = ', qfdopt
         WRITE(stdout,*) '    Minimal N quota for nano                 qnnmin    = ', qnnmin
         WRITE(stdout,*) '    Maximal N quota for nano                 qnnmax    = ', qnnmax
         WRITE(stdout,*) '    Minimal P quota for nano                 qpnmin    = ', qpnmin
         WRITE(stdout,*) '    Maximal P quota for nano                 qpnmax    = ', qpnmax
         WRITE(stdout,*) '    Minimal N quota for pico                 qnpmin    = ', qnpmin
         WRITE(stdout,*) '    Maximal N quota for pico                 qnpmax    = ', qnpmax
         WRITE(stdout,*) '    Minimal P quota for pico                 qppmin    = ', qppmin
         WRITE(stdout,*) '    Maximal P quota for pico                 qppmax    = ', qppmax
         WRITE(stdout,*) '    Minimal N quota for diatoms              qndmin    = ', qndmin
         WRITE(stdout,*) '    Maximal N quota for diatoms              qndmax    = ', qndmax
         WRITE(stdout,*) '    Minimal P quota for diatoms              qpdmin    = ', qpdmin
         WRITE(stdout,*) '    Maximal P quota for diatoms              qpdmax    = ', qpdmax
         WRITE(stdout,*) '    Maximal Fe quota for nanophyto.          qfnmax    = ', qfnmax
         WRITE(stdout,*) '    Maximal Fe quota for picophyto.          qfpmax    = ', qfpmax
         WRITE(stdout,*) '    Maximal Fe quota for diatoms             qfdmax    = ', qfdmax
      ENDIF
      !
      ! Metabolic cost of nitrate and ammonium utilisation
      xpsino3  = 2.3 * rno3
      xpsinh4  = 1.8 * rno3
      xpsiuptk = 1.0 / 6.625
      !
      xksi2_3 = xksi2 * xksi2 * xksi2
      !
      rlogfactnp = LOG( 4.0 / 0.8 )
      rlogfactdp = LOG( 6.0 / 0.8 )
      !
      xfracal (:,:,N) = 0._wp
      xlimphy (:,:,N) = 0._wp    ;   xlimdia (:,:,N) = 0._wp   ;   xlimpic (:,:,N) = 0._wp
      xlimnfe (:,:,N) = 0._wp    ;   xlimdfe (:,:,N) = 0._wp   ;   xlimpfe (:,:,N) = 0._wp
      xnanono3(:,:,N) = 0._wp    ;   xdiatno3(:,:,N) = 0._wp   ;   xpicono3(:,:,N) = 0._wp
      xnanonh4(:,:,N) = 0._wp    ;   xdiatnh4(:,:,N) = 0._wp   ;   xpiconh4(:,:,N) = 0._wp
      xnanofer(:,:,N) = 0._wp    ;   xdiatfer(:,:,N) = 0._wp   ;   xpicofer(:,:,N) = 0._wp
      xnanopo4(:,:,N) = 0._wp    ;   xdiatpo4(:,:,N) = 0._wp   ;   xpicopo4(:,:,N) = 0._wp
      xlimbac (:,:,N) = 0._wp    ;   xlimbacl(:,:,N) = 0._wp
      xqfuncfecn(:,:,N) = 0._wp  ;   xqfuncfecd(:,:,N) = 0._wp ;   xqfuncfecp(:,:,N) = 0._wp
      fvnuptk (:,:,N) = 0._wp    ;   fvduptk (:,:,N) = 0._wp   ;   fvpuptk(:,:,N)  = 0._wp
      xlimphys(:,:,N) = 0._wp    ;   xlimdias(:,:,N) = 0._wp
      xlimnpp (:,:,N) = 0._wp    ;   xlimnpn (:,:,N) = 0._wp   ;   xlimnpd (:,:,N) = 0._wp
      xlimpics(:,:,N) = 0._wp   
      !
   END SUBROUTINE p5z_lim_init


   INTEGER FUNCTION p5z_lim_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p5z_lim_alloc  ***
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_stop
      INTEGER ::   ierr(3)        ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      !*  Biological arrays for phytoplankton growth
      ALLOCATE( xpicono3(Istrp:Iendp,Jstrp:Jendp,N), xpiconh4(Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xpicopo4(Istrp:Iendp,Jstrp:Jendp,N), xpicodop(Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xnanodop(Istrp:Iendp,Jstrp:Jendp,N), xdiatdop(Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xpicofer(Istrp:Iendp,Jstrp:Jendp,N), xlimpfe (Istrp:Iendp,Jstrp:Jendp,N),       &
         &      fvnuptk (Istrp:Iendp,Jstrp:Jendp,N), fvduptk (Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xlimphys(Istrp:Iendp,Jstrp:Jendp,N), xlimdias(Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xlimnpp (Istrp:Iendp,Jstrp:Jendp,N), xlimnpn (Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xlimnpd (Istrp:Iendp,Jstrp:Jendp,N),                              &
         &      xlimpics(Istrp:Iendp,Jstrp:Jendp,N), xqfuncfecp(Istrp:Iendp,Jstrp:Jendp,N),     &
         &      fvpuptk (Istrp:Iendp,Jstrp:Jendp,N), xlimpic (Istrp:Iendp,Jstrp:Jendp,N),    STAT=ierr(1) )
         !
      !*  Minimum/maximum quotas of phytoplankton
      ALLOCATE( xqnnmin (Istrp:Iendp,Jstrp:Jendp,N), xqnnmax(Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xqpnmin (Istrp:Iendp,Jstrp:Jendp,N), xqpnmax(Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xqnpmin (Istrp:Iendp,Jstrp:Jendp,N), xqnpmax(Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xqppmin (Istrp:Iendp,Jstrp:Jendp,N), xqppmax(Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xqndmin (Istrp:Iendp,Jstrp:Jendp,N), xqndmax(Istrp:Iendp,Jstrp:Jendp,N),       &
         &      xqpdmin (Istrp:Iendp,Jstrp:Jendp,N), xqpdmax(Istrp:Iendp,Jstrp:Jendp,N),     STAT=ierr(2) )
         !
         !*  Chl/C in chloroplast
      ALLOCATE( ratchlp (Istrp:Iendp,Jstrp:Jendp,N),                          STAT=ierr(3) )
        !
      p5z_lim_alloc = MAXVAL( ierr )
       !
      IF( p5z_lim_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p5z_lim_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p5z_lim_alloc


   !!======================================================================
END MODULE p5zlim
