










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









MODULE sed
   !!======================================================================
   !!                        ***  sed  ***
   !! Sediment :   set sediment global variables
   !!======================================================================
   !! History :
   !!        !  06-12  (C. Ethe)  Orignal
   !!----------------------------------------------------------------------
   USE par_sed
   USE oce_sed
   USE in_out_manager

      !! * Substitutions






























   IMPLICIT NONE
   PUBLIC

   PUBLIC sed_alloc

   !! Namelist
   REAL(wp), PUBLIC               ::  reac_sil            !: reactivity of silicate in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_clay           !: reactivity of clay in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_ligc           !: reactivity of Ligands [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_poc1           !: reactivity of pocl in  [s-1]
   REAL(wp), PUBLIC               ::  reac_poc2           !: reactivity of pocl in  [s-1]
   REAL(wp), PUBLIC               ::  reac_poc3           !: reactivity of pocl in  [s-1]
   REAL(wp), PUBLIC               ::  reac_poc4           !: reactivity of pocl in  [s-1]
   REAL(wp), PUBLIC               ::  reac_poc5           !: reactivity of pocl in  [s-1]
   REAL(wp), PUBLIC               ::  reac_poc6           !: reactivity of pocl in  [s-1]
   REAL(wp), PUBLIC               ::  reac_nh4            !: reactivity of NH4 in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_h2s            !: reactivity of ODU in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_fe2            !: reactivity of Fe2+ in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_feh2s          !: reactivity of Fe2+ in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_feso           !: reactivity of FeS with O2 in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  reac_fesp           !: precipitation of FeS  [mol.l-1.s-1]
   REAL(wp), PUBLIC               ::  reac_fesd           !: Dissolution  of FeS  [s-1]
   REAL(wp), PUBLIC               ::  reac_cal            !: reactivity of cal in  [l.mol-1.s-1]
   REAL(wp), PUBLIC               ::  adsnh4              !: adsorption coefficient of NH4
   REAL(wp), PUBLIC               ::  adsfe2              !: adsorption coefficient of Fe2
   REAL(wp), PUBLIC               ::  ratligc             !: C/L ratio in POC
   REAL(wp), PUBLIC               ::  so2ut 
   REAL(wp), PUBLIC               ::  srno3 
   REAL(wp), PUBLIC               ::  spo4r 
   REAL(wp), PUBLIC               ::  srDnit 
   REAL(wp), PUBLIC               ::  dtsed               !: sedimentation time step
   INTEGER , PUBLIC               ::  nitsed000
   INTEGER , PUBLIC               ::  nitsedend
   LOGICAL , PUBLIC               ::  lrst_sed       !: logical to control the trc restart write
   LOGICAL , PUBLIC               ::  ln_rst_sed  = .TRUE.     !: initialisation from a restart file or not
   LOGICAL , PUBLIC               ::  ln_btbz     = .FALSE.    !: Depth variation of the bioturbation coefficient
   LOGICAL , PUBLIC               ::  ln_irrig    = .FALSE.    !: iActivation of the bioirrigation
   LOGICAL , PUBLIC               ::  ln_sed_2way = .FALSE.    !: 2 way coupling with 
   INTEGER             , PUBLIC   ::  nn_rstsed      !: control of the time step ( 0 or 1 ) for pass. tr.
   CHARACTER(len = 80) , PUBLIC   ::  cn_sedrst_in   !: suffix of pass. tracer restart name (input)
   CHARACTER(len = 256), PUBLIC   ::  cn_sedrst_indir  !: restart input directory
   CHARACTER(len = 80) , PUBLIC   ::  cn_sedrst_out  !: suffix of pass. tracer restart name (output)
   CHARACTER(len = 256), PUBLIC   ::  cn_sedrst_outdir  !: restart output directory
   INTEGER, PUBLIC                ::  nrosorder  !: order of the rosenbrock method
   REAL(wp), PUBLIC               ::  rosatol   !: Tolerance for absolute error
   REAL(wp), PUBLIC               ::  rosrtol   !: Tolerance for relative error

   !
   REAL(wp), PUBLIC, DIMENSION(:,:,:),  ALLOCATABLE ::  pwcp, pwcpa, pwcpaa
   REAL(wp), PUBLIC, DIMENSION(:,:,:),  ALLOCATABLE ::  solcp, solcpa
   REAL(wp), PUBLIC, DIMENSION(:,:,:),  ALLOCATABLE ::  seddiff

   !! * Shared module variables
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  pwcp_dta   !: pore water data at given time-step
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  rainrg     !: rain of each solid component in [g/(cm**2.s)]
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  fromsed    !:
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  tosed      !:
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  burial      !:
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  rearatpom  !: 
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  apluss, aminuss  !: 
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  xirrigtrd, xirrigtrdtmp  !: 
   !
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  temp       !: temperature
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  salt       !: salinity
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  fecratio   !: Fe/C ratio in falling particles to the sediments
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  dzdep, slatit, slongit   !: total thickness of solid material rained [cm] in each cell
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  zkbot      !: total thickness of solid material rained [cm] in each cell
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  wacc       !: total thickness of solid material rained [cm] in each cell
   !
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  hipor      !: [h+] in mol/kg*densSW 
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  co3por     !: [co3--]solid sediment at initial time
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  volw3d     !:  ???
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  vols3d     !:  ???
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE ::  volc
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  dens_sol   !: Density of each solid fraction


   !! Chemistry
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  densSW 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  borats 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  calcon2
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  akbs  
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  ak1s 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  ak2s   
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  akws  
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  ak12s  
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  ak1ps 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  ak2ps  
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  ak3ps 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  ak12ps 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  ak123ps
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  aksis 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  aknh3
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  aksps 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  akh2s
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  sieqs
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  aks3s
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  akf3s
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  sulfats
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  fluorids
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  co3sat

   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  mol_wgt    !: molecular weight of solid sediment data
 
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE ::  trc_data    !: tracer data to share with sediment model
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  sedmask     !: tracer data to share with sediment model
   !! Geometry
   INTEGER , PUBLIC, SAVE                          ::  jpoce       !: Ocean points ( number/indices )
   INTEGER , PUBLIC, DIMENSION(:, : ), ALLOCATABLE ::  jarr       !: Computation of 1D array of sediments points
   INTEGER , PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  jsvode, isvode   !: Computation of 1D array of sediments points

   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  epkbot        !: ocean bottom layer thickness
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  gdepbot       !: Depth of the sediment
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  dzkbot        !: ocean bottom layer thickness in meters
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  dz            !: sediment layers thickness
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  por           !: porosity profile     
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  por1          !: 1-por 
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  xtortuosity   !: Tortuosity
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  vols          !: volume of solid cell fraction
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  db            !: bioturbation ceofficient
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  irrig        !: bioturbation ceofficient
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  radssol, rads1sol
   REAL(wp), PUBLIC, DIMENSION(:,:  ), ALLOCATABLE ::  saturco3
   REAL(wp), PUBLIC, DIMENSION(:    ), ALLOCATABLE ::  rstepros      !:  Number of iteration of rosenbrock method
   REAL(wp)  ::   dens               !: density of solid material
   !! Inputs / Outputs
   CHARACTER( len = 80 ), DIMENSION(jptrased  ) ::  sedtrcl
   CHARACTER( len = 20 ), DIMENSION(jptrased  ) ::  sedtrcd , sedtrcu
   CHARACTER( len = 80 ), DIMENSION(jpdia2dsed) ::  seddia2l 
   CHARACTER( len = 20 ), DIMENSION(jpdia2dsed) ::  seddia2d, seddia2u
   !
   INTEGER, PUBLIC ::  numsed = 27    ! units
   INTEGER , PUBLIC               ::  numrsr, numrsw   !

   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE ::  pwcp_avg   !: pore water sediment data at initial time
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE ::  solcp_avg  !: solid sediment data at initial time
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE ::  trcsed_avg
   REAL(wp), PUBLIC, DIMENSION(:,:) , ALLOCATABLE  ::  flxsed_avg

   INTEGER, DIMENSION(jptrased)   ::  rstsed
   INTEGER, DIMENSION(jpsol)      ::  rstsol
   INTEGER                        ::  rstph, rstsedstep
   INTEGER                        ::  ncidwrised, nrecsedpis_avg
   INTEGER                        ::  nwrtsedpis_avg, ntssedpis_avg
   INTEGER                        ::  nrpfsedpis_avg, sedTsteppis_avg,sedTimepis_avg, sedTime2pis_avg
   REAL(wp)                       ::  timesedpis_avg
   LOGICAL                        ::  ldefsedpis_avg

CONTAINS

   INTEGER FUNCTION sed_alloc()
      !!-------------------------------------------------------------------
      !!                    *** ROUTINE sed_alloc ***
      !!-------------------------------------------------------------------
      USE lib_mpp, ONLY: ctl_stop
      !!-------------------------------------------------------------------
      !
      ALLOCATE( trc_data(Istrp:Iendp,Jstrp:Jendp,jpdta), sedmask(Istrp:Iendp,Jstrp:Jendp),  &
         &      epkbot(Istrp:Iendp,Jstrp:Jendp), gdepbot(Istrp:Iendp,Jstrp:Jendp),            &
         &      dz(jpksed)  , por(jpksed) , por1(jpksed) , vols(jpksed),     &
         &      xtortuosity(jpksed), mol_wgt(jpsol),                 STAT=sed_alloc )

      IF( sed_alloc /= 0 )   CALL ctl_stop( 'STOP', 'sed_alloc: failed to allocate arrays' )
      !
   END FUNCTION sed_alloc


END MODULE sed
