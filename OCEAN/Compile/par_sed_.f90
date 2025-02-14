










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









MODULE par_sed
   !!======================================================================
   !!                        ***  par_sed  ***
   !! Sediment :   set sediment parameter
   !!======================================================================
   !! History :
   !!        !  06-12  (C. Ethe)  Orignal
   !!----------------------------------------------------------------------
   !! $Id: par_sed.F90 15450 2021-10-27 14:32:08Z cetlod $

   !! Domain characteristics
!   USE par_kind
!   USE par_oce , ONLY :       &
!      jpi      =>   jpi   ,  & !: first  dimension of grid --> i
!      jpj      =>   jpj   ,  & !: second dimension of grid --> j
!      jpij     =>   jpij  ,  & !: jpi x jpj
!      jp_tem   =>   jp_tem,  & !: indice of temperature
!      jp_sal   =>   jp_sal     !: indice of salintity

   USE par_pisces

   INTEGER, PARAMETER :: jpdta = 17

   ! Vertical sediment geometry
   INTEGER, PUBLIC   ::      &
      jpksed   = 11

   ! sediment tracer species
   INTEGER, PARAMETER ::    &
      jpsol = 11,           &   !: number of solid component
      jpwat = 11,           &   !: number of pore water component
      jpads = 2                 !! number adsorbed species

   
   ! pore water components       
   INTEGER, PARAMETER :: &
      jwoxy  = 1,        & !: oxygen
      jwno3  = 2,        & !: nitrate
      jwpo4  = 3,        & !: phosphate
      jwnh4  = 4,        & !: Ammonium
      jwh2s  = 5,        & !: Sulfate
      jwso4  = 6,        & !: H2S
      jwfe2  = 7,        & !: Fe2+
      jwalk  = 8,        & !: Alkalinity
      jwlgw  = 9,        & !: Alkalinity
      jwdic  = 10,       & !: DIC
      jwsil  = 11          !: Silicate

   ! solid components       
   INTEGER, PARAMETER ::  &
      jsfeo   = 1,        & !: iron hydroxides
      jsfes   = 2,        & !: FeS
      jscal   = 3,        & !: Calcite
      jsopal  = 4,        & !: Opal
      jsclay  = 5,        & !: clay
      jspoc1  = 6,        &
      jspoc2  = 7,        &
      jspoc3  = 8,        &
      jspoc4  = 9,        &
      jspoc5  = 10,       &
      jspoc6  = 11


   INTEGER, PARAMETER ::  &
      jptrased   = jpsol + jpwat , &
      jpvode     = jptrased - 14  , &
      jpdia2dsed = 25

   INTEGER, PARAMETER ::  &
      r2dsed  = 0    ,    &
      r3dsol  = 17,    &
      r3dsed  = 20

!   REAL(wp), PUBLIC  :: rtrn  = 0.5 * EPSILON( 1.e0 )    !: truncation value


END MODULE par_sed
