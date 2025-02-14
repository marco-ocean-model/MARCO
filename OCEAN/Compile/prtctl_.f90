










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









MODULE prtctl


   USE trc

   IMPLICIT NONE
   PUBLIC

   PUBLIC prt_ctl_ini, prt_ctl_info, prt_ctl   

   !! * Substitutions





























   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p5zprod.F90 15459 2021-10-29 08:19:18Z cetlod $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

      SUBROUTINE prt_ctl_ini

 !     ALLOCATE( tra_ctl(jptra) )
 !     tra_ctl(:) = 0.e0           ! Initialization to zero

   END SUBROUTINE prt_ctl_ini

   SUBROUTINE prt_ctl_info (clinfo, cdcomp )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE prt_ctl_trc_info  ***
      !!
      !! ** Purpose : - print information without any computation
      !!----------------------------------------------------------------------
      CHARACTER(len=*),           INTENT(in) ::   clinfo
      CHARACTER(len=3), OPTIONAL, INTENT(in) ::   cdcomp   ! only 'top' is accepted
      CHARACTER(len=3) :: clcomp

      !!
      IF( PRESENT(cdcomp) ) THEN   ;   clcomp = cdcomp
      ELSE                         ;   clcomp = 'top'
      ENDIF
      !
   END SUBROUTINE prt_ctl_info

   SUBROUTINE prt_ctl( ptab, pmask, clinfo )

      REAL(wp), DIMENSION(:,:,:,:), INTENT(in) :: ptab
      CHARACTER(len=*),  INTENT(in)           :: clinfo   ! information about the tab3d array
      REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL ::   pmask

      INTEGER :: mi, mj, jk, jn
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N,jptra)        :: ztab
      REAL(wp)  :: zsum  

      DO jn = 1, SIZE( ptab, 4 )
         DO jk = 1, SIZE( ptab, 3 )
            DO mj = 1, SIZE( ptab, 2 )
               DO mi = 1, SIZE( ptab, 1 )
                  ztab(mi,mj,jk,jn) = ptab(mi,mj,jk,jn) * pmask(mi,mj,jk)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DO jn = 1, jptra
         zsum   = SUM( ztab(:,:,:,jn) )
         IF( .true. ) CALL mpp_sum( zsum )      ! min over the global domain
         IF( mynode .eq. 0 ) WRITE(stdout,FMT="(3x,a10,' : ',D23.16)") TRIM(ctrcnm(jn)), zsum
      END DO

   END SUBROUTINE prt_ctl      
! Empty module
END MODULE prtctl

