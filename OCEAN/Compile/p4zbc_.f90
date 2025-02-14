










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









MODULE p4zbc
   !!======================================================================
   !!                         ***  MODULE p4zbc  ***
   !! TOP :    surface boundary conditions of external inputs of nutrients
   !!======================================================================
   !! History :   3.5  !  2012-07 (O. Aumont, C. Ethe) Original code
   !!             3.*  !  2025-02  (S. Maishal, R. Person)  and optimization
   !!----------------------------------------------------------------------
   !!   p4z_bc        :  Read and interpolate time-varying nutrients fluxes
   !!   p4z_bc_init   :  Initialization of p4z_bc
   !!----------------------------------------------------------------------
   USE oce_trc
   USE sms_pisces      !   Source Minus Sink variables
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_bc
   PUBLIC   p4z_bc_init   

   !! * Module variables
   LOGICAL , PUBLIC ::   ln_dust      !: boolean for dust input from the atmosphere
   LOGICAL , PUBLIC ::   ln_ndepo     !: boolean for atmospheric deposition of N
   LOGICAL , PUBLIC ::   ln_ironsed   !: boolean for Fe input from sediments
   REAL(wp), PUBLIC ::   sedfeinput   !: Coastal release of Iron
   REAL(wp), PUBLIC ::   mfrac        !: Mineral Content of the dust
   REAL(wp), PUBLIC ::   wdust        !: Sinking speed of the dust 
   REAL(wp), PUBLIC ::   lgw_rath     !: Weak ligand ratio from sed hydro sources
   LOGICAL , PUBLIC ::   ll_bc
   LOGICAL , PUBLIC ::   ll_dust      !: boolean for dust input from the atmosphere

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   dust             !: dust fields
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ironsed          !: Coastal supply of iron
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   dustmo, ferdepmo !: 2 consecutive set of dust fields 
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   po4depmo, sildepmo !: 2 consecutive set of dust fields 
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   no3depmo, nh4depmo !: 2 consecutive set of dust fields 

   REAL(wp), PUBLIC :: sedsilfrac, sedcalfrac
   REAL(wp), PUBLIC :: year2daydta

   LOGICAL  :: l_dia_iron, l_dia_dust, l_dia_ndep

   !!* Substitution












































































   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zsbc.F90 10868 2019-04-15 10:32:56Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_bc( kt, Kbb, Kmm, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_bc  ***
      !!
      !! ** purpose :   read and interpolate the external sources of nutrients
      !!
      !! ** method  :   read the files and interpolate the appropriate variables
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt              ! ocean time step
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs  ! time level index     
      !
      INTEGER :: mi, mj, jk
      INTEGER, PARAMETER :: jpmois = 12
      INTEGER  :: irec1, irec2, i15
      REAL(wp) :: zpdtan, zpdtmo, zdemi, zt
      REAL(wp) :: zxy, zjulian, zsec
      !
      REAL(wp)   :: zdust, zwdust, zfact, ztra
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zno3dep, znh4dep
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zpo4dep, zsildep
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zirondep
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zw2d
      !!---------------------------------------------------------------------
      !
      ! Compute dust at ntstart or only if there is more than 1 time record in dust file
      IF( kt == ntstart ) THEN
        l_dia_iron   = iom_use( "Ironsed" ) 
        l_dia_dust   = iom_use( "Irondep" ) .OR. iom_use( "pdust" ) &
           &            .OR. iom_use( "Sildep" ) .OR. iom_use( "Po4dep" )
        l_dia_ndep   = iom_use( "No3dep" ) .OR. iom_use( "Nh4dep" )
        IF( mynode .eq. 0 ) THEN
           WRITE(stdout,*) ' '
           WRITE(stdout,*) ' Number of days per year in file year2daydta = ', year2daydta
           WRITE(stdout,*) ' '
        ENDIF
      ENDIF
      !
      zpdtan = ( year2daydta * day2sec ) / dt
      zpdtmo = zpdtan / float( jpmois )
      zdemi  = zpdtmo / 2.
      zt     = ( float( kt ) + zdemi) / zpdtmo


      !  recherche de l'indice des enregistrements
      !  du modele dynamique encadrant le pas de temps kt.
      !  --------------------------------------------------
      irec1 = int( zt )
      irec2 = irec1 + 1
      irec1 = MOD( irec1, jpmois )
      IF ( irec1 == 0 ) irec1 = jpmois
      irec2 = MOD( irec2, jpmois )
      IF ( irec2 == 0 ) irec2 = jpmois

      zxy = zt - float(int ( zt ) )
      !
      IF( ln_dust ) THEN
         ALLOCATE( zirondep(Istrp:Iendp,Jstrp:Jendp, N) )
         !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, zxy) &
      !$OMP SHARED(dust, dustmo, irec1, irec2)
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            dust(mi,mj) = ( 1. - zxy ) * dustmo(mi,mj,irec1) + zxy &
                           * dustmo(mi,mj,irec2)
         END DO   ;   END DO
      !$OMP END PARALLEL DO
         !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, zxy) &
      !$OMP SHARED(zirondep, ferdepmo, irec1, irec2, rfact2, e3t, Kmm)
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zirondep(mi,mj,1) =  ( ( 1. - zxy ) * ferdepmo(mi,mj,irec1) &
                    &                   + zxy   * ferdepmo(mi,mj,irec2) ) &
               &                   * rfact2 / Hz(mi,mj,N+1-1) 
         END DO   ;   END DO
      !$OMP END PARALLEL DO
         !                                              ! Iron solubilization of particles in the water column
         !                                              ! dust in kg/m2/s ---> 1/55.85 to put in mol/Fe ;  wdust in m/j
         zwdust = 0.03  / ( wdust / day2sec ) / ( 250. * day2sec )
!
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zwdust, rfact, mfrac, mMass_Fe) &
      !$OMP SHARED(dust, gdept, Kmm, wdust, zirondep)
         DO jk= 2, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zirondep(mi,mj,jk) = ( dust(mi,mj) * mfrac / mMass_Fe ) * zwdust &
                   &        * rfact * EXP( -((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) /( 250. * wdust ) )
         END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
         !                                              ! Iron solubilization of particles in the water column
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(tr, zirondep, jpfer, Krhs)
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) + zirondep(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
         !
         IF( .false. .AND. l_dia_dust ) THEN
            ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )   ;   zw3d(:,:,:) = 0._wp
            ALLOCATE( zw2d(-2:Lm+3+padd_X,-2:Mm+3+padd_E) )       ;   zw2d(:,:) = 0._wp
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, N+1-jk) &
      !$OMP SHARED(zirondep, rfact2r, e3t, tmask, Kmm, zw3d)
            DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zw3d(mi,mj,N+1-jk) = zirondep(mi,mj,jk) &
                    &         * 1.e+3 * rfact2r * Hz(mi,mj,N+1-jk) * tmask(mi,mj,jk)
            END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj) &
      !$OMP SHARED(dust, wdust, day2sec, tmask, zw2d)
            CALL iom_put( "Irondep", zw3d )  ! surface downward dust depo of iron
            !
            DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zw2d(mi,mj) = dust(mi,mj) / ( wdust * day2sec ) * tmask(mi,mj,1) ! dust concentration at surface
            END DO   ;   END DO
      !$OMP END PARALLEL DO
      !
            CALL iom_put( "pdust", zw2d ) ! dust concentration at surface
            DEALLOCATE( zw2d, zw3d )
            !
         ENDIF
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(zirondep, tmask, trc3d, Nirondep)
         DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            bioFlux(mi,mj,N+1-jk,Nirondep)  = zirondep(mi,mj,jk) * 1.e+3 * tmask(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
         DEALLOCATE( zirondep )
         !                                              
         !
      ENDIF
      !
      IF( ln_ndepo ) THEN
         ALLOCATE( zno3dep(Istrp:Iendp,Jstrp:Jendp) )
         !
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, ztra) &
      !$OMP SHARED(zxy, no3depmo, rfact2, e3t, tr, rno3)
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zno3dep(mi,mj) = ( 1. - zxy ) * no3depmo(mi,mj,irec1) &
               &                  + zxy   * no3depmo(mi,mj,irec2)
            ztra  = zno3dep(mi,mj) * rfact2 / Hz(mi,mj,N+1-1)
            t(mi,mj,N+1-1,Krhs,itemp+ntrc_salt+jpno3) = t(mi,mj,N+1-1,Krhs,itemp+ntrc_salt+jpno3) + ztra
            t(mi,mj,N+1-1,Krhs,itemp+ntrc_salt+jptal) = t(mi,mj,N+1-1,Krhs,itemp+ntrc_salt+jptal) - rno3 * ztra
         END DO   ;   END DO
      !$OMP END PARALLEL DO
         !
         IF( .false. .AND. l_dia_ndep ) THEN
             ALLOCATE( zw2d(-2:Lm+3+padd_X,-2:Mm+3+padd_E) )   ;   zw2d(:,:) = 0.
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj) &
      !$OMP SHARED(zno3dep, rno3, tmask, zw2d)
             DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
                zw2d(mi,mj) = zno3dep(mi,mj) * 1.e+3 * rno3 * tmask(mi,mj,1)
             END DO   ;   END DO
      !$OMP END PARALLEL DO
             CALL iom_put( "No3dep", zw2d )
             DEALLOCATE( zw2d )
         ENDIF
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj) &
      !$OMP SHARED(zno3dep, rno3, tmask, trc2d)
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            bioVSink(mi,mj,Nno3dep )  = zno3dep(mi,mj) * 1.e+3 * rno3 * tmask(mi,mj,1)
         END DO   ;   END DO
      !$OMP END PARALLEL DO
         DEALLOCATE( zno3dep )

      ENDIF
      ! Add the external input of iron from sediment mobilization
      ! ------------------------------------------------------
      IF( ln_ironsed ) THEN
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(tr, jpfer, Krhs, ironsed, rfact2)
         DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
           t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) = t(mi,mj,N+1-jk,Krhs,itemp+ntrc_salt+jpfer) &
               &                   + ironsed(mi,mj,jk) * rfact2
         END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
         !
         IF( .false. .AND. l_dia_iron ) THEN
            ALLOCATE( zw3d(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N) )   ;   zw3d(:,:,:) = 0.
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(zw3d, ironsed, tmask)
            DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
               zw3d(mi,mj,jk) = ironsed(mi,mj,jk) * 1.e+3 * tmask(mi,mj,jk)
            END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
            CALL iom_put( "Ironsed", zw3d )  ! iron inputs from sediments
            DEALLOCATE( zw3d )
         ENDIF
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(trc3d, ironsed, tmask)
         DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            bioFlux(mi,mj,N+1-jk,Nironsed ) = ironsed(mi,mj,jk) * 1e+3 * tmask(mi,mj,jk)  ! iron from  sediment
         END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
      ENDIF

   END SUBROUTINE p4z_bc


   SUBROUTINE p4z_bc_init( Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  routine p4z_bc_init  ***
      !!
      !! ** purpose :   initialization of the external sources of nutrients
      !!
      !! ** method  :   read the files and compute the budget
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   external netcdf files
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)  :: Kmm
!     NetCDF-3.
!
! netcdf version 3 fortran interface:
!

!
! external netcdf data types:
!
      integer nf_byte
      integer nf_int1
      integer nf_char
      integer nf_short
      integer nf_int2
      integer nf_int
      integer nf_float
      integer nf_real
      integer nf_double
      integer nf_ubyte
      integer nf_ushort
      integer nf_uint
      integer nf_int64
      integer nf_uint64

      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)

!
! default fill values:
!
      integer           nf_fill_byte
      integer           nf_fill_int1
      integer           nf_fill_char
      integer           nf_fill_short
      integer           nf_fill_int2
      integer           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double

      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690e+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690d+36)

!
! mode flags for opening and creating a netcdf dataset:
!
      integer nf_nowrite
      integer nf_write
      integer nf_clobber
      integer nf_noclobber
      integer nf_fill
      integer nf_nofill
      integer nf_lock
      integer nf_share
      integer nf_64bit_offset
      integer nf_64bit_data
      integer nf_cdf5
      integer nf_sizehint_default
      integer nf_align_chunk
      integer nf_format_classic
      integer nf_format_64bit
      integer nf_format_64bit_offset
      integer nf_format_64bit_data
      integer nf_format_cdf5
      integer nf_diskless
      integer nf_mmap

      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_64bit_data = 32)
      parameter (nf_cdf5 = nf_64bit_data)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      parameter (nf_format_64bit_offset = nf_format_64bit)
      parameter (nf_format_64bit_data = 5)
      parameter (nf_format_cdf5 = nf_format_64bit_data)
      parameter (nf_diskless = 8)
      parameter (nf_mmap = 16)

!
! size argument for defining an unlimited dimension:
!
      integer nf_unlimited
      parameter (nf_unlimited = 0)

!
! global attribute id:
!
      integer nf_global
      parameter (nf_global = 0)

!
! implementation limits:
!
      integer nf_max_dims
      integer nf_max_attrs
      integer nf_max_vars
      integer nf_max_name
      integer nf_max_var_dims

      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)

!
! error codes:
!
      integer nf_noerr
      integer nf_ebadid
      integer nf_eexist
      integer nf_einval
      integer nf_eperm
      integer nf_enotindefine
      integer nf_eindefine
      integer nf_einvalcoords
      integer nf_emaxdims
      integer nf_enameinuse
      integer nf_enotatt
      integer nf_emaxatts
      integer nf_ebadtype
      integer nf_ebaddim
      integer nf_eunlimpos
      integer nf_emaxvars
      integer nf_enotvar
      integer nf_eglobal
      integer nf_enotnc
      integer nf_ests
      integer nf_emaxname
      integer nf_eunlimit
      integer nf_enorecvars
      integer nf_echar
      integer nf_eedge
      integer nf_estride
      integer nf_ebadname
      integer nf_erange
      integer nf_enomem
      integer nf_evarsize
      integer nf_edimsize
      integer nf_etrunc

      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
!
! error handling modes:
!
      integer  nf_fatal
      integer nf_verbose

      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)

!
! miscellaneous routines:
!
      character*80   nf_inq_libvers
      external       nf_inq_libvers

      character*80   nf_strerror
!                         (integer             ncerr)
      external       nf_strerror

      logical        nf_issyserr
!                         (integer             ncerr)
      external       nf_issyserr

!
! control routines:
!
      integer         nf_inq_base_pe
!                         (integer             ncid,
!                          integer             pe)
      external        nf_inq_base_pe

      integer         nf_set_base_pe
!                         (integer             ncid,
!                          integer             pe)
      external        nf_set_base_pe

      integer         nf_create
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             ncid)
      external        nf_create

      integer         nf__create
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             initialsz,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__create

      integer         nf__create_mp
!                         (character*(*)       path,
!                          integer             cmode,
!                          integer             initialsz,
!                          integer             basepe,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__create_mp

      integer         nf_open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             ncid)
      external        nf_open

      integer         nf__open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__open

      integer         nf__open_mp
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             basepe,
!                          integer             chunksizehint,
!                          integer             ncid)
      external        nf__open_mp

      integer         nf_set_fill
!                         (integer             ncid,
!                          integer             fillmode,
!                          integer             old_mode)
      external        nf_set_fill

      integer         nf_set_default_format
!                          (integer             format,
!                          integer             old_format)
      external        nf_set_default_format

      integer         nf_redef
!                         (integer             ncid)
      external        nf_redef

      integer         nf_enddef
!                         (integer             ncid)
      external        nf_enddef

      integer         nf__enddef
!                         (integer             ncid,
!                          integer             h_minfree,
!                          integer             v_align,
!                          integer             v_minfree,
!                          integer             r_align)
      external        nf__enddef

      integer         nf_sync
!                         (integer             ncid)
      external        nf_sync

      integer         nf_abort
!                         (integer             ncid)
      external        nf_abort

      integer         nf_close
!                         (integer             ncid)
      external        nf_close

      integer         nf_delete
!                         (character*(*)       ncid)
      external        nf_delete

!
! general inquiry routines:
!

      integer         nf_inq
!                         (integer             ncid,
!                          integer             ndims,
!                          integer             nvars,
!                          integer             ngatts,
!                          integer             unlimdimid)
      external        nf_inq

! new inquire path

      integer nf_inq_path
      external nf_inq_path

      integer         nf_inq_ndims
!                         (integer             ncid,
!                          integer             ndims)
      external        nf_inq_ndims

      integer         nf_inq_nvars
!                         (integer             ncid,
!                          integer             nvars)
      external        nf_inq_nvars

      integer         nf_inq_natts
!                         (integer             ncid,
!                          integer             ngatts)
      external        nf_inq_natts

      integer         nf_inq_unlimdim
!                         (integer             ncid,
!                          integer             unlimdimid)
      external        nf_inq_unlimdim

      integer         nf_inq_format
!                         (integer             ncid,
!                          integer             format)
      external        nf_inq_format

!
! dimension routines:
!

      integer         nf_def_dim
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             len,
!                          integer             dimid)
      external        nf_def_dim

      integer         nf_inq_dimid
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             dimid)
      external        nf_inq_dimid

      integer         nf_inq_dim
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name,
!                          integer             len)
      external        nf_inq_dim

      integer         nf_inq_dimname
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name)
      external        nf_inq_dimname

      integer         nf_inq_dimlen
!                         (integer             ncid,
!                          integer             dimid,
!                          integer             len)
      external        nf_inq_dimlen

      integer         nf_rename_dim
!                         (integer             ncid,
!                          integer             dimid,
!                          character(*)        name)
      external        nf_rename_dim

!
! general attribute routines:
!

      integer         nf_inq_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len)
      external        nf_inq_att

      integer         nf_inq_attid
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             attnum)
      external        nf_inq_attid

      integer         nf_inq_atttype
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype)
      external        nf_inq_atttype

      integer         nf_inq_attlen
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             len)
      external        nf_inq_attlen

      integer         nf_inq_attname
!                         (integer             ncid,
!                          integer             varid,
!                          integer             attnum,
!                          character(*)        name)
      external        nf_inq_attname

      integer         nf_copy_att
!                         (integer             ncid_in,
!                          integer             varid_in,
!                          character(*)        name,
!                          integer             ncid_out,
!                          integer             varid_out)
      external        nf_copy_att

      integer         nf_rename_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        curname,
!                          character(*)        newname)
      external        nf_rename_att

      integer         nf_del_att
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_del_att

!
! attribute put/get routines:
!

      integer         nf_put_att_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             len,
!                          character(*)        text)
      external        nf_put_att_text

      integer         nf_get_att_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          character(*)        text)
      external        nf_get_att_text

      integer         nf_put_att_int1
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int1_t           i1vals(1))
      external        nf_put_att_int1

      integer         nf_get_att_int1
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int1_t           i1vals(1))
      external        nf_get_att_int1

      integer         nf_put_att_int2
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int2_t           i2vals(1))
      external        nf_put_att_int2

      integer         nf_get_att_int2
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int2_t           i2vals(1))
      external        nf_get_att_int2

      integer         nf_put_att_int
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          integer             ivals(1))
      external        nf_put_att_int

      integer         nf_get_att_int
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             ivals(1))
      external        nf_get_att_int

      integer         nf_put_att_int64
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          nf_int8_t           i8vals(1))
      external        nf_put_att_int64

      integer         nf_get_att_int64
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          nf_int8_t           i8vals(1))
      external        nf_get_att_int64

      integer         nf_put_att_real
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          real                rvals(1))
      external        nf_put_att_real

      integer         nf_get_att_real
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          real                rvals(1))
      external        nf_get_att_real

      integer         nf_put_att_double
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             xtype,
!                          integer             len,
!                          double              dvals(1))
      external        nf_put_att_double

      integer         nf_get_att_double
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          double              dvals(1))
      external        nf_get_att_double

!
! general variable routines:
!

      integer         nf_def_var
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             datatype,
!                          integer             ndims,
!                          integer             dimids(1),
!                          integer             varid)
      external        nf_def_var

      integer         nf_inq_var
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name,
!                          integer             datatype,
!                          integer             ndims,
!                          integer             dimids(1),
!                          integer             natts)
      external        nf_inq_var

      integer         nf_inq_varid
!                         (integer             ncid,
!                          character(*)        name,
!                          integer             varid)
      external        nf_inq_varid

      integer         nf_inq_varname
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_inq_varname

      integer         nf_inq_vartype
!                         (integer             ncid,
!                          integer             varid,
!                          integer             xtype)
      external        nf_inq_vartype

      integer         nf_inq_varndims
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ndims)
      external        nf_inq_varndims

      integer         nf_inq_vardimid
!                         (integer             ncid,
!                          integer             varid,
!                          integer             dimids(1))
      external        nf_inq_vardimid

      integer         nf_inq_varnatts
!                         (integer             ncid,
!                          integer             varid,
!                          integer             natts)
      external        nf_inq_varnatts

      integer         nf_rename_var
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        name)
      external        nf_rename_var

      integer         nf_copy_var
!                         (integer             ncid_in,
!                          integer             varid,
!                          integer             ncid_out)
      external        nf_copy_var

!
! entire variable put/get routines:
!

      integer         nf_put_var_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        text)
      external        nf_put_var_text

      integer         nf_get_var_text
!                         (integer             ncid,
!                          integer             varid,
!                          character(*)        text)
      external        nf_get_var_text

      integer         nf_put_var_int1
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int1_t           i1vals(1))
      external        nf_put_var_int1

      integer         nf_get_var_int1
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int1_t           i1vals(1))
      external        nf_get_var_int1

      integer         nf_put_var_int2
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int2_t           i2vals(1))
      external        nf_put_var_int2

      integer         nf_get_var_int2
!                         (integer             ncid,
!                          integer             varid,
!                          nf_int2_t           i2vals(1))
      external        nf_get_var_int2

      integer         nf_put_var_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ivals(1))
      external        nf_put_var_int

      integer         nf_get_var_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             ivals(1))
      external        nf_get_var_int

      integer         nf_put_var_real
!                         (integer             ncid,
!                          integer             varid,
!                          real                rvals(1))
      external        nf_put_var_real

      integer         nf_get_var_real
!                         (integer             ncid,
!                          integer             varid,
!                          real                rvals(1))
      external        nf_get_var_real

      integer         nf_put_var_double
!                         (integer             ncid,
!                          integer             varid,
!                          doubleprecision     dvals(1))
      external        nf_put_var_double

      integer         nf_get_var_double
!                         (integer             ncid,
!                          integer             varid,
!                          doubleprecision     dvals(1))
      external        nf_get_var_double

!
! single variable put/get routines:
!

      integer         nf_put_var1_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          character*1         text)
      external        nf_put_var1_text

      integer         nf_get_var1_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          character*1         text)
      external        nf_get_var1_text

      integer         nf_put_var1_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int1_t           i1val)
      external        nf_put_var1_int1

      integer         nf_get_var1_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int1_t           i1val)
      external        nf_get_var1_int1

      integer         nf_put_var1_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int2_t           i2val)
      external        nf_put_var1_int2

      integer         nf_get_var1_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          nf_int2_t           i2val)
      external        nf_get_var1_int2

      integer         nf_put_var1_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          integer             ival)
      external        nf_put_var1_int

      integer         nf_get_var1_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          integer             ival)
      external        nf_get_var1_int

      integer         nf_put_var1_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          real                rval)
      external        nf_put_var1_real

      integer         nf_get_var1_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          real                rval)
      external        nf_get_var1_real

      integer         nf_put_var1_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          doubleprecision     dval)
      external        nf_put_var1_double

      integer         nf_get_var1_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             index(1),
!                          doubleprecision     dval)
      external        nf_get_var1_double

!
! variable array put/get routines:
!

      integer         nf_put_vara_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          character(*)        text)
      external        nf_put_vara_text

      integer         nf_get_vara_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          character(*)        text)
      external        nf_get_vara_text

      integer         nf_put_vara_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_vara_int1

      integer         nf_get_vara_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_vara_int1

      integer         nf_put_vara_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_vara_int2

      integer         nf_get_vara_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_vara_int2

      integer         nf_put_vara_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             ivals(1))
      external        nf_put_vara_int

      integer         nf_get_vara_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             ivals(1))
      external        nf_get_vara_int

      integer         nf_put_vara_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          real                rvals(1))
      external        nf_put_vara_real

      integer         nf_get_vara_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          real                rvals(1))
      external        nf_get_vara_real

      integer         nf_put_vara_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          doubleprecision     dvals(1))
      external        nf_put_vara_double

      integer         nf_get_vara_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          doubleprecision     dvals(1))
      external        nf_get_vara_double

!
! strided variable put/get routines:
!

      integer         nf_put_vars_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          character(*)        text)
      external        nf_put_vars_text

      integer         nf_get_vars_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          character(*)        text)
      external        nf_get_vars_text

      integer         nf_put_vars_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_vars_int1

      integer         nf_get_vars_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_vars_int1

      integer         nf_put_vars_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_vars_int2

      integer         nf_get_vars_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_vars_int2

      integer         nf_put_vars_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             ivals(1))
      external        nf_put_vars_int

      integer         nf_get_vars_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             ivals(1))
      external        nf_get_vars_int

      integer         nf_put_vars_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          real                rvals(1))
      external        nf_put_vars_real

      integer         nf_get_vars_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          real                rvals(1))
      external        nf_get_vars_real

      integer         nf_put_vars_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          doubleprecision     dvals(1))
      external        nf_put_vars_double

      integer         nf_get_vars_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          doubleprecision     dvals(1))
      external        nf_get_vars_double

!
! mapped variable put/get routines:
!

      integer         nf_put_varm_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          character(*)        text)
      external        nf_put_varm_text

      integer         nf_get_varm_text
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          character(*)        text)
      external        nf_get_varm_text

      integer         nf_put_varm_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int1_t           i1vals(1))
      external        nf_put_varm_int1

      integer         nf_get_varm_int1
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int1_t           i1vals(1))
      external        nf_get_varm_int1

      integer         nf_put_varm_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int2_t           i2vals(1))
      external        nf_put_varm_int2

      integer         nf_get_varm_int2
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          nf_int2_t           i2vals(1))
      external        nf_get_varm_int2

      integer         nf_put_varm_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          integer             ivals(1))
      external        nf_put_varm_int

      integer         nf_get_varm_int
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          integer             ivals(1))
      external        nf_get_varm_int

      integer         nf_put_varm_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          real                rvals(1))
      external        nf_put_varm_real

      integer         nf_get_varm_real
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          real                rvals(1))
      external        nf_get_varm_real

      integer         nf_put_varm_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          doubleprecision     dvals(1))
      external        nf_put_varm_double

      integer         nf_get_varm_double
!                         (integer             ncid,
!                          integer             varid,
!                          integer             start(1),
!                          integer             count(1),
!                          integer             stride(1),
!                          integer             imap(1),
!                          doubleprecision     dvals(1))
      external        nf_get_varm_double

!     64-bit int functions.
      integer nf_put_var1_int64
      external nf_put_var1_int64
      integer nf_put_vara_int64
      external nf_put_vara_int64
      integer nf_put_vars_int64
      external nf_put_vars_int64
      integer nf_put_varm_int64
      external nf_put_varm_int64
      integer nf_put_var_int64
      external nf_put_var_int64
      integer nf_get_var1_int64
      external nf_get_var1_int64
      integer nf_get_vara_int64
      external nf_get_vara_int64
      integer nf_get_vars_int64
      external nf_get_vars_int64
      integer nf_get_varm_int64
      external nf_get_varm_int64
      integer nf_get_var_int64
      external nf_get_var_int64


!     NetCDF-4.
!     This is part of netCDF-4. Copyright 2006, UCAR, See COPYRIGHT
!     file for distribution information.

!     Netcdf version 4 fortran interface.

!     $Id: netcdf4.inc,v 1.28 2010/05/25 13:53:02 ed Exp $

!     New netCDF-4 types.
      integer nf_string
      integer nf_vlen
      integer nf_opaque
      integer nf_enum
      integer nf_compound

      parameter (nf_string = 12)
      parameter (nf_vlen = 13)
      parameter (nf_opaque = 14)
      parameter (nf_enum = 15)
      parameter (nf_compound = 16)

!     New netCDF-4 fill values.
      integer           nf_fill_ubyte
      integer           nf_fill_ushort
!      real              nf_fill_uint
!      real              nf_fill_int64
!      real              nf_fill_uint64
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)

!     New constants.
      integer nf_format_netcdf4
      parameter (nf_format_netcdf4 = 3)

      integer nf_format_netcdf4_classic
      parameter (nf_format_netcdf4_classic = 4)

      integer nf_netcdf4
      parameter (nf_netcdf4 = 4096)

      integer nf_classic_model
      parameter (nf_classic_model = 256)

      integer nf_chunk_seq
      parameter (nf_chunk_seq = 0)
      integer nf_chunk_sub
      parameter (nf_chunk_sub = 1)
      integer nf_chunk_sizes
      parameter (nf_chunk_sizes = 2)

      integer nf_endian_native
      parameter (nf_endian_native = 0)
      integer nf_endian_little
      parameter (nf_endian_little = 1)
      integer nf_endian_big
      parameter (nf_endian_big = 2)

!     For NF_DEF_VAR_CHUNKING
      integer nf_chunked
      parameter (nf_chunked = 0)
      integer nf_contiguous
      parameter (nf_contiguous = 1)
      integer nf_compact
      parameter (nf_compact = 2)

!     For NF_DEF_VAR_FLETCHER32
      integer nf_nochecksum
      parameter (nf_nochecksum = 0)
      integer nf_fletcher32
      parameter (nf_fletcher32 = 1)

!     For NF_DEF_VAR_DEFLATE
      integer nf_noshuffle
      parameter (nf_noshuffle = 0)
      integer nf_shuffle
      parameter (nf_shuffle = 1)

!     For NF_DEF_VAR_SZIP
      integer nf_szip_ec_option_mask
      parameter (nf_szip_ec_option_mask = 4)
      integer nf_szip_nn_option_mask
      parameter (nf_szip_nn_option_mask = 32)

!     For parallel I/O.
      integer nf_mpiio      
      parameter (nf_mpiio = 8192)
      integer nf_mpiposix
      parameter (nf_mpiposix = 16384)
      integer nf_pnetcdf
      parameter (nf_pnetcdf = 32768)

!     For NF_VAR_PAR_ACCESS.
      integer nf_independent
      parameter (nf_independent = 0)
      integer nf_collective
      parameter (nf_collective = 1)

!     For NF_DEF_VAR_QUANTIZE.
      integer nf_noquantize
      parameter (nf_noquantize = 0)
      integer nf_quantize_bitgroom
      parameter (nf_quantize_bitgroom = 1)

!     New error codes.
      integer nf_ehdferr        ! Error at HDF5 layer. 
      parameter (nf_ehdferr = -101)
      integer nf_ecantread      ! Can't read. 
      parameter (nf_ecantread = -102)
      integer nf_ecantwrite     ! Can't write. 
      parameter (nf_ecantwrite = -103)
      integer nf_ecantcreate    ! Can't create. 
      parameter (nf_ecantcreate = -104)
      integer nf_efilemeta      ! Problem with file metadata. 
      parameter (nf_efilemeta = -105)
      integer nf_edimmeta       ! Problem with dimension metadata. 
      parameter (nf_edimmeta = -106)
      integer nf_eattmeta       ! Problem with attribute metadata. 
      parameter (nf_eattmeta = -107)
      integer nf_evarmeta       ! Problem with variable metadata. 
      parameter (nf_evarmeta = -108)
      integer nf_enocompound    ! Not a compound type. 
      parameter (nf_enocompound = -109)
      integer nf_eattexists     ! Attribute already exists. 
      parameter (nf_eattexists = -110)
      integer nf_enotnc4        ! Attempting netcdf-4 operation on netcdf-3 file.   
      parameter (nf_enotnc4 = -111)
      integer nf_estrictnc3     ! Attempting netcdf-4 operation on strict nc3 netcdf-4 file.   
      parameter (nf_estrictnc3 = -112)
      integer nf_enotnc3        ! Attempting netcdf-3 operation on netcdf-4 file.   
      parameter (nf_enotnc3 = -113)
      integer nf_enopar         ! Parallel operation on file opened for non-parallel access.   
      parameter (nf_enopar = -114)
      integer nf_eparinit       ! Error initializing for parallel access.   
      parameter (nf_eparinit = -115)
      integer nf_ebadgrpid      ! Bad group ID.   
      parameter (nf_ebadgrpid = -116)
      integer nf_ebadtypid      ! Bad type ID.   
      parameter (nf_ebadtypid = -117)
      integer nf_etypdefined    ! Type has already been defined and may not be edited. 
      parameter (nf_etypdefined = -118)
      integer nf_ebadfield      ! Bad field ID.   
      parameter (nf_ebadfield = -119)
      integer nf_ebadclass      ! Bad class.   
      parameter (nf_ebadclass = -120)
      integer nf_emaptype       ! Mapped access for atomic types only.   
      parameter (nf_emaptype = -121)
      integer nf_elatefill      ! Attempt to define fill value when data already exists. 
      parameter (nf_elatefill = -122)
      integer nf_elatedef       ! Attempt to define var properties, like deflate, after enddef. 
      parameter (nf_elatedef = -123)
      integer nf_edimscale      ! Probem with HDF5 dimscales. 
      parameter (nf_edimscale = -124)
      integer nf_enogrp       ! No group found.
      parameter (nf_enogrp = -125)
      integer nf_estorage ! Can't specify both contiguous and chunking. 
      parameter (nf_estorage = -126)    
      integer nf_ebadchunk ! Bad chunksize. 
      parameter (nf_ebadchunk = -127)    
      integer nf_enotbuilt       ! NetCDF feature not built.
      parameter (nf_enotbuilt = -128)
      integer nf_ediskless ! Error in using diskless  access. 
      parameter (nf_ediskless = -129)    
      integer nf_ecantextend ! Attempt to extend dataset during ind. I/O operation. 
      parameter (nf_ecantextend = -130)    
      integer nf_empi !  operation failed. 
      parameter (nf_empi = -131)    
      integer nf_efilter ! Filter operation failed. 
      parameter (nf_efilter = -132)    
      integer nf_ercfile ! RC file failure 
      parameter (nf_ercfile = -133)    
      integer nf_enullpad ! Header Bytes not Null-Byte padded 
      parameter (nf_enullpad = -134)    
      integer nf_einmemory ! In-memory file error 
      parameter (nf_einmemory = -135)    
      integer nf_enofilter ! Filter not defined on variable. 
      parameter (nf_enofilter = -136)    
      integer nf_enczarr ! Error at NCZarr layer. 
      parameter (nf_enczarr = -137)    
      integer nf_es3 ! Generic S3 error 
      parameter (nf_es3 = -138)    
      integer nf_eempty ! Attempt to read empty NCZarr map key 
      parameter (nf_eempty = -139)    
      integer nf_eobject ! Some object exists when it should not 
      parameter (nf_eobject = -140)    
      integer nf_enoobject ! Some object not found 
      parameter (nf_enoobject = -141)    
      integer nf_eplugin ! Unclassified failure in accessing a dynamically loaded plugin> 
      parameter (nf_eplugin = -142)    


!     New functions.

!     Parallel I/O.
      integer nf_create_par
      external nf_create_par

      integer nf_open_par
      external nf_open_par

      integer nf_var_par_access
      external nf_var_par_access

!     Functions to handle groups.
      integer nf_inq_ncid
      external nf_inq_ncid

      integer nf_inq_grps
      external nf_inq_grps

      integer nf_inq_grpname
      external nf_inq_grpname

      integer nf_inq_grpname_full
      external nf_inq_grpname_full

      integer nf_inq_grpname_len
      external nf_inq_grpname_len

      integer nf_inq_grp_parent
      external nf_inq_grp_parent

      integer nf_inq_grp_ncid
      external nf_inq_grp_ncid

      integer nf_inq_grp_full_ncid
      external nf_inq_grp_full_ncid

      integer nf_inq_varids
      external nf_inq_varids

      integer nf_inq_dimids
      external nf_inq_dimids

      integer nf_def_grp
      external nf_def_grp

!     New rename grp function

      integer nf_rename_grp
      external nf_rename_grp

!     New options for netCDF variables.
      integer nf_def_var_deflate
      external nf_def_var_deflate

      integer nf_inq_var_deflate
      external nf_inq_var_deflate

      integer nf_def_var_zstandard
      external nf_def_var_zstandard

      integer nf_inq_var_zstandard
      external nf_inq_var_zstandard

      integer nf_def_var_szip
      external nf_def_var_szip

      integer nf_inq_var_szip
      external nf_inq_var_szip

      integer nf_def_var_quantize
      external nf_def_var_quantize

      integer nf_inq_var_quantize
      external nf_inq_var_quantize

      integer nf_def_var_fletcher32
      external nf_def_var_fletcher32

      integer nf_inq_var_fletcher32
      external nf_inq_var_fletcher32

      integer nf_def_var_chunking
      external nf_def_var_chunking

      integer nf_inq_var_chunking
      external nf_inq_var_chunking

      integer nf_def_var_fill
      external nf_def_var_fill

      integer nf_inq_var_fill
      external nf_inq_var_fill

      integer nf_def_var_endian
      external nf_def_var_endian

      integer nf_inq_var_endian
      external nf_inq_var_endian

      integer nf_def_var_filter
      external nf_def_var_filter

      integer nf_inq_var_filter
      external nf_inq_var_filter

!     User defined types.
      integer nf_inq_typeids
      external nf_inq_typeids

      integer nf_inq_typeid
      external nf_inq_typeid

      integer nf_inq_type
      external nf_inq_type

      integer nf_inq_user_type
      external nf_inq_user_type

!     User defined types - compound types.
      integer nf_def_compound
      external nf_def_compound

      integer nf_insert_compound
      external nf_insert_compound

      integer nf_insert_array_compound
      external nf_insert_array_compound

      integer nf_inq_compound
      external nf_inq_compound

      integer nf_inq_compound_name
      external nf_inq_compound_name

      integer nf_inq_compound_size
      external nf_inq_compound_size

      integer nf_inq_compound_nfields
      external nf_inq_compound_nfields

      integer nf_inq_compound_field
      external nf_inq_compound_field

      integer nf_inq_compound_fieldname
      external nf_inq_compound_fieldname

      integer nf_inq_compound_fieldindex
      external nf_inq_compound_fieldindex

      integer nf_inq_compound_fieldoffset
      external nf_inq_compound_fieldoffset

      integer nf_inq_compound_fieldtype
      external nf_inq_compound_fieldtype

      integer nf_inq_compound_fieldndims
      external nf_inq_compound_fieldndims

      integer nf_inq_compound_fielddim_sizes
      external nf_inq_compound_fielddim_sizes

!     User defined types - variable length arrays.
      integer nf_def_vlen
      external nf_def_vlen

      integer nf_inq_vlen
      external nf_inq_vlen

      integer nf_free_vlen
      external nf_free_vlen

!     User defined types - enums.
      integer nf_def_enum
      external nf_def_enum

      integer nf_insert_enum
      external nf_insert_enum

      integer nf_inq_enum
      external nf_inq_enum

      integer nf_inq_enum_member
      external nf_inq_enum_member

      integer nf_inq_enum_ident
      external nf_inq_enum_ident

!     User defined types - opaque.
      integer nf_def_opaque
      external nf_def_opaque

      integer nf_inq_opaque
      external nf_inq_opaque

!     Write and read attributes of any type, including user defined
!     types.
      integer nf_put_att
      external nf_put_att
      integer nf_get_att
      external nf_get_att

!     Write and read variables of any type, including user defined
!     types.
      integer nf_put_var
      external nf_put_var
      integer nf_put_var1
      external nf_put_var1
      integer nf_put_vara
      external nf_put_vara
      integer nf_put_vars
      external nf_put_vars
      integer nf_get_var
      external nf_get_var
      integer nf_get_var1
      external nf_get_var1
      integer nf_get_vara
      external nf_get_vara
      integer nf_get_vars
      external nf_get_vars

!     For helping F77 users with VLENs.
      integer nf_get_vlen_element
      external nf_get_vlen_element
      integer nf_put_vlen_element
      external nf_put_vlen_element

!     For dealing with file level chunk cache.
      integer nf_set_chunk_cache
      external nf_set_chunk_cache
      integer nf_get_chunk_cache
      external nf_get_chunk_cache

!     For dealing with per variable chunk cache.
      integer nf_set_var_chunk_cache
      external nf_set_var_chunk_cache
      integer nf_get_var_chunk_cache
      external nf_get_var_chunk_cache

!     NetCDF-2.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin netcdf 2.4 backward compatibility:
!

!      
! functions in the fortran interface
!
      integer nccre
      integer ncopn
      integer ncddef
      integer ncdid
      integer ncvdef
      integer ncvid
      integer nctlen
      integer ncsfil

      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil


      integer ncrdwr
      integer nccreat
      integer ncexcl
      integer ncindef
      integer ncnsync
      integer nchsync
      integer ncndirty
      integer nchdirty
      integer nclink
      integer ncnowrit
      integer ncwrite
      integer ncclob
      integer ncnoclob
      integer ncglobal
      integer ncfill
      integer ncnofill
      integer maxncop
      integer maxncdim
      integer maxncatt
      integer maxncvar
      integer maxncnam
      integer maxvdims
      integer ncnoerr
      integer ncebadid
      integer ncenfile
      integer nceexist
      integer nceinval
      integer nceperm
      integer ncenotin
      integer nceindef
      integer ncecoord
      integer ncemaxds
      integer ncename
      integer ncenoatt
      integer ncemaxat
      integer ncebadty
      integer ncebadd
      integer ncests
      integer nceunlim
      integer ncemaxvs
      integer ncenotvr
      integer nceglob
      integer ncenotnc
      integer ncfoobar
      integer ncsyserr
      integer ncfatal
      integer ncverbos
      integer ncentool


!
! netcdf data types:
!
      integer ncbyte
      integer ncchar
      integer ncshort
      integer nclong
      integer ncfloat
      integer ncdouble

      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)

!     
!     masks for the struct nc flag field; passed in as 'mode' arg to
!     nccreate and ncopen.
!     

!     read/write, 0 => readonly 
      parameter(ncrdwr = 1)
!     in create phase, cleared by ncendef 
      parameter(nccreat = 2)
!     on create destroy existing file 
      parameter(ncexcl = 4)
!     in define mode, cleared by ncendef 
      parameter(ncindef = 8)
!     synchronise numrecs on change (x'10')
      parameter(ncnsync = 16)
!     synchronise whole header on change (x'20')
      parameter(nchsync = 32)
!     numrecs has changed (x'40')
      parameter(ncndirty = 64)  
!     header info has changed (x'80')
      parameter(nchdirty = 128)
!     prefill vars on endef and increase of record, the default behavior
      parameter(ncfill = 0)
!     do not fill vars on endef and increase of record (x'100')
      parameter(ncnofill = 256)
!     isa link (x'8000')
      parameter(nclink = 32768)

!     
!     'mode' arguments for nccreate and ncopen
!     
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)

!     
!     'size' argument to ncdimdef for an unlimited dimension
!     
      integer ncunlim
      parameter(ncunlim = 0)

!     
!     attribute id to put/get a global attribute
!     
      parameter(ncglobal  = 0)

!     
!     advisory maximums:
!     
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
!     not enforced 
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)

!     
!     global netcdf error status variable
!     initialized in error.c
!     

!     no error 
      parameter(ncnoerr = nf_noerr)
!     not a netcdf id 
      parameter(ncebadid = nf_ebadid)
!     too many netcdfs open 
      parameter(ncenfile = -31)   ! nc_syserr
!     netcdf file exists && ncnoclob
      parameter(nceexist = nf_eexist)
!     invalid argument 
      parameter(nceinval = nf_einval)
!     write to read only 
      parameter(nceperm = nf_eperm)
!     operation not allowed in data mode 
      parameter(ncenotin = nf_enotindefine )   
!     operation not allowed in define mode 
      parameter(nceindef = nf_eindefine)   
!     coordinates out of domain 
      parameter(ncecoord = nf_einvalcoords)
!     maxncdims exceeded 
      parameter(ncemaxds = nf_emaxdims)
!     string match to name in use 
      parameter(ncename = nf_enameinuse)   
!     attribute not found 
      parameter(ncenoatt = nf_enotatt)
!     maxncattrs exceeded 
      parameter(ncemaxat = nf_emaxatts)
!     not a netcdf data type 
      parameter(ncebadty = nf_ebadtype)
!     invalid dimension id 
      parameter(ncebadd = nf_ebaddim)
!     ncunlimited in the wrong index 
      parameter(nceunlim = nf_eunlimpos)
!     maxncvars exceeded 
      parameter(ncemaxvs = nf_emaxvars)
!     variable not found 
      parameter(ncenotvr = nf_enotvar)
!     action prohibited on ncglobal varid 
      parameter(nceglob = nf_eglobal)
!     not a netcdf file 
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname) 
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)

!     
!     global options variable. used to determine behavior of error handler.
!     initialized in lerror.c
!     
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)

!
!     default fill values.  these must be the same as in the c interface.
!
      integer filbyte
      integer filchar
      integer filshort
      integer fillong
      real filfloat
      doubleprecision fildoub

      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690e+36)
      parameter (fildoub = 9.9692099683868690e+36)

!     This is to turn on netCDF internal logging.
      integer nf_set_log_level
      external nf_set_log_level
      INTEGER  :: mi, mj, jk, irec
      INTEGER :: ncid, varid, dimid, ierr, lstr, lenstr, nf_fread, nrec_dust, nrec_ndep
      INTEGER :: vartype, nvatts, latt, nvdims
      INTEGER :: vdims(5)
      INTEGER  :: ios                 ! Local integer output status for namelist read
      CHARACTER(len=16) :: varname, dimname, attname
      REAL(wp) :: zexpide, zdenitide, zmaskt

      REAL     ::  cycle_length
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  dustmp
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::  zcmask
      !
      NAMELIST/nampisbc/ln_dust, ln_ndepo,  ln_ironsed, &
        &                sedfeinput, wdust, mfrac, lgw_rath
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*)
         WRITE(stdout,*) 'p4z_bc_init : initialization of the external sources of nutrients '
         WRITE(stdout,*) '~~~~~~~~~~~~ '
      ENDIF
      !                            !* set file information
      REWIND(numnatp_ref);READ(numnatp_ref,nampisbc,IOSTAT=ios);CALL ctl_nam(ios,"nampisbc (ref)",.TRUE.)
      REWIND(numnatp_cfg);READ(numnatp_cfg,nampisbc,IOSTAT=ios);CALL ctl_nam(ios,"nampisbc (cfg)",.FALSE.)
      IF(mynode .eq. 0) WRITE ( numonp, nampisbc )      
 
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*) '   Namelist : nampissbc '
         WRITE(stdout,*) '      dust input from the atmosphere           ln_dust     = ', ln_dust
         WRITE(stdout,*) '      atmospheric deposition of n              ln_ndepo    = ', ln_ndepo
         WRITE(stdout,*) '      Fe input from sediments                  ln_ironsed  = ', ln_ironsed
         IF( ln_ironsed ) THEN
            WRITE(stdout,*) '      coastal release of iron                  sedfeinput  = ', sedfeinput
         ENDIF
         IF( ln_ligand ) THEN
            WRITE(stdout,*) '      Weak ligand ratio from sed hydro sources  lgw_rath   = ', lgw_rath
         ENDIF
         IF( ln_dust ) THEN
            WRITE(stdout,*) '      Mineral Fe content of the dust           mfrac       = ', mfrac
            WRITE(stdout,*) '      sinking speed of the dust                wdust       = ', wdust
         ENDIF
      ENDIF

      IF( ln_dust .OR. ln_ndepo .OR. ln_ironsed ) THEN   ;   ll_bc = .TRUE.
      ELSE                                               ;   ll_bc = .FALSE.
      ENDIF

      ll_dust = ln_dust .OR. ln_sediment

      ! coastal and island masks
      ! ------------------------
      IF( ln_ironsed ) THEN     
         ALLOCATE( zcmask(Istrp:Iendp,Jstrp:Jendp,N), ironsed(Istrp:Iendp,Jstrp:Jendp,N) )
         zcmask(:,:,:) = 0.0
         DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zcmask(mi,mj,N) = 1
         END DO   ;   END DO
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk, zmaskt) &
      !$OMP SHARED(tmask_i, zcmask)
         DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            IF( rmask(mi,mj) /= 0. ) THEN
               zmaskt = rmask(mi+1,mj) * rmask(mi-1,mj  ) * rmask(mi,mj+1)    &
                 &                      * rmask(mi  ,mj-1) * rmask(mi,mj  )
               IF( zmaskt == 0. )   zcmask(mi,mj,jk ) = 0.1
            ENDIF
         END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
         !
         DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            zexpide   = MIN( 8.,( ((-1)*(z_r(mi,mj,N+1-jk)-z_w(mi,mj,N))) / 500. )**(-1.5) )
            zdenitide = -0.9543 + 0.7662 * LOG( zexpide ) - 0.235 * LOG( zexpide )**2
            zcmask(mi,mj,jk) = zcmask(mi,mj,jk) * MIN( 1., EXP( zdenitide ) / 0.5 )
         END DO   ;   END DO   ;   END DO
         ! Coastal supply of iron
         ! -------------------------
         ironsed(:,:,N) = 0.
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(mi, mj, jk) &
      !$OMP SHARED(sedfeinput, zcmask, e3t, day2sec, ironsed)
         DO jk= 1, N  ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            ironsed(mi,mj,jk) = sedfeinput * zcmask(mi,mj,jk) / ( Hz(mi,mj,N+1-jk) * day2sec )
         END DO   ;   END DO   ;   END DO
      !$OMP END PARALLEL DO
         DEALLOCATE( zcmask)
      ENDIF
      !
      !
      !    READ DUST INPUT FROM ATMOSPHERE
      !    -------------------------------------
      IF( ln_dust .OR. ln_ndepo ) THEN
         lstr = lenstr(bioname)
         ierr = nf_open (bioname(1:lstr), nf_nowrite, ncid)
         IF (ierr .eq. nf_noerr ) THEN
            IF( mynode .eq. 0) WRITE(stdout,*) ' Read atmospheric deposition in file = ', TRIM(bioname)
         ELSE
            IF( mynode .eq. 0) CALL ctl_stop( 'STOP', 'p4z_bc_ini : Needed input file for atmospheric deposition &
                & or put ln_dust and ln_ndepo to false' )
         ENDIF
         ierr = nf_inq_varid(ncid,"dust_time",varid)
! bug if compilation with gfortran
!         ierr =nf_inq_var (ncid, varid, varname, vartype, nvdims,  vdims,  nvatts) 
         ierr =nf_inq_varnatts (ncid, varid, nvatts) 
!        sm -need to check again ..... (1./365.25)
         year2daydta = year2day
         DO mi = 1, nvatts
            ierr = nf_inq_attname (ncid, varid, mi, attname)
            IF (ierr == nf_noerr) THEN
               latt = lenstr(attname)
               IF (attname(1:latt) == 'cycle_length') THEN
                  ierr = nf_get_att_double (ncid, varid, attname(1:latt), cycle_length)
                  IF (ierr == nf_noerr) THEN
                     year2daydta = cycle_length
                  ELSE
                     IF (mynode .eq. 0) write(stdout,'(/1x,4A/)') 'SET_CYCLE ERROR while ', &
                     &        'reading attribute ''', attname(1:latt), '''.'
                  ENDIF
               ENDIF
            ENDIF
         END DO
      ENDIF

      IF ( ln_dust ) THEN
         lstr = lenstr(bioname)
         ierr = nf_open (bioname(1:lstr), nf_nowrite, ncid)
         IF ( ierr .NE. nf_noerr .AND. mynode .eq. 0 ) THEN
            WRITE(stdout,4) bioname
         ENDIF
         ierr = nf_inq_varid (ncid,"dust",varid)
         IF (ierr .NE. nf_noerr .AND. mynode .eq. 0 ) THEN
            WRITE(stdout,5) "dust", bioname
         ENDIF
         ierr = nf_inq_dimid(ncid,"dust_time",dimid)
         ierr = nf_inq_dimlen(ncid,dimid,nrec_dust)
         !
         IF (mynode .eq. 0) WRITE(stdout,*) ' nrec_dust = ', nrec_dust
         ALLOCATE( dustmp(-2:Lm+3+padd_X,-2:Mm+3+padd_E,nrec_dust), dustmo(-2:Lm+3+padd_X,-2:Mm+3+padd_E,nrec_dust) )
         ALLOCATE( ferdepmo(-2:Lm+3+padd_X,-2:Mm+3+padd_E,nrec_dust), dust(Istrp:Iendp,Jstrp:Jendp) )
         !
         DO irec = 1, nrec_dust
            ierr = nf_fread(dustmp(-2,-2,irec), ncid, varid, irec, r2dvar)
            IF (ierr .NE. nf_noerr .AND. mynode .eq. 0 ) THEN
               WRITE(stdout,6) "dust", irec
            ENDIF
         END DO
         !
         IF (mynode .eq. 0) WRITE(stdout,*)
         IF (mynode .eq. 0) WRITE(stdout,'(6x,A,1x,I4)') &
         &                   'TRCINI_PISCES -- Read dust deposition ', mynode
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(irec, mj, mi) &
      !$OMP SHARED(dustmp, dustmo, nrec_dust, Mmmpi, Lmmpi)
         DO irec = 1, nrec_dust
            DO mj = 1, Mmmpi
               DO mi = 1, Lmmpi
                  dustmo(mi,mj,irec) = MAX( 0., dustmp(mi,mj,irec) )
               ENDDO
            ENDDO
         ENDDO
      !$OMP END PARALLEL DO
         !
         ierr = nf_inq_varid (ncid,"dustfer",varid)
         IF (ierr .NE. nf_noerr .AND. mynode .eq. 0 ) THEN
            WRITE(stdout,5) "dustfer", bioname
         ENDIF
         DO irec = 1, nrec_dust
            ierr = nf_fread(dustmp(-2,-2,irec), ncid, varid, irec, r2dvar)
            IF (ierr .NE. nf_noerr .AND. mynode .eq. 0 ) THEN
               WRITE(stdout,6) "dustfer", irec
            ENDIF
         END DO
         !
         IF (mynode .eq. 0) WRITE(stdout,*)
         IF (mynode .eq. 0) WRITE(stdout,'(6x,A,1x,I4)') &
         &                   'TRCINI_PISCES -- Read iron deposition ', mynode
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(irec, mj, mi) &
      !$OMP SHARED(dustmp, ferdepmo, nrec_dust, Mmmpi, Lmmpi, mfrac, mMass_Fe)
         DO irec = 1, nrec_dust
            DO mj = 1, Mmmpi
               DO mi = 1, Lmmpi
                  ferdepmo(mi,mj,irec) = MAX( 0., dustmp(mi,mj,irec) ) * mfrac / mMass_Fe
               ENDDO
            ENDDO
         ENDDO
      !$OMP END PARALLEL DO
         !
         ierr = nf_close(ncid)
         DEALLOCATE(dustmp)
         !
      ELSE
         ALLOCATE( dust(Istrp:Iendp,Jstrp:Jendp) )
         dust(:,:) = 0.
      ENDIF
!
!    READ N DEPOSITION FROM ATMOSPHERE (use dust_time for time)
!    -------------------------------------
!
      IF (ln_ndepo) THEN
         lstr = lenstr(bioname)
         ierr = nf_open (bioname(1:lstr), nf_nowrite, ncid)
         IF (ierr .NE. nf_noerr .AND. mynode .eq. 0) THEN
            WRITE(stdout,4) bioname
         ENDIF
         ierr = nf_inq_dimid(ncid,"ndep_time",dimid)
         ierr = nf_inq_dimlen(ncid,dimid,nrec_ndep)
         ALLOCATE( dustmp(-2:Lm+3+padd_X,-2:Mm+3+padd_E,nrec_ndep) )
         ALLOCATE( no3depmo(-2:Lm+3+padd_X,-2:Mm+3+padd_E,nrec_ndep) )
         ierr = nf_inq_varid (ncid,"ndep",varid)
         IF (ierr .NE. nf_noerr .AND. mynode .eq. 0 ) THEN
            WRITE(stdout,5) "ndep", bioname
         ENDIF
         !
         DO irec = 1, nrec_ndep
            ierr = nf_fread(dustmp(-2,-2,irec), ncid, varid, irec, r2dvar)
            IF (ierr .NE. nf_noerr .AND. mynode .eq. 0 ) THEN
               WRITE(stdout,6) "ndep", irec
            ENDIF
         END DO
         !
         IF (mynode .eq. 0) WRITE(stdout,*)
         IF (mynode .eq. 0) WRITE(stdout,'(6x,A,1x,I4)') &
         &                   'TRCINI_PISCES -- Read Nitrate deposition ', mynode
      !$OMP PARALLEL DO &
      !$OMP PRIVATE(irec, mj, mi) &
      !$OMP SHARED(dustmp, no3depmo, nrec_ndep, Mmmpi, Lmmpi, rno3, mMass_N)
         DO irec = 1, nrec_ndep
            DO mj = 1, Mmmpi
               DO mi = 1, Lmmpi
                  no3depmo(mi,mj,irec) = MAX( 0., dustmp(mi,mj,irec) ) / rno3 / mMass_N 
               END DO
            END DO
         END DO
      !$OMP END PARALLEL DO
         !
         ierr = nf_close(ncid)
         DEALLOCATE( dustmp )

      ENDIF

  4   FORMAT(/,' TRCINI_PISCES - unable to open forcing netCDF ',1x,A)
  5   FORMAT(/,' TRCINI_PISCES - unable to find forcing variable: ',A, &
     &                               /,14x,'in forcing netCDF  ',A)
  6   FORMAT(/,' TRCINI_PISCES - error while reading variable: ',A,2x, &
     &                                           ' at TIME index = ',i4)

   END SUBROUTINE p4z_bc_init



   !!======================================================================
END MODULE p4zbc
