










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








      program partit 
!
! Generic netCDF partitioning tool: reads netCDF files corresonding
! to the whole physical grid and prepares multiple files which hold
! data corresponding to different subdomains. These files can be
! then read in parallel by different  processes. 
!
! Usage:  partit NP_XI NP_ETA ncname1 ... ncnameN
! ------  where NP_XI  number of subdomains along XI-direction
!               NP_ETA number of subdomains along ETA-direction
!               ncname1 ... ncnameN  names of netCDF files
!
! Non-partitionable objects of netCDF files, such as scalar variables
! and attributes (both global and attributes to variables) are copied
! redundantly into the partitioned files, while partitionable array
! data is subdivided into subdomains and distributed among the
! partitioned files in such a manner that all files contain
! individual data without any overlap or redundantly stored data.
!
! The partitioning algorithm works as follows: The partitionable
! dimensions ('xi_rho', 'xi_u', 'eta_rho' and 'eta_v') are identified
! by name, then then their values are read and compared in pairs to
! detect if any of the directions have periodicity. It is assumed
! that ghost points corresponding to physical boundaries are stored
! in the file, but computational margins (including periodic margins) 
! are not. Consequently, if xi_rho and xi_u are equal to each other, 
! XI-direction is periodic, and if they differ by one, it is not.
! ETA-direction is treated similarly. Once periodicity type is
! determined, the internal number of internal points in each
! direction (i.e. excluding ghost points corresponding to physical
! boundaries) id divided by the number of subdomains in that
! direction and then physical boundary points are attached to
! subdomains which are adjacent to the boundaries. This results in
! slightly different dimension sizes of netCDF files corresponding
! to diffeent subdomains.
!
! Once all dimensions are sorted out, data corresponding to
! subdomains is extracted from the source file and copied into
! partial files.
!
      implicit none
      integer stdout, max_buff_size, maxdims, maxvars, maxnodes
      parameter (stdout=6,           max_buff_size=3000*2000*100,
     &           maxdims=40,         maxvars=128,     maxnodes=4000)
      character*120 ncname0, ncname(0:maxnodes-1), string 
      character*32 dimname(maxdims), varname(maxvars)
      logical XiPeriodic, EtaPeriodic, part_switch(maxvars),
     &                                           series(maxvars)
      real*8 buff(max_buff_size)
      integer narg,  NP_XI, NNODES,  xi_rho,  id_xi_rho,  id_xi_psi,
     &        arg,  NP_ETA,   node,  xi_u,    id_xi_u,    id_xi_v,
     &        iargc, subXI,     ii,  eta_rho, id_eta_rho, id_eta_psi,
     &        ierr, subETA,     jj,  eta_v,   id_eta_v,   id_eta_u,
     &        ndims, nvars, ngatts,  tsize,   unlimdimid, varatts,
     &        i,j,k,  lstr,   lvar,  lenstr,  rec,  size, ncid0, 
     &  ncid(0:maxnodes-1), dimid(maxdims),   dimsize(maxdims),
     &      varid(maxvars), vartype(maxvars), vardims(maxvars),
     &      start(maxdims), count(maxdims),   start1(maxdims),
     &      dimids(maxdims,maxvars),          ibuff(maxdims)
      common /partit_main/  buff,  ncid,  dimid,  dimsize, varid,
     &    vartype, vardims, start, count, start1, dimids,  ibuff
      integer LLm, MMm
      integer chunk_size_X, margin_X
      integer chunk_size_E, margin_E      
      integer istrmpi, jstrmpi
      integer iendmpi, jendmpi
      integer Lmmpi, Mmmpi
      
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
!
! Check how many arguments are given, complain about the error,
! if too few, otherwise extract NP_X and NP_E from the first two
! arguments.
!
      narg=iargc()
      if (narg .lt. 3) then
        write(stdout,'(/1x,A,1x,A/32x,A/)') 'Usage of partit',
     &    'should be:', 'partit NP_X NP_E ncname1 ... ncnameN'
        stop
      endif

      call getarg(1,string)
      lvar=lenstr(string)
      NP_XI=0
      do i=1,lvar
        j=ichar(string(i:i))-48
        if (j.ge.0 .and. j.le.9) then
          NP_XI=10*NP_XI+j
        else
          write(stdout,'(/8x,2(A,1x),2A/)') 'ERROR: illegal first',
     &     'argument', string(1:lvar), ', must be an integer number.'
          stop
        endif
      enddo

      call getarg(2,string)
      lvar=lenstr(string)
      NP_ETA=0
      do i=1,lvar
        j=ichar(string(i:i))-48
        if (j.ge.0 .and. j.le.9) then
          NP_ETA=10*NP_ETA+j
        else
          write(stdout,'(/8x,2(A,1x),2A/)') 'ERROR: illegal second',
     &      'argument', string(1:lvar), ', must be an integer number.'
          stop
        endif
      enddo
      NNODES=NP_XI*NP_ETA
      if (NNODES.gt.maxnodes) then
        write(stdout,'(/8x,A,1x,A,I3,1x,A/15x,A/)') 'ERROR:',
     &      'requested number of nodes',NNODES,'exceeds limit.',
     &      'Increase parameter maxnodes in file "partit.F".'
        stop
      endif
      write(stdout,'(/4x,2(4x,A,I3)/)') 'NP_XI =',  NP_XI,
     &                                  'NP_ETA =', NP_ETA
!
! Process netCDF files: open, determine if it is already
! a partitioned file, then make general inquiry. Complain
! about error if it already partitioned, or if number of
! variables and/or dimensions exceeds specified limits.
!
      do arg=3,narg
        call getarg(arg,ncname0)
        lstr=lenstr(ncname0)
        ierr=nf_open (ncname0(1:lstr), nf_nowrite, ncid0)
        if (ierr .ne. nf_noerr) then
          write(stdout,'(/8x,A,1x,A,1x,A,A/)') 'ERROR: Cannot',
     &                  'open netCDF file', ncname0(1:lstr),'.'
          goto 97     !--> next file
        endif
        ierr=nf_inq_att (ncid0, nf_global, 'partition', i,j)
        if (ierr .eq. nf_noerr) then
          write(stdout,'(/8x,3A/17x,A/)') 'WARNING: netCDF file ''',
     &                 ncname0(1:lstr), ''' is already partitioned',
     &                 'file. It cannot be partitioned any further.'
          goto 97     !--> next file
        endif

        write(stdout,'(8x,3A)') 'Processing netCDF file ''',
     &                               ncname0(1:lstr), '''...'

        ierr=nf_inq (ncid0, ndims, nvars, ngatts, unlimdimid)
        if (ierr .ne. nf_noerr) then
           write(stdout,'(/8x,A,1x,A/15x,A,1x,A,A/)') 'ERROR:',
     &        'Cannot determine number of dimensions, variables',
     &        'and attributes in netCDF file',ncname0(1:lstr),'.'
          goto 97     !--> next file
        elseif (ndims .gt. maxdims) then
          write(stdout,'(/8x,A,I4,1x,4A/15x,A,1x,A/)')
     &        'ERROR: number of dimensions', ndims,  'in netCDF',
     &        'file ''', ncname0(1:lstr), '''', 'exceeds limit.',
     &        'Increase parameter maxdims in file "partit.F".'
          goto 97     !--> next file
         elseif (nvars .gt. maxvars) then
          write(stdout,'(/8x,A,I4,1x,4A/15x,A,1x,A/)')
     &        'ERROR: number of variables',  nvars,  'in netCDF',
     &        'file ''', ncname0(1:lstr), '''', 'exceeds limit.',
     &        'Increase parameter maxvars in file "partit.F".'
          goto 97     !--> next file
        endif
!
! Sort out dimensions: For each dimension find and save its name and
! size. Then check whether all partitionable dimensions (identified
! by names 'xi_rho', 'xi_u', 'eta_rho' and 'eta_v')  are present and
! save their IDs and sizes. 
!
        tsize=1      ! <-- default value.
        do i=1,ndims
          ierr=nf_inq_dimname (ncid0, i, dimname(i))
          if (ierr .ne. nf_noerr) then
             write(stdout,'(/8x,A,I3/15x,3A/)')
     &           'ERROR: Cannot determine name for dimension ID =',
     &            i,  'in netCDF file ''',  ncname0(1:lstr),  '''.' 
             goto 97     !--> next file
          endif
          ierr=nf_inq_dimlen  (ncid0, i, dimsize(i))
          if (ierr .ne. nf_noerr) then
             lvar=lenstr(dimname(i))
             write(stdout,'(/8x,A,A,A/15x,3A/)')
     &            'ERROR: Cannot determine length of dimension ''',
     &             dimname(i)(1:lvar),  '''',  'in netCDF file ''',
     &                                       ncname0(1:lstr), '''.'
             goto 97     !--> next file
          endif
          if (i.eq. unlimdimid) then
            tsize=dimsize(i)
            dimsize(i)=nf_unlimited
          endif
        enddo
!
! Determine IDs and sizes of partitionable dimensions, 'xi_rho',
! 'xi_u', 'eta_rho' and 'eta_v'. Also save IDs of obsolete dimensions
! 'xi_psi', 'xi_v', 'eta_psi' and and 'eta_u'. These are used to
! readress obsolete dimensions according to the rules:
!
        id_xi_rho=0                     ! xi_psi  --> xi_u
        id_xi_u=0                       ! xi_v    --> xi_rho
        id_eta_rho=0                    ! eta_psi --> eta_v
        id_eta_v=0                      ! eta_u   --> eta_rho 
        id_xi_psi=0
        id_xi_v=0
        id_eta_psi=0
        id_eta_u=0
        do i=1,ndims
          lvar=lenstr(dimname(i))
          if (lvar.eq.6 .and. dimname(i)(1:lvar).eq.'xi_rho') then
            id_xi_rho=i
            xi_rho=dimsize(i)
          elseif (lvar.eq.4 .and. dimname(i)(1:lvar).eq.'xi_u') then
            id_xi_u=i
            xi_u=dimsize(i)
          elseif (lvar.eq.7.and.dimname(i)(1:lvar).eq.'eta_rho') then
            id_eta_rho=i
            eta_rho=dimsize(i)
          elseif (lvar.eq.5 .and. dimname(i)(1:lvar).eq.'eta_v') then
            id_eta_v=i
            eta_v=dimsize(i)
          elseif (lvar.eq.6 .and.dimname(i)(1:lvar).eq.'xi_psi') then
            id_xi_psi=i
          elseif (lvar.eq.4 .and. dimname(i)(1:lvar).eq.'xi_v') then
            id_xi_v=i
          elseif (lvar.eq.7.and.dimname(i)(1:lvar).eq.'eta_psi') then
            id_eta_psi=i
          elseif (lvar.eq.5 .and. dimname(i)(1:lvar).eq.'eta_u') then
            id_eta_u=i
          endif
c**       write(*,'(I3,1x,A,T16,I3)') i,dimname(i)(1:lvar),dimsize(i)
        enddo
        if (id_xi_rho.eq.0  .or. id_xi_u.eq.0 .or.
     &      id_eta_rho.eq.0 .or. id_eta_v.eq.0) then
          write(stdout,'(/8x,2A/15x,3A/)') 'ERROR: not all ',
     &            'partitionable dimensions are found',
     &            'in netCDF file ''', ncname0(1:lstr), '''.'
          goto 97     !--> next file 
        endif
!
! Determine subdomain dimensions. Here "subXI" and "subETA" are
! the nimbers of internal grid points within each subdomains (that
! is, excluding physical boundary points and computational margins). 
! The number of internal points in either direction for the the
! whole computational domain must be divisible by NP_XI and NP_ETA  
! respectively. If it cannot be divided, complain about the error
! and exit.
!
        LLm = xi_rho - 2
        MMm = eta_rho - 2
        if (xi_rho .eq. xi_u) then
          XiPeriodic=.true.
          subXI=xi_rho/NP_XI
          if (subXI*NP_XI .ne. xi_rho) then
            write(stdout,'(/8x,A,1x,A,I4,1x,A/15x,A,I3,1x,A/)')
     &         'ERROR: Cannot partition XI-direction:', 'xi_rho =',
     &          xi_rho,  'is',  'not divisible by NP_XI =',  NP_XI,
     &         'in XI-periodic case.'

            goto 97
          endif
        elseif (xi_rho .eq. xi_u+1) then 
          XiPeriodic=.false.
          subXI=(xi_rho-2)/NP_XI
          chunk_size_X=(LLm+NP_XI-1)/NP_XI
          margin_X=(NP_XI*chunk_size_X-LLm)/2
          
C          if (subXI*NP_XI .ne. xi_rho-2) then
C            write(stdout,'(/8x,A,1x,A,I4,1x,A/15x,A,I3,1x,A/)')
C     &         'ERROR: Cannot partition XI-direction:', 'xi_rho-2 =',
C     &          xi_rho-2,  'is',  'not divisible by NP_XI =',  NP_XI,
C     &         'in nonperiodic XI case.'
C            goto 97
C          endif
        else
          write(stdout,'(/8x,2A/)') 'ERROR: inconsistent ',
     &                'dimensions ''xi_rho'' and ''xi_u''.' 
          goto 97     !--> next file
        endif

        if (eta_rho .eq. eta_v) then
          EtaPeriodic=.true.
          subETA=eta_rho/NP_ETA
      
          if (subETA*NP_ETA .ne. eta_rho) then
            write(stdout,'(/8x,A,1x,A,I4,1x,A/15x,A,I3,1x,A/)')
     &         'ERROR: Cannot partition ETA-direction:', 'eta_rho =',
     &          eta_rho,  'is',  'not divisible by NP_ETA =', NP_ETA,
     &         'in ETA-periodic case.'

            goto 97
          endif
        elseif (eta_rho .eq. eta_v+1) then
          EtaPeriodic=.false.
          subETA=(eta_rho-2)/NP_ETA
          
          chunk_size_E=(MMm+NP_ETA-1)/NP_ETA
          margin_E=(NP_ETA*chunk_size_E-MMm)/2          

C          if (subETA*NP_ETA .ne. eta_rho-2) then
C            write(stdout,'(/8x,A,1x,A,I4,1x,A/15x,A,I3,1x,A/)')
C     &        'ERROR: Cannot partition ETA-direction:','eta_rho-2 =',
C     &         eta_rho-2, 'is',  'not divisible by NP_ETA =', NP_ETA,
C     &        'in nonperiodic ETA case.'
C            goto 97
C          endif
        else
          write(stdout,'(/8x,A,1x,A/)') 'ERROR: inconsistent',
     &                 'dimensions ''eta_rho'' and ''eta_v''.'
          goto 97     !--> next file
        endif
!
! Create partitioned files.
! ====== =========== ======
!
        do node=0,NNODES-1
          lstr=lenstr(ncname0)
          ncname(node)=ncname0
          ierr=0
          call insert_node (ncname(node), lstr, node, NNODES, ierr)
          if (ierr. ne. 0) goto 97     !--> next file
!          ierr=nf_create (ncname(node)(1:lstr),nf_clobber,ncid(node))
          ierr=nf_create(ncname(node)(1:lstr),
     &                   nf_64bit_offset,ncid(node))
          if (ierr .eq. nf_noerr) then
            write(stdout,'(12x,3A)') 'Created partitioned file ''',
     &                                  ncname(node)(1:lstr), '''.' 
          else
            write(stdout,'(/8x,A,1x,3A/)') 'ERROR: cannot create',
     &              'netCDF file ''', ncname(node)(1:lstr), '''.'
            goto 97     !--> next file
          endif
!
! Define dimensions of partitioned files.
!
          jj=node/NP_XI
          ii=node-jj*NP_XI
          
      istrmpi=1+ii*chunk_size_X-margin_X
      iendmpi=istrmpi+chunk_size_X-1
      istrmpi=max(istrmpi,1)
      iendmpi=min(iendmpi,LLm)

      jstrmpi=1+jj*chunk_size_E-margin_E
      jendmpi=jstrmpi+chunk_size_E-1
      jstrmpi=max(jstrmpi,1)
      jendmpi=min(jendmpi,MMm)     
 
      Lmmpi=iendmpi-istrmpi+1
      Mmmpi=jendmpi-jstrmpi+1
      
          do i=1,ndims
            size=dimsize(i)
            if (i .eq. id_xi_rho) then
              size=Lmmpi
              if (.not.XiPeriodic) then
                if (ii.eq.0      ) size=size+1
                if (ii.eq.NP_XI-1) size=size+1
              endif
            elseif (i .eq. id_xi_u) then
              size=Lmmpi
              if (.not.XiPeriodic) then
                if (ii.eq.NP_XI-1) size=size+1 
              endif
            elseif (i .eq. id_eta_rho) then
              size=Mmmpi
              if (.not.EtaPeriodic) then
                if (jj.eq.0       ) size=size+1 
                if (jj.eq.NP_ETA-1) size=size+1 
              endif
            elseif (i .eq. id_eta_v) then
              size=Mmmpi
              if (.not.EtaPeriodic) then
                if (jj.eq.NP_ETA-1) size=size+1 
              endif
            endif
            if (i.ne.id_xi_psi   .and.  i.ne.id_xi_v  .and. 
     &          i.ne.id_eta_psi  .and.  i.ne.id_eta_u)  then
              lvar=lenstr(dimname(i))
              ierr=nf_def_dim (ncid(node), dimname(i)(1:lvar),
     &                                         size, dimid(i))
              if (ierr .ne. nf_noerr) then
                write(stdout,'(/8x,4A/15x,A,I4,A)') 'ERROR: ',
     &            'Cannot define dimension ''', dimname(i)(1:lvar),
     &            '''.',    'netCDF ettor status =',       ierr, '.'
              endif
c**           write(*,'(2I3,4x,2I3,I4,1x,A)') i, dimid(i), ii,
c**  &                            jj, size, dimname(i)(1:lvar)
            else
              dimid(i)=0
            endif
          enddo
!
! WARNING!!! ...After this moment array dimid(1:ndims) contains
! the set of NEW dimension IDs. Since the four dimensions, 'xi_psi',
! 'eta_psi', 'xi_v' and 'eta_u' have been eliminated, dimid(i) does
! not correspond to the set of dimension IDs of the original file
! [which would be just dimid(i)=i], but it is rather different.
! Array dimid(1:ndims) will be used later to remap old dimension
! IDs into new ones, see the remapping procedure approximately 80
! lines below.
!
! Put global attribute 'partition' which identifies subdomain
! within the processor grid individually for each file.
!
          ibuff(1)=ii
          ibuff(2)=jj
          ibuff(3)=NP_XI
          ibuff(4)=NP_ETA
          ierr=nf_put_att_int (ncid(node), nf_global, 'partition',
     &                                           nf_int, 4, ibuff)
        enddo
!
! Copy global attributes
!
        do i=1,ngatts 
          ierr=nf_inq_attname (ncid0, nf_global, i, string) 
          if (ierr. eq. nf_noerr) then
            lvar=lenstr(string)
            do node=0,NNODES-1
              ierr=nf_copy_att (ncid0, nf_global, string(1:lvar),
     &                                     ncid(node), nf_global)
              if (ierr. ne. nf_noerr) then
                lstr=lenstr(ncname(node))
                write(stdout,'(/8x,4A/15x,3A/)')  'ERROR: Cannot ',
     &            'copy global attribute ''', string(1:lvar), '''',
     &            'into netCDF file ''', ncname(node)(1:lstr),'''.'
                goto 97
              endif
            enddo
          else
            lstr=lenstr(ncname(0))
            write(stdout,'(/8x,2A,I3/15x,3A/)') 'ERROR: Cannot ',
     &         'determine mame of global attribute with ID =', i,
     &         'from netCDF file ''',    ncname0(1:lstr),   '''.'
            goto 97
          endif
        enddo
!
! Define variables and their attributes.
!
        do i=1,nvars
          ierr=nf_inq_var (ncid0,   i, varname(i),  vartype(i),
     &                      vardims(i), dimids(1,i),   varatts)
!
! Readress obsolete dimensions, if any:
!
          do j=1,vardims(i)
            if (dimids(j,i).eq.id_xi_psi) then
              dimids(j,i)=id_xi_u
            elseif (dimids(j,i).eq.id_xi_v) then
              dimids(j,i)=id_xi_rho
            elseif (dimids(j,i).eq.id_eta_psi) then
              dimids(j,i)=id_eta_v
            elseif (dimids(j,i).eq.id_eta_u) then
              dimids(j,i)=id_eta_rho
            endif
          enddo
!
! Determine whether partitionable dimensions or unlimited dimension
! are present for this variable. 
!
          series(i)=.false.
          part_switch(i)=.false.
          do j=1,vardims(i)
            if (dimids(j,i).eq.id_xi_rho .or.
     &          dimids(j,i).eq.id_xi_u    .or.
     &          dimids(j,i).eq.id_eta_rho .or.
     &          dimids(j,i).eq.id_eta_v) then
              part_switch(i)=.true. 
            elseif (dimids(j,i).eq.unlimdimid) then
              series(i)=.true. 
            endif
          enddo
!
! WARNING: Since dimids(1:vardims(i),i) contains dimension IDs
! corresponding to the set of IDs of the ORIGINAL file, and since
! some of the original dimensions were eliminated (merged), the
! set of dimension IDs in the NEW definitions is obtained by 
! inverse mapping of dimids(j,i) onto ibuff(j) using dimid(k) as
! a mapping array.  
!
          do j=1,vardims(i)
            do k=1,ndims
              if (dimids(j,i).eq.k) ibuff(j)=dimid(k)
            enddo
          enddo
c**       write(*,*) 'old_dimids:', (dimids(j,i),j=1,vardims(i))
c**       write(*,*) 'new_dimids:',    (ibuff(j),j=1,vardims(i))

          lvar=lenstr(varname(i))
          do node=0,NNODES-1
            ierr=nf_def_var (ncid(node), varname(i)(1:lvar),
     &              vartype(i), vardims(i), ibuff, varid(i))
          enddo
c**       write(stdout,'(I3,1x,A,T20,I3,1x,L1,1x,L1)') i,
c**  &                        varname(i)(1:lvar), vardims(i),
c**  &                            part_switch(i),  series(i)
          do j=1,varatts
            ierr=nf_inq_attname (ncid0, varid(i), j, string)
            lvar=lenstr(string)
            do node=0,NNODES-1
            ierr=nf_copy_att (ncid0, i, string(1:lvar),
     &                            ncid(node), varid(i))
            enddo
          enddo
        enddo
!
! Leave definition mode
!
        do node=0,NNODES-1
c          ierr=nf_set_fill (ncid(node), nf_nofill, i)
          ierr=nf_enddef(ncid(node))
        enddo
!
! Transfer variables into newly created files. 
!
        do rec=1,tsize
          if (tsize.gt.1) write(stdout,'(16x,A,I4,A)')
     &                 'Processing record', rec, '...'
          do i=1,nvars
            if (series(i) .or. rec.eq.1) then
              if (.not.part_switch(i) .and. .not.series(i)) then
!
! Scalar (zero-dimensional) variables:
!
                if (vartype(i) .eq. nf_char) then
                  ierr=nf_get_var_text (ncid0, i, buff)
                elseif (vartype(i) .eq. nf_int) then
                  ierr=nf_get_var_int    (ncid0, i, buff)
                elseif (vartype(i) .eq. nf_real) then
                  ierr=nf_get_var_real   (ncid0, i, buff)
                elseif (vartype(i) .eq. nf_double) then
                  ierr=nf_get_var_double (ncid0, i, buff)
                else
                  lvar=lenstr(varname(i))
                  write(stdout,'(/8x,4A/)') 'ERROR: scalar variable',
     &              ' ''', varname(i)(1:lvar), ''' has unknown type.'
                  goto 97
                endif
                if (ierr .eq. nf_noerr) then
                  do node=0,NNODES-1
                    if (vartype(i) .eq. nf_char) then
                     ierr=nf_put_var_text (ncid(node),varid(i),buff)
                    elseif (vartype(i) .eq. nf_int) then
                     ierr=nf_put_var_int   (ncid(node),varid(i),buff)
                    elseif (vartype(i) .eq. nf_real) then
                     ierr=nf_put_var_real  (ncid(node),varid(i),buff)
                    elseif (vartype(i) .eq. nf_double) then
                     ierr=nf_put_var_double(ncid(node),varid(i),buff)
                    endif
                    if (ierr .ne. nf_noerr) then
                      lvar=lenstr(varname(i))
                      lstr=lenstr(ncname(node))
                      write(stdout,'(/8x,3A/15x,3A,I4,A/)')
     &                    'ERROR: Cannot write scalar variable ''',
     &                     varname(i)(1:lvar), ''' into netCDF file',
     &                    '''',  ncname(node)(1:lstr), 
     &                    '''.     netCDF error code =',  ierr,   '.'
                      goto 97
                    endif
                  enddo
                else
                  lvar=lenstr(varname(i))
                  write(stdout,'(/8x,4A/)') 'ERROR: Cannot read ',
     &             'scalar variable ''', varname(i)(1:lvar), '''.'
                  goto 97
                endif
              elseif (.not.part_switch(i)) then
!
! Non-partitionable array.
!
                size=1
                do j=1,vardims(i) 
                  if (dimids(j,i).eq.unlimdimid) then
                    start(j)=rec
                    count(j)=1
                  else
                    start(j)=1
                    count(j)=dimsize(dimids(j,i))
                  endif 
                  size=size*count(j)
                enddo
                if (vartype(i) .eq. nf_char) then
                  size=size*1
                elseif (vartype(i) .eq. nf_int) then
                  size=size*4
                elseif (vartype(i) .eq. nf_real) then
                  size=size*4
                elseif (vartype(i) .eq. nf_double) then
                  size=size*8
                else
                  lvar=lenstr(varname(i))
                  write(stdout,'(/8x,3A/)') 'ERROR: variable ''',
     &                 varname(i)(1:lvar), ''' has unknown type.'
                  goto 97
                endif
                if (size .gt. 8*max_buff_size) then
                  write(stdout,'(/8x,A,3(/15x,A,I10,1x,A)/)')
     &              'ERROR: unsufficient buffer size in "partit.F":',
     &              'requested:',         size,      'Bytes,',
     &              'available:',   8*max_buff_size, 'Bytes.',
     &              'Increase parameter max_buff_size and recompile.'
                  goto 97
                endif

                if (vartype(i) .eq. nf_char) then
                  ierr=nf_get_vara_text   (ncid0, i, start,
     &                                          count, buff)
                elseif (vartype(i) .eq. nf_int) then
                  ierr=nf_get_vara_int    (ncid0, i, start,
     &                                         count, buff)
                elseif (vartype(i) .eq. nf_real) then
                  ierr=nf_get_vara_real   (ncid0, i, start,
     &                                         count, buff)
                elseif (vartype(i) .eq. nf_double) then
                  ierr=nf_get_vara_double (ncid0, i, start,
     &                                         count, buff)
                endif
                if (ierr .eq. nf_noerr) then
                  do node=0,NNODES-1
                    if (vartype(i) .eq. nf_char) then
                      ierr=nf_put_vara_text   (ncid(node), varid(i),
     &                                           start, count, buff)
                    elseif (vartype(i) .eq. nf_int) then
                      ierr=nf_put_vara_int    (ncid(node), varid(i),
     &                                           start, count, buff)
                    elseif (vartype(i) .eq. nf_real) then
                      ierr=nf_put_vara_real   (ncid(node), varid(i),
     &                                           start, count, buff)
                    elseif (vartype(i) .eq. nf_double) then
                      ierr=nf_put_vara_double (ncid(node), varid(i),
     &                                           start, count, buff)
                    endif
                    if (ierr .ne. nf_noerr) then
                      lvar=lenstr(varname(i))
                      lstr=lenstr(ncname(node))
                      write(stdout,'(/8x,3A,I3/15x,3A,I4,A/)')
     &                  'ERROR: Cannot write variable ''',
     &                   varname(i)(1:lvar),''' for time record',rec,
     &                  'into netCDF file ''',  ncname(node)(1:lstr),
     &                  '''. netCDF error code =', ierr, '.'
                      goto 97
                    endif
                  enddo
                else
                  lstr=lenstr(ncname0)
                  lvar=lenstr(varname(i))
                  write(stdout,'(/8x,4A,I3/15x,3A,I4/)') 'ERROR: ',
     &              'Cannot read variable ''',  varname(i)(1:lvar),
     &              ''' for time record',rec,'from netCDF file ''',
     &              ncname0(1:lstr),'''. netCDF error code =',ierr
                  goto 97
                endif
              elseif (part_switch(i)) then
!
! Partitioned array: 
!
                do node=0,NNODES-1
                  jj=node/NP_XI
                  ii=node-jj*NP_XI
                  
      istrmpi=1+ii*chunk_size_X-margin_X
      iendmpi=istrmpi+chunk_size_X-1
      istrmpi=max(istrmpi,1)
      iendmpi=min(iendmpi,LLm)

      jstrmpi=1+jj*chunk_size_E-margin_E
      jendmpi=jstrmpi+chunk_size_E-1
      jstrmpi=max(jstrmpi,1)
      jendmpi=min(jendmpi,MMm)     
 
      Lmmpi=iendmpi-istrmpi+1
      Mmmpi=jendmpi-jstrmpi+1
      
                  size=1
                  do j=1,vardims(i)
                    if (dimids(j,i).eq.id_xi_rho) then
                      start(j)=istrmpi
                      count(j)=Lmmpi
                      if (.not.XiPeriodic) then
                        if (ii.gt.0      ) start(j)=start(j)+1
                        if (ii.eq.0      ) count(j)=count(j)+1
                        if (ii.eq.NP_XI-1) count(j)=count(j)+1
                      endif
                      start1(j)=1

                    elseif (dimids(j,i).eq.id_xi_u) then
                      start(j)=istrmpi
                      count(j)=Lmmpi
                      if (.not.XiPeriodic) then
                        if (ii.eq.NP_XI-1) count(j)=count(j)+1
                      endif
                      start1(j)=1

                    elseif (dimids(j,i).eq.id_eta_rho) then
                      start(j)=jstrmpi
                      count(j)=Mmmpi
                      if (.not.EtaPeriodic) then
                        if (jj.gt.0       ) start(j)=start(j)+1
                        if (jj.eq.0       ) count(j)=count(j)+1
                        if (jj.eq.NP_ETA-1) count(j)=count(j)+1
                      endif
                      start1(j)=1

                    elseif (dimids(j,i).eq.id_eta_v) then
                      start(j)=jstrmpi
                      count(j)=Mmmpi
                      if (.not.EtaPeriodic) then
                        if (jj.eq.NP_ETA-1) count(j)=count(j)+1
                      endif
                      start1(j)=1

                    elseif (dimids(j,i).eq.unlimdimid) then
                      start(j)=rec
                      count(j)=1
                      start1(j)=rec
                    else
                      start(j)=1
                      count(j)=dimsize(dimids(j,i))
                      start1(j)=1
                    endif
                    size=size*count(j)
                  enddo
c**               write(*,*) 'dimids:', (dimids(j,i),j=1,vardims(i))
c**               write(*,*) ' start:',    (start(j),j=1,vardims(i))
c**               write(*,*) ' count:',    (count(j),j=1,vardims(i))


                  if (vartype(i) .eq. nf_char) then
                    size=size*1
                  elseif (vartype(i) .eq. nf_int) then
                    size=size*4
                  elseif (vartype(i) .eq. nf_real) then
                    size=size*4
                  elseif (vartype(i) .eq. nf_double) then
                    size=size*8
                  else
                    lvar=lenstr(varname(i))
                    write(stdout,'(/8x,4A/)') 'ERROR: variable ''',
     &                   varname(i)(1:lvar), ''' has unknown type.'
                    goto 97
                  endif
                  if (size .gt. 8*max_buff_size) then
                    write(stdout,'(/8x,A,3(/15x,A,I10,1x,A)/)')
     &              'ERROR: unsufficient buffer size in "partit.F":',
     &              'requested:',         size,      'Bytes,',
     &              'available:',   8*max_buff_size, 'Bytes.',
     &              'Increase parameter max_buff_size and recompile.'
                    goto 97
                  endif

                  if (vartype(i) .eq. nf_char) then
                    ierr=nf_get_vara_text   (ncid0, i, start,
     &                                           count, buff)
                  elseif (vartype(i) .eq. nf_int) then
                    ierr=nf_get_vara_int    (ncid0, i, start,
     &                                           count, buff)
                  elseif (vartype(i) .eq. nf_real) then
                    ierr=nf_get_vara_real   (ncid0, i, start,
     &                                           count, buff)
                  elseif (vartype(i) .eq. nf_double) then
                    ierr=nf_get_vara_double (ncid0, i, start,
     &                                           count, buff)
                  endif

                  if (ierr .eq. nf_noerr) then
                    if (vartype(i) .eq. nf_char) then
                      ierr=nf_put_vara_text   (ncid(node), varid(i),
     &                                          start1, count, buff)
                    elseif (vartype(i) .eq. nf_int) then
                      ierr=nf_put_vara_int    (ncid(node), varid(i),
     &                                          start1, count, buff)
                    elseif (vartype(i) .eq. nf_real) then
                      ierr=nf_put_vara_real   (ncid(node), varid(i),
     &                                          start1, count, buff)
                    elseif (vartype(i) .eq. nf_double) then
                      ierr=nf_put_vara_double (ncid(node), varid(i),
     &                                          start1, count, buff)
                    endif
                    if (ierr .ne. nf_noerr) then
                      lvar=lenstr(varname(i))
                      lstr=lenstr(ncname(node))
                      write(stdout,'(/8x,3A,I3/15x,3A,I4,A/)')
     &                  'ERROR: Cannot write partitioned array ''',
     &                   varname(i)(1:lvar),''' for time record',rec,
     &                  'into netCDF file ''',  ncname(node)(1:lstr),
     &                  '''. netCDF error code =', ierr, '.'
                      goto 97
                    endif
                  else
                    lstr=lenstr(ncname0)
                    lvar=lenstr(varname(i))
                    write(stdout,'(/8x,3A,I3/15x,3A,I4,A/)')
     &                  'ERROR: Cannot read partitioned array ''',
     &                   varname(i)(1:lvar),   ''' for time record',
     &                   rec, 'from netCDF file ''',ncname0(1:lstr),
     &                  '''. netCDF error code =',   ierr,    '.'

                    goto 97
                  endif
                enddo       ! <-- node=0,NNODES-1
              endif
            endif       ! <--series(i) .or. rec.eq.1
          enddo       ! <-- i=1,nvars
        enddo       ! <-- rec=1,tsize
!
! Close all netCDF files
!
  97    ierr=nf_close (ncid0)
        do node=0,NNODES-1
         ierr=nf_close (ncid(node))
        enddo
      enddo
      stop
      end


