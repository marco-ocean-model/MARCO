      subroutine t3dmix (tile)
      implicit none
      integer*4 tile, itrc, trd, omp_get_thread_num
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=41,   MMm0=42,   N=32)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 N_sl
      parameter (N_sl=0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=1,  NP_ETA=4,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 stdout, Np, NpHz, padd_X,padd_E
      parameter (stdout=6)
      parameter (Np=N+1)
      parameter (NpHz=(N+1+N_sl))
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d,N3dHz, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      parameter (N3dHz=size_XI*size_ETA*NpHz)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, NTA, itemp, NTot
      integer*4   ntrc_temp, ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      integer*4   ntrc_subs, ntrc_substot
      parameter (itemp=1)
      parameter (ntrc_temp=1)
      parameter (ntrc_salt=1)
      parameter (ntrc_pas=0)
      parameter (ntrc_bio=0)
      parameter (ntrc_subs=0, ntrc_substot=0)
      parameter (ntrc_sed=0)
      parameter (NTA=itemp+ntrc_salt)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed)
      parameter (NTot=NT)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_diapv
      integer*4   ntrc_diaeddy, ntrc_surf
     &          , isalt
      parameter (isalt=itemp+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_diapv=0)
      parameter (ntrc_diaeddy=0)
      parameter (ntrc_surf=0)
      real A2d(N2d,NSA,0:NPP-1), A3d(N3d,9,0:NPP-1)
     &    ,A3dHz(N3dHz,4,0:NPP-1)
      integer*4 B2d(N2d,0:NPP-1)
      common/private_scratch/ A2d,A3d,A3dHz
      common/private_scratch_bis/ B2d
      integer*4 chunk_size_X,margin_X,chunk_size_E,margin_E
      integer*4 Istr,Iend,Jstr,Jend, i_X,j_E
      chunk_size_X=(Lmmpi+NSUB_X-1)/NSUB_X
      margin_X=(NSUB_X*chunk_size_X-Lmmpi)/2
      chunk_size_E=(Mmmpi+NSUB_E-1)/NSUB_E
      margin_E=(NSUB_E*chunk_size_E-Mmmpi)/2
      j_E=tile/NSUB_X
      i_X=tile-j_E*NSUB_X
      Istr=1+i_X*chunk_size_X-margin_X
      Iend=Istr+chunk_size_X-1
      Istr=max(Istr,1)
      Iend=min(Iend,Lmmpi)
      Jstr=1+j_E*chunk_size_E-margin_E
      Jend=Jstr+chunk_size_E-1
      Jstr=max(Jstr,1)
      Jend=min(Jend,Mmmpi)
      trd=omp_get_thread_num()
      do itrc=1,NT
        call t3dmix_tile (istr,iend,jstr,jend, itrc,   A3d(1,1,trd),
     &                                                 A3d(1,2,trd),
     &                                    A3d(1,3,trd),A3d(1,4,trd),
     &                                    A2d(1,1,trd), A2d(1,2,trd)
     &                   ,A2d(1,10,trd),A2d(1,11,trd),A2d(1,12,trd),
     &                                   A2d(1,13,trd),A2d(1,14,trd)
     &          ,A3d(1,5,trd),A3d(1,6,trd),A3d(1,7,trd),A3d(1,8,trd)
     &                    )
      enddo
      return
      end
      subroutine t3dmix_tile (istr,iend,jstr,jend, itrc, LapT,
     &                                                          Akz,
     &                                                diff3u,diff3v,
     &                                      FX,FE
     &                                              ,FFC,CF,BC,CD,DC
     &    ,dTdr, dTdx, dTde, FC                    )
      implicit none
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=41,   MMm0=42,   N=32)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 N_sl
      parameter (N_sl=0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=1,  NP_ETA=4,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 stdout, Np, NpHz, padd_X,padd_E
      parameter (stdout=6)
      parameter (Np=N+1)
      parameter (NpHz=(N+1+N_sl))
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d,N3dHz, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      parameter (N3dHz=size_XI*size_ETA*NpHz)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, NTA, itemp, NTot
      integer*4   ntrc_temp, ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      integer*4   ntrc_subs, ntrc_substot
      parameter (itemp=1)
      parameter (ntrc_temp=1)
      parameter (ntrc_salt=1)
      parameter (ntrc_pas=0)
      parameter (ntrc_bio=0)
      parameter (ntrc_subs=0, ntrc_substot=0)
      parameter (ntrc_sed=0)
      parameter (NTA=itemp+ntrc_salt)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed)
      parameter (NTot=NT)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_diapv
      integer*4   ntrc_diaeddy, ntrc_surf
     &          , isalt
      parameter (isalt=itemp+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_diapv=0)
      parameter (ntrc_diaeddy=0)
      parameter (ntrc_surf=0)
      integer*4 istr,iend,jstr,jend, itrc, i,j,k,k1,k2, kmld,
     &        imin,imax,jmin,jmax, indx, idx,ide
      real LapT (Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &   diff3u (Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &   diff3v (Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &      Akz (Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &       CF (Istr-2:Iend+2,0:N),
     &       DC (Istr-2:Iend+2,0:N),
     &       CD (Istr-2:Iend+2,0:N),
     &       BC (Istr-2:Iend+2,0:N),
     &       FFC(Istr-2:Iend+2,0:N),
     &       FX (Istr-2:Iend+2,Jstr-2:Jend+2),
     &       FE (Istr-2:Iend+2,Jstr-2:Jend+2),
     &       FC (Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &     dTdr (Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &     dTdx (Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &     dTde (Istr-2:Iend+2,Jstr-2:Jend+2,0:N)
       real cff,cff1,
     &      TRIADS1,TRIADS2,TRIADS3,TRIADS4,sumX,sumE,sig,
     &      SLOPEXQ1,SLOPEXQ2,SLOPEXQ3,SLOPEXQ4,
     &      SLOPEYQ1,SLOPEYQ2,SLOPEYQ3,SLOPEYQ4
       real, dimension(0:4), parameter :: wgt
     &   =(/0.D0,1.D0,0.5D0,0.33333333333D0,0.25D0/)
      real h(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real hinv(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real f(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real fomn(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grid_h/h /grid_hinv/hinv /grid_f/f /grid_fomn/fomn
      real angler(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grid_angler/angler
      real latr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real lonr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real latu(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real lonu(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real latv(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real lonv(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grid_latr/latr /grid_lonr/lonr
      common /grid_latu/latu /grid_lonu/lonu
      common /grid_latv/latv /grid_lonv/lonv
      real pm(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pn(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real om_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real on_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real om_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real on_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real om_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real on_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real om_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real on_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pn_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pm_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pm_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pn_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /metrics_pm/pm    /metrics_pn/pn
      common /metrics_omr/om_r /metrics_on_r/on_r
      common /metrics_omu/om_u /metrics_on_u/on_u
      common /metrics_omv/om_v /metrics_on_v/on_v
      common /metrics_omp/om_p /metrics_on_p/on_p
      common /metrics_pnu/pn_u /metrics_pmv/pm_v
      common /metrics_pmu/pm_u /metrics_pnv/pn_v
      real dmde(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real dndx(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /metrics_dmde/dmde    /metrics_dndx/dndx
      real pmon_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pmon_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pmon_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pnom_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pnom_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pnom_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real grdscl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /metrics_pmon_p/pmon_p /metrics_pnom_p/pnom_p
      common /metrics_pmon_r/pmon_r /metrics_pnom_r/pnom_r
      common /metrics_pmon_u/pmon_u /metrics_pnom_v/pnom_v
      common /metrics_grdscl/grdscl
      real rmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real umask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real vmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pmask2(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mask_r/rmask
      common /mask_p/pmask
      common /mask_u/umask
      common /mask_v/vmask
      common /mask_p2/pmask2
      real zob(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /Z0B_VAR/zob
      real u(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3)
      real v(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3)
      real t(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t
      real Hz(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real Hz_bak(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real z_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real z_w(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real Huon(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real Hvom(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /grid_Hz_bak/Hz_bak /grid_zw/z_w /grid_Huon/Huon
      common /grid_Hvom/Hvom
      real We(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /grid_Hz/Hz /grid_zr/z_r /grid_We/We
      real rho1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real rho(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /ocean_rho1/rho1 /ocean_rho/rho
      real qp1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /ocean_qp1/qp1
      real qp2
      parameter (qp2=0.0000172D0)
      real visc2_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_sponge_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_sponge_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mixing_visc2_r/visc2_r /mixing_visc2_p/visc2_p
      common /mixing_visc2_sponge_r/visc2_sponge_r
      common /mixing_visc2_sponge_p/visc2_sponge_p
      real diff2_sponge(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real diff2(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /mixing_diff2_sponge/diff2_sponge
      common /mixing_diff2/diff2
      real diff4_sponge(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real diff4(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /mixing_diff4_sponge/diff4_sponge
      common /mixing_diff4/diff4
      real diff3d_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real diff3d_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /mixing_diff3d_u/diff3d_u
      common /mixing_diff3d_v/diff3d_v
      real dRdx(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real dRde(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real idRz(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /mixing_dRdx/dRdx
      common /mixing_dRde/dRde
      common /mixing_idRz/idRz
      real Rslope_max,Gslope_max
      parameter (Gslope_max=1.D0, Rslope_max=0.05D0)
      integer*4 ismooth
      real csmooth
      common /mixing_csmooth/ csmooth
      common /mixing_ismooth/ ismooth
      real Akv(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real Akt(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N,2)
      common /mixing_Akv/Akv /mixing_Akt/Akt
      real bvf(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /mixing_bvf/ bvf
      real ustar(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /lmd_kpp_ustar/ustar
      integer*4 kbl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      integer*4 kbbl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real hbbl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /lmd_kpp_kbl/ kbl
      common /lmd_kpp_hbbl/ hbbl
      common /lmd_kpp_kbbl/ kbbl
      real hbls(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /lmd_kpp_hbl/ hbls
      real ghats(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /lmd_kpp_ghats/ghats
      real dt, dtfast, time, time2, time_start, tdays, start_time
      integer*4 ndtfast, iic, kstp, krhs, knew, next_kstp
     &      , iif, nstp, nrhs, nnew, nbstep3d
      logical PREDICTOR_2D_STEP
      common /time_indices/  dt,dtfast, time, time2,time_start, tdays,
     &     ndtfast, iic, kstp, krhs, knew, next_kstp,
     &     start_time,
     &                       iif, nstp, nrhs, nnew, nbstep3d,
     &                       PREDICTOR_2D_STEP
      real time_avg, time2_avg, rho0
     &               , rdrg, rdrg2, Cdb_min, Cdb_max, Zobt
     &               , xl, el, visc2, visc4, gamma2
      real  theta_s,   theta_b,   Tcline,  hc
      real  sc_w(0:N), Cs_w(0:N), sc_r(N), Cs_r(N)
      real  rx0, rx1
      real  tnu2(NT),tnu4(NT)
      real weight(6,0:NWEIGHT)
      real  x_sponge,   v_sponge
       real  tauT_in, tauT_out, tauM_in, tauM_out
      integer*4 numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
      logical ldefhis
      logical got_tini(NT)
      common /scalars_main/
     &             time_avg, time2_avg,  rho0,      rdrg,    rdrg2
     &           , Zobt,       Cdb_min,   Cdb_max
     &           , xl, el,    visc2,     visc4,   gamma2
     &           , theta_s,   theta_b,   Tcline,  hc
     &           , sc_w,      Cs_w,      sc_r,    Cs_r
     &           , rx0,       rx1
     &           ,       tnu2,    tnu4
     &                      , weight
     &                      , x_sponge,   v_sponge
     &                      , tauT_in, tauT_out, tauM_in, tauM_out
     &      , numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
     &                      , got_tini
     &                      , ldefhis
      logical synchro_flag
      common /sync_flag/ synchro_flag
      integer*4 may_day_flag
      integer*4 tile_count, first_time, bc_count
      common /communicators_i/
     &        may_day_flag, tile_count, first_time, bc_count
      real hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      common /communicators_r/
     &     hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      real lonmin, lonmax, latmin, latmax
      common /communicators_lonlat/
     &     lonmin, lonmax, latmin, latmax
      real*8 Cu_Adv3d,  Cu_W, Cu_Nbq_X, Cu_Nbq_Y, Cu_Nbq_Z
      integer*4 i_cx_max, j_cx_max, k_cx_max
      common /diag_vars/ Cu_Adv3d,  Cu_W,
     &        i_cx_max, j_cx_max, k_cx_max
      real*8 volume, avgke, avgpe, avgkp, bc_crss
      common /communicators_rq/
     &          volume, avgke, avgpe, avgkp, bc_crss
      real*4 CPU_time(0:31,0:NPP)
      integer*4 proc(0:31,0:NPP),trd_count
      common /timers_roms/CPU_time,proc,trd_count
      logical EAST_INTER2, WEST_INTER2, NORTH_INTER2, SOUTH_INTER2
      logical EAST_INTER, WEST_INTER, NORTH_INTER, SOUTH_INTER
      logical CORNER_SW,CORNER_NW,CORNER_NE,CORNER_SE
      integer*4 mynode, mynode2, ii,jj, p_W,p_E,p_S,p_N, p_SW,p_SE,
     & p_NW,p_NE,NNODES2
      common /comm_setup/ mynode, mynode2, ii,jj, p_W,p_E,p_S,p_N,
     & p_SW,p_SE, p_NW,p_NE, EAST_INTER, WEST_INTER, NORTH_INTER,
     & SOUTH_INTER, EAST_INTER2, WEST_INTER2, NORTH_INTER2, 
     &                                                     SOUTH_INTER2,
     & CORNER_SW,CORNER_NW,CORNER_NE,CORNER_SE,NNODES2
      real pi, deg2rad, rad2deg
      parameter (pi=3.14159265358979323846D0, deg2rad=pi/180.D0,
     &                                      rad2deg=180.D0/pi)
      real Eradius, Erotation, g, day2sec,sec2day, jul_off,
     &     year2day,day2year
      parameter (Eradius=6371315.0D0,  Erotation=7.292115090D-5,
     &           day2sec=86400.D0, sec2day=1.D0/86400.D0,
     &           year2day=365.25D0, day2year=1.D0/365.25D0,
     &           jul_off=2440000.D0)
      parameter (g=9.81D0)
      real Cp
      parameter (Cp=3985.0D0)
      real vonKar
      parameter (vonKar=0.41D0)
      real spval
      parameter (spval=-999.0D0)
      logical mask_val
      parameter (mask_val = .true.)
      integer*4 IstrR,IendR,JstrR,JendR
      integer*4 IstrU
      integer*4 JstrV
      if ((istr.eq.1 .and. .not.WEST_INTER)) then
        IstrR=Istr-1
        IstrU=Istr+1
      else
        IstrR=Istr
        IstrU=Istr
      endif
      if ((iend.eq.Lmmpi .and. .not.EAST_INTER)) then
        IendR=Iend+1
      else
        IendR=Iend
      endif
      if ((jstr.eq.1.and. .not.SOUTH_INTER)) then
        JstrR=Jstr-1
        JstrV=Jstr+1
      else
        JstrR=Jstr
        JstrV=Jstr
      endif
      if ((jend.eq.Mmmpi .and. .not.NORTH_INTER)) then
        JendR=Jend+1
      else
        JendR=Jend
      endif
      if ((istr.eq.1 .and. .not.WEST_INTER)) then
        imin=istr
      else
        imin=istr-1
      endif
      if ((iend.eq.Lmmpi .and. .not.EAST_INTER)) then
        imax=iend
      else
        imax=iend+1
      endif
      if ((jstr.eq.1.and. .not.SOUTH_INTER)) then
        jmin=jstr
      else
        jmin=jstr-1
      endif
      if ((jend.eq.Mmmpi .and. .not.NORTH_INTER)) then
        jmax=jend
      else
        jmax=jend+1
      endif
!$acc kernels if(compute_on_device) default(present)
      do k=1,N
        do j=jmin,jmax
          do i=imin,imax+1
            diff3u(i,j,k)=
     &                     +sqrt(
     &                      0.5D0*(diff4(i,j,itrc)+diff4(i-1,j,itrc))
     &                     +diff3d_u(i,j,k)
     &                           )
          enddo
        enddo
        do j=jmin,jmax+1
          do i=imin,imax
            diff3v(i,j,k)=
     &                     +sqrt(
     &                      0.5D0*(diff4(i,j,itrc)+diff4(i,j-1,itrc))
     &                     +diff3d_v(i,j,k)
     &                           )
          enddo
        enddo
      enddo
      do k=1,N
        do j=jmin,jmax
          do i=imin,imax+1
            cff=0.5D0*(pm(i,j)+pm(i-1,j)) * umask(i,j)
            dTdx(i,j,k)=cff*( t(i  ,j,k,nstp,itrc)
     &                          -t(i-1,j,k,nstp,itrc)
     &                         )
          enddo
        enddo
        do j=jmin,jmax+1
          do i=imin,imax
            cff=0.5D0*(pn(i,j)+pn(i,j-1)) * vmask(i,j)
            dTde(i,j,k)=cff*( t(i,j  ,k,nstp,itrc)
     &                          -t(i,j-1,k,nstp,itrc)
     &                         )
          enddo
        enddo
      enddo
      do k=1,N-1
        do j=jmin-1,jmax+1
          do i=imin-1,imax+1
            dTdr(i,j,k)= idRz(i,j,k)*( t(i,j,k+1,nstp,itrc)
     &                                  - t(i,j,k  ,nstp,itrc)
     &                                  )
          enddo
        enddo
      enddo
      do j=jmin-1,jmax+1
        do i=imin-1,imax+1
          dTdr(i,j,0)= dTdr(i,j,1)
          dTdr(i,j,N)= dTdr(i,j,N-1)
        enddo
      enddo
      do k=1,N-1
        do j=jmin,jmax
          do i=imin,imax
            cff  = idRz(i,j,k)*( z_r (i,j,k+1)-z_r (i,j,k) )
                TRIADS1=dRdx(i  ,j,k  )*dTdr(i,j,k)-dTdx(i  ,j,k)
                TRIADS2=dRdx(i  ,j,k+1)*dTdr(i,j,k)-dTdx(i  ,j,k+1)
                TRIADS3=dRdx(i+1,j,k+1)*dTdr(i,j,k)-dTdx(i+1,j,k+1)
                TRIADS4=dRdx(i+1,j,k  )*dTdr(i,j,k)-dTdx(i+1,j,k)
                sumX = diff3u(i  ,j,k  )*dRdx(i  ,j,k  )*TRIADS1
     &               + diff3u(i  ,j,k+1)*dRdx(i  ,j,k+1)*TRIADS2
     &               + diff3u(i+1,j,k+1)*dRdx(i+1,j,k+1)*TRIADS3
     &               + diff3u(i+1,j,k  )*dRdx(i+1,j,k  )*TRIADS4
                idx = 4
                TRIADS1=dRde(i,j  ,k  )*dTdr(i,j,k)-dTde(i,j  ,k)
                TRIADS2=dRde(i,j  ,k+1)*dTdr(i,j,k)-dTde(i,j  ,k+1)
                TRIADS3=dRde(i,j+1,k+1)*dTdr(i,j,k)-dTde(i,j+1,k+1)
                TRIADS4=dRde(i,j+1,k  )*dTdr(i,j,k)-dTde(i,j+1,k)
                sumE = diff3v(i,j  ,k  )*dRde(i,j  ,k  )*TRIADS1
     &               + diff3v(i,j  ,k+1)*dRde(i,j  ,k+1)*TRIADS2
     &               + diff3v(i,j+1,k+1)*dRde(i,j+1,k+1)*TRIADS3
     &               + diff3v(i,j+1,k  )*dRde(i,j+1,k  )*TRIADS4
                ide = 4
                FC(i,j,k)=(sumX*wgt(idx)+sumE*wgt(ide))*cff
          enddo
        enddo
      enddo
      do j=jmin,jmax
        do i=imin,imax
          FC(i,j,0)  = 0.D0
          FC(i,j,N)  = 0.D0
          Akz (i,j,0 )= 0.0D0
          Akz (i,j,N )= 0.0D0
        enddo
      enddo
!$acc loop seq
      do k=1,N
          cff=0.5D0
          do j=jmin,jmax
            do i=imin,imax+1
              FX(i,j)=cff*diff3u(i,j,k)*(Hz(i,j,k)+Hz(i-1,j,k))
     &               *on_u(i,j)*(   dTdx(i,j,k) -
     &                0.25D0*dRdx(i,j,k)*( dTdr(i-1,j,k-1)+dTdr(i,j,k)
     &                                 + dTdr(i-1,j,k)+dTdr(i,j,k-1))
     &                                                              )
            enddo
          enddo
          do j=jmin,jmax+1
            do i=imin,imax
              FE(i,j)=cff*diff3v(i,j,k)*(Hz(i,j,k)+Hz(i,j-1,k))
     &               *om_v(i,j)*(   dTde(i,j,k) -
     &              0.25D0*dRde(i,j,k)*(   dTdr(i,j-1,k-1)+dTdr(i,j,k)
     &                              +    dTdr(i,j-1,k)+dTdr(i,j,k-1))
     &                                                              )
            enddo
          enddo
          do j=jmin,jmax
            do i=imin,imax
              LapT(i,j,k)=( pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                       +FE(i,j+1)-FE(i,j))
     &                     +FC(i,j,k)-FC(i,j,k-1)    )/Hz(i,j,k)
            enddo
          enddo
      enddo
        if ((istr.eq.1 .and. .not.WEST_INTER)) then
          do k=1,N
            do j=jmin,jmax
              LapT(istr-1,j,k)=LapT(istr,j,k)
            enddo
          enddo
        endif
        if ((iend.eq.Lmmpi .and. .not.EAST_INTER)) then
          do k=1,N
            do j=jmin,jmax
              LapT(iend+1,j,k)=LapT(iend,j,k)
            enddo
          enddo
        endif
        if ((jstr.eq.1.and. .not.SOUTH_INTER)) then
          do k=1,N
            do i=imin,imax
              LapT(i,jstr-1,k)=LapT(i,jstr,k)
            enddo
          enddo
        endif
        if ((jend.eq.Mmmpi .and. .not.NORTH_INTER)) then
          do k=1,N
            do i=imin,imax
              LapT(i,jend+1,k)=LapT(i,jend,k)
            enddo
          enddo
        endif
      do k=1,N
        do j=jstr,jend
          do i=istr,iend+1
            cff=0.5D0*(pm(i,j)+pm(i-1,j)) * umask(i,j)
            dTdx(i,j,k)=cff*(LapT(i,j,k)-LapT(i-1,j,k))
          enddo
        enddo
        do j=jstr,jend+1
          do i=istr,iend
            cff=0.5D0*(pn(i,j)+pn(i,j-1)) * vmask(i,j)
            dTde(i,j,k)=cff*(LapT(i,j,k)-LapT(i,j-1,k))
          enddo
        enddo
      enddo
      do k=1,N-1
        do j=jstr-1,jend+1
          do i=istr-1,iend+1
            dTdr(i,j,k)= idRz(i,j,k)*( LapT(i,j,k+1)-LapT(i,j,k) )
          enddo
        enddo
      enddo
      do j=jstr-1,jend+1
        do i=istr-1,iend+1
          dTdr(i,j,0)= dTdr(i,j,1)
          dTdr(i,j,N)=0.D0
        enddo
      enddo
      do k=1,N-1
        do j=jstr,jend
          do i=istr,iend
            FC(i,j,k)  = idRz(i,j,k)*( z_r (i,j,k+1)-z_r (i,j,k) )
                TRIADS1=dRdx(i  ,j,k  )*dTdr(i,j,k)-dTdx(i  ,j,k)
                TRIADS2=dRdx(i  ,j,k+1)*dTdr(i,j,k)-dTdx(i  ,j,k+1)
                TRIADS3=dRdx(i+1,j,k+1)*dTdr(i,j,k)-dTdx(i+1,j,k+1)
                TRIADS4=dRdx(i+1,j,k  )*dTdr(i,j,k)-dTdx(i+1,j,k)
                sumX = diff3u(i  ,j,k  )*dRdx(i  ,j,k  )*TRIADS1
     &               + diff3u(i  ,j,k+1)*dRdx(i  ,j,k+1)*TRIADS2
     &               + diff3u(i+1,j,k+1)*dRdx(i+1,j,k+1)*TRIADS3
     &               + diff3u(i+1,j,k  )*dRdx(i+1,j,k  )*TRIADS4
                idx = 4
                TRIADS1=dRde(i,j  ,k  )*dTdr(i,j,k)-dTde(i,j  ,k)
                TRIADS2=dRde(i,j  ,k+1)*dTdr(i,j,k)-dTde(i,j  ,k+1)
                TRIADS3=dRde(i,j+1,k+1)*dTdr(i,j,k)-dTde(i,j+1,k+1)
                TRIADS4=dRde(i,j+1,k  )*dTdr(i,j,k)-dTde(i,j+1,k)
                sumE = diff3v(i,j  ,k  )*dRde(i,j  ,k  )*TRIADS1
     &               + diff3v(i,j  ,k+1)*dRde(i,j  ,k+1)*TRIADS2
     &               + diff3v(i,j+1,k+1)*dRde(i,j+1,k+1)*TRIADS3
     &               + diff3v(i,j+1,k  )*dRde(i,j+1,k  )*TRIADS4
                ide = 4
                SLOPEXQ1=(FC(i,j,k)*dRdx(i  ,j,k  ))**2
                SLOPEXQ2=(FC(i,j,k)*dRdx(i  ,j,k+1))**2
                SLOPEXQ3=(FC(i,j,k)*dRdx(i+1,j,k+1))**2
                SLOPEXQ4=(FC(i,j,k)*dRdx(i+1,j,k  ))**2
                SLOPEYQ1=(FC(i,j,k)*dRde(i,j  ,k  ))**2
                SLOPEYQ2=(FC(i,j,k)*dRde(i,j  ,k+1))**2
                SLOPEYQ3=(FC(i,j,k)*dRde(i,j+1,k+1))**2
                SLOPEYQ4=(FC(i,j,k)*dRde(i,j+1,k  ))**2
                cff = 1.D0/(z_r(i,j,k+1)-z_r(i,j,k))
                Akz(i,j,k) = 8.D0*( max(
     &                       diff3u(i  ,j,k  )*SLOPEXQ1,
     &                       diff3u(i  ,j,k+1)*SLOPEXQ2,
     &                       diff3u(i+1,j,k+1)*SLOPEXQ3,
     &                       diff3u(i+1,j,k  )*SLOPEXQ4)
     &                           +max(
     &                       diff3v(i,j  ,k  )*SLOPEYQ1,
     &                       diff3v(i,j  ,k+1)*SLOPEYQ2,
     &                       diff3v(i,j+1,k+1)*SLOPEYQ3,
     &                       diff3v(i,j+1,k  )*SLOPEYQ4)
     &                          )*( max(
     &                diff3u(i  ,j,k  )*(pm(i  ,j)**2+SLOPEXQ1*cff**2),
     &                diff3u(i  ,j,k+1)*(pm(i  ,j)**2+SLOPEXQ2*cff**2),
     &                diff3u(i+1,j,k+1)*(pm(i+1,j)**2+SLOPEXQ3*cff**2),
     &                diff3u(i+1,j,k  )*(pm(i+1,j)**2+SLOPEXQ4*cff**2))
     &                             +max(
     &                diff3v(i,j  ,k  )*(pn(i,j  )**2+SLOPEYQ1*cff**2),
     &                diff3v(i,j  ,k+1)*(pn(i,j  )**2+SLOPEYQ2*cff**2),
     &                diff3v(i,j+1,k+1)*(pn(i,j+1)**2+SLOPEYQ3*cff**2),
     &                diff3v(i,j+1,k  )*(pn(i,j+1)**2+SLOPEYQ4*cff**2))
     &                            )
                FC(i,j,k)=-(sumX*wgt(idx)+sumE*wgt(ide))*FC(i,j,k)
          enddo
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          FC(i,j,0)  = 0.D0
          FC(i,j,N)  = 0.D0
          Akz (i,j,0 )= 0.0D0
          Akz (i,j,N )= 0.0D0
        enddo
      enddo
!$acc loop seq
      do k=1,N
          cff=0.5D0
          do j=jstr,jend
            do i=istr,iend+1
              FX(i,j)=-cff*diff3u(i,j,k)*(Hz(i,j,k)+Hz(i-1,j,k))
     &         *on_u(i,j)*(   dTdx(i  ,j,k) -
     &                0.25D0*dRdx(i,j,k)*( dTdr(i-1,j,k-1)+dTdr(i,j,k)
     &                                +  dTdr(i-1,j,k)+dTdr(i,j,k-1))
     &                                                              )
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              FE(i,j)=-cff*diff3v(i,j,k)*(Hz(i,j,k)+Hz(i,j-1,k))
     &        *om_v(i,j)*(  dTde(i,j  ,k) -
     &              0.25D0*dRde(i,j,k)*(   dTdr(i,j-1,k-1)+dTdr(i,j,k)
     &                              +    dTdr(i,j-1,k)+dTdr(i,j,k-1))
     &                                                              )
            enddo
          enddo
          do j=jstr,jend
            do i=istr,iend
              t(i,j,k,nnew,itrc)=Hz(i,j,k)*t(i,j,k,nnew,itrc)
     &         + dt*(
     &                   pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                    +FE(i,j+1)-FE(i,j))
     &                  +FC(i,j,k)-FC(i,j,k-1)    )
            enddo
          enddo
      enddo
       do j=Jstr,Jend
          indx=min(itrc,isalt)
          do i=istr,iend
            do k=1,N-1
              CD(i,k) = Akz(i,j,k)*(
     &           t(i,j,k+1,nstp,itrc)-t(i,j,k,nstp,itrc)
     &           )  / ( z_r(i,j,k+1)-z_r(i,j,k) )
            enddo
            CD(i,0) = 0.D0
            CD(i,N) = 0.D0
          enddo
          do i=istr,iend
            FFC(i,1)=dt*(Akt(i,j,1,indx)+Akz(i,j,1))
     &                               /( z_r(i,j,2)-z_r(i,j,1) )
            cff=1.D0/(Hz(i,j,1)+FFC(i,1))
            CF(i,1)= cff*FFC(i,1)
            DC(i,1)= cff*(t(i,j,1,nnew,itrc)-dt*(CD(i,1)-CD(i,0)))
          enddo
          do k=2,N-1,+1
            do i=istr,iend
              FFC(i,k)=dt*(Akt(i,j,k,indx)+Akz(i,j,k))
     &                              /( z_r(i,j,k+1)-z_r(i,j,k) )
              cff=1.D0/(Hz(i,j,k)+FFC(i,k)+FFC(i,k-1)*(1.D0-CF(i,k-1)))
              CF(i,k)=cff*FFC(i,k)
              DC(i,k)=cff*(t(i,j,k,nnew,itrc)+FFC(i,k-1)*DC(i,k-1)
     &                                    -dt*(CD(i,k)-CD(i,k-1))
     &                                                          )
            enddo
          enddo
          do i=istr,iend
             t(i,j,N,nnew,itrc)=( t(i,j,N,nnew,itrc)
     &                           -dt*(CD(i,N)-CD(i,N-1))
     &                           +FFC(i,N-1)*DC(i,N-1) )
     &                        /(Hz(i,j,N)+FFC(i,N-1)*(1.D0-CF(i,N-1)))
        enddo
          do k=N-1,1,-1
            do i=istr,iend
              t(i,j,k,nnew,itrc)=DC(i,k)+CF(i,k)*t(i,j,k+1,nnew,itrc)
            enddo
          enddo
      enddo
!$acc end kernels
        call exchange_r3d_tile (Istr,Iend,Jstr,Jend,
     &                          t(-1,-1,1,nnew,itrc))
      return
      end
