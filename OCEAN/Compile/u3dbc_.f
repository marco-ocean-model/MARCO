      subroutine u3dbc_tile (Istr,Iend,Jstr,Jend,grad)
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
      real rhoA(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real rhoS(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /coup_rhoA/rhoA           /coup_rhoS/rhoS
      real rufrc(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real rvfrc(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real rufrc_bak(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real rvfrc_bak(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /coup_rufrc/rufrc
      common /coup_rvfrc/rvfrc
      common /coup_rufrc_bak/rufrc_bak
      common /coup_rvfrc_bak/rvfrc_bak
      real Zt_avg1(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real DU_avg1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,5)
      real DV_avg1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,5)
      real DU_avg2(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real DV_avg2(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /ocean_Zt_avg1/Zt_avg1
      common /coup_DU_avg1/DU_avg1
      common /coup_DV_avg1/DV_avg1
      common /coup_DU_avg2/DU_avg2
      common /coup_DV_avg2/DV_avg2
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
      logical got_tbry(NT)
      common /bry_logical/ got_tbry
      real bry_time(2)
      REAL(kind=8) :: bry_origin_date_in_sec
      common /bry_indices_array/ bry_time,
     &        bry_origin_date_in_sec
      real bry_cycle
      common /bry_indices_real/ bry_cycle
      integer*4 bry_id, bry_time_id, bry_ncycle, bry_rec, itbry, ntbry
      common /bry_indices_integer/ bry_id, bry_time_id, bry_ncycle,
     &                             bry_rec, itbry, ntbry
      integer*4 zetabry_west_id
      common /zeta_west_id/ zetabry_west_id
      integer*4 ubarbry_west_id, vbarbry_west_id
      common /ubar_west_id/ ubarbry_west_id, vbarbry_west_id
      integer*4 ubry_west_id, vbry_west_id
      common /u_west_id/ ubry_west_id, vbry_west_id
      integer*4 tbry_west_id(NT)
      common /t_west_id/ tbry_west_id
      integer*4 zetabry_east_id
      common /zeta_east_id/ zetabry_east_id
      integer*4 ubarbry_east_id, vbarbry_east_id
      common /ubar_east_id/ ubarbry_east_id, vbarbry_east_id
      integer*4 ubry_east_id, vbry_east_id
      common /u_east_id/ ubry_east_id, vbry_east_id
      integer*4 tbry_east_id(NT)
      common /t_east_id/ tbry_east_id
      integer*4 zetabry_south_id
      common /zeta_south_id/ zetabry_south_id
      integer*4 ubarbry_south_id, vbarbry_south_id
      common /ubar_south_id/ ubarbry_south_id, vbarbry_south_id
      integer*4 ubry_south_id, vbry_south_id
      common /u_south_id/ ubry_south_id, vbry_south_id
      integer*4 tbry_south_id(NT)
      common /t_south_id/ tbry_south_id
      integer*4 zetabry_north_id
      common /zeta_north_id/ zetabry_north_id
      integer*4 ubarbry_north_id, vbarbry_north_id
      common /ubar_north_id/ ubarbry_north_id, vbarbry_north_id
      integer*4 ubry_north_id, vbry_north_id
      common /u_north_id/ ubry_north_id, vbry_north_id
      integer*4 tbry_north_id(NT)
      common /t_north_id/ tbry_north_id
      real zetabry_west(-1:Mm+2+padd_E),
     &    zetabry_west_dt(-1:Mm+2+padd_E,2)
      common /bry_zeta_west/ zetabry_west, zetabry_west_dt
      real ubarbry_west(-1:Mm+2+padd_E),
     &    ubarbry_west_dt(-1:Mm+2+padd_E,2)
     &    ,vbarbry_west(-1:Mm+2+padd_E),
     &    vbarbry_west_dt(-1:Mm+2+padd_E,2)
      common /bry_ubar_west/ ubarbry_west, ubarbry_west_dt,
     &                       vbarbry_west, vbarbry_west_dt
      real ubry_west(-1:Mm+2+padd_E,N),
     &    ubry_west_dt(-1:Mm+2+padd_E,N,2)
     &    ,vbry_west(-1:Mm+2+padd_E,N),
     &    vbry_west_dt(-1:Mm+2+padd_E,N,2)
      common /bry_u_west/ ubry_west, ubry_west_dt,
     &                    vbry_west, vbry_west_dt
      real tbry_west(-1:Mm+2+padd_E,N,NT),
     &     tbry_west_dt(-1:Mm+2+padd_E,N,2,NT)
      common /bry_t_west/ tbry_west, tbry_west_dt
      real zetabry_east(-1:Mm+2+padd_E),
     &    zetabry_east_dt(-1:Mm+2+padd_E,2)
      common /bry_zeta_east/ zetabry_east, zetabry_east_dt
      real ubarbry_east(-1:Mm+2+padd_E),
     &    ubarbry_east_dt(-1:Mm+2+padd_E,2)
     &    ,vbarbry_east(-1:Mm+2+padd_E),
     &    vbarbry_east_dt(-1:Mm+2+padd_E,2)
      common /bry_ubar_east/ ubarbry_east, ubarbry_east_dt,
     &                       vbarbry_east, vbarbry_east_dt
      real ubry_east(-1:Mm+2+padd_E,N),
     &    ubry_east_dt(-1:Mm+2+padd_E,N,2)
     &    ,vbry_east(-1:Mm+2+padd_E,N),
     &    vbry_east_dt(-1:Mm+2+padd_E,N,2)
      common /bry_u_east/ ubry_east, ubry_east_dt,
     &                    vbry_east, vbry_east_dt
      real tbry_east(-1:Mm+2+padd_E,N,NT),
     &    tbry_east_dt(-1:Mm+2+padd_E,N,2,NT)
      common /bry_t_east/ tbry_east, tbry_east_dt
      real zetabry_south(-1:Lm+2+padd_X),
     &    zetabry_south_dt(-1:Lm+2+padd_X,2)
      common /bry_zeta_south/ zetabry_south, zetabry_south_dt
      real ubarbry_south(-1:Lm+2+padd_X),
     &    ubarbry_south_dt(-1:Lm+2+padd_X,2)
     &    ,vbarbry_south(-1:Lm+2+padd_X),
     &    vbarbry_south_dt(-1:Lm+2+padd_X,2)
      common /bry_ubar_south/ ubarbry_south, ubarbry_south_dt,
     &                        vbarbry_south, vbarbry_south_dt
      real ubry_south(-1:Lm+2+padd_X,N),
     &    ubry_south_dt(-1:Lm+2+padd_X,N,2)
     &    ,vbry_south(-1:Lm+2+padd_X,N),
     &    vbry_south_dt(-1:Lm+2+padd_X,N,2)
      common /bry_u_south/ ubry_south, ubry_south_dt,
     &                     vbry_south, vbry_south_dt
      real tbry_south(-1:Lm+2+padd_X,N,NT),
     &    tbry_south_dt(-1:Lm+2+padd_X,N,2,NT)
      common /bry_t_south/ tbry_south, tbry_south_dt
      real zetabry_north(-1:Lm+2+padd_X),
     &    zetabry_north_dt(-1:Lm+2+padd_X,2)
      common /bry_zeta_north/ zetabry_north, zetabry_north_dt
      real ubarbry_north(-1:Lm+2+padd_X),
     &    ubarbry_north_dt(-1:Lm+2+padd_X,2)
     &    ,vbarbry_north(-1:Lm+2+padd_X),
     &    vbarbry_north_dt(-1:Lm+2+padd_X,2)
      common /bry_ubar_north/ ubarbry_north, ubarbry_north_dt,
     &                        vbarbry_north, vbarbry_north_dt
      real ubry_north(-1:Lm+2+padd_X,N),
     &    ubry_north_dt(-1:Lm+2+padd_X,N,2)
     &    ,vbry_north(-1:Lm+2+padd_X,N),
     &    vbry_north_dt(-1:Lm+2+padd_X,N,2)
      common /bry_u_north/ ubry_north, ubry_north_dt,
     &                     vbry_north, vbry_north_dt
      real tbry_north(-1:Lm+2+padd_X,N,NT),
     &    tbry_north_dt(-1:Lm+2+padd_X,N,2,NT)
      common /bry_t_north/ tbry_north, tbry_north_dt
      integer*4 Istr,Iend,Jstr,Jend, i,j,k
      real    grad(Istr-2:Iend+2,Jstr-2:Jend+2)
      real    cff,eps,
     &        cx,cy, dft,dfx,dfy, tau,tau_in,tau_out
      parameter (eps=1.D-20)
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
      tau_in=dt*tauM_in
      tau_out=dt*tauM_out
!$acc kernels if(compute_on_device) default(present)
      if ((istr.eq.1 .and. .not.WEST_INTER)) then
        do k=1,N
          do j=Jstr,Jend+1
            grad(Istr  ,j)=(u(Istr  ,j,k,nstp)-u(Istr  ,j-1,k,nstp))
     &                                                *pmask(Istr,j)
            grad(Istr+1,j)=(u(Istr+1,j,k,nstp)-u(Istr+1,j-1,k,nstp))
     &                                              *pmask(Istr+1,j)
          enddo
          do j=Jstr,Jend
            dft=u(Istr+1,j,k,nstp)-u(Istr+1,j,k,nnew)
            dfx=u(Istr+1,j,k,nnew)-u(Istr+2,j,k,nnew)
            if (dfx*dft .lt. 0.D0) then
              dft=0.D0
              tau=tau_in
            else
              tau=tau_out
            endif
            if (dft*(grad(Istr+1,j)+grad(Istr+1,j+1)) .gt. 0.D0) then
              dfy=grad(Istr+1,j)
            else
              dfy=grad(Istr+1,j+1)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
            cy=min(cff,max(dft*dfy,-cff))
            u(Istr,j,k,nnew)=(  cff*u(Istr,j,k,nstp)
     &                        +cx*u(Istr+1,j,k,nnew)
     &                    -max(cy,0.D0)*grad(Istr,j  )
     &                    -min(cy,0.D0)*grad(Istr,j+1)
     &                                    )/(cff+cx)
            u(Istr,j,k,nnew)=(1.D0-tau)*u(Istr,j,k,nnew)+tau*
     &                                      ubry_west(j,k)
            u(Istr,j,k,nnew)=u(Istr,j,k,nnew)*umask(Istr,j)
          enddo
        enddo
      endif
      if ((iend.eq.Lmmpi .and. .not.EAST_INTER)) then
        do k=1,N
          do j=Jstr,Jend+1
            grad(Iend  ,j)=(u(Iend  ,j,k,nstp)-u(Iend  ,j-1,k,nstp))
     &                                                *pmask(Iend,j)
            grad(Iend+1,j)=(u(Iend+1,j,k,nstp)-u(Iend+1,j-1,k,nstp))
     &                                              *pmask(Iend+1,j)
          enddo
          do j=Jstr,Jend
            dft=u(Iend,j,k,nstp)-u(Iend  ,j,k,nnew)
            dfx=u(Iend,j,k,nnew)-u(Iend-1,j,k,nnew)
            if (dfx*dft .lt. 0.D0) then
              dft=0.D0
              tau=tau_in
            else
              tau=tau_out
            endif
            if (dft*(grad(Iend,j)+grad(Iend,j+1)) .gt. 0.D0) then
              dfy=grad(Iend,j)
            else
              dfy=grad(Iend,j+1)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
            cy=min(cff,max(dft*dfy,-cff))
            u(Iend+1,j,k,nnew)=(  cff*u(Iend+1,j,k,nstp)
     &                              +cx*u(Iend,j,k,nnew)
     &                      -max(cy,0.D0)*grad(Iend+1,j  )
     &                      -min(cy,0.D0)*grad(Iend+1,j+1)
     &                                        )/(cff+cx)
            u(Iend+1,j,k,nnew)=(1.D0-tau)*u(Iend+1,j,k,nnew)+tau*
     &                                          ubry_east(j,k)
            u(Iend+1,j,k,nnew)=u(Iend+1,j,k,nnew)*umask(Iend+1,j)
          enddo
        enddo
      endif
      if ((jstr.eq.1.and. .not.SOUTH_INTER)) then
        do k=1,N
          do i=IstrU-1,Iend
            grad(i,Jstr-1)=u(i+1,Jstr-1,k,nstp)-u(i,Jstr-1,k,nstp)
            grad(i,Jstr  )=u(i+1,Jstr  ,k,nstp)-u(i,Jstr  ,k,nstp)
          enddo
          do i=IstrU,Iend
            dft=u(i,Jstr,k,nstp)-u(i,Jstr  ,k,nnew)
            dfx=u(i,Jstr,k,nnew)-u(i,Jstr+1,k,nnew)
            if (dfx*dft .lt. 0.D0) then
              dft=0.D0
              tau=tau_in
            else
              tau=tau_out
            endif
            if (dft*(grad(i-1,Jstr)+grad(i,Jstr)) .gt. 0.D0) then
              dfy=grad(i-1,Jstr)
            else
              dfy=grad(i  ,Jstr)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
            cy=min(cff,max(dft*dfy,-cff))
            u(i,Jstr-1,k,nnew)=(  cff*u(i,Jstr-1,k,nstp)
     &                              +cx*u(i,Jstr,k,nnew)
     &                      -max(cy,0.D0)*grad(i-1,Jstr-1)
     &                      -min(cy,0.D0)*grad(i  ,Jstr-1)
     &                                        )/(cff+cx)
           u(i,Jstr-1,k,nnew)=(1.D0-tau)*u(i,Jstr-1,k,nnew)
     &                               +tau*ubry_south(i,k)
            u(i,Jstr-1,k,nnew)=u(i,Jstr-1,k,nnew)*umask(i,Jstr-1)
          enddo
        enddo
      endif
      if ((jend.eq.Mmmpi .and. .not.NORTH_INTER)) then
        do k=1,N
          do i=IstrU-1,Iend
            grad(i,Jend  )=u(i+1,Jend  ,k,nstp)-u(i,Jend  ,k,nstp)
            grad(i,Jend+1)=u(i+1,Jend+1,k,nstp)-u(i,Jend+1,k,nstp)
          enddo
          do i=IstrU,Iend
            dft=u(i,Jend,k,nstp)-u(i,Jend  ,k,nnew)
            dfx=u(i,Jend,k,nnew)-u(i,Jend-1,k,nnew)
            if (dfx*dft .lt. 0.D0) then
              dft=0.D0
              tau=tau_in
            else
              tau=tau_out
            endif
            if (dft*(grad(i-1,Jend)+grad(i,Jend)) .gt. 0.D0) then
              dfy=grad(i-1,Jend)
            else
              dfy=grad(i  ,Jend)
            endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
            cy=min(cff,max(dft*dfy,-cff))
            u(i,Jend+1,k,nnew)=(  cff*u(i,Jend+1,k,nstp)
     &                              +cx*u(i,Jend,k,nnew)
     &                      -max(cy,0.D0)*grad(i-1,Jend+1)
     &                      -min(cy,0.D0)*grad(i  ,Jend+1)
     &                                        )/(cff+cx)
            u(i,Jend+1,k,nnew)=(1.D0-tau)*u(i,Jend+1,k,nnew)
     &                                +tau*ubry_north(i,k)
            u(i,Jend+1,k,nnew)=u(i,Jend+1,k,nnew)*umask(i,Jend+1)
          enddo
        enddo
      endif
      if ((istr.eq.1 .and. .not.WEST_INTER) .and. (jstr.eq.1.and. 
     &                                           .not.SOUTH_INTER)) then
        do k=1,N
          u(Istr,Jstr-1,k,nnew)=0.5D0*( u(Istr+1,Jstr-1,k,nnew)
     &                               +u(Istr  ,Jstr  ,k,nnew))
     &                          *umask(Istr,Jstr-1)
        enddo
      endif
      if ((iend.eq.Lmmpi .and. .not.EAST_INTER) .and. (jstr.eq.1.and. 
     &                                           .not.SOUTH_INTER)) then
        do k=1,N
          u(Iend+1,Jstr-1,k,nnew)=0.5D0*( u(Iend,Jstr-1,k,nnew)
     &                                 +u(Iend+1,Jstr,k,nnew))
     &                            *umask(Iend+1,Jstr-1)
        enddo
      endif
      if ((istr.eq.1 .and. .not.WEST_INTER) .and. (jend.eq.Mmmpi .and. 
     &                                           .not.NORTH_INTER)) then
        do k=1,N
          u(Istr,Jend+1,k,nnew)=0.5D0*( u(Istr+1,Jend+1,k,nnew)
     &                               +u(Istr  ,Jend  ,k,nnew))
     &                          *umask(Istr,Jend+1)
        enddo
      endif
      if ((iend.eq.Lmmpi .and. .not.EAST_INTER) .and. (jend.eq.Mmmpi 
     &                                     .and. .not.NORTH_INTER)) then
        do k=1,N
          u(Iend+1,Jend+1,k,nnew)=0.5D0*( u(Iend,Jend+1,k,nnew)
     &                                 +u(Iend+1,Jend,k,nnew))
     &                            *umask(Iend+1,Jend+1)
        enddo
      endif
!$acc end kernels
      return
      end
