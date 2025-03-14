      subroutine u2dbc_tile(Istr,Iend,Jstr,Jend,grad)
      implicit none
      integer*4 Istr,Iend,Jstr,Jend, i,j
      real    grad(Istr-2:Iend+2,Jstr-2:Jend+2)
      real    eps,cff, cx,cy,u_str,
     &        dft,dfx,dfy, tau,tau_in,tau_out,hx,zx
      parameter (eps=1.D-20)
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
      real zeta(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      real ubar(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      real vbar(-1:Lm+2+padd_X,-1:Mm+2+padd_E,4)
      common /ocean_zeta/zeta
      common /ocean_ubar/ubar
      common /ocean_vbar/vbar
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
      real sustr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real svstr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /forces_sustr/sustr /forces_svstr/svstr
      real sustrg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real svstrg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /smsdat_sustrg/sustrg /smsdat_svstrg/svstrg
      real    sustrp(2), svstrp(2), sms_time(2)
      real    sms_cycle, sms_scale
      integer*4 itsms, sms_ncycle, sms_rec, lsusgrd
      integer*4 lsvsgrd,sms_tid, susid, svsid
      real    sms_origin_date_in_sec
      common /smsdat1/ sustrp, svstrp, sms_time
      common /smsdat2/ sms_origin_date_in_sec
      common /smsdat3/ sms_cycle, sms_scale
      common /smsdat4/ itsms, sms_ncycle, sms_rec, lsusgrd
      common /smsdat5/ lsvsgrd,sms_tid, susid, svsid
      integer*4 lwgrd, wid
      common /smsdat5/ lwgrd, wid
      real bustr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real bvstr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /forces_bustr/bustr /forces_bvstr/bvstr
      real bustrg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real bvstrg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /bmsdat_bustrg/bustrg /bmsdat_bvstrg/bvstrg
      real bms_tintrp(2), bustrp(2),    bvstrp(2), tbms(2)
      real bmsclen, bms_tstart, bms_tend,  tsbms, sclbms
      integer*4 itbms,      bmstid,busid, bvsid,     tbmsindx
      logical bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      common /bmsdat1/bms_tintrp, bustrp,       bvstrp,    tbms
      common /bmsdat2/bmsclen, bms_tstart, bms_tend, tsbms, sclbms
      common /bmsdat3/itbms,      bmstid,busid, bvsid,     tbmsindx
      common /bmsdat4/bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      real stflx(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /forces_stflx/stflx
      real shflx_rsw(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_rsw/shflx_rsw
      real shflx_rlw(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_rlw/shflx_rlw
      real shflx_lat(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_lat/shflx_lat
      real shflx_sen(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /frc_shflx_sen/shflx_sen
      real stflxg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2,NT)
      common /stfdat_stflxg/stflxg
      real stflxp(2,NT), stf_time(2,NT)
      real stf_cycle(NT), stf_scale(NT)
      integer*4 itstf(NT), stf_ncycle(NT), stf_rec(NT)
      integer*4 lstfgrd(NT), stf_tid(NT), stf_id(NT)
      REAL(kind=8) :: stf_origin_date_in_sec
      common /stfdat1/ stflxp,  stf_time, stf_cycle, stf_scale
      common /stfdat2/ stf_origin_date_in_sec
      common /stfdat3/ itstf, stf_ncycle, stf_rec, lstfgrd
      common /stfdat4/  stf_tid, stf_id
      real btflx(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /forces_btflx/btflx
      real tair(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real rhum(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real prate(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real radlw(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real radsw(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real wspd(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real uwnd(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real vwnd(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /bulk_tair/ tair
      common /bulk_rhum/ rhum
      common /bulk_prate/ prate
      common /bulk_radlw/ radlw
      common /bulk_radsw/ radsw
      common /bulk_wspd/ wspd
      common /bulk_uwnd/ uwnd
      common /bulk_vwnd/ vwnd
      real tairg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real rhumg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real prateg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real radlwg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real radswg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real uwndg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real vwndg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real wspdg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /bulkdat_tairg/tairg
      common /bulkdat_rhumg/rhumg
      common /bulkdat_prateg/prateg
      common /bulkdat_radlwg/radlwg
      common /bulkdat_radswg/radswg
      common /bulk_uwndg/uwndg
      common /bulk_vwndg/vwndg
      common /bulkdat_wspdg/wspdg
      real    tairp(2),rhump(2),pratep(2),radlwp(2),radswp(2)
      real    uwndp(2),vwndp(2)
      real    bulk_time(2), bulk_cycle
      integer*4 tair_id,rhum_id,prate_id,radlw_id,radsw_id
      integer*4 ltairgrd,lrhumgrd,lprategrd,lradlwgrd,lradswgrd
      REAL(kind=8) :: blk_origin_date_in_sec
      integer*4 uwnd_id,vwnd_id,luwndgrd,lvwndgrd
      integer*4 itbulk,bulk_ncycle,bulk_rec,bulk_tid
      integer*4 bulkunused
      common /bulkdat1_for/ tair_id,rhum_id,prate_id,radlw_id,radsw_id
      common /bulkdat1_grd
     &                 / ltairgrd,lrhumgrd,lprategrd,lradlwgrd,lradswgrd
      common /bulkdat1_tim/ itbulk, bulk_ncycle, bulk_rec, bulk_tid
      common /bulkdat1_uns/ bulkunused
      common /bulkdat1_wnd/ uwnd_id,vwnd_id,luwndgrd,lvwndgrd
      common /bulkdat2_for/ tairp,rhump,pratep,radlwp,radswp
      common /bulkdat2_tim/ bulk_time, bulk_cycle, 
     &                                            blk_origin_date_in_sec
      common /bulkdat2_wnd/ uwndp,vwndp
      real srflx(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /forces_srflx/srflx
      real srflxg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /srfdat_srflxg/srflxg
      real srflxp(2),srf_time(2)
      real srf_cycle, srf_scale
      integer*4 itsrf, srf_ncycle, srf_rec
      integer*4 lsrfgrd, srf_tid, srf_id
      REAL(kind=8) :: srf_origin_date_in_sec
      common /srfdat1/ srflxp, srf_time, srf_cycle, srf_scale
      common /srfdat2/ srf_origin_date_in_sec
      common /srfdat3/ itsrf,srf_ncycle,srf_rec,lsrfgrd,srf_tid,srf_id
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
      tau_in=dtfast*tauM_in
      tau_out=dtfast*tauM_out
!$acc kernels if(compute_on_device) default(present) async(1)
      grad = 1.D0
      if ((istr.eq.1 .and. .not.WEST_INTER) .and. (jstr.eq.1.and. 
     &                                           .not.SOUTH_INTER)) then
        grad(Istr,Jstr) = 0.5D0
      endif
      if ((iend.eq.Lmmpi .and. .not.EAST_INTER) .and. (jstr.eq.1.and. 
     &                                           .not.SOUTH_INTER)) then
        grad(Iend+1,Jstr) = 0.5D0
      endif
      if ((istr.eq.1 .and. .not.WEST_INTER) .and. (jend.eq.Mmmpi .and. 
     &                                           .not.NORTH_INTER)) then
        grad(Istr,Jend) = 0.5D0
      endif
      if ((iend.eq.Lmmpi .and. .not.EAST_INTER) .and. (jend.eq.Mmmpi 
     &                                     .and. .not.NORTH_INTER)) then
        grad(Iend+1,Jend) = 0.5D0
      endif
      if ((istr.eq.1 .and. .not.WEST_INTER)) then
        do j=Jstr,Jend
          cff=0.5D0*(h(Istr-1,j)+zeta(Istr-1,j,kstp)+
     &             h(Istr  ,j)+zeta(Istr  ,j,kstp))
          hx=sqrt(g/cff)
          cx=dtfast*cff*hx*0.5D0*(pm(Istr-1,j)+pm(Istr,j))
          zx=(0.5D0+cx)*zeta(istr,j,kstp)+(0.5D0-cx)*zeta(istr-1,j,kstp)
          if (cx .gt. 0.292893218813452D0) then
            zx=zx + ( zeta(istr,j,knew) +cx*zeta(istr-1,j,kstp)
     &                             -(1.D0+cx)*zeta(istr  ,j,kstp)
     &                           )*(1.D0-0.292893218813452D0/cx)**2
          endif
          ubar(Istr,j,knew)= 0.5D0*( (1.D0-cx)*ubar(Istr  ,j,kstp)
     &                                 +cx*ubar(Istr+1,j,kstp)
     &                                     +ubarbry_west(j)
     &                                     -hx*( zx
     &                                        -zetabry_west(j)
     &                           ))
          ubar(Istr,j,knew)=ubar(Istr,j,knew)*umask(Istr,j)
        enddo
      endif
      if ((iend.eq.Lmmpi .and. .not.EAST_INTER)) then
        do j=Jstr,Jend
          cff=0.5D0*(h(Iend  ,j)+zeta(Iend  ,j,kstp)+
     &             h(Iend+1,j)+zeta(Iend+1,j,kstp))
          hx=sqrt(g/cff)
          cx=dtfast*cff*hx*0.5D0*(pm(Iend,j)+pm(Iend+1,j))
          zx=(0.5D0+cx)*zeta(iend,j,kstp)+(0.5D0-cx)*zeta(iend+1,j,kstp)
          if (cx .gt. 0.292893218813452D0) then
            zx=zx + ( zeta(iend,j,knew) +cx*zeta(iend+1,j,kstp)
     &                             -(1.D0+cx)*zeta(iend  ,j,kstp)
     &                           )*(1.D0-0.292893218813452D0/cx)**2
          endif
          ubar(Iend+1,j,knew)=ubar(Iend+1,j,knew)*umask(Iend+1,j)
          ubar(Iend+1,j,knew)= 0.5D0*( (1.D0-cx)*ubar(Iend+1,j,kstp)
     &                                   +cx*ubar(Iend  ,j,kstp)
     &                                       +ubarbry_east(j)
     &                                       +hx*( zx
     &                                         -zetabry_east(j)
     &                             ))
        enddo
      endif
      if ((jstr.eq.1.and. .not.SOUTH_INTER)) then
        do i=IstrU,Iend
          cff=sqrt(0.5D0*g*(h(i-1,Jstr-1)+zeta(i-1,Jstr-1,kstp)+
     &                    h(i  ,Jstr-1)+zeta(i  ,Jstr-1,kstp)))
          cx=dtfast*0.5D0*cff*(pn(i-1,Jstr-1)+pn(i,Jstr-1))
          ubar(i,Jstr-1,knew)=(ubar(i,Jstr-1,kstp)
     &                     +cx*ubar(i,Jstr  ,knew) )/(1.D0+cx)
     &                         *umask(i,Jstr-1)
        enddo
      endif
      if ((jend.eq.Mmmpi .and. .not.NORTH_INTER)) then
        do i=IstrU,Iend
          cff=sqrt(0.5D0*g*(h(i-1,Jend+1)+zeta(i-1,Jend+1,kstp)+
     &                    h(i  ,Jend+1)+zeta(i  ,Jend+1,kstp)))
          cx=dtfast*0.5D0*cff*(pn(i-1,Jend+1)+pn(i,Jend+1))
          ubar(i,Jend+1,knew)=(ubar(i,Jend+1,kstp)
     &                     +cx*ubar(i,Jend  ,knew))/(1.D0+cx)
     &                        *umask(i,Jend+1)
        enddo
      endif
      if ((istr.eq.1 .and. .not.WEST_INTER) .and. (jstr.eq.1.and. 
     &                                           .not.SOUTH_INTER)) then
        ubar(Istr,Jstr-1,knew)=0.5D0*( ubar(Istr+1,Jstr-1,knew)
     &                                  +ubar(Istr,Jstr,knew))
     &                        *umask(Istr,Jstr-1)
      endif
      if ((iend.eq.Lmmpi .and. .not.EAST_INTER) .and. (jstr.eq.1.and. 
     &                                           .not.SOUTH_INTER)) then
        ubar(Iend+1,Jstr-1,knew)=0.5D0*( ubar(Iend,Jstr-1,knew)
     &                                +ubar(Iend+1,Jstr,knew))
     &                        *umask(Iend+1,Jstr-1)
      endif
      if ((istr.eq.1 .and. .not.WEST_INTER) .and. (jend.eq.Mmmpi .and. 
     &                                           .not.NORTH_INTER)) then
        ubar(Istr,Jend+1,knew)=0.5D0*( ubar(Istr+1,Jend+1,knew)
     &                                  +ubar(Istr,Jend,knew))
     &                        *umask(Istr,Jend+1)
      endif
      if ((iend.eq.Lmmpi .and. .not.EAST_INTER) .and. (jend.eq.Mmmpi 
     &                                     .and. .not.NORTH_INTER)) then
        ubar(Iend+1,Jend+1,knew)=0.5D0*( ubar(Iend,Jend+1,knew)
     &                                +ubar(Iend+1,Jend,knew))
     &                        *umask(Iend+1,Jend+1)
      endif
!$acc end kernels
      return
      end
