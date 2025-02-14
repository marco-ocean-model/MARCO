










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









MODULE oce_trc


   use scalars
   use ncscrum
   USE par_pisces
   USE in_out_manager
   USE lib_mpp

   IMPLICIT NONE
   PUBLIC

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
! This is include file "grid.h": Environmental two-dimensional
! arrays associated with curvilinear horizontal coordinate system.
!
! h       Model topography (bottom depth [m] at RHO-points.)
! dh      Topograhy increment in case of moving bathymetry
! f       Coriolis parameter [1/s].
! fomn    Compound term, f/[pm*pn] at RHO points.
!
! angler  Angle [radians] between XI-axis and the direction
!             to the EAST at RHO-points.
!
! latr    Latitude (degrees_north) at RHO-, U-, and V-points.
! latu
! latv
! lonr    Longitude (degrees_east) at RHO-, U-, and V-points.
! lonu
! lonv
!
! xp      XI-coordinates [m] at PSI-points.
! xr      XI-coordinates (m] at RHO-points.
! yp      ETA-coordinates [m] at PSI-points.
! yr      ETA-coordinates [m] at RHO-points.
!
! pm      Coordinate transformation metric "m" [1/meters]
!              associated with the differential distances in XI.
! pn      Coordinate transformation metric "n" [1/meters]
!               associated with the differential distances in ETA.
! om_u    Grid spacing [meters] in the XI -direction at U-points.
! om_v    Grid spacing [meters] in the XI -direction at V-points.
! on_u    Grid spacing [meters] in the ETA-direction at U-points.
! on_v    Grid spacing [meters] in the ETA-direction at V-points.
!
! dmde    ETA-derivative of inverse metric factor "m", d(1/M)/d(ETA).
! dndx     XI-derivative  of inverse metric factor "n", d(1/N)/d(XI).
!
! pmon_p  Compound term, pm/pn at PSI-points.
! pmon_r  Compound term, pm/pn at RHO-points.
! pmon_u  Compound term, pm/pn at U-points.
! pnom_p  Compound term, pn/pm at PSI-points.
! pnom_r  Compound term, pn/pm at RHO-points.
! pnom_v  Compound term, pn/pm at V-points.
!
! rmask   Land-sea masking arrays at RHO-,U-,V- and PSI-points.
! umask   (rmask,umask,vmask) = (0=Land, 1=Sea);
! vmask
! pmask    pmask=(0=Land, 1=Sea, 1-gamma2 =boundary).
!
! reducu  reduction coefficient along x-axis for rivers sections
! reducv  reduction coefficient along y-axis for rivers sections

      real h(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real hinv(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real f(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real fomn(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /grid_h/h /grid_hinv/hinv /grid_f/f /grid_fomn/fomn

      real angler(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /grid_angler/angler

      real latr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real lonr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real latu(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real lonu(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real latv(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real lonv(-2:Lm+3+padd_X,-2:Mm+3+padd_E)

      common /grid_latr/latr /grid_lonr/lonr
      common /grid_latu/latu /grid_lonu/lonu
      common /grid_latv/latv /grid_lonv/lonv

      real pm(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pn(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pn_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pm_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pm_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pn_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /metrics_pm/pm    /metrics_pn/pn
      common /metrics_omr/om_r /metrics_on_r/on_r
      common /metrics_omu/om_u /metrics_on_u/on_u
      common /metrics_omv/om_v /metrics_on_v/on_v
      common /metrics_omp/om_p /metrics_on_p/on_p
      common /metrics_pnu/pn_u /metrics_pmv/pm_v
      common /metrics_pmu/pm_u /metrics_pnv/pn_v

      real dmde(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real dndx(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /metrics_dmde/dmde    /metrics_dndx/dndx

      real pmon_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pmon_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pmon_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pnom_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pnom_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pnom_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real grdscl(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /metrics_pmon_p/pmon_p /metrics_pnom_p/pnom_p
      common /metrics_pmon_r/pmon_r /metrics_pnom_r/pnom_r
      common /metrics_pmon_u/pmon_u /metrics_pnom_v/pnom_v
      common /metrics_grdscl/grdscl

      real rmask(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pmask(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real umask(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real vmask(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pmask2(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /mask_r/rmask
      common /mask_p/pmask
      common /mask_u/umask
      common /mask_v/vmask
      common /mask_p2/pmask2



      real zob(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /Z0B_VAR/zob

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

      real u(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3)
      real v(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3)
      real t(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t

      real Hz(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real Hz_bak(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real z_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real z_w(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Huon(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real Hvom(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /grid_Hz_bak/Hz_bak /grid_zw/z_w /grid_Huon/Huon
      common /grid_Hvom/Hvom

      real We(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /grid_Hz/Hz /grid_zr/z_r /grid_We/We

      real wz(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,3)
      common /ocean_wz/wz


      real rho1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real rho(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /ocean_rho1/rho1 /ocean_rho/rho
      real qp1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /ocean_qp1/qp1
      real qp2
      parameter (qp2=0.0000172)




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
!  This is include file "forces.h"
!--------------------------------------------------------------------
!  SURFACE MOMENTUM FLUX (WIND STRESS):
!--------------------------------------------------------------------
!  sustr |  XI- and ETA-components of kinematic surface momentum flux
!  svstr |  (wind stresses) defined at horizontal U- and V-points.
!            dimensioned as [m^2/s^2].
!
      real sustr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real svstr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /forces_sustr/sustr /forces_svstr/svstr

!
!  tsms      Time of surface momentum stresses.
!
!  sustrg |  Two-time level gridded data for XI- and ETA-componets
!  svstrg |  of kinematic surface momentum flux (wind stess).
!
!  sustrp |  Two-time level point data for XI- and ETA-componets
!  svstrp |  of kinematic surface momentum flux (wind stess).
!
      real sustrg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      real svstrg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      common /smsdat_sustrg/sustrg /smsdat_svstrg/svstrg

      real    sustrp(2), svstrp(2), sms_time(2)
      real    sms_cycle, sms_scale
      integer itsms, sms_ncycle, sms_rec, lsusgrd
      integer lsvsgrd,sms_tid, susid, svsid
      real    sms_origin_date_in_sec
      common /smsdat1/ sustrp, svstrp, sms_time
      common /smsdat2/ sms_origin_date_in_sec
      common /smsdat3/ sms_cycle, sms_scale
      common /smsdat4/ itsms, sms_ncycle, sms_rec, lsusgrd
      common /smsdat5/ lsvsgrd,sms_tid, susid, svsid

      integer lwgrd, wid
      common /smsdat5/ lwgrd, wid

!
!  BOTTOM MOMENTUM FLUX:
!--------------------------------------------------------------------
!  bustr |  XI- and ETA-components of kinematic bottom momentum flux
!  bvstr |  (drag) defined at horizontal U- and V-points [m^2/s^2].
      real bustr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real bvstr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /forces_bustr/bustr /forces_bvstr/bvstr
!
!  tbms      Time of surface momentum stresses.
!
!  bustrg |  Two-time level gridded data for XI- and ETA-componets
!  bvstrg |  of kinematic bottom momentum flux.
!
!  bustrp |  Two-time level point data for XI- and ETA-componets
!  bvstrp |  of kinematic bottom momentum flux.
!
      real bustrg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      real bvstrg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      common /bmsdat_bustrg/bustrg /bmsdat_bvstrg/bvstrg

      real bms_tintrp(2), bustrp(2),    bvstrp(2), tbms(2)
      real bmsclen, bms_tstart, bms_tend,  tsbms, sclbms
      integer itbms,      bmstid,busid, bvsid,     tbmsindx
      logical bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      common /bmsdat1/bms_tintrp, bustrp,       bvstrp,    tbms
      common /bmsdat2/bmsclen, bms_tstart, bms_tend, tsbms, sclbms
      common /bmsdat3/itbms,      bmstid,busid, bvsid,     tbmsindx
      common /bmsdat4/bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd

!
!  SURFACE TRACER FLUXES:
!--------------------------------------------------------------------
!  stflx   Kinematic surface fluxes of tracer type variables at
!          horizontal RHO-points. Physical dimensions [degC m/s] -
!          temperature; [PSU m/s] - salinity.
!
      real stflx(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /forces_stflx/stflx
      real shflx_rsw(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /frc_shflx_rsw/shflx_rsw
      real shflx_rlw(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /frc_shflx_rlw/shflx_rlw
      real shflx_lat(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /frc_shflx_lat/shflx_lat
      real shflx_sen(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /frc_shflx_sen/shflx_sen
!
!  stflxg   Two-time level surface tracer flux grided data.
!  stflxp   Two-time level surface tracer flux point  data.
!  tstflx   Time of surface tracer flux.
!
      real stflxg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2,NT)
      common /stfdat_stflxg/stflxg

      real stflxp(2,NT), stf_time(2,NT)
      real stf_cycle(NT), stf_scale(NT)
      integer itstf(NT), stf_ncycle(NT), stf_rec(NT)
      integer lstfgrd(NT), stf_tid(NT), stf_id(NT)
      REAL(kind=8) :: stf_origin_date_in_sec
      common /stfdat1/ stflxp,  stf_time, stf_cycle, stf_scale
      common /stfdat2/ stf_origin_date_in_sec
      common /stfdat3/ itstf, stf_ncycle, stf_rec, lstfgrd
      common /stfdat4/  stf_tid, stf_id
!
!  BOTTOM TRACER FLUXES:
!--------------------------------------------------------------------
!  btflx  Kinematic bottom fluxes of tracer type variables at
!         horizontal RHO-points. Physical dimensions [degC m/s] -
!         temperature; [PSU m/s] - salinity.
!
      real btflx(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /forces_btflx/btflx



!
!
!
!  HEAT FLUX BULK FORMULATION
!--------------------------------------------------------------------
!  tair     surface air temperature at 2m [degree Celsius].
!  wspd     wind speed at 10m [m s-1].
!  rhum     surface air relative humidity 2m [fraction]
!  prate    surface precipitation rate [cm day-1]
!  radlw    net terrestrial longwave radiation [Watts meter-2]
!  radsw    net solar shortwave radiation [Watts meter-2]
!  patm2d   atmospheric pressure above mean seal level
!  paref     reference pressure to compute inverse barometer effect
      real tair(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real rhum(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real prate(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real radlw(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real radsw(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real wspd(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real uwnd(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real vwnd(-2:Lm+3+padd_X,-2:Mm+3+padd_E)

      common /bulk_tair/ tair
      common /bulk_rhum/ rhum
      common /bulk_prate/ prate
      common /bulk_radlw/ radlw
      common /bulk_radsw/ radsw
      common /bulk_wspd/ wspd
      common /bulk_uwnd/ uwnd
      common /bulk_vwnd/ vwnd

      real tairg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      real rhumg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      real prateg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      real radlwg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      real radswg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      real uwndg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      real vwndg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      real wspdg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)

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
      integer tair_id,rhum_id,prate_id,radlw_id,radsw_id
      integer ltairgrd,lrhumgrd,lprategrd,lradlwgrd,lradswgrd
      REAL(kind=8) :: blk_origin_date_in_sec
      integer uwnd_id,vwnd_id,luwndgrd,lvwndgrd
      integer itbulk,bulk_ncycle,bulk_rec,bulk_tid
      integer bulkunused

      common /bulkdat1_for/ tair_id,rhum_id,prate_id,radlw_id,radsw_id
      common /bulkdat1_grd/ ltairgrd,lrhumgrd,lprategrd,lradlwgrd,lradswgrd
      common /bulkdat1_tim/ itbulk, bulk_ncycle, bulk_rec, bulk_tid
      common /bulkdat1_uns/ bulkunused
      common /bulkdat1_wnd/ uwnd_id,vwnd_id,luwndgrd,lvwndgrd

      common /bulkdat2_for/ tairp,rhump,pratep,radlwp,radswp
      common /bulkdat2_tim/ bulk_time, bulk_cycle, blk_origin_date_in_sec
      common /bulkdat2_wnd/ uwndp,vwndp
!
!  SOLAR SHORT WAVE RADIATION FLUX.
!--------------------------------------------------------------------
!  srflx  Kinematic surface shortwave solar radiation flux
!         [degC m/s] at horizontal RHO-points
!
      real srflx(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /forces_srflx/srflx
!
!  srflxg | Two-time-level grided and point data for surface
!  srflxp |      solar shortwave radiation flux grided data.
!  tsrflx   Time of solar shortwave radiation flux.
!
      real srflxg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      common /srfdat_srflxg/srflxg

      real srflxp(2),srf_time(2)
      real srf_cycle, srf_scale
      integer itsrf, srf_ncycle, srf_rec
      integer lsrfgrd, srf_tid, srf_id
      REAL(kind=8) :: srf_origin_date_in_sec
      common /srfdat1/ srflxp, srf_time, srf_cycle, srf_scale
      common /srfdat2/ srf_origin_date_in_sec
      common /srfdat3/ itsrf,srf_ncycle,srf_rec,lsrfgrd,srf_tid,srf_id




!--------------------------------------------------------------------
!  WIND INDUCED WAVES: everything is defined at rho-point
!--------------------------------------------------------------------
! wfrq | BBL/MRL | wind-induced wave frequency [rad/s]
! uorb | BBL     | xi-component  of wave-induced bed orbital velocity [m/s]
! vorb | BBL     | eta-component of wave-induced bed orbital velocity [m/s]
! wdrx | MRL     | cosine of wave direction [non dimension]
! wdre | MRL     | sine of   wave direction [non dimension]
! whrm | MRL     | (RMS) wave height (twice the wave amplitude) [m]
! wepb | MRL     | breaking dissipation rate (\epsilon_b term) [m3/s3]
! wepd | MRL     | frictional dissipation rate (\epsilon_d term) [m3/s3]
! wlm  | MRL     | mean length wave from input data (coupling or forcing)
! wepr | ROLLER  | roller dissipation rate (\epsilon_r term) [m3/s3]
! wbst | MRL/BKPP| frictional dissipation stress (e_d k/sigma) [m2/s2]
!--------------------------------------------------------------------





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
! This is include file "mixing.h"
!  ==== == ======= ==== ==========
!
      real visc2_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real visc2_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real visc2_sponge_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real visc2_sponge_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /mixing_visc2_r/visc2_r /mixing_visc2_p/visc2_p
      common /mixing_visc2_sponge_r/visc2_sponge_r
      common /mixing_visc2_sponge_p/visc2_sponge_p
      real diff2_sponge(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real diff2(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /mixing_diff2_sponge/diff2_sponge
      common /mixing_diff2/diff2
      real diff4_sponge(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real diff4(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /mixing_diff4_sponge/diff4_sponge
      common /mixing_diff4/diff4
      real diff3d_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real diff3d_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /mixing_diff3d_u/diff3d_u
      common /mixing_diff3d_v/diff3d_v
      real dRdx(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real dRde(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real idRz(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /mixing_dRdx/dRdx
      common /mixing_dRde/dRde
      common /mixing_idRz/idRz
      real Rslope_max,Gslope_max
      parameter (Gslope_max=1., Rslope_max=0.05)
      integer ismooth
      real csmooth
      common /mixing_csmooth/ csmooth
      common /mixing_ismooth/ ismooth

      real Akv(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Akt(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,2)
      common /mixing_Akv/Akv /mixing_Akt/Akt
      real Akv_old(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Akt_old(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /mixing_Akvold/Akv_old /mixing_Aktold/Akt_old

      real bvf(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /mixing_bvf/ bvf

      real hel(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /lmd_hel/hel

!
! Generic Length Scale
!
      real trb(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,2,NGLS)
      common /gls_trb/trb
      real Lscale(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /gls_lsc/Lscale
      real Eps_gls(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /gls_eps/Eps_gls
      integer kbl(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /gls_kbl/ kbl
      real hbl(-2:Lm+3+padd_X,-2:Mm+3+padd_E  )
      common /gls_hbl/ hbl
      real cm0
      common /gls_cm0/ cm0



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
! This is include file "diagnostics.h": tracer equation terms
! for output purposes:
!
!
!
      real bioFlux(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,NumFluxTerms)
      real bioVSink(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NumVSinkTerms)
      real bioFlux_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,NumFluxTerms)
      real bioVSink_avg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NumVSinkTerms)
      real timediabio_avg
      common /diag_bioFlux/bioFlux
      common /diag_bioVSink/bioVSink
      common /diag_bioFlux_avg/bioFlux_avg
      common /diag_bioVSink_avg/bioVSink_avg
      common /diag_timediabio_avg/timediabio_avg


   !! * Substitutions































   PUBLIC   trc_oce_rgb        ! routine called by traqsr.F90
   PUBLIC   trc_oce_rgb_read   ! routine called by traqsr.F90
   PUBLIC   trc_oce_ext_lev    ! function called by traqsr.F90 at least
   PUBLIC   trc_oce_alloc      ! function called by nemogcm.F90
   PUBLIC   ocean_2_pisces 

   PUBLIC   glob_sum
   PUBLIC   tracer_stat
   INTERFACE glob_sum
      MODULE PROCEDURE  glob_sum_2d, glob_sum_3d
   END INTERFACE   

   LOGICAL , PUBLIC ::   l_co2cpl  = .false.   !: atmospheric pco2 recieved from oasis
   LOGICAL , PUBLIC ::   l_offline = .false.   !: offline passive tracers flag
   REAL(wp), PUBLIC ::   r_si2                 !: largest depth of extinction (blue & 0.01 mg.m-3)  (RGB)
   LOGICAL , PUBLIC ::   ln_trcdc2dm           !: Diurnal cycle for TOP
   INTEGER, PUBLIC  ::   nksr    !: =nkV, i.e. maximum level of light extinction
   !
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:,:) ::   tra   !: traceur concentration for next time step
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:)   ::   tmask
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:)   ::   etot3     !: light absortion coefficient
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)     ::   oce_co2   !: ocean carbon flux
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)     ::   atm_co2   !: ocean carbon flux
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)     ::   fr_i      !: ocean carbon flux
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)     ::   e1e2t
   REAL(wp) , PUBLIC, DIMENSION(3,61)   ::   rkrgb    ! tabulated attenuation coefficients for RGB absorption


  REAL, DIMENSION(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N+1,3) :: gdepw   ! W-depht
  REAL, DIMENSION(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N+1,3) :: e3w     ! W-vertical scale factor

  INTEGER, PUBLIC :: Istrp,Iendp,Jstrp,Jendp
!$OMP threadprivate(Istrp,Iendp)
!$OMP threadprivate(Jstrp,Jendp)

  LOGICAL :: ln_ctl     = .false.
  LOGICAL :: ln_qsr_bio = .false.
  REAL(wp)    :: rn_si0     = 0.35
!  INTEGER :: stdout
  INTEGER :: jpdom_data = 1
  INTEGER :: jpdom_global = 1
  CHARACTER(len = 8)  ::   cn_cfg = "BENGUELA"      ! current name of the NetCDF mask file acting as a key
   

CONTAINS

   SUBROUTINE ocean_2_pisces( kt, Istr, Iend, Jstr, Jend )
   
      INTEGER, intent(in) :: kt, Istr, Iend, Jstr, Jend
      INTEGER :: mi, mj, jk, jl

      IF( kt == ntstart ) THEN
         Istrp = Istr
         Iendp = Iend
         Jstrp = Jstr
         Jendp = Jend
      ENDIF

      DO jl = 1, 3
         DO jk = 1, N+1
            DO mj =  Jstrp, Jendp 
               DO mi =  Istrp, Iendp 
                  gdepw(mi,mj,N+2-jk,jl) = -(z_w(mi,mj,jk-1)-z_w(mi,mj,N))
               END DO
            END DO
         END DO

         DO jk = 2, N
            DO mj =  Jstrp, Jendp 
               DO mi =  Istrp, Iendp 
                  e3w(mi,mj,jk,jl) = -z_r(mi,mj,N+1-jk) + z_r(mi,mj,N+2-jk)
              END DO
            END DO
         END DO

         DO mj =  Jstrp, Jendp 
            DO mi =  Istrp, Iendp 
               e3w(mi,mj,1  ,jl) = -2 * z_r(mi,mj,N)
               e3w(mi,mj,N+1,jl) = 2 * ( -z_w(mi,mj,0) + z_r(mi,mj,1) )
            END DO
         END DO
         !
      ENDDO
      !
   END SUBROUTINE ocean_2_pisces

   INTEGER FUNCTION trc_oce_alloc()
      !!----------------------------------------------------------------------
      !!                  ***  trc_oce_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( tra(Istrp:Iendp,Jstrp:Jendp,N,jptra),tmask(Istrp:Iendp,Jstrp:Jendp,N), e1e2t(Istrp:Iendp,Jstrp:Jendp), &
         &      etot3(Istrp:Iendp,Jstrp:Jendp,N), oce_co2(Istrp:Iendp,Jstrp:Jendp), &
         &     atm_co2(Istrp:Iendp,Jstrp:Jendp), fr_i(Istrp:Iendp,Jstrp:Jendp),STAT=trc_oce_alloc )

      IF( trc_oce_alloc /= 0 )   CALL ctl_warn('trc_oce_alloc: failed to allocate etot3 array')
      !
   END FUNCTION trc_oce_alloc


   SUBROUTINE trc_oce_rgb( prgb )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_opt_init  ***
      !!
      !! ** Purpose :   Initialization of of the optical scheme
      !!
      !! ** Method  :   Set a look up table for the optical coefficients
      !!                i.e. the attenuation coefficient for R-G-B light
      !!                tabulated in Chlorophyll class (from JM Andre)
      !!
      !! ** Action  :   prgb(3,61) tabulated R-G-B attenuation coef.
      !!
      !! Reference  : Lengaigne et al. 2007, Clim. Dyn., V28, 5, 503-516.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(3,61), INTENT(out) ::   prgb   ! tabulated attenuation coefficient
      !
      INTEGER  ::   jc     ! dummy loop indice
      INTEGER  ::   irgb   ! temporary integer
      REAL(wp) ::   zchl   ! temporary scalar
      REAL(wp), DIMENSION(4,61) ::   zrgb   ! tabulated attenuation coefficient (formerly read in 'kRGB61.txt')
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN
         WRITE(stdout,*)
         WRITE(stdout,*) '   trc_oce_rgb : Initialisation of the optical look-up table'
         WRITE(stdout,*) '   ~~~~~~~~~~~ '
      ENDIF
      !
      !  Chlorophyll        !     Blue attenuation     !     Green attenuation    !     Red attenuation      !
      zrgb(1, 1) =  0.010   ;   zrgb(2, 1) = 0.01618   ;   zrgb(3, 1) = 0.07464   ;   zrgb(4, 1) = 0.37807
      zrgb(1, 2) =  0.011   ;   zrgb(2, 2) = 0.01654   ;   zrgb(3, 2) = 0.07480   ;   zrgb(4, 2) = 0.37823
      zrgb(1, 3) =  0.013   ;   zrgb(2, 3) = 0.01693   ;   zrgb(3, 3) = 0.07499   ;   zrgb(4, 3) = 0.37840
      zrgb(1, 4) =  0.014   ;   zrgb(2, 4) = 0.01736   ;   zrgb(3, 4) = 0.07518   ;   zrgb(4, 4) = 0.37859
      zrgb(1, 5) =  0.016   ;   zrgb(2, 5) = 0.01782   ;   zrgb(3, 5) = 0.07539   ;   zrgb(4, 5) = 0.37879
      zrgb(1, 6) =  0.018   ;   zrgb(2, 6) = 0.01831   ;   zrgb(3, 6) = 0.07562   ;   zrgb(4, 6) = 0.37900
      zrgb(1, 7) =  0.020   ;   zrgb(2, 7) = 0.01885   ;   zrgb(3, 7) = 0.07586   ;   zrgb(4, 7) = 0.37923
      zrgb(1, 8) =  0.022   ;   zrgb(2, 8) = 0.01943   ;   zrgb(3, 8) = 0.07613   ;   zrgb(4, 8) = 0.37948
      zrgb(1, 9) =  0.025   ;   zrgb(2, 9) = 0.02005   ;   zrgb(3, 9) = 0.07641   ;   zrgb(4, 9) = 0.37976
      zrgb(1,10) =  0.028   ;   zrgb(2,10) = 0.02073   ;   zrgb(3,10) = 0.07672   ;   zrgb(4,10) = 0.38005
      zrgb(1,11) =  0.032   ;   zrgb(2,11) = 0.02146   ;   zrgb(3,11) = 0.07705   ;   zrgb(4,11) = 0.38036
      zrgb(1,12) =  0.035   ;   zrgb(2,12) = 0.02224   ;   zrgb(3,12) = 0.07741   ;   zrgb(4,12) = 0.38070
      zrgb(1,13) =  0.040   ;   zrgb(2,13) = 0.02310   ;   zrgb(3,13) = 0.07780   ;   zrgb(4,13) = 0.38107
      zrgb(1,14) =  0.045   ;   zrgb(2,14) = 0.02402   ;   zrgb(3,14) = 0.07821   ;   zrgb(4,14) = 0.38146
      zrgb(1,15) =  0.050   ;   zrgb(2,15) = 0.02501   ;   zrgb(3,15) = 0.07866   ;   zrgb(4,15) = 0.38189
      zrgb(1,16) =  0.056   ;   zrgb(2,16) = 0.02608   ;   zrgb(3,16) = 0.07914   ;   zrgb(4,16) = 0.38235
      zrgb(1,17) =  0.063   ;   zrgb(2,17) = 0.02724   ;   zrgb(3,17) = 0.07967   ;   zrgb(4,17) = 0.38285
      zrgb(1,18) =  0.071   ;   zrgb(2,18) = 0.02849   ;   zrgb(3,18) = 0.08023   ;   zrgb(4,18) = 0.38338
      zrgb(1,19) =  0.079   ;   zrgb(2,19) = 0.02984   ;   zrgb(3,19) = 0.08083   ;   zrgb(4,19) = 0.38396
      zrgb(1,20) =  0.089   ;   zrgb(2,20) = 0.03131   ;   zrgb(3,20) = 0.08149   ;   zrgb(4,20) = 0.38458
      zrgb(1,21) =  0.100   ;   zrgb(2,21) = 0.03288   ;   zrgb(3,21) = 0.08219   ;   zrgb(4,21) = 0.38526
      zrgb(1,22) =  0.112   ;   zrgb(2,22) = 0.03459   ;   zrgb(3,22) = 0.08295   ;   zrgb(4,22) = 0.38598
      zrgb(1,23) =  0.126   ;   zrgb(2,23) = 0.03643   ;   zrgb(3,23) = 0.08377   ;   zrgb(4,23) = 0.38676
      zrgb(1,24) =  0.141   ;   zrgb(2,24) = 0.03842   ;   zrgb(3,24) = 0.08466   ;   zrgb(4,24) = 0.38761
      zrgb(1,25) =  0.158   ;   zrgb(2,25) = 0.04057   ;   zrgb(3,25) = 0.08561   ;   zrgb(4,25) = 0.38852
      zrgb(1,26) =  0.178   ;   zrgb(2,26) = 0.04289   ;   zrgb(3,26) = 0.08664   ;   zrgb(4,26) = 0.38950
      zrgb(1,27) =  0.200   ;   zrgb(2,27) = 0.04540   ;   zrgb(3,27) = 0.08775   ;   zrgb(4,27) = 0.39056
      zrgb(1,28) =  0.224   ;   zrgb(2,28) = 0.04811   ;   zrgb(3,28) = 0.08894   ;   zrgb(4,28) = 0.39171
      zrgb(1,29) =  0.251   ;   zrgb(2,29) = 0.05103   ;   zrgb(3,29) = 0.09023   ;   zrgb(4,29) = 0.39294
      zrgb(1,30) =  0.282   ;   zrgb(2,30) = 0.05420   ;   zrgb(3,30) = 0.09162   ;   zrgb(4,30) = 0.39428
      zrgb(1,31) =  0.316   ;   zrgb(2,31) = 0.05761   ;   zrgb(3,31) = 0.09312   ;   zrgb(4,31) = 0.39572
      zrgb(1,32) =  0.355   ;   zrgb(2,32) = 0.06130   ;   zrgb(3,32) = 0.09474   ;   zrgb(4,32) = 0.39727
      zrgb(1,33) =  0.398   ;   zrgb(2,33) = 0.06529   ;   zrgb(3,33) = 0.09649   ;   zrgb(4,33) = 0.39894
      zrgb(1,34) =  0.447   ;   zrgb(2,34) = 0.06959   ;   zrgb(3,34) = 0.09837   ;   zrgb(4,34) = 0.40075
      zrgb(1,35) =  0.501   ;   zrgb(2,35) = 0.07424   ;   zrgb(3,35) = 0.10040   ;   zrgb(4,35) = 0.40270
      zrgb(1,36) =  0.562   ;   zrgb(2,36) = 0.07927   ;   zrgb(3,36) = 0.10259   ;   zrgb(4,36) = 0.40480
      zrgb(1,37) =  0.631   ;   zrgb(2,37) = 0.08470   ;   zrgb(3,37) = 0.10495   ;   zrgb(4,37) = 0.40707
      zrgb(1,38) =  0.708   ;   zrgb(2,38) = 0.09056   ;   zrgb(3,38) = 0.10749   ;   zrgb(4,38) = 0.40952
      zrgb(1,39) =  0.794   ;   zrgb(2,39) = 0.09690   ;   zrgb(3,39) = 0.11024   ;   zrgb(4,39) = 0.41216
      zrgb(1,40) =  0.891   ;   zrgb(2,40) = 0.10374   ;   zrgb(3,40) = 0.11320   ;   zrgb(4,40) = 0.41502
      zrgb(1,41) =  1.000   ;   zrgb(2,41) = 0.11114   ;   zrgb(3,41) = 0.11639   ;   zrgb(4,41) = 0.41809
      zrgb(1,42) =  1.122   ;   zrgb(2,42) = 0.11912   ;   zrgb(3,42) = 0.11984   ;   zrgb(4,42) = 0.42142
      zrgb(1,43) =  1.259   ;   zrgb(2,43) = 0.12775   ;   zrgb(3,43) = 0.12356   ;   zrgb(4,43) = 0.42500
      zrgb(1,44) =  1.413   ;   zrgb(2,44) = 0.13707   ;   zrgb(3,44) = 0.12757   ;   zrgb(4,44) = 0.42887
      zrgb(1,45) =  1.585   ;   zrgb(2,45) = 0.14715   ;   zrgb(3,45) = 0.13189   ;   zrgb(4,45) = 0.43304
      zrgb(1,46) =  1.778   ;   zrgb(2,46) = 0.15803   ;   zrgb(3,46) = 0.13655   ;   zrgb(4,46) = 0.43754
      zrgb(1,47) =  1.995   ;   zrgb(2,47) = 0.16978   ;   zrgb(3,47) = 0.14158   ;   zrgb(4,47) = 0.44240
      zrgb(1,48) =  2.239   ;   zrgb(2,48) = 0.18248   ;   zrgb(3,48) = 0.14701   ;   zrgb(4,48) = 0.44765
      zrgb(1,49) =  2.512   ;   zrgb(2,49) = 0.19620   ;   zrgb(3,49) = 0.15286   ;   zrgb(4,49) = 0.45331
      zrgb(1,50) =  2.818   ;   zrgb(2,50) = 0.21102   ;   zrgb(3,50) = 0.15918   ;   zrgb(4,50) = 0.45942
      zrgb(1,51) =  3.162   ;   zrgb(2,51) = 0.22703   ;   zrgb(3,51) = 0.16599   ;   zrgb(4,51) = 0.46601
      zrgb(1,52) =  3.548   ;   zrgb(2,52) = 0.24433   ;   zrgb(3,52) = 0.17334   ;   zrgb(4,52) = 0.47313
      zrgb(1,53) =  3.981   ;   zrgb(2,53) = 0.26301   ;   zrgb(3,53) = 0.18126   ;   zrgb(4,53) = 0.48080
      zrgb(1,54) =  4.467   ;   zrgb(2,54) = 0.28320   ;   zrgb(3,54) = 0.18981   ;   zrgb(4,54) = 0.48909
      zrgb(1,55) =  5.012   ;   zrgb(2,55) = 0.30502   ;   zrgb(3,55) = 0.19903   ;   zrgb(4,55) = 0.49803
      zrgb(1,56) =  5.623   ;   zrgb(2,56) = 0.32858   ;   zrgb(3,56) = 0.20898   ;   zrgb(4,56) = 0.50768
      zrgb(1,57) =  6.310   ;   zrgb(2,57) = 0.35404   ;   zrgb(3,57) = 0.21971   ;   zrgb(4,57) = 0.51810
      zrgb(1,58) =  7.079   ;   zrgb(2,58) = 0.38154   ;   zrgb(3,58) = 0.23129   ;   zrgb(4,58) = 0.52934
      zrgb(1,59) =  7.943   ;   zrgb(2,59) = 0.41125   ;   zrgb(3,59) = 0.24378   ;   zrgb(4,59) = 0.54147
      zrgb(1,60) =  8.912   ;   zrgb(2,60) = 0.44336   ;   zrgb(3,60) = 0.25725   ;   zrgb(4,60) = 0.55457
      zrgb(1,61) = 10.000   ;   zrgb(2,61) = 0.47804   ;   zrgb(3,61) = 0.27178   ;   zrgb(4,61) = 0.56870
      !
      prgb(:,:) = zrgb(2:4,:)
      !
      r_si2 = 1.e0 / zrgb(2, 1)        ! blue with the smallest chlorophyll concentration)
      IF(mynode .eq. 0) WRITE(stdout,*) '      RGB longest depth of extinction    r_si2 = ', r_si2
      !
      DO jc = 1, 61                         ! check
         zchl = zrgb(1,jc)
         irgb = NINT( 41 + 20.* LOG10( zchl ) + 1.e-15 )
         IF( irgb /= jc ) THEN
            IF(mynode .eq. 0) WRITE(stdout,*) '    jc =', jc, '  Chl = ', zchl, '  Chl class = ', irgb
            CALL ctl_stop( 'STOP', 'trc_oce_rgb : inconsistency in Chl tabulated attenuation coeff.' )
         ENDIF
      END DO
      !
   END SUBROUTINE trc_oce_rgb


   SUBROUTINE trc_oce_rgb_read( prgb )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_opt_init  ***
      !!
      !! ** Purpose :   Initialization of of the optical scheme
      !!
      !! ** Method  :   read the look up table for the optical coefficients
      !!
      !! ** input   :   xkrgb(61) precomputed array corresponding to the
      !!                          attenuation coefficient (from JM Andre)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(3,61), INTENT(out) ::   prgb   ! tabulated attenuation coefficient
      !
      INTEGER  ::   jc, jb ! dummy loop indice
      INTEGER  ::   irgb   ! temporary integer
      REAL(wp) ::   zchl   ! temporary scalar
      INTEGER  ::   numlight
      !!----------------------------------------------------------------------
      !
      IF(mynode .eq. 0) THEN                         ! control print
         WRITE(stdout,*)
         WRITE(stdout,*) ' trc_oce_rgb_read : optical look-up table read in kRGB61.txt file'
         WRITE(stdout,*) ' ~~~~~~~~~~~~~~~~'
         WRITE(stdout,*)
      ENDIF
      !
      CALL ctl_opn( numlight, 'kRGB61.txt', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, stdout, mynode .eq. 0 )
      DO jc = 1, 61
         READ(numlight,*) zchl, ( prgb(jb,jc), jb = 1, 3 )
         irgb = NINT( 41 + 20.* LOG10( zchl ) + 1.e-15 )
         IF(mynode .eq. 0) WRITE(stdout,*) '    jc =', jc, '  Chl = ', zchl, '  irgb = ', irgb
         IF( irgb /= jc ) THEN
            IF(mynode .eq. 0) WRITE(stdout,*) '    jc =', jc, '  Chl = ', zchl, '  Chl class = ', irgb
            CALL ctl_stop( 'STOP','trc_oce_rgb_read : inconsistency in Chl tabulated attenuation coeff.' )
         ENDIF
      END DO
      CLOSE( numlight )
      !
      r_si2 = 1.e0 / prgb(1, 1)      ! blue with the smallest chlorophyll concentration)
      IF(mynode .eq. 0) WRITE(stdout,*) '      RGB longest depth of extinction    r_si2 = ', r_si2
      !
   END SUBROUTINE trc_oce_rgb_read


   FUNCTION trc_oce_ext_lev( prldex, pqsr_frc ) RESULT( pjl )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE trc_oce_ext_lev  ***
      !!
      !! ** Purpose :   compute max. level for light penetration
      !!
      !! ** Method  :   the function provides the level at which irradiance
      !!                becomes negligible (i.e. = 1.e-15 W/m2) for 3 or 2 bands light
      !!                penetration: I(z) = pqsr_frc * EXP(hext/prldex) = 1.e-15 W/m2
      !!                # prldex is the longest depth of extinction:
      !!                   - prldex = 23 m (2 bands case)
      !!                   - prldex = 62 m (3 bands case: blue waveband & 0.01 mg/m2 for the chlorophyll)
      !!                # pqsr_frc is the fraction of solar radiation which penetrates,
      !!                considering Qsr=240 W/m2 and rn_abs = 0.58:
      !!                   - pqsr_frc = Qsr * (1-rn_abs)   = 1.00e2 W/m2 (2 bands case)
      !!                   - pqsr_frc = Qsr * (1-rn_abs)/3 = 0.33e2 W/m2 (3 bands case & equi-partition)
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::   prldex    ! longest depth of extinction
      REAL(wp), INTENT(in) ::   pqsr_frc  ! frac. solar radiation which penetrates
      !
      INTEGER  ::   jk, pjl            ! levels
      REAL(wp) ::   zhext              ! deepest level till which light penetrates
      REAL(wp) ::   zprec = 15._wp     ! precision to reach -LOG10(1.e-15)
      REAL(wp) ::   zem                ! temporary scalar
      !!----------------------------------------------------------------------
      !
      ! It is not necessary to compute anything below the following depth
      zhext = prldex * ( LOG(10._wp) * zprec + LOG(pqsr_frc) )
      !
      ! Level of light extinction
      pjl = N
      DO jk = N, 1, -1
         IF(SUM(tmask(Istrp:Iendp,Jstrp:Jendp,jk)) > 0 ) THEN
!            zem = MAXVAL( gdepw_1d(jk+1) * tmask(:,:,jk) )
            zem = MAXVAL( gdepw(Istrp:Iendp,Jstrp:Jendp,jk+1,1) * tmask(Istrp:Iendp,Jstrp:Jendp,jk) )
            IF( zem >= zhext )   pjl = jk                       ! last T-level reached by Qsr
         ELSE
            pjl = jk                                            ! or regional sea-bed depth
         ENDIF
      END DO
      !
   END FUNCTION trc_oce_ext_lev

      FUNCTION glob_sum_2d( cdname, ptab )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_sum_2d ***
      !!
      !! ** Purpose : perform a sum in calling DDPDD routine
      !!----------------------------------------------------------------------
      CHARACTER(len=*),  INTENT(in   ) ::   cdname  ! name of the calling subroutine
      REAL(wp), INTENT(in), DIMENSION(Istrp:Iendp,Jstrp:Jendp) ::   ptab
      REAL(wp)                             ::   glob_sum_2d   ! global masked sum
      !!
      REAL(wp)   ::   ztmp
      INTEGER    ::   mi, mj   ! dummy loop indices
      !!-----------------------------------------------------------------------
      !
      ztmp = 0.e0
      DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ztmp =  ptab(mi,mj) * rmask(mi,mj)
      END DO   ;   END DO
      IF( .true. )   CALL mpp_sum( ztmp )   ! sum over the global domain
      glob_sum_2d = REAL(ztmp,wp)
      !
   END FUNCTION glob_sum_2d

      FUNCTION glob_sum_3d( cdname, ptab )
      !!----------------------------------------------------------------------
      !!                  ***  FUNCTION  glob_sum_3d ***
      !!
      !! ** Purpose : perform a sum on a 3D array in calling DDPDD routine
      !!----------------------------------------------------------------------
      CHARACTER(len=*),  INTENT(in   ) ::   cdname  ! name of the calling subroutine
      REAL(wp), INTENT(in), DIMENSION(Istrp:Iendp,Jstrp:Jendp,N) ::   ptab
      REAL(wp)                               ::   glob_sum_3d   ! global masked sum
      !!
      REAL(wp)   ::   ztmp
      INTEGER    ::   mi, mj, jk   ! dummy loop indices
      !!-----------------------------------------------------------------------
      !
      !
      ztmp = 0.e0
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
         ztmp =  ptab(mi,mj,jk) * rmask(mi,mj)
      END DO   ;   END DO   ;   END DO
      IF( .true. )   CALL mpp_sum( ztmp )   ! sum over the global domain
      glob_sum_3d = REAL(ztmp,wp)
      !
   END FUNCTION glob_sum_3d

   SUBROUTINE tracer_stat( kt, clname )
      !!----------------------------------------------------------------------
      !!                    ***  trc_rst_stat  ***
      !!
      !! ** purpose  :   Compute tracers statistics
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)  :: kt
      CHARACTER(len=20), DIMENSION(jptra), INTENT(in) ::   clname
      !
      INTEGER  :: mi, mj, jk, jn
      REAL(wp) :: ztra, zmin, zmax, zmean, areatot, zcoef
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,1:N,jptra)  :: ptra
      REAL(wp), DIMENSION(Istrp:Iendp,Jstrp:Jendp,1:N)        :: zmask, zvol
      !!----------------------------------------------------------------------

      IF( mynode .eq. 0 ) THEN
         WRITE(stdout,*)
         WRITE(stdout,*) ' TRACER STAT at time-step kt = ', kt
         WRITE(stdout,*)
      ENDIF
      !
! to have coherent units when calling tracer_stat
      IF( kt .eq. ntstart ) THEN
        zcoef = 1.e-6
      ELSE
        zcoef = 1.
      ENDIF

      DO jn = 1, jptra
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            ptra(mi,mj,jk,jn) = t(mi,mj,N+1-jk,nnew,itemp+ntrc_salt+jn) * zcoef
         END DO   ;   END DO   ;   END DO
      ENDDO
      areatot = 0.                                                           ! total volume of the ocean
      DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
          zvol(mi,mj,jk)  = om_r(mi,mj) * on_r(mi,mj) * Hz(mi,mj,N+1-jk) * tmask(mi,mj,jk)
          zmask(mi,mj,jk) = tmask(mi,mj,jk) * rmask(mi,mj)
          areatot         = areatot + zvol(mi,mj,jk)
     END DO   ;   END DO   ;   END DO
     IF( .true. )   CALL mpp_sum( areatot )     ! sum over the global domain

     DO jn = 1, jptra
         ztra = 0.
         DO jk= 1, N ; DO mj=Jstrp,Jendp ; DO mi=Istrp,Iendp
            ztra  = ztra + ptra(mi,mj,jk,jn) * zvol(mi,mj,jk)
         END DO   ;   END DO   ;   END DO
         zmin  = MINVAL( ptra(:,:,:,jn), mask= ( zmask(:,:,:) /= 0. ) )
         zmax  = MAXVAL( ptra(:,:,:,jn), mask= ( zmask(:,:,:) /= 0. ) )
         IF( .true. ) THEN
            CALL mpp_sum( ztra )      ! min over the global domain
            CALL mpp_min( zmin )      ! min over the global domain
            CALL mpp_max( zmax )      ! max over the global domain
         END IF
         zmean  = ztra / areatot
         IF(mynode .eq. 0) WRITE(stdout,9000) jn, TRIM( clname(jn) ), zmean, zmin, zmax
      END DO
      WRITE(stdout,*)
9000  FORMAT(' tracer nb :',i2,'    name :',a10,'    mean :',e18.10,'    min :',e18.10, '    max :',e18.10 )
      !
   END SUBROUTINE tracer_stat

END MODULE oce_trc

