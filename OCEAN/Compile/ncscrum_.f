      module ncscrum
      use param
      integer*4 filetype_his, filetype_avg
     &       ,filetype_dia, filetype_dia_avg
     &       ,filetype_diaM, filetype_diaM_avg
     &       ,filetype_diags_vrt, filetype_diags_vrt_avg
     &       ,filetype_diags_ek, filetype_diags_ek_avg
     &       ,filetype_diags_pv, filetype_diags_pv_avg
     &       ,filetype_diags_eddy_avg
     &       ,filetype_surf, filetype_surf_avg
     &       ,filetype_diabio, filetype_diabio_avg
     &       ,filetype_abl, filetype_abl_avg
      parameter (filetype_his=1, filetype_avg=2,
     &           filetype_dia=3, filetype_dia_avg=4,
     &           filetype_diaM=5, filetype_diaM_avg=6,
     &           filetype_diags_vrt=7, filetype_diags_vrt_avg=8,
     &           filetype_diags_ek=9, filetype_diags_ek_avg=10,
     &           filetype_diags_pv=11, filetype_diags_pv_avg=12,
     &           filetype_diags_eddy_avg=17,
     &           filetype_surf=13, filetype_surf_avg=14,
     &           filetype_diabio=15,filetype_diabio_avg=16,
     &           filetype_abl=18, filetype_abl_avg=19)
      integer*4 iloop, indextemp
      integer*4 indxTime, indxZ, indxUb, indxVb
      parameter (indxTime=1, indxZ=2, indxUb=3, indxVb=4)
      integer*4 indxU, indxV
      parameter (indxU=6, indxV=7)
      integer*4 indxT
      parameter (indxT=indxV+1)
      integer*4 indxS
      parameter (indxS=indxV+ntrc_temp+1)
      integer*4 ncidpisrst, nrecpisrst
      integer*4 ncidsedrst, nrecsedrst
      integer*4 indxDIC, indxTAL, indxOXY, indxCAL, indxPO4,
     &        indxPOC, indxSIL, indxPHY, indxZOO, indxDOC,
     &        indxDIA, indxMES, indxDSI, indxFER, indxBFE,
     &        indxGOC, indxSFE, indxDFE, indxGSI, indxNFE,
     &        indxNCH, indxDCH, indxNO3, indxNH4
      parameter (indxDIC =indxV+ntrc_temp+ntrc_salt+ntrc_pas+1,
     &           indxTAL =indxDIC+1, indxOXY=indxDIC+2)
      parameter (indxPOC=indxDIC+3, indxPHY =indxDIC+4,
     &           indxZOO=indxDIC+5, indxDOC =indxDIC+6,
     &           indxNO3=indxDIC+7, indxFER =indxDIC+8)
      integer*4 indxLGW
      parameter (indxLGW=indxDIC+24)
      integer*4 indxDON, indxDOP, indxPON, indxPOP, indxNPH,
     &        indxPPH, indxNDI, indxPDI, indxPIC, indxNPI,
     &        indxPPI, indxPFE, indxPCH, indxGON, indxGOP
      parameter (indxDON=indxDIC+25, indxDOP=indxDIC+26,
     &           indxPON=indxDIC+27, indxPOP=indxDIC+28,
     &           indxNPH=indxDIC+29, indxPPH=indxDIC+30,
     &           indxNDI=indxDIC+31, indxPDI=indxDIC+32,
     &           indxPIC=indxDIC+33, indxNPI=indxDIC+34,
     &           indxPPI=indxDIC+35, indxPFE=indxDIC+36,
     &           indxPCH=indxDIC+37, indxGON=indxDIC+38,
     &           indxGOP=indxDIC+39)
      integer*4 indxBSD, indxBSS
      parameter (indxBSD=indxV+ntrc_temp+ntrc_salt+ntrc_pas+ntrc_bio+1,
     &           indxBSS=101)
      integer*4 indxbioFlux, indxbioVSink
      parameter (indxbioFlux=indxV+ntrc_temp+ntrc_salt
     &                           +ntrc_pas+ntrc_bio+ntrc_sed
     &                           +ntrc_diats+ntrc_diauv+ntrc_diavrt
     &                           +ntrc_diaek+ntrc_diapv+ntrc_diaeddy
     &                                               +ntrc_surf+400)
      parameter (indxbioVSink=indxbioFlux+NumFluxTerms)
      integer*4 indxO, indxW, indxR, indxVisc, indxDiff, indxAkv, 
     &                                                           indxAkt
      parameter (indxO=indxV+ntrc_temp+ntrc_salt+ntrc_pas+ntrc_bio
     &                      +ntrc_sed+ntrc_substot
     &           +ntrc_diats+ntrc_diauv+ntrc_diavrt+ntrc_diaek
     &           +ntrc_diapv+ntrc_diaeddy+ntrc_surf+ntrc_diabio+1,
     &           indxW=indxO+1, indxR=indxO+2, indxVisc=indxO+3,
     &           indxDiff=indxO+4,indxAkv=indxO+5, indxAkt=indxO+6)
      integer*4 indxAks
      parameter (indxAks=indxAkv+ntrc_temp+4)
      integer*4 indxHbl
      parameter (indxHbl=indxAkv+ntrc_temp+5)
      integer*4 indxTke
      parameter (indxTke=indxAkv+ntrc_temp+7)
      integer*4 indxGls
      parameter (indxGls=indxAkv+ntrc_temp+8)
      integer*4 indxLsc
      parameter (indxLsc=indxAkv+ntrc_temp+9)
      integer*4 indxAkk
      parameter (indxAkk=indxAkv+ntrc_temp+10)
      integer*4 indxAkp
      parameter (indxAkp=indxAkv+ntrc_temp+11)
      integer*4 indxSSH
      parameter (indxSSH=indxAkv+ntrc_temp+12)
      integer*4 indxbvf
      parameter (indxbvf=indxSSH+1)
      integer*4 indxrufrc
      parameter (indxrufrc=indxSSH+300)
      integer*4 indxrvfrc
      parameter (indxrvfrc=indxrufrc+1)
      integer*4 indxSUSTR, indxSVSTR
      integer*4 indxdRdx,indxdRde
      integer*4 indxru_nbq,indxrv_nbq
      integer*4 indxru_nbq_avg2,indxrv_nbq_avg2
      integer*4 indxqdmu_nbq,indxqdmv_nbq
      parameter (indxru_nbq=indxrvfrc+1,
     &           indxrv_nbq=indxru_nbq+1,
     &           indxru_nbq_avg2=indxrv_nbq+1,
     &           indxrv_nbq_avg2=indxru_nbq_avg2+1,
     &           indxqdmu_nbq=indxrv_nbq_avg2+1,
     &           indxqdmv_nbq=indxqdmu_nbq+1)
      parameter (indxdRdx=indxqdmv_nbq+1,
     &           indxdRde=indxdRdx+1)
      parameter (indxSUSTR=indxdRde+1,
     &           indxSVSTR=indxdRde+2)
      integer*4 indxTime2
      parameter (indxTime2=indxSSH+4)
      integer*4 indxShflx, indxShflx_rsw
      parameter (indxShflx=indxSUSTR+5)
      integer*4 indxSwflx
      parameter (indxSwflx=indxShflx+1, indxShflx_rsw=indxShflx+2)
      integer*4 indxSST, indxdQdSST
      parameter (indxSST=indxShflx_rsw+1, indxdQdSST=indxShflx_rsw+2)
      integer*4 indxWSPD,indxTAIR,indxRHUM,indxRADLW,indxRADSW,
     &        indxPRATE,indxUWND,indxVWND,indxPATM
      parameter (indxWSPD=indxSST+3,  indxTAIR=indxSST+4,
     &           indxRHUM=indxSST+5,  indxRADLW=indxSST+6,
     &           indxRADSW=indxSST+7, indxPRATE=indxSST+8,
     &           indxUWND=indxSST+9,  indxVWND=indxSST+10,
     &           indxPATM=indxSST+11)
      integer*4 indxShflx_rlw,indxShflx_lat,indxShflx_sen
      parameter (indxShflx_rlw=indxSST+12,
     &           indxShflx_lat=indxSST+13, indxShflx_sen=indxSST+14)
      integer*4 indxWstr
      parameter (indxWstr=indxSUSTR+23)
      integer*4 indxUWstr
      parameter (indxUWstr=indxSUSTR+24)
      integer*4 indxVWstr
      parameter (indxVWstr=indxSUSTR+25)
      integer*4 indxBostr
      parameter (indxBostr=indxSUSTR+26)
      integer*4 indxBustr, indxBvstr
      parameter (indxBustr=indxSUSTR+27,  indxBvstr=indxBustr+1)
      integer*4 indxWWA,indxWWD,indxWWP,indxWEB,indxWED,indxWER
      parameter (indxWWA=indxSUSTR+42, indxWWD=indxWWA+1,
     &           indxWWP=indxWWA+2
     &                             )
      integer*4 r2dvar, u2dvar, v2dvar, p2dvar, r3dvar,
     &                u3dvar, v3dvar, p3dvar, w3dvar,
     &                pw3dvar, b3dvar
      parameter (r2dvar=0, u2dvar=1, v2dvar=2, p2dvar=3,
     & r3dvar=4, u3dvar=5, v3dvar=6, p3dvar=7, w3dvar=8,
     & pw3dvar=11, b3dvar=12)
      integer*4 xi_rho,xi_u, eta_rho,eta_v
      integer*4 ncidfrc, ncidbulk, ncidclm,  ntsms
     &     , ncidqbar, ncidbtf
     &     , ntsrf,  ntssh,  ntsst, ntsss, ntuclm
     &     , ntbulk, ntqbar, ntww
      integer*4 nttclm(NT), ntstf(NT), nttsrc(NT)
     &       , ntbtf(NT)
      integer*4 ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
      integer*4 rstT(NT)
      integer*4 rstAkv,rstAkt
      integer*4 rstAks
      integer*4 rstTke,rstGls
      integer*4 rstBustr, rstBvstr
      integer*4 rstrufrc,rstrvfrc
      integer*4 rstru_nbq,rstrv_nbq
      integer*4 rstru_nbq_avg2,rstrv_nbq_avg2
      integer*4 rstqdmu_nbq,rstqdmv_nbq
      integer*4 rstdRdx,rstdRde
      integer*4  ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisBustr, hisBvstr
     &      , hisShflx, hisSwflx, hisShflx_rsw, hisBhflx, hisBwflx
     &      , hisUnbq, hisVnbq, hisWnbq, hisRnbq, hisCnbq
     &      , hisU,   hisV,   hisR,    hisHbl, hisHbbl
     &      , hisO,   hisW,   hisVisc, hisDiff
     &      , hisAkv, hisAkt, hisAks
     &      , hisbvf
     &      , hisTke, hisGls, hisLsc
     &      , hisShflx_rlw
     &      , hisShflx_lat,   hisShflx_sen
     &      , hisHel
      integer*4 hisT(NT)
      integer*4 nciddiabio, nrecdiabio, nrpfdiabio
     &      , diaTimebio, diaTime2bio, diaTstepbio
     &      , diabioFlux(NumFluxTerms)
     &      , diabioVSink(NumVSinkTerms)
     &      , diabioGasExc(NumGasExcTerms)
      integer*4 ncidavg, nrecavg,  nrpfavg
     &      , avgTime, avgTime2, avgTstep, avgZ, avgUb,  avgVb
     &      , avgBostr, avgWstr, avgUwstr, avgVwstr
     &      , avgBustr, avgBvstr
     &      , avgShflx, avgSwflx, avgShflx_rsw, avgBhflx, avgBwflx
     &      , avgU,   avgV,   avgR,    avgHbl, avgHbbl
     &      , avgO,   avgW,   avgVisc, avgDiff
     &      , avgAkv, avgAkt, avgAks
     &      , avgbvf
     &      , avgTke, avgGls, avgLsc
     &      , avgHel
      integer*4 avgT(NT)
      integer*4 avgShflx_rlw
     &      , avgShflx_lat,   avgShflx_sen
      integer*4 nciddiabio_avg, nrecdiabio_avg, nrpfdiabio_avg
     &      , diaTimebio_avg, diaTime2bio_avg, diaTstepbio_avg
     &      , diabioFlux_avg(NumFluxTerms)
     &      , diabioVSink_avg(NumVSinkTerms)
     &      , diabioGasExc_avg(NumGasExcTerms)
      logical wrthis(1000+NT)
     &      , wrtavg(1000+NT)
     &      , wrtdiabioFlux(NumFluxTerms+1)
     &      , wrtdiabioVSink(NumVSinkTerms+1)
     &      , wrtdiabioGasExc(NumGasExcTerms+1)
     &      , wrtdiabioFlux_avg(NumFluxTerms+1)
     &      , wrtdiabioVSink_avg(NumVSinkTerms+1)
     &      , wrtdiabioGasExc_avg(NumGasExcTerms+1)
      common/incscrum/
     &     ncidfrc, ncidbulk,ncidclm, ncidqbar, ncidbtf
     &     , ntsms, ntsrf, ntssh, ntsst
     &     , ntuclm, ntsss, ntbulk, ntqbar, ntww
     &      , xi_rho,  xi_u
     &      , eta_rho, eta_v
     &     ,  nttclm, ntstf, nttsrc, ntbtf
     &      , ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
     & ,   rstT
     &      , rstAkv,rstAkt
     &      , rstAks
     &      , rstTke,rstGls
     &      , rstBustr,rstBvstr
     &      , rstrufrc,rstrvfrc
     &      , rstru_nbq,rstrv_nbq
     &      , rstru_nbq_avg2,rstrv_nbq_avg2
     &      , rstqdmu_nbq,rstqdmv_nbq
     &      , rstdRdx,rstdRde
     &      , ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisBustr, hisBvstr
     &      , hisShflx, hisSwflx, hisShflx_rsw
     &      , hisBhflx, hisBwflx
     &      , hisUnbq, hisVnbq, hisWnbq, hisRnbq, hisCnbq
     &      , hisU,    hisV,     hisT,    hisR
     &      , hisO,    hisW,     hisVisc, hisDiff
     &      , hisAkv,  hisAkt,   hisAks
     &      , hisHbl,  hisHbbl
     &      , hisbvf
     &      , hisTke, hisGls, hisLsc
     &      , hisShflx_rlw
     &      , hisShflx_lat, hisShflx_sen
     &      , hisHel
     &      , nciddiabio, nrecdiabio, nrpfdiabio
     &      , diaTimebio, diaTime2bio, diaTstepbio, diabioFlux
     &      , diabioVSink
     &      , diabioGasExc
     &      , nciddiabio_avg, nrecdiabio_avg, nrpfdiabio_avg
     &      , diaTimebio_avg, diaTime2bio_avg, diaTstepbio_avg
     &      , diabioFlux_avg
     &      , diabioVSink_avg
     &      , diabioGasExc_avg
     &      , ncidavg,  nrecavg,  nrpfavg
     &      , avgTime, avgTime2, avgTstep, avgZ,    avgUb,  avgVb
     &      , avgBostr, avgWstr, avgUWstr, avgVWstr
     &      , avgBustr, avgBvstr
     &      , avgShflx, avgSwflx, avgShflx_rsw
     &      , avgBhflx, avgBwflx
     &      , avgU,    avgV
     &      ,     avgT
     &      ,     avgR
     &      , avgO,    avgW,     avgVisc,  avgDiff
     &      , avgAkv,  avgAkt,   avgAks
     &      , avgHbl,  avgHbbl
     &      , avgbvf
     &      , avgTke, avgGls, avgLsc
     &      , avgHel
     &      , avgShflx_rlw
     &      , avgShflx_lat, avgShflx_sen
     &      , wrthis
     &      , wrtavg
     &      , wrtdiabioFlux
     &      , wrtdiabioVSink
     &      , wrtdiabioGasExc
     &      , wrtdiabioFlux_avg
     &      , wrtdiabioVSink_avg
     &      , wrtdiabioGasExc_avg
      character*80 date_str, title
      character*80 origin_date, start_date_run, xios_origin_date
      integer*4      start_day, start_month, start_year
     &         ,   start_hour, start_minute, start_second
     &         ,   origin_day, origin_month, origin_year
     &         ,   origin_hour, origin_minute, origin_second
      REAL(kind=8) :: origin_date_in_sec, xios_origin_date_in_sec
      character*180 ininame,  grdname,  hisname
     &         ,   rstname,  frcname,  bulkname,  usrname
     &         ,   qbarname, tsrcname
     &         ,   btfname
     &                                ,  avgname
     &                                ,  dianamebio
     &                                ,  dianamebio_avg
     &                                ,   clmname
     &                                ,   bioname
      character*75  vname(20, 1000)
      common /cncscrum/   date_str,   title
     &         ,   origin_date, start_date_run
     &         ,   xios_origin_date
     &         ,   ininame,  grdname, hisname
     &         ,   rstname,  frcname, bulkname,  usrname
     &         ,   qbarname, tsrcname
     &         ,   btfname, origin_date_in_sec
     &         ,   xios_origin_date_in_sec
     &         ,   start_day, start_month, start_year
     &         ,   start_hour, start_minute, start_second
     &         ,   origin_day, origin_month, origin_year
     &         ,   origin_hour, origin_minute, origin_second
     &                                ,  avgname
     &                                ,  dianamebio
     &                                ,  dianamebio_avg
     &                                ,   clmname
     &                                ,   bioname
     &                                ,   vname
      end module ncscrum
