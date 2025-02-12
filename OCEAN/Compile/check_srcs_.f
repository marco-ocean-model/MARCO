      subroutine check_srcs
      implicit none
      integer*4 i
      integer*4 max_opt_size
      parameter (max_opt_size=5000)
      character*5000 Coptions,srcs
      common /strings/ Coptions,srcs
      do i=1,max_opt_size
        srcs(i:i)=' '
      enddo
      srcs(   1:   6)='main.F'
      srcs(   7:  13)=' step.F'
      srcs(  14:  24)=' read_inp.F'
      srcs(  25:  38)=' timers_roms.F'
      srcs(  39:  53)=' init_scalars.F'
      srcs(  54:  67)=' init_arrays.F'
      srcs(  68:  81)=' set_weights.F'
      srcs(  82:  94)=' set_scoord.F'
      srcs(  95: 105)=' ana_grid.F'
      srcs( 106: 119)=' setup_grid1.F'
      srcs( 120: 133)=' setup_grid2.F'
      srcs( 134: 147)=' set_nudgcof.F'
      srcs( 148: 161)=' ana_initial.F'
      srcs( 162: 174)=' analytical.F'
      srcs( 175: 183)=' zonavg.F'
      srcs( 184: 192)=' step2d.F'
      srcs( 193: 200)=' u2dbc.F'
      srcs( 201: 208)=' v2dbc.F'
      srcs( 209: 217)=' zetabc.F'
      srcs( 218: 231)=' obc_volcons.F'
      srcs( 232: 244)=' pre_step3d.F'
      srcs( 245: 255)=' step3d_t.F'
      srcs( 256: 268)=' step3d_uv1.F'
      srcs( 269: 281)=' step3d_uv2.F'
      srcs( 282: 290)=' prsgrd.F'
      srcs( 291: 298)=' rhs3d.F'
      srcs( 299: 310)=' set_depth.F'
      srcs( 311: 318)=' omega.F'
      srcs( 319: 328)=' uv3dmix.F'
      srcs( 329: 342)=' uv3dmix_spg.F'
      srcs( 343: 351)=' t3dmix.F'
      srcs( 352: 364)=' t3dmix_spg.F'
      srcs( 365: 376)=' hmix_coef.F'
      srcs( 377: 385)=' wetdry.F'
      srcs( 386: 395)=' wavedry.F'
      srcs( 396: 403)=' u3dbc.F'
      srcs( 404: 411)=' v3dbc.F'
      srcs( 412: 419)=' t3dbc.F'
      srcs( 420: 433)=' step3d_fast.F'
      srcs( 434: 444)=' step3d_w.F'
      srcs( 445: 457)=' rhs3d_w_nh.F'
      srcs( 458: 471)=' initial_nbq.F'
      srcs( 472: 482)=' grid_nbq.F'
      srcs( 483: 492)=' unbq_bc.F'
      srcs( 493: 502)=' vnbq_bc.F'
      srcs( 503: 512)=' wnbq_bc.F'
      srcs( 513: 522)=' rnbq_bc.F'
      srcs( 523: 530)=' w3dbc.F'
      srcs( 531: 546)=' nbq_bry_store.F'
      srcs( 547: 556)=' rho_eos.F'
      srcs( 557: 567)=' ab_ratio.F'
      srcs( 568: 578)=' alfabeta.F'
      srcs( 579: 591)=' alfabeta_k.F'
      srcs( 592: 602)=' ana_vmix.F'
      srcs( 603: 612)=' bvf_mix.F'
      srcs( 613: 623)=' lmd_vmix.F'
      srcs( 624: 636)=' gls_mixing.F'
      srcs( 637: 647)=' lmd_skpp.F'
      srcs( 648: 658)=' lmd_bkpp.F'
      srcs( 659: 671)=' lmd_swfrac.F'
      srcs( 672: 684)=' lmd_wscale.F'
      srcs( 685: 691)=' diag.F'
      srcs( 692: 700)=' wvlcty.F'
      srcs( 701: 712)=' checkdims.F'
      srcs( 713: 729)=' grid_stiffness.F'
      srcs( 730: 740)=' bio_diag.F'
      srcs( 741: 753)=' setup_kwds.F'
      srcs( 754: 766)=' check_kwds.F'
      srcs( 767: 779)=' check_srcs.F'
      srcs( 780: 797)=' check_switches1.F'
      srcs( 798: 815)=' check_switches2.F'
      srcs( 816: 823)=' debug.F'
      srcs( 824: 831)=' param.F'
      srcs( 832: 841)=' ncscrum.F'
      srcs( 842: 851)=' scalars.F'
      srcs( 852: 860)=' output.F'
      srcs( 861: 878)=' put_global_atts.F'
      srcs( 879: 889)=' nf_fread.F'
      srcs( 890: 902)=' nf_fread_x.F'
      srcs( 903: 915)=' nf_fread_y.F'
      srcs( 916: 929)=' nf_read_bry.F'
      srcs( 930: 940)=' get_date.F'
      srcs( 941: 949)=' lenstr.F'
      srcs( 950: 960)=' closecdf.F'
      srcs( 961: 974)=' insert_node.F'
      srcs( 975: 986)=' fillvalue.F'
      srcs( 987:1005)=' nf_add_attribute.F'
      srcs(1006:1017)=' set_cycle.F'
      srcs(1018:1031)=' def_grid_2d.F'
      srcs(1032:1045)=' def_grid_3d.F'
      srcs(1046:1055)=' def_his.F'
      srcs(1056:1065)=' def_rst.F'
      srcs(1066:1077)=' def_diags.F'
      srcs(1078:1090)=' def_diagsM.F'
      srcs(1091:1106)=' def_bio_diags.F'
      srcs(1107:1117)=' wrt_grid.F'
      srcs(1118:1127)=' wrt_his.F'
      srcs(1128:1137)=' wrt_avg.F'
      srcs(1138:1147)=' wrt_rst.F'
      srcs(1148:1159)=' wrt_diags.F'
      srcs(1160:1175)=' wrt_diags_avg.F'
      srcs(1176:1188)=' wrt_diagsM.F'
      srcs(1189:1205)=' wrt_diagsM_avg.F'
      srcs(1206:1221)=' wrt_bio_diags.F'
      srcs(1222:1241)=' wrt_bio_diags_avg.F'
      srcs(1242:1251)=' set_avg.F'
      srcs(1252:1267)=' set_diags_avg.F'
      srcs(1268:1284)=' set_diagsM_avg.F'
      srcs(1285:1304)=' set_bio_diags_avg.F'
      srcs(1305:1320)=' def_diags_vrt.F'
      srcs(1321:1336)=' wrt_diags_vrt.F'
      srcs(1337:1352)=' set_diags_vrt.F'
      srcs(1353:1372)=' set_diags_vrt_avg.F'
      srcs(1373:1392)=' wrt_diags_vrt_avg.F'
      srcs(1393:1407)=' def_diags_ek.F'
      srcs(1408:1422)=' wrt_diags_ek.F'
      srcs(1423:1437)=' set_diags_ek.F'
      srcs(1438:1456)=' set_diags_ek_avg.F'
      srcs(1457:1475)=' wrt_diags_ek_avg.F'
      srcs(1476:1490)=' def_diags_pv.F'
      srcs(1491:1505)=' wrt_diags_pv.F'
      srcs(1506:1520)=' set_diags_pv.F'
      srcs(1521:1539)=' set_diags_pv_avg.F'
      srcs(1540:1558)=' wrt_diags_pv_avg.F'
      srcs(1559:1575)=' def_diags_eddy.F'
      srcs(1576:1596)=' set_diags_eddy_avg.F'
      srcs(1597:1617)=' wrt_diags_eddy_avg.F'
      srcs(1618:1628)=' def_surf.F'
      srcs(1629:1639)=' wrt_surf.F'
      srcs(1640:1654)=' set_surf_avg.F'
      srcs(1655:1669)=' wrt_surf_avg.F'
      srcs(1670:1680)=' get_grid.F'
      srcs(1681:1694)=' get_initial.F'
      srcs(1695:1704)=' get_vbc.F'
      srcs(1705:1716)=' get_wwave.F'
      srcs(1717:1729)=' get_tclima.F'
      srcs(1730:1742)=' get_uclima.F'
      srcs(1743:1752)=' get_ssh.F'
      srcs(1753:1762)=' get_sss.F'
      srcs(1763:1775)=' get_smflux.F'
      srcs(1776:1788)=' get_stflux.F'
      srcs(1789:1801)=' get_srflux.F'
      srcs(1802:1811)=' get_sst.F'
      srcs(1812:1827)=' mod_tides_mas.F'
      srcs(1828:1838)=' tidedata.F'
      srcs(1839:1844)=' mas.F'
      srcs(1845:1856)=' get_tides.F'
      srcs(1857:1868)=' clm_tides.F'
      srcs(1869:1879)=' get_bulk.F'
      srcs(1880:1891)=' bulk_flux.F'
      srcs(1892:1901)=' get_bry.F'
      srcs(1902:1915)=' get_bry_bio.F'
      srcs(1916:1925)=' sstskin.F'
      srcs(1926:1939)=' get_psource.F'
      srcs(1940:1956)=' get_psource_ts.F'
      srcs(1957:1969)=' get_btflux.F'
      srcs(1970:1979)=' mrl_wci.F'
      srcs(1980:1991)=' wkb_wwave.F'
      srcs(1992:1999)=' wkbbc.F'
      srcs(2000:2013)=' get_bry_wkb.F'
      srcs(2014:2031)=' online_bulk_var.F'
      srcs(2032:2049)=' online_get_bulk.F'
      srcs(2050:2065)=' online_interp.F'
      srcs(2066:2091)=' online_interpolate_bulk.F'
      srcs(2092:2109)=' online_set_bulk.F'
      srcs(2110:2138)=' online_interpolate_3d_bulk.F'
      srcs(2139:2152)=' init_floats.F'
      srcs(2153:2165)=' wrt_floats.F'
      srcs(2166:2179)=' step_floats.F'
      srcs(2180:2192)=' rhs_floats.F'
      srcs(2193:2205)=' interp_rho.F'
      srcs(2206:2218)=' def_floats.F'
      srcs(2219:2239)=' init_arrays_floats.F'
      srcs(2240:2253)=' random_walk.F'
      srcs(2254:2274)=' get_initial_floats.F'
      srcs(2275:2285)=' init_sta.F'
      srcs(2286:2295)=' wrt_sta.F'
      srcs(2296:2306)=' step_sta.F'
      srcs(2307:2319)=' interp_sta.F'
      srcs(2320:2329)=' def_sta.F'
      srcs(2330:2347)=' init_arrays_sta.F'
      srcs(2348:2357)=' biology.F'
      srcs(2358:2366)=' o2sato.F'
      srcs(2367:2377)=' sediment.F'
      srcs(2378:2383)=' bbl.F'
      srcs(2384:2395)=' MPI_Setup.F'
      srcs(2396:2408)=' MessPass2D.F'
      srcs(2409:2421)=' MessPass3D.F'
      srcs(2422:2432)=' exchange.F'
      srcs(2433:2445)=' autotiling.F'
      srcs(2446:2452)=' zoom.F'
      srcs(2453:2463)=' update2D.F'
      srcs(2464:2482)=' set_nudgcof_fine.F'
      srcs(2483:2494)=' zoombc_2D.F'
      srcs(2495:2506)=' zoombc_3D.F'
      srcs(2507:2519)=' uv3dpremix.F'
      srcs(2520:2531)=' t3dpremix.F'
      srcs(2532:2542)=' update3D.F'
      srcs(2543:2558)=' zoombc_3Dfast.F'
      srcs(2559:2572)=' Agrif2Model.F'
      srcs(2573:2590)=' send_xios_diags.F'
      srcs(2591:2609)=' cpl_prism_define.F'
      srcs(2610:2625)=' cpl_prism_put.F'
      srcs(2626:2642)=' cpl_prism_init.F'
      srcs(2643:2658)=' cpl_prism_get.F'
      srcs(2659:2677)=' cpl_prism_getvar.F'
      srcs(2678:2694)=' cpl_prism_grid.F'
      srcs(2695:2713)=' module_substance.F'
      srcs(2714:2730)=' module_MUSTANG.F'
      srcs(2731:2752)=' module_OBSTRUCTIONS.F'
      srcs(2753:2763)=' abl_step.F'
      srcs(2764:2773)=' abl_ini.F'
      srcs(2774:2783)=' abl_tke.F'
      srcs(2784:2797)=' wrt_abl_his.F'
      srcs(2798:2811)=' def_abl_his.F'
      srcs(2812:2825)=' wrt_abl_avg.F'
      srcs(2826:2839)=' set_abl_avg.F'
      srcs(2840:2863)=' online_spectral_diags.F'
      srcs(2864:2878)=' var3d_oa_out.F'
      srcs(2879:2893)=' var2d_oa_out.F'
      srcs(2894:2913)=' scal0d_oa_out_loc.F'
      srcs(2914:2934)=' scal0d_oa_out_full.F'
      srcs(2935:2946)=' output_oa.F'
      srcs(2947:2964)=' copy_to_devices.F'
      srcs(2965:2987)=' exchange_device_host.F'
      srcs(2988:2998)=' init_acc.F'
      srcs(2999:3013)=' par_pisces.F90'
      srcs(3014:3025)=' oce_trc.F90'
      srcs(3026:3033)=' trc.F90'
      srcs(3034:3048)=' sms_pisces.F90'
      srcs(3049:3067)=' in_out_manager.F90'
      srcs(3068:3075)=' iom.F90'
      srcs(3076:3087)=' lib_mpp.F90'
      srcs(3088:3098)=' prtctl.F90'
      srcs(3099:3109)=' p4zche.F90'
      srcs(3110:3120)=' p4zint.F90'
      srcs(3121:3131)=' p4zlys.F90'
      srcs(3132:3142)=' p4zflx.F90'
      srcs(3143:3153)=' p4zlim.F90'
      srcs(3154:3165)=' p4zsink.F90'
      srcs(3166:3178)=' p4zmicro.F90'
      srcs(3179:3190)=' p4zmeso.F90'
      srcs(3191:3202)=' p4zmort.F90'
      srcs(3203:3213)=' p4zopt.F90'
      srcs(3214:3225)=' p4zprod.F90'
      srcs(3226:3236)=' p4zrem.F90'
      srcs(3237:3246)=' p4zbc.F90'
      srcs(3247:3257)=' p4zsed.F90'
      srcs(3258:3269)=' p4zdiaz.F90'
      srcs(3270:3280)=' p4zagg.F90'
      srcs(3281:3294)=' p4zfechem.F90'
      srcs(3295:3308)=' p4zligand.F90'
      srcs(3309:3319)=' p4zpoc.F90'
      srcs(3320:3330)=' p2zlim.F90'
      srcs(3331:3343)=' p2zmicro.F90'
      srcs(3344:3355)=' p2zmort.F90'
      srcs(3356:3367)=' p2zprod.F90'
      srcs(3368:3380)=' p5zmicro.F90'
      srcs(3381:3392)=' p5zmort.F90'
      srcs(3393:3404)=' p5zprod.F90'
      srcs(3405:3415)=' p5zlim.F90'
      srcs(3416:3427)=' p5zmeso.F90'
      srcs(3428:3438)=' p4zbio.F90'
      srcs(3439:3449)=' p4zsms.F90'
      srcs(3450:3461)=' trcsink.F90'
      srcs(3462:3479)=' trcwri_pisces.F90'
      srcs(3480:3497)=' trcnam_pisces.F90'
      srcs(3498:3515)=' trcsms_pisces.F90'
      srcs(3516:3533)=' trcini_pisces.F90'
      srcs(3534:3548)=' pisces_ini.F90'
      srcs(3549:3560)=' oce_sed.F90'
      srcs(3561:3572)=' par_sed.F90'
      srcs(3573:3583)=' sedadv.F90'
      srcs(3584:3594)=' sedbtb.F90'
      srcs(3595:3606)=' sedchem.F90'
      srcs(3607:3617)=' sedco3.F90'
      srcs(3618:3628)=' seddsr.F90'
      srcs(3629:3639)=' seddta.F90'
      srcs(3640:3647)=' sed.F90'
      srcs(3648:3661)=' seddsrjac.F90'
      srcs(3662:3673)=' sedfunc.F90'
      srcs(3674:3684)=' sedjac.F90'
      srcs(3685:3695)=' sedsol.F90'
      srcs(3696:3706)=' sedini.F90'
      srcs(3707:3720)=' sedinitrc.F90'
      srcs(3721:3733)=' sedinorg.F90'
      srcs(3734:3744)=' sedmat.F90'
      srcs(3745:3757)=' sedmodel.F90'
      srcs(3758:3769)=' sed_oce.F90'
      srcs(3770:3780)=' sedorg.F90'
      srcs(3781:3791)=' sedrst.F90'
      srcs(3792:3802)=' sedsfc.F90'
      srcs(3803:3813)=' sedstp.F90'
      srcs(3814:3824)=' sedwri.F90'
      srcs(3825:3834)=' trosk.F90'
      srcs(3835:3849)=' setavg_sed.F90'
      srcs(3850:3873)=' module_parameter_oa.F90'
      srcs(3874:3891)=' module_grd_oa.F90'
      srcs(3892:3910)=' module_oa_time.F90'
      srcs(3911:3930)=' module_oa_space.F90'
      srcs(3931:3952)=' module_oa_periode.F90'
      srcs(3953:3976)=' module_oa_variables.F90'
      srcs(3977:3995)=' module_oa_type.F90'
      srcs(3996:4014)=' module_tile_oa.F90'
      srcs(4015:4034)=' module_oa_level.F90'
      srcs(4035:4058)=' module_oa_interface.F90'
      srcs(4059:4069)=' var_oa.F90'
      srcs(4070:4085)=' tooldatosec.F90'
      srcs(4086:4102)=' toolsectodat.F90'
      srcs(4103:4120)=' tooldecompdat.F90'
      srcs(4121:4138)=' tooldatetosec.F90'
      srcs(4139:4157)=' toolorigindate.F90'
      srcs(4158:4169)=' flocmod.F90'
      srcs(4170:4186)=' comsubstance.F90'
      srcs(4187:4205)=' submassbalance.F90'
      srcs(4206:4219)=' substance.F90'
      srcs(4220:4239)=' OBSTRUCTIONS1DV.F90'
      srcs(4240:4260)=' com_OBSTRUCTIONS.F90'
      srcs(4261:4282)=' init_OBSTRUCTIONS.F90'
      srcs(4283:4299)=' OBSTRUCTIONS.F90'
      srcs(4300:4321)=' plug_OBSTRUCTIONS.F90'
      srcs(4322:4336)=' comMUSTANG.F90'
      srcs(4337:4356)=' coupler_MUSTANG.F90'
      srcs(4357:4378)=' sed_MUSTANG_CROCO.F90'
      srcs(4379:4394)=' sed_MUSTANG.F90'
      srcs(4395:4410)=' initMUSTANG.F90'
      srcs(4411:4433)=' plug_MUSTANG_CROCO.F90'
      return
      end
