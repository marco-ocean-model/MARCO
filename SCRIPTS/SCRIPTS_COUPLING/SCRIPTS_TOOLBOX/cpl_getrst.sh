#-------------------------------------------------------------------------------
#                                                                      Restart
#-------------------------------------------------------------------------------

if [[ ${RESTART_FLAG} == "FALSE" ]]; then

    module load $ncomod
#
    if [ ${USE_ATM} -eq 1 ] ; then
        maxatmdom=$( echo $wrfcpldom | wc -w )
        ocemaxdom=0
        for domatm in $wrfcpldom ; do
            varlist=""
            if [[ ${onlinecplmask} == "TRUE" ]]; then
                [[ $(( $AGRIFZ +1 )) > $maxatmdom ]] && { loopoce=`seq 0 $AGRIFZ` ;} || { loopoce=`seq 0 $(( $maxatmdom - 1 ))` ;}
            else
                loopoce=`seq 0 $(( $maxatmdom - 1 ))`
            fi

            atm_rst_file_dom=$(echo "$atm_rst_file" | sed "s/_dXX_/_${domatm}_/")
            if [ ! -f $atm_rst_file_dom ] ; then
                printf "\n ERROR: $atm_rst_file_dom not found. Exit"
                exit 1
            else
                if [[ ${CPL_restart} == "TRUE" ]] ; then
                    printf "Create atmosphere restart file for oasis from preexisting file ${atm_rst_file_dom}" #|tee cpl_getrst.log   
                    printf "\n  -> See log file ${EXEDIR}/cpl_getrst.log\n"

                    for ocedom in $loopoce; do
                        domoce="d0$(( ${ocedom} + 1 ))"
                        if [ ${domatm} == "d01" ]; then
                            . ${SCRIPTDIR}/OASIS_SCRIPTS/create_oasis_restart_from_preexisting_output_files.sh "${atm_rst_file_dom}" atm.nc wrf "WRF_${domatm}_EXT_${domoce}" >> cpl_getrst.log
                        else
                            . ${SCRIPTDIR}/OASIS_SCRIPTS/create_oasis_restart_from_preexisting_output_files.sh "${atm_rst_file_dom}" atm${domatm}.nc wrf "WRF_${domatm}_EXT_${domoce}" >> cpl_getrst.log
                        fi
                    done

                else
                    printf "Create restart file for oasis from calm conditions for variables: ${varlist}"  #|tee cpl_getrst.log
                    printf "\n  -> See log file ${EXEDIR}/cpl_getrst.log\n"

                    for ocedom in $loopoce; do
                        domoce="d0$(( ${ocedom} + 1 ))"
                        varlist="${varlist}WRF_${domatm}_EXT_${domoce}_SURF_NET_SOLAR WRF_${domatm}_EXT_${domoce}_EVAP-PRECIP WRF_${domatm}_EXT_${domoce}_SURF_NET_NON-SOLAR WRF_${domatm}_EXT_${domoce}_TAUE WRF_${domatm}_EXT_${domoce}_TAUN WRF_${domatm}_EXT_${domoce}_TAUMOD WRF_${domatm}_EXT_${domoce}_PSFC WRF_${domatm}_EXT_${domoce}_WND_E_01 WRF_${domatm}_EXT_${domoce}_WND_N_01 "
                    done

                    if [ ${domatm} == "d01" ]; then
                        . ${SCRIPTDIR}/OASIS_SCRIPTS/create_oasis_restart_from_calm_conditions.sh "${atm_rst_file_dom}" atm.nc wrf "${varlist}"  >> cpl_getrst.log
                    else 
                        . ${SCRIPTDIR}/OASIS_SCRIPTS/create_oasis_restart_from_calm_conditions.sh "${atm_rst_file_dom}" atm${domatm}.nc wrf "${varlist}"  >> cpl_getrst.log
                    fi
                fi
            fi
        done
    fi
#
    if [ ${USE_OCE} -eq 1 ] ; then
        for nn in `seq 0 ${AGRIFZ}`; do
            varlist=""
            if [ ${nn} -gt 0 ]; then
                agrif_ext=".${nn}"
            else
                agrif_ext=""
            fi
	    [[ ${AGRIFZ} > 0 ]] && { mm="_${nn}" ;} || { mm="" ;}
            varlist="${varlist}CROCO_SST${mm} CROCO_SSH${mm} CROCO_NOCE${mm} CROCO_EOCE${mm} "
            if [ ! -f ${oce_rst_file}${agrif_ext} ] ; then
                printf "\n ERROR: ${oce_rst_file}${agrif_ext} not found. Exit"
                exit 1
            else
                if [[ ${CPL_restart} == "TRUE" ]] ; then
                    printf "Create ocean restart file for oasis from preexisting file ${oce_rst_file}${agrif_ext}"  #|tee cpl_getrst.log     
                    printf "\n  -> See log file ${EXEDIR}/cpl_getrst.log\n"
                    . ${SCRIPTDIR}/OASIS_SCRIPTS/create_oasis_restart_from_preexisting_output_files.sh "${oce_rst_file}${agrif_ext}" oce.nc${agrif_ext} croco ${mm}  >> cpl_getrst.log
                else
                    printf "create restart file for oasis from calm conditions for variables: ${varlist}"  #|tee cpl_getrst.log     
                    printf "\n  -> See log file ${EXEDIR}/cpl_getrst.log\n"
                    . ${SCRIPTDIR}/OASIS_SCRIPTS/create_oasis_restart_from_calm_conditions.sh ${oce_rst_file}${agrif_ext} oce.nc${agrif_ext} croco "${varlist}"  >> cpl_getrst.log
                fi 
            fi           
        done
    fi
#
    if [ ${USE_WAV} -eq 1 ] ; then
        varlist='WW3_T0M1 WW3__OHS WW3__DIR WW3_ACHA WW3_TAWX WW3_TAWY WW3_TWOX WW3_TWOY WW3_USSX WW3_USSY WW3__FOC WW3___LM WW3__BHD WW3_UBRX WW3_UBRY'
        if [ ! -f ${wav_rst_file} ] ; then
            printf "\n ERROR: ${wav_rst_file} not found. Exit"
            exit 1
        else
            if [[ ${CPL_restart} == "TRUE" ]] ; then
                printf "Create wave restart file for oasis from preexisting file ${wav_rst_file}" #|tee cpl_getrst.log
                printf "\n  -> See log file ${EXEDIR}/cpl_getrst.log\n"
                . ${SCRIPTDIR}/OASIS_SCRIPTS/create_oasis_restart_from_preexisting_output_files.sh "${wav_rst_file}" wav.nc ww3  >> cpl_getrst.log
            else
                printf "Create restart file for oasis from calm conditions for variables: ${varlist}" #|tee cpl_getrst.log
                printf "\n  -> See log file ${EXEDIR}/cpl_getrst.log\n"
                . ${SCRIPTDIR}/OASIS_SCRIPTS/create_oasis_restart_from_calm_conditions.sh ${wav_rst_file} wav.nc ww3 "${varlist}"  >> cpl_getrst.log
            fi
        fi
    fi
#
    if [ ${USE_TOY} -ge 1 ] ; then
        varlist='TOY_V_01 TOY_U_01 TOY_TAUX TOY_TAUY TOY_TAUM TOYSRFLX TOYSTFLX TOY__EMP TOY_UOCE TOY_VOCE TOY_PSFC TOY__SST TOY__SSH TOY_T0M1 TOY___HS TOY__DIR TOY_TWOX TOY_TWOY TOY_TAWX TOY_TAWY TOY__CHA'

        printf "Create restart file for oasis from calm conditions for variables: ${varlist}" #|tee cpl_getrst.log
        printf "\n  -> See log file ${EXEDIR}/cpl_getrst.log\n"

        for k in `seq 0 $(( ${nbtoy} - 1 ))`; do
            . ${SCRIPTDIR}/OASIS_SCRIPTS/create_oasis_restart_from_calm_conditions.sh ${toyfile[$k]} ${toytype[$k]}.nc ${model_to_toy[$k]} "$varlist"  >> cpl_getrst.log
        done
    fi
    module unload $ncomod

else   

    if [ ${USE_ATM} -eq 1 ]; then
        for dom in $wrfcpldom; do
            if [ $dom == "d01" ]; then
                cpfile ${RESTDIR_IN}/atm_${CEXPER}_${DATE_END_JOBm1}.nc atm.nc
            else
                cpfile ${RESTDIR_IN}/atm${dom}_${CEXPER}_${DATE_END_JOBm1}.nc atm${dom}.nc
            fi
        done
    fi
    if [ ${USE_OCE} -eq 1 ]; then
        for nn in $( seq 0 ${AGRIFZ} ); do
            if [ ${nn} -gt 0 ];    then
                agrif_ext=".${nn}"
            else
                agrif_ext=""
            fi
            cpfile ${RESTDIR_IN}/oce_${CEXPER}_${DATE_END_JOBm1}.nc${agrif_ext} oce.nc${agrif_ext}
        done
    fi

    [ ${USE_WAV} -eq 1 ] && cpfile ${RESTDIR_IN}/wav_${CEXPER}_${DATE_END_JOBm1}.nc wav.nc 

    [[ ${USE_ATM} -eq 1 && ${USE_OCE} -eq 1 ]] && cp ${RESTDIR_IN}/*atmt_to_ocn* . && cp ${RESTDIR_IN}/*ocn*_to_atmt* . 
    [[ ${USE_ATM} -eq 1 && ${USE_WAV} -eq 1 ]] && cp ${RESTDIR_IN}/*atmt_to_ww3t* . && cp ${RESTDIR_IN}/*ww3t_to_atmt* .
    [[ ${USE_OCE} -eq 1 && ${USE_WAV} -eq 1 ]] && cp ${RESTDIR_IN}/*ocn*_to_ww3t* . && cp ${RESTDIR_IN}/*ww3t_to_ocn* .  

    if [ ${USE_TOY} -ge 1 ] ; then
        for k in `seq 0 $(( ${nbtoy} - 1 ))`; do
            cpfile ${RESTDIR_IN}/${toytype[$k]}_${CEXPER}_${DATE_END_JOBm1}.nc ${toytype[$k]}.nc
        done
	[ ${USE_OCE} -eq 1 ] && cp ${RESTDIR_IN}/*toy*_to_ocn* . && cp ${RESTDIR_IN}/*ocn*_to_toy* . 
	[ ${USE_WAV} -eq 1 ] && cp ${RESTDIR_IN}/*toy*_to_ww3t* . && cp ${RESTDIR_IN}/*ww3t_to_toy* .
        [ ${USE_ATM} -eq 1 ] && cp ${RESTDIR_IN}/*toy*_to_atmt* . && cp ${RESTDIR_IN}/*atmt_to_toy* .
        if [ ${nbtoy} -gt 1 ]; then
            cp ${RESTDIR_IN}/*toy*_to_toy* .
        fi
    fi
    

    cpfile2 ${RESTDIR_IN}/grids.nc ./ ; cpfile2 ${RESTDIR_IN}/masks.nc  ./ ; cpfile2 ${RESTDIR_IN}/areas.nc ./
    

fi
