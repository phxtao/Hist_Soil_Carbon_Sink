function fun_save_proda_transiesnt_summary(save_pathway, ...
    cesm2_simu_soc_stock, cesm2_simu_tot_eco_c, cesm2_simu_tot_veg_c, cesm2_simu_tot_lit_c, ...
    cpool_12c_stock_transient, cpool_12c_resp_transient, cpool_14c_stock_transient, coop_14c_resp_transient, ...
    cpool_potential_transient, cpol_capacity_transient, ...
    cpool_residence_time_transient, eco_process_transient)

transient_summary.cesm2_simu_soc_stock = cesm2_simu_soc_stock;
transient_summary.cesm2_simu_tot_eco_c = cesm2_simu_tot_eco_c;
transient_summary.cesm2_simu_tot_veg_c = cesm2_simu_tot_veg_c;
transient_summary.cesm2_simu_tot_lit_c = cesm2_simu_tot_lit_c;

transient_summary.cpool_12c_stock_transient = cpool_12c_stock_transient;
transient_summary.cpool_12c_resp_transient = cpool_12c_resp_transient;
transient_summary.cpool_14c_stock_transient = cpool_14c_stock_transient;
transient_summary.coop_14c_resp_transient = coop_14c_resp_transient;
transient_summary.cpool_potential_transient = cpool_potential_transient;
transient_summary.cpol_capacity_transient = cpol_capacity_transient;
transient_summary.cpool_residence_time_transient = cpool_residence_time_transient;
transient_summary.eco_process_transient = eco_process_transient;


save([save_pathway, '.mat'], 'transient_summary');

end
