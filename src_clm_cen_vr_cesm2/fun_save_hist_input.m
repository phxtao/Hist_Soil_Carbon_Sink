function fun_save_hist_input(save_pathway, ...
    cesm2_simu_nbedrock, cesm2_simu_altmax_origin, cesm2_simu_altmax_last_year_origin, cesm2_simu_cellsand, ...
    cesm2_simu_npp, cesm2_simu_soil_water_potnetial_origin, cesm2_simu_soil_temperature_origin, ...
    cesm2_simu_w_scalar_origin, cesm2_simu_t_scalar_origin, cesm2_simu_o_scalar_origin, cesm2_simu_n_scalar_origin, ...
    cesm2_simu_input_vector_litter1_origin, cesm2_simu_input_vector_litter2_origin, cesm2_simu_input_vector_litter3_origin, cesm2_simu_input_vector_cwd_origin, ...
    cesm2_simu_soc_stock, cesm2_simu_tot_eco_c, cesm2_simu_tot_veg_c, cesm2_simu_tot_lit_c)

grid_hist_input.cesm2_simu_nbedrock = cesm2_simu_nbedrock;
grid_hist_input.cesm2_simu_altmax_origin = cesm2_simu_altmax_origin;
grid_hist_input.cesm2_simu_altmax_last_year_origin = cesm2_simu_altmax_last_year_origin;
grid_hist_input.cesm2_simu_cellsand = cesm2_simu_cellsand;

grid_hist_input.cesm2_simu_npp = cesm2_simu_npp;
grid_hist_input.cesm2_simu_soil_water_potnetial_origin = cesm2_simu_soil_water_potnetial_origin;
grid_hist_input.cesm2_simu_soil_temperature_origin = cesm2_simu_soil_temperature_origin;

grid_hist_input.cesm2_simu_w_scalar_origin = cesm2_simu_w_scalar_origin;
grid_hist_input.cesm2_simu_t_scalar_origin = cesm2_simu_t_scalar_origin;
grid_hist_input.cesm2_simu_o_scalar_origin = cesm2_simu_o_scalar_origin;
grid_hist_input.cesm2_simu_n_scalar_origin = cesm2_simu_n_scalar_origin;

grid_hist_input.cesm2_simu_input_vector_litter1_origin = cesm2_simu_input_vector_litter1_origin;
grid_hist_input.cesm2_simu_input_vector_litter2_origin = cesm2_simu_input_vector_litter2_origin;
grid_hist_input.cesm2_simu_input_vector_litter3_origin = cesm2_simu_input_vector_litter3_origin;
grid_hist_input.cesm2_simu_input_vector_cwd_origin = cesm2_simu_input_vector_cwd_origin;

grid_hist_input.cesm2_simu_soc_stock = cesm2_simu_soc_stock;
grid_hist_input.cesm2_simu_tot_eco_c = cesm2_simu_tot_eco_c;
grid_hist_input.cesm2_simu_tot_veg_c = cesm2_simu_tot_veg_c;
grid_hist_input.cesm2_simu_tot_lit_c = cesm2_simu_tot_lit_c;

save([save_pathway, '.mat'], 'grid_hist_input');

end
