% clear;
% clc;
% %
% cross_valid_id = 7;
% worker_num = 100;
% worker_start = 1;
% worker_end = 100;
% 
% is_server = 0;
% 
% parallel setting
% delete(gcp('nocreate'));
% parpool(10);
% random_num = 1;

%% cesm2 Settings
cesm2_case_name_hat = 'hist_f05_g16_checked'; % 'sasu_f05_g16_checked_step4';
start_year = 1900;
end_year = 2010;
% start_year = 1901;
% end_year = 1905;
days_per_month = [31, 30, 31, 28, 31, 30, 31, 31, 30, 31, 30, 31];


year_num = length(start_year:end_year);

model_name = 'cesm2_clm5_cen_vr_v2';

nn_exp_name = 'exp_pc_cesm2_23';

time_domain = 'whole_time'; % 'whole_time', 'before_1985', 'after_1985'

is_direct_impact = 1;

for icross_valid = cross_valid_id
    % paths
    
    if is_server == 0
        % mac
        data_path = '/Users/ft254/DATAHUB/ENSEMBLE/';
        cd('/Users/ft254/Github/ENSEMBLE/SRC_TRANSIENT_SIMU/');
    else
        % server
        data_path = '/GFPS8p/cess11/taof/datahub/ensemble/';
        cd('/GFPS8p/cess11/taof/ensemble/src_transient_simu/');
    end
    
    nn_exp_name = ['exp_pc_cesm2_23_cross_valid_0_', num2str(icross_valid)];
    
    warning('off');
    format long e;
    
    
    %% set vertical soil pools
    month_num = 12;
    soil_cpool_num = 7;
    soil_decom_num = 20;
    %% load wosis data
    env_info = load([data_path, 'input_data/wosis_2019_snap_shot/wosis_2019_snapshot_hugelius_mishra_env_info.mat'], 'EnvInfo'); % NPP in this file will be used
    env_info = env_info.EnvInfo;
    % layer_info: "profile_id, date, upper_depth, lower_depth, node_depth, soc_layer_weight, soc_stock, bulk_denstiy, is_pedo"
    wosis_profile_info = ncread([data_path, 'input_data/wosis_2019_snap_shot/soc_profile_wosis_2019_snapshot_hugelius_mishra.nc'], 'soc_profile_info'); % wosis profile info
    wosis_soc_info = ncread([data_path, 'input_data/wosis_2019_snap_shot/soc_data_integrate_wosis_2019_snapshot_hugelius_mishra.nc'], 'data_soc_integrate'); % wosis SOC info
    
    %% 14c data
    % atmos 14c from Levin et al 2010
    atmo_14c = load([data_path, '../radiocarbon/input_data/atmo_14c_levin_2010.mat'], 'atmo_14c_levin_2010');
    atmo_14c = atmo_14c.atmo_14c_levin_2010; % 1900  - 2010
    
    %% soil depth information
    % width between two interfaces
    dz = [2.000000000000000E-002, 4.000000000000000E-002, 6.000000000000000E-002, ...
        8.000000000000000E-002, 0.120000000000000, 0.160000000000000, ...
        0.200000000000000, 0.240000000000000, 0.280000000000000, ...
        0.320000000000000, 0.360000000000000, 0.400000000000000, ...
        0.440000000000000, 0.540000000000000, 0.640000000000000, ...
        0.740000000000000, 0.840000000000000, 0.940000000000000, ...
        1.04000000000000, 1.14000000000000, 2.39000000000000, ...
        4.67553390593274, 7.63519052838329, 11.1400000000000, ...
        15.1154248593737]';
    
    % depth of the interface
    zisoi = [2.000000000000000E-002, 6.000000000000000E-002, ...
        0.120000000000000, 0.200000000000000, 0.320000000000000, ...
        0.480000000000000, 0.680000000000000, 0.920000000000000, ...
        1.20000000000000, 1.52000000000000, 1.88000000000000, ...
        2.28000000000000, 2.72000000000000, 3.26000000000000, ...
        3.90000000000000, 4.64000000000000, 5.48000000000000, ...
        6.42000000000000, 7.46000000000000, 8.60000000000000, ...
        10.9900000000000, 15.6655339059327, 23.3007244343160, ...
        34.4407244343160, 49.5561492936897]';
    
    % depth of the node
    zsoi = [1.000000000000000E-002, 4.000000000000000E-002, 9.000000000000000E-002, ...
        0.160000000000000, 0.260000000000000, 0.400000000000000, ...
        0.580000000000000, 0.800000000000000, 1.06000000000000, ...
        1.36000000000000, 1.70000000000000, 2.08000000000000, ...
        2.50000000000000, 2.99000000000000, 3.58000000000000, ...
        4.27000000000000, 5.06000000000000, 5.95000000000000, ...
        6.94000000000000, 8.03000000000000, 9.79500000000000, ...
        13.3277669529664, 19.4831291701244, 28.8707244343160, ...
        41.9984368640029]';
    
    % depth between two node
    dz_node = zsoi - [0; zsoi(1:end-1)];
    
    %% input from cesm2
    % cesm2 resolution
    cesm2_resolution_lat = 180/384;
    cesm2_resolution_lon = 360/576;
    
    cesm2_lon_grid = [(0 + cesm2_resolution_lon/2 : cesm2_resolution_lon : 180 - cesm2_resolution_lon/2), ...
        (-180 + cesm2_resolution_lon/2 : cesm2_resolution_lon : 0 - cesm2_resolution_lon/2)]';
    cesm2_lat_grid = (-90 + cesm2_resolution_lat/2 : cesm2_resolution_lat : 90 - cesm2_resolution_lat/2)';
    
    % lon and lat map
    cesm_lon_map = nan(length(cesm2_lon_grid), length(cesm2_lat_grid));
    cesm_lat_map = nan(length(cesm2_lon_grid), length(cesm2_lat_grid));
    
    for ilon = 1:length(cesm2_lon_grid)
        cesm_lon_map(ilon, :) = cesm2_lon_grid(ilon);
    end
    
    for ilat = 1:length(cesm2_lat_grid)
        cesm_lat_map(:, ilat) = cesm2_lat_grid(ilat);
    end
    
    % load cesm2 input
    var_name_list = {'nbedrock', 'ALTMAX', 'ALTMAX_LASTYEAR', 'CELLSAND', 'NPP',...
        'SOILPSI', 'TSOI', ...
        'W_SCALAR', 'T_SCALAR', 'O_SCALAR', 'FPI_vr', ...
        'LITR1_INPUT_ACC_VECTOR', 'LITR2_INPUT_ACC_VECTOR', 'LITR3_INPUT_ACC_VECTOR', 'CWD_INPUT_ACC_VECTOR', ...
        'TOTSOMC'};
    
    var_name_list_rename =  {'cesm2_simu_nbedrock', 'cesm2_simu_altmax', 'cesm2_simu_altmax_last_year', 'cesm2_simu_cellsand', 'cesm2_simu_npp',...
        'cesm2_simu_soil_water_potnetial', 'cesm2_simu_soil_temperature', ...
        'cesm2_simu_w_scalar', 'cesm2_simu_t_scalar', 'cesm2_simu_o_scalar', 'cesm2_simu_n_scalar', ...
        'cesm2_simu_input_vector_litter1', 'cesm2_simu_input_vector_litter2', 'cesm2_simu_input_vector_litter3', 'cesm2_simu_input_vector_cwd', ...
        'cesm2_simu_soc_stock'};
    
    
    %% grid information
    %--------------------------------- 0.5 deg grid info
    GlobalGrid = load([data_path, 'input_data/data4nn/world_grid_envinfo_present.mat']);
    GlobalGrid = GlobalGrid.EnvInfo;

    grid_para_result = csvread([data_path, 'output_data/neural_networking/grid_para_result_', model_name, '_', time_domain, '_', nn_exp_name, '.csv']);
    valid_grid_loc = csvread([data_path, 'output_data/neural_networking/valid_grid_loc_', model_name, '_', time_domain, '_', nn_exp_name, '.csv']);

    GlobalGrid = GlobalGrid(valid_grid_loc, :);
    GridInfo = GlobalGrid(:, [1:2, 12]);
    
    grid_interval = ceil(length(valid_grid_loc)/worker_num);
    grid_start = (worker_start-1)*grid_interval+1;
    grid_end = (worker_end)*grid_interval;
    
    %------------------------------- grid area map
    radius = 6371008.8; % unit m
    resolution_lon = 0.5;
    resolution_lat = 0.5;
    
    lon_grid_mask = (-180 + resolution_lon/2 : resolution_lon : 180 - resolution_lon/2)';
    lat_grid_mask = (90 - resolution_lat/2: -resolution_lat : -90 + resolution_lat/2)';
    
    length_top = (2*pi*radius*cos(abs(lat_grid_mask+resolution_lat/2)/180*pi)/360)*resolution_lon;
    length_down = (2*pi*radius*cos(abs(lat_grid_mask-resolution_lat/2)/180*pi)/360)*resolution_lon;
    height = (pi*radius/180)*resolution_lat;
    area_grid = (length_top + length_down)*height/2;
    
    area_map = nan(length(lat_grid_mask), length(lon_grid_mask));
    
    for ilat = 1: length(lat_grid_mask)
        area_map(ilat, :) = area_grid(ilat);
    end
    
    %----------------------------- relocation of forcing
    soc_land_mask = load([data_path, 'input_data/cesm2_simu/spinup_ss/hist_f05_g16_checked_1851_1900_cesm2_spin_up_ss_1851_1870_TOTSOMC.mat']);
    soc_land_mask = soc_land_mask.var_record;
    soc_land_mask = mean(soc_land_mask, 4, 'omitnan');
    soc_land_mask = mean(soc_land_mask, 3, 'omitnan');
    
    soc_land_mask(soc_land_mask == soc_land_mask(end, end) | soc_land_mask == 0) = nan;
    soc_land_mask_index = find(isnan(soc_land_mask) == 0);
    
    cesm_land_mask_lon_lat = [cesm_lon_map(soc_land_mask_index), cesm_lat_map(soc_land_mask_index)];
    
    grid_relocate = nan(length(valid_grid_loc), 1);
    
    for igrid  = 1:length(valid_grid_loc)
        grid_distance = (GridInfo(igrid, 1) - cesm_land_mask_lon_lat(:, 1)).^2 + (GridInfo(igrid, 2) - cesm_land_mask_lon_lat(:, 2)).^2;
        min_distance_loc = find(grid_distance == min(grid_distance, [], 'omitnan') & grid_distance <= 2);

        if isempty(min_distance_loc) == 0
            min_distance_loc = min_distance_loc(1);
            grid_relocate(igrid) = min_distance_loc;
        end
    end

    %-------------------------- carbon input allocation product
    npp_allocation_product = load([data_path, 'input_data/trendy_hist/input_allocation_data_product.mat']);
    npp_allocation_product = npp_allocation_product.summary_npp_dist;

    %% hist simulation
    if is_server == 0
        cd('/Users/ft254/Github/ENSEMBLE/SRC_POST_DA/src_clm_cen_vr_cesm2/');
    else
        cd('/GFPS8p/cess11/taof/ensemble/src_post_da/src_clm_cen_vr_cesm2/');
    end
    
    % component list
    component_list = {'input_allocation', 'vertical_transport'};
    
    iscenario = 3;
        
    for icomponent = 1:length(component_list)                
        disp([datestr(now,'HH:MM:SS'), ' processing grids from ',  num2str(grid_start), ' to ', num2str(grid_end), ' started']);
        
        soc_stock_hist_summary = nan(length(valid_grid_loc), length(1:10:111), soil_cpool_num*4); % dim: grid, time, pool and depth
        eco_process_hist_summary = nan(length(valid_grid_loc), length(1:10:111), 6); % dim: grid, time, pool and depth
        
        parfor igrid = min(grid_start, length(valid_grid_loc)) : min(grid_end, length(valid_grid_loc))
        % for igrid = min(grid_start, length(valid_grid_loc)) : min(grid_end, length(valid_grid_loc))
            disp([datestr(now,'HH:MM:SS'), ' component ', num2str(icomponent), ' processing grid ', num2str(igrid)])
            warning off
            
            grid_para = reshape(grid_para_result(igrid, :), [21, 1]);
            
            if isnan(grid_para(1)) == 1
                continue
            end
            
            grid_lon = GridInfo(igrid, 1);
            grid_lat = GridInfo(igrid, 2);

            grid_lon_coord = length((-180 + resolution_lon/2 : resolution_lon : grid_lon));
            grid_lat_coord = length((90 - resolution_lat/2: -resolution_lat : grid_lat));

            % carbon input allocation from a data product and changes of vertical transport
            if icomponent == 1
                soil_allocation_depth = [0, 10, 30, 50, 70, 90, 125, 175]';
                grid_npp_dist = reshape(npp_allocation_product(grid_lat_coord, grid_lon_coord, :), [size(npp_allocation_product, 3), 1]);
                grid_npp_allocation = [(1-grid_npp_dist(1)/100); grid_npp_dist(2:end)/sum(grid_npp_dist(2:end))*grid_npp_dist(1)/100];

                interp_allo = interp1(soil_allocation_depth/100, grid_npp_allocation, zsoi(1:soil_decom_num, 1), 'linear');

                interp_allo(interp_allo < 0) = 0;
                interp_allo(isnan(interp_allo)) = 0;
                vertical_input_product = interp_allo/sum(interp_allo);
                vertical_transport_scalar = 1;
            elseif icomponent == 2
                vertical_input_product = nan(soil_decom_num, 1);
                vertical_transport_scalar = 2;
            end

            grid_para(1:2) = grid_para(1:2)*vertical_transport_scalar;

            % scatter(soil_allocation_depth, grid_npp_allocation)
            % hold on
            % scatter(zsoi(1:soil_decom_num, 1)*100, interp_allo)

            %----------------------------------------------------
            % prepare transient forcing for each grid
            %----------------------------------------------------
            
            grid_hist_forcing = load([data_path, 'input_data/cesm2_simu/hist_simu/cesm2_hist_input_hist_f05_g16_checked_', ...
                num2str(start_year), '_', num2str(end_year), '_grid_', num2str(grid_relocate(igrid)), '.mat']);
            grid_hist_forcing = grid_hist_forcing.grid_hist_input;
            
            % load forcings from CESM2
            cesm2_simu_input_vector_cwd_origin = grid_hist_forcing.cesm2_simu_input_vector_cwd_origin;
            cesm2_simu_input_vector_litter1_origin = grid_hist_forcing.cesm2_simu_input_vector_litter1_origin;
            cesm2_simu_input_vector_litter2_origin = grid_hist_forcing.cesm2_simu_input_vector_litter2_origin;
            cesm2_simu_input_vector_litter3_origin = grid_hist_forcing.cesm2_simu_input_vector_litter3_origin;
            
            cesm2_simu_npp = grid_hist_forcing.cesm2_simu_npp;
            cesm2_simu_altmax_origin = grid_hist_forcing.cesm2_simu_altmax_origin;
            cesm2_simu_altmax_last_year_origin = grid_hist_forcing.cesm2_simu_altmax_last_year_origin;
            cesm2_simu_cellsand = grid_hist_forcing.cesm2_simu_cellsand;
            cesm2_simu_soil_water_potnetial_origin = grid_hist_forcing.cesm2_simu_soil_water_potnetial_origin;
            cesm2_simu_soil_temperature_origin = grid_hist_forcing.cesm2_simu_soil_temperature_origin;
            
            cesm2_simu_o_scalar_origin = grid_hist_forcing.cesm2_simu_o_scalar_origin;
            cesm2_simu_n_scalar_origin = grid_hist_forcing.cesm2_simu_n_scalar_origin;
            cesm2_simu_t_scalar_origin = grid_hist_forcing.cesm2_simu_t_scalar_origin;
            cesm2_simu_w_scalar_origin = grid_hist_forcing.cesm2_simu_w_scalar_origin;
            
            cesm2_simu_nbedrock = grid_hist_forcing.cesm2_simu_nbedrock;
            
            cesm2_simu_soc_stock = grid_hist_forcing.cesm2_simu_soc_stock;
            cesm2_simu_tot_eco_c = grid_hist_forcing.cesm2_simu_tot_eco_c;
            cesm2_simu_tot_veg_c = grid_hist_forcing.cesm2_simu_tot_veg_c;
            cesm2_simu_tot_lit_c = grid_hist_forcing.cesm2_simu_tot_lit_c;            
            
            ss_forcing_year_num = 20;
            
            % all changing
            if iscenario == 3
                % altmax current
                cesm2_simu_altmax = cesm2_simu_altmax_origin;
                % altmax last year
                cesm2_simu_altmax_last_year = cesm2_simu_altmax_last_year_origin;
                % temp
                cesm2_simu_soil_temperature = cesm2_simu_soil_temperature_origin;
                % water
                cesm2_simu_w_scalar = cesm2_simu_w_scalar_origin;
                % nitrogen
                cesm2_simu_n_scalar = cesm2_simu_n_scalar_origin;
                % oxygen
                cesm2_simu_o_scalar = cesm2_simu_o_scalar_origin;
                % input cwd
                cesm2_simu_input_vector_cwd = cesm2_simu_input_vector_cwd_origin;
                % input l1
                cesm2_simu_input_vector_litter1 = cesm2_simu_input_vector_litter1_origin;
                % input l2
                cesm2_simu_input_vector_litter2 = cesm2_simu_input_vector_litter2_origin;
                % input l3
                cesm2_simu_input_vector_litter3 = cesm2_simu_input_vector_litter3_origin;
            end
            
            % fully varying variable
            [cpool_12c_stock_transient, cpool_12c_resp_transient, ...
                cpool_14c_stock_transient, coop_14c_resp_transient, ...
                cpool_potential_transient, cpol_capacity_transient, ...
                cpool_residence_time_transient, eco_process_transient] = ...
                fun_hist_simu_seasonal(grid_para, atmo_14c, start_year, end_year, vertical_input_product, ...
                cesm2_simu_input_vector_cwd, cesm2_simu_input_vector_litter1, cesm2_simu_input_vector_litter2, cesm2_simu_input_vector_litter3, ...
                cesm2_simu_npp, cesm2_simu_altmax, cesm2_simu_altmax_last_year, cesm2_simu_nbedrock, ...
                cesm2_simu_o_scalar, cesm2_simu_n_scalar, cesm2_simu_cellsand, cesm2_simu_soil_temperature, cesm2_simu_w_scalar);
            
            soc_stock_hist_summary(igrid, :, :) = cpool_12c_stock_transient(1:10:111, :);
            eco_process_hist_summary(igrid, :, :) = eco_process_transient(1:10:111, :);
        end
        
        % save simulation record
        if is_server == 0
            save_pathway = ['/Users/ft254/Desktop/cesm2_hist_transient_summary_', ...
                '_mu_', num2str(idiff_mu), '_random_', num2str(irandom), '_', ...
                cesm2_case_name_hat, '_', num2str(start_year), '_',  num2str(end_year), '_cross_valid_', num2str(icross_valid), '_scenario_', num2str(iscenario)];
        else
            save_pathway = [data_path, 'output_data/hist_simu_proda/process_control/cesm2_hist_transient_summary_component_', component_list{icomponent}, '_', ...
                cesm2_case_name_hat, '_', num2str(start_year), '_',  num2str(end_year), '_cross_valid_', num2str(icross_valid), '_scenario_', num2str(iscenario), '_worker_', num2str(worker_start), '_', num2str(worker_end)];
        end
        
        hist_summary.soc_stock_hist_summary = soc_stock_hist_summary;
        hist_summary.eco_process_hist_summary = eco_process_hist_summary;
        
        save([save_pathway, '.mat'], 'hist_summary');
    end
end


disp('Global Projection Finished');
