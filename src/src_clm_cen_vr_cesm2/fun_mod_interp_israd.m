function [interp_14c_layer] = ...
    fun_mod_interp_israd(zsoi, soil_decom_num, ...
    mod_14c_layer, obs_14c_layer_time, obs_14c_layer_depth)

delta_14c_time_series = 1900:2010;
%----------------------------------------
% 14c layer interpolation
%----------------------------------------
% soil layer time loc
soil_layer_time_loc = nan(length(obs_14c_layer_depth), 1);
for ilayer = 1:length(obs_14c_layer_depth)
    time_diff_abs = abs(delta_14c_time_series - obs_14c_layer_time(ilayer));
    loc_min = find(time_diff_abs == min(time_diff_abs));
    soil_layer_time_loc(ilayer) = loc_min(1);
end

% soc layer interplation
interp_14c_layer = nan(length(obs_14c_layer_depth), 1);
for ilayer = 1:length(obs_14c_layer_depth)
    layer_time_loc = soil_layer_time_loc(ilayer);
    mod_14c_layer_middle = mod_14c_layer(layer_time_loc, :)';
    interp_14c_layer(ilayer) = interp1(zsoi(1:soil_decom_num, 1), mod_14c_layer_middle, obs_14c_layer_depth(ilayer), 'pchip');
end


end
