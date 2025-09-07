function [capacity_baseline, capacity_present, capacity_future] = ...
    fun_scenario_capacity(kp, atmo_14c, simu_year_start, simu_year_end, vertical_input_product, ...
    cesm2_simu_input_vector_cwd, cesm2_simu_input_vector_litter1, cesm2_simu_input_vector_litter2, cesm2_simu_input_vector_litter3, ...
    cesm2_simu_npp, cesm2_simu_altmax, cesm2_simu_altmax_last_year, cesm2_simu_nbedrock, ...
    cesm2_simu_o_scalar, cesm2_simu_n_scalar, cesm2_simu_cellsand, cesm2_simu_soil_temperature, cesm2_simu_w_scalar)

% constants
ss_forcing_year_num = 20;
warmup_year_size = 5000;

year_size = size(cesm2_simu_npp, 2); % length(simu_year_start : simu_year_end);
% timestep_num = 4;
depth_num = 4;
global kelvin_to_celsius
kelvin_to_celsius = 273.15;

global use_vertsoilc npool npool_vr n_soil_layer days_per_year secspday
use_vertsoilc = 1 ; % whether or not use vertical maxing part
npool = 7;  % number of pools if no vertical
npool_vr = 140; % number of pools if vertical
n_soil_layer = 20;  % number of soil layers
days_per_year = 365;
secspday = 24*60*60;

% time resolution
timestep_num = 12;
days_per_timestep = [31, 30, 31, 28, 31, 30, 31, 31, 30, 31, 30, 31];
% days_per_timestep = [92, 89, 92, 92];
month_num_per_timestep = 12/length(days_per_timestep);



decay_lamda_matrix = diag(0.0001209681*ones(npool_vr, 1)); % unit: yr-1, corresponding half life 5730 yr

% if the simulation is default CLM5
is_default = 0;

use_beta = 1;

normalize_q10_to_century_tfunc = false;

global max_altdepth_cryoturbation max_depth_cryoturb
max_altdepth_cryoturbation = 2;
max_depth_cryoturb = 3;

global dz dz_node zisoi zsoi
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
zisoi_0 = 0;
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

% define parameter names
% diffusion (bioturbation) 10^(-4) (m2/yr)
bio = kp(1)*(5*10^(-4) - 3*10^(-5)) + 3*10^(-5);
% cryoturbation 5*10^(-4) (m2/yr)
cryo = kp(2)*(16*10^(-4) - 3*10^(-5)) + 3*10^(-5);
% Q10 (unitless) 1.5
q10 = kp(3)*(3 - 1.2) + 1.2;
% Q10 when forzen (unitless) 1.5
fq10 = q10;
% parameters used in vertical discretization of carbon inputs 10 (metre)
efolding = kp(4)*(1 - 0.1) + 0.1; %proda: (1 - 0) + 0;
% turnover time of CWD (yr) 3.3333
tau4cwd = kp(5)*(6 - 1) + 1;
% tau for metabolic litter (yr) 0.0541
tau4l1 = kp(6)*(0.11 - 0.0001) + 0.0001; % proda *(0.11 - 0) + 0;
% tau for cellulose litter (yr) 0.2041
tau4l2 = kp(7)*(0.3 - 0.1) + 0.1;
% tau for lignin litter (yr)
tau4l3 = tau4l2;
% tau for fast SOC (yr) 0.1370
tau4s1 = kp(8)*(0.5 - 0.0001) + 0.0001; % proda: *(1 - 0) + 0;
% tau for slow SOC (yr) 5
tau4s2 = kp(9)*(10 - 1) + 1; % proda: *(50 - 1) + 1;
% tau for passive SOC (yr) 222.222
tau4s3 = kp(10)*(400 - 20) + 20; % proda *(1000 - 200) + 200;

% fraction from l1 to s2, 0.45
fl1s1 = kp(11)*(0.8 - 0.1) + 0.1;
% fraction from l2 to s1, 0.5
fl2s1 = kp(12)*(0.8 - 0.2) + 0.2;
% fraction from l3 to s2, 0.5
fl3s2 = kp(13)*(0.8 - 0.2) + 0.2;
% fraction from s1 to s2, sand dependeted
fs1s2 = kp(14)*(0.4 - 0.0001) + 0.0001; %proda *(0.4 - 0) + 0;
% fraction from s1 to s3, sand dependeted
fs1s3 = kp(15)*(0.1 - 0.0001) + 0.0001;% proda *(0.05 - 0) + 0;
% fraction from s2 to s1, 0.42
fs2s1 = kp(16)*(0.74 - 0.1) + 0.1;
% fraction from s2 to s3, 0.03
fs2s3 = kp(17)*(0.1 - 0.0001) + 0.0001;% proda *(0.1 - 0) + 0;
% fraction from s3 to s1, 0.45
fs3s1 = kp(18)*(0.9 - 0.0001) + 0.0001;% proda *(0.9 - 0) + 0;
% fraction from cwd to l2, 0.76
fcwdl2 = kp(19)*(1 - 0.5) + 0.5;

% water scaling factor
w_scaling = kp(20)*(5 - 0.0001) + 0.0001;% proda *(5 - 0) + 0;

% beta to describe the shape of vertical profile
% beta = 0.95;
beta = kp(21)*(0.9999 - 0.5) + 0.5;

% maximum and minimum water potential (MPa)
maxpsi= -0.0020;

minpsi= -2; % minimum water potential (MPa)

adv = 0; % parameter for advection (m/yr)

% lag_14c = kp(22)*(10 - 0) + 0;

%% steady state of SOC with baseline forcing

% steady state simu results
nbedrock = nan(timestep_num, 1);
sand_vector = nan(n_soil_layer, timestep_num);
npp_mean = nan(timestep_num, 1);
input_vector_cwd = nan(n_soil_layer, timestep_num);
input_vector_litter1 = nan(n_soil_layer, timestep_num);
input_vector_litter2 = nan(n_soil_layer, timestep_num);
input_vector_litter3 = nan(n_soil_layer, timestep_num);
altmax_current_profile = nan(timestep_num, 1);
altmax_lastyear_profile = nan(timestep_num, 1);
soil_temp_profile = nan(n_soil_layer, timestep_num);
soil_water_profile = nan(n_soil_layer, timestep_num);
xio = nan(n_soil_layer, timestep_num);
xin = nan(n_soil_layer, timestep_num);

forcing_loc_baseline = 1:ss_forcing_year_num;
nbedrock_monthly = mean(cesm2_simu_nbedrock(:, forcing_loc_baseline), 2);
sand_vector_monthly = mean(cesm2_simu_cellsand(:, :, forcing_loc_baseline), 3);
npp_mean_monthly = mean(cesm2_simu_npp(:, forcing_loc_baseline), 2);
input_vector_cwd_monthly = mean(cesm2_simu_input_vector_cwd(:, :, forcing_loc_baseline), 3);
input_vector_litter1_monthly = mean(cesm2_simu_input_vector_litter1(:, :, forcing_loc_baseline), 3);
input_vector_litter2_monthly = mean(cesm2_simu_input_vector_litter2(:, :, forcing_loc_baseline), 3);
input_vector_litter3_monthly = mean(cesm2_simu_input_vector_litter3(:, :, forcing_loc_baseline), 3);
altmax_current_profile_monthly = mean(cesm2_simu_altmax(:, forcing_loc_baseline), 2);
altmax_lastyear_profile_monthly = mean(cesm2_simu_altmax_last_year(:, forcing_loc_baseline), 2);
soil_temp_profile_monthly = mean(cesm2_simu_soil_temperature(:, :, forcing_loc_baseline), 3);
soil_water_profile_monthly = mean(cesm2_simu_w_scalar(:, :, forcing_loc_baseline), 3);
xio_monthly = mean(cesm2_simu_o_scalar(:, :, forcing_loc_baseline), 3);
xin_monthly = mean(cesm2_simu_n_scalar(:, :, forcing_loc_baseline), 3);

for iseason = 1:timestep_num
    nbedrock(iseason) = mean(nbedrock_monthly((iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep));
    sand_vector(:, iseason) = mean(sand_vector_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    npp_mean(iseason) = mean(npp_mean_monthly((iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep));
    input_vector_cwd(:, iseason) = sum(input_vector_cwd_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    input_vector_litter1(:, iseason) = sum(input_vector_litter1_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    input_vector_litter2(:, iseason) = sum(input_vector_litter2_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    input_vector_litter3(:, iseason) = sum(input_vector_litter3_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    altmax_current_profile(iseason) = mean(altmax_current_profile_monthly((iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep));
    altmax_lastyear_profile(iseason) = mean(altmax_lastyear_profile_monthly((iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep));
    soil_temp_profile(:, iseason) = mean(soil_temp_profile_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    soil_water_profile(:, iseason) = mean(soil_water_profile_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    xio(:, iseason) = mean(xio_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    xin(:, iseason) = mean(xin_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
end

npp_mean_baseline = npp_mean;
%----------------------------------------
% Environmental Scalar (Xi)
%----------------------------------------
xit = nan(n_soil_layer, timestep_num);
xiw = nan(n_soil_layer, timestep_num);

for iseason = 1:timestep_num
    % temperature related function xit
    % calculate rate constant scalar for soil temperature
    % assuming that the base rate constants are assigned for non-moisture
    % limiting conditions at 25 C.
    for ilayer = 1 : n_soil_layer
        if soil_temp_profile(ilayer, iseason) >= 0 + kelvin_to_celsius
            xit(ilayer, iseason) = q10^((soil_temp_profile(ilayer, iseason) - (kelvin_to_celsius + 25))/10);
        else
            xit(ilayer, iseason) = q10^((273.15 - 298.15)/10)*(fq10^((soil_temp_profile(ilayer, iseason) - (0 + kelvin_to_celsius))/10));
        end
    end

    catanf_30 = catanf(30);
    normalization_tref = 15;
    if normalize_q10_to_century_tfunc == true
        % scale all decomposition rates by a constant to compensate for offset between original CENTURY temp func and Q10
        normalization_factor = (catanf(normalization_tref)/catanf_30) / (q10^((normalization_tref-25)/10));
        xit = xit*normalization_factor;
    end

    % water related function xiw
    % calculate the rate constant scalar for soil water content.
    % Uses the log relationship with water potential given in
    % Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
    % a comparison of models. Ecology, 68(5):1190-1200.
    % and supported by data in
    % Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
    % and soil moisture. Soil Biol. Biochem., 15(4):447-453.

    %     for ilayer = 1 : n_soil_layer
    %         psi = min(maxpsi, soil_water_profile(ilayer, imonth));
    %
    %         if psi > minpsi
    %             xiw(ilayer, imonth) = (log(minpsi/psi))/(log(minpsi/maxpsi));
    %         else
    %             xiw(ilayer, imonth) = 0;
    %         end
    %         xiw(ilayer, imonth) = max(xiw(ilayer, imonth), 0.05);
    %     end

end

% xit = soil_temp_profile;

xiw = soil_water_profile*w_scaling;
xiw(xiw > 1) = 1;

xiw_baseline = xiw;
%----------------------------------------
% Triangle Matrix, A Matrix and K Matrix
%----------------------------------------
sand_vector = mean(sand_vector, 2, 'omitnan');
% Allocation matrix
a_ma = a_matrix(is_default, fl1s1, fl2s1, fl3s2, fs1s2, fs1s3, fs2s1, fs2s3, fs3s1, fcwdl2, sand_vector);

kk_ma_middle = nan(npool_vr, npool_vr, timestep_num);
tri_ma_middle = nan(npool_vr, npool_vr, timestep_num);
for iseason = 1:timestep_num
    % decomposition matrix
    monthly_xit = xit(:, iseason);
    monthly_xiw = xiw(:, iseason);
    monthly_xio = xio(:, iseason);
    monthly_xin = xin(:, iseason);
    kk_ma_middle(:, :, iseason) = kk_matrix(monthly_xit, monthly_xiw, monthly_xio, monthly_xin, efolding, tau4cwd, tau4l1, tau4l2, tau4l3, tau4s1, tau4s2, tau4s3);
    % tri matrix
    monthly_nbedrock = nbedrock(iseason);
    monthly_altmax_current_profile = altmax_current_profile(iseason);
    monthly_altmax_lastyear_profile = altmax_lastyear_profile(iseason);
    tri_ma_middle(:, :, iseason) = tri_matrix(monthly_nbedrock, monthly_altmax_current_profile, monthly_altmax_lastyear_profile, bio, adv, cryo);
end
tri_ma = mean(tri_ma_middle, 3, 'omitnan');
kk_ma = mean(kk_ma_middle, 3, 'omitnan');

%----------------------------------------
% Vertical Profile
%----------------------------------------
% in the original beta model in Jackson et al 1996, the unit for the depth
% of the soil is cm (dmax*100)
if isnan(vertical_input_product(1)) == 1
    m_to_cm = 100;
    vertical_prof = nan(n_soil_layer, 1);

    if mean(altmax_lastyear_profile) > 0
        for j = 1:n_soil_layer
            if j == 1
                vertical_prof(j) = (beta^((zisoi_0)*m_to_cm) - beta^(zisoi(j)*m_to_cm))/dz(j);
            else
                vertical_prof(j) = (beta^((zisoi(j - 1))*m_to_cm) - beta^(zisoi(j)*m_to_cm))/dz(j);
            end
        end
    else
        vertical_prof(1) = 1/dz(1);
        vertical_prof(2:end) = 0;
    end

    vertical_input = dz(1:n_soil_layer).*vertical_prof/sum(vertical_prof.*dz(1:n_soil_layer));
else
    vertical_input = vertical_input_product;
end
%----------------------------------------
% Analytical Solution of SOC
%----------------------------------------
matrix_in=nan(140,1);

% calculate annual carbon input, sum monthly input (gc/m3/month) as annual input (gc/m3/year)
input_tot_cwd = input_vector_cwd;
input_tot_litter1 = input_vector_litter1;
input_tot_litter2 = input_vector_litter2;
input_tot_litter3 = input_vector_litter3;

for iseason = 1:timestep_num
    input_tot_litter1(:, iseason) = input_vector_litter1(:, iseason).*dz(1:n_soil_layer);
    input_tot_litter2(:, iseason) = input_vector_litter2(:, iseason).*dz(1:n_soil_layer);
    input_tot_litter3(:, iseason) = input_vector_litter3(:, iseason).*dz(1:n_soil_layer);
    input_tot_cwd(:, iseason) = input_vector_cwd(:, iseason).*dz(1:n_soil_layer);
end

input_tot_litter1 = sum(sum(input_tot_litter1));
input_tot_litter2 = sum(sum(input_tot_litter2));
input_tot_litter3 = sum(sum(input_tot_litter3));
input_tot_cwd = sum(sum(input_tot_cwd));

matrix_in(1:20,1) = input_tot_cwd*vertical_input./(dz(1:n_soil_layer)*days_per_year); % litter input gc/m3/day
matrix_in(21:40,1) = input_tot_litter1*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(41:60,1) = input_tot_litter2*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(61:80,1) = input_tot_litter3*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(81:140,1) = 0;

cpool = (a_ma*kk_ma-tri_ma)\(-matrix_in);
cpool(cpool < 0) = 0;

cpool = [cpool(1:20), cpool(21:40), cpool(41:60), cpool(61:80), cpool(81:100), cpool(101:120), cpool(121:140)];

capacity_baseline_30cm = sum(cpool(1:4, :).*repmat(dz(1:4, 1), [1, npool]), 1) ...
    + cpool(5, :)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4));
capacity_baseline_100cm = sum(cpool(1:8, :).*repmat(dz(1:8, 1), [1, npool]), 1) ...
    + cpool(9, :)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8));
capacity_baseline = [capacity_baseline_30cm; capacity_baseline_100cm];

%% steady state of SOC with hist forcing

% steady state simu results
nbedrock = nan(timestep_num, 1);
sand_vector = nan(n_soil_layer, timestep_num);
npp_mean = nan(timestep_num, 1);
input_vector_cwd = nan(n_soil_layer, timestep_num);
input_vector_litter1 = nan(n_soil_layer, timestep_num);
input_vector_litter2 = nan(n_soil_layer, timestep_num);
input_vector_litter3 = nan(n_soil_layer, timestep_num);
altmax_current_profile = nan(timestep_num, 1);
altmax_lastyear_profile = nan(timestep_num, 1);
soil_temp_profile = nan(n_soil_layer, timestep_num);
soil_water_profile = nan(n_soil_layer, timestep_num);
xio = nan(n_soil_layer, timestep_num);
xin = nan(n_soil_layer, timestep_num);

forcing_loc_present = 92:111; % year_size-ss_forcing_year_num+1:year_size;
nbedrock_monthly = mean(cesm2_simu_nbedrock(:, forcing_loc_present), 2);
sand_vector_monthly = mean(cesm2_simu_cellsand(:, :, forcing_loc_present), 3);
npp_mean_monthly = mean(cesm2_simu_npp(:, forcing_loc_present), 2);
input_vector_cwd_monthly = mean(cesm2_simu_input_vector_cwd(:, :, forcing_loc_present), 3);
input_vector_litter1_monthly = mean(cesm2_simu_input_vector_litter1(:, :, forcing_loc_present), 3);
input_vector_litter2_monthly = mean(cesm2_simu_input_vector_litter2(:, :, forcing_loc_present), 3);
input_vector_litter3_monthly = mean(cesm2_simu_input_vector_litter3(:, :, forcing_loc_present), 3);
altmax_current_profile_monthly = mean(cesm2_simu_altmax(:, forcing_loc_present), 2);
altmax_lastyear_profile_monthly = mean(cesm2_simu_altmax_last_year(:, forcing_loc_present), 2);
soil_temp_profile_monthly = mean(cesm2_simu_soil_temperature(:, :, forcing_loc_present), 3);
soil_water_profile_monthly = mean(cesm2_simu_w_scalar(:, :, forcing_loc_present), 3);
xio_monthly = mean(cesm2_simu_o_scalar(:, :, forcing_loc_present), 3);
xin_monthly = mean(cesm2_simu_n_scalar(:, :, forcing_loc_present), 3);

for iseason = 1:timestep_num
    nbedrock(iseason) = mean(nbedrock_monthly((iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep));
    sand_vector(:, iseason) = mean(sand_vector_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    npp_mean(iseason) = mean(npp_mean_monthly((iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep));
    input_vector_cwd(:, iseason) = sum(input_vector_cwd_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    input_vector_litter1(:, iseason) = sum(input_vector_litter1_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    input_vector_litter2(:, iseason) = sum(input_vector_litter2_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    input_vector_litter3(:, iseason) = sum(input_vector_litter3_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    altmax_current_profile(iseason) = mean(altmax_current_profile_monthly((iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep));
    altmax_lastyear_profile(iseason) = mean(altmax_lastyear_profile_monthly((iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep));
    soil_temp_profile(:, iseason) = mean(soil_temp_profile_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    soil_water_profile(:, iseason) = mean(soil_water_profile_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    xio(:, iseason) = mean(xio_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    xin(:, iseason) = mean(xin_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
end

%----------------------------------------
% Environmental Scalar (Xi)
%----------------------------------------
xit = nan(n_soil_layer, timestep_num);
xiw = nan(n_soil_layer, timestep_num);

for iseason = 1:timestep_num
    % temperature related function xit
    % calculate rate constant scalar for soil temperature
    % assuming that the base rate constants are assigned for non-moisture
    % limiting conditions at 25 C.
    for ilayer = 1 : n_soil_layer
        if soil_temp_profile(ilayer, iseason) >= 0 + kelvin_to_celsius
            xit(ilayer, iseason) = q10^((soil_temp_profile(ilayer, iseason) - (kelvin_to_celsius + 25))/10);
        else
            xit(ilayer, iseason) = q10^((273.15 - 298.15)/10)*(fq10^((soil_temp_profile(ilayer, iseason) - (0 + kelvin_to_celsius))/10));
        end
    end

    catanf_30 = catanf(30);
    normalization_tref = 15;
    if normalize_q10_to_century_tfunc == true
        % scale all decomposition rates by a constant to compensate for offset between original CENTURY temp func and Q10
        normalization_factor = (catanf(normalization_tref)/catanf_30) / (q10^((normalization_tref-25)/10));
        xit = xit*normalization_factor;
    end

    % water related function xiw
    % calculate the rate constant scalar for soil water content.
    % Uses the log relationship with water potential given in
    % Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
    % a comparison of models. Ecology, 68(5):1190-1200.
    % and supported by data in
    % Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
    % and soil moisture. Soil Biol. Biochem., 15(4):447-453.

    %     for ilayer = 1 : n_soil_layer
    %         psi = min(maxpsi, soil_water_profile(ilayer, imonth));
    %
    %         if psi > minpsi
    %             xiw(ilayer, imonth) = (log(minpsi/psi))/(log(minpsi/maxpsi));
    %         else
    %             xiw(ilayer, imonth) = 0;
    %         end
    %         xiw(ilayer, imonth) = max(xiw(ilayer, imonth), 0.05);
    %     end

end

% xit = soil_temp_profile;

xiw = soil_water_profile*w_scaling;
xiw(xiw > 1) = 1;
xiw_present = xiw;
%----------------------------------------
% Triangle Matrix, A Matrix and K Matrix
%----------------------------------------
sand_vector = mean(sand_vector, 2, 'omitnan');
% Allocation matrix
a_ma = a_matrix(is_default, fl1s1, fl2s1, fl3s2, fs1s2, fs1s3, fs2s1, fs2s3, fs3s1, fcwdl2, sand_vector);

kk_ma_middle = nan(npool_vr, npool_vr, timestep_num);
tri_ma_middle = nan(npool_vr, npool_vr, timestep_num);
for iseason = 1:timestep_num
    % decomposition matrix
    monthly_xit = xit(:, iseason);
    monthly_xiw = xiw(:, iseason);
    monthly_xio = xio(:, iseason);
    monthly_xin = xin(:, iseason);
    kk_ma_middle(:, :, iseason) = kk_matrix(monthly_xit, monthly_xiw, monthly_xio, monthly_xin, efolding, tau4cwd, tau4l1, tau4l2, tau4l3, tau4s1, tau4s2, tau4s3);
    % tri matrix
    monthly_nbedrock = nbedrock(iseason);
    monthly_altmax_current_profile = altmax_current_profile(iseason);
    monthly_altmax_lastyear_profile = altmax_lastyear_profile(iseason);
    tri_ma_middle(:, :, iseason) = tri_matrix(monthly_nbedrock, monthly_altmax_current_profile, monthly_altmax_lastyear_profile, bio, adv, cryo);
end
tri_ma = mean(tri_ma_middle, 3, 'omitnan');
kk_ma = mean(kk_ma_middle, 3, 'omitnan');

%----------------------------------------
% Vertical Profile
%----------------------------------------
% in the original beta model in Jackson et al 1996, the unit for the depth
% of the soil is cm (dmax*100)
if isnan(vertical_input_product(1)) == 1
    m_to_cm = 100;
    vertical_prof = nan(n_soil_layer, 1);

    if mean(altmax_lastyear_profile) > 0
        for j = 1:n_soil_layer
            if j == 1
                vertical_prof(j) = (beta^((zisoi_0)*m_to_cm) - beta^(zisoi(j)*m_to_cm))/dz(j);
            else
                vertical_prof(j) = (beta^((zisoi(j - 1))*m_to_cm) - beta^(zisoi(j)*m_to_cm))/dz(j);
            end
        end
    else
        vertical_prof(1) = 1/dz(1);
        vertical_prof(2:end) = 0;
    end

    vertical_input = dz(1:n_soil_layer).*vertical_prof/sum(vertical_prof.*dz(1:n_soil_layer));
else
    vertical_input = vertical_input_product;
end
%----------------------------------------
% Analytical Solution of SOC
%----------------------------------------
matrix_in=nan(140,1);

% calculate annual carbon input, sum monthly input (gc/m3/month) as annual input (gc/m3/year)
input_tot_cwd = input_vector_cwd;
input_tot_litter1 = input_vector_litter1;
input_tot_litter2 = input_vector_litter2;
input_tot_litter3 = input_vector_litter3;

for iseason = 1:timestep_num
    input_tot_litter1(:, iseason) = input_vector_litter1(:, iseason).*dz(1:n_soil_layer);
    input_tot_litter2(:, iseason) = input_vector_litter2(:, iseason).*dz(1:n_soil_layer);
    input_tot_litter3(:, iseason) = input_vector_litter3(:, iseason).*dz(1:n_soil_layer);
    input_tot_cwd(:, iseason) = input_vector_cwd(:, iseason).*dz(1:n_soil_layer);
end

input_tot_litter1 = sum(sum(input_tot_litter1));
input_tot_litter2 = sum(sum(input_tot_litter2));
input_tot_litter3 = sum(sum(input_tot_litter3));
input_tot_cwd = sum(sum(input_tot_cwd));

matrix_in(1:20,1) = input_tot_cwd*vertical_input./(dz(1:n_soil_layer)*days_per_year); % litter input gc/m3/day
matrix_in(21:40,1) = input_tot_litter1*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(41:60,1) = input_tot_litter2*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(61:80,1) = input_tot_litter3*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(81:140,1) = 0;

cpool = (a_ma*kk_ma-tri_ma)\(-matrix_in);
cpool(cpool < 0) = 0;

cpool = [cpool(1:20), cpool(21:40), cpool(41:60), cpool(61:80), cpool(81:100), cpool(101:120), cpool(121:140)];
capacity_present_30cm = sum(cpool(1:4, :).*repmat(dz(1:4, 1), [1, npool]), 1) ...
    + cpool(5, :)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4));
capacity_present_100cm = sum(cpool(1:8, :).*repmat(dz(1:8, 1), [1, npool]), 1) ...
    + cpool(9, :)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8));

capacity_present = [capacity_present_30cm; capacity_present_100cm];


%% steady state of SOC with future forcing
% steady state simu results
nbedrock = nan(timestep_num, 1);
sand_vector = nan(n_soil_layer, timestep_num);
npp_mean = nan(timestep_num, 1);
input_vector_cwd = nan(n_soil_layer, timestep_num);
input_vector_litter1 = nan(n_soil_layer, timestep_num);
input_vector_litter2 = nan(n_soil_layer, timestep_num);
input_vector_litter3 = nan(n_soil_layer, timestep_num);
altmax_current_profile = nan(timestep_num, 1);
altmax_lastyear_profile = nan(timestep_num, 1);
soil_temp_profile = nan(n_soil_layer, timestep_num);
soil_water_profile = nan(n_soil_layer, timestep_num);
xio = nan(n_soil_layer, timestep_num);
xin = nan(n_soil_layer, timestep_num);

forcing_loc_future = year_size-ss_forcing_year_num+1:year_size;
nbedrock_monthly = mean(cesm2_simu_nbedrock(:, forcing_loc_future), 2);
sand_vector_monthly = mean(cesm2_simu_cellsand(:, :, forcing_loc_future), 3);
npp_mean_monthly = mean(cesm2_simu_npp(:, forcing_loc_future), 2);
input_vector_cwd_monthly = mean(cesm2_simu_input_vector_cwd(:, :, forcing_loc_baseline), 3);
input_vector_litter1_monthly = mean(cesm2_simu_input_vector_litter1(:, :, forcing_loc_baseline), 3);
input_vector_litter2_monthly = mean(cesm2_simu_input_vector_litter2(:, :, forcing_loc_baseline), 3);
input_vector_litter3_monthly = mean(cesm2_simu_input_vector_litter3(:, :, forcing_loc_baseline), 3);
altmax_current_profile_monthly = mean(cesm2_simu_altmax(:, forcing_loc_future), 2);
altmax_lastyear_profile_monthly = mean(cesm2_simu_altmax_last_year(:, forcing_loc_future), 2);
soil_temp_profile_monthly = mean(cesm2_simu_soil_temperature(:, :, forcing_loc_future), 3);
soil_water_profile_monthly = mean(cesm2_simu_w_scalar(:, :, forcing_loc_future), 3);
xio_monthly = mean(cesm2_simu_o_scalar(:, :, forcing_loc_future), 3);
xin_monthly = mean(cesm2_simu_n_scalar(:, :, forcing_loc_future), 3);

for iseason = 1:timestep_num
    nbedrock(iseason) = mean(nbedrock_monthly((iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep));
    sand_vector(:, iseason) = mean(sand_vector_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    npp_mean(iseason) = mean(npp_mean_monthly((iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep));
    input_vector_cwd(:, iseason) = sum(input_vector_cwd_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    input_vector_litter1(:, iseason) = sum(input_vector_litter1_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    input_vector_litter2(:, iseason) = sum(input_vector_litter2_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    input_vector_litter3(:, iseason) = sum(input_vector_litter3_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    altmax_current_profile(iseason) = mean(altmax_current_profile_monthly((iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep));
    altmax_lastyear_profile(iseason) = mean(altmax_lastyear_profile_monthly((iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep));
    soil_temp_profile(:, iseason) = mean(soil_temp_profile_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    soil_water_profile(:, iseason) = mean(soil_water_profile_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    xio(:, iseason) = mean(xio_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
    xin(:, iseason) = mean(xin_monthly(:, (iseason-1)*month_num_per_timestep+1:iseason*month_num_per_timestep), 2);
end

input_scaling_factor = sum(npp_mean'.*days_per_timestep*3600*24)/sum(npp_mean_baseline'.*days_per_timestep*3600*24);

%----------------------------------------
% Environmental Scalar (Xi)
%----------------------------------------
xit = nan(n_soil_layer, timestep_num);
xiw = nan(n_soil_layer, timestep_num);

for iseason = 1:timestep_num
    % temperature related function xit
    % calculate rate constant scalar for soil temperature
    % assuming that the base rate constants are assigned for non-moisture
    % limiting conditions at 25 C.
    for ilayer = 1 : n_soil_layer
        if soil_temp_profile(ilayer, iseason) >= 0 + kelvin_to_celsius
            xit(ilayer, iseason) = q10^((soil_temp_profile(ilayer, iseason) - (kelvin_to_celsius + 25))/10);
        else
            xit(ilayer, iseason) = q10^((273.15 - 298.15)/10)*(fq10^((soil_temp_profile(ilayer, iseason) - (0 + kelvin_to_celsius))/10));
        end
    end

    catanf_30 = catanf(30);
    normalization_tref = 15;
    if normalize_q10_to_century_tfunc == true
        % scale all decomposition rates by a constant to compensate for offset between original CENTURY temp func and Q10
        normalization_factor = (catanf(normalization_tref)/catanf_30) / (q10^((normalization_tref-25)/10));
        xit = xit*normalization_factor;
    end

    % water related function xiw
    % calculate the rate constant scalar for soil water content.
    % Uses the log relationship with water potential given in
    % Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
    % a comparison of models. Ecology, 68(5):1190-1200.
    % and supported by data in
    % Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
    % and soil moisture. Soil Biol. Biochem., 15(4):447-453.

    %     for ilayer = 1 : n_soil_layer
    %         psi = min(maxpsi, soil_water_profile(ilayer, imonth));
    %
    %         if psi > minpsi
    %             xiw(ilayer, imonth) = (log(minpsi/psi))/(log(minpsi/maxpsi));
    %         else
    %             xiw(ilayer, imonth) = 0;
    %         end
    %         xiw(ilayer, imonth) = max(xiw(ilayer, imonth), 0.05);
    %     end

end

% xit = soil_temp_profile;

% xiw = soil_water_profile*w_scaling;
xiw = (xiw_present./xiw_baseline).*xiw_baseline;
xiw(xiw > 1) = 1;

%----------------------------------------
% Triangle Matrix, A Matrix and K Matrix
%----------------------------------------
sand_vector = mean(sand_vector, 2, 'omitnan');
% Allocation matrix
a_ma = a_matrix(is_default, fl1s1, fl2s1, fl3s2, fs1s2, fs1s3, fs2s1, fs2s3, fs3s1, fcwdl2, sand_vector);

kk_ma_middle = nan(npool_vr, npool_vr, timestep_num);
tri_ma_middle = nan(npool_vr, npool_vr, timestep_num);
for iseason = 1:timestep_num
    % decomposition matrix
    monthly_xit = xit(:, iseason);
    monthly_xiw = xiw(:, iseason);
    monthly_xio = xio(:, iseason);
    monthly_xin = xin(:, iseason);
    kk_ma_middle(:, :, iseason) = kk_matrix(monthly_xit, monthly_xiw, monthly_xio, monthly_xin, efolding, tau4cwd, tau4l1, tau4l2, tau4l3, tau4s1, tau4s2, tau4s3);
    % tri matrix
    monthly_nbedrock = nbedrock(iseason);
    monthly_altmax_current_profile = altmax_current_profile(iseason);
    monthly_altmax_lastyear_profile = altmax_lastyear_profile(iseason);
    tri_ma_middle(:, :, iseason) = tri_matrix(monthly_nbedrock, monthly_altmax_current_profile, monthly_altmax_lastyear_profile, bio, adv, cryo);
end
tri_ma = mean(tri_ma_middle, 3, 'omitnan');
kk_ma = mean(kk_ma_middle, 3, 'omitnan');

%----------------------------------------
% Vertical Profile
%----------------------------------------
% in the original beta model in Jackson et al 1996, the unit for the depth
% of the soil is cm (dmax*100)
if isnan(vertical_input_product(1)) == 1
    m_to_cm = 100;
    vertical_prof = nan(n_soil_layer, 1);

    if mean(altmax_lastyear_profile) > 0
        for j = 1:n_soil_layer
            if j == 1
                vertical_prof(j) = (beta^((zisoi_0)*m_to_cm) - beta^(zisoi(j)*m_to_cm))/dz(j);
            else
                vertical_prof(j) = (beta^((zisoi(j - 1))*m_to_cm) - beta^(zisoi(j)*m_to_cm))/dz(j);
            end
        end
    else
        vertical_prof(1) = 1/dz(1);
        vertical_prof(2:end) = 0;
    end

    vertical_input = dz(1:n_soil_layer).*vertical_prof/sum(vertical_prof.*dz(1:n_soil_layer));
else
    vertical_input = vertical_input_product;
end
%----------------------------------------
% Analytical Solution of SOC
%----------------------------------------
matrix_in=nan(140,1);

% calculate annual carbon input, sum monthly input (gc/m3/month) as annual input (gc/m3/year)
input_tot_cwd = input_vector_cwd;
input_tot_litter1 = input_vector_litter1;
input_tot_litter2 = input_vector_litter2;
input_tot_litter3 = input_vector_litter3;

for iseason = 1:timestep_num
    input_tot_litter1(:, iseason) = input_vector_litter1(:, iseason).*dz(1:n_soil_layer);
    input_tot_litter2(:, iseason) = input_vector_litter2(:, iseason).*dz(1:n_soil_layer);
    input_tot_litter3(:, iseason) = input_vector_litter3(:, iseason).*dz(1:n_soil_layer);
    input_tot_cwd(:, iseason) = input_vector_cwd(:, iseason).*dz(1:n_soil_layer);
end

input_tot_litter1 = sum(sum(input_tot_litter1));
input_tot_litter2 = sum(sum(input_tot_litter2));
input_tot_litter3 = sum(sum(input_tot_litter3));
input_tot_cwd = sum(sum(input_tot_cwd));

matrix_in(1:20,1) = input_tot_cwd*vertical_input./(dz(1:n_soil_layer)*days_per_year); % litter input gc/m3/day
matrix_in(21:40,1) = input_tot_litter1*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(41:60,1) = input_tot_litter2*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(61:80,1) = input_tot_litter3*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(81:140,1) = 0;

cpool = (a_ma*kk_ma-tri_ma)\(-matrix_in*input_scaling_factor);
cpool(cpool < 0) = 0;

cpool = [cpool(1:20), cpool(21:40), cpool(41:60), cpool(61:80), cpool(81:100), cpool(101:120), cpool(121:140)];
capacity_future_30cm = sum(cpool(1:4, :).*repmat(dz(1:4, 1), [1, npool]), 1) ...
    + cpool(5, :)*dz(5)*(0.3 - zisoi(4))/(zisoi(5) - zisoi(4));
capacity_future_100cm = sum(cpool(1:8, :).*repmat(dz(1:8, 1), [1, npool]), 1) ...
    + cpool(9, :)*dz(9)*(1 - zisoi(8))/(zisoi(9) - zisoi(8));

capacity_future = [capacity_future_30cm; capacity_future_100cm];

end  % function

% ----- CENTURY T response function
function catanf_results = catanf(t1)
catanf_results = 11.75 +(29.7 / pi) * atan( pi * 0.031  * ( t1 - 15.4 ));
end

