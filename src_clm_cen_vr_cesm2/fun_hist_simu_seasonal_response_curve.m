function [cpool_12c_stock_transient, cpool_12c_resp_transient, ...
    cpool_14c_stock_transient, coop_14c_resp_transient] = ...
    fun_hist_simu_seasonal_response_curve(kp, atmo_14c, simu_year_start, simu_year_end, grid_npp_allocation, ...
    cesm2_simu_input_vector_cwd, cesm2_simu_input_vector_litter1, cesm2_simu_input_vector_litter2, cesm2_simu_input_vector_litter3, ...
    cesm2_simu_npp, cesm2_simu_altmax, cesm2_simu_altmax_last_year, cesm2_simu_nbedrock, ...
    cesm2_simu_o_scalar, cesm2_simu_n_scalar, cesm2_simu_cellsand, cesm2_simu_soil_temperature, cesm2_simu_w_scalar)

% constants
ss_forcing_year_num = 20;
warmup_year_size = 5000;

year_size = length(simu_year_start : simu_year_end);
month_num = 12;
season_num = 4;
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
% dt = secspday*30;
days_per_month = [31, 30, 31, 28, 31, 30, 31, 31, 30, 31, 30, 31];
days_per_season = [92, 89, 92, 92];

decay_lamda_matrix = diag(0.0001209681*ones(npool_vr, 1)); % unit: yr-1, corresponding half life 5730 yr

% steady state simu results
nbedrock = nan(season_num, 1);
sand_vector = nan(n_soil_layer, season_num);
npp_mean = nan(season_num, 1);
input_vector_cwd = nan(n_soil_layer, season_num);
input_vector_litter1 = nan(n_soil_layer, season_num);
input_vector_litter2 = nan(n_soil_layer, season_num);
input_vector_litter3 = nan(n_soil_layer, season_num);
altmax_current_profile = nan(season_num, 1);
altmax_lastyear_profile = nan(season_num, 1);
soil_temp_profile = nan(n_soil_layer, season_num);
soil_water_profile = nan(n_soil_layer, season_num);
xio = nan(n_soil_layer, season_num);
xin = nan(n_soil_layer, season_num);

nbedrock_monthly = mean(cesm2_simu_nbedrock(:, 1:ss_forcing_year_num), 2);
sand_vector_monthly = mean(cesm2_simu_cellsand(:, :, 1:ss_forcing_year_num), 3);
npp_mean_monthly = mean(cesm2_simu_npp(:, 1:ss_forcing_year_num), 2);
input_vector_cwd_monthly = mean(cesm2_simu_input_vector_cwd(:, :, 1:ss_forcing_year_num), 3);
input_vector_litter1_monthly = mean(cesm2_simu_input_vector_litter1(:, :, 1:ss_forcing_year_num), 3);
input_vector_litter2_monthly = mean(cesm2_simu_input_vector_litter2(:, :, 1:ss_forcing_year_num), 3);
input_vector_litter3_monthly = mean(cesm2_simu_input_vector_litter3(:, :, 1:ss_forcing_year_num), 3);
altmax_current_profile_monthly = mean(cesm2_simu_altmax(:, 1:ss_forcing_year_num), 2);
altmax_lastyear_profile_monthly = mean(cesm2_simu_altmax_last_year(:, 1:ss_forcing_year_num), 2);
soil_temp_profile_monthly = mean(cesm2_simu_soil_temperature(:, :, 1:ss_forcing_year_num), 3);
soil_water_profile_monthly = mean(cesm2_simu_w_scalar(:, :, 1:ss_forcing_year_num), 3);
xio_monthly = mean(cesm2_simu_o_scalar(:, :, 1:ss_forcing_year_num), 3);
xin_monthly = mean(cesm2_simu_n_scalar(:, :, 1:ss_forcing_year_num), 3);

month_num_per_season = 3;
for iseason = 1:season_num
    nbedrock(iseason) = mean(nbedrock_monthly((iseason-1)*month_num_per_season+1:iseason*month_num_per_season));
    sand_vector(:, iseason) = mean(sand_vector_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season), 2);
    npp_mean(iseason) = mean(npp_mean_monthly((iseason-1)*month_num_per_season+1:iseason*month_num_per_season));
    input_vector_cwd(:, iseason) = sum(input_vector_cwd_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season), 2);
    input_vector_litter1(:, iseason) = sum(input_vector_litter1_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season), 2);
    input_vector_litter2(:, iseason) = sum(input_vector_litter2_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season), 2);
    input_vector_litter3(:, iseason) = sum(input_vector_litter3_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season), 2);
    altmax_current_profile(iseason) = mean(altmax_current_profile_monthly((iseason-1)*month_num_per_season+1:iseason*month_num_per_season));
    altmax_lastyear_profile(iseason) = mean(altmax_lastyear_profile_monthly((iseason-1)*month_num_per_season+1:iseason*month_num_per_season));
    soil_temp_profile(:, iseason) = mean(soil_temp_profile_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season), 2);
    soil_water_profile(:, iseason) = mean(soil_water_profile_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season), 2);
    xio(:, iseason) = mean(xio_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season), 2);
    xin(:, iseason) = mean(xin_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season), 2);
end
% transient simu results
nbedrock_transient = nan(season_num, year_size);
sand_vector_transient = nan(n_soil_layer, season_num, year_size);
npp_mean_transient = nan(season_num, year_size);
input_vector_cwd_transient = nan(n_soil_layer, season_num, year_size);
input_vector_litter1_transient = nan(n_soil_layer, season_num, year_size);
input_vector_litter2_transient = nan(n_soil_layer, season_num, year_size);
input_vector_litter3_transient = nan(n_soil_layer, season_num, year_size);
altmax_current_profile_transient = nan(season_num, year_size);
altmax_lastyear_profile_transient = nan(season_num, year_size);
soil_temp_profile_transient = nan(n_soil_layer, season_num, year_size);
soil_water_profile_transient = nan(n_soil_layer, season_num, year_size);
xio_transient = nan(n_soil_layer, season_num, year_size);
xin_transient = nan(n_soil_layer, season_num, year_size);

nbedrock_transient_monthly = cesm2_simu_nbedrock(:, length(1900:simu_year_start):end);
sand_vector_transient_monthly = cesm2_simu_cellsand(:, :, length(1900:simu_year_start):end);
npp_mean_transient_monthly = cesm2_simu_npp(:, length(1900:simu_year_start):end);
input_vector_cwd_transient_monthly = cesm2_simu_input_vector_cwd(:, :, length(1900:simu_year_start):end);
input_vector_litter1_transient_monthly = cesm2_simu_input_vector_litter1(:, :, length(1900:simu_year_start):end);
input_vector_litter2_transient_monthly = cesm2_simu_input_vector_litter2(:, :, length(1900:simu_year_start):end);
input_vector_litter3_transient_monthly = cesm2_simu_input_vector_litter3(:, :, length(1900:simu_year_start):end);
altmax_current_profile_transient_monthly = cesm2_simu_altmax(:, length(1900:simu_year_start):end);
altmax_lastyear_profile_transient_monthly = cesm2_simu_altmax_last_year(:, length(1900:simu_year_start):end);
soil_temp_profile_transient_monthly = cesm2_simu_soil_temperature(:, :, length(1900:simu_year_start):end);
soil_water_profile_transient_monthly = cesm2_simu_w_scalar(:, :, length(1900:simu_year_start):end);
xio_transient_monthly = cesm2_simu_o_scalar(:, :, length(1900:simu_year_start):end);
xin_transient_monthly = cesm2_simu_n_scalar(:, :, length(1900:simu_year_start):end);

for iseason = 1:season_num
    nbedrock_transient(iseason, :) = mean(nbedrock_transient_monthly((iseason-1)*month_num_per_season+1:iseason*month_num_per_season, :), 1);
    sand_vector_transient(:, iseason, :) = mean(sand_vector_transient_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season, :), 2);
    npp_mean_transient(iseason, :) = mean(npp_mean_transient_monthly((iseason-1)*month_num_per_season+1:iseason*month_num_per_season, :), 1);
    input_vector_cwd_transient(:, iseason, :) = sum(input_vector_cwd_transient_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season, :), 2);
    input_vector_litter1_transient(:, iseason, :) = sum(input_vector_litter1_transient_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season, :), 2);
    input_vector_litter2_transient(:, iseason, :) = sum(input_vector_litter2_transient_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season, :), 2);
    input_vector_litter3_transient(:, iseason, :) = sum(input_vector_litter3_transient_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season, :), 2);
    altmax_current_profile_transient(iseason, :) = mean(altmax_current_profile_transient_monthly((iseason-1)*month_num_per_season+1:iseason*month_num_per_season, :), 1);
    altmax_lastyear_profile_transient(iseason, :) = mean(altmax_lastyear_profile_transient_monthly((iseason-1)*month_num_per_season+1:iseason*month_num_per_season, :), 1);
    soil_temp_profile_transient(:, iseason, :) = mean(soil_temp_profile_transient_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season, :), 2);
    soil_water_profile_transient(:, iseason, :) = mean(soil_water_profile_transient_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season, :), 2);
    xio_transient(:, iseason, :) = mean(xio_transient_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season, :), 2);
    xin_transient(:, iseason, :) = mean(xin_transient_monthly(:, (iseason-1)*month_num_per_season+1:iseason*month_num_per_season, :), 2);
end

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
efolding = kp(4)*(1 - 0) + 0;
% turnover time of CWD (yr) 3.3333
tau4cwd = kp(5)*(6 - 1) + 1;
% tau for metabolic litter (yr) 0.0541
tau4l1 = kp(6)*(0.11 - 0) + 0;
% tau for cellulose litter (yr) 0.2041
tau4l2 = kp(7)*(0.3 - 0.1) + 0.1;
% tau for lignin litter (yr)
tau4l3 = tau4l2;
% tau for fast SOC (yr) 0.1370
tau4s1 = kp(8)*(1 - 0) + 0;
% tau for slow SOC (yr) 5
tau4s2 = kp(9)*(50 - 1) + 1;
% tau for passive SOC (yr) 222.222
tau4s3 = kp(10)*(1000 - 200) + 200;

% fraction from l1 to s2, 0.45
fl1s1 = kp(11)*(0.8 - 0.1) + 0.1;
% fraction from l2 to s1, 0.5
fl2s1 = kp(12)*(0.8 - 0.2) + 0.2;
% fraction from l3 to s2, 0.5
fl3s2 = kp(13)*(0.8 - 0.2) + 0.2;
% fraction from s1 to s2, sand dependeted
fs1s2 = kp(14)*(0.4 - 0) + 0;
% fraction from s1 to s3, sand dependeted
fs1s3 = kp(15)*(0.05 - 0) + 0;
% fraction from s2 to s1, 0.42
fs2s1 = kp(16)*(0.74 - 0.1) + 0.1;
% fraction from s2 to s3, 0.03
fs2s3 = kp(17)*(0.1 - 0) + 0;
% fraction from s3 to s1, 0.45
fs3s1 = kp(18)*(0.9 - 0) + 0;
% fraction from cwd to l2, 0.76
fcwdl2 = kp(19)*(1 - 0.5) + 0.5;

% water scaling factor
w_scaling = kp(20)*(5 - 0) + 0;

% beta to describe the shape of vertical profile
% beta = 0.95;
beta = kp(21)*(0.9999 - 0.5) + 0.5;

% maximum and minimum water potential (MPa)
maxpsi= -0.0020;

minpsi= -2; % minimum water potential (MPa)

adv = 0; % parameter for advection (m/yr)

% lag_14c = kp(22)*(10 - 0) + 0;


% plant carbon input allocation
soil_allocation_depth = [0, 10, 30, 50, 70, 90, 125, 175]';

interp_allo = interp1(soil_allocation_depth/100, grid_npp_allocation, zsoi(1:n_soil_layer, 1), 'linear');

interp_allo(interp_allo < 0) = 0;
interp_allo(isnan(interp_allo)) = 0;
vertical_input_product = interp_allo/sum(interp_allo);

%% steady state of SOC
%----------------------------------------
% Environmental Scalar (Xi)
%----------------------------------------
xit = nan(n_soil_layer, season_num);
xiw = nan(n_soil_layer, season_num);

for iseason = 1:season_num
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

%----------------------------------------
% Triangle Matrix, A Matrix and K Matrix
%----------------------------------------
sand_vector = mean(sand_vector, 2, 'omitnan');
% Allocation matrix
a_ma = a_matrix(is_default, fl1s1, fl2s1, fl3s2, fs1s2, fs1s3, fs2s1, fs2s3, fs3s1, fcwdl2, sand_vector);

kk_ma_middle = nan(npool_vr, npool_vr, season_num);
tri_ma_middle = nan(npool_vr, npool_vr, season_num);
for iseason = 1:season_num
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
if isnan(grid_npp_allocation(1)) == 1
    m_to_cm = 100;
    vertical_prof = nan(n_soil_layer, 1);

    if altmax_lastyear_profile > 0
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

for iseason = 1:season_num
    input_tot_litter1(:, iseason) = input_vector_litter1(:, iseason).*dz(1:n_soil_layer);
    input_tot_litter2(:, iseason) = input_vector_litter2(:, iseason).*dz(1:n_soil_layer);
    input_tot_litter3(:, iseason) = input_vector_litter3(:, iseason).*dz(1:n_soil_layer);
    input_tot_cwd(:, iseason) = input_vector_cwd(:, iseason).*dz(1:n_soil_layer);
end

input_tot_litter1 = sum(sum(input_tot_litter1));
input_tot_litter2 = sum(sum(input_tot_litter2));
input_tot_litter3 = sum(sum(input_tot_litter1));
input_tot_cwd = sum(sum(input_tot_cwd));

matrix_in(1:20,1) = input_tot_cwd*vertical_input./(dz(1:n_soil_layer)*days_per_year); % litter input gc/m3/day
matrix_in(21:40,1) = input_tot_litter1*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(41:60,1) = input_tot_litter2*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(61:80,1) = input_tot_litter3*vertical_input./(dz(1:n_soil_layer)*days_per_year);
matrix_in(81:140,1) = 0;

carbon_input = input_tot_cwd + input_tot_litter1 + input_tot_litter2 + input_tot_litter3;

cpool = (a_ma*kk_ma-tri_ma)\(-matrix_in);
cpool(cpool < 0) = 0;
initial_cpool = cpool;

cpool = [cpool(1:20), cpool(21:40), cpool(41:60), cpool(61:80), cpool(81:100), cpool(101:120), cpool(121:140)];
cpool_layer = sum(cpool, 2);
total_c_stock_initial = sum(cpool_layer.*dz(1:n_soil_layer)); % unit gC/m2

%----------------------------------------
% Analytical Solution of 14C
%----------------------------------------
% ss_atmo_14c_fraction = mean(atmo_14c_nh(2:21, 2)/1000 + 1);
ss_atmo_modern_fraction = atmo_14c(1, 2)/1000 + 1;

cpool_14c = (a_ma*kk_ma-tri_ma-(decay_lamda_matrix/days_per_year))\(-matrix_in*ss_atmo_modern_fraction);
cpool_14c(cpool_14c < 0) = 0;

initial_delta_14c_pool = (cpool_14c./initial_cpool - 1)*1000;
initial_delta_14c_pool(initial_delta_14c_pool == Inf) = (0 - 1)*1000;

%% warmup forward simu before formal simu
ss_threshold = 0.5; % unit 1gc/m2/year
total_c_stock_record = total_c_stock_initial;
for iyear = 1:warmup_year_size
    year_count = mod(iyear, ss_forcing_year_num);
    if year_count == 0
        year_count = ss_forcing_year_num;
    end

    if year_count == 1 && iyear ~= 1
        total_c_change_rate = abs((total_c_stock - total_c_stock_record)/ss_forcing_year_num);
        if total_c_change_rate < ss_threshold % unit 1gc/m2/year
            break
        else
            total_c_stock_record = total_c_stock;
        end
    end
    % env. scalars
    xit = nan(n_soil_layer, 1);
    xiw = nan(n_soil_layer, 1);

    for iseason = 1 : season_num
        if iseason == 1 && iyear == 1
            current_cpool = initial_cpool;
            current_delta_14c_pool = initial_delta_14c_pool;
        end
        %----------------------------------------
        % Environmental Scalar (Xi)
        %----------------------------------------
        % temperature related function xit
        % calculate rate constant scalar for soil temperature
        % assuming that the base rate constants are assigned for non-moisture
        % limiting conditions at 25 C.
        for ilayer = 1 : n_soil_layer
            if soil_temp_profile_transient(ilayer, iseason, year_count) >= 0 + kelvin_to_celsius
                xit(ilayer) = q10^((soil_temp_profile_transient(ilayer, iseason, year_count) - (kelvin_to_celsius + 25))/10);
            else
                xit(ilayer) = q10^((273.15 - 298.15)/10)*(fq10^((soil_temp_profile_transient(ilayer, iseason, year_count) - (0 + kelvin_to_celsius))/10));
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

        % for ilayer = 1 : n_soil_layer
        %    psi = min(maxpsi, soil_water_profile(ilayer, imonth));
        %
        %    if psi > minpsi
        %        xiw(ilayer, imonth) = (log(minpsi/psi))/(log(minpsi/maxpsi));
        %    else
        %        xiw(ilayer, imonth) = 0;
        %    end
        %    xiw(ilayer, imonth) = max(xiw(ilayer, imonth), 0.05);
        % end

        % xit = soil_temp_profile(:, imonth, iyear);

        if is_default == 0
            xiw = w_scaling*soil_water_profile_transient(:, iseason, year_count);
            xiw(xiw > 1) = 1;
        else
            xiw = soil_water_profile_transient(:, iseason, year_count);
        end

        %----------------------------------------
        % Triangle Matrix, A Matrix and K Matrix
        %----------------------------------------
        % Allocation matrix
        a_ma = a_matrix(is_default, fl1s1, fl2s1, fl3s2, fs1s2, fs1s3, fs2s1, fs2s3, fs3s1, fcwdl2, sand_vector_transient(:, iseason, year_count));

        kk_ma = kk_matrix(xit, xiw, xio_transient(:, iseason, year_count), xin_transient(:, iseason, year_count), efolding, tau4cwd, tau4l1, tau4l2, tau4l3, tau4s1, tau4s2, tau4s3);
        % tri matrix
        tri_ma = tri_matrix(nbedrock_transient(iseason, year_count), altmax_current_profile_transient(iseason, year_count), altmax_lastyear_profile_transient(iseason, year_count), bio, adv, cryo);

        %----------------------------------------
        % Vertical Profile
        %----------------------------------------
        % in the original beta model in Jackson et al 1996, the unit for the depth
        % of the soil is cm (dmax*100)
        if is_default == 0
            if isnan(grid_npp_allocation(1)) == 1
                m_to_cm = 100;
                vertical_prof = nan(n_soil_layer, 1);

                if altmax_lastyear_profile > 0
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
        end
        %----------------------------------------
        % input vector
        %----------------------------------------

        matrix_in = nan(140, 1);
        input_matrix = nan(140, 1);

        % monthly input
        input_vector_cwd_month = input_vector_cwd_transient(:, iseason, year_count);
        input_vector_litter1_month = input_vector_litter1_transient(:, iseason, year_count);
        input_vector_litter2_month = input_vector_litter2_transient(:, iseason, year_count);
        input_vector_litter3_month = input_vector_litter3_transient(:, iseason, year_count);
        % total input amount
        input_tot_cwd = sum(input_vector_cwd_month.*dz(1:n_soil_layer)); % (gc/m2/season)
        input_tot_litter1 = sum(input_vector_litter1_month.*dz(1:n_soil_layer));
        input_tot_litter2 = sum(input_vector_litter2_month.*dz(1:n_soil_layer));
        input_tot_litter3 = sum(input_vector_litter3_month.*dz(1:n_soil_layer));


        if is_default == 0
            % redistribution by beta
            matrix_in(1:20,1) = input_tot_cwd*vertical_input./(dz(1:n_soil_layer)*days_per_season(iseason)); % litter input gc/m3/day
            matrix_in(21:40,1) = input_tot_litter1*vertical_input./(dz(1:n_soil_layer)*days_per_season(iseason));
            matrix_in(41:60,1) = input_tot_litter2*vertical_input./(dz(1:n_soil_layer)*days_per_season(iseason));
            matrix_in(61:80,1) = input_tot_litter3*vertical_input./(dz(1:n_soil_layer)*days_per_season(iseason));
            matrix_in(81:140,1) = 0;
        else
            matrix_in(1:20,1) = input_vector_cwd_month./days_per_season(iseason); % litter input gc/m3/day
            matrix_in(21:40,1) = input_vector_litter1_month./days_per_season(iseason);
            matrix_in(41:60,1) = input_vector_litter2_month./days_per_season(iseason);
            matrix_in(61:80,1) = input_vector_litter3_month./days_per_season(iseason);
            matrix_in(81:140,1) = 0;
        end


        %----------------------------------------
        % soc stepwise iteration
        %----------------------------------------
        next_cpool = (matrix_in + (a_ma*kk_ma-tri_ma)*current_cpool)*days_per_season(iseason) + current_cpool;
        next_cpool(next_cpool < 0) = 0;
        %----------------------------------------
        % 14c stepwise iteration
        %----------------------------------------
        % convert soil delta 14c to fraction
        current_fraction_modern = current_delta_14c_pool/1000 + 1;
        % convert atmospheric delta 14c to fraction
        fraction_modern_atmo = atmo_14c(1, 2)/1000 + 1;
        % calculate fraction for the next time step, unit: g 14c/m3
        next_14c_pool = (matrix_in.*fraction_modern_atmo + (a_ma*kk_ma-tri_ma)*(current_fraction_modern.*current_cpool))*days_per_season(iseason)...
            - decay_lamda_matrix/days_per_year*days_per_season(iseason)*(current_fraction_modern.*current_cpool)...
            + (current_fraction_modern.*current_cpool);
        next_14c_pool(next_14c_pool < 0) = 0;

        next_delta_14c_pool = (next_14c_pool./next_cpool - 1)*1000;

        %----------------------------------------
        % carbon pool update
        %----------------------------------------
        current_cpool = next_cpool;
        current_delta_14c_pool = next_delta_14c_pool;

        cpool_layer = sum([current_cpool(1:20), current_cpool(21:40), current_cpool(41:60), current_cpool(61:80), ...
            current_cpool(81:100), current_cpool(101:120), current_cpool(121:140)], 2);
        total_c_stock = sum(cpool_layer.*dz(1:n_soil_layer)); % unit gC/m2
    end % imonth = ...
end

initial_cpool = current_cpool;
initial_delta_14c_pool = current_delta_14c_pool;

%% transient simulation
gradient_series = (0.0:0.1:1)';

% 12c related var
cpool_12c_stock_transient = nan(year_size, npool*depth_num, length(gradient_series), 2);
cpool_12c_resp_transient = nan(year_size, npool*depth_num, length(gradient_series), 2);
% 14c related var
cpool_14c_stock_transient = nan(year_size, npool*depth_num, length(gradient_series), 2);
coop_14c_resp_transient = nan(year_size, npool*depth_num, length(gradient_series), 2);

cpool_transient = nan(year_size, npool*n_soil_layer, length(gradient_series), 2);

anpp_origin = grid_npp_allocation(1);

for icomponent = 1:2
    for igradient = 1:length(gradient_series)
        if icomponent == 1 % change the input allocation
            anpp_remain = anpp_origin*(1-gradient_series(igradient));
            anpp_reallo = anpp_origin*gradient_series(igradient);
            bnpp_allo_new = anpp_reallo*grid_npp_allocation(2:end)/sum(grid_npp_allocation(2:end)) + grid_npp_allocation(2:end);
            grid_npp_allocation_new = [anpp_remain; bnpp_allo_new];

            interp_allo = interp1(soil_allocation_depth/100, grid_npp_allocation_new, zsoi(1:n_soil_layer, 1), 'linear');

            interp_allo(interp_allo < 0) = 0;
            interp_allo(isnan(interp_allo)) = 0;
            vertical_input_product = interp_allo/sum(interp_allo);
            % not change the v transport
            bio_new = bio;
            cryo_new = cryo;
        else % change the v transport
            % but restore the initial input allocation
            interp_allo = interp1(soil_allocation_depth/100, grid_npp_allocation, zsoi(1:n_soil_layer, 1), 'linear');
            interp_allo(interp_allo < 0) = 0;
            interp_allo(isnan(interp_allo)) = 0;
            vertical_input_product = interp_allo/sum(interp_allo);
            % change the v transport
            bio_new = bio*(1+gradient_series(igradient));
            cryo_new = cryo*(1+gradient_series(igradient));
        end

        if total_c_change_rate > ss_threshold*10 % which means here the threshold is 1gc/m2/yr
            % disp('no steady state reached before hist simu')
        else
            % env. scalars
            xit = nan(n_soil_layer, 1);
            xiw = nan(n_soil_layer, 1);

            for iyear = 1 : year_size
                matrix_in_annual = nan(140, season_num);
                input_matrix_annual = nan(140, season_num);
                % allocation matrix
                a_ma_annual = nan(npool*n_soil_layer, npool*n_soil_layer, season_num);
                kk_ma_annual = nan(npool*n_soil_layer, npool*n_soil_layer, season_num);
                % tri matrix
                tri_ma_annual = nan(npool*n_soil_layer, npool*n_soil_layer, season_num);

                env_scalar_annual = nan(season_num, n_soil_layer);
                carbon_input_annual =  nan(season_num, 1);
                soil_temp_annual =  nan(season_num, 1);
                soil_water_annual =  nan(season_num, 1);

                cpool_12c_stock_annual = nan(season_num, npool*4);
                cpool_12c_resp_annual = nan(season_num, npool*4);
                cpool_14c_stock_annual = nan(season_num, npool*4);
                cpool_14c_resp_annual = nan(season_num, npool*4);

                cpool_annual = nan(season_num, npool*n_soil_layer);
                for iseason = 1 : season_num
                    if iseason == 1 && iyear == 1
                        current_cpool = initial_cpool;
                        current_delta_14c_pool = initial_delta_14c_pool;
                    end
                    %----------------------------------------
                    % Environmental Scalar (Xi)
                    %----------------------------------------
                    % temperature related function xit
                    % calculate rate constant scalar for soil temperature
                    % assuming that the base rate constants are assigned for non-moisture
                    % limiting conditions at 25 C.
                    for ilayer = 1 : n_soil_layer
                        if soil_temp_profile_transient(ilayer, iseason, iyear) >= 0 + kelvin_to_celsius
                            xit(ilayer) = q10^((soil_temp_profile_transient(ilayer, iseason, iyear) - (kelvin_to_celsius + 25))/10);
                        else
                            xit(ilayer) = q10^((273.15 - 298.15)/10)*(fq10^((soil_temp_profile_transient(ilayer, iseason, iyear) - (0 + kelvin_to_celsius))/10));
                        end
                    end

                    catanf_30 = catanf(30);
                    normalization_tref = 15;
                    if normalize_q10_to_century_tfunc == true
                        % scale all decomposition rates by a constant to compensate for offset between original CENTURY temp func and Q10
                        normalization_factor = (catanf(normalization_tref)/catanf_30) / (q10^((normalization_tref-25)/10));
                        xit = xit*normalization_factor;
                    end

                    soil_temp_annual(iseason) = sum(soil_temp_profile_transient(1:n_soil_layer, iseason, iyear).*dz(1:n_soil_layer))/sum(dz(1:n_soil_layer));

                    soil_water_annual_middle = minpsi./(exp(soil_water_profile_transient(1:n_soil_layer, iseason, iyear)*(log(minpsi/maxpsi))));
                    soil_water_annual(iseason) = sum(soil_water_annual_middle.*dz(1:n_soil_layer))/sum(dz(1:n_soil_layer));

                    % water related function xiw
                    % calculate the rate constant scalar for soil water content.
                    % Uses the log relationship with water potential given in
                    % Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
                    % a comparison of models. Ecology, 68(5):1190-1200.
                    % and supported by data in
                    % Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
                    % and soil moisture. Soil Biol. Biochem., 15(4):447-453.

                    % for ilayer = 1 : n_soil_layer
                    %    psi = min(maxpsi, soil_water_profile(ilayer, imonth));
                    %
                    %    if psi > minpsi
                    %        xiw(ilayer, imonth) = (log(minpsi/psi))/(log(minpsi/maxpsi));
                    %    else
                    %        xiw(ilayer, imonth) = 0;
                    %    end
                    %    xiw(ilayer, imonth) = max(xiw(ilayer, imonth), 0.05);
                    % end

                    % xit = soil_temp_profile(:, imonth, iyear);

                    if is_default == 0
                        xiw = w_scaling*soil_water_profile_transient(:, iseason, iyear);
                        xiw(xiw > 1) = 1;
                    else
                        xiw = soil_water_profile_transient(:, iseason, iyear);
                    end

                    env_scalar_annual(iseason, :) = xit.*xiw.*exp(-zsoi(1:n_soil_layer)/efolding);

                    %----------------------------------------
                    % Triangle Matrix, A Matrix and K Matrix
                    %----------------------------------------
                    % Allocation matrix
                    a_ma = a_matrix(is_default, fl1s1, fl2s1, fl3s2, fs1s2, fs1s3, fs2s1, fs2s3, fs3s1, fcwdl2, sand_vector_transient(:, iseason, iyear));

                    kk_ma = kk_matrix(xit, xiw, xio_transient(:, iseason, iyear), xin_transient(:, iseason, iyear), efolding, tau4cwd, tau4l1, tau4l2, tau4l3, tau4s1, tau4s2, tau4s3);
                    % tri matrix
                    tri_ma = tri_matrix(nbedrock_transient(iseason, iyear), altmax_current_profile_transient(iseason, iyear), altmax_lastyear_profile_transient(iseason, iyear), bio_new, adv, cryo_new);

                    a_ma_annual(:, :, iseason) = a_ma;
                    kk_ma_annual(:, :, iseason) = kk_ma;
                    tri_ma_annual(:, :, iseason) = tri_ma;

                    %----------------------------------------
                    % Vertical Profile
                    %----------------------------------------
                    % in the original beta model in Jackson et al 1996, the unit for the depth
                    % of the soil is cm (dmax*100)
                    if is_default == 0
                        if isnan(vertical_input_product(1)) == 1
                            m_to_cm = 100;
                            vertical_prof = nan(n_soil_layer, 1);

                            if altmax_lastyear_profile > 0
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
                    end
                    %----------------------------------------
                    % input vector
                    %----------------------------------------

                    matrix_in = nan(140, 1);
                    input_matrix = nan(140, 1);

                    % monthly input
                    input_vector_cwd_month = input_vector_cwd_transient(:, iseason, iyear);
                    input_vector_litter1_month = input_vector_litter1_transient(:, iseason, iyear);
                    input_vector_litter2_month = input_vector_litter2_transient(:, iseason, iyear);
                    input_vector_litter3_month = input_vector_litter3_transient(:, iseason, iyear);
                    % total input amount
                    input_tot_cwd = sum(input_vector_cwd_month.*dz(1:n_soil_layer)); % (gc/m2/season)
                    input_tot_litter1 = sum(input_vector_litter1_month.*dz(1:n_soil_layer));
                    input_tot_litter2 = sum(input_vector_litter2_month.*dz(1:n_soil_layer));
                    input_tot_litter3 = sum(input_vector_litter3_month.*dz(1:n_soil_layer));
                    carbon_input = input_tot_cwd + input_tot_litter1 + input_tot_litter2 + input_tot_litter3;

                    carbon_input_annual(iseason) = carbon_input;

                    if is_default == 0
                        % redistribution by beta
                        matrix_in(1:20,1) = input_tot_cwd*vertical_input./(dz(1:n_soil_layer)*days_per_season(iseason)); % litter input gc/m3/day
                        matrix_in(21:40,1) = input_tot_litter1*vertical_input./(dz(1:n_soil_layer)*days_per_season(iseason));
                        matrix_in(41:60,1) = input_tot_litter2*vertical_input./(dz(1:n_soil_layer)*days_per_season(iseason));
                        matrix_in(61:80,1) = input_tot_litter3*vertical_input./(dz(1:n_soil_layer)*days_per_season(iseason));
                        matrix_in(81:140,1) = 0;

                        matrix_in_annual(:, iseason) = matrix_in;
                    else
                        matrix_in(1:20,1) = input_vector_cwd_month./days_per_season(iseason); % litter input gc/m3/day
                        matrix_in(21:40,1) = input_vector_litter1_month./days_per_season(iseason);
                        matrix_in(41:60,1) = input_vector_litter2_month./days_per_season(iseason);
                        matrix_in(61:80,1) = input_vector_litter3_month./days_per_season(iseason);
                        matrix_in(81:140,1) = 0;
                    end

                    if is_default == 0
                        input_matrix(1:20,1) = input_tot_cwd/carbon_input.*vertical_input./dz(1:n_soil_layer);
                        input_matrix(21:40,1) = input_tot_litter1/carbon_input.*vertical_input./dz(1:n_soil_layer);
                        input_matrix(41:60,1) = input_tot_litter2/carbon_input.*vertical_input./dz(1:n_soil_layer);
                        input_matrix(61:80,1) = input_tot_litter3/carbon_input.*vertical_input./dz(1:n_soil_layer);
                        input_matrix(81:140,1) = 0;

                        input_matrix_annual(:, iseason) = input_matrix;
                    else
                        vertical_input = (input_vector_cwd + input_vector_litter1 + input_vector_litter2 + input_vector_litter3).*dz(1:n_soil_layer)./(input_tot_cwd + input_tot_litter1 + input_tot_litter2 + input_tot_litter3);
                        input_matrix(1:20,1) = input_tot_cwd/carbon_input.*input_vector_cwd.*dz(1:n_soil_layer)./input_tot_cwd./dz(1:n_soil_layer);
                        input_matrix(21:40,1) = input_tot_litter1/carbon_input.*input_vector_litter1.*dz(1:n_soil_layer)./input_tot_litter1./dz(1:n_soil_layer);
                        input_matrix(41:60,1) = input_tot_litter2/carbon_input.*input_vector_litter2.*dz(1:n_soil_layer)./input_tot_litter2./dz(1:n_soil_layer);
                        input_matrix(61:80,1) = input_tot_litter3/carbon_input.*input_vector_litter3.*dz(1:n_soil_layer)./input_tot_litter3./dz(1:n_soil_layer);
                        input_matrix(81:140,1) = 0;
                    end

                    %----------------------------------------
                    % soc stepwise iteration
                    %----------------------------------------
                    next_cpool = (matrix_in + (a_ma*kk_ma-tri_ma)*current_cpool)*days_per_season(iseason) + current_cpool; % unit: g c/m3
                    next_cpool(next_cpool < 0) = 0;
                    %----------------------------------------
                    % 14c stepwise iteration
                    %----------------------------------------
                    % convert soil delta 14c to fraction
                    current_fraction_modern = current_delta_14c_pool/1000 + 1;
                    % convert atmospheric delta 14c to fraction
                    fraction_modern_atmo = atmo_14c(iyear, 2)/1000 + 1;
                    % calculate fraction for the next time step, unit: g 14c/m3
                    next_14c_pool = (matrix_in.*fraction_modern_atmo + (a_ma*kk_ma-tri_ma)*(current_fraction_modern.*current_cpool))*days_per_season(iseason)...
                        - decay_lamda_matrix/days_per_year*days_per_season(iseason)*(current_fraction_modern.*current_cpool)...
                        + (current_fraction_modern.*current_cpool);
                    next_14c_pool(next_14c_pool < 0) = 0;


                    next_delta_14c_pool = (next_14c_pool./next_cpool - 1)*1000;


                    %----------------------------------------
                    % respired flux 14c calculation
                    %----------------------------------------
                    r_vector = -sum(a_ma*kk_ma, 1); % respiration rate, unit: day-1

                    respired_cpool = r_vector'.*(current_cpool)*days_per_season(iseason); % unit: gc m-3 month-1

                    % 14 c
                    respired_14c_pool = r_vector'.*(current_fraction_modern.*current_cpool)*days_per_season(iseason); % unit: g 14c m-3 month-1

                    %----------------------------------------
                    % annual record
                    %----------------------------------------
                    cpool_annual(iseason, :) = next_cpool;
                    for ipool = 1:npool
                        cpool_12c_stock_annual(iseason, ((ipool-1)*depth_num+1):ipool*depth_num) = fun_layer_summary(next_cpool(((ipool-1)*n_soil_layer+1):ipool*n_soil_layer));
                        cpool_14c_stock_annual(iseason, ((ipool-1)*depth_num+1):ipool*depth_num) = fun_layer_summary(next_14c_pool(((ipool-1)*n_soil_layer+1):ipool*n_soil_layer));
                        cpool_12c_resp_annual(iseason, ((ipool-1)*depth_num+1):ipool*depth_num) = fun_layer_summary(respired_cpool(((ipool-1)*n_soil_layer+1):ipool*n_soil_layer));
                        cpool_14c_resp_annual(iseason, ((ipool-1)*depth_num+1):ipool*depth_num) = fun_layer_summary(respired_14c_pool(((ipool-1)*n_soil_layer+1):ipool*n_soil_layer));

                    end

                    %----------------------------------------
                    % transient record
                    %----------------------------------------
                    if iseason == season_num
                        % annual mean
                        cpool_annual_mean = mean(cpool_annual, 1, 'omitnan')';

                        cpool_transient(iyear, :, igradient, icomponent) = cpool_annual_mean;

                        cpool_12c_stock_transient(iyear, :, igradient, icomponent) = mean(cpool_12c_stock_annual, 1, 'omitnan')';
                        cpool_14c_stock_transient(iyear, :, igradient, icomponent) = mean(cpool_14c_stock_annual, 1, 'omitnan')';

                        cpool_12c_resp_transient(iyear, :, igradient, icomponent) = sum(cpool_12c_resp_annual, 1, 'omitnan')';
                        coop_14c_resp_transient(iyear, :, igradient, icomponent) = sum(cpool_14c_resp_annual, 1, 'omitnan')';

                    end

                    %----------------------------------------
                    % carbon pool update
                    %----------------------------------------
                    current_cpool = next_cpool;

                    current_delta_14c_pool = next_delta_14c_pool;
                end % imonth = ...
            end % iyear = ...

            %     plot(1900:(2014-1900)/(length(flux_14c_transient)-1):2014, flux_14c_transient)
            %     hold on
            %     plot(atmo_14c_nh(:, 1), atmo_14c_nh(:, 2))
        end
    end
end
% 
% for igradient = 1:11
%     plot(sum(cpool_12c_stock_transient(:, 4:depth_num:end, igradient, 1)/1000, 2))
%     hold on
% end
% ylim([78, 82])
% 
% for igradient = 1:10
%     plot(sum(cpool_12c_stock_transient(:, 4:depth_num:end, igradient, 2)/1000, 2))
%     hold on
% end
% ylim([78, 82])

end  % function

% ----- CENTURY T response function
function catanf_results = catanf(t1)
catanf_results = 11.75 +(29.7 / pi) * atan( pi * 0.031  * ( t1 - 15.4 ));
end

