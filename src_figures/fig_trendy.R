## Packages
library(R.matlab)
library(ggplot2)
library(cowplot)
library(viridis)
# library(matlab)
library(GGally)
library(raster)

library(sf)
library(sp)
library(proj4)


# dev.off()
##
rm(list = ls())

setwd('/Users/ft254/Google_Drive/R')

#############################################################################
# Data Path
#############################################################################
model_name = 'cesm2_clm5_cen_vr_v2'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

data_dir_output = '/Users/ft254/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/ft254/DATAHUB/ENSEMBLE/INPUT_DATA/'

#############################################################################
# function to increase vertical spacing between legend keys
#############################################################################
# @clauswilke
draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.6, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}

# register new key drawing function, 
# the effect is global & persistent throughout the R session
GeomBar$draw_key = draw_key_polygon3

#################################################################################
# CESM2 valid grid in simulation
#################################################################################
valid_grid_loc = read.csv(paste(data_dir_output, 'nn_soil_carbon_sink/global_grid_valid_id.csv', sep = ''), header = FALSE)
valid_grid_loc = valid_grid_loc$V1

grid_var_names = c('Lon', 'Lat', 'Date', 
                   'Rmean', 'Rmax', 'Rmin', 
                   'ESA_Land_Cover', 
                   'ET',
                   'IGBP', 'Climate', 'Soil_Type', 'NPPmean', 'NPPmax', 'NPPmin',
                   'Veg_Cover', 
                   'Annual Mean Temperature', 'Mean Diurnal Range', 'Isothermality', 'Temperature Seasonality', 'Max Temperature of Warmest Month', 'Min Temperature of Coldest Month', 'Temperature Annual Range', 'Mean Temperature of Wettest Quarter', 'Mean Temperature of Driest Quarter', 'Mean Temperature of Warmest Quarter', 'Mean Temperature of Coldest Quarter', 'Annual Precipitation', 'Precipitation of Wettest Month', 'Precipitation of Driest Month', 'Precipitation Seasonality', 'Precipitation of Wettest Quarter', 'Precipitation of Driest Quarter', 'Precipitation of Warmest Quarter', 'Precipitation of Coldest Quarter', 
                   'Abs_Depth_to_Bedrock',
                   'Bulk_Density_0cm', 'Bulk_Density_30cm', 'Bulk_Density_100cm',
                   'CEC_0cm', 'CEC_30cm', 'CEC_100cm',
                   'Clay_Content_0cm', 'Clay_Content_30cm', 'Clay_Content_100cm',
                   'Coarse_Fragments_v_0cm', 'Coarse_Fragments_v_30cm', 'Coarse_Fragments_v_100cm', 
                   'Depth_Bedrock_R', 
                   'Garde_Acid', 
                   'Occurrence_R_Horizon', 
                   'pH_Water_0cm', 'pH_Water_30cm', 'pH_Water_100cm', 
                   'Sand_Content_0cm', 'Sand_Content_30cm', 'Sand_Content_100cm',
                   'Silt_Content_0cm', 'Silt_Content_30cm', 'Silt_Content_100cm', 
                   'SWC_v_Wilting_Point_0cm', 'SWC_v_Wilting_Point_30cm', 'SWC_v_Wilting_Point_100cm', 
                   'Texture_USDA_0cm', 'Texture_USDA_30cm', 'Texture_USDA_100cm', 
                   'USDA_Suborder', 
                   'WRB_Subgroup', 
                   'Drought',
                   'Elevation',
                   'Max_Depth', 
                   'Koppen_Climate_2018', 
                   'cesm2_npp', 'cesm2_npp_std',
                   'cesm2_gpp', 'cesm2_gpp_std',
                   'cesm2_vegc',
                   'nbedrock')

grid_env_info = readMat(paste(data_dir_input, 'data4nn/world_grid_envinfo_present_2024.mat', sep = ''))
grid_env_info = grid_env_info$EnvInfo[valid_grid_loc, ]
colnames(grid_env_info) = grid_var_names

global_lat_lon = grid_env_info[ , c('Lon', 'Lat')]

#################################################################################
# trendy original simulation
#################################################################################
# model list in trendy
trendy_model_list = c('CABLE-POP', 'CLASSIC', 'CLM5.0-CRU', 'DLEM', 'IBIS', 'ISAM',
                      'ISBA-CTRIP', 'JSBACH', 'JULES-ES-1p0', 'LPJ-GUESS', 'LPX-Bern', 'OCN',
                      'ORCHIDEE', 'ORCHIDEE-CNP', 'ORCHIDEEv3', 'SDGVM', 'VISIT', 'YIBs', 'CLM5.0-GSWP3v1')

lon_grid_mask = seq(from = -180 + 0.5/2, to = 180 - 0.5/2, by = 0.5)
lat_grid_mask = seq(from = 90 - 0.5/2, to = -90 + 0.5/2, by = -0.5)

grid_lat_map = matlab::repmat(lat_grid_mask, c(1, 720))
grid_lon_map = matlab::repmat(lon_grid_mask, c(360, 1))


land_mask = array(NA, dim = c(360, 720))
for (igrid in 1:length(valid_grid_loc)) {
  target_lon_coord = length(seq(from = -180 + 0.5/2, to = global_lat_lon[igrid, 'Lon'], by = 0.5))
  target_lat_coord = length(seq(from = 90 - 0.5/2, to = global_lat_lon[igrid, 'Lat'], by = -0.5))
  land_mask[target_lat_coord, target_lon_coord] = 1
}


trendy_lon_lat = cbind(grid_lon_map[which(land_mask == 1)], 
                       grid_lat_map[which(land_mask == 1)])

# trendy soc hist
trendy_soc_hist = readMat(paste( data_dir_input, 'trendy_hist/trendy_soc_tot_trend_s_2.mat', sep = ''))
trendy_soc_hist = trendy_soc_hist$summary.soc.time.series
trendy_soc_hist = cbind(trendy_soc_hist, NA)
# trendy npp hist
trendy_npp_hist = readMat(paste( data_dir_input, 'trendy_hist/trendy_npp_tot_trend_s_2.mat', sep = ''))
trendy_npp_hist = trendy_npp_hist$summary.npp.time.series
trendy_npp_hist = cbind(trendy_npp_hist, NA)

# trendy veg hist
trendy_veg_hist = readMat(paste( data_dir_input, 'trendy_hist/trendy_veg_tot_trend_s_2.mat', sep = ''))
trendy_veg_hist = trendy_veg_hist$summary.veg.time.series
trendy_veg_hist = cbind(trendy_veg_hist, NA)

trendy_sink_tot = array(NA, dim = c(nrow(global_lat_lon), 19))
trendy_sink_rate = array(NA, dim = c(nrow(global_lat_lon), 19))

# imodel = 1
# for (imodel in 1:length(trendy_model_list)) {
#   print(paste('process model ', imodel, sep = ''))
#   data_path = paste(data_dir_input, 'trendy_hist/trendy_soc_1900_2010_transient_s_2_',  trendy_model_list[imodel], '.mat', sep = '')
#   if (file.exists(data_path) == TRUE) {
#     trendy_simu_middle = readMat(data_path)
#     trendy_simu_middle = trendy_simu_middle$summary.soc.transient
# 
#     trendy_simu_tot_middle = (trendy_simu_middle[ , , 111] - trendy_simu_middle[ , , 1])*1000
#     trendy_simu_rate_middle = (trendy_simu_middle[ , , 111] - trendy_simu_middle[ , , 101])*1000
# 
#     trendy_sink_tot[ , imodel] =  trendy_simu_tot_middle[which(land_mask == 1)]
#     trendy_sink_rate[ , imodel] =  trendy_simu_rate_middle[which(land_mask == 1)]
#   }
# }
# 
# writeMat(paste( data_dir_output, 'hist_simu_proda/trendy_input/model_ensemble_trendy_sink_tot.mat', sep = ''), data = trendy_sink_tot)
# writeMat(paste( data_dir_output, 'hist_simu_proda/trendy_input/model_ensemble_trendy_sink_rate.mat', sep = ''), data = trendy_sink_rate)


trendy_sink_tot = readMat(paste( data_dir_output, 'hist_simu_proda/trendy_input/model_ensemble_trendy_sink_tot.mat', sep = ''), data = trendy_sink_tot)
trendy_sink_tot = trendy_sink_tot$data

trendy_sink_rate = readMat(paste( data_dir_output, 'hist_simu_proda/trendy_input/model_ensemble_trendy_sink_rate.mat', sep = ''), data = trendy_sink_rate)
trendy_sink_rate = trendy_sink_rate$data

#################################################################################
# cesm2 hist simulation with trendy input
#################################################################################
year_start = 1900
year_end = 2010

year_series = seq(from = 1, to = 111, by = 10)
year_num = length(year_series)
iscenario = 3
idepth = 4

# trendy npp
treny_npp_mean = readMat(paste( data_dir_input, 'trendy_hist/global_grid_trendy_npp_mean.mat', sep = ''))
treny_npp_mean = cbind(treny_npp_mean$grid.trendy.npp.mean[valid_grid_loc, ], grid_env_info[ , 'cesm2_npp'])

## grid_area    
radius = 6371008.8
resolution_lon = 0.5
resolution_lat = 0.5

length_top = (2*pi*radius*cos(abs(global_lat_lon[ , 2]+resolution_lat/2)/180*pi)/360)*resolution_lon
length_down = (2*pi*radius*cos(abs(global_lat_lon[ , 2]-resolution_lat/2)/180*pi)/360)*resolution_lon
height = (pi*radius/180)*resolution_lat
grid_area = (length_top + length_down)*height/2

##------------------------- clm5 simu
# clm5_carbon_input = array(NA, dim = c(45230, 12, 10))
# clm5_soc_stock = array(NA, dim = c(45230, 12, 10))
#
# cross_valid_num = 8
# for (cross_valid_num in 1:10) {
#   print(paste('processing cross valid ', as.character(cross_valid_num), sep = ''))
#   # clm5 input
#   process_hist_raw = readMat(paste( data_dir_output, 'hist_simu_proda/scenario_simu/cesm2_hist_process_proda_default_cross_valid_', as.character(cross_valid_num), '_scenario_simu.mat', sep = ''))
#   clm5_carbon_input_middle = process_hist_raw$process.map.cesm2[ , , 21, 1]
#
#   clm5_carbon_input[ , , cross_valid_num] = clm5_carbon_input_middle
#
#   # clm5 soc
#   soc_hist_raw = readMat(paste(data_dir_output, 'hist_simu_proda/scenario_simu/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_scenario_simu.mat', sep = ''))
#   clm5_soc_stock_middle = apply(soc_hist_raw$soc.stock.map.cesm2[ , , 17:28, 4], c(1, 2), sum, na.rm = TRUE)
#
#   clm5_soc_stock[ , , cross_valid_num] = clm5_soc_stock_middle
# }
#
# clm5_carbon_input = apply(clm5_carbon_input, c(1, 2), mean, na.rm = TRUE)
# clm5_soc_stock = apply(clm5_soc_stock, c(1, 2), mean, na.rm = TRUE)
# 
##------------------------- sink for total soil carbon pool
# mip_sink_tot = array(NA, dim = c(length(grid_area), 19, 10))
# mip_sink_rate = array(NA, dim = c(length(grid_area), 19, 10))
# mip_input_change = array(NA, dim = c(length(grid_area), 19, 10))
# mip_input_mean = array(NA, dim = c(length(grid_area), 19, 10))
#
#
# for (cross_valid_num in 1:10) {
#   imodel = 1
#   for (imodel in 1:19) {
#     print(paste('cross valid ', as.character(cross_valid_num), ' processing model ', imodel, sep = ''))
#
#     if (imodel == 19) {
#       global_carbon_stock_tot = clm5_soc_stock
#       global_carbon_input = clm5_carbon_input
#     }else {
#
#       data_path = paste( data_dir_output, 'hist_simu_proda/trendy_input/cesm2_hist_transient_summary_trendy_model_',
#                          as.character(imodel), '_hist_f05_g16_checked_', as.character(year_start), '_', as.character(year_end),
#                          '_cross_valid_', as.character(cross_valid_num), '_scenario_', as.character(iscenario), '.mat', sep = '')
#
#       # soc stock simu
#       if (file.exists(data_path) == TRUE) {
#         his_summary = readMat(data_path)
#         global_carbon_stock = his_summary$hist.summary
#         # global_carbon_input = his_summary$hist.summary[[2]][ , , 6]
#         # total soc storage for the whole depth
#         global_carbon_stock_tot = apply(global_carbon_stock[ , , 17:28], c(1, 2), sum, na.rm = FALSE)
#
#       } else {
#         global_carbon_stock_tot = array(NA, dim = c(length(valid_grid_loc), year_num))
#         # global_carbon_input = array(NA, dim = c(length(valid_grid_loc), year_num))
#       }
#     }
#
#     accu_change_carbon_stock_tot = global_carbon_stock_tot - matlab::repmat(global_carbon_stock_tot[ , 1], c(1, year_num))
#     invalid_loc = which(global_carbon_stock_tot[ , 1] > 1000000 | global_carbon_stock_tot[ , 1] < 0 |
#                           accu_change_carbon_stock_tot[ , year_num] < -10000 |
#                           accu_change_carbon_stock_tot[ , year_num] > 10000)
#
#     accu_change_carbon_stock_tot[invalid_loc, ] = NA
#     global_carbon_stock_tot[invalid_loc, ] = NA
#
#     mip_sink_tot[ , imodel, cross_valid_num] = accu_change_carbon_stock_tot[ , year_num]
#     mip_sink_rate[ , imodel, cross_valid_num] = accu_change_carbon_stock_tot[ , year_num] - accu_change_carbon_stock_tot[ , (year_num - 1)]
#     mip_input_mean[ , imodel, cross_valid_num] = treny_npp_mean[ , imodel]
#     # mip_input_change[ , imodel] = global_carbon_input[ , year_num] - global_carbon_input[ , 1]
#
#   }
# }
#
# mip_sink_tot = apply(mip_sink_tot, c(1, 2), mean, na.rm = TRUE)
# mip_sink_rate = apply(mip_sink_rate, c(1, 2), mean, na.rm = TRUE)
# mip_input_mean = apply(mip_input_mean, c(1, 2), mean, na.rm = TRUE)
# mip_input_change = apply(mip_input_change, c(1, 2), mean, na.rm = TRUE)
#
# writeMat(paste( data_dir_output, 'hist_simu_proda/trendy_input/model_ensemble_mip_sink_tot.mat', sep = ''), data = mip_sink_tot)
# writeMat(paste( data_dir_output, 'hist_simu_proda/trendy_input/model_ensemble_mip_sink_rate.mat', sep = ''), data = mip_sink_rate)
# writeMat(paste( data_dir_output, 'hist_simu_proda/trendy_input/model_ensemble_mip_input_mean.mat', sep = ''), data = mip_input_mean)
# writeMat(paste( data_dir_output, 'hist_simu_proda/trendy_input/model_ensemble_mip_input_change.mat', sep = ''), data = mip_input_change)

mip_sink_tot = readMat(paste( data_dir_output, 'hist_simu_proda/trendy_input/model_ensemble_mip_sink_tot.mat', sep = ''), data = mip_sink_tot)
mip_sink_tot = mip_sink_tot$data

mip_sink_rate = readMat(paste( data_dir_output, 'hist_simu_proda/trendy_input/model_ensemble_mip_sink_rate.mat', sep = ''), data = mip_sink_rate)
mip_sink_rate = mip_sink_rate$data

mip_input_mean = readMat(paste( data_dir_output, 'hist_simu_proda/trendy_input/model_ensemble_mip_input_mean.mat', sep = ''), data = mip_input_mean)
mip_input_mean = mip_input_mean$data

mip_input_change = readMat(paste( data_dir_input, '/trendy_hist/trendy_npp_change_map.mat', sep = ''))
mip_input_change = mip_input_change$trendy.npp.change
mip_input_change = mip_input_change[ , 2, ] - mip_input_change[ , 1, ]

#########################################################
# figures - boxes
#########################################################
color_scheme = c('#08519c', '#6baed6', '#c6dbef',
                 '#006d2c', '#74c476', '#c7e9c0', 
                 '#252525', '#969696', '#d9d9d9', 
                 '#d94801', '#fd8d3c', '#fdd0a2', 
                 '#54278f', '#9e9ac8', '#dadaeb', 
                 '#a50f15', '#fb6a4a', '#fcbba1',
                 '#D41159')


sink_tot_model_summary_proda = apply(mip_sink_tot*matlab::repmat(grid_area, c(1, 19)), 2, sum, na.rm = TRUE)/10**15
sink_tot_model_summary_proda[sink_tot_model_summary_proda == 0] = NA
sink_tot_model_summary_trendy = apply(trendy_sink_tot*matlab::repmat(grid_area, c(1, 19)), 2, sum, na.rm = TRUE)/10**15
sink_tot_model_summary_trendy[sink_tot_model_summary_trendy == 0] = NA

rate_model_summary_proda = apply(mip_sink_rate*matlab::repmat(grid_area, c(1, 19)), 2, sum, na.rm = TRUE)/10**15/10
rate_model_summary_proda[rate_model_summary_proda == 0] = NA
rate_model_summary_trendy = apply(trendy_sink_rate*matlab::repmat(grid_area, c(1, 19)), 2, sum, na.rm = TRUE)/10**15/10
rate_model_summary_trendy[rate_model_summary_trendy == 0] = NA

current_data = rbind(cbind(sink_tot_model_summary_trendy, rate_model_summary_trendy, 'source' = 1, c(1:19)),
                     cbind(sink_tot_model_summary_proda, rate_model_summary_proda, 'source' = 2, c(1:19)))

current_data = data.frame(current_data)

colnames(current_data) = c('accu', 'rate', 'source', 'model')
current_data[current_data == 0] = NA

var(current_data$accu[current_data$source == 1], na.rm = TRUE)
var(current_data$accu[current_data$source == 2], na.rm = TRUE)
var(current_data$accu[current_data$source == 2], na.rm = TRUE)/var(current_data$accu[current_data$source == 1], na.rm = TRUE)

var(current_data$rate[current_data$source == 1], na.rm = TRUE)
var(current_data$rate[current_data$source == 2], na.rm = TRUE)
var(current_data$rate[current_data$source == 2], na.rm = TRUE)/var(current_data$rate[current_data$source == 1], na.rm = TRUE)

box_hinge = rbind(quantile(current_data$accu[current_data$source == 1], probs = c(0.025, 0.16, 0.5, 0.84, 0.975), na.rm = TRUE) ,
                  quantile(current_data$accu[current_data$source == 2], probs = c(0.025, 0.16, 0.5, 0.84, 0.975), na.rm = TRUE) )


p_box_accu =
  ggplot() +
  geom_boxplot(data = current_data, aes(x = as.factor(source), y = accu), size = 2, outlier.colour = NA, ymin = box_hinge[ , 1], lower = box_hinge[ , 2], middle = box_hinge[ , 3], upper = box_hinge[ , 4], ymax = box_hinge[ , 5]) + 
  # geom_violin(data = current_data, aes(x = as.factor(source), y = accu), size = 2, outlier.colour = NA, trim = FALSE) +
  geom_jitter(data = current_data, aes(x = as.factor(source), y = accu, color = as.factor(model)), position = position_jitter(0.1), shape = 16, size = 9) +
  scale_x_discrete(labels = c('TRENDY', 'PRODA')) +
  scale_color_manual(name = '', values = color_scheme, label = trendy_model_list) +   
  theme_classic() +
  # add title
  labs(title = '', x = '', y = expression(paste('Total sink (Pg C)', sep = ''))) +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'vertical')  +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the margin
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 35), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))

p_box_accu_proda =
  ggplot() +
  geom_boxplot(data = current_data[current_data$source == 2, ], aes(x = as.factor(source), y = accu), size = 2, outlier.colour = NA, ymin = box_hinge[2, 1], lower = box_hinge[2, 2], middle = box_hinge[2, 3], upper = box_hinge[2, 4], ymax = box_hinge[2, 5]) + 
  # geom_violin(data = current_data, aes(x = as.factor(source), y = accu), size = 2, outlier.colour = NA, trim = FALSE) +
  geom_jitter(data = current_data[current_data$source == 2, ], aes(x = as.factor(source), y = accu, color = as.factor(model)), position = position_jitter(0.1), shape = 16, size = 9) +
  scale_x_discrete(labels = 'PRODA') +
    scale_y_continuous(limits = c(0, 50)) + 
  scale_color_manual(name = '', values = color_scheme, label = trendy_model_list) +   
  theme_classic() +
  # add title
  labs(title = '', x = '', y = expression(paste('Total sink (Pg C)', sep = ''))) +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'vertical')  +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the margin
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 35), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))


box_hinge = rbind(quantile(current_data$rate[current_data$source == 1], probs = c(0.025, 0.16, 0.5, 0.84, 0.975), na.rm = TRUE) ,
                  quantile(current_data$rate[current_data$source == 2], probs = c(0.025, 0.16, 0.5, 0.84, 0.975), na.rm = TRUE) )

p_box_rate =
  ggplot() +
  # geom_violin(data = current_data, aes(x = as.factor(source), y = rate), size = 2, trim = FALSE) +
  geom_boxplot(data = current_data, aes(x = as.factor(source), y = rate), size = 2, outlier.colour = NA, ymin = box_hinge[ , 1], lower = box_hinge[ , 2], middle = box_hinge[ , 3], upper = box_hinge[ , 4], ymax = box_hinge[ , 5]) + 
  geom_jitter(data = current_data, aes(x = as.factor(source), y = rate, color = as.factor(model)), position = position_jitter(0.1), shape = 16, size = 9) +
  scale_x_discrete(labels = c('TRENDY', 'PRODA')) +
  scale_color_manual(name = '', values = color_scheme, label = trendy_model_list) +   
  theme_classic() +
  # add title
  labs(title = '', x = '', y = expression(paste('Sink rate (Pg C yr'^'-1', ')', sep = ''))) +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'vertical')  +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the margin
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 35), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))


p_box_rate = p_box_rate + 
  theme(legend.justification = c(0, 1), legend.position = c(1.1, 1.1), legend.background = element_rect(fill = NA)) 

p_box_rate_proda =
  ggplot() +
  # geom_violin(data = current_data, aes(x = as.factor(source), y = rate), size = 2, trim = FALSE) +
  geom_boxplot(data = current_data[current_data$source == 2, ], aes(x = as.factor(source), y = rate), size = 2, outlier.colour = NA, ymin = box_hinge[2, 1], lower = box_hinge[2, 2], middle = box_hinge[2, 3], upper = box_hinge[2, 4], ymax = box_hinge[2, 5]) + 
  geom_jitter(data = current_data[current_data$source == 2, ], aes(x = as.factor(source), y = rate, color = as.factor(model)), position = position_jitter(0.1), shape = 16, size = 9) +
  scale_x_discrete(labels = c('PRODA')) +
  scale_color_manual(name = '', values = color_scheme, label = trendy_model_list) +   
    scale_y_continuous(limits = c(0, 1.0)) + 
  theme_classic() +
  # add title
  labs(title = '', x = '', y = expression(paste('Sink rate (Pg C yr'^'-1', ')', sep = ''))) +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'vertical')  +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the margin
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 35), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))

#########################################################
# figures - maps
#########################################################
color_scheme = c('#b2182b', '#ef8a62', '#fddbc7', '#f7f7f7', '#d1e5f0', '#67a9cf', '#2166ac')

world_coastline = st_read('/Users/ft254/Google_Drive/Tsinghua_Luo/World_Vector_Shape/ne110m/ne_110m_land.shp', layer = 'ne_110m_land')
world_coastline <- st_transform(world_coastline, CRS('+proj=robin'))

ocean_left = cbind(rep(-180, 100), seq(from = 80, to = -56, by = -(80 + 56)/(100 -1)))
ocean_right = cbind(rep(180, 100), seq(from = -56, to = 80, by = (80 + 56)/(100 -1)))
ocean_top = cbind(seq(from = 180, to = -180, by = -(360)/(100 -1)), rep(80, 100))
ocean_bottom = cbind(seq(from = -180, to = 180, by = (360)/(100 -1)), rep(-56, 100))

world_ocean = rbind(ocean_left, ocean_bottom, ocean_right, ocean_top)
world_ocean = as.matrix(world_ocean)

world_ocean <- project(xy = world_ocean, proj = '+proj=robin')

world_ocean = data.frame(world_ocean)
colnames(world_ocean) = c('long', 'lat')


#----------------------------------------------
# TRENDY soil carbon sink, total
#----------------------------------------------
current_data = data.frame(cbind(trendy_lon_lat, 
                                apply(trendy_sink_tot, 1, median, na.rm = TRUE),
                                apply(trendy_sink_tot, 1, sd, na.rm = TRUE),
                                apply(trendy_sink_tot, 1, quantile, prob = 0.25, na.rm = TRUE),
                                apply(trendy_sink_tot, 1, quantile, prob = 0.75, na.rm = TRUE), 
                                0))


colnames(current_data) = c('lon', 'lat', 'sink', 'std', 'lower_sigma', 'upper_sigma', 'sig_index')
current_data$sig_index[current_data$sink > 0 & current_data$lower_sigma > 0] = 1
current_data$sig_index[current_data$sink < 0 & current_data$upper_sigma < 0] = 1
# current_data$sig_index[(abs(current_data$sink) -  current_data$std) > 0] = 1


lon_lat_transfer = project(xy = as.matrix(current_data[ , c('lon', 'lat')]), proj = '+proj=robin') 
current_data[ , c('lon', 'lat')] = lon_lat_transfer

lat_limits = rbind(c(0, -56), c(0, 80))
lat_limits_robin = project(xy = as.matrix(lat_limits), proj = '+proj=robin') 

p_trendy_sink_total =
  ggplot() +
  geom_tile(data = current_data, aes(x = lon, y = lat, fill = sink), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradient2(expression(paste('Total sink \n(1900 - 2010) g C m'^'-2', sep = '')), low = '#d73027', mid = 'white', high = '#3288bd', na.value="transparent", limits = c(-500, 500), trans = 'identity', oob = scales::squish) +
  # geom_point(data = current_data[current_data$sig_index == 0, ], aes(x = lon, y = lat), shape = 4, size = 0.3, stroke = 0.1, color = 'black') +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
  # change the background to black and white
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  theme(legend.justification = c(0, 0), legend.position = c(0.1, -0.2), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 30, barheight = 2.5, title.position = 'right', title.hjust = 1, label.hjust = 0.5, label.vjust = 2.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 30, ), legend.title = element_text(size = 35)) +
  # add title
  labs(title = 'TRENDY ensemble', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))

p_trendy_sink_total_std =
  ggplot() +
  geom_tile(data = current_data, aes(x = lon, y = lat, fill = std), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradientn(expression(paste('Standard deviation \n(1900 - 2010) g C m'^'-2', sep = '')), colours = c('#fef0d9','#fdcc8a','#fc8d59','#e34a33','#b30000'), na.value="transparent", limits = c(0, 500), trans = 'identity', oob = scales::squish) +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
  # change the background to black and white
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  theme(legend.justification = c(0, 0), legend.position = c(0.1, -0.2), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 30, barheight = 2.5, title.position = 'right', title.hjust = 1, label.hjust = 0.5, label.vjust = 2.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 30, ), legend.title = element_text(size = 35)) +
  # add title
  labs(title = 'TRENDY ensemble', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))



#----------------------------------------------
# TRENDY soil carbon sink, rate
#----------------------------------------------

current_data = data.frame(cbind(trendy_lon_lat, 
                                apply(trendy_sink_rate/10, 1, median, na.rm = TRUE),
                                apply(trendy_sink_rate/10, 1, sd, na.rm = TRUE),
                                apply(trendy_sink_rate/10, 1, quantile, prob = 0.25, na.rm = TRUE),
                                apply(trendy_sink_rate/10, 1, quantile, prob = 0.75, na.rm = TRUE), 
                                0))


colnames(current_data) = c('lon', 'lat', 'sink', 'std', 'lower_sigma', 'upper_sigma', 'sig_index')
current_data$sig_index[current_data$sink > 0 & current_data$lower_sigma > 0] = 1
current_data$sig_index[current_data$sink < 0 & current_data$upper_sigma < 0] = 1


lon_lat_transfer = project(xy = as.matrix(current_data[ , c('lon', 'lat')]), proj = '+proj=robin') 
current_data[ , c('lon', 'lat')] = lon_lat_transfer

lat_limits = rbind(c(0, -56), c(0, 80))
lat_limits_robin = project(xy = as.matrix(lat_limits), proj = '+proj=robin') 

p_trendy_sink_rate =
  ggplot() +
  geom_tile(data = current_data, aes(x = lon, y = lat, fill = sink), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradient2(expression(paste('Sink rate \n(2000 - 2010) g C m'^'-2', 'yr'^'-1', sep = '')), low = '#d73027', mid = 'white', high = '#3288bd', na.value="transparent", limits = c(-20, 20), trans = 'identity', oob = scales::squish) +
  # geom_point(data = current_data[current_data$sig_index == 0, ], aes(x = lon, y = lat), shape = 4, size = 0.3, stroke = 0.1, color = 'black') +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
  # change the background to black and white
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  theme(legend.justification = c(0, 0), legend.position = c(0.1, -0.2), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 30, barheight = 2.5, title.position = 'right', title.hjust = 1, label.hjust = 0.5, label.vjust = 2.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 30, ), legend.title = element_text(size = 35)) +
  # add title
  labs(title = 'TRENDY ensemble', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))

p_trendy_sink_rate_std =
  ggplot() +
  geom_tile(data = current_data, aes(x = lon, y = lat, fill = std), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradientn(expression(paste('Standard deviation \n(2000 - 2010) g C m'^'-2', 'yr'^'-1', sep = '')), colours = c('#fef0d9','#fdcc8a','#fc8d59','#e34a33','#b30000'), na.value="transparent", limits = c(0, 20), trans = 'identity', oob = scales::squish) +
  # geom_point(data = current_data[current_data$sig_index == 0, ], aes(x = lon, y = lat), shape = 4, size = 0.3, stroke = 0.1, color = 'black') +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
  # change the background to black and white
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  theme(legend.justification = c(0, 0), legend.position = c(0.1, -0.2), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 30, barheight = 2.5, title.position = 'right', title.hjust = 1, label.hjust = 0.5, label.vjust = 2.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 30, ), legend.title = element_text(size = 35)) +
  # add title
  labs(title = 'TRENDY ensemble', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))

#----------------------------------------------
# PRODA soil carbon sink, total
#----------------------------------------------
current_data = data.frame(cbind(global_lat_lon, 
                                apply(mip_sink_tot, 1, median, na.rm = TRUE),
                                apply(mip_sink_tot, 1, sd, na.rm = TRUE),
                                apply(mip_sink_tot, 1, quantile, prob = 0.25, na.rm = TRUE),
                                apply(mip_sink_tot, 1, quantile, prob = 0.75, na.rm = TRUE), 
                                0))


colnames(current_data) = c('lon', 'lat', 'sink', 'std', 'lower_sigma', 'upper_sigma', 'sig_index')
current_data$sig_index[current_data$sink > 0 & current_data$lower_sigma > 0] = 1
current_data$sig_index[current_data$sink < 0 & current_data$upper_sigma < 0] = 1
# current_data$sig_index[(abs(current_data$sink) -  current_data$std) > 0] = 1


lon_lat_transfer = project(xy = as.matrix(current_data[ , c('lon', 'lat')]), proj = '+proj=robin') 
current_data[ , c('lon', 'lat')] = lon_lat_transfer

lat_limits = rbind(c(0, -56), c(0, 80))
lat_limits_robin = project(xy = as.matrix(lat_limits), proj = '+proj=robin') 

p_sink_total =
  ggplot() +
  geom_tile(data = current_data, aes(x = lon, y = lat, fill = sink), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradient2(expression(paste('Total sink \n(1900 - 2010) g C m'^'-2', sep = '')), low = '#d73027', mid = 'white', high = '#3288bd', na.value="transparent", limits = c(-500, 500), trans = 'identity', oob = scales::squish) +
  # geom_point(data = current_data[current_data$sig_index == 0, ], aes(x = lon, y = lat), shape = 4, size = 0.3, stroke = 0.1, color = 'black') +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
  # change the background to black and white
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  theme(legend.justification = c(0, 0), legend.position = c(0.1, -0.2), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 30, barheight = 2.5, title.position = 'right', title.hjust = 1, label.hjust = 0.5, label.vjust = 2.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 30, ), legend.title = element_text(size = 35)) +
  # add title
  labs(title = 'PRODA ensemble', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))

p_sink_total_std =
  ggplot() +
  geom_tile(data = current_data, aes(x = lon, y = lat, fill = std), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradientn(expression(paste('Standard deviation \n(1900 - 2010) g C m'^'-2', sep = '')), colours = c('#fef0d9','#fdcc8a','#fc8d59','#e34a33','#b30000'), na.value="transparent", limits = c(0, 500), trans = 'identity', oob = scales::squish) +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
  # change the background to black and white
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  theme(legend.justification = c(0, 0), legend.position = c(0.1, -0.2), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 30, barheight = 2.5, title.position = 'right', title.hjust = 1, label.hjust = 0.5, label.vjust = 2.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 30, ), legend.title = element_text(size = 35)) +
  # add title
  labs(title = 'PRODA ensemble', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))

#----------------------------------------------
# PRODA soil carbon sink, rate
#----------------------------------------------
current_data = data.frame(cbind(global_lat_lon, 
                                apply(mip_sink_rate/10, 1, median, na.rm = TRUE),
                                apply(mip_sink_rate/10, 1, sd, na.rm = TRUE),
                                apply(mip_sink_rate/10, 1, quantile, prob = 0.25, na.rm = TRUE),
                                apply(mip_sink_rate/10, 1, quantile, prob = 0.75, na.rm = TRUE), 
                                0))


colnames(current_data) = c('lon', 'lat', 'sink', 'std', 'lower_sigma', 'upper_sigma', 'sig_index')
current_data$sig_index[current_data$sink > 0 & current_data$lower_sigma > 0] = 1
current_data$sig_index[current_data$sink < 0 & current_data$upper_sigma < 0] = 1


lon_lat_transfer = project(xy = as.matrix(current_data[ , c('lon', 'lat')]), proj = '+proj=robin') 
current_data[ , c('lon', 'lat')] = lon_lat_transfer

lat_limits = rbind(c(0, -56), c(0, 80))
lat_limits_robin = project(xy = as.matrix(lat_limits), proj = '+proj=robin') 

p_sink_rate =
  ggplot() +
  geom_tile(data = current_data, aes(x = lon, y = lat, fill = sink), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradient2(expression(paste('Sink rate \n(2000 - 2010) g C m'^'-2', 'yr'^'-1', sep = '')), low = '#d73027', mid = 'white', high = '#3288bd', na.value="transparent", limits = c(-20, 20), trans = 'identity', oob = scales::squish) +
  # geom_point(data = current_data[current_data$sig_index == 0, ], aes(x = lon, y = lat), shape = 4, size = 0.3, stroke = 0.1, color = 'black') +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
  # change the background to black and white
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  # theme(legend.position = 'none') +
  theme(legend.justification = c(0, 0), legend.position = c(0.1, -0.2), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 30, barheight = 2.5, title.position = 'right', title.hjust = 1, label.hjust = 0.5, label.vjust = 2.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 30, ), legend.title = element_text(size = 35)) +
  # add title
  labs(title = 'PRODA ensemble', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))

p_sink_rate_std =
  ggplot() +
  geom_tile(data = current_data, aes(x = lon, y = lat, fill = std), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradientn(expression(paste('Standard deviation \n(2000 - 2010) g C m'^'-2', 'yr'^'-1', sep = '')), colours = c('#fef0d9','#fdcc8a','#fc8d59','#e34a33','#b30000'), na.value="transparent", limits = c(0, 20), trans = 'identity', oob = scales::squish) +
  # geom_point(data = current_data[current_data$sig_index == 0, ], aes(x = lon, y = lat), shape = 4, size = 0.3, stroke = 0.1, color = 'black') +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
  # change the background to black and white
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  theme(legend.justification = c(0, 0), legend.position = c(0.1, -0.2), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 30, barheight = 2.5, title.position = 'right', title.hjust = 1, label.hjust = 0.5, label.vjust = 2.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 30, ), legend.title = element_text(size = 35)) +
  # add title
  labs(title = 'PRODA ensemble', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))

# ----------------------------------------------
# TRENDY plant carbon input changes
#----------------------------------------------
current_data = data.frame(cbind(global_lat_lon, 
                                apply(mip_input_change, 1, median, na.rm = TRUE),
                                apply(mip_input_change, 1, sd, na.rm = TRUE),
                                apply(mip_input_change, 1, quantile, prob = 0.25, na.rm = TRUE),
                                apply(mip_input_change, 1, quantile, prob = 0.75, na.rm = TRUE), 
                                0))


colnames(current_data) = c('lon', 'lat', 'input', 'std', 'lower_sigma', 'upper_sigma', 'sig_index')
current_data$sig_index[current_data$sink > 0 & current_data$lower_sigma > 0] = 1
current_data$sig_index[current_data$sink < 0 & current_data$upper_sigma < 0] = 1


lon_lat_transfer = project(xy = as.matrix(current_data[ , c('lon', 'lat')]), proj = '+proj=robin') 
current_data[ , c('lon', 'lat')] = lon_lat_transfer

lat_limits = rbind(c(0, -56), c(0, 80))
lat_limits_robin = project(xy = as.matrix(lat_limits), proj = '+proj=robin') 

p_input_change =
  ggplot() +
  geom_tile(data = current_data, aes(x = lon, y = lat, fill = input), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradient2(expression(paste('NPP change (1900 - 2010) ', ' g C m'^'-2', 'yr'^'-1', sep = '')), low = '#d73027', mid = 'white', high = '#3288bd', na.value="transparent", limits = c(-200, 200), trans = 'identity', oob = scales::squish) +
  # geom_point(data = current_data[current_data$sig_index == 0, ], aes(x = lon, y = lat), shape = 4, size = 0.3, stroke = 0.1, color = 'black') +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
  # change the background to black and white
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  # theme(legend.position = 'none') +
  theme(legend.justification = c(0, 0), legend.position = c(0.1, -0.2), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 30, barheight = 2.5, title.position = 'right', title.hjust = 1, label.hjust = 0.5, label.vjust = 2.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 30, ), legend.title = element_text(size = 35)) +
  # add title
  labs(title = 'TRENDY ensemble', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))

p_input_change_std =
  ggplot() +
  geom_tile(data = current_data, aes(x = lon, y = lat, fill = std), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradientn(expression(paste('Standard deviation \n(1900 - 2010) g C m'^'-2', 'yr'^'-1', sep = '')), colours = c('#fef0d9','#fdcc8a','#fc8d59','#e34a33','#b30000'), na.value="transparent", limits = c(0, 200), trans = 'identity', oob = scales::squish) +
  # geom_point(data = current_data[current_data$sig_index == 0, ], aes(x = lon, y = lat), shape = 4, size = 0.3, stroke = 0.1, color = 'black') +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
  # change the background to black and white
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  theme(legend.justification = c(0, 0), legend.position = c(0.1, -0.2), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 30, barheight = 2.5, title.position = 'right', title.hjust = 1, label.hjust = 0.5, label.vjust = 2.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 30, ), legend.title = element_text(size = 35)) +
  # add title
  labs(title = 'TRENDY ensemble', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))


jpeg(paste('./Soil_Sink/proda_global_seq_maps.jpeg', sep = ''), width = 22, height = 18, units = 'in', res = 300)
plot_grid(NULL, p_sink_total, NULL, p_box_accu_proda, NULL, 
          NULL, p_sink_rate, NULL,  p_box_rate_proda, NULL,
          nrow = 2,
          ncol = 5,
          rel_widths = c(0.0, 1, 0.0, 0.5, 0.0),
          labels = c(' ', 'a', ' ', 'b', ' ',
                     ' ', 'c', ' ', 'd', ' '), 
          label_size = 50,
          label_x = 0.01, label_y = 0.92,
          label_fontfamily = 'Arial',
          label_fontface = 'bold'
)
dev.off()


jpeg(paste('./Soil_Sink/proda_and_trendy_global_seq_maps.jpeg', sep = ''), width = 38, height = 18, units = 'in', res = 300)
plot_grid(NULL, p_trendy_sink_total, NULL, p_sink_total, NULL, p_box_accu, NULL, 
          NULL, p_trendy_sink_rate, NULL, p_sink_rate, NULL,  p_box_rate, NULL,
          nrow = 2,
          ncol = 7,
          rel_widths = c(0.0, 1, 0.0, 1, 0.0, 0.5, 0.0),
          labels = c(' ', 'a', ' ', 'b', ' ', 'c', ' ',
                     ' ', 'd', ' ', 'e', ' ', 'f', ' '), 
          label_size = 50,
          label_x = 0.01, label_y = 0.92,
          label_fontfamily = 'Arial',
          label_fontface = 'bold'
)
dev.off()


jpeg(paste('./Soil_Sink/trendy_npp_change.jpeg', sep = ''), width = 18, height = 9, units = 'in', res = 300)
p_input_change + 
  theme(legend.justification = c(0, 0), legend.position = c(0.1, -0.13), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 30, barheight = 2.5, title.position = 'right', title.vjust = 1, label.hjust = 0.5, label.vjust = 3, frame.linewidth = 0), reverse = FALSE) +
  labs(title = '', x = '', y = '') 
dev.off()

