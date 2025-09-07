## Packages
library(R.matlab)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(viridis)
library(raster)
library(scales)
# dev.off()
##
rm(list = ls())

setwd('/Users/ft254/Google_Drive/R')
# Custom log transformation
sym_log_trans <- function(base = exp(1), thresh = 1, scale = 1) {
  trans <- function(x) sign(x) * log(abs(x) / thresh + 1, base = base) * scale
  inv <- function(x) sign(x) * (base^abs(x / scale) - 1) * thresh
  
  trans_new(paste0("sym_log-", format(base)), trans, inv, 
            domain = c(-Inf, Inf))
}

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
# Load Projected SOC PRODA
#################################################################################
vert_component = 1
x_axis_label = c('Reallocated aboveground input (%)', 'Enhanced vertical transport (%)')
cross_valid_num = 7 
cpool_num = 7
year_start = 1900
year_end = 2010
year_series = seq(from = year_start, to = year_end, by = 10)
year_num = length(year_series)
control_gradient = seq(from = 0, to = 1, by = 0.1)
#-----------------scenario simulation
soc_hist_raw = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_response_curve_soc_tot.mat', sep = ''))
soc_hist_raw = soc_hist_raw$global.soc.sink.tot

resp_hist_raw = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_response_curve_resp_tot.mat', sep = ''))
resp_hist_raw = resp_hist_raw$global.soc.resp.tot

tau_hist_raw = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_response_curve_tau_tot.mat', sep = ''))
tau_hist_raw = tau_hist_raw$global.soc.tau.tot

process_hist_raw = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_response_curve_process_tot.mat', sep = ''))
process_hist_raw = process_hist_raw$global.process.tot


soil_depth_list = c('0 - 30cm', '30 - 100cm', '100 - 200cm', '200cm - bedrock')
depth_color_scheme = c('#8c510a','#dfc27d','#80cdc1','#01665e')

#################################################################################
# depth resolved  process change
#################################################################################
component_name_list = c('Carbon input', 'Soil respiration', 'Carbon transfer', 'Baseline Decomposition', 'Vertical transport', 'Environmental modification')

current_data = c()
idepth = 4
for (idepth in 1:4) {
  soc_sink = soc_hist_raw[ , seq(from = (16+idepth), to = 28, by = 4), 1, vert_component]
  soc_sink = apply(soc_sink, c(1), sum, na.rm = TRUE)
  soc_sink_rate = c(NA, diff(soc_sink))/10
  # soc_sink_rate[1] = soc_sink_rate[2]
  soc_sink_rate[1] = 0
  
  resp_change = apply(resp_hist_raw[1:12, seq(from = (16+idepth), to = 28, by = 4), 1, vert_component], c(1), sum, na.rm = TRUE)

  input_change = soc_sink_rate + resp_change
  current_data = rbind(current_data, 
                       cbind(year_series, 
                             input_change,
                             (input_change - input_change[1])/input_change[1], 
                             0.1,
                             idepth),
                       cbind(year_series, 
                             resp_change,
                             (resp_change - resp_change[1])/resp_change[1], 
                             0.2,
                             idepth)
                       )
  
  
  
  process_change = process_hist_raw[ , seq(from = idepth, to = 20, by = 4), 1, vert_component]
  iprocess = 1
  for (iprocess in c(1, 3, 4, 5)){
    current_data = rbind(current_data, cbind(year_series, 
                                             process_change[ , iprocess],
                                             (process_change[ , iprocess] - process_change[1, iprocess])/process_change[1, iprocess], 
                                             iprocess, 
                                             idepth))
  }
}

current_data = data.frame(current_data)
colnames(current_data) = c('year', 'process', 'process_change', 'process_name', 'depth')

#################################################################################
# env modifier different depth
#################################################################################
current_data_plot = current_data[current_data$process_name == 5, ]
soil_temp_mean = readMat(paste(data_dir_output, 'hist_simu_proda/global_env_change_soil_temp_mean.mat', sep = ''))
soil_temp_mean = soil_temp_mean$global.soil.temp.mean

soil_temp_std = readMat(paste(data_dir_output, 'hist_simu_proda/global_env_change_soil_temp_std.mat', sep = ''))
soil_temp_std = soil_temp_std$global.soil.temp.std

soil_water_mean = readMat(paste(data_dir_output, 'hist_simu_proda/global_env_change_soil_water_mean.mat', sep = ''))
soil_water_mean = soil_water_mean$global.soil.water.mean

soil_water_std = readMat(paste(data_dir_output, 'hist_simu_proda/global_env_change_soil_water_std.mat', sep = ''))
soil_water_std = soil_water_std$global.soil.water.std

soil_depth_scalar = readMat(paste(data_dir_output, 'hist_simu_proda/global_env_change_soil_depth_scalar.mat', sep = ''))
soil_depth_scalar = soil_depth_scalar$global.soil.depth.scalar

env_change_plot = c()
for (idepth in 1:4) {
  env_change_plot = rbind(env_change_plot, 
                          cbind(seq(from = 1900, to = 2010), 
                          soil_temp_mean[idepth, ], 
                          soil_temp_std[idepth, ],
                          soil_water_mean[idepth, ],
                          soil_water_std[idepth, ],
                          soil_depth_scalar[idepth, ], 
                          idepth)
  )
}

colnames(env_change_plot) = c('year', 'temp_mean', 'temp_std', 'water_mean', 'water_std', 'depth_scalar', 'depth')
env_change_plot = data.frame(env_change_plot)

p_temp_mean = 
  ggplot() +
  geom_line(data = env_change_plot, aes(x = year, y = temp_mean, color = as.factor(depth), group = as.factor(depth)), size = 2.5) +
  scale_color_manual(name = '', values = depth_color_scheme, label = soil_depth_list) +
  
  # scale_color_gradientn(name = c('Soil carbon sink'), limits = c(-20, 20), colors = color_gradient_scheme, oob = scales::squish) +
  scale_y_continuous(limits = c(NA, NA), trans = 'identity') +
  # scale_shape_manual(name = 'Components', values = c(15, 16, 17, 18, 8), label = component_name_list) +
  coord_cartesian() + 
  theme_classic() +
  # add title
  labs(title = '', x = 'Year', y = 'Soil temperature (K)') +
  # change the legend properties
  guides(color = guide_legend(ncol = 2)) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'horizontal', legend.box = 'horizontal')  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))

p_temp_cv = 
ggplot() +
  geom_line(data = env_change_plot, aes(x = year, y = temp_std/temp_mean, color = as.factor(depth), group = as.factor(depth)), size = 2.5) +
  scale_color_manual(name = '', values = depth_color_scheme, label = soil_depth_list) +
  
  # scale_color_gradientn(name = c('Soil carbon sink'), limits = c(-20, 20), colors = color_gradient_scheme, oob = scales::squish) +
  scale_y_continuous(limits = c(NA, NA), trans = 'identity') +
  # scale_shape_manual(name = 'Components', values = c(15, 16, 17, 18, 8), label = component_name_list) +
  coord_cartesian() + 
  theme_classic() +
  # add title
  labs(title = '', x = 'Year', y = 'Absolute CV of soil temperature') +
  # change the legend properties
  guides(color = guide_legend(ncol = 2)) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'horizontal', legend.box = 'horizontal')  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))

p_water_mean = 
ggplot() +
  geom_line(data = env_change_plot, aes(x = year, y = water_mean, color = as.factor(depth), group = as.factor(depth)), size = 2.5) +
  scale_color_manual(name = '', values = depth_color_scheme, label = soil_depth_list) +
  
  # scale_color_gradientn(name = c('Soil carbon sink'), limits = c(-20, 20), colors = color_gradient_scheme, oob = scales::squish) +
  scale_y_continuous(limits = c(NA, NA), trans = 'identity') +
  # scale_shape_manual(name = 'Components', values = c(15, 16, 17, 18, 8), label = component_name_list) +
  coord_cartesian() + 
  theme_classic() +
  # add title
  labs(title = '', x = 'Year', y = 'Soil water potential (kPa)') +
  # change the legend properties
  guides(color = guide_legend(ncol = 2)) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'horizontal', legend.box = 'horizontal')  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))

p_water_cv = 
  ggplot() +
  geom_line(data = env_change_plot, aes(x = year, y = abs(water_std/water_mean), color = as.factor(depth), group = as.factor(depth)), size = 2.5) +
  scale_color_manual(name = '', values = depth_color_scheme, label = soil_depth_list) +
  
  # scale_color_gradientn(name = c('Soil carbon sink'), limits = c(-20, 20), colors = color_gradient_scheme, oob = scales::squish) +
  scale_y_continuous(limits = c(NA, NA), trans = 'identity') +
  # scale_shape_manual(name = 'Components', values = c(15, 16, 17, 18, 8), label = component_name_list) +
  coord_cartesian() + 
  theme_classic() +
  # add title
  labs(title = '', x = 'Year', y = 'Absolute CV of soil water potential') +
  # change the legend properties
  guides(color = guide_legend(ncol = 2)) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'horizontal', legend.box = 'horizontal')  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))


p_depth_scalar = 
  ggplot() +
  geom_line(data = env_change_plot, aes(x = year, y = depth_scalar, color = as.factor(depth), group = as.factor(depth)), size = 2.5) +
  scale_color_manual(name = '', values = depth_color_scheme, label = soil_depth_list) +
  
  # scale_color_gradientn(name = c('Soil carbon sink'), limits = c(-20, 20), colors = color_gradient_scheme, oob = scales::squish) +
  scale_y_continuous(limits = c(NA, NA), trans = 'identity') +
  # scale_shape_manual(name = 'Components', values = c(15, 16, 17, 18, 8), label = component_name_list) +
  coord_cartesian() + 
  theme_classic() +
  # add title
  labs(title = '', x = 'Year', y = 'Soil depth scalar (unitless)') +
  # change the legend properties
  guides(color = guide_legend(ncol = 2)) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'horizontal', legend.box = 'horizontal')  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))


p_total_scalar = 
ggplot() +
  geom_hline(yintercept = 0, color = 'grey', size = 2, linetype = 'dashed') +
  geom_line(data = current_data_plot, aes(x = year, y = process, color = as.factor(depth), group = as.factor(depth)), size = 2.5) +
  geom_point(data = current_data_plot, aes(x = year, y = process, color = as.factor(depth), group = as.factor(depth)), size = 5, stroke = 5) +
  scale_color_manual(name = '', values = depth_color_scheme, label = soil_depth_list) +
  
  # scale_color_gradientn(name = c('Soil carbon sink'), limits = c(-20, 20), colors = color_gradient_scheme, oob = scales::squish) +
  scale_y_continuous(limits = c(0.0, 0.7), trans = 'identity') +
  # scale_shape_manual(name = 'Components', values = c(15, 16, 17, 18, 8), label = component_name_list) +
  coord_cartesian() + 
  theme_classic() +
  # add title
  labs(title = '', x = 'Year', y = 'Environmental stress') +
  # change the legend properties
  guides(color = guide_legend(ncol = 2)) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'horizontal', legend.box = 'horizontal')  +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))

jpeg(paste('./Soil_Sink/env_modifier_across_depth.jpeg', sep = ''), width = 25, height = 30, units = 'in', res = 300)
plot_grid(
  NULL, p_total_scalar, p_depth_scalar, NULL, 
  NULL, p_temp_mean, p_temp_cv, NULL,
  NULL, p_water_mean, p_water_cv, NULL,
  nrow = 3,
  labels = c(' ', 'a', 'b', ' ',
             ' ', 'c', 'd', ' ',
             ' ', 'e', 'f', ' '
  ),
  rel_widths = c(0.02, 1, 1, 0.02),
  label_size = 50,
  label_x = 0., label_y = 1.01,
  label_fontfamily = 'Arial',
  label_fontface = 'bold'
)
dev.off()

