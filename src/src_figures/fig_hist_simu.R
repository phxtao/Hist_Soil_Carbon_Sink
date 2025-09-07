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

setwd('/Users/phxtao/Google_Drive/R')
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

data_dir_output = '/Users/phxtao/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/phxtao/DATAHUB/ENSEMBLE/INPUT_DATA/'

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
cross_valid_num = 1
cpool_num = 7
year_start = 1900
year_end = 2010
year_series = seq(from = year_start, to = year_end, by = 10)
year_num = length(year_series)
control_gradient = seq(from = 0, to = 1, by = 0.1)

soc_hist_raw = array(NA, dim = c(12, 28, 10))
resp_hist_raw = array(NA, dim = c(12, 28, 10))
tau_hist_raw = array(NA, dim = c(12, 28, 10))
process_hist_raw = array(NA, dim = c(12, 21, 10))

for (cross_valid_num in 1:10) {
#-----------------scenario simulation
soc_hist_raw_middle = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_scenario_simu_soc_tot.mat', sep = ''))
soc_hist_raw[ , , cross_valid_num] = soc_hist_raw_middle$global.soc.sink.tot[ , , 4]

resp_hist_raw_middle = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_scenario_simu_resp_tot.mat', sep = ''))
resp_hist_raw[ , , cross_valid_num] = resp_hist_raw_middle$global.soc.resp.tot[ , , 4]

tau_hist_raw_middle = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_scenario_simu_tau_tot.mat', sep = ''))
tau_hist_raw[ , , cross_valid_num] = tau_hist_raw_middle$global.soc.tau.tot[ , , 4]

process_hist_raw_middle = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_scenario_simu_process_tot.mat', sep = ''))
process_hist_raw[ , , cross_valid_num] = process_hist_raw_middle$global.process.tot[ , , 4]
}

soil_depth_list = c('0 - 30cm', '30 - 100cm', '100 - 200cm', '200cm - bedrock')
depth_color_scheme = c('#8c510a','#dfc27d','#80cdc1','#01665e')

#################################################################################
# depth resolved soil carbon sink
#################################################################################
current_data = c()

idepth = 1
for (idepth in 1:4) {
  soc_change = apply(soc_hist_raw[ , seq(from = (16+idepth), to = 28, by = 4), ], c(1, 3), sum, na.rm = TRUE)
  
  tau_change = tau_hist_raw[ , seq(from = (16+idepth), to = 28, by = 4), ]
  
    soc_weight = soc_hist_raw[ , seq(from = (16+idepth), to = 28, by = 4), ]/aperm(matlab::repmat(soc_change, c(1, 1, 3)), c(1, 3, 2))
    tau_change_control = apply(tau_change*soc_weight, c(1, 3), sum, na.rm = TRUE)

        # tau_change_control = apply(tau_change[ , , imanagement], c(1), sum, na.rm = TRUE)
    current_data = rbind(current_data, cbind(apply(soc_change - matlab::repmat(soc_change[1, ], c(12, 1)), 1, mean, na.rm = TRUE), 
                                             apply(soc_change - matlab::repmat(soc_change[1, ], c(12, 1)), 1, sd, na.rm = TRUE), 
                                             apply(tau_change_control, 1, mean, na.rm = TRUE), 
                                             apply(tau_change_control, 1, sd, na.rm = TRUE), 
                                             year_series, idepth))
}


colnames(current_data) = c('sink', 'sink_std', 'tau', 'tau_std', 'year', 'depth')
current_data = data.frame(current_data)

current_data$sink_error_base = NA


icontrol = 1
idepth = 1
for (iyear in 1:length(year_series)) {
  for (idepth in 1:4) {
    data_loc = which(current_data$depth == idepth & current_data$year == year_series[iyear])
    if (current_data$sink[data_loc] > 0){
      positive_loc = which(current_data$depth >= idepth 
                           & current_data$year == year_series[iyear]
                           & current_data$sink > 0)
      
      current_data$sink_error_base[data_loc] = sum(current_data$sink[positive_loc])
    } else {
      negative_loc = which(current_data$depth >= idepth 
                           & current_data$year == year_series[iyear]
                           & current_data$sink <= 0)
      current_data$sink_error_base[data_loc] = sum(current_data$sink[negative_loc])
    }
  }
}

#----------------------------------
# figure
#----------------------------------
current_data_base = current_data
current_data_year = apply(soc_hist_raw[ , 17:28, ], c(1, 3), sum, na.rm = TRUE)
current_data_year = cbind(year_series, 
                          apply(current_data_year - matlab::repmat(current_data_year[1, ], c(12, 1)), 1, mean, na.rm = TRUE),
                          apply(current_data_year - matlab::repmat(current_data_year[1, ], c(12, 1)), 1, sd, na.rm = TRUE)
                          )
colnames(current_data_year) = c('year', 'sink', 'sink_std')
current_data_year = data.frame(current_data_year)

p_sink_depth =
ggplot() +
  geom_bar(data = current_data, aes(x = year, y = sink, fill = as.factor(depth)), color = 'white', width = 9.5, stat = 'identity', position = 'stack', size = 0.5) +
  geom_errorbar(data = current_data, aes(x = year, y = sink_error_base, ymin = sink_error_base - 2*sink_std, ymax = sink_error_base + 2*sink_std), stat = 'identity', width = 2, size = 1) +
  geom_ribbon(data = current_data_year, aes(x = year, y = sink, ymin = sink - 2*sink_std, ymax = sink + 2*sink_std), fill = "black", alpha = 0.5) +
  geom_line(data = current_data_year, aes(x = year, y = sink), size = 1.5) +
  geom_hline(yintercept = 0, color = 'black', size = 1) +
  scale_fill_manual(name = c(''), labels = soil_depth_list, values = depth_color_scheme) +
  # scale_x_continuous(n.breaks = 7, sec.axis = sec_axis(trans = ~.))
  # scale_y_continuous(limits = c(-30, 100), n.breaks = 7) +
  theme_classic() +
  # add title
  labs(title = '', x = '', y = expression(paste('Stock change (Pg C)', sep = ''))) +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 0), legend.direction = 'vertical')  +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1), legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.02, vjust = 0, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch')) +
  theme(panel.background = element_rect(fill = 'transparent', color = NA),  # plot panel
        plot.background  = element_rect(fill = 'transparent', color = NA),  # entire plot
        legend.background = element_rect(fill = 'transparent', color = NA), # legend bg
        legend.box.background = element_rect(fill = 'transparent', color = NA)) # legend box




ilayer = 1
for (ilayer in 1:4) {
  if (ilayer == 1) {
    p =
      ggplot() +
      geom_ribbon(data = current_data_base[current_data_base$depth == ilayer, ], aes(x = year, y = tau, ymin = tau - 2*tau_std, ymax = tau + 2*tau_std), fill = depth_color_scheme[ilayer], alpha = 0.5) +
      geom_line(data = current_data_base[current_data_base$depth == ilayer, ], aes(x = year, y = tau), color = depth_color_scheme[ilayer], linetype = 'solid', size = 2.5) +
      geom_point(data = current_data_base[current_data_base$depth == ilayer, ], aes(x = year, y = tau), color = depth_color_scheme[ilayer], shape = 18, size = 10) +
      # geom_point(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = year, y = net_flux), color = depth_color_scheme[ilayer], size = 5, stroke = 5) +
      scale_y_continuous(limits = c(NA, NA), n.breaks = 4, trans = 'identity') +
      coord_cartesian(clip = 'off') + 
      theme_classic() +
      # add title
      labs(title = '', x = '', y = expression(paste('Residence time (year)'))) +
      # change the legend properties
      theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'vertical', legend.box = "horizontal")  +
      theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
      theme(legend.justification = c(0, 0), legend.position = 'None', legend.background = element_rect(fill = NA)) +
      theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
      # modify the position of title
      theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 40)) +
      # modify the font size
      # modify the margin
      theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
      theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
      theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch')) +
      theme(axis.title.y = element_text(hjust = -0.6))
  } else {
    p = 
      ggplot() +
      geom_ribbon(data = current_data_base[current_data_base$depth == ilayer, ], aes(x = year, y = tau, ymin = tau - 2*tau_std, ymax = tau + 2*tau_std), fill = depth_color_scheme[ilayer], alpha = 0.5) +
      geom_line(data = current_data_base[current_data_base$depth == ilayer, ], aes(x = year, y = tau), color = depth_color_scheme[ilayer], linetype = 'solid', size = 2.5) +
      geom_point(data = current_data_base[current_data_base$depth == ilayer, ], aes(x = year, y = tau), color = depth_color_scheme[ilayer], shape = 18, size = 10) +
      # geom_point(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = year, y = net_flux), color = depth_color_scheme[ilayer], size = 5, stroke = 5) +
      scale_y_continuous(limits = c(NA, NA), n.breaks = 4, trans = 'identity') +
      coord_cartesian() + 
      theme_classic() +
      # add title
      labs(title = '', x = '', y = '') +
      # change the legend properties
      theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 40)) +
      # modify the font size
      # modify the margin
      theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
      theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
      theme(axis.text.y = element_text(size = 30, color = 'black'), axis.text.x = element_text(size = 0, color = NULL),
            axis.title.y = element_text(size = 35), axis.title.x = element_text(size = 0), 
            axis.line.y = element_line(size = 1), axis.line.x = element_line(size = 0), 
            axis.ticks.y = element_line(size = 1, color = 'black'), axis.ticks.x = element_line(size = 1, color = NULL), 
            axis.ticks.length.y = unit(0.12, 'inch'), axis.ticks.length.x = unit(0, 'inch'))
  }
  
  eval(parse(text = paste('p', ilayer, ' = p', sep = '')))
  
}

p_tau_depth =
  plot_grid(
    p4, NULL,
    p3, NULL, 
    p2, NULL, 
    p1, NULL,
    ncol = 2,
    nrow = 4, 
    rel_heights = c(1, 1, 1, 1.4),
    rel_widths = c(1, 0.1),
    greedy = FALSE,
    align = 'v'
  )


#################################################################################
# depth resolved  process change
#################################################################################
component_name_list = c('Carbon input', 'Soil respiration', 'Carbon transfer', 'Baseline Decomposition', 'Vertical transport', 'Environmental modification')

current_data = c()
idepth = 1
for (idepth in 1:4) {
  soc_sink = soc_hist_raw[ , seq(from = (16+idepth), to = 28, by = 4), ]
  soc_sink = apply(soc_sink, c(1, 3), sum, na.rm = TRUE)
  soc_sink_rate = rbind(0, apply(soc_sink, 2, diff))/10

  resp_change = apply(resp_hist_raw[1:12, seq(from = (16+idepth), to = 28, by = 4), ], c(1, 3), sum, na.rm = TRUE)
  
  input_change = soc_sink_rate + resp_change
  current_data = rbind(current_data, 
                       cbind(year_series, 
                             apply(input_change, 1, mean, na.rm = TRUE),
                             apply(input_change, 1, sd, na.rm = TRUE),
                             0.1,
                             idepth),
                       cbind(year_series, 
                             apply(resp_change, 1, mean, na.rm = TRUE),
                             apply(resp_change, 1, sd, na.rm = TRUE),
                             0.2,
                             idepth),
                       cbind(year_series, 
                             apply((input_change - resp_change), 1, mean, na.rm = TRUE),
                             apply((input_change - resp_change), 1, sd, na.rm = TRUE),
                             0.3,
                             idepth)
                       
  )
}

current_data = data.frame(current_data)
colnames(current_data) = c('year', 'process', 'process_std', 'process_name', 'depth')

#---------------------------------------------------
current_data_plot = cbind(current_data[current_data$process_name == 0.1, c('year', 'depth')], 
                          current_data[current_data$process_name == 0.1, 'process'],
                          current_data[current_data$process_name == 0.1, 'process_std'],
                          current_data[current_data$process_name == 0.2, 'process'],
                          current_data[current_data$process_name == 0.2, 'process_std'],
                          current_data[current_data$process_name == 0.3, 'process'],
                          current_data[current_data$process_name == 0.3, 'process_std']
)
colnames(current_data_plot) = c('year', 'depth', 'influx', 'influx_std', 'outflux', 'outflux_std', 'net_flux', 'net_flux_std')

ilayer = 4
for (ilayer in 1:4) {
  if (ilayer == 4) {
    p =
      ggplot() +
      geom_ribbon(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = year, y = influx*100, ymin = (influx - 2*influx_std)*100, ymax = (influx + 2*influx_std)*100), fill = depth_color_scheme[ilayer], alpha = 0.5) +
      geom_ribbon(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = year, y = outflux*100, ymin = (outflux - 2*outflux_std)*100, ymax = (outflux + 2*outflux_std)*100), fill = depth_color_scheme[ilayer], alpha = 0.5) +
      geom_line(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = year, y = influx*100), color = depth_color_scheme[ilayer], linetype = 'solid', size = 2.5) +
      geom_line(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = year, y = outflux*100), color = depth_color_scheme[ilayer], linetype = 'dashed', size = 2.5)  +
    # geom_point(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = year, y = net_flux), color = depth_color_scheme[ilayer], size = 5, stroke = 5) +
      scale_y_continuous(limits = c(NA, NA), n.breaks = 4, trans = 'identity') +
      coord_cartesian(clip = 'off') + 
      theme_classic() +
      annotate('text', label = expression(paste('×10'^'-2', sep = '')), x = 1900, y = 1.2, hjust = 0.5, vjust = -0.5,  size = 10) + 
      # add title
      labs(title = '', x = '', y = expression(paste('In- and out-fluxes (Pg C yr'^-1, ')'))) +
      # change the legend properties
      theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'vertical', legend.box = "horizontal")  +
      theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
      theme(legend.justification = c(0, 0), legend.position = 'None', legend.background = element_rect(fill = NA)) +
      theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
      # modify the position of title
      theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 40)) +
      # modify the font size
      # modify the margin
      theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
      theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
      theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch')) +
      theme(axis.title.y = element_text(hjust = -0.3))
  } else {
    p = 
      ggplot() +
      geom_ribbon(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = year, y = influx, ymin = (influx - 2*influx_std), ymax = (influx + 2*influx_std)), fill = depth_color_scheme[ilayer], alpha = 0.5) +
      geom_ribbon(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = year, y = outflux, ymin = (outflux - 2*outflux_std), ymax = (outflux + 2*outflux_std)), fill = depth_color_scheme[ilayer], alpha = 0.5) +
      geom_line(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = year, y = influx), color = depth_color_scheme[ilayer], linetype = 'solid', size = 2.5) +
      geom_line(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = year, y = outflux), color = depth_color_scheme[ilayer], linetype = 'dashed', size = 2.5) +
      # geom_point(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = year, y = net_flux), color = depth_color_scheme[ilayer], size = 5, stroke = 5) +
      scale_y_continuous(limits = c(NA, NA), n.breaks = 4, trans = 'identity') +
      coord_cartesian() + 
      theme_classic() +
      # add title
      labs(title = '', x = '', y = '') +
      # change the legend properties
      theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 40)) +
      # modify the font size
      # modify the margin
      theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
      theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
      theme(axis.text.y = element_text(size = 30, color = 'black'), axis.text.x = element_text(size = 0, color = NULL),
            axis.title.y = element_text(size = 35), axis.title.x = element_text(size = 0), 
            axis.line.y = element_line(size = 1), axis.line.x = element_line(size = 0), 
            axis.ticks.y = element_line(size = 1, color = 'black'), axis.ticks.x = element_line(size = 1, color = NULL), 
            axis.ticks.length.y = unit(0.12, 'inch'), axis.ticks.length.x = unit(0, 'inch'))
  }
  
  eval(parse(text = paste('p', ilayer, ' = p', sep = '')))
  
}

p_layer_flux =
  plot_grid(
    p1, NULL,
    p2, NULL, 
    p3, NULL, 
    p4, NULL,
    ncol = 2,
    nrow = 4, 
    rel_heights = c(1, 1, 1, 1.4),
    rel_widths = c(1, 0.1),
    greedy = FALSE,
    align = 'v'
  )

p_flux_change =
  ggplot() +
    geom_ribbon(data = current_data_plot, aes(x = year, y = net_flux, ymin = (net_flux - 2*net_flux_std), ymax = (net_flux + 2*net_flux_std), fill = as.factor(depth)), alpha = 0.5) +
  geom_hline(yintercept = 0, color = 'grey', size = 2, linetype = 'dashed') +
  geom_line(data = current_data_plot, aes(x = year, y = net_flux, color = as.factor(depth)), size = 2.5) +
  geom_point(data = current_data_plot, aes(x = year, y = net_flux, color = as.factor(depth)), size = 5, stroke = 5) +
  scale_y_continuous(limits = c(NA, NA), n.breaks = 7) +
  scale_shape_manual(name = '', values = c(15, 16), label = c('Carbon input', 'Soil Respiration')) + 
    scale_color_manual(name = '', values = depth_color_scheme, label = soil_depth_list) + 
    scale_fill_manual(name = '', values = depth_color_scheme, label = soil_depth_list) + 
    coord_cartesian() + 
  theme_classic() +
  # add title
  labs(title = '', x = '', y = expression(paste('Net flux changes (Pg C yr'^-1, ')'))) +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'vertical', legend.box = "horizontal")  +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 0), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))


jpeg(paste('./Soil_Sink/proda_global_sink.jpeg', sep = ''), width = 25, height = 20, units = 'in', res = 300)
p_flux_change = p_flux_change + labs(x = 'Year')
p_tau_depth = p_tau_depth + labs(x = 'Year')

plot_grid(
  NULL, p_sink_depth, p_layer_flux, NULL,
  NULL, p_flux_change, p_tau_depth, NULL, 
  nrow = 2, 
  labels = c(' ', 'a', 'b', ' ',
             ' ', 'c', 'd', ' '), 
  rel_widths = c(0.03, 1, 1, 0.03),
  label_size = 50,
  label_x = 0., label_y = 1.01,
  label_fontfamily = 'Arial',
  label_fontface = 'bold',
  align = 'None'
)
dev.off()

