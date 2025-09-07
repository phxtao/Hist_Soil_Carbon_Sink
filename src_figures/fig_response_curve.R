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
# response curve
#################################################################################
vert_component = 1

for (vert_component in 1) {
  x_axis_label = c('Reallocated aboveground input (%)', 'Enhanced vertical transport (%)')
  cpool_num = 7
  year_start = 1900
  year_end = 2010
  year_series = seq(from = year_start, to = year_end, by = 10)
  year_num = length(year_series)
  control_gradient = seq(from = 0, to = 1, by = 0.1)
  
  #################################################################################
  # depth resolved soil carbon sink
  #################################################################################
  current_data_net_control_raw = array(0, dim = c(length(control_gradient), 10))
  current_data_raw = array(0, dim = c(528, 5, 10))
  current_data_flux_raw  = array(NA, dim = c(220, 7, 10))
  
  cross_valid_num = 1
  for (cross_valid_num in 1:10){
    # base simu
    soc_hist_raw_base = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_scenario_simu_soc_tot.mat', sep = ''))
    soc_hist_raw_base = soc_hist_raw_base$global.soc.sink.tot[ , , 4]
    
    resp_hist_raw_base = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_scenario_simu_resp_tot.mat', sep = ''))
    resp_hist_raw_base = resp_hist_raw_base$global.soc.resp.tot[ , , 4]
    
    tau_hist_raw_base = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_scenario_simu_tau_tot.mat', sep = ''))
    tau_hist_raw_base = tau_hist_raw_base$global.soc.tau.tot[ , , 4]
    
    process_hist_raw_base = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_scenario_simu_process_tot.mat', sep = ''))
    process_hist_raw_base = process_hist_raw_base$global.process.tot[ , , 4]
    
    
    # response curve
    soc_hist_raw = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_response_curve_soc_tot.mat', sep = ''))
    soc_hist_raw = soc_hist_raw$global.soc.sink.tot
    soc_hist_raw[ , , 1, 1] = soc_hist_raw_base
    soc_hist_raw[ , , 1, 2] = soc_hist_raw_base
    
    resp_hist_raw = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_response_curve_resp_tot.mat', sep = ''))
    resp_hist_raw = resp_hist_raw$global.soc.resp.tot
    resp_hist_raw[ , , 1, 1] = resp_hist_raw_base
    resp_hist_raw[ , , 1, 2] = resp_hist_raw_base
    
    tau_hist_raw = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_response_curve_tau_tot.mat', sep = ''))
    tau_hist_raw = tau_hist_raw$global.soc.tau.tot
    tau_hist_raw[ , , 1, 1] = tau_hist_raw_base
    tau_hist_raw[ , , 1, 2] = tau_hist_raw_base
    
    process_hist_raw = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_response_curve_process_tot.mat', sep = ''))
    process_hist_raw = process_hist_raw$global.process.tot
    process_hist_raw[ , , 1, 1] = process_hist_raw_base
    process_hist_raw[ , , 1, 2] = process_hist_raw_base
    
    #---------------flux changes
    current_data_flux_middle = c()
    idepth = 2
    for (idepth in 1:4) {
      soc_sink = soc_hist_raw[12, seq(from = (16+idepth), to = 28, by = 4), , vert_component] - soc_hist_raw[1, seq(from = (16+idepth), to = 28, by = 4), , vert_component]
      soc_sink = apply(soc_sink, c(2), sum, na.rm = TRUE)
      
      resp_change = apply(apply(resp_hist_raw[1:12, seq(from = (16+idepth), to = 28, by = 4), , vert_component], c(2, 3), mean, na.rm = TRUE), c(2), sum, na.rm = TRUE)
      # resp_change = resp_change/resp_change[1]
      
      input_change = soc_sink/110 + resp_change
      # input_change = apply(process_hist_raw[1:12, (4+idepth), , 1], c(2), mean, na.rm = TRUE)
      # input_change = input_change/input_change[1]
      
      # process_hist_raw[ , (5+idepth-1), , ] = process_hist_raw[ , (5+idepth-1), , ]*process_hist_raw[ , 20, , ]
      process_change = (process_hist_raw[ , seq(from = idepth, to = 20, by = 4), , vert_component] - matlab::repmat(process_hist_raw[ , seq(from = idepth, to = 20, by = 4), 1, vert_component], c(1, 1, 11)))/matlab::repmat(process_hist_raw[ , seq(from = idepth, to = 20, by = 4), 1, vert_component], c(1, 1, 11))
      
      process_change = apply(process_change, c(2, 3), mean, na.rm = TRUE)
      
      iprocess = 1
      for (iprocess in 1:5){
        current_data_flux_middle = rbind(current_data_flux_middle, cbind(control_gradient, soc_sink, resp_change, input_change, process_change[iprocess, ], iprocess, idepth))
      }
    }
    current_data_flux_raw[ , , cross_valid_num] = current_data_flux_middle
    
    # data summary
    current_data_middle = c()
    
    idepth = 1
    for (idepth in 1:4) {
      soc_change = apply(soc_hist_raw[ , seq(from = (16+idepth), to = 28, by = 4), , vert_component], c(1, 3), sum, na.rm = TRUE)
      
      current_data_net_control_raw[ , cross_valid_num] = current_data_net_control_raw[ , cross_valid_num] + (soc_change[year_num, ] - soc_change[1, ])
      
      tau_change = tau_hist_raw[ , seq(from = (16+idepth), to = 28, by = 4), , vert_component]
      
      imanagement = 1
      for (imanagement in 1:11) {
        soc_weight = soc_hist_raw[ , seq(from = (16+idepth), to = 28, by = 4), imanagement, vert_component]/matlab::repmat(soc_change[ , imanagement], c(1, 3))
        tau_change_control = apply(tau_change[ , , imanagement]*soc_weight, c(1), sum, na.rm = TRUE)
        # tau_change_control = apply(tau_change[ , , imanagement], c(1), sum, na.rm = TRUE)
        current_data_middle = rbind(current_data_middle, cbind(soc_change[ , imanagement] - soc_change[1, imanagement], tau_change_control, year_series, idepth, control_gradient[imanagement]))
      }
    }
    current_data_raw[ , , cross_valid_num] = current_data_middle
  }
  
  current_data_net_control = data.frame(cbind(
    control_gradient,
    apply(current_data_net_control_raw, 1, mean, na.rm = TRUE),
    apply(current_data_net_control_raw, 1, sd, na.rm = TRUE)
  ))
  
  colnames(current_data_net_control) = c('control', 'sink', 'sink_std')
  
  
  current_data = cbind(apply(current_data_raw[ , 1, ], 1, mean, na.rm = TRUE),
                       apply(current_data_raw[ , 1, ], 1, sd, na.rm = TRUE),
                       apply(current_data_raw[ , 2, ], 1, mean, na.rm = TRUE),
                       apply(current_data_raw[ , 2, ], 1, sd, na.rm = TRUE),
                       current_data_raw[ , 3, 1],
                       current_data_raw[ , 4, 1],
                       current_data_raw[ , 5, 1]
  )
  colnames(current_data) = c('sink', 'sink_std', 'tau', 'tau_std', 'year', 'depth', 'control')
  current_data = data.frame(current_data)
  
  
  sum(current_data[current_data$year == 2010 & current_data$control == 1, 1]) - 
    sum(current_data[current_data$year == 2000 & current_data$control == 1, 1])
  #----------------------------------
  
  soil_depth_list = c('0 - 30cm', '30 - 100cm', '100 - 200cm', '200cm - bedrock')
  depth_color_scheme = c('#8c510a','#dfc27d','#80cdc1','#01665e')
  
  
  current_data_base = current_data[current_data$year == 2010, ]
  current_data_base$sink_error_base = NA
  
  
  icontrol = 1
  idepth = 1
  for (icontrol in 1:11) {
    for (idepth in 1:4) {
      data_loc = which(current_data_base$depth == idepth & current_data_base$control == control_gradient[icontrol])
      if (current_data_base$sink[data_loc] > 0){
        positive_loc = which(current_data_base$depth >= idepth 
                             & current_data_base$control == control_gradient[icontrol]
                             & current_data_base$sink > 0)
        
        current_data_base$sink_error_base[data_loc] = sum(current_data_base$sink[positive_loc])
      } else {
        negative_loc = which(current_data_base$depth >= idepth 
                             & current_data_base$control == control_gradient[icontrol]
                             & current_data_base$sink <= 0)
        current_data_base$sink_error_base[data_loc] = sum(current_data_base$sink[negative_loc])
      }
    }
  }
  
  p_soc_control_exp =
    ggplot() +
    geom_bar(data = current_data_base, aes(x = control*100, y = sink, fill = as.factor(depth)), color = 'white', width = 9.5, stat = 'identity', position = 'stack', size = 0.5) +
    geom_errorbar(data = current_data_base, aes(x = control*100, y = sink_error_base, ymin = sink_error_base - 2*sink_std, ymax = sink_error_base + 2*sink_std), stat = 'identity', width = 2, size = 1) +
    geom_ribbon(data = current_data_net_control, aes(x = control*100, y = sink, ymin = sink - 2*sink_std, ymax = sink + 2*sink_std), fill = "black", alpha = 0.5) +
    geom_line(data = current_data_net_control, aes(x = control*100, y = sink), size = 1.5) +
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
  
  
  
  ilayer = 2
  for (ilayer in 1:4) {
    if (ilayer == 1) {
      p =
        ggplot() +
        geom_ribbon(data = current_data_base[current_data_base$depth == ilayer, ], aes(x = control*100, y = tau, ymin = tau - 2*tau_std, ymax = tau + 2*tau_std), alpha = 0.5, fill = depth_color_scheme[ilayer]) +
        geom_line(data = current_data_base[current_data_base$depth == ilayer, ], aes(x = control*100, y = tau), color = depth_color_scheme[ilayer], linetype = 'solid', size = 2.5) +
        geom_point(data = current_data_base[current_data_base$depth == ilayer, ], aes(x = control*100, y = tau), color = depth_color_scheme[ilayer], shape = 18, size = 10) +
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
        theme(axis.title.y = element_text(hjust = -0.6)) +
        theme(panel.background = element_rect(fill = 'transparent', color = NA),  # plot panel
              plot.background  = element_rect(fill = 'transparent', color = NA),  # entire plot
              legend.background = element_rect(fill = 'transparent', color = NA), # legend bg
              legend.box.background = element_rect(fill = 'transparent', color = NA)) # legend box
      
    } else {
      p =
        ggplot() +
        geom_ribbon(data = current_data_base[current_data_base$depth == ilayer, ], aes(x = control*100, y = tau, ymin = tau - 2*tau_std, ymax = tau + 2*tau_std), alpha = 0.5, fill = depth_color_scheme[ilayer]) +
        geom_line(data = current_data_base[current_data_base$depth == ilayer, ], aes(x = control*100, y = tau), color = depth_color_scheme[ilayer], linetype = 'solid', size = 2.5) +
        geom_point(data = current_data_base[current_data_base$depth == ilayer, ], aes(x = control*100, y = tau), color = depth_color_scheme[ilayer], shape = 18, size = 10) +
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
              axis.ticks.length.y = unit(0.12, 'inch'), axis.ticks.length.x = unit(0, 'inch')) + 
        theme(panel.background = element_rect(fill = 'transparent', color = NA),  # plot panel
              plot.background  = element_rect(fill = 'transparent', color = NA),  # entire plot
              legend.background = element_rect(fill = 'transparent', color = NA), # legend bg
              legend.box.background = element_rect(fill = 'transparent', color = NA)) # legend box
      
    }
    
    eval(parse(text = paste('p', ilayer, ' = p', sep = '')))
    
  }
  
  p_tau_control_exp =
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
  # fluxes changes with reallocation
  #################################################################################
  current_data_flux = cbind(current_data_flux_raw[ , 1, 1], 
                            apply(current_data_flux_raw[ , 2, ], 1, mean, na.rm = TRUE), # soc_sink
                            apply(current_data_flux_raw[ , 2, ], 1, sd, na.rm = TRUE),
                            apply(current_data_flux_raw[ , 3, ], 1, mean, na.rm = TRUE), # resp_change
                            apply(current_data_flux_raw[ , 3, ], 1, sd, na.rm = TRUE),
                            apply(current_data_flux_raw[ , 4, ], 1, mean, na.rm = TRUE), # input_change
                            apply(current_data_flux_raw[ , 4, ], 1, sd, na.rm = TRUE),
                            apply(current_data_flux_raw[ , 5, ], 1, mean, na.rm = TRUE), # process_change
                            apply(current_data_flux_raw[ , 5, ], 1, sd, na.rm = TRUE),
                            apply(current_data_flux_raw[ , 4, ] - current_data_flux_raw[ , 3, ], 1, mean, na.rm = TRUE), # net flux
                            apply(current_data_flux_raw[ , 4, ] - current_data_flux_raw[ , 3, ], 1, sd, na.rm = TRUE), # net flux
                            current_data_flux_raw[ , 6, 1],
                            current_data_flux_raw[ , 7, 1] 
  )
  
  current_data_flux = data.frame(current_data_flux)
  colnames(current_data_flux) = c('gradient', 'soc_sink', 'soc_sink_std',
                                  'resp_change', 'resp_change_std',
                                  'input_change', 'input_change_std',
                                  'process_change', 'process_change_std',
                                  'net_flux', 'net_flux_std',
                                  'process_name', 'depth')
  
  #-----------------------------------------
  current_data_plot = current_data_flux[current_data_flux$process_name == 2, c('gradient', 'resp_change', 'resp_change_std',
                                                                               'input_change', 'input_change_std',
                                                                               'net_flux', 'net_flux_std', 'depth')]
  
  p_flux =
    ggplot() +
    geom_ribbon(data = current_data_plot, aes(x = gradient*100, y = net_flux, ymin = net_flux - 2*net_flux_std, ymax = net_flux + 2*net_flux_std, fill = as.factor(depth)), alpha = 0.5) +
    
    geom_hline(yintercept = 0, color = 'grey', size = 2, linetype = 'dashed') +
    geom_line(data = current_data_plot, aes(x = gradient*100, y = net_flux, color = as.factor(depth)), size = 2) +
    geom_point(data = current_data_plot, aes(x = gradient*100, y = net_flux, color = as.factor(depth)), shape = 16, size = 9) +
    geom_ribbon(data = current_data_net_control, aes(x = control*100, y = sink/111, ymin = sink/111 - 2*sink_std/111, ymax = sink/111 + 2*sink_std/111), fill = 'black', alpha = 0.5) +
    geom_line(data = current_data_net_control, aes(x = control*100, y = sink/111), color = 'black', size = 2) +
  scale_color_manual(name = '', values = depth_color_scheme, label = soil_depth_list) +
    scale_fill_manual(name = '', values = depth_color_scheme, label = soil_depth_list) +
    scale_shape_manual(name = '', values = c(15, 16), label = c('Soil respiration', 'Carbon input')) +
    scale_y_continuous(n.breaks = 7, limits = c(-0.4, 0.6)) +
    coord_cartesian() +
    theme_classic() +
    # add title
    labs(title = '', x = '', y = expression(paste('Net carbon fluxes (Pg C yr'^'-1', ')'))) +
    # change the legend properties
    guides(color = guide_legend(ncol = 2)) +
    guides(shape = guide_legend(ncol = 2)) +
    theme(legend.text = element_text(size = 25), legend.title = element_text(size = 0), legend.direction = 'vertical', legend.box = "vertical")  +
    theme(legend.justification = c(0.0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
    theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.5, vjust = 0, size = 40)) +
    # modify the font size
    # modify the margin
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
    theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
    theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
    theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch')) +
    theme(panel.background = element_rect(fill = "transparent", color = NA),  # plot panel
          plot.background  = element_rect(fill = "transparent", color = NA),  # entire plot
          legend.background = element_rect(fill = "transparent", color = NA), # legend bg
          legend.box.background = element_rect(fill = "transparent", color = NA)) # legend box
  
  
  #################################################################################
  # process time series changes
  #################################################################################
  ilayer = 1
  for (ilayer in 1:4) {
    if (ilayer == 4) {
      p =
        ggplot() +
        geom_ribbon(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = gradient*100, y = input_change*100, ymin = (input_change - 2*input_change_std)*100, ymax = (input_change + 2*input_change_std)*100), fill = depth_color_scheme[ilayer], alpha = 0.5) +
        geom_ribbon(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = gradient*100, y = resp_change*100, ymin = (resp_change - 2*resp_change_std)*100, ymax = (resp_change + 2*resp_change_std)*100), fill = depth_color_scheme[ilayer], alpha = 0.5) +
        geom_line(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = gradient*100, y = input_change*100), color = depth_color_scheme[ilayer], linetype = 'solid', size = 2.5) +
        geom_line(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = gradient*100, y = resp_change*100), color = depth_color_scheme[ilayer], linetype = 'dashed', size = 2.5) +
        # geom_point(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = year, y = net_flux), color = depth_color_scheme[ilayer], size = 5, stroke = 5) +
        scale_y_continuous(limits = c(NA, NA), n.breaks = 4, trans = 'identity') +
        coord_cartesian(clip = 'off') + 
        theme_classic() +
        annotate('text', label = expression(paste('×10'^'-2', sep = '')), x = 0, y = 2.3, hjust = 0.5, vjust = -0.5,  size = 10) + 
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
        theme(axis.title.y = element_text(hjust = -0.3)) +
        theme(panel.background = element_rect(fill = "transparent", color = NA),  # plot panel
              plot.background  = element_rect(fill = "transparent", color = NA),  # entire plot
              legend.background = element_rect(fill = "transparent", color = NA), # legend bg
              legend.box.background = element_rect(fill = "transparent", color = NA)) # legend box
      
    } else {
      p =
        ggplot() +
        geom_ribbon(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = gradient*100, y = input_change, ymin = (input_change - 2*input_change_std), ymax = (input_change + 2*input_change_std)), fill = depth_color_scheme[ilayer], alpha = 0.5) +
        geom_ribbon(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = gradient*100, y = resp_change, ymin = (resp_change - 2*resp_change_std), ymax = (resp_change + 2*resp_change_std)), fill = depth_color_scheme[ilayer], alpha = 0.5) +
        geom_line(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = gradient*100, y = input_change), color = depth_color_scheme[ilayer], linetype = 'solid', size = 2.5) +
        geom_line(data = current_data_plot[current_data_plot$depth == ilayer, ], aes(x = gradient*100, y = resp_change), color = depth_color_scheme[ilayer], linetype = 'dashed', size = 2.5) +
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
              axis.ticks.length.y = unit(0.12, 'inch'), axis.ticks.length.x = unit(0, 'inch')) +
        theme(panel.background = element_rect(fill = "transparent", color = NA),  # plot panel
              plot.background  = element_rect(fill = "transparent", color = NA),  # entire plot
              legend.background = element_rect(fill = "transparent", color = NA), # legend bg
              legend.box.background = element_rect(fill = "transparent", color = NA)) # legend box
      
    }
    
    eval(parse(text = paste('p', ilayer, ' = p', sep = '')))
    
  }
  
  p_input_change =
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
  
  eval(parse(text = paste('p_input_change_', vert_component, '= p_input_change', sep = '')))
  eval(parse(text = paste('p_flux_', vert_component, '= p_flux', sep = '')))
  eval(parse(text = paste('p_soc_control_exp_', vert_component, '= p_soc_control_exp', sep = '')))
  eval(parse(text = paste('p_tau_control_exp_', vert_component, '= p_tau_control_exp', sep = '')))
  
}

jpeg(paste('./Soil_Sink/proda_global_gradient_sink_flux_1.jpeg', sep = ''), width = 23, height = 20, units = 'in', res = 300)
p_flux_1 = p_flux_1 + labs(x = x_axis_label[1]) + 
  theme(axis.title.x = element_text(hjust = 3.5))

plot_grid(
  NULL, p_soc_control_exp_1, p_input_change_1, NULL, 
  NULL, p_flux_1, p_tau_control_exp_1, NULL, 
  nrow = 2, 
  labels = c(' ', 'a', 'b', ' ',
             ' ', 'c', 'd', ' '), 
  rel_widths = c(0.03, 1, 1, 0.03),
  label_size = 50,
  label_x = 0., label_y = 1.01,
  label_fontfamily = 'Arial',
  label_fontface = 'bold'
)
dev.off()


jpeg(paste('./Soil_Sink/proda_global_gradient_sink_flux_2.jpeg', sep = ''), width = 23, height = 20, units = 'in', res = 300)
p_flux_2 = p_flux_2 + labs(x = x_axis_label[2]) + 
  theme(axis.title.x = element_text(hjust = 3.5))

plot_grid(
  NULL, p_soc_control_exp_2, p_input_change_2, NULL, 
  NULL, p_flux_2, p_tau_control_exp_2, NULL, 
  nrow = 2, 
  labels = c(' ', 'a', 'b', ' ',
             ' ', 'c', 'd', ' '), 
  rel_widths = c(0.03, 1, 1, 0.03),
  label_size = 50,
  label_x = 0., label_y = 1.01,
  label_fontfamily = 'Arial',
  label_fontface = 'bold'
)
dev.off()

# 
# jpeg(paste('./Soil_Sink/proda_global_gradient_sink_flux_1.jpeg', sep = ''), width = 30, height = 10, units = 'in', res = 300)
# p_flux_1 = p_flux_1 + labs(x = x_axis_label[1]) + 
#   theme(axis.title.x = element_text(hjust = 0.5))
# 
# plot_grid(
#   NULL, p_soc_control_exp_1, p_flux_1, p_tau_control_exp_1, NULL, 
#   nrow = 1, 
#   labels = c(' ', 'a', 'b', 'c', ' '), 
#   rel_widths = c(0.03, 1, 1, 1, 0.03),
#   label_size = 50,
#   label_x = 0., label_y = 1.01,
#   label_fontfamily = 'Arial',
#   label_fontface = 'bold'
# )
# dev.off()
# 
# 
# jpeg(paste('./Soil_Sink/proda_global_gradient_sink_flux_2.jpeg', sep = ''), width = 30, height = 10, units = 'in', res = 300)
# p_flux_2 = p_flux_2 + labs(x = x_axis_label[2]) + 
#   theme(axis.title.x = element_text(hjust = 0.5))
# 
# plot_grid(
#   NULL, p_soc_control_exp_2, p_flux_2, p_tau_control_exp_2, NULL, 
#   nrow = 1, 
#   labels = c(' ', 'a', 'b', 'c', ' '), 
#   rel_widths = c(0.03, 1, 1,  1, 0.03),
#   label_size = 50,
#   label_x = 0., label_y = 1.01,
#   label_fontfamily = 'Arial',
#   label_fontface = 'bold'
# )
# dev.off()


