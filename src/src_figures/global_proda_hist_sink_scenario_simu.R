## Packages
library(R.matlab)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(viridis)
library(raster)
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
# Load Projected SOC PRODA
#################################################################################
depth_num = 4

cross_valid_num = 7 
cpool_num = 7
year_start = 1900
year_end = 2010
year_series = seq(from = year_start, to = year_end, by = 10)
year_num = length(year_series)
scenario_name = c('Static forcing', 'Input change only', 'Climate change only', 'Input + climate change')


#-----------------scenario simulation
soc_pool_hist_raw = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_scenario_simu_tot.mat', sep = ''))
soc_pool_hist_raw = soc_pool_hist_raw$global.soc.sink.tot

soil_sink_summary = c()

iscenario = 1
for (iscenario in 1:length(scenario_name)) {
  ipool = 5
  for (ipool in 5:7){
    cpool_middle_accu = apply(soc_pool_hist_raw[ , ((ipool-1)*4+1):((ipool-1)*4+4), iscenario], c(1), sum, na.rm = TRUE)
    cpool_middle_accu = cpool_middle_accu - cpool_middle_accu[1]
    cpool_middle_rate = c(0, diff(cpool_middle_accu))/10
    soil_sink_summary = rbind(soil_sink_summary, 
                              cbind(year_series, cpool_middle_accu, cpool_middle_rate, ipool, iscenario))
  }
}


colnames(soil_sink_summary) = c('year', 'sink', 'rate', 'pool', 'scenario')
soil_sink_summary = data.frame(soil_sink_summary)


#################################################################################
# global time series changes
#################################################################################
color_scheme = c('#E69F00', '#009E73', '#0072B2')
depth_list = c('0 - 30cm', '30 - 100cm', '100 - 200cm', '200cm - bedrock')

iscenario = 4
for (iscenario in 1:4) {
  valid_loc = which(soil_sink_summary[ , 'scenario'] == iscenario)
  current_data = data.frame(soil_sink_summary[valid_loc, ])
  
  p_sink =
    ggplot() +
    geom_bar(data = current_data, aes(x = year, y = sink, fill = as.factor(pool)), color = 'white', size = 1, width = 9.5, stat = 'identity', position = 'stack', size = 2) +
    geom_hline(yintercept = 0, color = 'black', size = 1) +
    scale_fill_manual(name = c(''), labels = c('Labile', 'Moderate', 'Persistent'), values = color_scheme) +
    # scale_x_continuous(n.breaks = 7, sec.axis = sec_axis(trans = ~.))
    scale_y_continuous(limits = c(-25, 45), n.breaks = 7) +
    theme_classic() +
    # add title
    labs(title = scenario_name[iscenario], x = '', y = expression(paste('Stock change (Pg C)', sep = ''))) +
    # change the legend properties
    theme(legend.text = element_text(size = 35), legend.title = element_text(size = 35), legend.direction = 'horizontal')  +
    theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
    theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(1, 'inch')) +
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.02, vjust = 0, size = 40)) +
    # modify the font size
    # modify the margin
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
    theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
    theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
    theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))
  
  
  p_rate =
    ggplot() +
    geom_bar(data = current_data, aes(x = year, y = rate, fill = as.factor(pool)), color = 'white', size = 1, stat = 'identity', position = 'stack', width = 9.5) + 
    geom_hline(yintercept = 0, color = 'black', size = 1) +
    scale_fill_manual(name = c(''), labels = c('Labile', 'Moderate', 'Persistent'), values = color_scheme) +
    # scale_x_continuous(n.breaks = 7, sec.axis = sec_axis(trans = ~.)) +
    scale_y_continuous(limits = c(-0.5, 1.0), n.breaks = 7) +
    theme_classic() +
    # add title
    labs(title = scenario_name[iscenario], x = '', y = expression(paste('Stock change rate (Pg C yr'^'-1', ')', sep = ''))) +
    # change the legend properties
    theme(legend.text = element_text(size = 35), legend.title = element_text(size = 35), legend.direction = 'horizontal')  +
    theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
    theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(1, 'inch')) +
    # modify the position of title
    theme(plot.title = element_text(hjust = 0.02, vjust = 0, size = 40)) +
    # modify the font size
    # modify the margin
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
    theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
    theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
    theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))
  
  
  eval(parse(text = paste('p_sink_', iscenario, ' = p_sink', sep = '')))
  eval(parse(text = paste('p_rate_', iscenario, ' = p_rate', sep = '')))
  
}

########################## time change
jpeg(paste('./Soil_Sink/proda_global_seq_time_seris_static_forcing.jpeg', sep = ''), width = 16, height = 16, units = 'in', res = 300)
plot_grid(
  NULL, p_sink_1, NULL,
  NULL, p_rate_1, NULL, 
  nrow = 2, 
  labels = c(' ', 'a', ' ', 
             ' ', 'b', ' '), 
  rel_widths = c(0.02, 1, 0.02),
  label_size = 50,
  label_x = 0., label_y = 1.01,
  label_fontfamily = 'Arial',
  label_fontface = 'bold'
)
dev.off()


p_sink_4 = p_sink_4 +
  theme(legend.justification = c(0, 1), legend.position = c(0, 0.95), legend.background = element_rect(fill = NA)) 

jpeg(paste('./Soil_Sink/proda_global_seq_time_seris.jpeg', sep = ''), width = 30, height = 25, units = 'in', res = 300)
plot_grid(
  NULL, p_sink_4, NULL, p_rate_4, NULL,
  NULL, p_sink_2, NULL, p_rate_2, NULL,
  NULL, p_sink_3, NULL, p_rate_3, NULL, 
  nrow = 3, 
  labels = c(' ', 'a', ' ', 'b', ' ',
             ' ', 'c', ' ', 'd', ' ', 
             ' ', 'e', ' ', 'f', ' '), 
  rel_widths = c(0.02, 1, 0.02, 1, 0.02),
  label_size = 50,
  label_x = 0., label_y = 1.01,
  label_fontfamily = 'Arial',
  label_fontface = 'bold'
)
dev.off()

