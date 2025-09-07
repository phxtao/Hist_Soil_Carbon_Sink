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


#-----------------scenario simulation
soc_pool_hist_raw = readMat(paste( data_dir_output, 'hist_simu_proda/process_control/cesm2_hist_soc_proda_default_cross_valid_', as.character(cross_valid_num), '_response_curve_tot.mat', sep = ''))
soc_pool_hist_raw = soc_pool_hist_raw$global.soc.sink.tot

response_curve_summary = c()

idepth = 1
for (idepth in 1:4) {
  icomponent = 1
  for (icomponent in 1:2) {
    if (icomponent == 1) {
      response_ratio = seq(0, 1, 0.1)
    } else {
      response_ratio = seq(0, 1, 0.1)
    }
    
    if (idepth == 1) {
      cpool_middle = soc_pool_hist_raw[ , seq(idepth+16, 28, 4), , icomponent]
    } else {
      cpool_middle = soc_pool_hist_raw[ , seq(idepth+16, 28, 4), , icomponent] - soc_pool_hist_raw[ , seq((idepth - 1)+16, 28, 4), , icomponent]
    }
    
    response_curve_middle = apply(cpool_middle, c(1, 3), sum, na.rm = TRUE)
    
    if (idepth == 4) {
      cpool_net = soc_pool_hist_raw[ , seq(idepth+16, 28, 4), , icomponent]
      response_curve_net = apply(cpool_net, c(1, 3), sum, na.rm = TRUE)
    }
    
    iresponse = 1
    for (iresponse in 1:ncol(response_curve_middle)) {
      response_curve_summary = rbind(response_curve_summary, 
                                     cbind(year_series, 
                                           response_curve_middle[ , iresponse] - response_curve_middle[1, iresponse],
                                           response_ratio[iresponse],
                                           icomponent,
                                           idepth))
    }
    
    if (idepth == 4) {
      for (iresponse in 1:ncol(response_curve_middle)) {
        response_curve_summary = rbind(response_curve_summary, 
                                       cbind(year_series, 
                                             response_curve_net[ , iresponse] - response_curve_net[1, iresponse],
                                             response_ratio[iresponse],
                                             icomponent,
                                             idepth+1))
      }
    }
  }
}

response_curve_summary

colnames(response_curve_summary) = c('year', 'sink', 'response_ratio', 'component', 'depth')
response_curve_summary = data.frame(response_curve_summary)

#################################################################################
# global time series changes
#################################################################################
color_scheme = c('#009E73', '#004D40', '#1E88E5', '#08519c')

current_data = response_curve_summary[response_curve_summary$component == 1 & response_curve_summary$year == 2010 & response_curve_summary$depth < 5, ]
current_data_net = response_curve_summary[response_curve_summary$component == 1 & response_curve_summary$year == 2010 & response_curve_summary$depth == 5, ]
p_npp =
  ggplot() +
  geom_bar(data = current_data, aes(x = response_ratio*100, y = sink, fill = as.factor(depth)), color = 'white', size = 0.5, stat = 'identity', position = 'stack', width = 10) +
  geom_line(data = current_data_net, aes(x = response_ratio*100, y = sink), color = 'darkred', size = 3) +
  geom_hline(yintercept = 0, color = 'black', size = 1) +
  scale_fill_manual(name = c(''), labels = c('0 - 30cm', '30 - 100cm', '100 - 200cm', '> 200cm'), values = color_scheme) +
  scale_x_continuous(n.breaks = 7) +
  scale_y_continuous(limits = c(-40, 100), n.breaks = 7) +
  theme_classic() +
  # add title
  labs(title = '', x = 'Reallocated aboveground input (%)', y = expression(paste('Stock change (Pg C)', sep = ''))) +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'horizontal')  +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 1), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.02, vjust = -7, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))



current_data = response_curve_summary[response_curve_summary$component == 2 & response_curve_summary$year == 2010 & response_curve_summary$depth < 5, ]
current_data_net = response_curve_summary[response_curve_summary$component == 2 & response_curve_summary$year == 2010 & response_curve_summary$depth == 5, ]

p_vtransport =
  ggplot() +
  geom_bar(data = current_data, aes(x = response_ratio*100, y = sink, fill = as.factor(depth)), color = 'white', size = 0.5, stat = 'identity', position = 'stack', width = 10) +
  geom_line(data = current_data_net, aes(x = response_ratio*100, y = sink), color = 'darkred', size = 3) +
  geom_hline(yintercept = 0, color = 'black', size = 1) +
  scale_fill_manual(name = c(''), labels = c('0 - 30cm', '30 - 100cm', '100 - 200cm', '> 200cm'), values = color_scheme) +
  scale_x_continuous(n.breaks = 7) +
  scale_y_continuous(limits = c(-40, 100), n.breaks = 7) +
  theme_classic() +
  # add title
  labs(title = '', x = 'Enhanced vertical transport (%)', y = '') +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'vertical')  +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0.1, 1), legend.position = c(0.1, 1), legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.02, vjust = -7, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30), axis.title = element_text(size = 30), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.12, 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))


jpeg(paste('./Soil_Sink/proda_global_seq_response_curve.jpeg', sep = ''), width = 20, height = 10, units = 'in', res = 300)
plot_grid(NULL, p_npp, p_vtransport, NULL,
          nrow = 1, 
          labels = c(' ', 'a', 'b', ' '), 
          rel_widths = c(0.02, 1, 1, 0.02),
          label_size = 50,
          label_x = 0., label_y = 1.01,
          label_fontfamily = 'Arial',
          label_fontface = 'bold'
)
dev.off()
