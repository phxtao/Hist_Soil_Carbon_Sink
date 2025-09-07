## Packages
library(R.matlab)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(viridis)
library(raster)
library(scales)


library(tidyverse)
library(ggExtra)
library(ggpointdensity)


library(sf)
library(sp)
library(proj4)


# dev.off()
##
rm(list = ls())

setwd('/Users/ft254/Google_Drive/R')



#############################################################################
# site distribution
#############################################################################
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


current_data = read.table('/Users/ft254/DATAHUB/ENSEMBLE/OUTPUT_DATA/nn_soil_carbon_sink/valid_sites.csv', sep = ',')
colnames(current_data) = c('lon', 'lat', 'type')

#-----------------current data projection
lon_lat_transfer = project(xy = as.matrix(current_data[ , c('lon', 'lat')]), proj = '+proj=robin') 
current_data[ , c('lon', 'lat')] = lon_lat_transfer

lat_limits = rbind(c(0, -56), c(0, 80))
lat_limits_robin = project(xy = as.matrix(lat_limits), proj = '+proj=robin') 


jpeg(paste('./Soil_Sink/proda_loc.jpeg', sep = ''), width = 20, height = 10, units = 'in', res = 300)
ggplot() +
  geom_point(data = current_data, aes(x = lon, y = lat, size = as.character(type), color = as.character(type)), shape = 19, alpha = 1, na.rm = TRUE) +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
  scale_color_manual(name = '', labels = c('SOC', '14C'), values = c('#994F00', '#006CD1')) +
  scale_size_manual(name = 'Sample size', values = c(0.5, 3), guide = 'none') +
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # theme(legend.position = 'none') +
  theme(legend.justification = c(0.1, 0.3), legend.position = c(0.1, 0.3), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  theme(legend.text = element_text(size = 35), legend.title = element_text(size = 35))  +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  theme(legend.text = element_text(size = 35), legend.title = element_text(size = 35)) +
  # add title
  labs(title = '', x = '', y = '') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 40)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  # modify the margin
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'inch')) +
  theme(axis.text=element_text(size = 30))

dev.off()

#############################################################################
# agreement
#############################################################################
model_name = 'cesm2_clm5_cen_vr_v2'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

data_dir_output = '/Users/ft254/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/ft254/DATAHUB/ENSEMBLE/INPUT_DATA/'


current_data_12c = readMat('/Users/ft254/DATAHUB/ENSEMBLE/OUTPUT_DATA/hist_simu_proda/12c_summary.mat')
current_data_12c = current_data_12c$prediction.summary
colnames(current_data_12c) = c('proda', 'obs')
current_data_12c = data.frame(current_data_12c)

pred_nse_12c = 1 - sum((current_data_12c$proda - current_data_12c$obs)**2)/sum((current_data_12c$obs - mean(current_data_12c$obs))**2)

p_12c=
  ggplot(data = current_data_12c) + 
  stat_bin_hex(aes(x = obs/1000, y = proda/1000), bins = 100) +
  geom_abline(slope = 1, intercept = 0, size = 2, color = 'black') +
  # geom_pointdensity(aes(x = obs, y = proda), shape = 19, size = 1) +
  scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 500), oob = scales::squish) +
  scale_x_continuous(limits = c(0.1, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
  scale_y_continuous(limits = c(0.1, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
  geom_text(aes(x = 0.1, y = 1000, label = paste('Explained variation = ', round(pred_nse_12c*100, 1), '%', sep = '')), vjust = 1, hjust = 0, size = 9, show.legend = FALSE) + 
  theme_classic() + 
  # add title
  labs(title = '', x = expression(paste('Observed SOC (kg C m'^'-3', ')', sep = '')), y = expression(paste('Modelled SOC (kg C m'^'-3', ')', sep = ''))) + 
  # change the legend properties
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'right', title.hjust = 0, title.vjust = 0.8, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.position = 'None', legend.justification = c(0, 1), legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch')) 


current_data_14c = readMat('/Users/ft254/DATAHUB/ENSEMBLE/OUTPUT_DATA/hist_simu_proda/14c_summary.mat')
current_data_14c = current_data_14c$prediction.summary
colnames(current_data_14c) = c('proda', 'obs')
current_data_14c = data.frame(current_data_14c)
pred_nse_14c = 1 - sum((current_data_14c$proda - current_data_14c$obs)**2)/sum((current_data_14c$obs - mean(current_data_14c$obs))**2)

p_14c=
  ggplot(data = current_data_14c) +
  geom_point(aes(x = obs, y = proda), color = '#005AB5', shape = 19, size = 2) +
  geom_abline(slope = 1, intercept = 0, size = 2, color = 'black') +
  scale_x_continuous(limits = c(-1000, 250), n.breaks = 7, trans = 'identity') +
  scale_y_continuous(limits = c(-1000, 250), n.breaks = 7, trans = 'identity') +
  geom_text(aes(x = -1000, y = 250, label = paste('Explained variation = ', round(pred_nse_14c*100, 1), '%', sep = '')), vjust = 1, hjust = 0, size = 9, show.legend = FALSE) + 
  
  theme_classic() +
  # add title
  labs(title = '', x = expression(paste('Observed ', Delta, '14C', sep = '')), y = expression(paste('Modelled ', Delta, '14C', sep = ''))) +
  # change the legend properties
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'right', title.hjust = 0, title.vjust = 0.8, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.position = 'None', legend.justification = c(0, 1), legend.background = element_rect(fill = NA)) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) +
  # modify the font size
  # modify the margin
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))



jpeg(paste('./Soil_Sink/simu_agreement.jpeg', sep = ''), width = 16, height = 8, units = 'in', res = 300)

plot_grid(
  NULL, p_12c, p_14c, NULL,
  nrow = 1, 
  labels = c(' ', 'a', 'b', ' '), 
  rel_widths = c(0.03, 1, 1, 0.03),
  label_size = 40,
  label_x = 0., label_y = 1.01,
  label_fontfamily = 'Arial',
  label_fontface = 'bold',
  align = 'None'
)
dev.off()

