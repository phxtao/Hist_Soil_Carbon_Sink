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

library(lme4)
library(lmerTest)
library(quantreg)
library(ggpubr)

library(sf)
library(sp)
library(proj4)

# dev.off()
##
rm(list = ls())

setwd('/Users/phxtao/Google_Drive/R')

#############################################################################
# Data Path
#############################################################################
model_name = 'cesm2_clm5_cen_vr_v2'
nn_exp_name = 'exp_pc_cesm2_23'
time_domain = 'whole_time'

data_dir_output = '/Users/phxtao/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/phxtao/DATAHUB/ENSEMBLE/INPUT_DATA/'

#############################################################################
# plant allocation product
#############################################################################

npp_allocation_product = readMat(paste(data_dir_input, 'trendy_hist/input_allocation_data_product.mat', sep = ''))
npp_allocation_product = npp_allocation_product$summary.npp.dist

soil_allocation_depth = c(0, 10, 30, 50, 70, 90, 125, 175)

resolution_lon = 0.5
resolution_lat = 0.5

grid_lon_list = seq(from =-180 + resolution_lon/2, to = 180 - resolution_lon/2, by = resolution_lon)
grid_lat_list = seq(from = 90 - resolution_lat/2, to = -90 + resolution_lon/2, by = -resolution_lat)

subsoil_allocation_fraction = apply(npp_allocation_product[ , , 4:8], c(1, 2), sum, na.rm = FALSE)/apply(npp_allocation_product[ , , 2:8], c(1, 2), sum, na.rm = FALSE)*(npp_allocation_product[ , , 1])/100

# note here the unit for npp is Mg ha-1 = 100g m-2
subsoil_npp = apply(npp_allocation_product[ , , 4:8], c(1, 2), sum, na.rm = FALSE)
topsoil_npp = apply(npp_allocation_product[ , , 2:3], c(1, 2), sum, na.rm = FALSE)

abovesoil_npp = (subsoil_npp + topsoil_npp)/(npp_allocation_product[ , , 1]/100)*(1 - npp_allocation_product[ , , 1]/100)
# ANPP 
# summary(as.array(1 - npp_allocation_product[ , , 1]))
#############################################################################
# global peatland
#############################################################################
global_peatland = readMat(paste('/Users/phxtao/DATAHUB/global_peatland/global_peatland_half_deg.mat', sep = ''))
global_peatland = global_peatland$peatland.map.half.deg

#############################################################################
# soc products
#############################################################################
## WISE
WISE_SOC_Original = readMat('/Users/phxtao/Google_Drive/Tsinghua_Luo/Projects/ORCHIDEE/Tao_ORCHIDEE_MICT/WISE_Half_Deg.mat')
WISE_SOC_Original = WISE_SOC_Original$wise.half.deg
colnames(WISE_SOC_Original) = c('Lon', 'Lat', 'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7')

WISE_SOC = cbind(WISE_SOC_Original[ , c('Lon', 'Lat')], 
                 rowSums(cbind(WISE_SOC_Original[ , 'D1'], WISE_SOC_Original[ , 'D2']*(40 - 30)/(40 - 20)), na.rm = TRUE),
                 rowSums(WISE_SOC_Original[ , 3:7], na.rm = TRUE), 
                 rowSums(WISE_SOC_Original[ , 3:9], na.rm = TRUE))
colnames(WISE_SOC) = c('Lon', 'Lat', 'Top_SOC', 'Sub_SOC', 'Bottom_SOC')

## HWSD
HWSD_SOC = readMat('/Users/phxtao/Google_Drive/Tsinghua_Luo/Projects/CLM/Tao_CLM/HWSD_Half_Deg.mat')
HWSD_SOC = cbind(HWSD_SOC$HWSD.Half.Deg, NA)
colnames(HWSD_SOC) = c('Lon', 'Lat', 'Top_SOC', 'Sub_SOC', 'Bottom_SOC')
HWSD_SOC[ , 'Sub_SOC'] = HWSD_SOC[ , 'Sub_SOC'] + HWSD_SOC[ , 'Top_SOC']

## NCSCD
NCSCD_SOC = readMat('/Users/phxtao/Google_Drive/Tsinghua_Luo/Projects/ORCHIDEE/Tao_ORCHIDEE_MICT/ncscd_half_deg.mat')
NCSCD_SOC = NCSCD_SOC$ncscd.half.deg
colnames(NCSCD_SOC) = c('Lon', 'Lat', 'Top_SOC', 'Sub_SOC', 'Bottom_SOC')
NCSCD_SOC[ , 'Bottom_SOC'] = NCSCD_SOC[ , 'Sub_SOC'] + NCSCD_SOC[ , 'Bottom_SOC']

NCSCD_SOC[which(NCSCD_SOC[ , 'Top_SOC'] == 'NaN'), 'Top_SOC'] = WISE_SOC[which(NCSCD_SOC[ , 'Top_SOC'] == 'NaN'), 'Top_SOC']
NCSCD_SOC[which(NCSCD_SOC[ , 'Sub_SOC'] == 'NaN'), 'Sub_SOC'] = WISE_SOC[which(NCSCD_SOC[ , 'Sub_SOC'] == 'NaN'), 'Sub_SOC']
NCSCD_SOC[which(NCSCD_SOC[ , 'Bottom_SOC'] == 'NaN'), 'Bottom_SOC'] = WISE_SOC[which(NCSCD_SOC[ , 'Bottom_SOC'] == 'NaN'), 'Bottom_SOC']

invalid_loc = which(NCSCD_SOC[ , 'Top_SOC'] > NCSCD_SOC[ , 'Sub_SOC'] 
                    | NCSCD_SOC[ , 'Sub_SOC'] > NCSCD_SOC[ , 'Bottom_SOC']
                    | NCSCD_SOC[ , 'Top_SOC'] > NCSCD_SOC[ , 'Bottom_SOC'])
NCSCD_SOC[invalid_loc, c('Top_SOC', 'Sub_SOC', 'Bottom_SOC')] = NA


##############################################
# global envinfo 
##############################################
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


grid_env_info = readMat(paste(data_dir_input, 'data4nn/world_grid_envinfo_present.mat', sep = ''))
grid_env_info = grid_env_info$EnvInfo
colnames(grid_env_info) = grid_var_names

# load climate info
global_climate = grid_env_info[ , 'Koppen_Climate_2018']
var_climate = array(NA, dim = c(length(global_climate), 1))
# koppen climatte class 2
var_climate = global_climate
var_climate[which(global_climate == 1)] = 101 # Af
var_climate[which(global_climate == 2)] = 102 # Am
var_climate[which(global_climate == 3)] = 103 # Aw
var_climate[which(global_climate >= 4 & global_climate <= 5)] = 104 # BwX
var_climate[which(global_climate >= 6 & global_climate <= 7)] = 105 # BsW
var_climate[which(global_climate >= 8 & global_climate <= 10)] = 106 # CsX
var_climate[which(global_climate >= 11 & global_climate <= 13)] = 107 # CwX
var_climate[which(global_climate >= 14 & global_climate <= 16)] = 108 # CfX
var_climate[which(global_climate >= 17 & global_climate <= 20)] = 109 # DsX
var_climate[which(global_climate >= 21 & global_climate <= 24)] = 110 # DwX
var_climate[which(global_climate >= 25 & global_climate <= 28)] = 111 # DfX
var_climate[which(global_climate >= 29 & global_climate <= 30)] = 112 # E

# load soil order
global_soilorder = grid_env_info[ , 'USDA_Suborder']
var_soilorder = array(NA, dim = c(length(global_soilorder), 1))

var_soilorder[] = 113 # others
var_soilorder[which(global_soilorder >= 5 & global_soilorder <= 7)] = 101 # Gelisols
var_soilorder[which(global_soilorder >= 10 & global_soilorder <= 13)] = 102 # Histosols
var_soilorder[which(global_soilorder >= 15 & global_soilorder <= 19)] = 103 # Spodosols
var_soilorder[which(global_soilorder >= 20 & global_soilorder <= 27)] = 104 # Andisols
var_soilorder[which(global_soilorder >= 30 & global_soilorder <= 34)] = 105 # Oxisols
var_soilorder[which(global_soilorder >= 40 & global_soilorder <= 45)] = 106 # Vertisols
var_soilorder[which(global_soilorder >= 50 & global_soilorder <= 56)] = 107 # Aridisols
var_soilorder[which(global_soilorder >= 60 & global_soilorder <= 64)] = 108 # Ultisols
var_soilorder[which(global_soilorder >= 69 & global_soilorder <= 77)] = 109 # Mollisols
var_soilorder[which(global_soilorder >= 80 & global_soilorder <= 84)] = 110 # Alfisols
var_soilorder[which(global_soilorder >= 85 & global_soilorder <= 86)] = 111 # Inceptisols
var_soilorder[which(global_soilorder >= 89 & global_soilorder <= 94)] = 111 # Inceptisols
var_soilorder[which(global_soilorder >= 95 & global_soilorder <= 99)] = 112 # Entisols

# load ESA land cover
global_landcover = grid_env_info[ , 'ESA_Land_Cover']
var_landcover = array(NA, dim = c(length(global_landcover), 1))

var_landcover[which(global_landcover >= 1 & global_landcover <= 2)] = 101 # agriculture
var_landcover[which(global_landcover >= 3 & global_landcover <= 4)] = 102 # mosaic agriculture
var_landcover[which(global_landcover >= 5 & global_landcover <= 6)] = 103 # broadleaved forest 
var_landcover[which(global_landcover >= 7 & global_landcover <= 8)] = 104 # needleleaved forest 
var_landcover[which(global_landcover >= 9 & global_landcover <= 9)] = 105 # mixed forest 
var_landcover[which(global_landcover >= 10 & global_landcover <= 11)] = 106 # Mosaic tree and shrub
var_landcover[which(global_landcover >= 12 & global_landcover <= 12)] = 107 # shrub
var_landcover[which(global_landcover >= 13 & global_landcover <= 13)] = 108 # grassland
var_landcover[which(global_landcover >= 14 & global_landcover <= 15)] = 109 # Lichen & mosses, sparse vegatation
var_landcover[which(global_landcover >= 16 & global_landcover <= 18)] = 110 # wetland
var_landcover[which(global_landcover >=19)] = 111 # urban and other

#############################################################################
# soc products
#############################################################################
current_data = array(NA, dim = c(nrow(WISE_SOC), 8))

igrid = 1
for (igrid in 1:nrow(WISE_SOC)) {
  grid_lon = WISE_SOC[igrid, 'Lon']
  grid_lat = WISE_SOC[igrid, 'Lat']
  
  grid_lon_coord = length(seq(from = -180 + resolution_lon/2, to = grid_lon, by = resolution_lon))
  grid_lat_coord = length(seq(from = 90 - resolution_lat/2, to = grid_lat, by = -resolution_lat ))
  
  current_data[igrid, ] = c(grid_lon, grid_lat, 
                            1 - npp_allocation_product[grid_lat_coord, grid_lon_coord, 1]/100,
                            topsoil_npp[grid_lat_coord, grid_lon_coord],
                            subsoil_npp[grid_lat_coord, grid_lon_coord],
                            abovesoil_npp[grid_lat_coord, grid_lon_coord],
                            subsoil_allocation_fraction[grid_lat_coord, grid_lon_coord],
                            global_peatland[grid_lat_coord, grid_lon_coord]
  )
}

current_data = data.frame(current_data)
colnames(current_data) = c('lon', 'lat', 'f_anpp', 'topsoil_npp', 'subsoil_npp', 'abovesoil_npp', 'subsoil_allocation', 'peatland')

current_data$wise_subsoc_fraction = 1 - WISE_SOC[ , 'Top_SOC']/WISE_SOC[ , 'Sub_SOC']
current_data$hwsd_subsoc_fraction = 1 - HWSD_SOC[ , 'Top_SOC']/HWSD_SOC[ , 'Sub_SOC']
current_data$ncscd_subsoc_fraction = 1 - NCSCD_SOC[ , 'Top_SOC']/NCSCD_SOC[ , 'Sub_SOC']
current_data$wise_subsoc = WISE_SOC[ , 'Sub_SOC'] - WISE_SOC[ , 'Top_SOC']
current_data$hwsd_subsoc = HWSD_SOC[ , 'Sub_SOC'] - HWSD_SOC[ , 'Top_SOC']
current_data$ncscd_subsoc = NCSCD_SOC[ , 'Sub_SOC'] - NCSCD_SOC[ , 'Top_SOC']

current_data$wise_totsoc = WISE_SOC[ , 'Sub_SOC']
current_data$hwsd_totsoc = HWSD_SOC[ , 'Sub_SOC']
current_data$ncscd_totsoc = NCSCD_SOC[ , 'Sub_SOC']

current_data$climate = var_climate
current_data$soil_order = var_soilorder
current_data$land_cover = var_landcover
current_data$mat = grid_env_info[ , 'Annual Mean Temperature']
current_data$map = grid_env_info[ , 'Annual Precipitation']
current_data$clay = grid_env_info[ , 'Clay_Content_30cm']
current_data$npp = grid_env_info[ , 'NPPmean']

# current_data$subsoil_allocation[current_data$subsoil_allocation == 0] = NA
# current_data$subsoil_allocation[current_data$wise_subsoc == 0] = NA
# current_data$subsoil_allocation[current_data$hwsd_subsoc == 0] = NA
# current_data$subsoil_allocation[current_data$ncscd_subsoc == 0] = NA

current_data$subsoc_mean = (current_data$wise_subsoc + current_data$hwsd_subsoc + current_data$ncscd_subsoc)/3
current_data$totsoc_mean = (current_data$wise_totsoc + current_data$hwsd_totsoc + current_data$ncscd_totsoc)/3
current_data$topsoc_mean = current_data$totsoc_mean - current_data$subsoc_mean
current_data$subsoc_fraction_mean = (current_data$wise_subsoc_fraction + current_data$hwsd_subsoc_fraction + current_data$ncscd_subsoc_fraction)/3

current_data$topsoc_fraction_mean = 1-current_data$subsoc_fraction_mean
current_data$topsoil_allocation = current_data$topsoil_npp/(current_data$topsoil_npp + current_data$subsoil_npp + current_data$abovesoil_npp)

##-------------------statistical summary
summary(current_data$f_anpp)
summary(current_data$subsoil_allocation)
summary(current_data$subsoil_allocation/(1-current_data$f_anpp))

#############################################################################
# suballocation vs subsoc fraction 
#############################################################################
cor.test(current_data$subsoc_fraction_mean, current_data$subsoil_allocation)
cor.test(current_data$topsoc_fraction_mean, current_data$topsoil_allocation)

regression_model_data = data.frame(current_data[ , c('subsoil_allocation', 'topsoil_allocation', 'npp', 'topsoil_npp', 'subsoil_npp', 'mat', 'map', 'clay', 'totsoc_mean', 'topsoc_mean', 'subsoc_mean', 'subsoc_fraction_mean', 'topsoc_fraction_mean', 'climate', 'soil_order', 'land_cover')])
regression_model_data = regression_model_data[is.na(apply(regression_model_data, 1, sum, na.rm = FALSE)) == 0, ]

#--------------------------------------
# mixed model on sub
#--------------------------------------
mix_model_sub = lmer(subsoc_fraction_mean ~ subsoil_allocation + subsoil_npp + subsoc_mean + (1 | land_cover) + (1 | climate) + (1 | soil_order), data = regression_model_data)
mix_model_summary_sub = summary(mix_model_sub)
mix_model_summary_sub

r_squared_mixed_sub = summary(lm(regression_model_data$subsoc_fraction_mean ~ predict(mix_model_sub)))
r_squared_mixed_sub

fit_function_meta = function(x) {100*mix_model_summary_sub$coefficients[1, 1] + x*(mix_model_summary_sub$coefficients[2, 1])}

if (mix_model_summary_sub$coefficients[1, 5] < 0.001) {
  p_intercept = paste(' < 0.001', sep = '')
} else {
  p_intercept = paste(' = ', round(mix_model_summary_sub$coefficients[1, 5], 3), sep = '')
}

if (mix_model_summary_sub$coefficients[2, 5] < 0.001) {
  p_slope = paste(' < 0.001', sep = '')
} else {
  p_slope = paste(' = ', round(mix_model_summary_sub$coefficients[2, 5], 3), sep = '')
}

text_data = data.frame('x_axis' = c(0.01), 
                       'y_axis' = c(20), 
                       'equation' = paste('fSOC = ', round(mix_model_summary_sub$coefficients[1, 1], 2), ' + ', round(mix_model_summary_sub$coefficients[2, 1], 2), '*fNPP', '\nP(Intercept)', p_intercept, ', P(fNPP)', p_slope, '\nExplained Variation = ', round(r_squared_mixed_sub$r.squared*100, 2), '%', sep = ''))


scatter_plot_subsoil =
  ggplot(data = current_data[current_data$subsoil_npp > 0, ]) + 
  geom_pointdensity(aes(x = subsoil_allocation*100, y = (subsoc_fraction_mean)*100), shape = 19, size = 1, method='neighbors') +
  scale_color_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 600), oob = scales::squish) +
  scale_x_continuous(limits = c(NA, NA), n.breaks = 7, trans = 'identity') +
  scale_y_continuous(limits = c(NA, NA), n.breaks = 7, trans = 'identity') +
  geom_function(fun = fit_function_meta, size = 2, color = 'black') + 
  geom_label(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 9, show.legend = FALSE, fill = 'white', alpha = 0.6, label.size = NA) + 
  theme_classic() + 
  # add title
  labs(title = '', x = 'Deep soil plant input fraction (%)', y = 'Deep SOC fraction (%)') + 
  # change the legend properties
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'right', title.hjust = 0, title.vjust = 0.8, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.position = 'None', legend.justification = c(0, 1), legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch')) 

scatter_plot_subsoil =
  ggMarginal(scatter_plot_subsoil, type = 'histogram')

#--------------------------------------
# mixed model on top
#--------------------------------------

mix_model_top = lmer(topsoc_fraction_mean ~ topsoil_allocation + topsoil_npp + topsoc_mean + (1 | land_cover) + (1 | climate) + (1 | soil_order), data = regression_model_data)
mix_model_summary_top = summary(mix_model_top)
mix_model_summary_top

r_squared_mixed_top = summary(lm(regression_model_data$topsoc_fraction_mean ~ predict(mix_model_top)))
r_squared_mixed_top

fit_function_meta = function(x) {100*mix_model_summary_top$coefficients[1, 1] + x*(mix_model_summary_top$coefficients[2, 1])}

if (mix_model_summary_top$coefficients[1, 5] < 0.001) {
  p_intercept = paste(' < 0.001', sep = '')
} else {
  p_intercept = paste(' = ', round(mix_model_summary_top$coefficients[1, 5], 3), sep = '')
}

if (mix_model_summary_top$coefficients[2, 5] < 0.001) {
  p_slope = paste(' < 0.001', sep = '')
} else {
  p_slope = paste(' = ', round(mix_model_summary_top$coefficients[2, 5], 3), sep = '')
}

text_data = data.frame('x_axis' = c(8), 
                       'y_axis' = c(100), 
                       'equation' = paste('fSOC = ', round(mix_model_summary_top$coefficients[1, 1], 2), ' + ', round(mix_model_summary_top$coefficients[2, 1], 2), '*fNPP', '\nP(Intercept)', p_intercept, ', P(fNPP)', p_slope, '\nExplained Variation = ', round(r_squared_mixed_top$r.squared*100, 2), '%', sep = ''))


# jpeg(paste('/Users/phxtao/Desktop/fig.jpeg', sep = ''), width = 9, height = 9, units = 'in', res = 300)
# 
#   ggplot(data = current_data) + 
#     stat_bin_hex(aes(x = npp, y = topsoc_mean*1000/0.3/grid_env_info[ , 'Bulk_Density_0cm']), bins = 100) +
#     scale_fill_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 50), oob = scales::squish) +
#       scale_x_continuous(limits = c(NA, NA), n.breaks = 7, trans = 'identity') +
#   scale_y_continuous(limits = c(0.1, 1000), trans = 'log10', labels = trans_format('log10', math_format(10^.x))) + 
#   theme_classic() + 
#   # add title
#   labs(title = '', x = expression(paste('Plant carbon input (g C m'^'-2', ' yr'^'-1', ')', sep = '')), y = expression(paste('SOC (g C kg'^'-1', ')', sep = ''))) + 
#   # change the legend properties
#   guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'right', title.hjust = 0, title.vjust = 0.8, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
#   theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
#   theme(legend.position = 'None', legend.justification = c(0, 1), legend.background = element_rect(fill = NA)) + 
#   # modify the position of title
#   theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
#   # modify the font size
#   # modify the margin
#   theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
#   theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch')) 
# dev.off()


scatter_plot_topsoil =
  ggplot(data = current_data) + 
  geom_pointdensity(aes(x = (topsoil_npp/(abovesoil_npp + subsoil_npp + topsoil_npp)*100), y = (1 - subsoc_fraction_mean)*100), shape = 19, size = 1, method='neighbors') +
  scale_color_gradientn(name = 'Count', colors = viridis(7), trans = 'identity', limits = c(1, 600), oob = scales::squish) +
  scale_x_continuous(limits = c(NA, NA), n.breaks = 7, trans = 'identity') +
  scale_y_continuous(limits = c(NA, NA), n.breaks = 7, trans = 'identity') +
  geom_function(fun = fit_function_meta, size = 2, color = 'black') + 
  geom_label(data = text_data, aes(x = x_axis, y = y_axis, label = equation), vjust = 1, hjust = 0, size = 9, show.legend = FALSE, fill = 'white', alpha = 0.6, label.size = NA) + 
  theme_classic() + 
  # add title
  labs(title = '', x = 'Shallow soil plant input fraction (%)', y = 'Shallow SOC fraction (%)') + 
  # change the legend properties
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 15, barheight = 2.5, title.position = 'right', title.hjust = 0, title.vjust = 0.8, label.hjust = 0.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.position = 'None', legend.justification = c(0, 1), legend.background = element_rect(fill = NA)) + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, size = 50)) + 
  # modify the font size
  # modify the margin
  theme(plot.margin = unit(c(0., 0.2, 0.2, 0.2), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch')) 

scatter_plot_topsoil = 
  ggMarginal(scatter_plot_topsoil, type = 'histogram')





#############################################################################
# global maps 
#############################################################################
color_scheme = c('#b2182b', '#ef8a62', '#fddbc7', '#f7f7f7', '#d1e5f0', '#67a9cf', '#2166ac')

world_coastline = st_read('/Users/phxtao/Google_Drive/Tsinghua_Luo/World_Vector_Shape/ne110m/ne_110m_land.shp', layer = 'ne_110m_land')
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


#-----------------current data projection
lon_lat_transfer = project(xy = as.matrix(current_data[ , c('lon', 'lat')]), proj = '+proj=robin') 
current_data[ , c('lon', 'lat')] = lon_lat_transfer

lat_limits = rbind(c(0, -56), c(0, 80))
lat_limits_robin = project(xy = as.matrix(lat_limits), proj = '+proj=robin') 

p_soc_fraction =
  ggplot() +
  geom_tile(data = current_data, aes(x = lon, y = lat, fill = subsoc_fraction_mean*100), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradientn(expression(paste('Deep soil SOC fraction (%)', sep = '')), colours = viridis(7), na.value="transparent", limits = c(30, 60), trans = 'identity', oob = scales::squish) +
  # geom_point(data = current_data[current_data$sig_index == 0, ], aes(x = lon, y = lat), shape = 4, size = 0.3, stroke = 0.1, color = 'black') +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
  # change the background to black and white
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  theme(legend.justification = c(0, 0), legend.position = c(0.1, -0.2), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # theme(legend.justification = c(0, 0), legend.position = 'None', legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 30, barheight = 2.5, title.position = 'right', title.hjust = 1, title.vjust = 0.85, label.hjust = 0.5, label.vjust = 2.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 30, ), legend.title = element_text(size = 35)) +
  # add title
  labs(title = '', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))


p_allocation =
  ggplot() +
  geom_tile(data = current_data, aes(x = lon, y = lat, fill = subsoil_allocation*100), height = 60000, width = 60000, na.rm = TRUE) +
  scale_fill_gradientn(expression(paste('Deep soil plant input fraction (%)', sep = '')), colours = viridis(7), na.value="transparent", limits = c(0, 20), trans = 'identity', oob = scales::squish) +
  # geom_point(data = current_data[current_data$sig_index == 0, ], aes(x = lon, y = lat), shape = 4, size = 0.3, stroke = 0.1, color = 'black') +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
  # change the background to black and white
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  theme(legend.justification = c(0, 0), legend.position = c(0.1, -0.2), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # theme(legend.justification = c(0, 0), legend.position = 'None', legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 30, barheight = 2.5, title.position = 'right', title.hjust = 1, title.vjust = 0.85, label.hjust = 0.5, label.vjust = 2.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 30, ), legend.title = element_text(size = 35)) +
  # add title
  labs(title = '', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))


jpeg(paste('./Soil_Sink/subsoil_soc_fraction_vs_allocation.jpeg', sep = ''), width = 18, height = 17, units = 'in', res = 300)
plot_grid(NULL, p_allocation, NULL,  
          NULL, p_soc_fraction, NULL, 
          nrow = 2,
          ncol = 3,
          rel_widths = c(0.1, 1, 0.1),
          labels = c(' ', 'a', ' ', 
                     ' ', 'b', ' '), 
          label_size = 50,
          label_x = 0.01, label_y = 0.92,
          label_fontfamily = 'Arial',
          label_fontface = 'bold'
)
dev.off()

#############################################################################
# deep rooting
#############################################################################
data_dir_output = '/Users/phxtao/DATAHUB/ENSEMBLE/OUTPUT_DATA/'
data_dir_input = '/Users/phxtao/DATAHUB/ENSEMBLE/INPUT_DATA/'

deep_root_data = read.csv(paste('/Users/phxtao/DATAHUB/SOIL_CARBON_SINK/metadata/PEC-data-202408.csv', sep = ''), header = TRUE)

deep_root_data$depth_node = (deep_root_data$Upper_Layer + deep_root_data$Lower_Layer)/2
deep_root_data$rr_soc = log((deep_root_data$SOC_t)/(deep_root_data$SOC_c))


deep_root_data$depth_type[which(deep_root_data$depth_node <= 30)] = 1
deep_root_data$depth_type[which(deep_root_data$depth_node > 30 & deep_root_data$depth_node <= 100)] = 2
deep_root_data$depth_type[which(deep_root_data$depth_node > 100)] = 3

deep_root_data$depth_type = as.character(deep_root_data$depth_type)
# Perform one-way ANOVA
anova_result = aov(rr_soc ~ depth_type, data = deep_root_data)

# Summary of the ANOVA test
summary(anova_result)

tukey_result = TukeyHSD(anova_result)

print(tukey_result)

# Perform one-sample t-test (e.g., test if mean is different from 0)
t_test_result = c(t.test(deep_root_data[deep_root_data$depth_type == 1, 'rr_soc'], mu = 0)$p.value,
                  t.test(deep_root_data[deep_root_data$depth_type == 2, 'rr_soc'], mu = 0)$p.value,
                  t.test(deep_root_data[deep_root_data$depth_type == 3, 'rr_soc'], mu = 0)$p.value
)

data_count = rbind(c(length(which(deep_root_data$depth_type == 1)), 1),
                   c(length(which(deep_root_data$depth_type == 2)), 2),
                   c(length(which(deep_root_data$depth_type == 3)), 3)
)
data_count = data.frame(data_count)
colnames(data_count) = c('count', 'type')
data_count$type = as.character(data_count$type)
#-------------------------------------------------------
# figures
#-------------------------------------------------------
depth_label = c('0-30cm', '30-100cm', '>100cm')

# jpeg(paste('./Soil_Sink/meta_deep_root.jpeg', sep = ''), width = 8, height = 8, units = 'in', res = 300)
p_root =
  ggplot(data = deep_root_data, aes(x = depth_type, y = rr_soc)) + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1, color = 'grey') + 
  geom_jitter(width = 0.1, alpha = 0.5, size = 2.5, shape = 16, color = '#005AB5') + 
  geom_violin(color = 'black', width = 0.8, fill = NA, size = 2) + 
  geom_boxplot(outlier.shape = NA, width = 0.3, size = 2., alpha = 1, color = 'black', fill = NA) + 
  geom_text(data = data_count, aes(x = type, y = -1, label = count), size = 10) + 
  stat_compare_means(comparisons = list(c(1, 2), c(2, 3), c(1, 3)), tip.length = 0.02, bracket.size = 1.5, size = 10,
                     method = 'wilcox.test', label = 'p.format') +
  scale_x_discrete(labels = depth_label) + 
  theme_classic() +
  # add title
  labs(title = '', x = 'Soil depths', y = expression(paste('log(Response Ratio)'['SOC'], sep = ''))) +
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
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch'))

# dev.off()

# ggplot() + 
# geom_boxplot(data = deep_root_data, aes(x = as.character(depth_type), y = rr_soc, fill = as.character(depth_type)))



#############################################################################
# deep injection
#############################################################################
deep_inject_data = read.csv(paste('/Users/phxtao/DATAHUB/SOIL_CARBON_SINK/metadata/deep_injection_SOC_stock.csv', sep = ''), header = TRUE)
colnames(deep_inject_data) = c('site', 'year', 'replicates', 'layer', 'normal', 'deep_injection')

current_data = rbind(cbind(deep_inject_data[ , c('site', 'year', 'replicates', 'layer')], 'soc' = as.array(deep_inject_data$normal), 'tech' = 'A_shallow'),
                     cbind(deep_inject_data[ , c('site', 'year', 'replicates', 'layer')], 'soc' = as.array(deep_inject_data$deep_injection), 'tech' = 'B_deep')
)


anova_result = aov(soc ~ tech, data = current_data[current_data$layer == '0-20', ])
summary(anova_result)
anova_result = aov(soc ~ tech, data = current_data[current_data$layer == '20-40', ])
summary(anova_result)

mean(current_data[current_data$layer == '20-40' & current_data$tech == 'shallow', 'soc'])
mean(current_data[current_data$layer == '20-40' & current_data$tech == 'deep', 'soc'])

p_top =
  ggplot(data = current_data[current_data$layer == '0-20', ], aes(x = as.factor(year), y = soc/10, color = tech)) + 
  geom_boxplot(position = position_dodge(width = 0.7), outlier.shape = NA, width = 0.6, size = 1.5, alpha = 1,  fill = NA) +
  geom_point(position=position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1), alpha = 0.8, size = 2.5, shape = 16) +
  scale_color_manual(name = '', label = c('Normal Management', 'Deep injection'), values = c('#D55E00', '#009E73')) + 
  coord_cartesian() + 
  theme_classic() +
  # add title
  labs(title = '0 - 20cm', x = '', y = '') +
  # change the legend properties
  theme(plot.title = element_text(hjust = 0.1, vjust = -3, size = 40)) +
  
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'vertical', legend.box = "horizontal")  +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(1, 0), legend.position = c(1, 0), legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text.y = element_text(size = 30, color = 'black'), axis.text.x = element_text(size = 0, color = NULL),
        axis.title.y = element_text(size = 35), axis.title.x = element_text(size = 0), 
        axis.line.y = element_line(size = 1), axis.line.x = element_line(size = 0), 
        axis.ticks.y = element_line(size = 1, color = 'black'), axis.ticks.x = element_line(size = 1, color = NULL), 
        axis.ticks.length.y = unit(0.12, 'inch'), axis.ticks.length.x = unit(0, 'inch'))


p_deep =
  ggplot(data = current_data[current_data$layer == '20-40', ], aes(x = as.factor(year), y = soc/10, color = tech)) + 
  geom_boxplot(position = position_dodge(width = 0.7), outlier.shape = NA, width = 0.6, size = 1.5, alpha = 1,  fill = NA) +
  geom_point(position=position_jitterdodge(dodge.width = 0.9, jitter.width = 0.1), alpha = 0.8, size = 2.5, shape = 16) +
  scale_color_manual(values = c('#D55E00', '#009E73')) + 
  coord_cartesian(clip = 'off') + 
  theme_classic() +
  # add title
  labs(title = '20 - 40cm', x = 'Year', y = expression(paste('SOC stock (kg m'^'-2', ')'))) +
  # change the legend properties
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25), legend.direction = 'vertical', legend.box = "horizontal")  +
  theme(legend.text = element_text(size = 25), legend.title = element_text(size = 25))  +
  theme(legend.justification = c(0, 0), legend.position = 'None', legend.background = element_rect(fill = NA)) +
  theme(legend.key = element_rect(color = NA, fill = NA), legend.key.size = unit(0.5, 'inch')) +
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.1, vjust = -3, size = 40)) +
  # modify the font size
  # modify the margin
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 30, color = 'black'), axis.title = element_text(size = 35), axis.line = element_line(size = 1), axis.ticks = element_line(size = 1, color = 'black'), axis.ticks.length = unit(0.12, 'inch')) +
  theme(axis.title.y = element_text(hjust = 30))


p_injection =
  plot_grid(
    p_top, NULL,
    p_deep, NULL,
    ncol = 2,
    nrow = 2, 
    rel_heights = c(1, 1.4),
    rel_widths = c(1, 0.1),
    greedy = FALSE,
    align = 'v'
  )



jpeg(paste('./Soil_Sink/global_allocation_soc_relationship.jpeg', sep = ''), width = 20, height = 10, units = 'in', res = 300)

plot_grid(NULL, scatter_plot_topsoil, scatter_plot_subsoil, NULL, 
          nrow = 1,
          ncol = 4,
          rel_widths = c(0.1, 1, 1, 0.1),
          labels = c(' ', 'a', 'b', ' '), 
          label_size = 50,
          label_x = -0.02, label_y = 0.98,
          label_fontfamily = 'Arial',
          label_fontface = 'bold'
)

dev.off()


jpeg(paste('./Soil_Sink/deep_rooting_and_injection.jpeg', sep = ''), width = 20, height = 10, units = 'in', res = 300)

plot_grid(NULL, p_root, p_injection, NULL, 
          nrow = 1,
          ncol = 4,
          rel_widths = c(0.1, 1, 1, 0.1),
          labels = c(' ', 'a', 'b', ' '), 
          label_size = 50,
          label_x = -0.02, label_y = 0.98,
          label_fontfamily = 'Arial',
          label_fontface = 'bold'
)

dev.off()

jpeg(paste('./Soil_Sink/meta_and_global_synthesis.jpeg', sep = ''), width = 20, height = 20, units = 'in', res = 300)
plot_grid(NULL, scatter_plot_topsoil, scatter_plot_subsoil, NULL, 
          NULL, p_root, p_injection, NULL, 
          nrow = 2,
          ncol = 4,
          rel_widths = c(0.1, 1, 1, 0.1),
          labels = c(' ', 'a', 'b', ' ', 
                     ' ', 'c', 'd', ' '), 
          label_size = 50,
          label_x = 0.01, label_y = 0.92,
          label_fontfamily = 'Arial',
          label_fontface = 'bold'
)
dev.off()


#############################################################################
# world distribution
#############################################################################
loc_identifier = deep_root_data$Lat**2 + deep_root_data$Lon**2
unique_identifier = unique(loc_identifier)

data_dist_energy = c()

iloc = 1
for (iloc in 1:length(unique_identifier)) {
  unique_loc = which(loc_identifier == unique_identifier[iloc])
  data_dist_energy = rbind(data_dist_energy, 
                           cbind(deep_root_data[unique_loc[1], c('Lon', 'Lat')], length(unique_loc), 1)
  )
}
colnames(data_dist_energy) = c('lon', 'lat', 'number', 'case')

data_dist_injection = rbind(cbind(116.73, 38.03, 24, 2),
                            cbind(116.62, 37.72, 28, 2),
                            cbind(116.33, 37.60, 24, 2),
                            cbind(117.08, 37.10, 18, 2),
                            cbind(123.45, 41.82, 6, 2)
)
colnames(data_dist_injection) = c('lon', 'lat', 'number', 'case')

current_data = rbind(data_dist_energy, data_dist_injection)
current_data = data.frame(current_data)

lon_lat_transfer = project(xy = as.matrix(current_data[ , c('lon', 'lat')]), proj = '+proj=robin') 
current_data[ , c('lon', 'lat')] = lon_lat_transfer

lat_limits = rbind(c(0, -56), c(0, 80))
lat_limits_robin = project(xy = as.matrix(lat_limits), proj = '+proj=robin') 

jpeg(paste('./Soil_Sink/meta_loc.jpeg', sep = ''), width = 20, height = 10, units = 'in', res = 300)
ggplot() +
  geom_point(data = current_data, aes(x = lon, y = lat, size = number, color = as.character(case)), stroke = 2, shape = 1, height = 60000, width = 60000, na.rm = TRUE) +
  geom_sf(data = world_coastline, fill = NA, color = 'black', linewidth = 1) + 
  geom_polygon(data = world_ocean, aes(x = long, y = lat), fill = NA, color = 'black', size = 2) +
  scale_color_manual(name = 'Methods', labels = c('Deep rooting', 'Deep injection'), values = c('#994F00', '#006CD1')) +
  scale_size_continuous(name = 'Sample size', range = c(3, 20)) +
  # theme_map() +
  ylim(lat_limits_robin[ , 2]) +
  # change the legend properties
  theme(legend.justification = c(0.15, 0.1), legend.position = c(0.15, 0.1), legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # theme(legend.justification = c(0, 0), legend.position = 'None', legend.background = element_rect(fill = NA), legend.text.align = 0) +
  # change the size of colorbar
  guides(fill = guide_colorbar(direction = 'horizontal', barwidth = 30, barheight = 2.5, title.position = 'right', title.hjust = 1, title.vjust = 0.85, label.hjust = 0.5, label.vjust = 2.5, frame.linewidth = 0), reverse = FALSE) +
  theme(legend.text = element_text(size = 30, ), legend.title = element_text(size = 35)) +
  # add title
  labs(title = '', x = '', y = '') + 
  # modify the position of title
  theme(plot.title = element_text(hjust = 0.5, vjust = -1, size = 50)) + 
  # modify the font size
  theme(axis.title = element_text(size = 20)) + 
  theme(panel.background = element_rect(fill = NA, colour = NA)) +
  # modify the margin
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), 'inch')) +
  theme(axis.text=element_text(size = 35, color = 'black'))

dev.off()
