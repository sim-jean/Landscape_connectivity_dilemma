################################################################################
########## Generate large scale landscapes #####################################
################################################################################
rm(list = ls())
library(dplyr)
library(ggplot2)
library(magrittr)
library(GA)
library(landscapeR)
library(here)
library(raster)
library(terra)
library(landscapetools)
library(NLMR)

set.seed(123)

size = 100
budget = round(size^2/5) + 1

green_palette = c("Age 0 to 5" = 'cornsilk',
                  "Age 5 to 10" = 'lightgreen',
                  "Age > 10" = 'darkgreen')

show_landscape_own = function(landscape_raster, color_palette = green_palette){
  loc_dat  = as.data.frame(landscape_raster, xy=T)
  loc_dat['Class'] = factor(loc_dat$layer_Categories, levels = names(green_palette))
  loc_dat%>%
    ggplot(aes(x=x, y=y, fill = Class))+
    geom_tile()+
    scale_fill_manual(values = color_palette)+
    ylab(' ')+
    xlab(' ')+
    theme_void()+
    theme(
      axis.title.x = element_blank(),       # Remove x-axis title
      axis.text.x = element_blank(),        # Remove x-axis text (labels)
      axis.ticks.x = element_blank(),       # Remove x-axis ticks
      axis.line.x = element_blank(),        # Remove x-axis line
      axis.title.y = element_blank(),       # Remove y-axis title
      axis.text.y = element_blank(),        # Remove y-axis text (labels)
      axis.ticks.y = element_blank(),       # Remove y-axis ticks
      axis.line.y = element_blank(),
      legend.position = 'bottom',
      legend.direction = 'horizontal'
    )
}


distributions_ = matrix(nrow = 4, ncol = 3)
distributions_[1, ] = c(.33, .33, .34)
distributions_[2, ] = c(0.1, .45, .45)
distributions_[3, ] = c(.1, .3, .6)
distributions_[4, ] = c(.1, .6, .3)


autocorrelations = c(.5, .7, .9, 1.3)


i = 1
for(autocorr_ in autocorrelations){
  for(distrib_ in 1:nrow(distributions_)){
      
    tryer = nlm_fbm(size, size, 1, fract_dim = autocorr_)
        
    classified_landscape = util_classify(tryer,
                                         weighting = distributions_[distrib_, ],
                                         level_names = c("Age 0 to 5",
                                                         "Age 5 to 10",
                                                         "Age > 10"))
    #classified_landscape@data@values  = classified_landscape@data@values - 1
        
        #classified_landscape = util_classify(tryer,
        #                                     n = 3,
        #                                     level_names = c("Age 0 to 5", 
        #                                                     "Age 5 to 10",
        #                                                     "Age > 10"))
      
    plot_ = show_landscape_own(classified_landscape)
    ggsave(here('Large_scale_landscapes', paste0('large_land_autocorr_',autocorr_,'_distrib_', distrib_,".jpg")), 
           plot_, 
           height = 12,
           width = 10,
           units = 'cm')
    
    to_save_data = Matrix(classified_landscape@data@values, nrow = size^2, ncol = size^2, byrow=T, sparse = T)
    writeMM(to_save_data,here('Large_scale_landscapes', paste0('large_land_', i,'.mtx')))
    
    print(paste('Landscape ', i, 'is done'))
    i = i+1
  }
}

gc()
  
for(distrib_ in 1:nrow(distributions_)){
  for(p_ in c(.3, .5, .7)){
    tryer = nlm_randomcluster(ncol = size, nrow= size, p = p_, ai = distributions_[distrib_,], neighbourhood = 8, rescale = F )
    
    classified_landscape = util_classify(tryer,
                                         weighting = distributions_[distrib_, ],
                                         level_names = c("Age 0 to 5",
                                                         "Age 5 to 10",
                                                         "Age > 10"))
    }
}


