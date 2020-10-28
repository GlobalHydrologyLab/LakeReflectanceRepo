fui.hue <- function(R, G, B) {
  
  # Convert R,G, and B spectral reflectance to dominant wavelength based
  # on CIE chromaticity color space
  
  # see Wang et al 2015. MODIS-Based Radiometric Color Extraction and
  # Classification of Inland Water With the Forel-Ule
  # Scale: A Case Study of Lake Taihu
  
  require(colorscience)
  # chromaticity.diagram.color.fill()
  Xi <- 2.7689*R + 1.7517*G + 1.1302*B
  Yi <- 1.0000*R + 4.5907*G + 0.0601*B
  Zi <- 0.0565*G + 5.5943*B
  
  # calculate coordinates on chromaticity diagram
  x <-  Xi / (Xi + Yi +  Zi)
  y <-  Yi / (Xi + Yi +  Zi)
  z <-  Zi / (Xi + Yi +  Zi)
  
  # calculate hue angle
  alpha <- atan2((x - 0.33), (y - 0.33)) * 180/pi
  
  # make look up table for hue angle to wavelength conversion
  cie <- cccie31 %>%
    mutate(a = atan2( (x - 0.33), (y - 0.33)) * 180/pi) %>%
    dplyr::filter(wlnm <= 700) %>%
    dplyr::filter(wlnm >=380)
  
  # find nearest dominant wavelength to hue angle
  wl <- cie[as.vector(sapply(alpha,function(x) which.min(abs(x - cie$a)))), 'wlnm']
  
  #out <- cbind(as.data.frame(alpha), as.data.frame(wl))
  
  return(wl)
}

## Function for comparing sensors during coincident time periods
correctionPlot <- function(band, sat, dataPre, dataPost){
  if(sat == 'l8'){

    df <- tibble(l7 = quantile(dataPost[sat == 'l7' & year > 2012, ..band
                           ][[1]], seq(.01,.99,.01)),
    Original = quantile(dataPre[sat == sat & year > 2012, ..band
                                 ][[1]], seq(.01,.99,.01)),
    PostCorrection = quantile(dataPost[sat == sat & year > 2012, ..band
                                       ][[1]], seq(.01,.99,.01)))
  }else if(sat == 'l5'){

    df <- tibble(l7 = quantile(dataPost[sat == 'l7' & year < 2012 & year > 1999, ..band
                           ][[1]], seq(.01,.99,.01)),
    Original = quantile(dataPre[sat == sat & year < 2012 & year > 1999, ..band
                                ][[1]], seq(.01,.99,.01)),
    PostCorrection = quantile(dataPost[sat == sat & year < 2012 & year > 1999, ..band
                                       ][[1]], seq(.01,.99,.01)))
  }
  
  ogBias <- round(Metrics::bias(df$l7, df$Original), 3)
  CBias <- round(Metrics::bias(df$l7, df$PostCorrection), 3)
  
  df <- df %>% gather(Original, PostCorrection, key = "Correction", value = 'Reflectance')
  
  ggplot(df, aes(x = l7, y = Reflectance, color = Correction)) + geom_point(alpha = .8) + 
    geom_abline(color = 'red') + 
    scale_color_viridis_d(end = .7, labels = c('Original', 'Post\nCorrection')) +
    #stat_regline_equation(aes(label =  paste(..adj.rr.label..))) +
    annotate('text', x= Inf, y = Inf, vjust = 4.5, hjust = 1, 
             label = paste0('Original Bias: ', ogBias,'\nCorrected Bias: ', CBias)) +
    theme_bw() +
    theme(axis.title = element_blank()) +
    scale_y_continuous(trans = 'log10') +
    scale_x_continuous(trans = 'log10') #
    #labs(y = paste0(capitalize(sat), " SR"), x = 'L7 SR', title = paste0(capitalize(sat),' ', capitalize(band), " Correction"))
}
  

## Code for connnecting dWL to forel-ule index 
fui.lookup <- tibble(dWL = c(471:583), fui = NA)

fui.lookup$fui[fui.lookup$dWL <= 583] = 21
fui.lookup$fui[fui.lookup$dWL <= 581] = 20
fui.lookup$fui[fui.lookup$dWL <= 579] = 19
fui.lookup$fui[fui.lookup$dWL <= 577] = 18
fui.lookup$fui[fui.lookup$dWL <= 575] = 17
fui.lookup$fui[fui.lookup$dWL <= 573] = 16
fui.lookup$fui[fui.lookup$dWL <= 571] = 15
fui.lookup$fui[fui.lookup$dWL <= 570] = 14
fui.lookup$fui[fui.lookup$dWL <= 569] = 13
fui.lookup$fui[fui.lookup$dWL <= 568] = 12
fui.lookup$fui[fui.lookup$dWL <= 567] = 11
fui.lookup$fui[fui.lookup$dWL <= 564] = 10
fui.lookup$fui[fui.lookup$dWL <= 559] = 9
fui.lookup$fui[fui.lookup$dWL <= 549] = 8
fui.lookup$fui[fui.lookup$dWL <= 530] = 7
fui.lookup$fui[fui.lookup$dWL <= 509] = 6
fui.lookup$fui[fui.lookup$dWL <= 495] = 5
fui.lookup$fui[fui.lookup$dWL <= 489] = 4
fui.lookup$fui[fui.lookup$dWL <= 485] = 3
fui.lookup$fui[fui.lookup$dWL <= 480] = 2
fui.lookup$fui[fui.lookup$dWL <= 475 & fui.lookup$dWL >470] = 1


# Actual Forel-Ule Colors
fui.colors <- c(
  "#2158bc", "#316dc5", "#327cbb", "#4b80a0", "#568f96", "#6d9298", "#698c86", 
  "#759e72", "#7ba654", "#7dae38", "#94b660","#94b660", "#a5bc76", "#aab86d", 
  "#adb55f", "#a8a965", "#ae9f5c", "#b3a053", "#af8a44", "#a46905", "#9f4d04")

