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
  
  df <- df %>% gather(Original, PostCorrection, key = "Correction", value = 'Reflectance')
  
  ggplot(df, aes(x = l7, y = Reflectance, color = Correction)) + geom_point(alpha = .5) + 
    geom_abline(color = 'red') + 
    scale_color_viridis_d(begin = .2, end = .8) +
    stat_regline_equation(aes(label =  paste(..adj.rr.label..))) +
    theme_bw() +
    #scale_y_continuous(sec.axis = sec_axis(~.*(100/max.x), name = 'Quantile')) + 
    labs(y = paste0(capitalize(sat), " SR"), x = 'L7 SR', title = paste0(capitalize(sat),' ', capitalize(band), " Correction"))
}
  
 


## Old Correction plot function
# correctionPlot <- function(band, sat){
#   if(sat == 'l8'){
#     x <- refCompPre %>% filter(year > 2012)
#   }else if(sat == 'l5'){
#     x <- refCompPre %>% filter(year > 1999, year < 2012)
#   }
#   max.x <- max(x[,band])
#   df <- tibble(y = quantile(x %>% filter(sat == 'l7') %>% .[,band], 
#                             seq(.01,.99, .01)), 
#                x = quantile(x %>% filter(sat == sat) %>% .[,band], 
#                             seq(.01,.99, .01)))
#   lm <- lm(y~poly(x, 2), data = df)
#   r2 <- summary(lm)$adj.r.squared %>% round(4)
#   
#   df <- df %>% bind_cols(Corrected = lm$fitted.values) %>%
#     gather(x, Corrected, key = 'Reflectance', value = x) %>%
#     mutate(Reflectance = ifelse(Reflectance == 'x', 'Original', Reflectance))
#   
#   
#   ggplot(df, aes(x = x, y = y, color = Reflectance)) + geom_point(alpha = .5) + 
#     geom_abline(color = 'red') + 
#     annotate('text', label = paste0("italic(R)^2 ==", r2), parse = T, x = -Inf, y = Inf, vjust = 1.2, hjust = -.2) +
#     scale_color_viridis_d(begin = .2, end = .8) +
#     theme_bw() +
#     #scale_y_continuous(sec.axis = sec_axis(~.*(100/max.x), name = 'Quantile')) + 
#     labs(x = paste0(capitalize(sat), " SR"), y = 'L7 SR', title = paste0(capitalize(sat),' ', capitalize(band), " Correction")) +
#     labs(title = paste0(capitalize(band))) 
# }
