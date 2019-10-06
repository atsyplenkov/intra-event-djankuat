# Boxplot upgrade ---------------------------------------------------------------
# Source: https://owi.usgs.gov/blog/boxplots/

n_fun <- function(x, y_position = log10(10^4)){
  return(data.frame(y = 0.99*y_position,
                    label = paste0("n = ", length(x))))
}

prettyLogs <- function(x){
  pretty_range <- range(x[x > 0])
  pretty_logs <- 10^(-10:10)
  log_index <- which(pretty_logs < pretty_range[2] & 
                       pretty_logs > pretty_range[1])
  log_index <- c(log_index[1]-1,log_index, log_index[length(log_index)]+1)
  pretty_logs_new <-  pretty_logs[log_index] 
  return(pretty_logs_new)
}

fancyNumbers <- function(n){
  nNoNA <- n[!is.na(n)]
  x <-gsub(pattern = "1e",replacement = "10^",
           x = format(nNoNA, scientific = TRUE))
  exponents <- as.numeric(sapply(strsplit(x, "\\^"), function(j) j[2]))
  
  base <- ifelse(exponents == 0, "1", ifelse(exponents == 1, "10","10^"))
  exponents[base == "1" | base == "10"] <- ""
  textNums <- rep(NA, length(n))  
  textNums[!is.na(n)] <- paste0(base,exponents)
  
  textReturn <- parse(text=textNums)
  return(textReturn)
}

# Hydrological events ----------------------------------------------------------
# Defining Local Minimum from HYSEP
# https://github.com/USGS-R/DVstats/blob/master/R/hysep.R
locmin <- function(x, datetime, window = 1) {
  
  timestep <- as.double(signif(difftime(head(datetime)[5],
                                        head(datetime)[4],
                                        units = "hours"), 4))
  
  N2star <- round(window / timestep)
  N2star <- ifelse(N2star %% 2 == 0, N2star  + 1, N2star)
  Nobs <- length(x)
  Ngrp <- ceiling(Nobs / N2star)
  Nfil <- (N2star - 1L) / 2L
  Mid <- as.integer((N2star) / 2)
  LocMin <- sapply(seq(N2star, Nobs), function(i)
    min(x[seq(i - N2star + 1L, i)]) == x[i - Mid]
  )
  LocMin <- c(rep(FALSE, Nfil), LocMin, rep(FALSE, Nfil))
  return(LocMin)
  
}

# HYDROLOGICAL EVENTS -------------------------------------------------------
hydro_events <- function(dataframe,
                         q = NULL,
                         datetime = NULL,
                         window = 1) {
  
  dataframe %>% 
    mutate(q = zoo::na.approx(q, rule = 2)) %>% 
    mutate(LocMin = locmin(x = q,
                           datetime = datetime,
                           window = window),
           test = ifelse(LocMin == FALSE, NA, 1)) -> dataframe
  
  # Remove multiple local minimums
  dataframe$test <- lapply(seq(1, length(dataframe$test)), function(i)
    dataframe$test[i - 1L] == dataframe$test[i])
  
  dataframe$LocMin[dataframe$test == T] <- FALSE
  
  # Name the hydrological events
  dataframe$he <- NA
  dataframe$he[dataframe$LocMin == T] <- 2:(length(dataframe$he[dataframe$LocMin == T])+1)
  dataframe$he[1] <- 1
  
  # Fill the gaps of he's
  dataframe %>%
    mutate(he = zoo::na.locf(he),
           he = as.factor(he)) %>% 
    dplyr::select(-LocMin, - test) %>% 
    as_tibble() -> dataframe 
  
  return(dataframe)
}

# Hysteresis Index calculation -----------------------------------------------
hi <- function(i) {
  M <- lm(log(ssc) ~ I(log(q)), data = i)
  a <- exp(coef(M)[1])
  b <- coef(M)[2]
  res <- i$ssc - a*i$q^b
  r <- mean(res[1:which.max(i$q)])
  f <- mean(res[c((which.max(i$q)+1):length(i$q))])
  i$HI <- r-f
  return(i)
}

# GGPLOT2 THEME -----------------------------------------------------------
theme_clean <- function(base_font_family = "Ubuntu",
                        base_font_size = 12,
                        legend = "bottom") {
  
  ggpubr::theme_pubclean(base_family = base_font_family,
                         base_size = base_font_size) +
    theme(
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.position = legend,
      strip.background = element_blank()
    )
}

# Insert minor blank ticks to GGPLOT2 -------------------------------------
# Source:: https://stackoverflow.com/questions/34533472/insert-blanks-into-a-vector-for-e-g-minor-tick-labels-in-r
every_nth <- function(x, nth, empty = TRUE, inverse = TRUE) 
{
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if(empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}