#############################################################################
#
# Suspended sediment budget and intra-event sediment dynamics of a small
# glaciated mountainous catchment in the Northern Caucasus
#
# Tsyplenkov A., Vanmaercke M., Golosov V., Chalov S.
#
# Part 1. Data reading and pre-processing
#
#############################################################################
Sys.setlocale("LC_ALL", "Russian_Russia")
Sys.setenv(TZ = "UTC")


library(tidyverse)
library(readxl)
library(purrr)
library(lubridate)
library(extrafont)
library(gghighlight)
library(corrr)
library(zoo)
library(openair)

source("R/00_own-functions.R")

# DATA LOAD ---------------------------------------------------------------
# Read SSC-NTU table
ssc_ntu <- readxl::read_xlsx("data/raw/ssc-ntu_2015_2017.xlsx") %>% 
  rename(ntu = 6)

ssc_ntu %>% 
  filter(!is.na(ntu), ntu < 1000) %>% 
  ggplot(aes(x = ntu, y = SSC)) +
  geom_point(size = 1.8) +
  geom_smooth(method = "lm",
              formula = y ~ poly(x, 2),
              linetype = "dashed",
              size = 1.1,
              color = "dimgrey",
              se = F) +
  ggpmisc::stat_poly_eq(aes(label =  paste(stat(eq.label),
                                           stat(adj.rr.label),
                                           sep = "~~~")),
                        formula = y ~ poly(x, 2), 
                        eq.x.rhs = "italic(T)",
                        eq.with.lhs = "italic(SSC)~`=`~",
                        parse = TRUE) +
  scale_x_continuous(breaks = c(0, 250, 500, 750, 1000)) +
  scale_y_continuous(labels = function(...) prettyNum(..., big.mark = " ")) +
  expand_limits(x = c(0, 0),
                y = c(0, 0)) +
  labs(x = expression(italic("T")*", NTU"),
       y = expression(italic(SSC)*", g"%.%m^-3)) +
  theme_clean() -> ssc_ntu1

ssc_ntu %>% 
  filter(!is.na(ntu), ntu >= 1000) %>% 
  ggplot(aes(x = ntu, y = SSC)) +
  geom_point(size = 1.8) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              linetype = "dashed",
              size = 1.1,
              color = "dimgrey",
              se = F) +
  ggpmisc::stat_poly_eq(aes(label =  paste(stat(eq.label),
                                           stat(adj.rr.label),
                                           sep = "~~~")),
                        eq.x.rhs = "italic(T)",
                        eq.with.lhs = "italic(SSC)~`=`~",
                        formula = y ~ x, 
                        parse = TRUE) +
  scale_y_continuous(labels = function(...) prettyNum(..., big.mark = " ")) +
  scale_x_continuous(breaks = c(1, 10, 20, 30) * 10^3,
                     labels = function(...) prettyNum(..., big.mark = " ")) +
  expand_limits(
    y = c(0, 0)) +
  labs(x = expression(italic("T")*", NTU"),
       y = expression(italic(SSC)*", g"%.%m^-3)) +
  theme_clean() -> ssc_ntu2

ggsave(
  "figures/fig02_ssc-ntu.png",
  ggpubr::ggarrange(ssc_ntu1, ssc_ntu2, labels = c("a", "b")),
  dpi = 500, h = 5, w = 10
)

# Create SSC-T models
ssc_ntu %>% 
  filter(!is.na(ntu), ntu < 1000) %>% 
  lm(SSC ~ poly(ntu, 2), data = .) -> lm1

ssc_ntu %>% 
  filter(!is.na(ntu), ntu >= 1000) %>% 
  lm(SSC ~ ntu, data = .) -> lm2

# Database can be accessed via doi.org/10.5063/F1H1307Q

## See Rets, E. P., Popovnin, V. V., Toropov, P. A., Smirnov,
## A. M., Tokarev, I. V., Chizhova, J. N., … Kornilova, E. D. (2018).
## Djankuat Glacier Station in the North Caucasus, Russia:
## A Database of complex glaciological, hydrological, meteorological
## observations and stable isotopes sampling results during 2007-2017.
## Earth System Science Data Discussions, 1–31.
## https://doi.org/10.5194/essd-2018-124

# Read djankuat hydrology data
df2017 <- readxl::read_xlsx("data/raw/Hydrology_Djan_Gauge_2007-2017.xlsx",
                            sheet = "2017", na = "n/a") %>% 
  magrittr::set_colnames(c("datetime", "q", "cond", "sal",
                           "t", "ntu", "ssc",
                           "do18", "dd")) %>% 
  filter(datetime >= as.POSIXct("2017-06-06 09:00:00", tz = "UTC"),
         datetime <= as.POSIXct("2017-09-24 23:00:00", tz = "UTC")) %>% 
  mutate(#ssc = ifelse(ntu > 1000, 1.2109*ntu - 749.54,
    #            0.0023*(ntu^2) - 0.2175*ntu + 218.7),
    ssc = ifelse(ntu > 1000,
                 predict(lm2, .),
                 predict(lm1, .)))

# Read data from MID and GL stations
mid <- read_delim("data/raw/data_mid.csv", delim = ";") %>% 
  mutate(datetime = glue::glue("{date} {time}"),
         datetime = as.POSIXct(datetime)) %>% 
  dplyr::select(datetime, ssc.mid)

gl <- read_delim("data/raw/data_gl.csv", delim = ";") %>% 
  mutate(datetime = glue::glue("{date} {time}"),
         datetime = as.POSIXct(datetime)) %>% 
  dplyr::select(datetime, ssc.gl, h.gl)

# Read data from rain collector
rain <- read_delim("data/raw/rain_out.csv", delim = ";") %>% 
  mutate(start = as.POSIXct(start, tz = "UTC"), # rain begin
         end = as.POSIXct(end, tz = "UTC"), # rain end
         p = as.double(p), # rainfall amount [mm]
         duration = difftime(end, start, tz = "UTC"),
         duration = as.double(duration) / 60, # length of an event [h]
         intensity = signif(p / duration, 3))  # intensity [mm/h]

out <- df2017 %>% dplyr::select(datetime, q, ssc)

# Complete 2017s database
full_join(out, mid, by = "datetime") %>% 
  full_join(gl, by = "datetime") -> df17 

rain %>% 
  rename(datetime = start) %>% 
  select(-end) %>% 
  full_join(df17, by = "datetime") %>% 
  arrange(datetime) %>% 
  dplyr::select(datetime, 5:9, 2:4) %>% 
  slice(-c(1:13)) -> df17

# SAVE -----------------------------------------------------------------------
rm(mid, gl, out, df2017)
save("df17", "rain", file =  "data/raw/df17_raw.Rdata")