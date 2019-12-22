#############################################################################
#
# Suspended sediment budget and intra-event sediment dynamics of a small
# glaciated mountainous catchment in the Northern Caucasus
#
# Tsyplenkov A., Vanmaercke M., Golosov V., Chalov S.
#
# Part 2. Data skimming
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
load("data/raw/df17_raw.Rdata")

# Create a column where an effect of rain on events will be described:
# SOURCE: https://stackoverflow.com/a/51884231/9300556
# TRUE - if there was rainfall during measurement
# FALSE - if there was no rain at all
# Overlap timeseries and time interval
df17 <- sqldf::sqldf('select df17.*, case when rain.p is not null then "TRUE" else "FALSE" end as factor_rain 
      from df17
      left join rain
      on start <= datetime and
      datetime <= end') %>% 
  as_tibble() 

# Hydrological events -------------------------------------------------------
df17 <- hydro_events(df17, q, datetime, 9)

# Hydrological event database -----------------------------------------------
xts::xts(df17[, -1], df17$datetime) %>% 
  timetk::tk_tbl() %>% 
  mutate(he = factor(he)) %>% 
  group_by(he) %>% 
  summarise(start = first(index),
            end = last(index),
            length = as.double(signif(difftime(end, start, units = "hours"), 3))) %>% 
  rowwise %>%
  mutate(mean_date = mean.POSIXct(c(start, end))) %>% 
  arrange(start) %>% 
  as_tibble() -> df17_db

# Calculate delay time of rainfall start and runoff peak
df17 %>% 
  filter(datetime > as.POSIXct("2017-06-30 9:00:00"),
         datetime < as.POSIXct("2017-08-30 9:00:00")) %>% 
  filter(factor_rain == T) %>% 
  group_by(he) %>% 
  summarise(peak.time = datetime[which.max(q)],
            rain.start = datetime[which(!is.na("p"))]) %>% 
  mutate(delay = as.double((peak.time - rain.start) / 3600)) %>%  # convert it to hours
  summarise(lag = signif(mean(delay), 3),
            med = median(delay),
            max = max(delay),
            min = min(delay),
            sd = sd(delay)
  ) %>% 
  pull(1) -> rain.lag

# Expand the duration of the rain by mean delay time
rain[, 2] <- rain[, 2] + rain.lag * 3600

df17 <- dplyr::select(df17, -factor_rain) # remove previous factor.rain
df17 <- sqldf::sqldf('select df17.*, case when rain.p is not null then "TRUE" else "FALSE" end as factor_rain 
      from df17
      left join rain
      on start <= datetime and
      datetime <= end') %>% 
  as_tibble()

# SSC graph analysis --------------------------------------------------------
# Filter data: remove all points, which may be influenced by rain (factor_rain = TRUE)
df17 %>% 
  mutate(ssc = zoo::na.approx(ssc, rule = 2)) %>% 
  filter(factor_rain != "TRUE" & ssc.mid > 1) %>% 
  as_tibble() %>% 
  mutate(datetime = format.POSIXct(datetime,
                                   "%Y-%m-%d %H:%M")) -> df17_norain

lbl <- c(1:20, 22:39, 41:44, 46:50, 52:54, 56:66, 68:81, 84:87)

df17_norain %>% 
  ggplot(aes(x = ssc.gl, y = ssc)) + 
  geom_point(size = 2, color = "#868686FF") +
  geom_abline(intercept = 0,             # add 1:1 line
              slope = 1, 
              color = "#EFC000FF",
              size = 1,
              alpha = 0.8) +
  ggrepel::geom_text_repel(data = df17_norain[-lbl,], # name outliers
                           aes(label = as.character(datetime)),
                           point.padding = unit(0.5, "lines"),
                           segment.colour = "black") +
  geom_point(data = df17_norain[-lbl,], # highlight outliers
             color = "#CD534CFF") +
  geom_smooth(data = df17_norain[lbl,],
              method = "lm",
              formula = y ~ x - 1,
              se = F,
              color = "#CD534CFF",           # add linear model result
              size = 1,
              alpha = 0.8,) +
  ggpmisc::stat_poly_eq(data = df17_norain[lbl,],
                        formula = y ~ x - 1,
                        parse = T,
                        eq.with.lhs = "italic(SSC[OUT])~`=`~",
                        eq.x.rhs = "`·`~italic(SSC[GL])",
                        aes(label =  paste(stat(eq.label),
                                           stat(adj.rr.label),
                                           sep = "~~~"))) +
  lims(y = c(0, 1000), x = c(0, 1000)) +
  xlab(expression(italic(SSC["ДАЛ"])*","~~ г %.% м^-3)) +
  ylab(expression(italic(SSC["ДЖАН"])*","~~ г %.% м^-3)) +
  theme_clean() -> ssc_ssc

ggsave("figures/rus/fig03_ssc-ssc_rus.png",
       plot = ssc_ssc,
       dpi = 500, w = 8, h = 6)

gl.ratio <- as.numeric(format(coef(lm(ssc ~ ssc.gl - 1,
                                      data = df17_norain[lbl,]))[1], digits = 2))

# Water discharge reconstruction -----------------------------------------------
# Subset all non-rainy days at GL station
df17 %>% 
  mutate(q.gl = ifelse(factor_rain == F & ssc.gl > 1,
                       q * gl.ratio, NA)) -> df17

df17 %>% 
  filter(factor_rain == F & ssc.gl > 1) %>%
  ggplot(aes(x = h.gl, y = q.gl)) +
  geom_point(color = "#868686FF", size = 2) +
  geom_smooth(method = "lm", se = F,
              color='#CD534CFF',
              size = 1.4,
              alpha = 0.8) +
  ggpmisc::stat_poly_eq(formula = y ~ x,
                        parse = T,
                        eq.with.lhs = "italic(Q[ДАЛ])~`=`~",
                        eq.x.rhs = "`·`~italic(H[ДАЛ])",
                        aes(label =  paste(stat(eq.label),
                                           stat(adj.rr.label),
                                           sep = "~~~"))) +
  xlim(70, 83) +
  ylim(0, 2) +
  labs(x = expression(italic(H[ДАЛ])*","*~см),
       y = expression(italic(Q[ДАЛ])*","*~м^3%.%с^"-1")) +
  theme_clean() -> qfh

ggsave("figures/rus/fig04_qfh_rus.png", plot = qfh,
       dpi = 500, w = 5, h = 5)

# Build linear model q.gl = f(h.gl)
qfh <- lm(q.gl ~ h.gl, data = filter(df17, factor_rain == F & ssc.gl > 1))

df17 %>% 
  # Calculate water discharge at MID and GL
  dplyr::mutate(q.gl = ifelse(is.na(q.gl), predict(qfh,.), q.gl),
                q.mid = q) %>% 
  # Calculate sediment discharge
  dplyr::mutate(r = ssc * q,
                r.mid = ssc.mid * q.mid,
                r.gl = ssc.gl * q.gl) %>% 
  dplyr::group_by(he) %>% 
  dplyr::summarise(r = mean(r, na.rm = T),
                   r.mid = mean(r.mid, na.rm = T),
                   r.gl = mean(r.gl, na.rm = T),
                   ssc.out.mean = mean(ssc, na.rm = T),
                   ssc.out.max = max(ssc, na.rm = T),
                   q.out.mean = mean(q, na.rm = T),
                   q.out.max = max(q, na.rm = T),
                   p = sum(p, na.rm = T),
                   peak.q = datetime[which.max(q)],
                   peak.ssc = datetime[which.max(ssc)[1]],
                   q.lag = as.double(difftime(peak.q, peak.ssc, "hours")))  %>% 
  dplyr::select(-peak.q, -peak.ssc) %>% 
  dplyr::left_join(df17_db, by = "he") %>% 
  dplyr::mutate(sl.out = round((r * 3600 * length / 10^6), 2),
                sl.mid = round((r.mid * 3600 * length / 10^6), 2),
                sl.gl = round((r.gl * 3600 * length / 10^6), 2),
                # Edit Inf values in SSC's
                ssc.out.max = ifelse(is.infinite(ssc.out.max),
                                     NA, ssc.out.max)) %>% 
  dplyr::mutate_if(is.numeric, list(~signif(., 2))) -> df17_db

is.na(df17_db) <- sapply(df17_db, is.na)

# Data skimming --------------------------------------------------------
df17 %>% 
  dplyr::select(1:5, q.gl) %>% 
  gather(variable, value, -datetime) %>% 
  group_by(variable) %>% 
  summarise(Mean = mean(value, na.rm = T),
            SD = sd(value, na.rm = T),
            Median = median(value, na.rm = T),
            Max = max(value, na.rm = T),
            Max.date = datetime[which.max(value)],
            Min = min(value, na.rm = T),
            Min.date = datetime[which.min(value)]) %>%
  mutate_if(is.numeric, list(~signif(.,3))) %>% 
  mutate_if(is.numeric, list(~prettyNum(., " "))) %>% 
  mutate_if(is.Date, list(~format.POSIXct("%F %H:%M"))) %>% 
  mutate_at(vars(Max.date, Min.date),
            list(~as.character.POSIXt(.)))-> table2

df17 %>% 
  dplyr::select(1:5, q.gl) %>% 
  mutate(day = lubridate::date(datetime),
         sl.out = q * ssc / 1000, # sediment load at the OUT [kg/s]
         sl.mid = q * ssc.mid / 1000, # sediment load at the MID [kg/s]
         sl.gl = q.gl * ssc.gl / 1000) %>% # sediment load at the GL [kg/s]
  group_by(day) %>% 
  summarise(sl.out = mean(sl.out, na.rm = T) * 86400 / 1000,
            sl.mid = mean(sl.mid, na.rm = T) * 86400 / 1000,
            sl.gl = mean(sl.gl, na.rm = T) * 86400 / 1000) %>% 
  filter(day >= as.Date("2017-06-30"),
         day <= as.Date("2017-08-30")) %>% 
  gather(station, sl, -day) %>% 
  group_by(station) %>% 
  summarise(`Observation period` = sum(sl, na.rm = T),
            Total = 1.02 * `Observation period`,
            Mean = mean(sl, na.rm = T),
            SD = sd(sl, na.rm = T),
            Median = median(sl, na.rm = T),
            Max = max(sl, na.rm = T),
            Max.date = day[which.max(sl)],
            Min = min(sl, na.rm = T),
            Min.date = day[which.min(sl)]) %>%
  mutate_if(is.numeric, list(~signif(.,3))) %>% 
  mutate_if(is.numeric, list(~prettyNum(., " "))) %>% 
  mutate_at(vars(Max.date, Min.date),
            list(~as.character.POSIXt(.))) %>% 
  mutate(Total = ifelse(station == "sl.out", Total, NA)) -> table3

df17 %>% 
  dplyr::select(1:5, q.gl, p) %>% 
  mutate(day = lubridate::date(datetime),
         sl.out = q * ssc / 1000, # sediment load at the OUT [kg/s]
         sl.mid = q * ssc.mid / 1000, # sediment load at the MID [kg/s]
         sl.gl = q.gl * ssc.gl / 1000) %>% # sediment load at the GL [kg/s]
  group_by(day) %>% 
  summarise(sl.out = mean(sl.out, na.rm = T) * 86400 / 1000,
            sl.mid = mean(sl.mid, na.rm = T) * 86400 / 1000,
            sl.gl = mean(sl.gl, na.rm = T) * 86400 / 1000,
            p = sum(p, na.rm = T)) %>% 
  mutate(rain = ifelse(p != 0, "yes", "no")) %>% 
  filter(day >= as.Date("2017-06-30"),
         day <= as.Date("2017-08-30")) %>%
  gather(station, sl, -day, -rain) %>%
  filter(station != "p") %>% 
  group_by(rain, station) %>% 
  summarise(`Observation period` = sum(sl, na.rm = T),
            Total = 1.02 * `Observation period`,
            Mean = mean(sl, na.rm = T),
            SD = sd(sl, na.rm = T),
            Median = median(sl, na.rm = T),
            Max = max(sl, na.rm = T),
            Max.date = day[which.max(sl)],
            Min = min(sl, na.rm = T),
            Min.date = day[which.min(sl)]) %>% 
  mutate_if(is.numeric, list(~signif(.,3))) %>% 
  mutate_if(is.numeric, list(~prettyNum(., " "))) %>% 
  mutate_at(vars(Max.date, Min.date),
            list(~as.character.POSIXt(.))) %>% 
  mutate(Total = ifelse(station == "sl.out", Total, NA)) -> table35
  

# Hydrograph plot ------------------------------------------------------------
# Sys.setlocale("LC_TIME", "English")

df17 %>% 
  mutate_at(vars(q, ssc), list(~zoo::na.approx(., rule = 2))) %>%
  mutate(p = ifelse(is.na(p), 0, p),
         ssc = log10(ssc)) %>% 
  select(datetime, q, ssc, p) %>% 
  gather(var, val, -datetime) %>% 
  ggplot(aes(x = datetime, y = val, color = var)) +
  geom_line() +
  facet_wrap(~var, scales = "free_y", nrow = 3, strip.position = "right",
             labeller = labeller(var = c(
               p = "Precipitation",
               q = "Water discharge",
               ssc = "log10(SSC)"
             ))) +
  labs(x = "", y = "") +
  ggsci::scale_color_nejm() +
  scale_x_datetime(date_breaks = "7 days",
                   date_labels = "%b %d") +
  theme_clean(legend = "none") -> hydrograph

ggsave("figures/fig06_hydrograph.png", hydrograph,
       dpi = 500, w = 10, h = 6)
  
# SAVE -------------------------------------------------------------------------
save("df17", "df17_db", file =  "data/tidy/djan17.Rdata")

openxlsx::write.xlsx(list(table2, table3, table35),
                     "analysis/tables.xlsx")
