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
  xlab(expression(italic(SSC[GL])*","~~ g %.% m^-3)) +
  ylab(expression(italic(SSC[OUT])*","~~ g %.% m^-3)) +
  theme_clean() -> ssc_ssc

ggsave("figures/fig03_ssc-ssc.png", plot = ssc_ssc,
       dpi = 500, w = 8, h = 6)

gl.ratio <- as.numeric(format(coef(lm(ssc ~ ssc.gl - 1,
                                      data = df17_norain[lbl,]))[1], digits = 2))

# Water discharge reconstruction -----------------------------------------------
# Read q = f(H) data from OUT station
read_xlsx("data/raw/discharges_dja_2017_v3.xlsx", sheet = "Q(H)",
                     range = "A2:D62") %>% 
  select(-2) %>% 
  magrittr::set_colnames(c("date", "h", "q")) %>% 
  mutate(date = force_tz(date, "Europe/Moscow")) %>% 
  filter(!is.na(q), !str_detect(q, "-")) %>% 
  mutate(date = na.locf(date),
         Month = month(date, label = T, locale = "English-US"),
         q = as.double(q)) %>%  
  ggplot(aes(x = h, y = q, color = Month)) +
  geom_point() +
  geom_smooth(method = "lm",
              se = F,
              size = 1.4,
              alpha = 0.8) +
  ggpmisc::stat_poly_eq(formula = y ~ x,
                        parse = T, vstep = 0.05, 
                        eq.with.lhs = "italic(Q[OUT])~`=`~",
                        eq.x.rhs = "`·`~italic(H[OUT])",
                        aes(label =  paste(stat(eq.label),
                                           stat(adj.rr.label),
                                           sep = "~~~"))) +
  scale_y_continuous(breaks = seq(0, 3, 0.5),
                     lim = c(0, 3)) +
  labs(x = expression(italic(H[OUT])*","*~cm),
       y = expression(italic(Q[OUT])*","*~m^3%.%s^"-1")) +
  ggsci::scale_color_npg() +
  # scale_color_manual(values = nationalparkcolors::park_palette("SmokyMountains")) +
  theme_clean() -> qfh.out

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
                        eq.with.lhs = "italic(Q[GL])~`=`~",
                        eq.x.rhs = "`·`~italic(H[GL])",
                        aes(label =  paste(stat(eq.label),
                                           stat(adj.rr.label),
                                           sep = "~~~"))) +
  scale_y_continuous(breaks = seq(0, 3, 0.5),
                     lim = c(0, 3)) +
  xlim(65, 85) +
  labs(x = expression(italic(H[GL])*","*~cm),
       y = expression(italic(Q[GL])*","*~m^3%.%s^"-1")) +
  theme_clean() -> qfh.gl

ggpubr::ggarrange(qfh.out, qfh.gl, ncol = 2, align = "v",
                  common.legend = T, legend = "bottom",
                  labels = c("a", "b")) %>% 
ggsave(filename = "figures/fig04_qfh.png",
       plot = .,
       dpi = 1000,
       w = 10, h = 5)

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
                   intensity = max(intensity, na.rm = T),
                   peak.q = datetime[which.max(q)],
                   peak.ssc = datetime[which.max(ssc)[1]],
                   q.lag = as.double(difftime(peak.q, peak.ssc, "hours")))  %>% 
  dplyr::select(-peak.q, -peak.ssc) %>% 
  dplyr::left_join(df17_db, by = "he") %>% 
  dplyr::mutate(sl.out = round((r * 3600 * length / 10^6), 2),
                sl.mid = round((r.mid * 3600 * length / 10^6), 2),
                sl.gl = round((r.gl * 3600 * length / 10^6), 2),
                # Edit Inf values in SSC's
                intensity = ifelse(is.infinite(intensity),
                                     NA, intensity),
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
            list(~as.character.POSIXt(.))) -> table2

# Calculate suspended load statistics as event based
df17_db %>% 
  dplyr::select(he, sl.out, sl.mid, sl.gl) %>% 
  gather(station, sl, -he) %>% 
  group_by(station) %>% 
  summarise(n = sum(!is.na(sl)),
            `Observation period` = sum(sl, na.rm = T),
            Total = 1.02 * `Observation period`,
            Mean = mean(sl, na.rm = T),
            SD = sd(sl, na.rm = T),
            Median = median(sl, na.rm = T),
            Max = max(sl, na.rm = T),
            # Max.date = he[which.max(sl)],
            Max.date = paste0(df17_db$start[he[which.max(sl)]],
                              " — ",
                              df17_db$end[he[which.max(sl)]]),
            Min = min(sl, na.rm = T),
            # Min.date = he[which.min(sl)],
            Min.date = paste0(df17_db$start[he[which.min(sl)]],
                              " — ",
                              df17_db$end[he[which.min(sl)]])) %>% 
  mutate_if(is.numeric, list(~round(., 2))) %>% 
  mutate_if(is.numeric, list(~prettyNum(., " "))) %>% 
  mutate_at(vars(Max.date, Min.date),
            list(~as.character.POSIXt(.))) %>% 
  mutate(Total = ifelse(station == "sl.out", Total, NA)) -> table3_he

# Calculate suspended load statistics as daily based
# 
# df17 %>% 
#   dplyr::select(1:5, q.gl) %>% 
#   mutate(day = lubridate::date(datetime),
#          sl.out = q * ssc / 1000, # sediment load at the OUT [kg/s]
#          sl.mid = q * ssc.mid / 1000, # sediment load at the MID [kg/s]
#          sl.gl = q.gl * ssc.gl / 1000) %>% # sediment load at the GL [kg/s]
#   group_by(day) %>% 
#   summarise(sl.out = mean(sl.out, na.rm = T) * 86400 / 1000,
#             sl.mid = mean(sl.mid, na.rm = T) * 86400 / 1000,
#             sl.gl = mean(sl.gl, na.rm = T) * 86400 / 1000) %>% 
#   gather(station, sl, -day) %>% 
#   group_by(station) %>% 
#   summarise(n = sum(!is.na(sl)),
#             `Observation period` = sum(sl, na.rm = T),
#             Total = 1.02 * `Observation period`,
#             Mean = mean(sl, na.rm = T),
#             SD = sd(sl, na.rm = T),
#             Median = median(sl, na.rm = T),
#             Max = max(sl, na.rm = T),
#             Max.date = day[which.max(sl)],
#             Min = min(sl, na.rm = T),
#             Min.date = day[which.min(sl)]) %>%  
#   mutate_if(is.numeric, list(~signif(., 5))) 
#   mutate_if(is.numeric, list(~prettyNum(., " "))) %>% 
#   mutate_at(vars(Max.date, Min.date),
#             list(~as.character.POSIXt(.))) %>% 
#   mutate(Total = ifelse(station == "sl.out", Total, NA)) -> table3_day

# Hydrograph plot ------------------------------------------------------------
Sys.setlocale("LC_TIME", "English")

# Read and transform air temperature data
temp <- readxl::read_xlsx("data/raw/base_camp_AWS.xlsx",
                          sheet = "2017-daily") %>% 
  select(datetime = 1, t = 2) %>% 
  mutate(datetime = force_tz(datetime, "Europe/Moscow"),
         datetime = floor_date(datetime, unit = "minutes"),
         t = as.numeric(t)) %>% 
  filter(datetime != "2017-06-06 00:00:00") %>%
  mutate(t = zoo::na.approx(t, rule = 2)) %>% 
  rename(date = datetime) %>% 
  openair::timeAverage(mydata = ., avg.time = "hour") %>% 
  rename(datetime = date) %>% 
  arrange(datetime) %>% 
  mutate(t = zoo::na.approx(t, rule = 2))

# Plot timeseries of raindall, air temperature, water discharge and SSC
full_join(df17, temp, by = "datetime") %>% 
  arrange(datetime) %>% 
  mutate(t = zoo::na.approx(t),
         p = ifelse(is.na(p), NA, p)) %>%
         # ) %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = t, color = "Air temperature",
                fill = "Air temperature"), na.rm = T) +
  # geom_line(aes(x = datetime, y = p/10, color = "Precipitation")) +
  geom_col(aes(x = datetime, y = p/10, color = "Precipitation",
               fill = "Precipitation")) +
  scale_y_continuous(name = "Air temperature, °C",
                     sec.axis = sec_axis(~.*10, name = "Precipitation, mm")) +
  scale_x_datetime(name = "",
                   date_breaks = "7 days",
                   date_labels = "%b %d",
                   limits = c(as.POSIXct("2017-06-06 12:00:00"),
                              as.POSIXct("2017-09-24 00:00:00"))) +
  ggsci::scale_color_nejm(name = "") +
  ggsci::scale_fill_nejm(name = "") +
  theme_clean() -> temp_rain

df17 %>% 
  mutate_at(vars(q, ssc), list(~zoo::na.approx(., rule = 2))) %>%
  mutate(ssc = log10(ssc)) %>% 
  ggplot() +
  geom_line(aes(x = datetime, y = ssc-2, color = "Suspended sediment concentration")) +
  geom_line(aes(x = datetime, y = q, color = "Water discharge"), na.rm = T) +
  scale_y_continuous(name = expression(italic(Q[OUT])*","~m^3%.%s^-1),
                     sec.axis = sec_axis(~.+2,
                                         name = expression(log[10]*italic(SSC[OUT])*","~ g %.% m^-3))) +
  scale_x_datetime(name = "",
                   date_breaks = "7 days",
                   date_labels = "%b %d",
                   limits = c(as.POSIXct("2017-06-06 12:00:00"),
                              as.POSIXct("2017-09-24 00:00:00"))) +
  ggsci::scale_color_npg(name = "") +
  theme_clean() -> hydro_ssc

ggpubr::ggarrange(temp_rain, hydro_ssc,
                  ncol = 1,
                  align = "h",
                  legend = "bottom",
                  labels = c("a", "b")) %>% 
  ggsave(filename = "figures/fig06_hydrograph.png", plot = ., dpi = 500,
         w = 10, h = 6)

# SAVE -------------------------------------------------------------------------
save("df17", "df17_db", file =  "data/tidy/djan17.Rdata")

openxlsx::write.xlsx(list(table2, table3_he),
                     "analysis/tables.xlsx")