#############################################################################
#
# Suspended sediment budget and intra-event sediment dynamics of a small
# glaciated mountainous catchment in the Northern Caucasus
#
# Tsyplenkov A., Vanmaercke M., Golosov V., Chalov S.
#
# Part 2. Sediment budget reconstruction
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
  xlim(70, 83) +
  ylim(0, 2) +
  labs(x = expression(italic(H[GL])*","*~cm),
       y = expression(italic(Q[GL])*","*~m^3%.%s^"-1")) +
  theme_clean() -> qfh

ggsave("figures/fig04_qfh.png", plot = qfh,
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
                   p = sum(p, na.rm = T))  %>% 
  dplyr::left_join(df17_db, by = "he") %>% 
  dplyr::mutate(w.out = round((r * 3600 * length / 10^6), 2),
                w.mid = round((r.mid * 3600 * length / 10^6), 2),
                w.gl = round((r.gl * 3600 * length / 10^6), 2)) %>% 
  dplyr::mutate_if(is.numeric, list(~signif(., 2))) -> df17_db

is.na(df17_db) <- sapply(df17_db, is.na)

# Plot sediment budget
df17_db %>% 
  mutate(w.out = ifelse(w.out >= 1800, 720, w.out),
         w.out = na.approx(w.out),
         he = as.numeric(he),
         w.mid = ifelse(he > 25 & he < 108, na.approx(w.mid, rule = 2), NA),
         w.gl = ifelse(he > 25 & he < 108, na.approx(w.gl, rule = 2), NA)) %>% 
  dplyr::select(he, w.out, w.mid, w.gl, p) %>% 
  gather(station, sl, -he, -p) %>% 
  ggplot(aes(x = he, y = sl)) +
  geom_col(aes(y = p * 2), fill = "grey80") +
  geom_line(aes(group = factor(station), color = factor(station)),
            size = 1) +
  # geom_col(aes(fill = factor(station)),
  #          position = position_dodge(preserve = "total"),           
  #          alpha = .7) +
  coord_cartesian(ylim = 0:700) +
  geom_curve(data = tibble(x = 100, y = 620, xend = 111, yend = 720),
             aes(x = x, y = y, xend = xend, yend = yend),
             curvature = -.3,
             color = "black",
             inherit.aes = F,
             arrow = arrow(type = "closed", length = unit(.25, "cm"))) +
  annotate("text", x = 95, y = 615, hjust = 0, vjust = 1, lineheight = .9,
           label = paste0(round(max(df17_db$w.out, na.rm = T), 0), " ton"),
           size = 4, family = "Ubuntu", color = "black") +
  labs(y = expression(italic("Sediment load")*","~t),
       x = "Hydrological event") +
  scale_y_continuous(sec.axis = sec_axis(~./6,
                                         breaks = c(0, 200/6, 400/6, 100),
                                         labels = c(0, 30, 60, 100),
                                         name = expression(italic(p)*", mm"))) +
  scale_x_continuous(breaks = seq(0, nrow(df17_db), 5)) +
  scale_color_manual(name = "Sediment load at",
                     labels = c("GL", "MID", "OUT"),
                     values = c("w.gl" = "#0073C2FF",
                                "w.mid" = "#EFC000FF",
                                "w.out" = "#868686FF")) +
  theme_clean() -> sed_budget17

ggsave("figures/fig05_sed-budget17.png", plot = sed_budget17,
       dpi = 500, w = 10, h = 5)


# SAVE
save("df17", "df17_db", file =  "data/tidy/djan17.Rdata")
