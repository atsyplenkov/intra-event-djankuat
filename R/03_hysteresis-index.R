#############################################################################
#
# Suspended sediment budget and intra-event sediment dynamics of a small
# glaciated mountainous catchment in the Northern Caucasus
#
# Tsyplenkov A., Vanmaercke M., Golosov V., Chalov S.
#
# Part 3. Hysteresis index calculation
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
library(magrittr)

source("R/00_own-functions.R")
load("data/tidy/djan17.Rdata")

# Simple Hysteresis Index ------------------------------------------------------
df17 %>% 
  filter(!is.na(ssc),
         he == "75") %>%
  dplyr::select(1:3) %>% 
  mutate(ssc = ssc / 300) %>% 
  gather(variable, value, -datetime) %>% 
  ggplot(aes(x = datetime, y = value, color = variable)) +
  geom_line(size = 1.2) +
  scale_y_continuous(limits = c(0, 2.5),
                     sec.axis = sec_axis(~. * 300,
                                         breaks = seq(0, 750, 150),
                                         name = expression(italic(SSC)*","~ g %.% m^-3))) +
  labs(x = "", y = expression(italic(Q)*","~m^3%.%s^-1)) +
  ggsci::scale_color_nejm(name = "",
                          label = c("Q", "SSC")) +
  theme_clean() -> test_shi1

df17 %>% 
  filter(!is.na(ssc),
         he == "75") %>%
  dplyr::select(1:3) %>% 
  lm(ssc ~ q, data = .) -> hi_lm   

df17 %>% 
  filter(!is.na(ssc),
         he == "75") %>%
  dplyr::select(1:3) %>% 
  mutate(id = 1:n(),
         type = ifelse(id %in% c(1:which.max(.$q)), "rising", "falling"),
         ssc_pred = predict(hi_lm, .)) %>% 
  ggplot(aes(x = q, y = ssc)) +
  geom_segment(aes(x = q, y = ssc, xend = q, yend = ssc_pred),
               color = "grey60", linetype = "dashed") +
  geom_path(arrow = arrow(length = unit(3, "mm"), ends = "last")) +
  geom_point(aes(color = type), size = 2) +
  geom_smooth(method = "lm", se = F,
              color = "dimgrey",
              linetype = "longdash") +
  scale_y_continuous(limits = c(0, 750),
                     breaks = seq(0, 750, 150)) +
  labs(x = expression(italic(Q)*","~m^3%.%s^-1),
       y = expression(italic(SSC)*","~ g %.% m^-3))+
  ggsci::scale_color_jama(name = "",
                          label = c("Falling limb",
                                    "Rising limb")) +
  theme_clean() -> test_shi2

ggpubr::ggarrange(test_shi1, test_shi2,
                  labels = "AUTO",
                  nrow = 1,
                  align = "hv") %>% 
  ggsave(plot = .,
         filename = "figures/fig05_shi-example.png",
         dpi = 500, w = 10, h = 5)

# Apply SHI to the database
df17 %>% 
  filter(!is.na(ssc)) %>% 
  group_by(he = as.factor(he)) %>%
  mutate(n = n()) %>%
  filter(n >= 5) %>% 
  do(hi(.)) %>%
  summarise(SHI = round(mean(HI), 2) / max(ssc, na.rm = T)) %>% 
  ungroup() %>% 
  full_join(df17_db, by = "he") -> df17_db

# Visual identification of Hysteresis loops
# Plot Hysteresis events
# df17 %>% 
#   # mutate(ssc = na.approx(ssc, rule = 2)) %>%
#   filter(!is.na(ssc)) %>% 
#   ggplot(aes(x = q, y = ssc, group = he)) + 
#   geom_path(arrow = arrow(length = unit(3, "mm"), ends = "last")) +
#   # geom_point(shape = 1, colour = "black", fill = "white", stroke = 1, size = 2) +
#   labs(x = expression(italic(Q)*","~m^3%.%s^-1),
#        y = expression(italic(SSC)*","~ g %.% m^-3))+
#   theme_clean() -> hys2017
# 
# # Save result to .pdf book
# # Now you can make a visiual definition of loop type
# pdf("analysis/hysteresis_plots2017.pdf")
# ggplus::facet_multiple(plot=hys2017,
#                        facets="he",
#                        scales = "free",
#                        ncol = 2,
#                        nrow = 2)
# dev.off() 

# Write down loop types
df17_db %<>% 
  mutate(type = case_when(
    he %in% c(2:4, 5, 6:9, 12:14,
              17:24, 31, 36,
              40, 42:43, 45:48,
              52:56, 58:62,
              67:68, 70:71,
              73:75, 77:78, 80:81, 85,
              89:93, 95:96, 98, 101:103,
              106:109, 114:115, 117, 119:134) ~ "CW",
    he %in% c(26, 30, 34,
              41, 44, 
              49, 65, 69, 72,
              79, 83:84, 86, 94, 112, 116) ~ "AW",
    he %in% c(29, 50, 51,
              63:64, 66, 76) ~ "F8",
    TRUE ~ "NA"
  )) %>% 
  # Edit wrong assigned SHI values
  mutate(SHI = ifelse(SHI > 0 & type == "AW",
                      NA, SHI))

# Check the distribution of loop types and SHI
df17_db %>% 
  filter(!is.na(SHI)) %>% 
  mutate(outlier.high = SHI > quantile(SHI, .75) + 1.50*IQR(SHI),
         outlier.low = SHI < quantile(SHI, .25) - 1.50*IQR(SHI),
         outlier.color = case_when(outlier.high ~ "red",
                                   outlier.low ~ "red",
                                   outlier.low == F | outlier.high == F ~ "black")) %>% 
  ggplot(aes(x = type, y = SHI)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom ='errorbar',
               width = .25) +
  geom_jitter(aes(color = outlier.color),
              width = .1,
              alpha = .6,
              show.legend = F) +
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  ggsci::scale_color_lancet() +
  labs(x = "", y = expression(italic(SHI))) +
  theme_clean() -> shi_boxplot

# Create table with event characteristics
df17_db %>% 
  group_by(type) %>% 
  summarise(n = n(),
            n.pct = n*100/nrow(.),
            mean.duration = mean(length),
            mean.lag = mean(q.lag),
            sd.lag = sd(q.lag),
            shi.mean = mean(SHI, na.rm = T),
            shi.max = max(SHI, na.rm = T),
            shi.min = min(SHI, na.rm = T),
            q.mean = mean(q.out.mean, na.rm = T),
            q.max = max(q.out.max, na.rm = T),
            ssc.mean = mean(ssc.out.mean, na.rm = T),
            ssc.max = max(ssc.out.max, na.rm = T),
            sl.mean = mean(sl.out, na.rm = T),
            sl = sum(sl.out, na.rm = T),
            sl.pct = 100 * sum(sl.out, na.rm = T) /
              (df17_db %>% pull(sl.out) %>% sum(., na.rm = T))) %>% 
  mutate_if(is.numeric, list(~signif(., 3)))  %>% 
  mutate_if(is.numeric, list(~signif(.,3))) %>% 
  mutate_if(is.numeric, list(~prettyNum(., " "))) -> table4

# Transpose 
table4 %<>% 
  transpose() %>% 
  set_rownames(names(table4)) %>% 
  set_colnames(table4$type)

pairwise.wilcox.test(g = df17_db$type, x = df17_db$r)

# SAVE -------------------------------------------------------------------
openxlsx::loadWorkbook("analysis/tables.xlsx") -> wb
openxlsx::addWorksheet(wb, "Table 4")
openxlsx::writeData(wb, sheet = "Table 4", x = table4, rowNames = T)
openxlsx::saveWorkbook(wb, "analysis/tables.xlsx", overwrite = T)
