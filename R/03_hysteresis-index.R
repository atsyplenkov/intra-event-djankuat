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
  do(hi(.)) %>%
  summarise(SHI = round(mean(HI), 2) / max(ssc, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(df17_db, by = "he") -> df17_db

# Visual identification of Hysteresis loops
# Plot Hysteresis events
df17 %>% 
  # mutate(ssc = na.approx(ssc, rule = 2)) %>%
  filter(!is.na(ssc)) %>% 
  ggplot(aes(x = q, y = ssc, group = he)) + 
  geom_path(arrow = arrow(length = unit(3, "mm"), ends = "last")) +
  # geom_point(shape = 1, colour = "black", fill = "white", stroke = 1, size = 2) +
  labs(x = expression(italic(Q)*","~m^3%.%s^-1),
       y = expression(italic(SSC)*","~ g %.% m^-3))+
  theme_clean() -> hys2017

# Save result to .pdf book
# Now you can make a visiual definition of loop type
pdf("analysis/hysteresis_plots2017.pdf")
ggplus::facet_multiple(plot=hys2017,
                       facets="he",
                       scales = "free",
                       ncol = 2,
                       nrow = 2)
dev.off() 

# Write down loop types
df17_db %<>% 
  mutate(type = case_when(
    he %in% c(2:4, 6:9, 12:15, 17:24, 31, 36,
              40, 42:43, 45:47, 49,
              51:55, 58:60, 65,
              67:71, 73:74, 76:77,
              80, 82:87, 89,
              93:95, 97:100, 104,
              105:107, 109:120,
              122:124) ~ "CW",
    he %in% c(25, 26, 30, 34,
              41, 44, 48, 63,
              79, 88, 96) ~ "AW",
    he %in% c(5, 28, 29, 56, 61:62,
              64, 66, 72, 75, 78,
              121) ~ "F8",
    TRUE ~ "NA"
  ))
