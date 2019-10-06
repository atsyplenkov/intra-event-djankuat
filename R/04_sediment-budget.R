#############################################################################
#
# Suspended sediment budget and intra-event sediment dynamics of a small
# glaciated mountainous catchment in the Northern Caucasus
#
# Tsyplenkov A., Vanmaercke M., Golosov V., Chalov S.
#
# Part 4. Sediment budget analysis
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

# Sediment budget
df17_db %>% 
  mutate(dw = sl.gl - sl.out,
         erosion = ifelse(dw > 0, "deposition", "erosion")) %>%
  group_by(erosion) %>% 
  summarise(n = n(),
            n.pct = 100*n/58,
            w = sum(dw, na.rm = T))

df17_db %>% 
  mutate(dw = sl.gl - sl.out) %>% 
  ggplot(aes(x = he, y = dw, fill = type)) +
  geom_col() +
  labs(x = "Hydrological event",
       y = expression(SL[GL]~-~SL[OUT]~", t")) +
  # ggrepel::geom_text_repel(aes(label = he),
  #                          point.padding = unit(0.5, "lines"),
  #                          segment.colour = "black", nudge_y = 10) +
  scale_fill_manual(name = "",
                    label = c("Anticlockwise",
                              "Clockwise",
                              "Figure-eight",
                              "Not clear"),
                    values = c( "#4DBBD5FF", "#E64B35FF",
                                "#00A087FF", "gray")) +
  scale_x_discrete(breaks = seq(0, 136, 5)) +
  scale_y_continuous(breaks = seq(-50, 100, 25)) +
  theme_clean() -> sed_budget

# Glacier role
df17_db %>% 
  mutate(delta = sl.gl * 100 / sl.out) %>%
  filter(delta <= 100) %>% 
  summarise(role = median(delta)) %>% 
  # summarise(role = density(delta)$x[which.max(density(delta)$y)]) %>% 
  pull(role) %>% cat("Glacier role is", ., "%")

df17_db %>% 
  mutate(delta = sl.gl * 100 / sl.out) %>%
  filter(delta <= 100) %>%
  ggplot(aes(x = delta)) +
  geom_segment(data = tibble(x = df17_db %>% 
                               mutate(delta = sl.gl * 100 / sl.out) %>%
                               filter(delta <= 100) %>%
                               pull(delta) %>% median,
                             y = 0.02),
               aes(x = x, xend = x,
                   y = 0, yend = y),
               linetype = "dashed") +
  annotate("text", x = 77.6, y = 0.02, label = "Median",
           vjust = -.5) +
  geom_density(color = "grey80",
               fill = "#0288b7",
               alpha = .35) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  labs(x = expression(italic(SL[GL])*", %"),
       y = "Probability density") +
  theme_clean() -> glacier_pdf

# SHI vs MID and GL
df17_db %>% 
  mutate(delta = 100 - sl.gl * 100 / sl.out) %>%
  filter(p == 0, !he %in% c(49, 56, 96)) %>%
  ggplot(aes(x = SHI, y = delta)) +
  geom_point(aes(color = type)) +
  # ggrepel::geom_text_repel(aes(label = he)) +
  ggpmisc::stat_poly_eq(aes(label =  paste(stat(rr.label))),
                        formula = y ~ x,
                        parse = T) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              color = "black",
              linetype = "dashed",
              se = F) +
  scale_x_continuous(limits = c(-0.2, 0.5)) +
  scale_color_manual(name = "",
                     label = c("Anticlockwise",
                               "Clockwise",
                               "Figure-eight",
                               "Not clear"),
                     values = c( "#4DBBD5FF", "#E64B35FF",
                                 "#00A087FF", "gray")) +
  labs(x = expression(italic(SHI)),
       y = expression(italic(SL[GL])*", %")) +
  theme_clean() -> shi_gl

df17_db %>% 
  mutate(delta = 100 - sl.mid * 100 / sl.out) %>%
  filter(p > 0, he != 89) %>%
  ggplot(aes(x = SHI, y = delta)) +
  geom_point(aes(color = type)) +
  # ggrepel::geom_text_repel(aes(label = he)) +
  ggpmisc::stat_poly_eq(aes(label =  paste(stat(rr.label))),
                        formula = y~x,
                        parse = T) +
  geom_smooth(method = "lm",
              color = "black",
              linetype = "dashed",
              se = F) +
  scale_x_continuous(limits = c(-0.2, 0.5)) +
  scale_color_manual(name = "",
                     label = c("Anticlockwise",
                               "Clockwise",
                               "Figure-eight",
                               "Not clear"),
                     values = c( "#4DBBD5FF", "#E64B35FF",
                                 "#00A087FF", "gray")) +  
  labs(x = expression(italic(SHI)),
       y = expression(italic(SL[MID])*", %")) +
  theme_clean() -> shi_mid

df17_db %>% 
  # mutate(delta = 100 - sl.mid * 100 / sl.out) %>%
  # filter(p > 0, he != 89) %>%
  mutate(delta = 100 - sl.gl * 100 / sl.out) %>%
  filter(p == 0, !he %in% c(49, 56, 96)) %>%
  select_if(is.numeric) %>% 
  correlate() %>% 
  focus(SHI)

# SAVE -------------------------------------------------------------------
ggsave("figures/fig08_sed-budget.png", sed_budget,
       dpi = 500, w = 10, h = 6)

ggsave("figures/fig09_glacier-pdf.png", glacier_pdf,
       dpi = 500, w = 7, h = 5)

ggpubr::ggarrange(shi_gl, shi_mid, nrow = 1, align = "hv",
                  labels = "AUTO", common.legend = T,
                  legend = "bottom") %>% 
  ggsave(plot = .,
         filename = "figures/fig10_shi-pct.png",
         dpi = 500, w = 10, h = 5)

