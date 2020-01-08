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

# Pairwise Wilcox test
pairwise.wilcox.test(x = df17_db$q.out.mean,
                     g = as_factor(df17_db$type),
                     p.adjust.method = "bonferroni")

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
  summarise(n = n(),
            med = median(delta),
            mean = mean(delta),
            pdf = density(delta)$x[which.max(density(delta)$y)]) 

df17_db %>% 
  mutate(delta = sl.gl * 100 / sl.out) %>%
  filter(delta <= 100) %>%
  ggplot(aes(x = delta)) +
  geom_segment(data = tibble(x = df17_db %>% 
                               mutate(delta = sl.gl * 100 / sl.out) %>%
                               filter(delta <= 100) %>%
                               pull(delta) %>% median,
                             y = 14),
               aes(x = x, xend = x,
                   y = 0, yend = y),
               linetype = "dashed") +
  annotate("text", x = 77.6, y = 15, label = "Median",
           vjust = .5) +
  geom_histogram(binwidth = 10,
                 color = "grey80",
                 fill = "#0288b7",
                 alpha = .35) +
  # geom_density(color = "grey80",
  #              fill = "#0288b7",
  #              alpha = .35) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  labs(x = expression(italic(SL[GL])*", %"),
       y = "Number of events") +
  theme_clean() -> glacier_pdf

# SHI vs MID and GL
df17_db %>% 
  mutate(delta = 100 * (sl.gl / sl.out)) %>%
  filter(p == 0, !he %in% c(49, 56, 96)) %>%
  ggplot(aes(x = SHI, y = delta, color = type)) +
  geom_point() +
  # ggrepel::geom_text_repel(aes(label = he)) +
  ggpmisc::stat_poly_eq(data = . %>% filter(type == "CW"),
                        aes(label =  paste(stat(rr.label))),
                        formula = y ~ x,
                        parse = T) +
  geom_smooth(data = . %>% filter(type == "CW"),
              method = "lm",
              formula = y ~ x,
              linetype = "dashed",
              se = F) +
  scale_x_continuous(limits = c(-0.2, 0.5)) +
  # scale_y_continuous(limits = c(-50, 80)) +
  scale_color_manual(name = "",
                     label = c("Anticlockwise",
                               "Clockwise",
                               "Figure-eight",
                               "Not clear"),
                     values = c( "#4DBBD5FF", "#E64B35FF",
                                 "#00A087FF", "gray")) +
  labs(title = "During non-rainfall events",
       subtitle = "a",
       x = expression(italic(SHI)),
       y = expression(italic(SL["GL,REL"])*", %")) +
  theme_clean() -> shi_gl_norain

df17_db %>% 
  mutate(delta = 100 * sl.gl / sl.out) %>%
  filter(p > 8) %>%
  ggplot(aes(x = SHI, y = delta, color = type)) +
  geom_point() +
  # ggrepel::geom_text_repel(aes(label = he)) +
  ggpmisc::stat_poly_eq(data = . %>% filter(type == "CW"),
                        aes(label =  paste(stat(rr.label))),
                        formula = y ~ x,
                        parse = T) +
  geom_smooth(data = . %>% filter(type == "CW"),
              method = "lm",
              formula = y ~ x,
              linetype = "dashed",
              se = F) +
  scale_x_continuous(limits = c(-0.2, 0.5)) +
  # scale_y_continuous(limits = c(-50, 80)) +
  scale_color_manual(name = "",
                     label = c("Anticlockwise",
                               "Clockwise",
                               "Figure-eight",
                               "Not clear"),
                     values = c( "#4DBBD5FF", "#E64B35FF",
                                 "#00A087FF", "gray")) +
  labs(title = "During rainfall events > 8 mm",
       subtitle = "c",
       x = expression(italic(SHI)),
       y = expression(italic(SL["GL,REL"])*", %")) +
  theme_clean() -> shi_gl_rain

df17_db %>% 
  mutate(delta = 100 * ((sl.mid - sl.gl)/sl.out)) %>%
  # mutate(delta = 100 - sl.mid * 100 / sl.out) %>%
  filter(p > 8) %>%
  ggplot(aes(x = SHI, y = delta, color = type)) +
  geom_point() +
  # ggrepel::geom_text_repel(aes(label = he)) +
  ggpmisc::stat_poly_eq(data = . %>% filter(type == "CW"),
                        aes(label =  paste(stat(rr.label))),
                        formula = y~x,
                        parse = T) +
  geom_smooth(data = . %>% filter(type == "CW"),
              method = "lm",
              linetype = "dashed",
              se = F) +
  scale_x_continuous(limits = c(-0.2, 0.5)) +
  # scale_y_continuous(limits = c(-50, 80)) +
  scale_color_manual(name = "",
                     label = c("Anticlockwise",
                               "Clockwise",
                               "Figure-eight",
                               "Not clear"),
                     values = c( "#4DBBD5FF", "#E64B35FF",
                                 "#00A087FF", "gray")) +  
  labs(subtitle = "d",
       x = expression(italic(SHI)),
       y = expression(italic(SL["MID,REL"])*", %")) +
  theme_clean() -> shi_mid_rain

df17_db %>% 
  mutate(delta = 100 * ((sl.mid - sl.gl)/sl.out)) %>%
  # mutate(delta = 100 - sl.mid * 100 / sl.out) %>%
  filter(p == 0) %>%
  ggplot(aes(x = SHI,
             y = delta,
             color = type)) +
  geom_point() +
  # ggrepel::geom_text_repel(aes(label = he)) +
  ggpmisc::stat_poly_eq(data = . %>% filter(type == "CW"),
                        aes(label =  paste(stat(rr.label))),
                        formula = y~x,
                        parse = T) +
  geom_smooth(data = . %>% filter(type == "CW"),
              method = "lm",
              linetype = "dashed",
              se = F) +
  scale_x_continuous(limits = c(-0.2, 0.5)) +
  # scale_y_continuous(limits = c(-50, 80)) +
  scale_color_manual(name = "",
                     label = c("Anticlockwise",
                               "Clockwise",
                               "Figure-eight",
                               "Not clear"),
                     values = c( "#4DBBD5FF", "#E64B35FF",
                                 "#00A087FF", "gray")) +  
  labs(subtitle = "b",
       x = expression(italic(SHI)),
       y = expression(italic(SL["MID,REL"])*", %")) +
  theme_clean() -> shi_mid_norain

df17_db %>% 
  mutate(delta = 100 * ((sl.mid - sl.gl)/sl.out)) %>%
  # mutate(delta = 100 - sl.mid * 100 / sl.out) %>%
  # filter(p > 8) %>% 
  select_if(is.numeric) %>% 
  correlate(method = "spearman")  %>% 
  focus(intensity) %>%
  arrange(intensity)

# SAVE -------------------------------------------------------------------
ggsave("figures/fig08_sed-budget.png", sed_budget,
       dpi = 500, w = 10, h = 6)

ggpubr::ggarrange(shi_gl_norain, shi_mid_norain,
                  shi_gl_rain, shi_mid_rain,
                  nrow = 2, ncol = 2,
                  align = "hv",
                  common.legend = T,
                  legend = "bottom") %>% 
  ggsave(plot = .,
         filename = "figures/fig09_shi-pct.png",
         dpi = 500, w = 6, h = 6)

# SAVE database as APPENDIX
df17_db %>%
  arrange(he) %>%
  select(-mean_date) %>% 
  mutate_at(vars(start, end), ~as.character(.)) %>% 
  mutate_all(~str_replace_na(.)) %>% 
  set_colnames(c("Hydrological event",
                 "SHI [—]",
                 "R.OUT [g/s]",
                 "R.MID [g/s]",
                 "R.GL [g/s]",
                 "SSC.OUT.mean [g/cbm]",
                 "SSC.OUT.max [g/cbm]",
                 "Q.OUT.mean [cbm/s]",
                 "Q.OUT.max [cbm/s]",
                 "Precipitation [mm]",
                 "Rainfall intensity [mm/h]",
                 "Mean Q lag time [h]",
                 "Start of the event [—]",
                 "End of the event[—]",
                 "Event duration [h]",
                 "SL.OUT [t/event]",
                 "SL.MID [t/event]",
                 "SL.GL [t/event]",
                 "Loop type")) %>% 
  # Reorder columns
  select(c(1, 13, 14, 15, 10, 11, 6, 7, 8, 9,
           3, 4, 5, 12, 16, 17, 18, 2, 19)) %>% 
  openxlsx::write.xlsx(x = .,
                       "analysis/APPENDIX.xlsx", )
