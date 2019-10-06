# EXTRAS

## Explore time lag between Q and SSC maximums
# 1) Select clockwise loops with negative lag time (?WTF)
df17_db %>%
  filter(type == "CW", q.lag < 0) %>%
  pull(he) %>%
  as.character() %>%
  as.numeric() -> what

# 2) Create a table 
df17_db %>% 
  mutate(rain = ifelse(p > 0, "yes", "no")) %>% 
  filter(!he %in%  what) %>% 
  mutate(type = case_when(type == "AW" ~ "Отрицательный",
                          type == "CW" ~ "Положительный",
                          type == "F8" ~ "Сложный",
                          TRUE ~ "Не определено")) %>% 
  group_by(type, rain) %>% 
  mutate(q.lag = abs(q.lag)) %>% 
  summarise(n = n(),
            delay_med = median(q.lag),
            delay_mean = mean(q.lag),
            delay_sd = sd(q.lag)) %>% 
  mutate(rain = ifelse(rain == "yes", "с осадками", "без осадков")) %>% 
  arrange(rain, desc(type)) %>% 
  set_colnames(c("Тип зависимости", "Тип события",
                 "Кол-во",
                 "Медиана", "Среднее", "SD")) %>% 
  openxlsx::write.xlsx(file = "analysis/tables_rus.xlsx")

# 3) Perform pairwise Wilcoxon Rank Sum test
df17_db %>% 
  mutate(rain = ifelse(p > 0, "yes", "no")) %>% 
  filter(!he %in%  what) %>%
  mutate(t = forcats::fct_cross(as_factor(type), as_factor(rain))) %>% 
  summarise(pairwise.wilcox.test(q.lag, t), )
  
pairwise.wilcox.test(abs(tt$q.lag), tt$t)



