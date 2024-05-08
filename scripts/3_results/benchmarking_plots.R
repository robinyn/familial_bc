library(tidyverse)
library(ggplot2)

benchmark = read_tsv("~/Desktop/RESULTS/benchmark_forR.tsv")
all_runs = read_tsv("~/Desktop/RESULTS/benchmark_all_runs.tsv")

benchmark = benchmark %>% 
  pivot_longer(cols=c(run_time, cpu_time), names_to = "type", values_to = "time") %>% 
  mutate(type = str_remove(type, "_time"))

all_runs = all_runs %>% 
  pivot_longer(cols=c(run_time, cpu_time), names_to = "type", values_to = "time") %>% 
  mutate(type = str_remove(type, "_time")) 

df = all_runs %>% 
  group_by(pipeline, type) %>% 
  summarise(mean=mean(time), sd=sd(time))

plot_theme = theme(panel.background = element_blank(),
                   axis.line = element_line(color = "black", linewidth = 1),
                   axis.ticks = element_line(linewidth = 1),
                   axis.ticks.length = unit(0.2, "cm"),
                   axis.text.x = element_text(size = 18, vjust = 0.5),
                   axis.text.y = element_text(size = 18, hjust = 0.5),
                   axis.title = element_text(size = 20))


ggplot(df %>% 
         mutate(pipeline=str_to_sentence(pipeline)) %>% 
         mutate(pipeline=factor(pipeline, levels=c("Unoptimized", "Optimized"))) %>% 
         mutate(type=str_replace(type, "cpu", "CPU time")) %>% 
         mutate(type=str_replace(type, "run", "Run time")), 
       aes(x=type, y=mean, fill=pipeline)) +
  geom_bar(position="dodge", stat="identity") +
  plot_theme +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size=18),
        legend.box.background = element_blank(),
        legend.position = c(0.8, 0.9)) +
  scale_fill_manual(values=c("black", "red")) +
  labs(x="", y="Mean time (s)")


ggplot(benchmark %>% 
         filter(pipeline=="optimized") %>% 
         mutate(type=str_replace(type, "cpu", "CPU time")) %>% 
         mutate(type=str_replace(type, "run", "Run time")), aes(x=no_cores, y=time, shape=type, color=type)) +
  geom_point(size=5) +
  plot_theme +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0.17, 0.1),
        legend.text = element_text(size=18),
        legend.box.background = element_blank()) +
  scale_color_manual(values=c("black", "red")) + 
  labs(x="Number of cores", y="Mean time (s)")
