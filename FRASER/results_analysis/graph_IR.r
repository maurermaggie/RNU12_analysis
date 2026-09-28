library(readxl)
library(magrittr)
library(biomaRt)
library(tidyverse)
require(data.table)

UDN029422_fp <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/UDN029422_non_nygc/FRASER1/FRASER_output.csv"
UDN029422 <- read_csv(UDN029422_fp)

UDN029422 <- UDN029422 %>% filter(type == "theta") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3) %>%
        group_by(sampleID) %>% tally %>% arrange(desc(n))

colnames(UDN029422) <- c("sampleID", "n")

UDN029422$numbers <- rownames(UDN029422)
title_lab <- expression("Intron Retention Outliers in UDN Subset")
x_lab <- expression(atop("Samples ordered by Number of Junctions with",
                     "Significant Intron Retention Outliers"))
y_lab <- expression("Number of Junctions Significant Intron Retention Outliers")

theta_plot <- ggplot(UDN029422, aes(x = reorder(sampleID, n), y = n)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#7D5E62") +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title=title_lab)  +
  #geom_point(data=g1, aes(colour= numbers)) +
  labs(colour="Sample ID") +
  #geom_text_repel(data=(joined %>% filter(numbers %in% top_10)), size=3, aes(x = reorder(numbers, no_MIGs_theta_juncs), y = no_MIGs_theta_juncs, label = numbers), position=position_dodge(width=0.9),
  #                  hjust=1.5,colour="black", max.overlaps = Inf) +
  geom_hline(data=UDN029422, aes(yintercept=mean(n)+sd(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_hline(data=UDN029422, aes(yintercept=mean(n)+sd(n), color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = UDN029422, aes(label = "Mean + 2*SD", x = 25, y=mean(UDN029422$n) + sd(UDN029422$n) + sd(UDN029422$n) + 3),  size=8, hjust = 1, vjust=0.000025) +
  geom_text(data = UDN029422, aes(label = "Mean + 1*SD", x = 25, y=mean(UDN029422$n) + sd(UDN029422$n) + 3),  size=8, hjust = 1, vjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(filename="/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/results_analysis/theta_plot.pdf", plot=theta_plot,  limitsize = FALSE, units = "in", height=12, width=10)
