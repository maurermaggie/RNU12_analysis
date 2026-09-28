library(readxl)
library(magrittr)
library(biomaRt)
library(scales)
library(ComplexUpset)
library(tidyverse)
require(data.table)

#############################----Sample ID labels----#########################################
# UDN112958 = E1, NR_029422.2:n.85G>A, NR_029422.2:n.96G>A
# UDN029422 = F1, NR_029422.2:​n.89A>G,NR_029422.2:​n.77T>A
# UDN907359 = G1, NR_029422.2:n.141G>C, NR_029422.2:n.134del
# UDN629307 = G1_father
# UDN664300 = G1_mother
#
# File paths below still use UDN IDs; sampleID values and plot labels use family codes.

SAMPLE_LABELS <- c(
  "UDN112958" = "E1",
  "UDN029422" = "F1",
  "UDN907359" = "G1",
  "UDN629307" = "G1_father",
  "UDN664300" = "G1_mother"
)

relabel_sampleID <- function(df) {
  df %>% mutate(sampleID = recode(sampleID, !!!SAMPLE_LABELS, .default = sampleID))
}

#############################----Combine results FRASER 1 Theta----#########################################
F1_MIRE <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/UDN029422/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_F1' = 'n') %>%
                    relabel_sampleID()

E1_MIRE <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/UDN112958/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_E1' = 'n') %>%
                    relabel_sampleID()

G1_MIRE <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/UDN907359/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_G1' = 'n') %>%
                    relabel_sampleID()

G1_father_MIRE <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/UDN629307/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_G1_father' = 'n') %>%
                    relabel_sampleID()

G1_mother_MIRE <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Theta/UDN664300/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_G1_mother' = 'n') %>%
                    relabel_sampleID()


MIRE <- full_join(F1_MIRE, E1_MIRE) %>% full_join(G1_father_MIRE) %>% full_join(G1_mother_MIRE) %>% full_join(G1_MIRE)

#############################----Combine results FRASER 1 Jaccard----#########################################
F1_J <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/UDN029422/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_F1' = 'n') %>%
                    relabel_sampleID()

E1_J <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/UDN112958/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_E1' = 'n') %>%
                    relabel_sampleID()

G1_father_J <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/UDN629307/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_G1_father' = 'n') %>%
                    relabel_sampleID()

G1_mother_J <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/UDN664300/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_G1_mother' = 'n') %>%
                    relabel_sampleID()

G1_J <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER2/UDN907359/minor/junction_level/juncs_w_minor_intron_retention/unique_junction_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_G1' = 'n') %>%
                    relabel_sampleID()

Jaccard_FRASER2 <- full_join(F1_J, E1_J) %>% full_join(G1_father_J) %>% full_join(G1_mother_J) %>% full_join(G1_J)

#############################----G1's psi3 and psi5 results----#########################################
G1_psi3 <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi3/UDN907359/minor/junction_level/juncs_w_intron_retention_in_genes_w_minor_intron/unique_gene_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_G1' = 'n') %>%
                    relabel_sampleID()

G1_psi5 <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/compiled_results/FRASER1/Psi5/UDN907359/minor/junction_level/juncs_w_intron_retention_in_genes_w_minor_intron/unique_gene_counts.csv") %>% 
                    select(-pairedEnd, -Z_score, -bamFile) %>%
                    rename('n_G1' = 'n') %>%
                    relabel_sampleID()
                    

#############################----Graph FRASER1 theta results----#########################################
x_lab <- expression(atop(""))
y_lab <- expression("Number of MIREs")

upper_2sd_F1 <- mean(MIRE$n_F1, na.rm = TRUE) + 2 * sd(MIRE$n_F1, na.rm = TRUE)
upper_1sd_F1 <- mean(MIRE$n_F1, na.rm = TRUE) + sd(MIRE$n_F1, na.rm = TRUE)

max_row_F1 <- MIRE %>%
  filter(sampleID == "F1")

F1_plt <- ggplot(F1_MIRE, aes(x = reorder(sampleID, n_F1), y = n_F1)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#7D5E62") +
  xlab(x_lab) +
  ylab(y_lab) +
  geom_point(
    data = max_row_F1,
    color = "red",
    size = 4
  ) +
  geom_hline(data=MIRE, aes(yintercept=upper_2sd_F1, color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = MIRE, aes(label = "Mean + 2*SD", x = 25, upper_2sd_F1+.5),  size=8, hjust = -1, vjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=-1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
F1_plt

upper_2sd_E1 <- mean(MIRE$n_E1, na.rm = TRUE) + 2 * sd(MIRE$n_E1, na.rm = TRUE)
upper_1sd_E1 <- mean(MIRE$n_E1, na.rm = TRUE) + sd(MIRE$n_E1, na.rm = TRUE)

max_row_E1 <- MIRE %>%
  filter(sampleID == "E1")

E1_plt <- ggplot(MIRE, aes(x = reorder(sampleID, n_E1), y = n_E1)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#7D5E62") +
  xlab(x_lab) +
  ylab(y_lab) +
  geom_point(
    data = max_row_E1,
    color = "red",
    size = 4
  ) +
  geom_hline(data=MIRE, aes(yintercept=upper_2sd_E1, color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = MIRE, aes(label = "Mean + 2*SD", x = 25, upper_2sd_E1+.5),  size=8, hjust = -1, vjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=-1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
E1_plt

upper_2sd_G1 <- mean(MIRE$n_G1, na.rm = TRUE) + 2 * sd(MIRE$n_G1, na.rm = TRUE)
upper_1sd_G1 <- mean(MIRE$n_G1, na.rm = TRUE) + sd(MIRE$n_G1, na.rm = TRUE)

max_row_G1 <- MIRE %>%
  filter(sampleID == "G1")

G1_plt <- ggplot(MIRE, aes(x = reorder(sampleID, n_G1), y = n_G1)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#7D5E62") +
  xlab(x_lab) +
  ylab(y_lab) +
  geom_point(
    data = max_row_G1,
    color = "red",
    size = 4
  ) +
  geom_hline(data=MIRE, aes(yintercept=upper_2sd_G1, color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = MIRE, aes(label = "Mean + 2*SD", x = 25, upper_2sd_G1+.05), hjust=-1, size=8, vjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

G1_plt

upper_2sd_G1_mother <- mean(MIRE$n_G1_mother, na.rm = TRUE) + 2 * sd(MIRE$n_G1_mother, na.rm = TRUE)
upper_1sd_G1_mother <- mean(MIRE$n_G1_mother, na.rm = TRUE) + sd(MIRE$n_G1_mother, na.rm = TRUE)

max_row_G1_mother <- MIRE %>%
  filter(sampleID == "G1_mother")

G1_mother_plt <- ggplot(MIRE, aes(x = reorder(sampleID, n_G1_mother), y = n_G1_mother)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#7D5E62") +
  xlab(x_lab) +
  ylab(y_lab) +
  geom_point(
    data = max_row_G1_mother,
    color = "red",
    size = 4
  ) +
  geom_hline(data=MIRE, aes(yintercept=upper_2sd_G1_mother, color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = MIRE, aes(label = "Mean + 2*SD", x = 25, upper_2sd_G1_mother+.02), hjust=-1, size=8, vjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

G1_mother_plt

upper_2sd_G1_father <- mean(MIRE$n_G1_father, na.rm = TRUE) + 2 * sd(MIRE$n_G1_father, na.rm = TRUE)
upper_1sd_G1_father <- mean(MIRE$n_G1_father, na.rm = TRUE) + sd(MIRE$n_G1_father, na.rm = TRUE)

max_row_G1_father <- MIRE %>%
  filter(sampleID == "G1_father")

G1_father_plt <- ggplot(MIRE, aes(x = reorder(sampleID, n_G1_father), y = n_G1_father)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#7D5E62") +
  xlab(x_lab) +
  ylab(y_lab) +
  geom_point(
    data = max_row_G1_father,
    color = "red",
    size = 4
  ) +
  geom_hline(data=MIRE, aes(yintercept=upper_2sd_G1_father, color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = MIRE, aes(label = "Mean + 2*SD", x = 25, upper_2sd_G1_father+.05), hjust=-1, size=8, vjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

G1_father_plt

row1_title <- ggplot() +
  theme_void() +
  annotate(
    "text", x = 0.5, y = 0.5,
    label = "RNU12 Probands",
    size = 9,
    fontface = "bold"
  ) +
  xlim(0, 1) +
  ylim(0, 1)

row2_title <- ggplot() +
  theme_void() +
  annotate(
    "text", x = 0.5, y = 0.5,
    label = "Unaffected Family Members",
    size = 9,
    fontface = "bold"
  ) +
  xlim(0, 1) +
  ylim(0, 1)

layout <- "
AAAAAA
BBCCDD
EEEEEE
#FFGG#
"

F1_plt <- F1_plt + labs(title="F1")
E1_plt <- E1_plt + labs(title="E1")
G1_plt <- G1_plt + labs(title="G1")
G1_mother_plt <- G1_mother_plt + labs(title="G1_mother")
G1_father_plt <- G1_father_plt + labs(title="G1_father")

final_total_outliers <- patchwork::wrap_plots(
  A = row1_title,
  B = E1_plt,
  C = F1_plt,
  D = G1_plt,
  E = row2_title,
  F = G1_mother_plt,
  G = G1_father_plt,
  design = layout,
  heights = c(0.12, 1, 0.12, 1)
)
final_total_outliers

ggsave("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/results_analysis/graphs/final_total_outliers_labeled.png", final_total_outliers, width = 22, height = 10, dpi = 300)

#############################----Graph FRASER2 jaccard results----#########################################
x_lab <- expression(atop(""))
y_lab <- expression("Number of Jaccard MI outliers")

upper_2sd_F1 <- mean(Jaccard_FRASER2$n_F1, na.rm = TRUE) + 2 * sd(Jaccard_FRASER2$n_F1, na.rm = TRUE)
upper_1sd_F1 <- mean(Jaccard_FRASER2$n_F1, na.rm = TRUE) + sd(Jaccard_FRASER2$n_F1, na.rm = TRUE)

max_row_F1 <- Jaccard_FRASER2 %>%
  filter(sampleID == "F1")

F1_plt_J <- ggplot(F1_J, aes(x = reorder(sampleID, n_F1), y = n_F1)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#7D5E62") +
  xlab(x_lab) +
  ylab(y_lab) +
  geom_point(
    data = max_row_F1,
    color = "red",
    size = 4
  ) +
  geom_hline(data=F1_J, aes(yintercept=upper_2sd_F1, color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = F1_J, aes(label = "Mean + 2*SD", x = 25, upper_2sd_F1+.5),  size=8, hjust = -1, vjust=0.000025) +
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=-1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
F1_plt_J

upper_2sd_E1 <- mean(Jaccard_FRASER2$n_E1, na.rm = TRUE) + 2 * sd(Jaccard_FRASER2$n_E1, na.rm = TRUE)
upper_1sd_E1 <- mean(Jaccard_FRASER2$n_E1, na.rm = TRUE) + sd(Jaccard_FRASER2$n_E1, na.rm = TRUE)

max_row_E1 <- Jaccard_FRASER2 %>%
  filter(sampleID == "E1")

E1_plt_J <- ggplot(Jaccard_FRASER2, aes(x = reorder(sampleID, n_E1), y = n_E1)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#7D5E62") +
  xlab(x_lab) +
  ylab(y_lab) +
  geom_point(
    data = max_row_E1,
    color = "red",
    size = 4
  ) +
  geom_hline(data=Jaccard_FRASER2, aes(yintercept=upper_2sd_E1, color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = Jaccard_FRASER2, aes(label = "Mean + 2*SD", x = 25, upper_2sd_E1+.5),  size=8, hjust = -1, vjust=0.000025) +
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=-1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
E1_plt_J

upper_2sd_G1 <- mean(Jaccard_FRASER2$n_G1, na.rm = TRUE) + 2 * sd(Jaccard_FRASER2$n_G1, na.rm = TRUE)
upper_1sd_G1 <- mean(Jaccard_FRASER2$n_G1, na.rm = TRUE) + sd(Jaccard_FRASER2$n_G1, na.rm = TRUE)

max_row_G1 <- Jaccard_FRASER2 %>%
  filter(sampleID == "G1")

G1_plt_J <- ggplot(Jaccard_FRASER2, aes(x = reorder(sampleID, n_G1), y = n_G1)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#7D5E62") +
  xlab(x_lab) +
  ylab(y_lab) +
  geom_point(
    data = max_row_G1,
    color = "red",
    size = 4
  ) +
  geom_hline(data=Jaccard_FRASER2, aes(yintercept=upper_2sd_G1, color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = Jaccard_FRASER2, aes(label = "Mean + 2*SD", x = 25, upper_2sd_G1+.05), hjust=-1, size=8, vjust=0.000025) +
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

G1_plt_J

ggsave("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/results_analysis/G1_plt_J.png", G1_plt_J, dpi = 300)


upper_2sd_G1_mother <- mean(Jaccard_FRASER2$n_G1_mother, na.rm = TRUE) + 2 * sd(Jaccard_FRASER2$n_G1_mother, na.rm = TRUE)
upper_1sd_G1_mother <- mean(Jaccard_FRASER2$n_G1_mother, na.rm = TRUE) + sd(Jaccard_FRASER2$n_G1_mother, na.rm = TRUE)

max_row_G1_mother <- Jaccard_FRASER2 %>%
  filter(sampleID == "G1_mother")

G1_mother_plt_J <- ggplot(Jaccard_FRASER2, aes(x = reorder(sampleID, n_G1_mother), y = n_G1_mother)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#7D5E62") +
  xlab(x_lab) +
  ylab(y_lab) +
  geom_point(
    data = max_row_G1_mother,
    color = "red",
    size = 4
  ) +
  geom_hline(data=Jaccard_FRASER2, aes(yintercept=upper_2sd_G1_mother, color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = Jaccard_FRASER2, aes(label = "Mean + 2*SD", x = 25, upper_2sd_G1_mother+.02), hjust=-1, size=8, vjust=0.000025) +
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

G1_mother_plt_J

upper_2sd_G1_father <- mean(Jaccard_FRASER2$n_G1_father, na.rm = TRUE) + 2 * sd(Jaccard_FRASER2$n_G1_father, na.rm = TRUE)
upper_1sd_G1_father <- mean(Jaccard_FRASER2$n_G1_father, na.rm = TRUE) + sd(Jaccard_FRASER2$n_G1_father, na.rm = TRUE)

max_row_G1_father <- Jaccard_FRASER2 %>%
  filter(sampleID == "G1_father")

G1_father_plt_J <- ggplot(Jaccard_FRASER2, aes(x = reorder(sampleID, n_G1_father), y = n_G1_father)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#7D5E62") +
  xlab(x_lab) +
  ylab(y_lab) +
  geom_point(
    data = max_row_G1_father,
    color = "red",
    size = 4
  ) +
  geom_hline(data=Jaccard_FRASER2, aes(yintercept=upper_2sd_G1_father, color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = Jaccard_FRASER2, aes(label = "Mean + 2*SD", x = 25, upper_2sd_G1_father+.05), hjust=-1, size=8, vjust=0.000025) +
  theme_classic(base_size = 20)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

G1_father_plt_J

row1_title <- ggplot() +
  theme_void() +
  annotate(
    "text", x = 0.5, y = 0.5,
    label = "RNU12 Probands",
    size = 9,
    fontface = "bold"
  ) +
  xlim(0, 1) +
  ylim(0, 1)

row2_title <- ggplot() +
  theme_void() +
  annotate(
    "text", x = 0.5, y = 0.5,
    label = "Unaffected Family Members",
    size = 9,
    fontface = "bold"
  ) +
  xlim(0, 1) +
  ylim(0, 1)

layout <- "
AAAAAA
BBCCDD
EEEEEE
#FFGG#
"

F1_plt_J <- F1_plt_J + labs(title="F1")
E1_plt_J <- E1_plt_J + labs(title="E1")
G1_plt_J <- G1_plt_J + labs(title="G1")
G1_mother_plt_J <- G1_mother_plt_J + labs(title="G1_mother")
G1_father_plt_J <- G1_father_plt_J + labs(title="G1_father")

final_total_outliers_J <- patchwork::wrap_plots(
  A = row1_title,
  B = E1_plt_J,
  C = F1_plt_J,
  D = G1_plt_J,
  E = row2_title,
  F = G1_mother_plt_J,
  G = G1_father_plt_J,
  design = layout,
  heights = c(0.12, 1, 0.12, 1)
)
final_total_outliers_J

ggsave("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/results_analysis/graphs/final_jaccard_outliers_labeled.png", final_total_outliers_J, width = 22, height = 10, dpi = 300)

#############################----Graph COMBINED FRASER1 theta results----#########################################
MIRE <- MIRE %>% mutate(avg_MIRE = rowMeans(across(n_F1:n_G1), na.rm = TRUE))

MIRE_mean <- mean(MIRE$avg_MIRE)
MIRE_sd <- sd(MIRE$avg_MIRE)
MIRE$Z_score <- (MIRE$avg_MIRE - MIRE_mean)/MIRE_sd

max_row <- MIRE %>%
  filter(sampleID %in% c("E1", "F1", "G1"))

upper_2sd <- mean(MIRE$avg_MIRE, na.rm = TRUE) + 2 * sd(MIRE$avg_MIRE, na.rm = TRUE)
upper_1sd <- mean(MIRE$avg_MIRE, na.rm = TRUE) + sd(MIRE$avg_MIRE, na.rm = TRUE)

MIRE_plt <- ggplot(MIRE, aes(x = reorder(sampleID, avg_MIRE), y = avg_MIRE)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#7D5E62") +
  xlab(x_lab) +
  ylab(y_lab) +
  geom_point(
    data = max_row,
    color = "red",
    size = 4
  ) +
  geom_hline(data=MIRE, aes(yintercept=upper_2sd, color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = MIRE, aes(label = "Mean + 2*SD", x = 25, upper_2sd+.05), hjust=-1, size=8, vjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#############################----Graph COMBINED FRASER1 theta results----#########################################
Jaccard_FRASER2 <- Jaccard_FRASER2 %>% mutate(avg_Jaccard = rowMeans(across(n_F1:n_G1), na.rm = TRUE))

J_mean <- mean(Jaccard_FRASER2$avg_Jaccard)
J_sd <- sd(Jaccard_FRASER2$avg_Jaccard)
Jaccard_FRASER2$Z_score <- (Jaccard_FRASER2$avg_Jaccard - J_mean)/J_sd

max_row <- Jaccard_FRASER2 %>%
  filter(sampleID %in% c("E1", "F1", "G1"))

upper_2sd <- mean(Jaccard_FRASER2$avg_Jaccard, na.rm = TRUE) + 2 * sd(Jaccard_FRASER2$avg_Jaccard, na.rm = TRUE)
upper_1sd <- mean(Jaccard_FRASER2$avg_Jaccard, na.rm = TRUE) + sd(Jaccard_FRASER2$avg_Jaccard, na.rm = TRUE)

x_lab <- expression("Samples ordered by the number of Jaccard MI outliers")
y_lab <- expression("Number of Jaccard MI outliers")

Jaccard_plt <- ggplot(Jaccard_FRASER2, aes(x = reorder(sampleID, avg_Jaccard), y = avg_Jaccard)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#7D5E62") +
  xlab(x_lab) +
  ylab(y_lab) +
  labs(title="Number of Jaccard index outliers including minor introns from averaged pool")+
  geom_point(
    data = max_row,
    color = "red",
    size = 4
  ) +
  geom_hline(data=Jaccard_FRASER2, aes(yintercept=upper_2sd, color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = Jaccard_FRASER2, aes(label = "Mean + 2*SD", x = 25, upper_2sd+.2), hjust=1, size=8, vjust=-1) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

Jaccard_plt
ggsave("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/results_analysis/Jaccard_plt.png", Jaccard_plt, width = 22, height = 10, dpi = 300)

#############################----Graph psi3 results----#########################################
upper_2sd_psi3 <- mean(G1_psi3$n_G1, na.rm = TRUE) + 2 * sd(G1_psi3$n_G1, na.rm = TRUE)
upper_1sd_psi3 <- mean(G1_psi3$n_G1, na.rm = TRUE) + sd(G1_psi3$n_G1, na.rm = TRUE)

max_row_G1_psi3 <- G1_psi3 %>%
  filter(sampleID == "G1")

G1_plt_psi3 <- ggplot(G1_psi3, aes(x = reorder(sampleID, n_G1), y = n_G1)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#7D5E62") +
  ylab("Number of Psi3 outliers in MIGs") +
  xlab("Samples ordered by the \nnumber of Psi3 outliers in MIGs")+
  labs(title="Psi3")+
  geom_point(
    data = max_row_G1_psi3,
    color = "red",
    size = 4
  ) +
  geom_hline(data=G1_psi3, aes(yintercept=upper_2sd_psi3, color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = G1_psi3, aes(label = "Mean + 2*SD", x = 25, upper_2sd_psi3+.05), hjust=-1, size=8, vjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

G1_plt_psi3
ggsave("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/results_analysis/psi3_plt.png", G1_plt_psi3, dpi = 300)

#############################----Graph psi5 results----#########################################
upper_2sd_psi5 <- mean(G1_psi5$n_G1, na.rm = TRUE) + 2 * sd(G1_psi5$n_G1, na.rm = TRUE)
upper_1sd_psi5 <- mean(G1_psi5$n_G1, na.rm = TRUE) + sd(G1_psi5$n_G1, na.rm = TRUE)

max_row_G1_psi5 <- G1_psi5 %>%
  filter(sampleID == "G1")

G1_plt_psi5 <- ggplot(G1_psi5, aes(x = reorder(sampleID, n_G1), y = n_G1)) +
  geom_point(shape=21, alpha=.6, size=5, fill="#7D5E62") +
  ylab("Number of Psi5 outliers in MIGs") +
  xlab("Samples ordered by the \nnumber of Psi5 outliers in MIGs")+
  labs(title="Psi5")+
  geom_point(
    data = max_row_G1_psi5,
    color = "red",
    size = 4
  ) +
  geom_hline(data=G1_psi5, aes(yintercept=upper_2sd_psi5, color = group),
              color="blue", linetype="dashed", size=0.8) +
  geom_text(data = G1_psi5, aes(label = "Mean + 2*SD", x = 25, upper_2sd_psi5+.05), hjust=-1, size=8, vjust=0.000025) +
  theme_classic(base_size = 25)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

G1_plt_psi5
ggsave("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/results_analysis/psi5_plt.png", G1_plt_psi5, dpi = 300)

#############################----Compare genes with jaccard outliers including minor introns----#########################################
MIDB_all_fp <- "/oak/stanford/groups/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/Phenotype_Paper/data/Homo_sapiens_intron_new.csv" 
MIDB_all <- read_csv(MIDB_all_fp) %>% filter(intron_class == "minor")

G1 <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/UDN907359_nygc/FRASER2/FRASER2_output_filtered_p_lt_051.csv", col_names = FALSE)
colnames(G1) <- c("seqnames", "start", "end", "width", "strand", "sampleID", "hgncSymbol", "type", "pValue", "padjust", 
                            "psiValue", "deltaPsi", "meanCounts", "meanTotalCounts", "counts", "totalCounts", "nonsplitCounts",
                            "nonsplitProportion", "nonsplitProportion_99quantile")
G1 <- G1 %>% filter(sampleID == "UDN907359") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3)

E1 <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/UDN112958_nygc/FRASER2/FRASER2_output_filtered_p_lt_051.csv", col_names = FALSE)
colnames(E1) <- c("seqnames", "start", "end", "width", "strand", "sampleID", "hgncSymbol", "type", "pValue", "padjust", 
                            "psiValue", "deltaPsi", "meanCounts", "meanTotalCounts", "counts", "totalCounts", "nonsplitCounts",
                            "nonsplitProportion", "nonsplitProportion_99quantile")
E1 <- E1 %>% filter(sampleID == "UDN112958") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3)

F1 <- read_csv("/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/output/FRASER/UDN029422_non_nygc/FRASER2/FRASER2_output_filtered_p_lt_051.csv", col_names = FALSE)
colnames(F1) <- c("seqnames", "start", "end", "width", "strand", "sampleID", "hgncSymbol", "type", "pValue", "padjust", 
                            "psiValue", "deltaPsi", "meanCounts", "meanTotalCounts", "counts", "totalCounts", "nonsplitCounts",
                            "nonsplitProportion", "nonsplitProportion_99quantile")
F1 <- F1 %>% filter(sampleID == "UDN029422") %>% filter(padjust < 0.05) %>% filter(abs(deltaPsi) > 0.3)

get_genes_with_MIR <- function(data_frame) {
    pro_ensg <- mean(startsWith(data_frame$hgncSymbol, "ENSG"), na.rm = TRUE)

    if (pro_ensg > 0.5) {
      intron_list <- MIDB_all %>% pull(ensembl_gene_id) %>% unique()
      MIDB_for_join <- MIDB_all %>% rename(hgncSymbol = ensembl_gene_id)
      F1_theta_mig <- data_frame %>% filter(hgncSymbol %in% intron_list)
      MIDB_joined <- left_join(MIDB_for_join, F1_theta_mig, by = "hgncSymbol")
    } else {
      intron_list <- MIDB_all %>% pull(gene_symbol) %>% unique()
      F1_theta_mig <- data_frame %>%
        filter(hgncSymbol %in% intron_list) %>%
        rename(gene_symbol = hgncSymbol)
      MIDB_joined <- left_join(MIDB_all, F1_theta_mig, by = "gene_symbol")
    }

    theta_minor_junctions <- MIDB_joined %>%
      filter(intron_class == "minor") %>%
      filter(start == intron_start | start == intron_start - 1 | start == intron_start + 1 |
            end == intron_end | end == intron_end - 1 | end == intron_end + 1) %>%
      select(any_of(c("sampleID", "gene_symbol", "intron_name", "intron_start", "intron_end",
                      "start", "end", "intron_class", "pValue", "padjust", "zScore", "psiValue",
                      "deltaPsi", "meanCounts", "meanTotalCounts", "counts", "totalCounts"))) %>%
      unique()

    theta_minor_junctions %>% pull(gene_symbol) %>% unique()
}

F1_genes <-get_genes_with_MIR(F1)
E1_genes <-get_genes_with_MIR(E1)
G1_genes <-get_genes_with_MIR(G1)

gene_list <- c(G1_genes, F1_genes, E1_genes) %>% unique 
genes <- data.frame(gene_list)
genes$G1 <- ifelse(genes$gene_list %in% G1_genes, TRUE, FALSE)
genes$F1 <- ifelse(genes$gene_list %in% F1_genes, TRUE, FALSE)
genes$E1 <- ifelse(genes$gene_list %in% E1_genes, TRUE, FALSE)
rownames(genes) <- genes$gene_list

y_lab <- expression(" ")

upset <- upset(genes, colnames(genes)[2:4], name='gene_list', width_ratio=0.1,
    matrix=(
        intersection_matrix(geom=geom_point(shape='circle filled', size=3))
        + scale_color_manual(
            values=c('G1'='#DA786C', 'F1'='#88ACAA', 'E1'='#F2D177'),
            guide=guide_legend(override.aes=list(shape='circle'))
        )
    ),
    queries=list(
        upset_query(set='G1', fill='#DA786C'),
        upset_query(set='F1', fill='#88ACAA'),
        upset_query(set='E1', fill='#F2D177')
    )
) +
theme_classic(base_size=30)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(plot.margin = margin(1,1,1.5,1.5, "cm"))+
  xlab("Number of Shared Jaccard outliers including MIs")+
  ylab(y_lab)

upset
ggsave(filename= "/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis/FRASER/results_analysis/upset.png", plot=upset, limitsize = FALSE, units = "in", height=15, width=15)
