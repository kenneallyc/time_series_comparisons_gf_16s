setwd("~/Desktop/GF_16s_1")
getwd()
path <- "/Users/charlotte.kenneally/Desktop/GF_16s_1"

## needs
require(grid)
library(phyloseq)
library(ggplot2)
library(data.table)
library(tidyverse)

otu <-read.csv("OTU600_Jan8_no1Csamples.csv", header = TRUE, row.names = 1)
map <-read.csv("samdat600_Jan8_no1Csamples.csv", header = TRUE, row.names = 1)
tax <-read.csv("tax600_GF_Dec1_23.csv", header = TRUE, row.names = colnames(otu))
tax <- as.matrix(tax[2:8])

### useful info -----
"#290AD8FF","#D5D5FFFF","#72D9FFFF","#AAF7FFFF","#1E8E99FF", "#FFEE99FF",
"#FF2B00FF","#D9F0D3FF","#FFAD72FF","#F76D5EFF","#D82632FF","#FF7080FF"

"#5991FFFF","#8CB2FFFF","#BFD4FFFF","#E6EEFFFF","#F7FAFFFF","#FFFFCCFF","#FFFF99FF","#FFFF00FF","#FFCC00FF","#FF9900FF","#FF6600FF","#FF0000FF"))+
  "#8E0152FF","#C51B7DFF","#DE77AEFF","#F1B6DAFF","#FDE0EFFF","#F7F7F7FF","#E6F5D0FF","#B8E186FF","#7FBC41FF","#4D9221FF","#276419FF","#417839FF"))+ 
  "#A7D3D4FF", "#F8D564FF", "#EDD6D5FF", "#EA879CFF", "#FF3D7FFF", "#7FC7AFFF", "#E9E29CFF",
  "#E0E4CCFF", "#F38630FF", "#FA6900FF", "#BF9BDDFF"



## phyloseq objects ------
ps3 <- phyloseq(otu_table(otu, taxa_are_rows = FALSE), 
                tax_table(tax), sample_data(map)) 
rowSums(otu_table(ps3)) %>%
  sort()
ps3_cut3000 <- rarefy_even_depth(ps3, rngseed = 777, sample.size = 3000, replace = FALSE)
rowSums(otu_table(ps3_cut3000)) %>%
  sort()

## separate and make phyloseq objects for each round
round1 <- subset_samples(ps3_cut3000, round == "1")
round2 <- subset_samples(ps3_cut3000, round == "2")
## those still include all samples from those rounds 
df_round1 <- data.table(psmelt(round1))
df_round2 <- data.table(psmelt(round2))

### phyloseq objects for each donor
pdx_rec <- subset_samples(ps3_cut3000, donor == "Cr")
mr1_rec <- subset_samples(ps3_cut3000, donor == "mr1ko")
bl6_rec <- subset_samples(ps3_cut3000, donor == "BL6-1")
pooled1_rec <- subset_samples(ps3_cut3000, donor == "Pooled1")

## data frames for each
df_pdx_rec <- data.table(psmelt(pdx_rec))

## week 1 and week 8 separated data frames to plot separate ------
## specify to week 1 and then to week 8
round1_wk1 <- subset_samples(round1, days_post_fmt == "7")
round1_wk8 <- subset_samples(round1, days_post_fmt == "56")
## remove cecum and colon samples from week 8 
round1_wk8 <- subset_samples(round1_wk8, sample_type == "stool")
## make data frame of those to check 
df_round1_wk1 <- data.table(psmelt(round1_wk1))
df_round1_wk8 <- data.table(psmelt(round1_wk8))

## now break it up by cage 
#round1_wk1_cageA <- subset_samples(round1_wk1, donor == "Cr")
#round1_wk1_cageB <- subset_samples(round1_wk1, donor == "mr1ko")
#round1_wk1_cageD <- subset_samples(round1_wk1, donor == "Pooled1")
## turn into data frames
#df_round1_wk1_cageA <- data.table(psmelt(round1_wk1_cageA))
#df_round1_wk1_cageB <- data.table(psmelt(round1_wk1_cageB))
#df_round1_wk1_cageD <- data.table(psmelt(round1_wk1_cageD))

## week 8 turn
#round1_wk8_cageA <- subset_samples(round1_wk8, donor == "Cr")
#round1_wk8_cageB <- subset_samples(round1_wk8, donor == "mr1ko")
#round1_wk8_cageD <- subset_samples(round1_wk8, donor == "Pooled1")
## turn into data frames
#df_round1_wk8_cageA <- data.table(psmelt(round1_wk8_cageA))
#df_round1_wk8_cageB <- data.table(psmelt(round1_wk8_cageB))
#df_round1_wk8_cageD <- data.table(psmelt(round1_wk8_cageD))




### filter data frames for days and sample types I want ------
## by day
df_rd1_wk1_8_allsamples <- filter(df_round1, days_post_fmt == "7" | days_post_fmt == "56")
### filter out the colon and cecum samples
df_rd1_wk1_8 <- filter(df_rd1_wk1_8_allsamples, sample_type == "stool")
### separate by cage 
df_rd1_wk1_8_cageA <- filter(df_rd1_wk1_8, donor == "Cr")
df_rd1_wk1_8_cageB <- filter(df_rd1_wk1_8, donor == "mr1ko")
df_rd1_wk1_8_cageD <- filter(df_rd1_wk1_8, donor == "Pooled1")

## now round 2
df_rd2_wk1_8_allsamples <- filter(df_round2, days_post_fmt == "7" | days_post_fmt == "56")
df_rd2_wk1_8 <- filter(df_rd2_wk1_8_allsamples, sample_type == "stool")
df_rd2_wk1_8_cageA <- filter(df_rd2_wk1_8, donor == "Cr")
df_rd2_wk1_8_cageB <- filter(df_rd2_wk1_8, donor == "BL6-1")
df_rd2_wk1_8_cageC <- filter(df_rd2_wk1_8, donor == "mr1ko")
df_rd2_wk1_8_cageD <- filter(df_rd2_wk1_8, donor == "Pooled1")

### round 1 D7/D56 abundance --------

A = df_rd1_wk1_8_cageA %>%
  drop_na(Phylum)%>%
  ggplot(aes(days_post_fmt, Abundance))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  ggtitle("PDxCre(KO)")+
  labs(x = "Days Post FMT")+
  scale_x_continuous(breaks = seq(7, 56, by = 49))+
  theme_classic()+
  scale_fill_manual(values = c("#290AD8FF","#D5D5FFFF","#72D9FFFF","#AAF7FFFF","#1E8E99FF", "#FFEE99FF",
                               "#FF2B00FF","#D9F0D3FF","#FFAD72FF","#F76D5EFF","#D82632FF","#FF7080FF"))+
  theme(axis.text.x = element_text(color = "#000000", size = 12))+
  theme(axis.text.y = element_text(color = "#000000", size = 12))+
  theme(axis.title.x = element_text(size = 15, color = "#000000"))+
  theme(axis.title.y = element_text(size = 15, color = "#000000"))+
  theme(legend.title = element_text(size = 15))+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
B = df_rd1_wk1_8_cageB %>%
  drop_na(Phylum)%>%
  ggplot(aes(days_post_fmt, Abundance))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  ggtitle("Mr1KO")+
  labs(x = "Days Post FMT")+
  scale_x_continuous(breaks = seq(7, 56, by = 49))+
  theme_classic()+
  scale_fill_manual(values = c("#290AD8FF","#D5D5FFFF","#72D9FFFF","#AAF7FFFF","#1E8E99FF", "#FFEE99FF",
                               "#FF2B00FF","#D9F0D3FF","#FFAD72FF","#F76D5EFF","#D82632FF","#FF7080FF"))+
  theme(axis.text.x = element_text(color = "#000000", size = 12))+
  theme(axis.text.y = element_text(color = "#000000", size = 12))+
  theme(axis.title.x = element_text(size = 15, color = "#000000"))+
  theme(axis.title.y = element_text(size = 15, color = "#000000"))+
  theme(legend.title = element_text(size = 15))+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
D = df_rd1_wk1_8_cageD %>%
  drop_na(Phylum)%>%
  ggplot(aes(days_post_fmt, Abundance))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  ggtitle("Pooled 1")+
  labs(x = "Days Post FMT")+
  scale_x_continuous(breaks = seq(7, 56, by = 49))+
  theme_classic()+
  scale_fill_manual(values = c("#290AD8FF","#D5D5FFFF","#72D9FFFF","#AAF7FFFF","#1E8E99FF", "#FFEE99FF",
                               "#FF2B00FF","#D9F0D3FF","#FFAD72FF","#F76D5EFF","#D82632FF","#FF7080FF"))+
  theme(axis.text.x = element_text(color = "#000000", size = 12))+
  theme(axis.text.y = element_text(color = "#000000", size = 12))+
  theme(axis.title.x = element_text(size = 15, color = "#000000"))+
  theme(axis.title.y = element_text(size = 15, color = "#000000"))+
  theme(legend.title = element_text(size = 15))+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
wk1_8_abundance_rd1 <- ggarrange(A + rremove("xlab"), B + rremove("ylab") + rremove("xlab"), D + rremove("ylab") + rremove("xlab")
          , nrow = 1, common.legend = TRUE, legend = "right") 
## add back in the x axis title 
rd1_1_8_abdc <- annotate_figure(wk1_8_abundance_rd1, bottom =  textGrob("Days Post FMT", gp = gpar(cex = 1.3)))

### round 2 D7/D56 abundance ------

E = df_rd2_wk1_8_cageA %>%
  drop_na(Phylum)%>%
  ggplot(aes(days_post_fmt, Abundance))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  ggtitle("PDxCre(KO)")+
  labs(x = "Days Post FMT")+
  scale_x_continuous(breaks = seq(7, 56, by = 49))+
  theme_classic()+
  scale_fill_manual(values = c("#290AD8FF","#D5D5FFFF","#72D9FFFF","#AAF7FFFF","#1E8E99FF", "#FFEE99FF",
                               "#FF2B00FF","#D9F0D3FF","#FFAD72FF","#F76D5EFF","#D82632FF","#FF7080FF"))+
  theme(axis.text.x = element_text(color = "#000000", size = 12))+
  theme(axis.text.y = element_text(color = "#000000", size = 12))+
  theme(axis.title.x = element_text(size = 15, color = "#000000"))+
  theme(axis.title.y = element_text(size = 15, color = "#000000"))+
  theme(legend.title = element_text(size = 15))+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
G = df_rd2_wk1_8_cageB %>%
  drop_na(Phylum)%>%
  ggplot(aes(days_post_fmt, Abundance))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  ggtitle("C57BL6")+
  labs(x = "Days Post FMT")+
  scale_x_continuous(breaks = seq(7, 56, by = 49))+
  theme_classic()+
  scale_fill_manual(values = c("#290AD8FF","#D5D5FFFF","#72D9FFFF","#AAF7FFFF","#1E8E99FF", "#FFEE99FF",
                               "#FF2B00FF","#D9F0D3FF","#FFAD72FF","#F76D5EFF","#D82632FF","#FF7080FF"))+
  theme(axis.text.x = element_text(color = "#000000", size = 12))+
  theme(axis.text.y = element_text(color = "#000000", size = 12))+
  theme(axis.title.x = element_text(size = 15, color = "#000000"))+
  theme(axis.title.y = element_text(size = 15, color = "#000000"))+
  theme(legend.title = element_text(size = 15))+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
H = df_rd2_wk1_8_cageB %>%
  drop_na(Phylum)%>%
  ggplot(aes(days_post_fmt, Abundance))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  ggtitle("Mr1KO")+
  labs(x = "Days Post FMT")+
  scale_x_continuous(breaks = seq(7, 56, by = 49))+
  theme_classic()+
  scale_fill_manual(values = c("#290AD8FF","#D5D5FFFF","#72D9FFFF","#AAF7FFFF","#1E8E99FF", "#FFEE99FF",
                               "#FF2B00FF","#D9F0D3FF","#FFAD72FF","#F76D5EFF","#D82632FF","#FF7080FF"))+
  theme(axis.text.x = element_text(color = "#000000", size = 12))+
  theme(axis.text.y = element_text(color = "#000000", size = 12))+
  theme(axis.title.x = element_text(size = 15, color = "#000000"))+
  theme(axis.title.y = element_text(size = 15, color = "#000000"))+
  theme(legend.title = element_text(size = 15))+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
I = df_rd2_wk1_8_cageD %>%
  drop_na(Phylum)%>%
  ggplot(aes(days_post_fmt, Abundance))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  ggtitle("Pooled 1")+
  labs(x = "Days Post FMT")+
  scale_x_continuous(breaks = seq(7, 56, by = 49))+
  theme_classic()+
  scale_fill_manual(values = c("#290AD8FF","#D5D5FFFF","#72D9FFFF","#AAF7FFFF","#1E8E99FF", "#FFEE99FF",
                               "#FF2B00FF","#D9F0D3FF","#FFAD72FF","#F76D5EFF","#D82632FF","#FF7080FF"))+
  theme(axis.text.x = element_text(color = "#000000", size = 12))+
  theme(axis.text.y = element_text(color = "#000000", size = 12))+
  theme(axis.title.x = element_text(size = 15, color = "#000000"))+
  theme(axis.title.y = element_text(size = 15, color = "#000000"))+
  theme(legend.title = element_text(size = 15))+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
wk1_8_abundance_rd2 <- ggarrange(E + rremove("xlab"), G + rremove("ylab") + rremove("xlab"), H + rremove("ylab") + rremove("xlab"), I + rremove("ylab") + rremove("xlab")
                                 , nrow = 1, common.legend = TRUE, legend = "right") 
rd2_1_8_abdc <- annotate_figure(wk1_8_abundance_rd2, bottom =  textGrob("Days Post FMT", gp = gpar(cex = 1.3))) 

### bray curtis alternatives to ordination CC207 -----
library(tidyverse)
library(vegan)
install.packages("Hmisc")

## shared_tbl <- read_tsv("OTU_day_name_BC.txt") 

## days to look into -- week 1 and week 8
days_wanted <- c(7, 56)

shared_tbl <- read_tsv("OTU_DAY_DONOR_RD1_RD2_no_cecum_colon_no_rd3.txt") %>%
  mutate(days = as.numeric(str_replace(Group, ".*D", "")))%>%
  filter(days %in% days_wanted)%>%
  select(Group, starts_with("OTU"))
shared_df <- shared_tbl %>%
  column_to_rownames("Group")
## GENERATE DISTANCE MATRIX // bray curtis rarified to 18XX
mice_dist <- avgdist(shared_df, sample = 1828)
mice_dist %>%
  as.matrix()%>%
  as.tibble(rownames = "samples")%>%
  pivot_longer(-samples)%>%
  filter(samples < name)%>%
  separate(samples, into = c("mouse_a", "day_a"), "D", convert = TRUE)%>%
  separate(name, into = c("mouse_b", "day_b"), "D", convert = TRUE)%>%
  mutate(comparison = case_when(
    mouse_a != mouse_b & day_a == 7 & day_a == day_b ~ "week1",
    mouse_a != mouse_b & day_a == 56 & day_a == day_b ~ "week8",
    mouse_a == mouse_b & day_a == 56 & day_b == 7 ~ "same",
    TRUE ~ NA_character_),
    comparison = factor(comparison, levels = c("week1", "week8", "same"))
    ) %>% 
  drop_na()%>%
  ggplot(aes(x=comparison, y=value))+
  geom_jitter(width = 0.25, color = "#1E8E99FF")+
  stat_summary(fun.data = median_hilow, color = "red", size = 1, 
               fun.args = list(conf.int=0.50))+
  labs(x = NULL, y = "Bray-Curtis distances")+
  scale_x_discrete(breaks = c("same", "week1", "week8"),
                   labels = c("Intra-mouse\ndistances between 7 and\n56 days post-FMT", 
                              "Inter-mouse\ndistances at 7 days\npost-FMT", 
                              "Inter-mouse\ndistances at 56 days\npost-FMT"))+
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))+
  theme_classic()+
  ggtitle(label = "Bray Curtis Distances at Week 1 and Week 8 post-FMT")
### this is a bray curtis distance plot to compare
### A: distance differences between different mice at week 1 D7
### B: distance difference between different mice at week 8 D56
### C: distance difference between SAME mice between D7 and D56

### CC204 -----
days_wanted <- c(7:56)

shared_tbl <- read_tsv("OTU_DAY_DONOR_RD1_RD2_no_cecum_colon_no_rd3.txt") %>%
  mutate(days = as.numeric(str_replace(Group, ".*D", ""))) %>%
  filter(days %in% days_wanted) %>%
  select(Group, starts_with("OTU"))
shared_df <- shared_tbl %>%
  column_to_rownames("Group")
## GENERATE DISTANCE MATRIX // bray curtis rarified to 18XX
mice_dist <- avgdist(shared_df, sample = 1828)

mice_dist %>%
  as.matrix() %>%
  as.tibble(rownames = "samples") %>%
  pivot_longer(-samples) %>%
  filter(samples < name) %>%
  mutate(mouse_a = str_replace(samples, "D.*", ""),
         mouse_b = str_replace(name, "D.*", ""),
         day_a = as.numeric(str_replace(samples, ".*D", "")),
         day_b = as.numeric(str_replace(name, ".*D", "")),
         #cage_a = str_replace(samples, "-.*", ""),
         Donor = str_replace(name, "-.*", ""),
         diff = abs(day_a - day_b)) %>%
  filter(mouse_a == mouse_b & diff > 0) %>%
  group_by(diff, mouse_a, Donor) %>%
  summarize(median = median(value)) %>%
  ungroup() %>%
  print(n = 100)
  ggplot(aes(x = diff, y = median, color = Donor, group = mouse_a))+
  geom_line(size = 1)+
  labs(x = "Days", y = "Median Bray-Curtis distance")+
  scale_x_continuous(breaks = seq(7, 56, by = 7))+
  scale_color_manual(name = "Donor", breaks = c("1A", "1B", "1D", "2A", "2B", "2C", "2D"),
                     labels = c("PDxCre(KO)", "Mr1KO", "Pooled 1", "PDxCre(KO)", 
                                "C57BL6", "Mr1KO", "Pooled 1"),
                     values = c("#290AD8FF","#FFAD72FF","#1E8E99FF", "#290AD8FF",
                                "#FF2B00FF","#FFAD72FF","#1E8E99FF"))+
  theme_classic()
### bray curtis distances over the entire 8 week period (minus days 1-7) 
### each line represents a single mouse 










## PDX received mice OTU change thru 8 weeks -----
pdx_rec <- subset_samples(ps3_cut3000, donor == "Cr")
pdx_rec <- subset_samples(pdx_rec, sample_type == "stool")
pdx_rec <- subset_samples(pdx_rec, days_post_fmt != "0")
df_pdx_rec <- data.table(psmelt(pdx_rec))

### when piping into ggplot, need to specify to plot as factor, otherwise it will plot as integer and keep a numerical axis

df_pdx_rec %>%
  drop_na(Phylum)%>%
  ggplot(aes(x = as.factor(days_post_fmt), y = Abundance))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  ggtitle("PDxCre(KO) Received abundance")+
  labs(x = "Days Post FMT")+
  theme_classic()+
  scale_fill_manual(values = c("#290AD8FF","#D5D5FFFF","#72D9FFFF","#AAF7FFFF","#1E8E99FF", "#FFEE99FF",
                               "#FF2B00FF","#D9F0D3FF","#FFAD72FF","#F76D5EFF","#D82632FF","#FF7080FF"))+
  theme(axis.text.x = element_text(color = "#000000", size = 12))+
  theme(axis.text.y = element_text(color = "#000000", size = 12))+
  theme(axis.title.x = element_text(size = 15, color = "#000000"))+
  theme(axis.title.y = element_text(size = 15, color = "#000000"))+
  theme(legend.title = element_text(size = 15))+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

### MR1KO received OTUs thru 8 wks -----
mr1ko_rec <- subset_samples(ps3_cut3000, donor == "mr1ko")
mr1ko_rec <- subset_samples(mr1ko_rec, sample_type == "stool")
mr1ko_rec <- subset_samples(mr1ko_rec, days_post_fmt != "0")
df_mr1ko_rec <- data.table(psmelt(mr1ko_rec))

df_mr1ko_rec %>%
  drop_na(Phylum)%>%
  ggplot(aes(x = as.factor(days_post_fmt), y = Abundance))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  ggtitle("Mr1KO Received abundance")+
  labs(x = "Days Post FMT")+
  theme_classic()+
  scale_fill_manual(values = c("#290AD8FF","#D5D5FFFF","#72D9FFFF","#AAF7FFFF","#1E8E99FF", "#FFEE99FF",
                               "#FF2B00FF","#D9F0D3FF","#FFAD72FF","#F76D5EFF","#D82632FF","#FF7080FF"))+
  theme(axis.text.x = element_text(color = "#000000", size = 12))+
  theme(axis.text.y = element_text(color = "#000000", size = 12))+
  theme(axis.title.x = element_text(size = 15, color = "#000000"))+
  theme(axis.title.y = element_text(size = 15, color = "#000000"))+
  theme(legend.title = element_text(size = 15))+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))



### C57BL6 received thru 8 weeks -----
bl6_rec <- subset_samples(ps3_cut3000, donor == "BL6-1")
bl6_rec <- subset_samples(bl6_rec, sample_type == "stool")
bl6_rec <- subset_samples(bl6_rec, days_post_fmt != "0")
df_bl6_rec <- data.table(psmelt(bl6_rec))

df_bl6_rec %>%
  drop_na(Phylum)%>%
  ggplot(aes(x = as.factor(days_post_fmt), y = Abundance))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  ggtitle("C57BL6 Received abundance")+
  labs(x = "Days Post FMT")+
  theme_classic()+
  scale_fill_manual(values = c("#290AD8FF","#D5D5FFFF","#72D9FFFF","#AAF7FFFF","#1E8E99FF", "#FFEE99FF",
                               "#FF2B00FF","#D9F0D3FF","#FFAD72FF","#F76D5EFF","#D82632FF","#FF7080FF"))+
  theme(axis.text.x = element_text(color = "#000000", size = 12))+
  theme(axis.text.y = element_text(color = "#000000", size = 12))+
  theme(axis.title.x = element_text(size = 15, color = "#000000"))+
  theme(axis.title.y = element_text(size = 15, color = "#000000"))+
  theme(legend.title = element_text(size = 15))+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))




### pooled 1 received thru 8 wks -------
pooled1_rec <- subset_samples(ps3_cut3000, donor == "Pooled1")
pooled1_rec <- subset_samples(pooled1_rec, sample_type == "stool")
pooled1_rec <- subset_samples(pooled1_rec, days_post_fmt != "0")
pooled1_rec <- subset_samples(pooled1_rec, cage != "3B")
df_pooled1_rec <- data.table(psmelt(pooled1_rec))

df_pooled1_rec %>%
  drop_na(Phylum)%>%
  ggplot(aes(x = as.factor(days_post_fmt), y = Abundance))+
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill")+
  ggtitle("Pooled 1 Received abundance")+
  labs(x = "Days Post FMT")+
  theme_classic()+
  scale_fill_manual(values = c("#290AD8FF","#D5D5FFFF","#72D9FFFF","#AAF7FFFF","#1E8E99FF", "#FFEE99FF",
                               "#FF2B00FF","#D9F0D3FF","#FFAD72FF","#F76D5EFF","#D82632FF","#FF7080FF"))+
  theme(axis.text.x = element_text(color = "#000000", size = 12))+
  theme(axis.text.y = element_text(color = "#000000", size = 12))+
  theme(axis.title.x = element_text(size = 15, color = "#000000"))+
  theme(axis.title.y = element_text(size = 15, color = "#000000"))+
  theme(legend.title = element_text(size = 15))+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))

### time series line plot abundance -----
df_pooled1_rec %>%
  drop_na(Phylum)%>%
  ggplot(aes(x = days_post_fmt, y = Abundance, color = Phylum, group = paste0(Phylum, cage)))+
  geom_line()


mice_dist %>%
  as.matrix() %>%
  as.tibble(rownames = "samples") %>%
  pivot_longer(-samples) %>%
  filter(samples < name) %>%
  mutate(mouse_a = str_replace(samples, "D.*", ""),
         mouse_b = str_replace(name, "D.*", ""),
         day_a = as.numeric(str_replace(samples, ".*D", "")),
         day_b = as.numeric(str_replace(name, ".*D", "")),
         #cage_a = str_replace(samples, "-.*", ""),
         Donor = str_replace(name, "-.*", ""),
         diff = abs(day_a - day_b)) %>%
  filter(mouse_a == mouse_b) %>%
  group_by(diff, mouse_a, Donor) %>%
  summarize(median = median(value)) %>%
