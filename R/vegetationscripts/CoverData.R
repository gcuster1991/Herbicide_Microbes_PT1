library("tidyverse")
library("lmerTest")
library("multcomp")
library("emmeans")
library("cowplot")

cover.dat_1 <- read_csv("/Users/gordoncuster/Desktop/HerbPt1Analysis/AK_WeedCover/CoverClasses_GC_Jan_2024.csv") 

cover.dat_1$Total_Weed_Cover_Midpoint <- as.double(cover.dat_1$Total_Weed_Cover_Midpoint)

cover.dat0 <- cover.dat_1 %>% dplyr::select(Plot, SubSample, Time, Herbicide, ends_with("Midpoint")) %>%
  mutate(Herb.common = case_when(
    Herbicide == "Aatrex" ~ "atrazine",
    Herbicide == "Clarity" ~ "dicamba",
    Herbicide == "Roundup Powermax" ~ "glyphosate",
    Herbicide == "Hand" ~ "hand-weeded",
    Herbicide == "Non-Treated" ~ "non-treated")) %>%
  pivot_longer(cols = ends_with("Midpoint"),
               names_to = "Species",
               values_to = "Cover") %>%
  mutate(Species = str_remove(Species, "_Midpoint"),
         Species = case_when(
           Species == "Corn" ~ "corn",
           Species == "Pigweed" ~ "redroot pigweed",
           Species == "Lambsquarter" ~ "common lambsquarters",
           Species == "Nightshade" ~ "hairy nightshade",
           Species == "Foxtail" ~ "green foxtail",
           Species == "Unknown" ~ "other",
           Species == "Total_Weed_Cover" ~ "total weedy vegetation"))
glimpse(cover.dat0)

cover.dat0 %>%
  distinct(Species)

#cover.dat0 %>%
#  filter(Species %in% c("other", "green foxtail", "hairy nightshade")) %>%
#  group_by(Herb.common, Species) %>%
#  summarize(mean(Cover, na.rm = TRUE))

cover.dat <- cover.dat0 %>%
  filter(!Species %in% c("other", "green foxtail", "hairy nightshade")) %>%
  mutate(Species = factor(Species, 
                          levels = rev(c("corn", 
                                     "redroot pigweed", 
                                     "common lambsquarters", 
                                     "total weedy vegetation"))))



### Sample time = 1
lmer(Cover ~ Herb.common * Species + (1 | Plot),
     data = cover.dat %>% filter(Time == "T1")) -> cover.t1.lmer
VarCorr(cover.t1.lmer)
anova(cover.t1.lmer)
cover.t1.cld <- cld(emmeans(cover.t1.lmer, ~ Species), Letters = LETTERS)
cover.t1.cld

#t1.gg <- ggplot(cover.t1.cld,
#       aes(x = emmean, y = Species)) +
#  geom_bar(aes(fill = Species),
#           stat = "identity") +
#  theme_minimal_grid() +
#  xlab("Cover (%)") +
#  ylab(element_blank()) +
#  ggtitle("Pre-treatment") +
#  theme(legend.position = "none",
 #       plot.title.position = "plot") 


lmer(Cover ~ Herb.common * Species + (1 | Plot),
     data = cover.dat %>% filter(Time == "T1")) -> cover.t1.lmer
VarCorr(cover.t1.lmer)
anova(cover.t1.lmer)
cover.t1.cld <- cld(emmeans(cover.t1.lmer, ~ Herb.common | Species), 
                    Letters = LETTERS)

t1.gg <- ggplot(cover.t1.cld,
                aes(x = emmean, y = factor(Species, levels = c("common lambsquarters", "redroot pigweed", "total weedy vegetation", "corn")))) +
  geom_bar(aes(fill = Species),
           stat = "identity",
           position = "dodge") + xlim(0,100) +
  theme_minimal_grid() +
  xlab("Cover (%)") +
  ylab(element_blank()) +
  facet_wrap( ~ Herb.common,
              ncol = 1) +
  ggtitle("B. Pre-treatment") +
  theme(legend.position = "none",
        plot.title.position = "plot") 

  
### Sample time = 2
lmer(Cover ~ Herb.common * Species + (1 | Plot),
     data = cover.dat %>% filter(Time == "T2")) -> cover.t2.lmer
VarCorr(cover.t2.lmer)
anova(cover.t2.lmer)
cover.t2.cld <- cld(emmeans(cover.t2.lmer, ~ Herb.common | Species), 
                    Letters = LETTERS)
cover.t2.cld
t2.gg <- ggplot(cover.t2.cld,
       aes(x = emmean, y = reorder(Species, emmean))) +
  geom_bar(aes(fill = Species),
           stat = "identity",
           position = "dodge") + xlim(0,100) +
  theme_minimal_grid() +
  xlab("Cover (%)") +
  ylab(element_blank()) +
  facet_wrap( ~ Herb.common,
              ncol = 1) +
  theme(legend.position = "none") +
  ggtitle("A. 10 d post-treatment") +
  theme(legend.position = "none",
        plot.title.position = "plot") 

### Sample time = 3
lmer(Cover ~ Herb.common * Species + (1 | Plot),
     data = cover.dat %>% filter(Time == "T3")) -> cover.t3.lmer
VarCorr(cover.t3.lmer)
anova(cover.t3.lmer)
cover.t3.cld <- cld(emmeans(cover.t3.lmer, ~ Herb.common | Species), 
                    Letters = LETTERS)
cover.t3.cld
t3.gg <- ggplot(cover.t3.cld,
       aes(x = emmean, y = reorder(Species, emmean))) +
  geom_bar(aes(fill = Species),
           stat = "identity",
           position = "dodge") + xlim(0,100) +
  theme_minimal_grid() +
  xlab("Cover (%)") +
  ylab(element_blank()) +
  facet_wrap( ~ Herb.common,
              ncol = 1) +
  ggtitle("B. 20 d post-treatment") +
  theme(legend.position = "none",
        plot.title.position = "plot") 

t123.gg <- plot_grid(t1.gg, t2.gg, t3.gg, ncol = 3, align = "h")

ggsave("/Users/gordoncuster/Desktop/Git_projects/Herbicide_Microbes_PT1/Figures/FinalPlots/Individual Figures/Figure1_T123_vegetation.png", t123.gg, dpi = 600, bg = "white")
ggsave("/Users/gordoncuster/Desktop/Git_projects/Herbicide_Microbes_PT1/Figures/FinalPlots/Individual Figures/Figure1_T123_vegetation.eps", t123.gg, dpi = 600, bg = "white")


t1.gg
t23.gg


ggsave("T1.png", t1.gg, width = 4, height = 1.7, units = "in", dpi = 600)
ggsave("T2 and T3.png", t23.gg, width = 8, height = 5, units = "in", dpi = 600)



#

#remake bar for total cover
dat <- read_csv("/Users/gordoncuster/Desktop/HerbPt1Analysis/AK_WeedCover/CoverClasses_GC_Jan_2024.csv") 

index <- c("Clarity", "Roundup Powermax", "Aatrex", "Hand", "Non-Treated")
values <- c("Dicamba", "Glyphosate", "Atrazine-Mesotrione", "Handweeded", "Non-Treated")
dat$Herbicide <- as.factor(values[match(dat$Herbicide, index)])


dat$total_Weed_Veg_Midpoint <- dat$Corn_Midpoint + dat$Pigweed_Midpoint + dat$Nightshade_Midpoint + dat$Lambsquarter_Midpoint + dat$Foxtail_Midpoint + dat$Unknown_Midpoint

#ggplot(data = dat, aes(x = Time, y = Total_Weed_Veg, fill = Herbicide)) + geom_boxplot() + facet_grid(cols = vars(Herbicide)) + theme_classic() 

ggplot(data = dat, aes(x = Time, y = total_Weed_Veg_Midpoint, fill = Herbicide)) + geom_bar(aes(fill=Herbicide), stat="identity", position = "dodge") + facet_grid(cols = vars(Herbicide)) + 
  theme_minimal_grid() 
#ggplot(data = dat, aes(x = Time, y = total_Weed_Veg_Midpoint, fill = Herbicide)) + geom_boxplot() + facet_grid(cols = vars(Herbicide)) + theme_minimal_grid() + coord_flip()
#ggplot(data = dat, aes(x = Time, y = Pigweed_Midpoint, fill = Herbicide)) + geom_boxplot() + facet_grid(cols = vars(Herbicide)) + theme_classic() 

#ggplot(data = dat, aes(x = Time, y = Lambsquarter_Midpoint, fill = Herbicide)) + geom_boxplot() + facet_grid(cols = vars(Herbicide)) + theme_classic() 


