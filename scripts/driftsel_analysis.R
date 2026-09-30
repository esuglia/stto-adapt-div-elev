library(tidyverse)
# Source the necessary scripts from driftsel source code
source("rafm/RAFM.r")
source("driftsel/driftsel.r")
source("scripts/viz_trait_tidy.R")

geno1 = read.table('data/matrix_split.txt',head=FALSE)
geno = as.matrix(geno1[,2:ncol(geno1)])
pos = read.table('data/data.012.pos',head=FALSE)
ID = read.table('data/data.012.indv',head=FALSE)
pops = scan("data/pop.txt", what = character())
climate_data = read.table('data/climate_pc.csv', head=TRUE, sep = ",", row.names = 1)

pop_num = as.numeric(as.factor(pops))

pop_labels = data.frame(population = pops, code = pop_num) |> 
distinct()

n = dim(geno)[1] #get the number of individuals
p = dim(geno)[2]/2

geno2 = geno
colnames(geno2) = paste(paste('loc',sep='',rep(1:p, each=2)),sep='.',rep(1:2,p))

geno2 = cbind(pop_num, geno2)

colnames(geno2)[1] = 'subpop'

ptm = proc.time()
afm = do.all(geno2, 15000, 5000, 10)   # 15000 iterations, first 5000 as burn-in, the rest saved in every 10th iteractions
proc.time() - ptm


saveRDS(afm, "data/afm_15000_5000_10.RDS")
afm = readRDS("data/afm_15000_5000_10.RDS")


Y = read.table('data/merged_traits.csv', head=TRUE, sep = "\t")
ped_sto = as.matrix(data.frame(ID = 1:length(pop_num),
    sire = pop_num+1000,
    dam = pop_num+2000,
    sire.pop = pop_num,
    dam.pop = pop_num))

Y = left_join(Y, climate_data, by = c("pop.x" = "pop"))

covars_scaled = cbind(1:length(pop_num), 1)

traits = as.matrix(Y[c("sla",
"postvern_stem_diam",
 "postvern_height",
 "postvern_lngst_lf",
 "postvern_no_lvs"
 )])
traits_scaled = scale(traits)
traits_scaled = cbind(1:length(pop_num), traits_scaled)

library(MCMCpack)
samp = MH(afm$theta, ped_sto, 
    covars_scaled,
    traits_scaled,
    15000, 5000, 10, alt=T)

saveRDS(samp, "data/samp_new_scaled_no_covar.RDS")
# Elena recreated figure with new data on 7/30/26, so input the file "data/samp_new_scaled_covar_2026.RDS" that Rishav uploaded to GitHub:
samp = readRDS("data/samp_new_scaled_covar_2026.RDS")

# traits = 1: "sla",
# 2: "postvern_stem_diam",
# 3: "postvern_height",
# 4: "postvern_lngst_lf",
# 5: "postvern_no_lvs"

fixedpost = samp$fixed.ef
popefpost = samp$pop.ef
Gpost = samp$G
THpost = samp$theta
traits = c(1,5)
trait_names = c("SLA", "Stem Diameter", "Height", "Longest Leaf",  "No. of Leaves")
population_labels = c("BH", "CP2", "DPR", "IH", "KC2", "LV1", "LV2", "LV3", "LVTR", 
"SHA", "SQ1", "SQ3", "TM2", "WL1", "WL2", "WL3", "WV", "YO1", 
"YO10", "YO11")

# Colors matching the driftsel figures Rishav made:
c("#000",   "#CF5A6A",  "#7ACC5B",  "#428EDE",   "#6EDFE2",  "#BB27B3",  "#F2A93B",  "#9F9F9F",  "#A1FC4F", "#000",  "#CF596A",  "#7DCE62",  "#4995E1",  "#68DCE0",  "#BC2AB5",  "#F2A93B",  "#9F9F9F", "#88FD20", "#000", "#CC5364")

# "BH"   "CP2"  "DPR"  "IH"   "KC2"  "LV1"  "LV2"  "LV3"  "LVTR" "SHA"  "SQ1"  "SQ3"  "TM2"  "WL1"  "WL2"  "WL3"  "WV"   "YO1"  "YO10" "YO11"

############## Figure 5 #########################################
source("scripts/viz_trait_tidy.R")
# Visualize the trait relationship using the custom ggplot function
# Modify the function in the viz_trait_tidy.R script as needed to make the colors consistent
# Traits:
# 1: SLA
# 2: stem diam
# 3: height
# 4: longest leaf
# 5: # leaves

# Height, # leaves
fig5a = viz_traits_tidy(fixedpost, popefpost, Gpost, THpost, 
                      traits = c(3,5), 
                      size_param = 0.5, 
                      #plot_title = "Trait Relationships",
                      population_labels = population_labels,
                      trait_names = trait_names)
fig5a

# Height, stem diam
fig5b = viz_traits_tidy(fixedpost, popefpost, Gpost, THpost, 
                    traits = c(3,2), 
                    size_param = 0.5, 
                    #plot_title = "Trait Relationships",
                    population_labels = population_labels,
                    trait_names = trait_names)
fig5b

# stem diam, # leaves
fig5c = viz_traits_tidy(fixedpost, popefpost, Gpost, THpost, 
                        traits = c(2,5), 
                        size_param = 0.5, 
                        #plot_title = "Trait Relationships",
                        population_labels = population_labels,
                        trait_names = trait_names)
fig5c

# Height, SLA
fig5d = viz_traits_tidy(fixedpost, popefpost, Gpost, THpost, 
                        traits = c(3,1), 
                        size_param = 0.5, 
                        #plot_title = "Trait Relationships",
                        population_labels = population_labels,
                        trait_names = trait_names)
fig5d

library(cowplot)
title = ggdraw() + draw_label("Trait relationships from driftsel multivariate analysis", fontface = 'bold', size = 16)
plot_grid(title, plot_grid(fig5a, fig5b, fig5c, fig5d, labels = c("A", "B", "C", "D"), ncol = 2), ncol = 1, rel_heights = c(0.1, 1.5))

#ggsave("figures/fig5a-d_driftsel_7.28.25.png")

############## Suppl. Figure 2 #########################################

# Height, # leaves
suppfig2a = viz_traits_tidy(fixedpost, popefpost, Gpost, THpost, 
                        traits = c(1,5), 
                        size_param = 0.5, 
                        #plot_title = "Trait Relationships",
                        population_labels = population_labels,
                        trait_names = trait_names)
suppfig2a

# Height, stem diam
suppfig2b = viz_traits_tidy(fixedpost, popefpost, Gpost, THpost, 
                        traits = c(3,4), 
                        size_param = 0.5, 
                        #plot_title = "Trait Relationships",
                        population_labels = population_labels,
                        trait_names = trait_names)
suppfig2b

# stem diam, # leaves
suppfig2c = viz_traits_tidy(fixedpost, popefpost, Gpost, THpost, 
                        traits = c(4,5), 
                        size_param = 0.5, 
                        #plot_title = "Trait Relationships",
                        population_labels = population_labels,
                        trait_names = trait_names)
suppfig2c

library(cowplot)
title2 = ggdraw() + draw_label("Trait relationships from driftsel multivariate analysis", fontface = 'bold', size = 16)
plot_grid(title2, plot_grid(suppfig2a, suppfig2b, suppfig2c, labels = c("A", "B", "C"), ncol = 2), ncol = 1, rel_heights = c(0.3, 1.5))

#ggsave("figures/suppfig2a-c_driftsel_7.28.25.png", width = 10, height = 10)

# Caterpillar plot

# This shows each population's deviation from ancestral mean
# with credible intervals and neutral drift expectations

library(ggplot2)
samp = readRDS("data/samp_new_scaled_covar_2026.RDS")
# Extract data for one trait
trait_index <- 2  # choose your trait
traits = c(1,6)
trait_names = c("SLA", "Leaf Thickness", "Stem Diameter", "Height", "Longest Leaf",  "No. of Leaves")
population_labels = c("BH", "CP2", "DPR", "IH", "KC2", "LV1", "LV2", "LV3", "LVTR", 
"SHA", "SQ1", "SQ3", "TM2", "WL1", "WL2", "WL3", "WV", "YO1", 
"YO10", "YO11")
n_pops = length(population_labels)

plot_trait_univariate = function(trait_index = NULL) {
# Calculate posterior summaries
pop_means <- apply(samp$pop.ef[, trait_index, ], 1, mean)
pop_CI <- t(apply(samp$pop.ef[, trait_index, ], 1, quantile, 
                  probs=c(0.025, 0.975)))

# Calculate neutral expectations (as before)
neutral_dist <- matrix(NA, nrow=n_pops, ncol=dim(samp$pop.ef)[3])

for(iter in 1:dim(samp$pop.ef)[3]) {
  G_sd <- sqrt(samp$G[trait_index, trait_index, iter])
  
  for(pop in 1:n_pops) {
    theta_ii <- samp$theta[pop, pop, iter]
    neutral_dist[pop, iter] <- sqrt(theta_ii) * G_sd
  }
}

neutral_mean <- rowMeans(neutral_dist)
neutral_CI <- t(apply(neutral_dist, 1, quantile, probs=c(0.025, 0.975)))

# Create dataframe

plot_data <- data.frame(
  population = 1:n_pops,
  pop_labels = population_labels,
  observed = pop_means,
  obs_lower = pop_CI[,1],
  obs_upper = pop_CI[,2],
  neutral_sd = neutral_mean,
  neut_lower = -neutral_CI[,2],  # symmetric around 0
  neut_upper = neutral_CI[,2]
)
#plot_data$pop_labels = factor(plot_data$pop_labels, levels = population_labels)

# plot by elevation instead of alphabetically
plot_data$pop_labels = factor(plot_data$pop_labels, levels = c("TM2", "IH", "BH", "WV", "KC2", "DPR", "SHA", "WL1", "SQ1", "WL2", "WL3", "YO1", "CP2", "LV3", "SQ3", "LV2", "LV1", "LVTR", "YO11", "YO10"))

# Plot
ggplot(plot_data, aes(y=observed, x=pop_labels)) +
  # Neutral expectation ribbons
  # Observed credible intervals
  geom_errorbar(aes(ymin=obs_lower, ymax=obs_upper), 
                width=0.2, linewidth=1) +
  # Observed means
  geom_point(aes(y=observed), size=3, color="black") +
  # Reference line at ancestral mean
  geom_hline(yintercept=0, linetype="dashed", color="red") +
  geom_ribbon(aes(y=observed, x=population, ymin=neut_lower, ymax=neut_upper), 
              fill="gray80", alpha=0.3, inherit.aes = FALSE) +
  labs(x="Population", 
       y="Deviation from ancestral mean (genetic SD)",
       title=paste(trait_names[trait_index], " - Selection evidence")) +
  theme_bw(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(text=element_text(size=22))

ggsave(paste0("figures/Trait_means_", trait_names[trait_index],".png"), width = 8, height = 8)
}

plot_trait_univariate(1)
# Populations outside the gray ribbon show evidence of selection

# use cowplot to create grid of all plots
library(cowplot)
library(magick)

p1 <- ggdraw() + draw_image("figures/Trait_means_Longest leaf.png")
p2 <- ggdraw() + draw_image("figures/Trait_means_Height.png")
p3 <- ggdraw() + draw_image("figures/Trait_means_Stem Diameter.png")
p4 <- ggdraw() + draw_image("figures/Trait_means_Leaf Thickness.png")
p5 <- ggdraw() + draw_image("figures/Trait_means_SLA.png")
p6 <- ggdraw() + draw_image("figures/Trait_means_No. of Leaves.png")

plot_grid(p1, p2, p3, p4, p5, p6,
          ncol = 2, 
          labels = c("A", "B", "C", "D", "E", "F"))

ggsave("figures/Trait_means_combined.png", width = 16, height = 24)
