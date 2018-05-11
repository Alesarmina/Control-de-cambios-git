#######use phyloseq######
library(phyloseq)
library(ggplot2)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(vegan)
library(plyr)
####import biom otu table with sample metadata and taxonomy######
setwd("/Users/Beto/Desktop/BIOSFERA")

bios_otu_table<-("./BIOSFERA_otu-table_n64_rare21712_barrel_collapsed_w_envsmd_CN.biom")
bios_tree<-("./rep_set.tre")
bios_physeq_data<-import_biom(bios_otu_table,
                              bios_tree,
                              parseFunction = parse_taxonomy_greengenes)
bios_physeq_data


####make PCoA and ordination plot of unifrac distances####
theme_set(theme_bw())
bios_pcoa<-ordinate(bios_physeq_data, "PCoA", "unifrac", weighted=TRUE)
PCoA_unifrac_plot<-plot_ordination(bios_physeq_data,
                                   bios_pcoa,
                                   color = "Plant",
                                   shape = "Environment")
PCoA_unifrac_plot + geom_point(size=7, alpha=0.75) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "PCoA Axis 1 (38.8%)", y = "PCoA Axis 2 (16.5%)")


####CCA of Bray-Curtis distance and environmental variables as vectors####
####Using Phyloseq#####
summary(bios_physeq_data@sam_data)
bios_physeq_data@sam_data<-bios_physeq_data@sam_data[,c(-1,-6,-10:-14)]
bios_physeq_data@sam_data$C.N<-as.numeric(bios_physeq_data@sam_data$C.N)
bios_physeq_data@sam_data$Moisture<-as.numeric(bios_physeq_data@sam_data$Moisture)
bios_physeq_data@sam_data$Potential_C.mineralization<-as.numeric(bios_physeq_data@sam_data$Potential_C.mineralization)
bios_physeq_data@sam_data$SOM<-as.numeric(bios_physeq_data@sam_data$SOM)
bios_physeq_data@sam_data$pH<-as.numeric(bios_physeq_data@sam_data$pH)
bios_physeq_data@sam_data$Root_volume<-as.numeric(bios_physeq_data@sam_data$Root_volume)

summary(bios_physeq_data@sam_data)


bios_cca<-ordinate(bios_physeq_data,
                    method = "CCA",
                    distance = "bray",
                    formula = bios_physeq_data ~ C.N + pH +
                     Moisture + SOM + Potential_C.mineralization)
cca_bios_plot<-plot_ordination(bios_physeq_data,
                                        bios_cca,
                                        color = "Plant",
                                        shape = "Environment")

arrowmat<-vegan::scores(bios_cca, display = "bp")
arrowdf<-data.frame(labels = rownames(arrowmat), arrowmat)
library(plyr)
arrowdf$labels1<-revalue(arrowdf$labels, c("Potential_C.mineralization" = "C. Min",
                               "C.N" = "C:N"))
arrow_map<-aes(xend = CCA1,
               yend = CCA2,
               x = 0,
               y = 0,
               shape = NULL,
               color = NULL,
               label = labels1)
label_map<-aes(x = 1.2 * CCA1, y = 1.2 * CCA2,
               shape = NULL, color = NULL, label = labels1)
arrowhead<-arrow(length = unit(0.05, "npc"))

cca_bios_plot_vec<-cca_bios_plot +
  geom_segment(arrow_map,
               size = 0.5, 
               data = arrowdf,
               color = "gray",
               arrow = arrowhead) +
  geom_text(label_map, size = 5, data = arrowdf, show.legend = FALSE) +
  geom_point(size=3, alpha=0.75) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

print(cca_bios_plot_vec)

ggsave("cca_otus_plot_env_all.png",
       width = 25,
       height = 20,
       units = "cm",
       dpi = 600)


####NMDS of UniFrac distances pylum level####
bios_phylum<-tax_glom(bios_physeq_data, taxrank = "Phylum")
bios_phylum_nmds<-ordinate(bios_phylum,
                           method = "NMDS",
                           distance = "wunifrac")
NMDS_bios_phylum_unifrac_plot<-plot_ordination(bios_phylum,
                                               bios_phylum_nmds,
                                               color = "Plant",
                                               shape = "Environment") +
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  theme_bw()
print(NMDS_bios_phylum_unifrac_plot)
vec.sp<-envfit(bios_phylum_nmds, as.data.frame(t(bios_phylum@otu_table)), perm=1000)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$species<-rownames(vec.sp.df)
taxa_phylum<-as.data.frame(bios_phylum@tax_table)
vec.sp.df$phylum<-taxa_phylum$Phylum
vec.sp.df$pval<-vec.sp$vectors$pvals
vec.sp.df<-filter(vec.sp.df, pval < 0.01)



NMDS_bios_phylum_unifrac_plot_df<-plot_ordination(bios_phylum, bios_phylum_nmds, justDF = TRUE)

ggplot(data = NMDS_bios_phylum_unifrac_plot_df, aes(NMDS1, NMDS2)) +
  geom_point(data = NMDS_bios_phylum_unifrac_plot_df, aes(colour = Plant, shape = Environment), size = 3, alpha = 0.75) +
  geom_segment(data = vec.sp.df, aes(x=0,xend=NMDS1*0.25,y=0,yend=NMDS2*0.25),
               arrow = arrow(length = unit(0.3, "cm")),
               colour="grey", inherit.aes = FALSE) +
  geom_text(data=vec.sp.df, aes(x=NMDS1*0.25,y=NMDS2*0.25,label=phylum),size=3) +
  theme_bw() +
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t")


###Stacked barplots###
#melt to long fomat table(ggplot format)
#prune phlya below 1% in each sample
bios_phylum<-bios_physeq_data %>%   
  tax_glom(taxrank = "Phylum") %>%  #aglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #transform to relative abundance
  psmelt() %>% #melt to long format
  filter(Abundance > 0.01) %>% #filter low abundance taxa
  arrange(Phylum) #sort data frame alphabetically by phylum

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)
ggplot(bios_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  theme_bw() +
  facet_grid(~Plant, scale = "free", space = "free_x") +
  scale_y_continuous(labels = percent) + ylab("Abundace (%)")

#####Analisis de variables ambientales######

soil_data<-read.table(file = "BIOS_SoilData.txt", header = TRUE, sep = "\t", row.names = 1)
library(plyr)
soil_data<-rename(soil_data, c("Potential_C.mineralization" = "C. Min",
                                         "Soil_Moisture" = "Soil Moisture",
                                         "C.N" = "C:N",
                              "Root_volume" = "Root volume"))
soil_data_NMDS<-soil_data[c(1:30),c(5:7,11,12)]

summary(soil_data_NMDS)

anova_pH<-summary(aov(pH~Biome, soil_data))
anova_pH

anova_soilMoist<-summary(aov(Soil.Moisture..gravimetric.~Biome, soil_data))
anova_soilMoist
TukeyHSD(aov(Soil.Moisture..gravimetric.~Biome, soil_data), conf.level = 0.95)


nmds_soil<-metaMDS(soil_data_NMDS, distance = "euclidean")
nmds_soil_df<-nmds_soil$points
nmds_soil_df<-as.data.frame(nmds_soil_df)

nmds_soil_df$Biome<-soil_data$Biome
nmds_soil_df$Plant<-soil_data$Plant

vec.env<-envfit(nmds_soil, as.data.frame(soil_data_NMDS), perm=1000)
vec.env.df<-as.data.frame(vec.env$vectors$arrows*sqrt(vec.env$vectors$r))
vec.env.df$species<-rownames(vec.env.df)
samples<-as.data.frame(rownames(soil_data_NMDS))


soil_nmds_plot<-ggplot(nmds_soil_df, aes(MDS1, MDS2)) +
  geom_point(data = nmds_soil_df, aes(color = Plant, shape = Biome), size = 3, alpha = 0.75) +
  geom_segment(data = vec.env.df, aes(x=0,xend=NMDS1*0.25,y=0,yend=NMDS2*0.25),
             arrow = arrow(length = unit(0.3, "cm")),
             colour="grey", inherit.aes = FALSE) +
  geom_text(data=vec.env.df, aes(x=NMDS1*0.25,y=NMDS2*0.25,label=species),size=5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

 
print(soil_nmds_plot)

ggsave("nmds_plot_env_all_C:N.png",
       width = 25,
       height = 20,
       units = "cm",
       dpi = 600)

##environmental variables only Desert & Savannah

soil_data_NMDS_DS<-soil_data[c(7:30),c(5:8,11,12)]
soil_data_NMDS_DS$Root_volume<-as.numeric(soil_data_NMDS_DS$Root_volume)

nmds_soil_DS<-metaMDS(soil_data_NMDS_DS, distance = "euclidean")
nmds_soil_DS_df<-nmds_soil_DS$points
nmds_soil_DS_df<-as.data.frame(nmds_soil_DS_df)

nmds_soil_DS_df$Biome<-soil_data$Biome[7:30]
nmds_soil_DS_df$Plant<-soil_data$Plant[7:30]

vec.env_DS<-envfit(nmds_soil_DS, as.data.frame(soil_data_NMDS_DS), perm=1000)
vec.env.df_DS<-as.data.frame(vec.env_DS$vectors$arrows*sqrt(vec.env_DS$vectors$r))
vec.env.df_DS$species<-rownames(vec.env.df_DS)
samples<-as.data.frame(rownames(soil_data_NMDS))


soil_nmds_plot_DS<-ggplot(nmds_soil_DS_df, aes(MDS1, MDS2)) +
  geom_point(data = nmds_soil_DS_df, aes(color = Plant, shape = Biome), size = 3, alpha = 0.75) +
  geom_segment(data = vec.env.df_DS, aes(x=0,xend=NMDS1*0.25,y=0,yend=NMDS2*0.25),
               arrow = arrow(length = unit(0.3, "cm")),
               colour="grey", inherit.aes = FALSE) +
  geom_text(data=vec.env.df_DS, aes(x=NMDS1*0.25,y=NMDS2*0.25,label=species),size=5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


print(soil_nmds_plot_DS)

ggsave("nmds_plot_env_all_C:N_DesertSavannah.png",
       width = 38,
       height = 25,
       units = "cm",
       dpi = 600)

