library(vegan)
library(ggplot2)
datos<-read.csv("MICs.csv", header=TRUE, sep=",")
as.data.frame(datos)
row.names(datos)<-datos$KEY
head(datos)
datos<-datos[,-1]
# calculate distance for NMDS
sol <- metaMDS(datos)
# Create meta data for grouping
Phylogroup<-read.csv("meta_MICS.csv", header=TRUE, sep=",")
Phylogroup = data.frame(
  sites = Phylogroup$KEY,
  amt = Phylogroup$Isolation.Date,
  row.names = "sites")
# plot NMDS using basic plot function and color points by "amt" from MyMeta
plot(sol$points, col = Phylogroup$amt)
# draw dispersion ellipses around data points
ordiellipse(sol, Phylogroup$amt, display = "sites", kind = "sd", label = T)

# same in ggplot2
NMDS = data.frame(MDS1 = sol$points[,1], MDS2 = sol$points[,2])
ggplot(data = NMDS, aes(MDS1, MDS2)) + 
  geom_point(aes(data = Phylogroup, color = Phylogroup$amt))
