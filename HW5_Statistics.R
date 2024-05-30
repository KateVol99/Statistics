install.packages("ade4")
install.packages("vegan")
install.packages("ggplot2")
install.packages("igraph")
# Load the vegan package
library(ade4)
library(vegan)
library(ggplot2)
library(igraph)

data(doubs)

env <- doubs$env

spe <- doubs$fish

spa <- doubs$xy

head(spe)[, 1:8]

str(spe)

str(env)

summary(env)

range(env)

names(env)


ab <- table(unlist(spe))

barplot(ab, las = 1, 
        xlab = "Abundance class", ylab = "Frequency",
        col = grey(5:0/5))


heatmap(abs(cor(env)), 
        
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))


env.z <- decostand(env, method = "standardize")


round(apply(env.z, 2, mean), 1)


apply(env.z, 2, sd)
exists("spe")



spe.rda <- rda(spe.hel ~ ., data = env.z)

summary(spe.rda)


fwd.sel <- ordiR2step(rda(spe.hel ~ 1, data = env.z), 
                      scope = formula(spe.rda), 
                      direction = "forward",
                      R2scope = TRUE, 
                      pstep = 1000,
                      trace = FALSE) 


spe.hel <- decostand(spe, method = "hellinger")

spe.hel <- decostand(spe, method = "hellinger")

spe.rda <- rda(spe.hel ~ ., data = env.z)

fwd.sel <- ordiR2step(rda(spe.hel ~ 1, data = env.z), 
                      scope = formula(spe.rda), 
                      direction = "forward",
                      R2scope = TRUE, 
                      pstep = 1000,
                      trace = FALSE) 

fwd.sel$call


spe.rda.signif <- rda(spe.hel ~ dfs + oxy + bdo, data = env.z)

RsquareAdj(spe.rda.signif)


anova.cca(spe.rda.signif, step = 1000)

anova.cca(spe.rda.signif, step = 1000, by = "term")

ordiplot(spe.rda.signif, scaling = 1, type = "text")

ordiplot(spe.rda.signif, scaling = 2, type = "text")

perc <- round(100*(summary(spe.rda.signif)$cont$importance[2, 1:2]), 2)

sc_si <- scores(spe.rda.signif, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(spe.rda.signif, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(spe.rda.signif, display="bp", choices=c(1, 2), scaling=1)

plot(spe.rda.signif,
     scaling = 1, 
     type = "none",
     frame = FALSE,
     
     xlim = c(-1,1), 
     ylim = c(-1,1),
     main = "Triplot RDA - scaling 1",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)

points(sc_si, 
       pch = 21, 
       col = "black", 
       bg = "steelblue",
       cex = 1.2) 

points(sc_sp, 
       pch = 22, 
       col = "black",
       bg = "#f2bd33", 
       cex = 1.2)

text(sc_sp + c(0.03, 0.09), 
     labels = rownames(sc_sp), 
     col = "grey40", 
     font = 2, # bold
     cex = 0.6)

arrows(0,0, # start them from (0,0)
       sc_bp[,1], sc_bp[,2],
       col = "red", 
       lwd = 3)

text(x = sc_bp[,1] -0.1, 
     y = sc_bp[,2] - 0.03, 
     labels = rownames(sc_bp), 
     col = "red", 
     cex = 1, 
     font = 2)

env.topo <- subset(env.z, select = c(alt, dfs, slo, flo))
env.chem <- subset(env.z, select = c(pH, har, pho, nit, amm,
                                     oxy, bdo))


spe.partial.rda <- rda(spe.hel, env.chem, env.topo)
summary(spe.partial.rda)

ordiplot(spe.partial.rda, scaling = 2, main = "Doubs River partial RDA - Scaling 2")

diversity(spe, "shannon")

diversity(spe, "simpson")


ggplot(spa, aes(x = x, y = y)) +        
  geom_path(col = "lightblue", lwd = 5) +
  geom_label(label = c(1:30)) +
  theme_void() +
  labs(title = "Location of the sampling sites for the Doubs dataset") 

summary(spa)

spa$site <- 1:nrow(spa)
spa <- spa[-8, ]
which(spa$y < 82)

spa$group <- NA
spa$group[spa$y < 82] <- 1
spa$group[spa$y >= 82 & spa$y < 156] <- 2
spa$group[spa$y >= 156] <- 3


ggplot(data = spa) +
  geom_point(aes(x = x, 
                 y = y, 
                 col = as.factor(group)), 
             size = 4) +
  labs(color = "Groups", 
       x = "Longitude", 
       y = "Latitude") +
  scale_color_manual(values = c("#3b5896", "#e3548c", "#ffa600")) +
  theme_classic() + 
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18))

dist(spe)
class(dist(spe))
as.matrix(dist(spe))

spe.D.Euclid <- dist(x = spe, method = "euclidean")
is.euclid(spe.D.Euclid)

dim(as.matrix(dist(spe)))
adjacency_matrix <- as.matrix(dist(spe)) <= 0.01

graph <- graph.adjacency(adjacency_matrix, mode = "undirected")
plot(graph, main = "Graph of proximity between sampling points")




