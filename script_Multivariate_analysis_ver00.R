########################################################################
## Description: R code to the implementation of the Multivariate
## analysis...
##
## Maintainer: UNIVALI / EMCT / LEMA - iAtlantic | FAPESC Projects
## Author: Rodrigo Sant'Ana
## Created: seg ago 22 12:30:23 2022 (-0300)
## Version: 0.0.1
##
## URL: Communications Earth & Environment -
## https://www.nature.com/commsenv/
## Doc URL:
##
## Database info:
## a) https://doi.pangaea.de/10.1594/PANGAEA.946402
##
### Commentary:
##
### Code:
########################################################################

########################################################################
######@> Setup R...

######@> Loading R packages...
library(dplyr)
library(lubridate)
library(ggplot2)
library(ggforce)
library(ggrepel)
library(tidyr)
library(cowplot)
library(janitor)
library(patchwork)
library(grImport)
library(ade4)
library(adespatial)
library(FactoMineR)
library(factoextra)
library(mvpart)
library(vegan)
library(ecotraj)
library(readxl)
library(extrafont)
library(viridis)

######@> Importing fonts...
loadfonts(device = "postscript")

######@> Custom theme for ggplot...
rgb01 <- "black"
rgb02 <- "black"
seta <- grid::arrow(length = grid::unit(0.2, "cm"), type = "open")
seta2 <- grid::arrow(length = grid::unit(0.2, "cm"), type = "open",
                     ends = "both")
my_theme <- function(base_size = 18, base_family = "Helvetica") {
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(axis.ticks = element_line(colour = rgb01),
              axis.line = element_line(colour = rgb01, size = 0.2),
              axis.text = element_text(colour = rgb02, size = 14),
              axis.title = element_text(size = 18),
              legend.background = element_blank(),
              legend.key = element_blank(),
              panel.background = element_blank(),
              panel.grid = element_line(linetype = "solid",
                                        size = 0.2,
                                        colour = "gray90"),
              plot.background = element_blank(),
              complete = TRUE)
}

#####@> Testing ggplot theme...
df <- data.frame(x = rnorm(10), y = rnorm(10),
                 z = rep(c("A", "B"), each = 5))
ggplot(data = df, aes(x = x, y = y, fill = z)) +
    geom_point(pch = 21, size = 5) +
    my_theme()

######@> Standardizing decimal marks...
options(scipen = 10)

########################################################################
######@> Loading datasets...

######@> [Pangaea datasets] Catch time-series...
db00 <- read.delim("Data/Demersal_Fisheries_Catch_Time_Series_BMM.tab",
                   skip = 19)

######@> [Pangaea datasets] Thermal affinity...
thermal <- read.delim("Data/Thermal_Species_Affinities_BMM.tab",
                      skip = 32)

########################################################################
######@> Cleaning and Preparing dataset...

######@> Rename columns / variables...
names(db00) <- c("Year", "Species", "Catch.kg")
names(thermal) <- c("Species", "Species.Cod", "AlphaID_URI",
                    "AlphaID_Semantic", "ITIS_TSN_URI",
                    "ITIS_TSN_Semantic", "FAO_ISSCAAP", "FAO_TAXOCODE",
                    "FAO_ASFIS", "Thermal_Affinity")

######@> Pivot dataset...
db01 <- db00 %>%
    pivot_wider(names_from = "Species",
                values_from = "Catch.kg") %>%
    as.data.frame()

######@> Verifying species with 0 catch...
any(colSums(db01[, 2:ncol(db01)]) == 0)

######@> Transforming data using Hellinger's, as proposed by Legendre &
######@> Gauthier (2014)...
db01.hel <- decostand(db01[, 2:ncol(db01)],
                           method = "hellinger")
row.names(db01.hel) <- db01$Year

########################################################################
######@> Applying Multivariate Analysis...

######@> RDA...
out.rda <- rda(db01.hel)

######@> Analysing year compositions...
year <- data.frame(year = as.numeric(row.names(db01.hel)))
res <- mvpart(data.matrix(db01.hel) ~ year, data = year, xv = "pick",
              xvmult = 100)

#####@> Exporting groups...
gr.mrt <- data.frame(grupos = res$where)
gr.mrt$class <- ifelse(gr.mrt$grupos == 7, "I",
                ifelse(gr.mrt$grupos == 5, "II",
                ifelse(gr.mrt$grupos == 6, "III","IV")))

######@> PCA with Hellinger transformation...
out.pca <- PCA(db01.hel)

#####@> Percentage of explanation...
tab01 <- data.frame(out.pca$eig)
tab01$id <- 1:nrow(tab01)
tab01[nrow(tab01) + 1,] <- c(0, 0, 0, 0)
tab01 <- arrange(tab01, id)

####@> Looking to that...
p00 <- ggplot(data = tab01) +
    geom_bar(aes(x = id, y = percentage.of.variance), stat = "identity") +
    geom_line(aes(x = id, y = cumulative.percentage.of.variance)) +
    labs(x = "Principal components", y = "Percentage of variance") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100.5)) +
    my_theme()
p00

#####@> Looking to PCA...

####@> Extracting values...
var <- get_pca_var(out.pca)

###@> Contributions...
var$contrib

###@> Visualizing outputs...
p01 <- fviz_pca_var(out.pca, col.var = "contrib") +
    scale_color_gradient2(name = "% of Contribution",
                          low = "blue",
                          mid = "orange",
                          high = "red",
                          midpoint = 2.5) +
    labs(title = "") +
    my_theme() +
    theme(legend.position = "right")
p01
p02 <- fviz_pca_ind(out.pca, geom = c("text"),
                    habillage = factor(gr.mrt$class),
                    col.ind = "black", labelsize = 8,
                    repel = TRUE) +
    geom_mark_ellipse(aes(color = Groups),
                      alpha = 0, tol = .005) +
    scale_x_continuous(breaks = seq(-10, 10, 2), limits = c(-10, 10)) +
    scale_y_continuous(breaks = seq(-10, 10, 2), limits = c(-10, 10)) +
    labs(x = "PC01 (20.5%)", y = "PC02 (15.5%)",
         title = "", fill = "Groups", shape = "Groups",
         colour = "Groups") +
    scale_color_viridis(discrete = TRUE, end = 0.7) +
    my_theme() +
    theme(legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 18))
p02
temp01 <- data.frame(spp = row.names(var$contrib),
                     contrib.PC01 = var$contrib[, 1],
                     contrib.PC02 = var$contrib[, 2])
temp02 <- data.frame(spp = row.names(var$coord),
                     coord.PC01 = var$coord[, 1],
                     coord.PC02 = var$coord[, 2])
df01 <- full_join(temp01, temp02, by = "spp")
df01 <- filter(df01, contrib.PC01 >= mean(temp01$contrib.PC01))
df01$spp <- gsub("_", " ", df01$spp)

temp01 <- data.frame(spp = row.names(var$contrib),
                     contrib.PC01 = var$contrib[, 1],
                     contrib.PC02 = var$contrib[, 2])
temp02 <- data.frame(spp = row.names(var$coord),
                     coord.PC01 = var$coord[, 1],
                     coord.PC02 = var$coord[, 2])
df02 <- full_join(temp01, temp02, by = "spp")
df02 <- filter(df02, contrib.PC01 < mean(temp01$contrib.PC01))
df02$spp <- gsub("_", " ", df02$spp)

temp01 <- thermal %>%
    select(Species, Thermal_Affinity) %>%
    as.data.frame()
df01 <- df01 %>%
    left_join(temp01, by = c("spp" = "Species"))
df01$Classe <- ifelse(df01$Thermal_Affinity > 21.10635, "Warm water",
                      "Cold water")
df01$origin <- 0
df02 <- df02 %>%
    left_join(temp01, by = c("spp" = "Species"))
df02$Classe <- ifelse(df02$Thermal_Affinity > 21.10635, "Warm water",
                      "Cold water")
df02$origin <- 0

p03 <- ggplot(data = df01, aes(x = coord.PC01, y = coord.PC02)) +
    geom_label_repel(aes(label = spp, colour = Classe), size = 5,
                     box.padding = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual("Thermal affinity", values = c("blue", "red")) +
    scale_x_continuous(limits = c(-0.9, 0.9)) +
    scale_y_continuous(limits = c(-0.8, 0.8)) +
    labs(x = "PC01 (20.5%)", y = "PC02 (15.5%)",) +
    my_theme() +
    theme(legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 18))
p03

######@> Multiple comparisons...

#####@> Comparison matrix...
matC <- expand.grid(ano01 = 2000:2019, ano02 = 2000:2019)
matC$teste <- matC$ano02 > matC$ano01

####@> Reference...
matC$linha01 <- as.numeric(as.factor(matC$ano01))
matC$linha02 <- as.numeric(as.factor(matC$ano02))

####@> Cleaning unecessary cases...
matC <- filter(matC, teste == TRUE)

####@> Selecting data...
g04 <- db01.hel[matC$linha01, ]
g05 <- db01.hel[matC$linha02, ]

####@> TBI...
res.TBID <- TBI(g04, g05, method = "%difference", nperm = 999,
                test.t.perm = TRUE, save.BC = TRUE)

####@> Extracting values from TBI...
tab02 <- res.TBID$BCD.mat
tab02$TBI <- res.TBID$TBI
tab02$p.value <- res.TBID$p.TBI
tab02$ID <- paste0(substr(row.names(g04), 1, 4),
                   " - ",
                   substr(row.names(g05), 1, 4))
tab02$ano01 <- as.numeric(substr(tab02$ID, 1, 4))
tab02$ano02 <- as.numeric(substr(tab02$ID, 8, 11))
tab02$Change <- gsub(" ", "", tab02$Change)

####@> Standardizing changes...
tab02$TBI_Change <- ifelse(tab02$Change == "+", tab02$TBI,
                           tab02$TBI * (-1))

###@> Heatmap 01...
p00 <- ggplot(data = tab02, aes(x = ano01, y = ano02)) +
    geom_tile(aes(fill = TBI)) +
    geom_text(aes(label = paste0(round(p.value, 4),
                                 " (", Change, ") ")),
              size = 4, fontface = "bold") +
    scale_fill_viridis_c(begin = 0.2, option = "A") +
    scale_x_continuous(breaks = seq(2000, 2019, 1)) +
    scale_y_continuous(breaks = seq(2000, 2019, 1)) +
    labs(x = "Year (Past)", y = "Year (Future)") +
    my_theme() +
    theme(legend.key.size = unit(2, "cm"),
          legend.key.width = unit(0.8, "cm"),
          legend.box.background = element_rect(colour = "black"),
          legend.background = element_rect(fill = "white"),
          legend.position = c(0.98, 0.03),
          legend.justification = c(1, 0))
p00

###@> Heatmap 02...
p01 <- ggplot(data = tab02, aes(x = ano01, y = ano02)) +
    geom_tile(aes(fill = TBI_Change)) +
    geom_text(aes(label = round(p.value, 4)),
              size = 4) +
    scale_fill_distiller(palette = "Spectral", direction = 1) +
    scale_x_continuous(breaks = seq(2000, 2019, 1)) +
    scale_y_continuous(breaks = seq(2000, 2019, 1)) +
    labs(x = "Year (t - 1)", y = "Year (t)", fill = "TBI - Change") +
    my_theme() +
    theme(legend.key.size = unit(3, "cm"),
          legend.key.width = unit(0.8, "cm"),
          legend.box.background = element_rect(colour = "black"),
          legend.background = element_rect(fill = "white"),
          legend.position = c(0.98, 0.03),
          legend.justification = c(1, 0))
p01

#####@> Beta diversity analysis...
mat01 <- db01[, 2:ncol(db01)]
out01 <- beta.div(mat01, method = "hellinger", sqrt.D = FALSE,
                  samp = FALSE, nperm = 9999, save.D = FALSE)

tab02 <- data.frame(spp = names(out01$SCBD),
                    SCBD = out01$SCBD)
tab02 <- arrange(tab02, desc(SCBD))

########################################################################
##
##                  Creative Commons License 4.0
##                       (CC BY-NC-SA 4.0)
##
##  This is a humam-readable summary of (and not a substitute for) the
##  license (https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
##
##  You are free to:
##
##  Share - copy and redistribute the material in any medium or format.
##
##  The licensor cannot revoke these freedoms as long as you follow the
##  license terms.
##
##  Under the following terms:
##
##  Attribution - You must give appropriate credit, provide a link to
##  license, and indicate if changes were made. You may do so in any
##  reasonable manner, but not in any way that suggests the licensor
##  endorses you or your use.
##
##  NonCommercial - You may not use the material for commercial
##  purposes.
##
##  ShareAlike - If you remix, transform, or build upon the material,
##  you must distributive your contributions under the same license
##  as the  original.
##
##  No additional restrictions — You may not apply legal terms or
##  technological measures that legally restrict others from doing
##  anything the license permits.
##
########################################################################
