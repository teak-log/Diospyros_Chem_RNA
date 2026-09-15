library(ggtree)
library(ggplot2)
library(corHMM)
library(phytools)


### Choosing a Character Model ###
dio.tree <- read.tree("C:/Users/Admin/Documents/rna/cagee/trial2/astral2_unann.nwk")
dio.tree <- drop.tip(dio.tree, c("sandwicensis", "tireliae", "glans", "ptparviflora", "inexplorata"))  # dropping species with no soil data available
dio.data <- read.csv("soil_preference.txt", sep = "\t", header = TRUE)
soil.pref <- read.csv("soil_preference.txt", sep = "\t", row.names = 1,
                      stringsAsFactors = TRUE)
soil.pref <- setNames(soil.pref[,1], rownames(soil.pref))

cols=c("#005F60", "#FF8800")
ggtree(dio.tree) %<+% dio.data +
  geom_tippoint(aes(colour = soil_preference), size = 4) +
  guides(color = guide_legend(title = "native soil type")) +
  scale_color_manual(values = cols) + geom_tiplab(offset = 0.05, fontface=3)  # plot phylogeny showing species' soil preferences

fitER <- fitMk(dio.tree,soil.pref,model = "ER")  # equal rates
fitARD <- fitMk(dio.tree,soil.pref, model = "ARD")  # all rates different
fit01 <- fitMk(dio.tree,soil.pref,model = matrix(c(0,1,0,0),2,2,byrow = TRUE))  # only non-ultramafic to ultramafic transitions allowed
fit10 <- fitMk(dio.tree,soil.pref,model = matrix(c(0,0,1,0),2,2,byrow = TRUE))  # only ultramafic to non-ultramafic transitions allowed
fitER  # check the optimised rates
fitARD
fit01
fit10

aic <- c(AIC(fitER),AIC(fitARD),AIC(fit01),AIC(fit10))  # compare transition models using AIC scores
data.frame(model=c("ER","ARD","sed->ult","ult->sed"),
           logL=c(logLik(fitER),logLik(fitARD),
                  logLik(fit01),logLik(fit10)),AIC=aic,delta.AIC=aic-min(aic))  # print AIC score table


### Marginal Ancestral State Reconstruction ###
dio.data <- data.frame(sp=names(soil.pref),soil.pref=as.numeric(soil.pref)-1)
fit.marginal <- corHMM(dio.tree,dio.data,node.states = "marginal",
                       rate.cat = 1)

fit.marginal
head(fit.marginal$states)

soil.pref <- soil.pref[dio.tree$tip.label]  #reorders soil preference table
names(cols) <- levels(soil.pref)
tip.cols <- cols[soil.pref]

plotTree(dio.tree, ftype="i", offset=0.5)  # plot ancestral state reconstructions on phylogeny
tiplabels(
            pch = 21,           # filled circle
            bg = tip.cols,      # fill colour for each tip
            cex = 1.5,          # point size
            adj = 0.5           # position: 0.5 = on tip, 1 = left of label
        )
legend("topleft", legend = levels(soil.pref),pch = 21,pt.cex = 1.5,
       pt.bg = cols,bty = "n",cex = 0.8)
nodelabels(pie = fit.marginal$states,piecol = c("#005F60", "#FF8800"),cex=0.5)


### Stochastic Character Mapping ###
mtrees <- make.simmap(dio.tree,soil.pref,model = "ARD", nsim = 10000, Q="mcmc",
                      prior=list(use.empirical=TRUE),samplefreq=10)  # 10000 stochastic character maps to build posterior distribution of possible histories

par(mar=c(5.1,4.1,2.1,2.1))
plot(d<-density(sapply(mtrees,function(x) x[["Q"]][1,2]),bw=0.005),bty="n",main = "",
     xlab = "q",xlim=c(0,30),ylab="Posterior Density from MCMC",las=1,cex.axis=0.8)  # plot posterior density of non-ultramafic to ultramafic transitions from MCMC
polygon(d,col = make.transparent("darkgray",0.25))  # colour area under the curve
lines(d<-density(sapply(mtrees,function(x) x[["Q"]][2,1]),bw=0.005),bty="n",main = "",
      xlab = "q",xlim=c(0,30),ylab="Posterior Density from MCMC",las=1,cex.axis=0.8)  # on the same graph, plot posterior density of ultramafic to non-ultramafic transitions from MCMC
polygon(d,col = make.transparent("darkgray",0.25))

pd <- summary(mtrees)
pd

plot(pd,colors=cols,fsize=1,ftype="i",lwd=2,offset=0.4,ylim=c(-1,Ntip(dio.tree)), cex=c(0.5,0.3))  # plot ancestral state reconstructions on phylogeny
legend("topleft",legend = levels(soil.pref),pch = 22,pt.cex = 1.5,pt.bg = cols, bty = "n",cex = 0.8)
