###############################################################################
#
# Figures for PEER Covariate Testing
#
#
#
###############################################################################

library(ggplot2)


sigcis <- read.delim('ciseQTL/significant_cis.txt', sep = ' ', header = F, stringsAsFactors = F)
sigcis$time = 1
sigcis$cutoff = 1
sigcis$cov = 'ok'


for(i in 1:nrow(sigcis)){
  
  line = strsplit(sigcis[i,2], '_')
  line = line[[1]]
  sigcis[i,'time'] <- line[1]
  sigcis[i,'cutoff'] <- line[3]
  sigcis[i,'cov'] <- strsplit(line[4], '\\.')[[1]][1] 
  sigcis$cov <- factor(sigcis$cov, levels = c('nobatch', 'nocrc', '0', '1', '2', '3', '4', '5', '10', '15', '20', '25', '30'))
  sigcis$time <- factor(sigcis$time, levels = c('T4', 'T30', 'T44'))
  
}

siglist = unique(sigcis$cutoff)
subset = NULL
for(i in siglist){
  subset = sigcis[sigcis$cutoff == i,]
  mainlab = paste('Effect of Covariates on cis-eQTL Discovery')
  ylim = max(subset$V1)*1.2
  p <- ggplot(subset, aes(cov, V1)) + 
    geom_point(mapping = aes(cov, V1, color = time), size = 3) +
    labs(title = mainlab, x = 'Covariates Used', y = 'Number of significant cis-eQTL')
  print(p)
  
}


p <- ggplot(sigcis, aes(cov, V1)) +
  facet_grid(rows = vars(cutoff), scales = 'free')+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(mapping = aes(cov, V1, color = time), size = 3) +
  labs(title = mainlab, x = 'Covariates Used', y = 'Number of significant cis-eQTL')
print(p)



