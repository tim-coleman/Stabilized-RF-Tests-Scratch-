source('C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Pvalue_combine.R')
qsar_fish_toxicity <- read.csv("C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Data/qsar_fish_toxicity.csv", header=FALSE, sep=";")
names(qsar_fish_toxicity)[7] <- "LC50"

set.seed(1)
out <- stable_MSE_Test.default(X = qsar_fish_toxicity[,-7], y = qsar_fish_toxicity$LC50, B = 250,
                               var = names(qsar_fish_toxicity[,-7]), n_val = nrow(qsar_fish_toxicity))

pdf(file = 'C:/Users/drain/Box Sync/Mentch Coleman Shared/RF Permutation Tests/StabilizedRFPerm/Pval_Toxicity.pdf', width = 6, height = 4)
ggplot(data.frame(out)) + 
  geom_histogram(aes(x = Pvalues, y = ..density..), bins = 10, col = 'blue', fill = 'lightskyblue') + 
  ylab("") + xlab(expression(p[S])) + ggtitle("Fish Toxicity P-value Distribution") + 
  theme_bw() + theme(plot.title = element_text(size = 14, hjust = 0.5),
                     axis.title.x = element_text(size = 13))
dev.off()
ggplot(data.frame(out)) + stat_ecdf(aes(x = Pvalues)) + theme_bw()

### Comparing against knockoffs
library(knockoff)
train.qf <- qsar_fish_toxicity[,-7] %>% as.matrix()
result.qf <- knockoff.filter(X = train.qf, y = qsar_fish_toxicity$LC50, knockoffs = function(X) create.second_order(X = train.qf), fdr = .25,
                             statistic = stat.random_forest)
print(result.qf$selected)