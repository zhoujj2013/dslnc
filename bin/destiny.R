library(ggplot2)

df<-read.delim('all.gene.filter.score.txt',header = T,sep = '\t')

ggplot(df, aes(x =DSscore), alpha=0.4)+geom_density(aes(fill = Gene), alpha=0.4)

ggsave('destiny.pdf',width = 5,height = 4)
