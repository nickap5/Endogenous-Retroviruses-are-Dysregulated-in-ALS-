library(tidyverse)

# Fill this in with appropriate values
df <- data.frame(Balanced_Accuracy_Random = c(0.9,
                                              0.9,
                                              0.9,
                                              0.9,
                                              0.9,
                                              0.9,
                                              0.9,
                                              0.9,
                                              0.9,
                                              0.9),
                 Difference_Sig = c(0.2,
                                    0.2,
                                    0.2,
                                    0.2,
                                    0.2,
                                    0.2,
                                    0.2,
                                    0.2,
                                    0.2,
                                    0.2),
                 row.names = 1:10)
x_labs = factor(rownames(df), levels = unique(rownames(df)))

wilcox.test(df$Balanced_Accuracy_Random,
            mu = 0.7 # Sig Value
)
# Plot Differences
ggplot(df, aes(x=x_labs, y=Difference_Sig)) + 
  geom_point(size = 4) +
  geom_hline(yintercept = 0.09,color = "red") + # Replace with Median Diff
  ylim(0,0.12) +
  ggtitle("Differences in Balanced Accuracy Between\nBalanced Random and Balanced DEA Models") +
  xlab("Random Trial Number") + ylab("Difference in Balanced Accuracies") +
  theme_bw(base_size = 13) 
ggsave("Bootstrap_Res_CTX_MoR.jpeg")