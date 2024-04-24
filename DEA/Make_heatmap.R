library(stringr)
library(readxl)
library(tidyverse)


df <- as.data.frame(readxl::read_xls("IPA_Data_Canon.xls", col_names = TRUE))
coln = c("Canonical Pathways","All Patients", "Females", "Males")
colnames(df) <- coln
df<-df[-1,]


df2 = head(df, 10)
df_long <- pivot_longer(data = df2,
                        cols = -`Canonical Pathways`,
                        names_to = "DEA_Analysis", 
                        values_to = "Z_Score")
df_long$Z_Score = as.numeric(df_long$Z_Score)
df_long$DEA_Analysis = str_extract(df_long$DEA_Analysis, "[A-Z].*s")
df_long$DEA_Analysis = gsub("_", " ", df_long$DEA_Analysis)
df_long$`Canonical Pathways` = stringr::str_wrap(df_long$`Canonical Pathways`, width = 20)


# Create heatmap with ggplot2
library(ggplot2)
ggp <- ggplot(df_long, aes(DEA_Analysis, `Canonical Pathways`)) +                           
  geom_tile(aes(fill = Z_Score)) +
  labs(title = "Canonical Pathway Comparison Cortical DEAs", x = "", y = "") +
  scale_fill_gradient(name = "Z-Score (from IPA)",
                      high = "#FFFFFF",
                      low = "#1C57D6") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        plot.title = element_text(face = "bold",
                                  size = 16,
                                  hjust = 0.5))



ggp
ggsave(paste0("Heatmap_Canon_CTX.jpeg"),
       width = 12,
       height = 10)



df <- as.data.frame(readxl::read_xls("IPA_Data_UR.xls", col_names = TRUE))
coln = c("Upstream Regulators","All Patients", "Females", "Males")
colnames(df) <- coln
df<-df[-1,]


df2 = head(df, 10)
df_long <- pivot_longer(data = df2,
                        cols = -`Upstream Regulators`,
                        names_to = "DEA_Analysis", 
                        values_to = "Z_Score")
df_long$Z_Score = as.numeric(df_long$Z_Score)
df_long$DEA_Analysis = str_extract(df_long$DEA_Analysis, "[A-Z].*s")
df_long$DEA_Analysis = gsub("_", " ", df_long$DEA_Analysis)
df_long$`Upstream Regulators` = stringr::str_wrap(df_long$`Upstream Regulators`, width = 20)


# Create heatmap with ggplot2
library(ggplot2)
ggp <- ggplot(df_long, aes(DEA_Analysis, `Upstream Regulators`)) +                           
  geom_tile(aes(fill = Z_Score)) +
  labs(title = "Upstream Regulators Comparison Cortical DEAs", x = "", y = "") +
  scale_fill_gradient(name = "Z-Score (from IPA)",
                      low = "#6B33FF",
                      high = "#FF9033") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        plot.title = element_text(face = "bold",
                                  size = 16,
                                  hjust = 0.5))



ggp
ggsave(paste0("Heatmap_UR_CTX.jpeg"),
       width = 12,
       height = 10)