#############################################################
# Conidia & Blastospore Count Analysis
# 9 isolates, 2 replications, CRD design
#############################################################

# ---- 0. Load packages ----
# install.packages(c("agricolae", "dplyr", "ggplot2"))  # run once if not installed
install.packages("agricolae")
library(agricolae)
library(dplyr)
library(ggplot2)

# ---- 1. Read data ----
df <- read.csv("CountAssay.csv", stringsAsFactors = FALSE)
head(df)
str(df)   # check structure before proceeding

# ---- 2. Convert columns to correct type ----
df$Isolate     <- as.factor(df$Isolate)
df$Replication <- as.factor(df$Replication)


# check for NAs introduced during conversion
sum(is.na(df$Conidia_count))
sum(is.na(df$Blastospore_count))

str(df)   # confirm types are now correct

# ---- 4. Run ANOVA ----
aov_conidia <- aov(Conidia_count ~ Isolate, data = df)
summary(aov_conidia)

aov_blasto <- aov(Blastospore_count ~ Isolate, data = df)
summary(aov_blasto)

# ---- 5. Mean separation (Tukey HSD, letter grouping) ----
tukey_conidia <- HSD.test(aov_conidia, "Isolate", group = TRUE)
tukey_blasto  <- HSD.test(aov_blasto, "Isolate", group = TRUE)

tukey_conidia$groups
tukey_blasto$groups

# extract letters into data frames, keyed by Isolate
letters_conidia <- data.frame(
  Isolate = rownames(tukey_conidia$groups),
  Letter_Conidia = tukey_conidia$groups$groups
)

letters_blasto <- data.frame(
  Isolate = rownames(tukey_blasto$groups),
  Letter_Blastospore = tukey_blasto$groups$groups
)

# ---- 6. Summary table (back-transformed, actual means) ----
summary_table <- df %>%
  group_by(Isolate) %>%
  summarise(
    Mean_Conidia      = mean(Conidia_count),
    SE_Conidia        = sd(Conidia_count) / sqrt(n()),
    Mean_Blastospore  = mean(Blastospore_count),
    SE_Blastospore    = sd(Blastospore_count) / sqrt(n())
  ) %>%
  left_join(letters_conidia, by = "Isolate") %>%
  left_join(letters_blasto, by = "Isolate")

summary_table

# save summary table as CSV
write.csv(summary_table, "spore_summary_table.csv", row.names = FALSE)

# ---- 7. Correlation between conidia and blastospore counts ----
cor.test(df$Conidia_count, df$Blastospore_count, method = "pearson")

#############################################################
# ---- 8. Grouped bar chart with error bars + letters ----
#############################################################

# reshape summary table to long format for grouped bars
plot_data <- summary_table %>%
  select(Isolate, Mean_Conidia, SE_Conidia, Letter_Conidia,
         Mean_Blastospore, SE_Blastospore, Letter_Blastospore) %>%
  tidyr::pivot_longer(
    cols = -Isolate,
    names_to = c(".value", "Type"),
    names_pattern = "(Mean|SE|Letter)_(.*)"
  )

# plot_data now has columns: Isolate, Type, Mean, SE, Letter

ggplot(plot_data, aes(x = Isolate, y = Mean, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = Mean - SE, ymax = Mean + SE),
    position = position_dodge(0.8), width = 0.25
  ) +
  geom_text(
    aes(label = Letter, y = Mean + SE),
    position = position_dodge(0.8), vjust = -0.5, size = 3.5
  ) +
  scale_fill_manual(values = c("Conidia" = "#56B4E9", "Blastospore" = "#D55E00"),
                    labels = c("Conidia", "Blastospore")) +
  labs(
    x = "Isolate",
    y = "Mean Count (± SE) (10^8/mL)",
    fill = "Type",
    title = "Conidia and Blastospore Production Across Isolates"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

# save the plot
ggsave("spore_barplot.png", width = 8, height = 6, dpi = 300)
