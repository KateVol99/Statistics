install.packages("vegan")
library(vegan)
data(BCI)
library(tidyverse)
library(ggplot2)


subset_data <- BCI[, 1:15]
summary(subset_data)

gathered_data <- gather(subset_data)
hist_plot <- ggplot(gathered_data, aes(x = value)) +
   geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
   facet_wrap(~ key, scales = "free") +
   theme_minimal() +
   labs(title = "Histograms of BCI Data (Subset)")
print(hist_plot)

#2
nmds_result <- metaMDS(BCI)
plot(nmds_result)



