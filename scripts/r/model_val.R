if (!require(moments)) install.packages("moments")  # for kurtosis
library(moments)
library(terra)
library(dplyr)
library(ggplot2)
library(caret)

# Pull in full raster result
wf_haz <- rast("D:\\GIST909\\GoogleExports\\norm_wildfire_hazard2.tif")

# Compute and create full raster histogram
# ---- Extract values (remove NA) ----
vals <- values(wf_haz, na.rm = TRUE)

# Compute histogram without plotting
h <- hist(vals, breaks = 40, plot = FALSE)

# Normal curve
xfit <- seq(min(vals), max(vals), length = 100)
yfit <- dnorm(xfit, mean = mean(vals), sd = sd(vals))

# Scale to histogram
yfit <- yfit * length(vals) * diff(h$breaks)[1]

# Plot histogram with expanded y-axis
hist(vals,
     breaks = 50,
     col = "skyblue",
     border = "white",
     main = "Hazard Index Value Distribution",
     xlab = "Pixel Values",
     ylab = "Frequency",
     ylim = c(0, max(c(h$counts, yfit)) * 1.1))

# Add curve
lines(xfit, yfit, col = "darkblue", lwd = 2)

# ---- Compute statistics ----
median_val <- median(vals)
sd_val     <- sd(vals)
kurt_val   <- kurtosis(vals)

# ---- Print statistics ----
cat("Raster Statistics:\n")
cat("Median: ", median_val, "\n")
cat("Standard Deviation: ", sd_val, "\n")
cat("Kurtosis: ", kurt_val, "\n")

# Histogram with stat print out
hist(vals,
     main = paste0("Wildfire Hazard Index\nMedian=", round(median_val, 3),
                   " | SD=", round(sd_val, 2),
                   " | Kurtosis=", round(kurt_val, 2)),
     xlab = "Pixel Values",
     ylab = "Frequency",
     col = "skyblue",
     breaks = 50)

# Wilcox test here
# Set working directory for sample points (created in QGIS)
setwd('D:/GIST909/QGIS/SamplePoints')
dir()

inside <- read.csv("samplevaluesInside2.csv")
outside <- read.csv("samplevaluesOutside2.csv")

# Make sure the value column is labelled "WHI"
inside_whi <- inside$WRI
outside_whi <- outside$WRI

# wilcox.test(inside_whi, outside_whi, alternative = "greater")

# create object for Wilcox/Mann-Whitney test
test <- wilcox.test(inside_whi, outside_whi, alternative = "greater")

median(inside_whi)
median(outside_whi)

# plot for distribution and median comparison
boxplot(inside_whi, outside_whi,
        names = c("Inside Fire Perimeter", "Outside Perimeter"),
        main = "Random Sample Point Comparison",
        ylab = "Wildfire Hazard Index (WHI) Value",
        col = c("#FA807260", "#7FFFD480"),
        border = "black",
        notch = TRUE,
        outline = FALSE,
        cex.lab = 1.2,
        cex.axis = 1.1)

# print test result
test

# effect
W <- test$statistic
W # print result

# set up rank bi-serial equation
n1 <- length(inside_whi)
n2 <- length(outside_whi)

# rank bi-serial correlation
r_rb <- (2 * W) / (n1 * n2) - 1
r_rb


# Start of dNBR and WHI comparison
whi <- rast("D:\\GIST909\\QGIS\\wri_clip.tif")
dnbr <- rast("D:\\GIST909\\QGIS\\dNBR_Clipped2.tif")

plot(whi)
plot(dnbr)

# Make sure they are the same extent
ext(whi)
ext(dnbr)

# stack rasters
stacked <- c(whi, dnbr)
names(stacked) <- c("WHI", "dNBR")

# create random samples for each raster
samples <- spatSample(stacked, size = 1000, method = "random", na.rm = TRUE, as.df = TRUE)

# Using equal interval
n_classes <- 4

# WHI breaks Equal Interval
whi_breaks <- seq(min(samples$WHI), max(samples$WHI), length.out = n_classes + 1)

# dNBR breaks Equal Interval
dnbr_breaks <- seq(min(samples$dNBR), max(samples$dNBR), length.out = n_classes + 1)

samples$WHI_class <- cut(samples$WHI,
                         breaks = whi_breaks,
                         include.lowest = TRUE,
                         labels = c("Low", "Moderate", "High", "VeryHigh"))

samples$dNBR_class <- cut(samples$dNBR,
                          breaks = dnbr_breaks,
                          include.lowest = TRUE,
                          labels = c("Low", "Moderate", "High", "VeryHigh"))

# check samples
head(samples)

# scatter plot for trend
ggplot(samples, aes(x = WHI, y = dNBR)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Trend Between WHI and dNBR",
       x = "Wildfire Hazard Index (WHI)",
       y = "dNBR (Burn Severity)") +
  theme_minimal()

# scatter plot - dNBR distribution across WHI classes
ggplot(samples, aes(x = WHI_class, y = dNBR)) +
  geom_boxplot(outlier.alpha = 0.1) +
  geom_jitter(alpha = 0.2, width = 0.2) +
  labs(title = "dNBR Distribution Across WHI Classes",
       x = "WHI Class",
       y = "dNBR") +
  theme_minimal()

# alternative scatter - without points
ggplot(samples, aes(x = WHI_class, y = dNBR, fill = WHI_class)) +
  geom_boxplot(alpha = 0.6) +
  labs(title = "Burn Severity (dNBR) by WHI Equal Interval Class",
       x = "WHI Class",
       y = "dNBR") +
  theme_minimal() +
  theme(legend.position = "none")

# build confusion matrix
conf_mat <- table(samples$WHI_class, samples$dNBR_class)

# quick print
conf_mat

# calculate accuracy
overall_acc <- sum(diag(conf_mat)) / sum(conf_mat)
overall_acc

# Producer's accuracy (by column)
prod_acc <- diag(conf_mat) / colSums(conf_mat)

# User's accuracy (by row)
user_acc <- diag(conf_mat) / rowSums(conf_mat)

# check accuracies
prod_acc
user_acc

# Prints table and statistics
confusionMatrix(conf_mat)

# create a data frame for plotting
conf_df <- as.data.frame(conf_mat)

# confusion matrix heatmap plot
ggplot(conf_df, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 3) +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b80") +
  labs(x = "WHI Class (Predicted)",
       y = "dNBR Class (Observed)",
       fill = "Count",
       title = "Confusion Matrix: WHI vs dNBR Classes") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

