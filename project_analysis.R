###### Growth curve section ######
# download package
install.packages("ggplot2")
install.packages('minpack.lm')
install.packages('zoo')
install.packages('dplyr')
install.packages('nls.multstart')
install.packages('tidyverse')

library("ggplot2")
library("minpack.lm")
library("zoo")
library("dplyr")
library("nls.multstart")
library('tidyverse')
install.packages('pracma')
library('pracma')

# Load data 
setwd('/Users/shiyu/Desktop/project')

ex1 <- read.csv('temp.csv') ## the OD of monobacterial grow assay
ex2 <- read.csv('interaction1.csv') ## the OD of CFSM exposure experiment
ex3 <- read.csv('st_interaction.csv') ### swfit temperature experiment

#Define fuctions
# Define the logistic growth model function

logistic_model <- function(t, r_max, N_max, N_0) {
  return(N_0 * N_max * exp(r_max * t) / (N_max + N_0 * (exp(r_max * t) - 1)))
}

# Create an empty data frame to store the results
result_df1 <- data.frame()
result_df2 <- data.frame()
result_df3 <- data.frame()

# Get unique combinations of experimental conditions
conditions1 <- unique(ex1[, c( "temp.G.", "species")])
conditions2 <- unique(ex2[, c("spent.media", "temp.G.", "species")])
conditions3 <- unique(ex3[, c("spent.media","temp.SM." ,"temp.G.", "species")])

# Loop through each experimental condition
# Experiment 1
for (i in 1:nrow(conditions1)) {
  
  print(i)
  
  current_condition1 <- conditions1[i, ]
  
  # Subset data for the current condition
  current_data1 <- subset(ex1, 
                          temp.G. == current_condition1$"temp.G." & 
                            species == current_condition1$"species")
  
  # Fit the logistic growth model
  
  fit_logistic <- nls_multstart(od ~ logistic_model(t = time.point, r_max, N_max, N_0), current_data1,
                                start_lower = c(r_max=0, N_0 = -2, N_max = 0),
                                start_upper = c(r_max=2, N_0 = 1, N_max = 1),
                                lower = c(r_max=0, N_0 = -1, N_max = 0),
                                iter = 500,
                                supp_errors = "Y")
  
  rmax1 <- coef(fit_logistic)["r_max"]
  
  K1 <- coef(fit_logistic)['N_max']
  
  K01 <- coef(fit_logistic)['N_0']
  
  AUC1 <- trapz(current_data1$time.point, current_data1$od) 
  
  # Plot the microbial growth curve for the current condition
  
  preds_time <- seq(1, 120)
  preds_growth <- logistic_model(t = preds_time, r_max = rmax1, N_max = K1, N_0 = K01)
  pred.df <- data.frame(preds_time, preds_growth)
  
  p <- ggplot(current_data1, aes(x = time.point, y = od)) +
    geom_point() +
    geom_line(data = pred.df, aes(x = preds_time, y = preds_growth)) +
    labs(title = paste("Microbial Growth Curve for", current_condition1$"species", "at", current_condition1$"temp.G.", "degrees"),
         x = "Time (hours)", y = "OD") +
    theme_minimal()
  
  # Save the plot to a file
  
  plot_filename <- paste("plot_", current_condition1$"temp.G.", current_condition1$"species", ".png", sep = "")
  ggsave(plot_filename, plot = p)
  
  # Store the results in the result data frame
  result_df1 <- rbind(result_df1, data.frame(
    Temp = current_condition1$"temp.G.",
    Species = current_condition1$"species",
    Rmax = rmax1,
    K = K1,
    AUC = AUC1
  ))
}

# Save the results data frame to a CSV file
write.csv(result_df1, "growth_parameters1.csv", row.names = FALSE)

# Experiment 2

for (i in 1:nrow(conditions2)) {
  
  print(i)
  
  current_condition2 <- conditions2[i, ]
  
  # Subset data for the current condition
  current_data2 <- subset(ex2, 
                          spent.media == current_condition2$"spent.media" &
                            temp.G. == current_condition2$"temp.G." & 
                            species == current_condition2$"species")
  
  # Fit the logistic growth model
  
  fit_logistic2 <- nls_multstart(od ~ logistic_model(t = time.point, r_max, N_max, N_0), current_data2,
                                 start_lower = c(r_max=0, N_0 = -2, N_max = 0),
                                 start_upper = c(r_max=2, N_0 = 1, N_max = 1),
                                 lower = c(r_max=0, N_0 = -1, N_max = 0),
                                 iter = 500,
                                 supp_errors = "Y")
  
  rmax2 <- coef(fit_logistic2)["r_max"]
  
  K2 <- coef(fit_logistic2)['N_max']
  
  K02 <- coef(fit_logistic2)['N_0']
  
  AUC2 <- trapz(current_data2$time.point, current_data2$od) 
  
  # Plot the microbial growth curve for the current condition
  
  preds_time <- seq(1, 120)
  preds_growth2 <- logistic_model(t = preds_time, r_max = rmax2, N_max = K2, N_0 = K02)
  pred.df2 <- data.frame(preds_time, preds_growth2)
  
  p <- ggplot(current_data2, aes(x = time.point, y = od)) +
    geom_point() +
    geom_line(data = pred.df2, aes(x = preds_time, y = preds_growth2)) +
    labs(title = paste("Microbial Growth Curve for", current_condition2$"species", "at", current_condition2$"temp.G.", "degrees"),
         x = "Time (hours)", y = "OD") +
    theme_minimal()
  
  # Save the plot to a file
  
  plot_filename <- paste("plot_", current_condition2$"temp.G.", current_condition2$"species", ".png", sep = "")
  ggsave(plot_filename, plot = p)
  
  # Store the results in the result data frame
  result_df2 <- rbind(result_df2, data.frame(
    Spent_media = current_condition2$"spent.media",
    Temp = current_condition2$"temp.G.",
    Species = current_condition2$"species",
    Rmax = rmax2,
    K = K2,
    AUC = AUC2
  ))
}

# Save the results data frame to a CSV file
write.csv(result_df2, "growth_parameters2.csv", row.names = FALSE)

# Experiment 3

for (i in 1:nrow(conditions3)) {
  
  print(i)
  
  current_condition3 <- conditions3[i, ]
  
  # Subset data for the current condition
  current_data3 <- subset(ex3, 
                          spent.media == current_condition3$"spent.media" &
                            temp.G. == current_condition3$"temp.G." & 
                            species == current_condition3$"species")
  
  # Fit the logistic growth model
  
  fit_logistic3 <- nls_multstart(OD ~ logistic_model(t = time.point, r_max, N_max, N_0), current_data3,
                                 start_lower = c(r_max=0, N_0 = -2, N_max = 0),
                                 start_upper = c(r_max=2, N_0 = 1, N_max = 1),
                                 lower = c(r_max=0, N_0 = -1, N_max = 0),
                                 iter = 500,
                                 supp_errors = "Y")
  
  rmax3 <- coef(fit_logistic3)["r_max"]
  
  K3 <- coef(fit_logistic3)['N_max']
  
  K03 <- coef(fit_logistic3)['N_0']
  
  AUC3 <- trapz(current_data3$time.point, current_data3$OD) 
  
  # Plot the microbial growth curve for the current condition
  
  preds_time3 <- seq(1, 120)
  preds_growth3 <- logistic_model(t = preds_time3, r_max = rmax3, N_max = K3, N_0 = K03)
  pred.df3 <- data.frame(preds_time3, preds_growth3)
  
  p <- ggplot(current_data3, aes(x = time.point, y = OD)) +
    geom_point() +
    geom_line(data = pred.df3, aes(x = preds_time3, y = preds_growth3)) +
    labs(title = paste("Microbial Growth Curve for", current_condition3$"species", "at", current_condition3$"temp.G.", "degrees"),
         x = "Time (hours)", y = "OD") +
    theme_minimal()
  
  # Save the plot to a file
  
  plot_filename <- paste("plot_", current_condition3$"temp.G.", current_condition3$"species", ".png", sep = "")
  ggsave(plot_filename, plot = p)
  
  # Store the results in the result data frame
  result_df3 <- rbind(result_df3, data.frame(
    Spent_media = current_condition3$"spent.media",
    Temp.SM = current_condition3$"temp.SM.",
    Temp = current_condition3$"temp.G.",
    Species = current_condition3$"species",
    Rmax = rmax3,
    K = K3,
    AUC = AUC3
  ))
}

# Save the results data frame to a CSV file
write.csv(result_df3, "growth_parameters3.csv", row.names = FALSE)

######TPC AUC~Temp section ############
# install package from GitHub
install.packages('remotes')
remotes::install_github("padpadpadpad/rTPC")
# load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(ggplot2)
library(dplyr)

# load data
rm(list=ls())
setwd('/Users/shiyu/Desktop/project')
tpc1 <- read.csv('growth_parameters1.csv')
tpc2 <- read.csv('growth_parameters2.csv')
tpc3 <- read.csv('growth_parameters3.csv')
tpc4 <- read.csv('gp4.csv')

ggplot(tpc1, aes(x = Temp, y = AUC)) +
  geom_point() +
  facet_wrap(~Species, scales = "free")

ggplot(tpc1, aes(x = Rmax, y = K)) +
  geom_point(aes(col = Temp)) +
  facet_wrap(~Species, scales = "free")

ggplot(tpc2, aes(x = Temp, y =AUC , color = Spent_media)) +
  geom_point() +
  geom_line()+
  facet_wrap(~Species, scales = "free")

ggplot(tpc3, aes(x = Temp, y = AUC,color = Spent_media)) +
  geom_point() +
  geom_line()+
  facet_wrap(~Species, scales = "free")

####calculate optimal temperature & activation energy ##############

# Create an empty data frame to store the results
result.df <- data.frame()
result.df2 <- data.frame()
result.df3 <- data.frame()
result.df4 <- data.frame()
#Get unique combinations of experimental conditions

con1 <- unique(tpc1[, c( "Species")])
con2 <- unique(tpc2[, c("Spent_media", "Species")])
con3 <- unique(tpc3[, c("Spent_media", "Species")])
con4 <- unique(tpc4[, c("Spent_media", "Species")])

# choose model
mod = 'sharpschoolhigh_1981'


for (i in 1:length(con1)) {
  print(i)
  current_species <- con1[i]
  
  # Subset data for the current condition
  #current_data <- subset(tpc1, 
  #                         Species == current_con$"Species")
  current_data <- tpc1 %>% filter(Species == current_species)
  
  # get start vals
  start_vals <- get_start_vals(current_data$Temp, current_data$AUC, model_name = 'sharpeschoolhigh_1981')
  
  # get limits
  low_lims <- get_lower_lims(current_data$Temp, current_data$AUC, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(current_data$Temp, current_data$AUC, model_name = 'sharpeschoolhigh_1981')
  
  # fit model 
  fit <- nls_multstart(AUC~sharpeschoolhigh_1981(temp = Temp, r_tref,e,eh,th, tref = 15),
                       data = current_data,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  
  if(!is.null(fit)){
    ae <- coef(fit)['e']
    Topt <- get_topt(fit)
    thermal_tolerance <- get_thermaltolerance(fit)
    
    # Store the results in the result data frame
    result.df <- rbind(result.df, data.frame(
      Species = current_species,
      Activation_Energy = ae,
      Optimum_temp = Topt,
      Thermal_tolerance = thermal_tolerance
    ))
  } else{
    result.df <- rbind(result.df, data.frame(
      Species = current_species,
      Activation_Energy = NA,
      Optimum_temp = NA,
      Thermal_tolerance = NA))
  }
}
write.csv(result.df, "tpc_para1.csv", row.names = FALSE)
for (i in 1:nrow(con2)) {
  print(i)
  current_species <- con2[i, ]
  
  current_data <- tpc2 %>% filter(Species == current_species$'Species' &
                                    Spent_media == current_species$'Spent_media')
  
  # get start vals
  start_vals <- get_start_vals(current_data$Temp, current_data$AUC, model_name = 'sharpeschoolhigh_1981')
  
  # get limits
  low_lims <- get_lower_lims(current_data$Temp, current_data$AUC, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(current_data$Temp, current_data$AUC, model_name = 'sharpeschoolhigh_1981')
  
  # fit model 
  fit <- nls_multstart(AUC~sharpeschoolhigh_1981(temp = Temp, r_tref,e,eh,th, tref = 15),
                       data = current_data,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  
  if(!is.null(fit)){
    ae <- coef(fit)['e']
    Topt <- get_topt(fit)
    thermal_tolerance <- get_thermaltolerance(fit)
    
    # Store the results in the result data frame
    result.df2 <- rbind(result.df2, data.frame(
      Species = current_species$Species,
      spent_media = current_species$Spent_media,
      Activation_Energy = ae,
      Optimum_temp = Topt,
      Thermal_tolerance = thermal_tolerance
    ))
  } else{
    result.df2 <- rbind(result.df2, data.frame(
      Species = current_species$Species,
      spent_media = current_species$Spent_media,
      Activation_Energy = NA,
      Optimum_temp = NA,
      Thermal_tolerance = NA))
  }
}
write.csv(result.df2,"tpc_para2.csv", row.names = FALSE)
for (i in 1:nrow(con3)) {
  print(i)
  current_species <- con3[i,]
  
  current_data <- tpc3 %>% filter(Species == current_species$Species &
                                    Spent_media == current_species$Spent_media)
  
  # get start vals
  start_vals <- get_start_vals(current_data$Temp, current_data$AUC, model_name = 'sharpeschoolhigh_1981')
  
  # get limits
  low_lims <- get_lower_lims(current_data$Temp, current_data$AUC, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(current_data$Temp, current_data$AUC, model_name = 'sharpeschoolhigh_1981')
  
  # fit model 
  fit <- nls_multstart(AUC~sharpeschoolhigh_1981(temp = Temp, r_tref,e,eh,th, tref = 15),
                       data = current_data,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  
  if(!is.null(fit)){
    ae <- coef(fit)['e']
    Topt <- get_topt(fit)
    thermal_tolerance <- get_thermaltolerance(fit)
    
    # Store the results in the result data frame
    result.df3 <- rbind(result.df3, data.frame(
      Species = current_species$Species,
      spent_media = current_species$Spent_media,
      Activation_Energy = ae,
      Optimum_temp = Topt,
      Thermal_tolerance = thermal_tolerance
    ))
  } else{
    result.df3 <- rbind(result.df3, data.frame(
      Species = current_species$Species,
      spent_media = current_species$Spent_media,
      Activation_Energy = NA,
      Optimum_temp = NA,
      Thermal_tolerance = NA))
  }
}
write.csv(result.df3,"tpc_para3.csv", row.names = FALSE)
for (i in 1:nrow(con4)) {
  print(i)
  current_species <- con4[i,]
  
  current_data <- tpc4 %>% filter(Species == current_species$Species &
                                    Spent_media == current_species$Spent_media)
  
  # get start vals
  start_vals <- get_start_vals(current_data$Temp, current_data$AUC, model_name = 'sharpeschoolhigh_1981')
  
  # get limits
  low_lims <- get_lower_lims(current_data$Temp, current_data$AUC, model_name = 'sharpeschoolhigh_1981')
  upper_lims <- get_upper_lims(current_data$Temp, current_data$AUC, model_name = 'sharpeschoolhigh_1981')
  
  # fit model 
  fit <- nls_multstart(AUC~sharpeschoolhigh_1981(temp = Temp, r_tref,e,eh,th, tref = 15),
                       data = current_data,
                       iter = 500,
                       start_lower = start_vals - 10,
                       start_upper = start_vals + 10,
                       lower = low_lims,
                       upper = upper_lims,
                       supp_errors = 'Y')
  
  if(!is.null(fit)){
    ae <- coef(fit)['e']
    Topt <- get_topt(fit)
    thermal_tolerance <- get_thermaltolerance(fit)
    
    # Store the results in the result data frame
    result.df4 <- rbind(result.df4, data.frame(
      Species = current_species$Species,
      spent_media = current_species$Spent_media,
      Activation_Energy = ae,
      Optimum_temp = Topt,
      Thermal_tolerance = thermal_tolerance
    ))
  } else{
    result.df4 <- rbind(result.df4, data.frame(
      Species = current_species$Species,
      spent_media = current_species$Spent_media,
      Activation_Energy = NA,
      Optimum_temp = NA,
      Thermal_tolerance = NA))
  }
}
write.csv(result.df4,"tpc_para4.csv", row.names = FALSE)
## didn't fit a mechanistic model of thermal performance in the end 
## to the data because they're based on fundamental metabolic rates, 
## but I used AUC so maybe its not appropriate
###### interaction experiment section########
library(broom)
library(ggplot2)
#  load data
rm(list = ls())
data <- read.csv("interratio.csv")
data1 <- data.frame()
data1 <- subset(data, Spent_media == "fresh")

### plot AUC~Temp in fresh media###
p <- ggplot(data1, aes(x = Temp, y = AUC, color = 'pink')) +
  geom_smooth(method = "loess", se = TRUE, fill="grey") +
  geom_point(data = data1, aes(x = Temp, y = AUC), color = 'black', alpha = 0.5) +
  facet_wrap(~ Species, ncol = 3, nrow = 3) +
  labs(title = "Fitted Thermal Performance Curves (Fresh media)",
       x = "Temperature (ºC)",
       y = "AUC") +
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_blank())  # White background
ggsave("output.png", plot = p)

#### compare the differences of growth status in CSFMs and fresh media##
num_colors <- 10
colors <- color_palette(num_colors)

a <- ggplot(data, aes(x = Temp, y = AUC, color = Spent_media)) +
  geom_smooth(aes(size = ifelse(Spent_media == "fresh", "bold", "normal")), linetype = "solid", method = "loess", se = TRUE, fill="grey") +
  geom_point(data = data, aes(x = Temp, y = AUC), color = 'black', alpha = 0.4) +
  scale_size_manual(values = c("normal" = 0.8, "bold" = 1.5)) +
  scale_color_manual(values = colors)+
  facet_wrap(~ Species, ncol = 3, nrow = 3) +
  labs(title = "Fitted Thermal Performance Curves",
       x = "Temperature (ºC)",
       y = "AUC") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_blank()) 
ggsave("outputtpc.png", plot = a)


### AUC ratio~ Temp

ggplot(data, aes(x = Temp, y = aucratio, color = Spent_media)) +
  geom_smooth(aes(size = ifelse(Spent_media == "fresh", "bold", "normal")), 
              linetype = "solid", method = "loess", se = FALSE) +
  geom_point(data = data, aes(x = Temp, y = aucratio, color = Spent_media), alpha = 0.4) +
  scale_size_manual(values = c("normal" = 0.8, "bold" = 1.5)) +
  scale_color_manual(values = colors) +  
  facet_wrap(~ Species, ncol = 3, nrow = 3) +
  labs(title = "Temperature effects on interactions",
       x = "Temperature (ºC)",
       y = "Interaction (AUC ratio)") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_blank())  # White background

## fiting regression model

data2 <- subset(data, Spent_media != "fresh")

model <- lm(aucratio ~ poly(Temp, 2)*Spent_media + Species, data = data2)
summary(model)

model2 <- lm(aucratio ~ as.character(Temp)*Spent_media + Species, data = data)
summary(model2)

model3 <- lm(aucratio ~ Temp*Spent_media + Species, data = data)
summary(model3)
## choose a better one 

plot(model)
model_summary <- tidy(model)

write.csv(model_summary, file = "model_summary.csv", row.names = FALSE)




######something about anther experiment i did and the TPCmodel i had tried.######
####but not used in the end
#########preperation########
# load packages
library(rTPC)
library(nls.multstart)
library(purrr)
library(broom)
library(ggplot2)
require(Matrix)
#  load data
rm(list = ls())
data <- read.csv("interratio.csv")

#  choose model
chosen_model <- 'sharpeschoolhigh_1981'

############## model fuction#################
# 准备初始参数函数 prepare start parameters

prepare_start_params <- function(Temp, AUC) {
  get_start_vals(Temp, AUC, model_name = chosen_model)
}


# 在fit_curve函数中使用来传递多个参数Function to fit a curve
fit_curve <- function(curve_data, Species, Spent_media) {
  start_params <- prepare_start_params(curve_data$Temp, curve_data$AUC)
  
  tryCatch({
    fit <- nls_multstart(AUC ~ sharpeschoolhigh_1981(temp = Temp, r_tref, e, eh, th, tref = 15),
                         data = curve_data,
                         iter = c(3, 3, 3, 3),
                         start_lower = start_params - 10,
                         start_upper = start_params + 10,
                         lower = get_lower_lims(curve_data$Temp, curve_data$AUC, model_name = chosen_model),
                         upper = get_upper_lims(curve_data$Temp, curve_data$AUC, model_name = chosen_model),
                         supp_errors = 'Y',
                         convergence_count = FALSE)
    # 将Species和Spent_media信息添加到模型对象的属性
    fit$Species <- Species
    fit$Spent_media <- Spent_media
    
    return(fit)
  }, error = function(e) {
    message("Error fitting curve with Species: ", Species, ", Spent_media: ", Spent_media)
    return(NULL)
  })
}


########try to fit model & get predicted data#########
#  create new data frame for saving predicted data
predicted_data <- data.frame()

#  Fit models for each Species and Spent_media combination
for (species in unique(data$Species)) {
  for (spent_media in unique(data$Spent_media)) {
    # select date 
    subset_data <- data[data$Species == species & data$Spent_media == spent_media, ]
    
    # fit model
    fit <- fit_curve(subset_data, Species = species, Spent_media = spent_media)
    
    if (!is.null(fit)) {
      #  new data frame for preds 
      new_data <- data.frame(Temp = seq(min(subset_data$Temp), max(subset_data$Temp), length.out = 100))
      new_data$Species <- species
      new_data$Spent_media <- spent_media
      
      #  Predict and collect data
      predictions <- augment(fit, newdata = new_data)
      # Spent_media&Species
      predictions$Species <- species
      predictions$Spent_media <- spent_media
      # combine together
      predicted_data <- rbind(predicted_data, predictions)
    }
  }
}

#  check result 

print(predicted_data)

#####error fitting curve####
# Error fitting curve with Species: S18, Spent_media: fresh
# Error fitting curve with Species: S18, Spent_media: S18
# Error fitting curve with Species: S18, Spent_media: S56
# Error fitting curve with Species: S18, Spent_media: WYM25_03
# Error fitting curve with Species: S56, Spent_media: fresh
# Error fitting curve with Species: S56, Spent_media: S18
# Error fitting curve with Species: S56, Spent_media: BB19_16
# Error fitting curve with Species: BB19_16, Spent_media: S18
# Error fitting curve with Species: BB19_16, Spent_media: BB19_16
# Error fitting curve with Species: BB19_16, Spent_media: S56
# Error fitting curve with Species: BB19_16, Spent_media: WYC41_02
# Error fitting curve with Species: WYC41_02, Spent_media: BB19_16
# Error fitting curve with Species: WYM25_02, Spent_media: S18
# Error fitting curve with Species: WYM25_02, Spent_media: WYM25_03
# Error fitting curve with Species: WYM25_03, Spent_media: S18
# Error fitting curve with Species: BB13_09, Spent_media: BB19_16
# Error fitting curve with Species: SP03_02, Spent_media: fresh
# Error fitting curve with Species: SP03_02, Spent_media: S18
# Error fitting curve with Species: SP03_02, Spent_media: SP03_19
# Error fitting curve with Species: SP03_02, Spent_media: WYC41_02
# Error fitting curve with Species: SP03_19, Spent_media: S18
# Error fitting curve with Species: SP03_19, Spent_media: BB19_16
# Error fitting curve with Species: SP03_19, Spent_media: SP03_19
# Error fitting curve with Species: SP03_19, Spent_media: WYM25_02
# Error fitting curve with Species: SP03_19, Spent_media: WYM25_03
##### 绘制图形 plot TPCs ####

ggplot(predicted_data, aes(x = Temp, y = .fitted, color = Spent_media)) +
  geom_line(aes(size = ifelse(Spent_media == "fresh", "bold", "normal")), linetype = "solid") +
  geom_point(data = data, aes(x = Temp, y = AUC), color = 'black', alpha = 0.5) +
  scale_size_manual(values = c("normal" = 0.8, "bold" = 1.5)) +
  facet_wrap(~ Species, ncol = 3, nrow = 3) +
  labs(title = "Fitted Thermal Performance Curves",
       x = "Temperature (ºC)",
       y = "AUC") +
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"),
        strip.background = element_blank())  # White background
