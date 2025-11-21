#R version 4.5.1 (macOS)
#18th Nov 2025

Sys.setenv(LANGUAGE = "en")
dir<-""
setwd(dir)

#load packages
pacman::p_load("MASS")
pacman::p_load("ggplot2")
pacman::p_load("gridExtra")
pacman::p_load("grid")
pacman::p_load("gtable")
pacman::p_load("stats")
pacman::p_load("car")
pacman::p_load("tidyverse")
pacman::p_load("dplyr")
pacman::p_load("multcomp") 
pacman::p_load("rlang") 

rm(list=ls()) #clean the environment

# Define the parameters to iterate over
rho_values <- c(0.9, 0.6)
size_values <-c(2000,1000,500,250,100,50)
sce_values <- c("tirzepatide","mDiet", "MEDI0382") #"tirzepatide","mDiet", "MEDI0382"

# Loop through each combination of rho and sce
for (rho in rho_values) {
  for (sce in sce_values) {
    for(size in size_values){
 
    nsim <- 10000
    n <- size 
    mu1 <- ifelse(sce == "tirzepatide", 105.6, ifelse(sce == "mDiet", 89, 95.9))
    s1 <- ifelse(sce == "tirzepatide", 22.92, ifelse(sce == "mDiet", 12, 18.9))
    pwc <- ifelse(sce == "tirzepatide", -0.178, ifelse(sce == "mDiet", -0.029, -0.023)) #tirze pwc cut from -0.209 to -0.1, -0.05
    
    
#//////////////////////////////////
#1. Create a data frame to contain the simulation results----
    results <- data.frame(
      #simulation setting (1st row only)
      sce = integer(nsim + 1),
      pwc = integer(nsim + 1),
      rho = numeric(nsim + 1),
      nobs_all = integer(nsim + 1),
      nsim = integer(nsim + 1),
      theta = numeric(nsim + 1),
      
      #summary stats for dataset specific values (1st row only)
      alloc_mean = numeric(nsim + 1),
      alloc_sd = numeric(nsim + 1),
      ATE_set_mean = numeric(nsim + 1),
      ATE_set_sd = numeric(nsim + 1),
      ATE_t_mean = numeric(nsim + 1),
      ATE_t_sd = numeric(nsim + 1),
      ATE_ancova2c_mean = numeric(nsim + 1),
      ATE_ancova2c_sd = numeric(nsim + 1),
      ATE_ancova_mean = numeric(nsim + 1),
      ATE_ancova_sd = numeric(nsim + 1),
      
      #performance measures 1st row only
      bias_t = numeric(nsim + 1),
      mcerror_bias_t = numeric(nsim + 1),
      bias_ancova2c = numeric(nsim + 1),
      mcerror_bias_ancova2c = numeric(nsim + 1),
      bias_ancova = numeric(nsim + 1),
      mcerror_bias_ancova = numeric(nsim + 1),
      EmpSE_t = numeric(nsim + 1),
      mcerror_EmpSE_t = numeric(nsim + 1),
      EmpSE_ancova2c = numeric(nsim + 1),
      mcerror_EmpSE_ancova2c = numeric(nsim + 1),
      EmpSE_ancova = numeric(nsim + 1),
      mcerror_EmpSE_ancova = numeric(nsim + 1),
      MSE_t = numeric(nsim + 1),
      mcerror_MSE_t = numeric(nsim + 1),
      MSE_ancova2c = numeric(nsim + 1),
      mcerror_MSE_ancova2c = numeric(nsim + 1),
      MSE_ancova = numeric(nsim + 1),
      mcerror_MSE_ancova = numeric(nsim + 1),
      cover_t = numeric(nsim + 1),
      cover_ancova2c = numeric(nsim + 1),
      cover_ancova = numeric(nsim + 1),
      
      #dataset specific results from 2nd row and onwards
      setid = numeric(nsim + 1),
      alloc = numeric(nsim + 1), #allocation balance
      ATE_set = numeric(nsim + 1),
      ATE_t = numeric(nsim + 1),
      ll_t = numeric(nsim + 1),
      ul_t = numeric(nsim + 1),
      CIrange_t = numeric(nsim + 1),
      ATE_ancova2c = numeric(nsim + 1),
      ll_ancova2c = numeric(nsim + 1),
      ul_ancova2c = numeric(nsim + 1),
      CIrange_ancova2c = numeric(nsim + 1),
      R2_ancova2c = numeric(nsim + 1),
      ATE_ancova = numeric(nsim + 1),
      ll_ancova = numeric(nsim + 1),
      ul_ancova = numeric(nsim + 1),
      CIrange_ancova = numeric(nsim + 1),
      R2_ancova = numeric(nsim + 1)
      
    )
    
    #add simulation setting values to the results
    results$sce[1] <- sce
    results$pwc[1] <- pwc
    results$rho[1] <- rho
    results$nobs_all[1] <- n
    results$nsim[1] <- nsim
    theta <- mu1 * pwc
    results$theta[1] <- theta
    
    for (i in 2:(nsim + 1)) {
#2. Generate data------
      #Bivariate normal distribution
      data <- mvrnorm(n, mu = c(mu1, mu1), Sigma = matrix(c(s1^2, s1^2 * rho, s1^2 * rho, s1^2), 2, 2))
      data <- as.data.frame(data)
      names(data) <- c("x", "y0") # More descriptive than V1, V2
      
      # Calculate quartiles
      quartiles <- quantile(data$x, probs = seq(0.25, 1, 0.25))
      data <- data %>%
        #y1 calculation using case_when
        mutate(y1 = case_when(
          x <= quartiles[1] ~ y0 * (1 + pwc * 0.2),
          x <= quartiles[2] ~ y0 * (1 + pwc * 0.8),
          x <= quartiles[3] ~ y0 * (1 + pwc * 1.2),
          TRUE ~ y0 * (1 + pwc * 1.8)
        )) %>%
        #Random assignment of treatment (z)
        mutate(z = sample(0:1, n, replace = TRUE, prob = c(0.5, 0.5))) %>%
        #Observed outcome (y) based on treatment assignment
        mutate(y = ifelse(z == 0, y0, y1)) %>%
        #Calculate the individual observed change
        mutate(obschange = y - x) %>%
        #Select and arrange columns
        dplyr::select(z, x, y, obschange, y0, y1) %>%
        arrange(x)
      
      #dataset specific values to go into the results dataframe
      results$alloc[i] <- mean(data$z) # allocation balance (closet to 0.5)
      ATE_set <- mean(data$y1) - mean(data$y0) #counter factual: cannot be observed
      results$ATE_set[i] <- ATE_set
      
#3. T-test of follow up weight ----
      ttest_model <- t.test(y ~ factor(z), data = data) #t-test of y1 between trt0 and trt1
      ATE_t <- as.numeric(ttest_model$estimate[2] - ttest_model$estimate[1])
      ll_t <- -ttest_model$conf.int[2]
      ul_t <- -ttest_model$conf.int[1]
      CIrange_t <- ul_t - ll_t
      
      results$ATE_t[i] <- ATE_t
      results$ll_t[i] <- ll_t
      results$ul_t[i] <- ul_t
      results$CIrange_t[i] <- CIrange_t
      
#4. ANCOVA2_centering-------------------------------------------------------
      data$y_centred <- as.numeric(scale(data$y, center = TRUE, scale = FALSE)) #centering
      data$x_centred <- as.numeric(scale(data$x, center = TRUE, scale = FALSE)) #centering
      data$z_centred <- as.numeric(scale(data$z, center = TRUE, scale = FALSE))
      ancova2c <- lm(y_centred ~ x_centred * factor(z_centred), data = data) #model
      summary_ancova2c <- summary(ancova2c)
      coefficients2c <- summary_ancova2c$coefficients
      
      #diagnosis (off for simulation)
      if (0) {
        par(mfrow = c(2, 2)) # Set up the graphics layout
        plot(ancova2c) # Produce y0 diagnostic plots
        # Q-Q plot of residuals
        qqPlot(ancova2c$residuals, main = "Q-Q Plot")
        # Shapiro-Wilk test for normality
        shapiro.test(ancova2c$residuals)
        # Check Homogeneity of Variances
        # Levene's Test for Homogeneity of Variances
        leveneTest(y ~ factor(z), data = data)
      }
      
      #predict mean based on the model
      mean_x <- mean(data$x_centred[data$z == 1])
      hypothesis_matrix <- rbind(c(0, mean_x, 1, mean_x) - c(0, mean_x, 0, 0))
      lincom <- glht(ancova2c, linfct = hypothesis_matrix) # calc. difference between treatment level for mean y0_centred
      summary_lincom <- summary(lincom)
      confint_lincom <- confint(lincom)
      
      ATE_ancova2c <- as.numeric(summary_lincom$test$coefficients) #extract the stats
      ll_ancova2c <- confint_lincom$confint[, "lwr"]
      ul_ancova2c <- confint_lincom$confint[, "upr"]
      CIrange_ancova2c <- ul_ancova2c - ll_ancova2c
      R2_ancova2c <- summary_ancova2c$r.squared
      
      results$ATE_ancova2c[i] <- ATE_ancova2c
      results$ll_ancova2c[i] <- ll_ancova2c
      results$ul_ancova2c[i] <- ul_ancova2c
      results$CIrange_ancova2c[i] <- CIrange_ancova2c
      results$R2_ancova2c[i] <- R2_ancova2c
      
#5. ANCOVA_centering-------------------------------------------------------
      #data$y_centred <- as.numeric(scale(data$y, center = TRUE, scale = FALSE)) #centering
      #data$x_centred <- as.numeric(scale(data$x, center = TRUE, scale = FALSE)) #centering
      #data$z_centred <- as.numeric(scale(data$z, center = TRUE, scale = FALSE))
      ancova <- lm(y_centred ~ x_centred + z_centred, data = data) #model
      summary_ancova <- summary(ancova)
      coefficients <- summary_ancova$coefficients
      
      #predict mean based on the model
      #mean_x <- mean(data$x_centred[data$z == 1])
      hypothesis_matrix <- matrix(c(0, 0, 1), 1)
      rownames(hypothesis_matrix) <- "ATE: z1 vs z0"
      lincom <- glht(ancova, linfct = hypothesis_matrix) # calc. difference between treatment level for mean y0_centred
      summary_lincom <- summary(lincom)
      confint_lincom <- confint(lincom)
      
      ATE_ancova <- as.numeric(summary_lincom$test$coefficients) #extract the stats
      ll_ancova <- confint_lincom$confint[, "lwr"]
      ul_ancova <- confint_lincom$confint[, "upr"]
      CIrange_ancova <- ul_ancova - ll_ancova
      R2_ancova <- summary_ancova$r.squared
      
      results$ATE_ancova[i] <- ATE_ancova
      results$ll_ancova[i] <- ll_ancova
      results$ul_ancova[i] <- ul_ancova
      results$CIrange_ancova[i] <- CIrange_ancova
      results$R2_ancova[i] <- R2_ancova
      
#6. ----
      results$setid[i] <- (i - 1)
    } 
    
#7. ----
    #summary stat for dataset sepecific ATE values
    alloc_mean <- mean(results$alloc[2:(1 + nsim)])
    alloc_sd <- sd(results$alloc[2:(1 + nsim)])
    results$alloc_mean[1] <- alloc_mean
    results$alloc_sd[1] <- alloc_sd
    results$alloc[1] <- alloc_mean
    
    ATE_set_mean <- mean(results$ATE_set[2:(1 + nsim)])
    ATE_set_sd <- sd(results$ATE_set[2:(1 + nsim)])
    results$ATE_set_mean[1] <- ATE_set_mean
    results$ATE_set_sd[1] <- ATE_set_sd
    results$ATE_set[1] <- ATE_set_mean
    
    ATE_t_mean <- mean(results$ATE_t[2:(1 + nsim)])
    ATE_t_sd <- sd(results$ATE_t[2:(1 + nsim)])
    results$ATE_t_mean[1] <- ATE_t_mean
    results$ATE_t_sd[1] <- ATE_t_sd
    results$ATE_t[1] <- ATE_t_mean
    
    ATE_ancova2c_mean <- mean(results$ATE_ancova2c[2:(1 + nsim)])
    ATE_ancova2c_sd <- sd(results$ATE_ancova2c[2:(1 + nsim)])
    results$ATE_ancova2c_mean[1] <- ATE_ancova2c_mean
    results$ATE_ancova2c_sd[1] <- ATE_ancova2c_sd
    results$ATE_ancova2c[1] <- ATE_ancova2c_mean

    ATE_ancova_mean <- mean(results$ATE_ancova[2:(1 + nsim)])
    ATE_ancova_sd <- sd(results$ATE_ancova[2:(1 + nsim)])
    results$ATE_ancova_mean[1] <- ATE_ancova_mean
    results$ATE_ancova_sd[1] <- ATE_ancova_sd
    results$ATE_ancova[1] <- ATE_ancova_mean
    
    #performance measures
    results$EmpSE_t[1] <- sd(results$ATE_t[2:(1 + nsim)])
    results$mcerror_EmpSE_t[1] <- (2 * (nsim - 1))^-0.5 * results$EmpSE_t[1]
    results$EmpSE_ancova2c[1] <- sd(results$ATE_ancova2c[2:(1 + nsim)])
    results$mcerror_EmpSE_ancova2c[1] <- (2 * (nsim - 1))^-0.5 * results$EmpSE_ancova2c[1]
    results$EmpSE_ancova[1] <- sd(results$ATE_ancova[2:(1 + nsim)])
    results$mcerror_EmpSE_ancova[1] <- (2 * (nsim - 1))^-0.5 * results$EmpSE_ancova[1]
    
    results$bias_t[1] <- mean(results$ATE_t[2:(1 + nsim)]) - ATE_set_mean
    results$mcerror_bias_t[1] <- (nsim)^-0.5 * results$EmpSE_t[1] #mcerror_bias needs to be after empSE as it factors in empse
    results$bias_ancova2c[1] <- mean(results$ATE_ancova2c[2:(1 + nsim)]) - ATE_set_mean
    results$mcerror_bias_ancova2c[1] <- (nsim)^-0.5 * results$EmpSE_ancova2c[1]
    results$bias_ancova[1] <- mean(results$ATE_ancova[2:(1 + nsim)]) - ATE_set_mean
    results$mcerror_bias_ancova[1] <- (nsim)^-0.5 * results$EmpSE_ancova[1]
    
    results$MSE_t[1] <- mean((results$ATE_t[2:(1 + nsim)] - ATE_set)^2)
    results$MSE_ancova2c[1] <- mean((results$ATE_ancova2c[2:(1 + nsim)] - ATE_set)^2)
    results$MSE_ancova[1] <- mean((results$ATE_ancova[2:(1 + nsim)] - ATE_set)^2)
    #prep to calc mcerror for MSE
    theta_hat <- results$ATE_t[2:(1 + nsim)]
    numerator <- sum(((theta_hat - ATE_set)^2 - results$MSE_t[1])^2)
    denominator <- nsim * (nsim - 1)
    results$mcerror_MSE_t[1] <- sqrt(numerator / denominator)
    #prep to calc mcerror for MSE
    theta_hat <- results$ATE_ancova2c[2:(1 + nsim)]
    numerator <- sum(((theta_hat - ATE_set)^2 - results$MSE_ancova2c[1])^2)
    denominator <- nsim * (nsim - 1)
    results$mcerror_MSE_ancova2c[1] <- sqrt(numerator / denominator)
    #prep to calc mcerror for MSE
    theta_hat <- results$ATE_ancova[2:(1 + nsim)]
    numerator <- sum(((theta_hat - ATE_set)^2 - results$MSE_ancova[1])^2)
    denominator <- nsim * (nsim - 1)
    results$mcerror_MSE_ancova[1] <- sqrt(numerator / denominator)
    
    ATE_set2 <- ATE_set 
    results$cover_t[1] <- results %>%
      slice(-1) %>%
      filter(ll_t < ATE_set2 & ul_t > ATE_set2) %>%
      nrow() 
    results$cover_ancova2c[1] <- results %>%
      slice(-1) %>%
      filter(ll_ancova2c < ATE_set2 & ul_ancova2c > ATE_set2) %>%
      nrow() 
    results$cover_ancova[1] <- results %>%
      slice(-1) %>%
      filter(ll_ancova < ATE_set2 & ul_ancova > ATE_set2) %>%
      nrow() 
    
#8 Export results----
    #summary
    dir.create("log", showWarnings = FALSE)
    filename <- paste0("log/", sce, "_pwc(", pwc, ")", "rho(", rho, ")n(", n, ")nsim(", nsim, ").csv")
    write.csv(results, file = filename, row.names = FALSE)
    #dataset specific records
    dir.create("log2", showWarnings = FALSE)
    results_set <- results %>% dplyr::select(1:6, 38:54)
    filename2 <- paste0("log2/2", sce, "_pwc(", pwc, ")", "rho(", rho, ")n(", n, ")nsim(", nsim, ").csv")
    write.csv(results_set, file = filename2, row.names = FALSE)
    
    } 
  } 
} 

print("Simulation and saving complete for all parameter combinations.")