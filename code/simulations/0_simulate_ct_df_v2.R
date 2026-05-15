library(survival)

set.seed(12345)

# ---------------------------
# Global parameters
# ---------------------------
nsim <- 5000
ftime <- 1.5
n_values  <- c(30, 100)
ve_values <- c(10, 20, 30, 40, 50, 60, 80)
#n_values  <- c(30, 100, 500)
#ve_values <- c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90)
fc_values <- seq(0.1, 0.9, 0.1)
cens_rate_pred_values <- c(0.05, 0.1, 0.2)

shape_true <- 1.4

root_dir <- "code/simulations/sim_ct_df_with_censore_weibull_5000"
if (!dir.exists(root_dir)) dir.create(root_dir)

# ---------------------------
# Helper function
# ---------------------------
get_scale <- function(fc, t, k) {
  t / (-log(1 - fc))^(1 / k)
}

# ===========================================================
# Scenario A — vary fc, fixed VE, fixed pred censor rate
# ===========================================================
# fixed_ve <- 30
# fixed_hr <- 1 - fixed_ve/100
# cens_rate <- 0.1 # ALS censoring rate


# for (n in n_values) {
#   for (fc in fc_values) {
    
#     scale_true <- get_scale(fc, ftime, shape_true)
    
#     folder_name <- sprintf("sim_n%d_fc%.1f_ve%d_censpred%.1f", n, fc, fixed_ve, cens_rate)
#     folder_path <- file.path(root_dir, folder_name)
#     if (!dir.exists(folder_path)) dir.create(folder_path)
    
#     for (i in 1:nsim) {
#       if (i %% 1000 == 0) cat("Scenario A:", n, fc, "i =", i, "\n")
#       # ---------------------------
#       # Observed data
#       # ---------------------------
#       group <- c(rep(0, n), rep(1, n))
#       hr_vector <- ifelse(group == 0, 1, fixed_hr)
      
#       U <- runif(2*n)
      
#       T_obs <- scale_true * (-log(U))^(1/shape_true) /
#         (hr_vector)^(1/shape_true)
      
#       C_obs <- rexp(2*n, rate = cens_rate)
      
#       time_obs  <- pmin(T_obs, C_obs, ftime)
#       event_obs <- as.integer(T_obs <= C_obs & T_obs <= ftime)
      
#       # ---------------------------
#       # Noise levels
#       # ---------------------------
#       noise_levels <- seq(0, 2, by = 0.2) 
      
#       all_pred <- list()
      
#       for (k in seq_along(noise_levels)) {
        
#         noise_sd <- noise_levels[k]
        
#         # ---------------------------
#         # Predicted (noisy placebo)
#         # ---------------------------
#         #shape_pred <- shape_true + rnorm(2*n, 0, noise_sd)
#         #shape_pred <- pmax(shape_pred, 0.2)
        
#         #scale_pred <- get_scale(fc, ftime, shape_pred)
        
#         T_pred <- scale_true * (-log(U))^(1/shape_true)
#         T_pred_noisy <- T_pred * exp(rnorm(2*n, 0, noise_sd))
        
#         C_pred <- rexp(2*n, rate = cens_rate)
        
#         time_pred  <- pmin(T_pred_noisy, C_pred, ftime)
#         event_pred <- as.integer(T_pred_noisy <= C_pred & T_pred_noisy <= ftime)
        
#         # Store
#         all_pred[[k]] <- data.frame(
#           id = 1:(2*n),
#           group = group,
          
#           time_obs = time_obs,
#           event_obs = event_obs,
          
#           time_pred = time_pred,
#           event_pred = event_pred,
          
#           noise = noise_sd
#         )
#       }
      
#       # Combine all noise levels
#       sim_df <- do.call(rbind, all_pred)
      
#       fname <- sprintf("%s/sim_%04d.csv", folder_path, i)
#       write.csv(sim_df, fname, row.names = FALSE)
#     }
#   }
# }

# ===========================================================
# Scenario B — vary VE, fixed fc, fixed pred censor rate
# ===========================================================
fixed_fc <- 0.4
cens_rate <- 0.1 # ALS censoring rate


for (n in n_values) {
  
  scale_true <- get_scale(fixed_fc, ftime, shape_true)
  
  for (ve in ve_values) {
    
    hr <- 1 - ve/100
    
    folder_name <- sprintf("sim_n%d_fc%.1f_ve%d_censpred%.1f", n, fixed_fc, ve, cens_rate)
    folder_path <- file.path(root_dir, folder_name)
    if (!dir.exists(folder_path)) dir.create(folder_path)
    
    for (i in 1:nsim) {
      if (i %% 1000 == 0) cat("Scenario B:", n, ve, "i =", i, "\n")
      # ---------------------------
      # Observed data
      # ---------------------------
      group <- c(rep(0, n), rep(1, n))
      hr_vector <- ifelse(group == 0, 1, hr)
      
      U <- runif(2*n)
      
      T_obs <- scale_true * (-log(U))^(1/shape_true) /
        (hr_vector)^(1/shape_true)
      
      C_obs <- rexp(2*n, rate = cens_rate)
      
      time_obs  <- pmin(T_obs, C_obs, ftime)
      event_obs <- as.integer(T_obs <= C_obs & T_obs <= ftime)
      
      # ---------------------------
      # Noise levels
      # ---------------------------
      noise_levels <- seq(0, 2, by = 0.2) 
      
      all_pred <- list()
      
      for (k in seq_along(noise_levels)) {
        
        noise_sd <- noise_levels[k]
        
        # ---------------------------
        # Predicted (noisy placebo)
        # ---------------------------
        #shape_pred <- shape_true + rnorm(2*n, 0, noise_sd)
        #shape_pred <- pmax(shape_pred, 0.2)
        
       # scale_pred <- get_scale(fc, ftime, shape_pred)
        
        T_pred <- scale_true * (-log(U))^(1/shape_true)
        T_pred_noisy <- T_pred * exp(rnorm(2*n, 0, noise_sd))
        
        C_pred <- rexp(2*n, rate = cens_rate)
        
        time_pred  <- pmin(T_pred_noisy, C_pred, ftime)
        event_pred <- as.integer(T_pred_noisy <= C_pred & T_pred_noisy <= ftime)
        
        # Store
        all_pred[[k]] <- data.frame(
          id = 1:(2*n),
          group = group,
          
          time_obs = time_obs,
          event_obs = event_obs,
          
          time_pred = time_pred,
          event_pred = event_pred,
          
          noise = noise_sd
        )
      }
      
      # Combine all noise levels
      sim_df <- do.call(rbind, all_pred)
      
      fname <- sprintf("%s/sim_%04d.csv", folder_path, i)
      write.csv(sim_df, fname, row.names = FALSE)
    }
  }
}

# ===========================================================
# Scenario C — vary red censor rate, fixed fc, fixed VE
# ===========================================================
# fixed_fc <- 0.4
# fixed_ve <- 30

# for (n in n_values) {
  
#   scale_true <- get_scale(fixed_fc, ftime, shape_true)
  
#   for (cens_rate_pred in cens_rate_pred_values) {
    
#     hr <- 1 - fixed_ve/100
    
#     folder_name <- sprintf("sim_n%d_fc%.1f_ve%d_censpred%.2f", n, fixed_fc, fixed_ve, cens_rate_pred)
#     folder_path <- file.path(root_dir, folder_name)
#     if (!dir.exists(folder_path)) dir.create(folder_path)
    
#     for (i in 1:nsim) {
#       if (i %% 1000 == 0) cat("Scenario C:", n, cens_rate_pred, "i =", i, "\n")
#       # ---------------------------
#       # Observed data
#       # ---------------------------
#       group <- c(rep(0, n), rep(1, n))
#       hr_vector <- ifelse(group == 0, 1, hr)
      
#       U <- runif(2*n)
      
#       T_obs <- scale_true * (-log(U))^(1/shape_true) /
#         (hr_vector)^(1/shape_true)
      
#       C_obs <- rexp(2*n, rate = cens_rate)
      
#       time_obs  <- pmin(T_obs, C_obs, ftime)
#       event_obs <- as.integer(T_obs <= C_obs & T_obs <= ftime)
      
#       # ---------------------------
#       # Noise levels
#       # ---------------------------
#       noise_levels <- seq(0, 2, by = 0.2) 
      
#       all_pred <- list()
      
#       for (k in seq_along(noise_levels)) {
        
#         noise_sd <- noise_levels[k]
        
#         # ---------------------------
#         # Predicted (noisy placebo)
#         # ---------------------------
  
#         T_pred <- scale_true * (-log(U))^(1/shape_true)
#         T_pred_noisy <- T_pred * exp(rnorm(2*n, 0, noise_sd))
        
#         C_pred <- rexp(2*n, rate = cens_rate_pred)
        
#         time_pred  <- pmin(T_pred_noisy, C_pred, ftime)
#         event_pred <- as.integer(T_pred_noisy <= C_pred & T_pred_noisy <= ftime)
        
#         # Store
#         all_pred[[k]] <- data.frame(
#           id = 1:(2*n),
#           group = group,
          
#           time_obs = time_obs,
#           event_obs = event_obs,
          
#           time_pred = time_pred,
#           event_pred = event_pred,
          
#           noise = noise_sd
#         )
#       }
      
#       # Combine all noise levels
#       sim_df <- do.call(rbind, all_pred)
      
#       fname <- sprintf("%s/sim_%04d.csv", folder_path, i)
#       write.csv(sim_df, fname, row.names = FALSE)
#     }
#   }
# }

cat("✔ All Weibull simulations with correlated predictions created.\n")


