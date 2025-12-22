library(survobj)
library(survival)
set.seed(12345)

# ---------------------------
# Simulation parameters
# ---------------------------
nsim <- 1000
ftime <- 1.5

n_values  <- c(30, 100, 200, 1000)
ve_values <- c(5, 10, 20, 30,40,50,60,70,80, 90)
fc_values <- seq(0.1, 0.9, 0.1)

root_dir <- "sim_ct_df"
if (!dir.exists(root_dir)) dir.create(root_dir)


# -----------------------------------------------------------
# Scenario A — vary event probability in the control group (fc), fixed ve = 20
# -----------------------------------------------------------
fixed_ve <- 20
fixed_hr <- 1 - fixed_ve/100

for (n in n_values) {
  for (fc in fc_values) {
    
    # Folder for this parameter combination
    folder_name <- sprintf("sim_n%d_fc%.1f_ve%d", n, fc, fixed_ve)
    folder_path <- file.path(root_dir, folder_name)
    if (!dir.exists(folder_path)) dir.create(folder_path)
    
    # Distribution for control group
    s_events <- s_exponential(fail = fc, t = ftime)
    s_events_pred <- s_exponential(fail = fc, t = ftime)
    
    for (i in 1:nsim) {
      
      group <- c(rep(0, n), rep(1, n))
      hr_vector <- ifelse(group == 0, 1, fixed_hr)
      
      sim_time <- s_events$rsurvhr(hr_vector)
      cevent   <- censor_event(censor_time = ftime, time = sim_time, event = 1)
      ctime    <- censor_time(censor_time = ftime, time = sim_time)
      
      # Counterfactual time from control distribution
      counterfactual_time_full <- s_events$rsurvhr(rep(1, 2*n))   # hr=1 → untreated
      counter_time  <- censor_time(censor_time = ftime, time = counterfactual_time_full)
      counter_event <- censor_event(censor_time = ftime, time = counterfactual_time_full, event = 1)
      
      # Keep counterfactual values only for treated; NA for controls
      counter_time[group == 0]  <- NA
      counter_event[group == 0] <- NA
      
      sim_df <- data.frame(
        time = ctime,
        event = cevent,
        group = group,
        counter_time = counter_time,
        counter_event = counter_event
      )
      # File name
      fname <- sprintf("%s/sim_%04d.csv", folder_path, i)
      write.csv(sim_df, fname, row.names = FALSE)
    }
  }
}

# -----------------------------------------------------------
# Scenario B — vary VE, fixed fail_control = 0.4
# -----------------------------------------------------------
fixed_fc <- 0.4

for (n in n_values) {
  
  s_events <- s_exponential(fail = fixed_fc, t = ftime)
  
  for (ve in ve_values) {
    hr <- 1 - ve/100
    
    # Folder for this parameter combination
    folder_name <- sprintf("sim_n%d_fc%.1f_ve%d", n, fixed_fc, ve)
    folder_path <- file.path(root_dir, folder_name)
    if (!dir.exists(folder_path)) dir.create(folder_path)
    
    for (i in 1:nsim) {
      
      group <- c(rep(0, n), rep(1, n))
      hr_vector <- ifelse(group == 0, 1, hr)
      
      sim_time <- s_events$rsurvhr(hr_vector)
      cevent   <- censor_event(censor_time = ftime, time = sim_time, event = 1)
      ctime    <- censor_time(censor_time = ftime, time = sim_time)
      
      # Counterfactual time from control distribution
      counterfactual_time_full <- s_events$rsurvhr(rep(1, 2*n))   # hr=1 → untreated
      counter_time  <- censor_time(censor_time = ftime, time = counterfactual_time_full)
      counter_event <- censor_event(censor_time = ftime, time = counterfactual_time_full, event = 1)
      
      # Keep counterfactual values only for treated; NA for controls
      counter_time[group == 0]  <- NA
      counter_event[group == 0] <- NA
      
      sim_df <- data.frame(
        time = ctime,
        event = cevent,
        group = group,
        counter_time = counter_time,
        counter_event = counter_event
      )
      
      fname <- sprintf("%s/sim_%04d.csv", folder_path, i)
      write.csv(sim_df, fname, row.names = FALSE)
    }
  }
}

cat("✔ All simulation folders and CSV files successfully created.\n")