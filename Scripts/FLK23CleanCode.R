##############################################################################
# SKUA CASE STUDY -- CONSIDERATIONS FOR DISEASE/MOVEMENT ECOLOGY STUDY #
# Falkland Islands/Islas Malvinas Initial Analyses (Cleaning and HMM) #
# Last updated: 22Mar26 BS #
# R version 4.5.2 (2025-10-31) #

##############################################################################

# load packages
library(dplyr)
library(data.table)
library(readr)
library(readxl)
library(stringr)
library(forcats)
library(lubridate)
library(scales)
library(ggplot2)
library(ggpubr)
library(viridis)
library(purrr)

library(sf)
library(raster)
library(ggspatial)
library(rnaturalearth)

library(amt)
library(momentuHMM)
library(diveMove)
library(hmmTMB)
library(adehabitatHR)
library(bayesplot)
library(rstan)
library(posterior)

# call data
skua_fks23_gps <- read_csv("Raw Data/FLK23_Raw.csv")

# preliminary cleaning, gps data...
  ## create date-time column
  skua_fks23_gps$year <- paste0("20", skua_fks23_gps$year)
  
  skua_fks23_gps$datetime <- as.POSIXct(paste(skua_fks23_gps$year, skua_fks23_gps$month, skua_fks23_gps$day, skua_fks23_gps$hour, skua_fks23_gps$min, skua_fks23_gps$sec), format = "%Y %m %d %H %M %S", tz = "UTC")
  
  ## remove lat and long at 0.0000 (failed/uncalibrated)
  skua_fks23_gps <- skua_fks23_gps[skua_fks23_gps$latitude != 0,]
  
  ## remove any duplicates or near duplicates (occured within 10 s of another fix)
  d1new <- skua_fks23_gps %>%
    group_by(ID) %>%
    do(distinct(., datetime, .keep_all = TRUE)) %>%
    do(mutate(
      .,
      dup = difftime(datetime, lag(datetime), units = "secs") < 10)) %>%
    do(arrange(., order(datetime))) %>%
    dplyr:: filter(.,!dup) %>%
    dplyr:: select(-dup) %>%
    arrange(datetime)
  
  d1.rep.new <- nrow(skua_fks23_gps) - nrow(d1new)
  cat(sprintf("%d duplicate time &/or near duplicate location/time records removed\n", d1.rep.new)) 

  ## remove extreme travel rate locations (30 m/s) (consistent across literature)...
  vmax <- 30
  
  d5new <- d1new %>%
    do(mutate(., keep = grpSpeedFilter(cbind(datetime, longitude, latitude), speed.thr = vmax)))
  
  d5.rep.new <- nrow(d1new) - sum(d5new$keep)
  cat(sprintf(
    paste(
      "\n%d records with travel rates >",
      vmax,
      "m/s will be ignored by SSM filter\n"
    ),
    d5.rep.new
  )) 
  
  NewTracks <- d5new %>%
    filter(keep == "TRUE") %>%
    mutate(lc = "G") %>%
    ungroup()
  
  ## add capture type
  NewTracks$breeding_status <- NA
  allids <- unique(NewTracks$ID)

    onnest <- c("912", "916", "920", "927", "928", "929", "933", "934", "935", "941", "952", "954", "958")
    offnest <- setdiff(allids, onnest)
    
    for(i in 1:nrow(NewTracks)) {
      if(NewTracks$ID[i] %in% onnest) {
        NewTracks$breeding_status[i] <- "On-nest"
      }
      else if (NewTracks$ID[i] == 915) {
        NewTracks$breeding_status[i] <- "Off-nest?"
      }
      else {
        NewTracks$breeding_status[i] <- "Off-nest"
      }
    }

  ## save as new dataset, format datetime
  cleanedbreed <- NewTracks
  
  cleanedbreed$datetime <- as.POSIXct(cleanedbreed$datetime, format = "%Y %m %d %H %M %S", tz = "UTC")

  ## remove 915 (many unrealistic points)
  cleanedbreed <- cleanedbreed[cleanedbreed$ID != 915,]
  
  ## remove 912, 921, 931 with too few points
  cleanedbreed <- cleanedbreed[cleanedbreed$ID != 912,]
  cleanedbreed <- cleanedbreed[cleanedbreed$ID != 921,]
  cleanedbreed <- cleanedbreed[cleanedbreed$ID != 931,]    

  ## make subsets by island, add as column
  new <- c("930", "932", "933", "934", "935", "941", "942", "944", "956", "959")
  saunders <- c("911", "912", "914", "918", "927", "928", "929")
  bleaker <- setdiff(allids, union(new, saunders))
  
  cleanedbreed <- cleanedbreed %>%
    mutate(
      island = case_when(
        ID %in% new ~ "New Island",
        ID %in% saunders ~ "Saunders",
        TRUE ~ "Bleaker"
      ))  
  
  # count points per individual
  skua_wdl24_argos %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(count = n())
  
  # build GPS heatmap summarizing fixes 
  cleanedbreed$date <- as.Date(cleanedbreed$datetime)
  ARGOSdata_recording <- as.data.frame(table(cleanedbreed$ID, cleanedbreed$date))
  colnames(ARGOSdata_recording) <- c("ID", "date", "points")
  ARGOSdata_recording$date <- as.Date(ARGOSdata_recording$date)
  
  ARGOSfull_data <- expand_grid(
    date = seq(min(ARGOSdata_recording$date), max(ARGOSdata_recording$date), by = "day"),
    ID = unique(ARGOSdata_recording$ID)
  ) %>%
    left_join(ARGOSdata_recording, by = c("date", "ID")) %>%
    mutate(points = replace_na(points, 0))  # Fill NAs with 0
  
  heatmapFLK <- ggplot(data = ARGOSfull_data, aes(x = date, y = ID, fill = points)) +
    geom_tile() +
    scale_fill_viridis(limits = c(0, max(ARGOSdata_recording$points)), option = "D") +
    scale_x_date(date_breaks = "5 days", date_labels = "%b %d") +
    xlab("Date") + ylab("ID") +
    labs(fill = "GPS fixes per day") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

######################################################################################
# compute Hidden Markov Model (HMM, uses package 'momentuHMM')
  
  data_hmm <- cleanedbreed %>% 
    dplyr::rename(lon=longitude) %>%
    dplyr::rename(lat=latitude) %>%
    dplyr::rename(time=datetime) %>%
    mutate(month = substr(time, 6, 7))
  
  data_hmm <- data_hmm %>%
    dplyr::select(ID, lon, lat, time, breeding_status, island)
  
  # Check for individuals with too little data
  table(data_hmm$ID)

  ## project to UTM
  llcoord <- st_as_sf(data_hmm[, c("lon", "lat")], coords = c("lon", "lat"), crs = CRS("+proj=longlat +datum=WGS84"))
  
  utmcoord <- st_transform(llcoord, 32721)
  
  # extract UTM coordinates back to the dataframe
  utm_coords <- st_coordinates(utmcoord)
  data_hmm$lon <- utm_coords[,1]
  data_hmm$lat <- utm_coords[,2]
  
  ## table of time intervals in data
  plot(table(diff(data_hmm$time)),
       xlim = c(0, 1000), 
       xlab = "Time interval (seconds)", 
       ylab = "Count")
  
  data_hmm <- data_hmm %>%
    arrange(ID, time)
  
  ## prepare data for HMM (OPTIONAL: correlated random walk; REQUIRED: compute step lengths and turning angles)
  data_hmmCRWL <- crawlWrap(obsData = data_hmm, 
                            timeStep = "1 hour", 
                            Time.name = "time", 
                            coord = c("lon", "lat"),
                            attempts = 10000)
  
  ## retrieve step lengths and turning angles
  data_hmmPrepped <- momentuHMM::prepData(data_hmmCRWL, type = "UTM")
  
        # Build GPS heatmap summarizing regularized fixes
        data_hmmPrepped$date <- as.Date(data_hmmPrepped$time)
        ARGOSdata_recording <- as.data.frame(table(data_hmmPrepped$ID, data_hmmPrepped$date))
        colnames(ARGOSdata_recording) <- c("ID", "date", "points")
        ARGOSdata_recording$date <- as.Date(ARGOSdata_recording$date)
        
        ARGOSfull_data <- expand_grid(
          date = seq(min(ARGOSdata_recording$date), max(ARGOSdata_recording$date), by = "day"),
          ID = unique(ARGOSdata_recording$ID)
        ) %>%
          left_join(ARGOSdata_recording, by = c("date", "ID")) %>%
          mutate(points = replace_na(points, 0))  # Fill NAs with 0
        
        heatmapFLKregularized <- ggplot(data = ARGOSfull_data, aes(x = date, y = ID, fill = points)) +
          geom_tile() +
          scale_fill_viridis(limits = c(0, max(ARGOSdata_recording$points)), option = "D") +
          scale_x_date(date_breaks = "5 days", date_labels = "%b %d") +
          xlab("Date") + ylab("ID") +
          labs(fill = "GPS fixes per day") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        
  ## remove NAs for turning angle and step length, fix format and covariates lost in prepping
  data_hmmPrepped <- data_hmmPrepped %>%
    dplyr::filter(!is.na(step) & !is.na(angle))
  
    summary(data_hmmPrepped$step)
  
    ### checking the step length distribution
    hist(data_hmmPrepped$step, xlab = "step length", main = "") # heavy tail... 
    
    ### checking the turning angle distribution
    hist(data_hmmPrepped$angle, breaks = seq(-pi, pi, length = 15), xlab = "angle", main = "")
  
    ## fix capture type and island
    data_hmmPrepped$breeding_status <- NA
    
    for(i in 1:nrow(data_hmmPrepped)) {
      if(data_hmmPrepped$ID[i] %in% onnest) {
        data_hmmPrepped$breeding_status[i] <- "On-nest"
      }
      else if (data_hmmPrepped$ID[i] == 915) {
        data_hmmPrepped$breeding_status[i] <- "Off-nest?"
      }
      else {
        data_hmmPrepped$breeding_status[i] <- "Off-nest"
      }
    }  
    
    data_hmmPrepped$breeding_status <- factor(data_hmmPrepped$breeding_status,
                                              levels = c("On-nest", "Off-nest"))
    
    data_hmmPrepped <- data_hmmPrepped %>%
      mutate(
        island = case_when(
          ID %in% new ~ "New Island",
          ID %in% saunders ~ "Saunders",
          TRUE ~ "Bleaker"
        ))  
    
    data_hmmPrepped$island <- factor(data_hmmPrepped$island,
                                              levels = c("New Island", "Saunders", "Bleaker"))
    
    # make day of season column (start at 1)
    season_start <- as.Date("2023-12-06")
    
    data_hmmPrepped <- data_hmmPrepped %>%
      mutate(
        date = as.Date(time),
        day_of_season = as.integer(date - season_start) + 1
      )
    
    # check ACF
    acf(data_hmmPrepped$step[!is.na(data_hmmPrepped$step)],lag.max=100) 
    
  # initial values for HMM
    n_states = 3
    step_mean0 = c(10, 40, 500) 
    step_sd0<-c(5, 20, 250) 
    angle_mean0<- c(0, 0, 0) 
    angle_rho0 <- c(0.3,0.5,0.7) # angle concentration
  
  Par0<-list(step = list(mean = step_mean0, sd = step_sd0), 
             angle = list(mu = angle_mean0, rho = angle_rho0))
  
  # observation process
  obs <- Observation$new(data_hmmPrepped,
                       dists = list(step="gamma2",angle="wrpcauchy"),
                       par = Par0,
                       n_states = n_states)
  
  # transition matrix
  hid <- MarkovChain$new(n_states = n_states,
                         formula = ~breeding_status * island + day_of_season + s(ID, bs = "re"),
                       data = data_hmmPrepped,
                       initial_state = "stationary")
  
  hmm1<-HMM$new(obs=obs,hid=hid)
  
  
  # specification of priors...
      
      hmm1$priors()
      
      # OBSERVATION MODEL
        #mean step length...
            plot(density(exp(rnorm(1e4, mean = log(10), sd = 0.1))),
                 main = "state 1 mu SD = 0.1", xlab = NA)
            plot(density(exp(rnorm(1e4, mean = log(10), sd = 0.5))),
                 main = "state 1 mu SD = 0.5", xlab = NA)
            
            plot(density(exp(rnorm(1e4, mean = log(40), sd = 0.1))),
                 main = "state 2 mu SD = 0.1", xlab = NA)
            plot(density(exp(rnorm(1e4, mean = log(40), sd = 0.5))),
                 main = "state 2 mu SD = 0.5", xlab = NA)
            
            plot(density(exp(rnorm(1e4, mean = log(500), sd = 0.1))),
                 main = "state 3 mu SD = 0.1", xlab = NA)
            plot(density(exp(rnorm(1e4, mean = log(500), sd = 0.5))),
                 main = "state 3 mu SD = 0.5", xlab = NA)
            
        # specify based on observed patterns and initial values that worked in frequentist model, but keeping it loose/uninformed
            prior_obs <- matrix(c(log(10), 0.5,
                                  log(50), 1,
                                  log(500), 1.5,
                                  log(5), 1,
                                  log(25), 1,
                                  log(250), 1,
                                  0, 0.5,
                                  0, 0.5,
                                  0, 0.5,
                                  qlogis(0.8), 1,
                                  qlogis(0.8), 1,
                                  qlogis(0.8), 1),
                                ncol = 2, byrow = TRUE)
      
      # HIDDEN MODEL (transition probabilities)
            prior_hid <- matrix(c(-2, 1,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2, 
                                  -2, 1,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  -2, 1,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  -2, 1,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  -2, 1,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  -2, 1,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2,
                                  0, 2),
                                ncol = 2, byrow = TRUE)
      
        # HIDDEN LAMBDA/RE PRIORS
              prior_lambda <- matrix(c(
                0, 1,  
                0, 1,  
                0, 1,  
                0, 1,  
                0, 1,  
                0, 1),
                ncol = 2, byrow = TRUE)
        
      # update the priors...
         hmm1$set_priors(new_priors = list(coeff_fe_obs = prior_obs,
                                           coeff_fe_hid = prior_hid,
                                           log_lambda_hid = prior_lambda))
     
    # check that they're set!
     hmm1$priors()
      
  # fit the model...
  hmm1$fit_stan(chains = 4, iter = 2000)
  
  # check output...
  hmm1$out_stan()
  
  hmm1$plot_ts("x", "y") +
    geom_point(size = 0.5) + 
    labs(x = "Easting (UTM Zone 21S)", y = "Northing (UTM Zone 21S)") + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),  # tilt x-axis labels
      axis.text.y = element_text(angle = 45, hjust = 1)   # tilt y-axis labels
    ) +
    scale_x_continuous(labels = scales::comma_format()) +  # full numbers with commas
    scale_y_continuous(labels = scales::comma_format())  
    
  hmmstep <- hmm1$plot_dist("step")
  hmmangle <- hmm1$plot_dist("angle")
  
  traceplot(hmm1$out_stan(),pars = "coeff_fe_hid")
  stan_dens(hmm1$out_stan(),pars = "coeff_fe_hid")
  
  hmm1$plot("delta", var = "breeding_status")
  
  hmm1$par()
  hmm1$coeff_fe()
  hmm1$coeff_re()
  
  data_hmmPrepped$state<-hmm1$viterbi()
  data_hmmPrepped$state<-as.factor(data_hmmPrepped$state)
  data_hmmPrepped <- data_hmmPrepped %>%
    mutate(behavior = factor(state, 
                             levels = 1:3,
                             labels = c("Resting", "Foraging", "Commuting")))

  save(hmm1, file = "HMMfull3state.RData")
  
  # get iterations for visualization
  iters <- hmm1$iters()
  iters <- as.data.frame.table(iters)
  
  # now, visualize posteriors for transition probabilities and observation parameters
  # create a mapping for better parameter labels
  parameter_labels <- c(
    # step parameters
    "step.mean.state1" = "Step Mean\nResting",
    "step.mean.state2" = "Step Mean\nForaging", 
    "step.mean.state3" = "Step Mean\nCommuting",
    "step.sd.state1" = "Step SD\nResting",
    "step.sd.state2" = "Step SD\nForaging",
    "step.sd.state3" = "Step SD\nCommuting",
    
    # angle parameters
    "angle.mu.state1" = "Angle Mean\nResting",
    "angle.mu.state2" = "Angle Mean\nForaging",
    "angle.mu.state3" = "Angle Mean\nCommuting",
    "angle.rho.state1" = "Angle Concentration\nResting",
    "angle.rho.state2" = "Angle Concentration\nForaging",
    "angle.rho.state3" = "Angle Concentration\nCommuting",
    
    # transition probabilities
    "S1>S1" = "Resting → Resting",
    "S1>S2" = "Resting → Foraging",
    "S1>S3" = "Resting → Commuting",
    "S2>S1" = "Foraging → Resting",
    "S2>S2" = "Foraging → Foraging",
    "S2>S3" = "Foraging → Commuting",
    "S3>S1" = "Commuting → Resting",
    "S3>S2" = "Commuting → Foraging",
    "S3>S3" = "Commuting → Commuting"
  )
  
  # apply the labels
  iters$Parameter_Label <- parameter_labels[as.character(iters$Var2)]
  
  # create separate plots for movement parameters and transition probabilities
  # movement parameters (step and angle)
  movement_params <- iters[grepl("step|angle", iters$Var2), ]
  
  movement_plot <- ggplot(movement_params, aes(x = Freq)) +
    geom_histogram(bins = 30, fill = "lightblue", col = "white", alpha = 0.8) +
    facet_wrap(~ Parameter_Label, nrow = 3, scales = "free") +
    labs(x = "Parameter Value", 
         y = "Posterior Density") +
    theme_bw() +
    theme(strip.text = element_text(size = 9),
          axis.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 12))
  
  # transition probabilities
  transition_params <- iters[grepl("S[1-3]>S[1-3]", iters$Var2), ]
  
  transition_plot <- ggplot(transition_params, aes(x = Freq)) +
    geom_histogram(bins = 30, fill = "lightcoral", col = "white", alpha = 0.8) +
    facet_wrap(~ Parameter_Label, nrow = 3, scales = "free_y") +
    labs(x = "Transition Probability", 
         y = "Posterior Density") +
    theme_bw() +
    theme(strip.text = element_text(size = 9),
          axis.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 12))
  
  # display the plots
  print(movement_plot)
  print(transition_plot)
  
  # save
  ggsave("obs_posteriorsHMM.png", dpi = 600, bg = "white", units = "in", height = 5, width = 7)
  ggsave("tpm_posteriorsHMM.png", dpi = 600, bg = "white", units = "in", height = 4, width = 6)
  
  # visualize effects of covariates on transition probabilites
  draws<-as.array(hmm1$out_stan(),pars = "coeff_fe_hid")
  draws <- as_draws_df(draws)
  summary_df<-summarise_draws(
    draws,
    mean = ~ mean(.x),           
    median =~ median(.x),
    sd = ~ sd(.x),
    lower = ~ quantile(.x, probs = 0.025) ,
    upper = ~ quantile(.x, probs = 0.975) # 95% credible interval
  )%>%
    mutate(transition = rep(c("Pr (1->2)","Pr (1->3)","Pr (2->1)","Pr (2->3)","Pr (3->1)","Pr (3->2)"), each=7),
           parameter = rep(c("Intercept","Off-nest","Saunders","Bleaker","Day of season","Off-nest x Saunders","Off-nest x Bleaker"),6))%>%
    mutate(parameter = factor(parameter, levels = c("Off-nest x Saunders", "Off-nest x Bleaker", "Day of season","Saunders",
                                                    "Bleaker","Off-nest","Intercept")))
  
  ggplot(summary_df, aes(x = parameter, y = mean)) +
    geom_pointrange(aes(ymin = `2.5%`, ymax = `97.5%`)) +
    geom_hline(yintercept = 0,colour ="red",linetype =2)+
    facet_wrap(~factor(transition),scales = "free_x")+
    theme_bw() +
    ylab("Mean") +
    xlab("Parameter") +
    coord_flip()
  
  ggsave("tpmMedError.png", dpi = 600, bg = "white", units = "in", height = 4, width = 6)
  
  # to inspect patterns by covariates
  
      # create grid of ordinal day
      newdata <- data.frame(day_of_season = seq(min(data_hmmPrepped$day_of_season), max(data_hmmPrepped$day_of_season), length = 100),
                            island = "Bleaker",
                            breeding_status = "Off-nest",
                            ID = 911)
      
      # OR create grid of island
      newdata <- data.frame(day_of_season = 0,
                            island = as.factor(unique(data_hmmPrepped$island)),
                            breeding_status = "Off-nest",
                            ID = 911)
      
      # OR create grid of breeding status
      newdata <- data.frame(day_of_season = 0,
                            island = "Bleaker",
                            breeding_status = as.factor(unique(data_hmmPrepped$breeding_status)),
                            ID = 911)
      
      ## THEN randomly resample
      ind_post <- sort(sample(1:1000, size = 100))
      
      # loop over randomly selected posterior samples
      probs_list <- lapply(ind_post, function(i_post) {
        
        hmm1$update_par(iter = i_post)
        
        tpm_post <- hmm1$predict(what = "tpm", newdata = newdata)
        
        n_states <- dim(tpm_post)[1]
        n_cov    <- dim(tpm_post)[3]
        
        transitions <- which(row(tpm_post[,,1]) != col(tpm_post[,,1]), arr.ind = TRUE)
        
        probs <- do.call(rbind, lapply(seq_len(nrow(transitions)), function(i){
          
          from <- transitions[i,1]
          to   <- transitions[i,2]
          
          data.frame(
            prob = tpm_post[from, to, ],
            transition = paste0("Pr(", from, " -> ", to, ")")
          )
        }))
        
        probs$group <- paste0("iter ", i_post, " - ", probs$transition)
        
        return(probs)
      }) 
      
      
      # ... and combine data frames for plot
      
      # for each island...
      probs_df <- do.call(what = rbind, args = probs_list)
      probs_df$day_of_season <- newdata$day_of_season
      probs_df$island<-newdata$island
      island_plot <- ggplot(probs_df, aes(island, prob, group = group)) +
        geom_point(alpha = 0.5) +
        labs(x = "Island", y = "Transition Probabilities") +
        facet_wrap("transition") +theme_classic()
      
      # OR by capture type...
      probs_df <- do.call(what = rbind, args = probs_list)
      probs_df$day_of_season <- newdata$day_of_season
      probs_df$breeding_status<-newdata$breeding_status
      captureplot <- ggplot(probs_df, aes(breeding_status, prob, group = group)) +
        geom_point(alpha = 0.5) +
        labs(x = "Capture Type", y = "Transition Probabilities") +
        facet_wrap("transition") +theme_classic()
      
      # OR by day of season
      probs_df <- do.call(what = rbind, args = probs_list)
      probs_df$day_of_season <- newdata$day_of_season
      dayplot <- ggplot(probs_df, aes(day_of_season, prob, group = group)) +
        geom_line(linewidth = 0.1, alpha = 0.5) +
        labs(x = "Day of Season", y = "Transition Probabilities") +
        facet_wrap("transition") +theme_classic()
      
      ### for stationary probabilities... use newdata above then...
      
      probs_df <- data.frame()
      
      for(i in ind_post) {
        hmm1$update_par(iter = i)
        probs <- data.frame(state = rep(paste0("state ", 1:3), each = 100),
                            prob = as.vector(hmm1$predict(what = "delta",
                                                          newdata = newdata)))
        probs$group <- paste0("iter ", i, " - ", probs$state)
        probs_df <- rbind(probs_df, probs)
      }
      
      probs_df$day_of_season <- newdata$day_of_season
      probs_df$breeding_status <- newdata$breeding_status
      probs_df$island <- newdata$island
      
      # plot stationary state probs against day of season
      probs_day <- ggplot(probs_df, aes(day_of_season, prob, group = group, col = state)) +
        geom_line(size = 0.3, alpha = 0.5) +
        labs(x = "time of day", y = "stationary state probabilities", col = NULL) +
        scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
        guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1)))
      
      ggsave("probDay.png", probs_day, dpi = 600, bg = "white", units = "in", height = 7, width = 10)
      
      # plot stationary state probs against capture type...
      probs_cap <- ggplot(probs_df, aes(breeding_status, prob, group = group, col = state)) +
        geom_line(size = 0.1, alpha = 0.5) +
        labs(x = "capture type", y = "stationary state probabilities", col = NULL) +
        scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
        guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1)))
      
      # plot stationary state probs against island...
      ggplot(probs_df, aes(island, prob, group = group, col = state)) +
        geom_line(size = 0.1, alpha = 0.5) +
        labs(x = "island", y = "stationary state probabilities", col = NULL) +
        scale_color_manual(values = hmmTMB:::hmmTMB_cols) +
        guides(color = guide_legend(override.aes = list(size = 0.5, alpha = 1)))

  # after viterbi decoding... map !
      
      # convert back to Lat/Lon
        # create sf object from x, y coordinates (UTM 21S)
        utm_points <- data_hmmPrepped %>%
          filter(!is.na(x) & !is.na(y)) %>%
          st_as_sf(coords = c("x", "y"), crs = 32721)
        
        # transform to WGS84 (lat/lon degrees)
        wgs84_points <- st_transform(utm_points, crs = 4326)
      
        # extract coordinates
        coords_latlon <- st_coordinates(wgs84_points)
        
        # fill the lat/lon columns in your original dataframe
        valid_rows <- which(!is.na(data_hmmPrepped$x) & !is.na(data_hmmPrepped$y))
        data_hmmPrepped$lon[valid_rows] <- coords_latlon[, "X"]  # Longitude
        data_hmmPrepped$lat[valid_rows] <- coords_latlon[, "Y"]  # Latitude
        
        head(data_hmmPrepped[, c("x", "y", "lon", "lat")])
      
####################
    # plot with osm

        library(osmdata)
        
        # Query coastline polygons for Falklands bounding box
        bbox <- c(-62, -55, -55, -50) # xmin, ymin, xmax, ymax
        falklands_osm <- opq(bbox) %>% add_osm_feature(key = 'natural', value = 'coastline') %>% 
          osmdata_sf()
        
        islands_poly <- falklands_osm$osm_polygons
        main_coast_lines <- falklands_osm$osm_lines
        
        coast_union <- st_union(main_coast_lines)
        mainland_poly <- st_polygonize(coast_union)
        
        all_land <- bind_rows(
          st_sf(geometry = mainland_poly),
          islands_poly
        )
        
        tracks_sf <- data_hmmPrepped %>%
          st_as_sf(coords = c("lon", "lat"), crs = 4326) # WGS84 CRS for lon/lat
        
        tracks_lines <- tracks_sf %>%
          group_by(ID) %>%
          summarise(do_union = FALSE) %>%               # keep trajectories per ID
          st_cast("LINESTRING")
        
        # define the bounding boxes for the zoom areas
        zoom_boxes <- tibble::tribble(
          ~name,     ~xmin,   ~xmax,   ~ymin,   ~ymax,
          "Saunders", -60.5,  -60.1,  -51.4,  -51.2,
          "New",      -61.8,  -60.8,  -52.4,  -51.4,
          "Bleaker",  -59,    -58.7,  -52.3,  -52.1
        )
        
        # convert to sf polygons
        zoom_boxes_sf <- zoom_boxes %>%
          rowwise() %>%
          mutate(geometry = list(
            st_polygon(list(matrix(
              c(xmin, ymin,
                xmin, ymax,
                xmax, ymax,
                xmax, ymin,
                xmin, ymin),
              ncol = 2,
              byrow = TRUE
            )))
          )) %>%
          st_as_sf(crs = 4326)
        
        # whole plot...
        ggplot() +
          geom_sf(data = all_land, fill = "grey80", colour = "black") +
          geom_sf(data = tracks_lines, aes(colour = ID), linewidth = 0.5, alpha = 0.8) +
          geom_sf(data = zoom_boxes_sf, fill = NA, colour = "black", linewidth = 0.5) +
          scale_colour_viridis_d(option = "turbo") +
          coord_sf(xlim = c(-62, -55), ylim = c(-53.5, -50.5), expand = FALSE) +
          theme_void() +
          labs(
            colour = "Individual"
          )
          
          ggsave("FALKwholeMAPtracks.png", dpi = 600, units = "in", width = 10, height = 6, bg = "white")
        
          # saunders...
          ggplot() +
            geom_sf(data = all_land, fill = "grey80", colour = "black") +
            geom_sf(data = tracks_lines, aes(colour = ID), linewidth = 0.4, alpha = 0.8) +
            scale_colour_viridis_d(option = "turbo") +
            coord_sf(xlim = c(-60.5, -60.1), ylim = c(-51.4, -51.2), expand = FALSE) +
            theme_classic() +
            theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) + 
            labs(
              colour = "Individual"
            ) 
          
          ggsave("SAUwholeMAPtracks.png", dpi = 600, units = "in", width = 6, height = 6, bg = "white")
          
          # new...
          ggplot() +
            geom_sf(data = all_land, fill = "grey80", colour = "black") +
            geom_sf(data = tracks_lines, aes(colour = ID), linewidth = 0.4, alpha = 0.8) +
            scale_colour_viridis_d(option = "turbo") +
            coord_sf(xlim = c(-61.8, -60.8), ylim = c(-52.4, -51.4), expand = FALSE) +
            theme_classic() +
            theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) + 
            labs(
              colour = "Individual"
            )
          
          ggsave("NEWwholeMAPtracks.png", dpi = 600, units = "in", width = 3, height = 5, bg = "white")
          
          # bleaker...
          ggplot() +
            geom_sf(data = all_land, fill = "grey80", colour = "black") +
            geom_sf(data = tracks_lines, aes(colour = ID), linewidth = 0.4, alpha = 0.8) +
            scale_colour_viridis_d(option = "turbo") +
            coord_sf(xlim = c(-59, -58.7), ylim = c(-52.3, -52.1), expand = FALSE) +
            theme_classic() +
            theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) + 
            labs(
              colour = "Individual"
            )
          
          ggsave("BLKwholeMAPtracks.png", dpi = 600, units = "in", width = 5, height = 5, bg = "white")
          

###################################################################################
# NON-OVERLAPPING 3-DAY MOVING WINDOW COMPUTATIONS WITH MULTIPLE UD LEVELS AND STATE-SPECIFIC UDs
      
      # create amt track
      trk <- mk_track(data_hmmPrepped, .x = x, .y = y, .t = time, id = ID, state = state)
      cat("Track time range:", as.character(range(trk$t_, na.rm = TRUE)), "\n")
      
      # function to create non-overlapping 3-day windows
      create_3day_windows <- function(start_date, end_date) {
        windows <- list()
        current_date <- start_date
        
        while (current_date <= end_date) {
          window_start <- as.POSIXct(paste(current_date, "00:00:00"), tz = "UTC")
          window_end <- as.POSIXct(paste(current_date + 2, "23:59:59"), tz = "UTC")  # 3 days total
          
          if (window_end <= as.POSIXct(paste(end_date, "23:59:59"), tz = "UTC")) {
            windows <- append(windows, list(list(
              start = window_start,
              end = window_end,
              center_date = current_date + 1  # middle day of the 3-day window
            )))
          }
          
          current_date <- current_date + 3  # move by 3 days for non-overlapping windows
        }
        
        return(windows)
      }
      
      # function to compute KDE and extract areas
      compute_kde_areas <- function(clean_track_data, id, center_date, window_start, window_end) {
        results <- list()
        
        tryCatch({
          # use clean track data directly for overall KDE
          kde_overall <- clean_track_data %>% hr_kde()
          
          # extract 50% and 95% isopleths, then compute areas
          hr_50 <- kde_overall %>% hr_isopleths(levels = 0.50)
          hr_95 <- kde_overall %>% hr_isopleths(levels = 0.95)
          
          area_50_overall <- as.numeric(st_area(hr_50)) / 1000000  # Convert to km²
          area_95_overall <- as.numeric(st_area(hr_95)) / 1000000
          
          results$area_50_overall_km2 <- area_50_overall
          results$area_95_overall_km2 <- area_95_overall
          
          cat("  Overall UDs: 50% =", round(area_50_overall, 3), "km², 95% =", round(area_95_overall, 3), "km²\n")
          
          # state-specific KDEs for all 3 states from the HMM
          for (state in c(1, 2, 3)) {
            state_data <- clean_track_data %>% filter(state == !!state)
            
            if (nrow(state_data) >= 5) {  # minimum points for KDE
              tryCatch({
                kde_state <- state_data %>% hr_kde()
                hr_50_state <- kde_state %>% hr_isopleths(levels = 0.50)
                hr_95_state <- kde_state %>% hr_isopleths(levels = 0.95)
                
                area_50_state <- as.numeric(st_area(hr_50_state)) / 1000000
                area_95_state <- as.numeric(st_area(hr_95_state)) / 1000000
                
                results[[paste0("area_50_state_", state, "_km2")]] <- area_50_state
                results[[paste0("area_95_state_", state, "_km2")]] <- area_95_state
                results[[paste0("n_points_state_", state)]] <- nrow(state_data)
                
                cat("  State", state, "UDs: 50% =", round(area_50_state, 3), "km², 95% =", round(area_95_state, 3), "km² (n =", nrow(state_data), ")\n")
                
              }, error = function(e) {
                cat("  Error computing state", state, "KDE:", e$message, "\n")
                results[[paste0("area_50_state_", state, "_km2")]] <- NA
                results[[paste0("area_95_state_", state, "_km2")]] <- NA
                results[[paste0("n_points_state_", state)]] <- nrow(state_data)
              })
            } else {
              cat("  Insufficient points for state", state, "KDE (n =", nrow(state_data), ")\n")
              results[[paste0("area_50_state_", state, "_km2")]] <- NA
              results[[paste0("area_95_state_", state, "_km2")]] <- NA
              results[[paste0("n_points_state_", state)]] <- nrow(state_data)
            }
          }
          
          # add metadata
          results$ID <- id
          results$center_date <- center_date
          results$window_start <- window_start
          results$window_end <- window_end
          results$success <- TRUE
          
          return(results)
          
        }, error = function(e) {
          cat("  Error computing overall KDE:", e$message, "\n")
          return(list(
            ID = id,
            center_date = center_date,
            window_start = window_start,
            window_end = window_end,
            success = FALSE,
            error_message = e$message
          ))
        })
      }
      
      # create date column
      data_hmmPrepped$date <- as.Date(data_hmmPrepped$time)
      unique_days <- unique(data_hmmPrepped$date)
      print(unique_days)
      
      # moving window computation...
      MW_results <- list()
      processed_count <- 0
      total_combinations <- 0
      
      # calculate total number of combinations (ID x windows)
      for (i in unique(trk$id)) {
        indiv_data <- trk %>% filter(id == i) %>% arrange(t_)  # Keep as track object
        if (nrow(indiv_data) > 0) {
          start_date <- as.Date(min(indiv_data$t_))
          end_date <- as.Date(max(indiv_data$t_))
          windows <- create_3day_windows(start_date, end_date)
          total_combinations <- total_combinations + length(windows)
        }
      }
      cat("Total ID-window combinations to process:", total_combinations, "\n")
      
      for (i in unique(trk$id)) {    # computing metrics for each window; looping per individual    
        cat("Processing ID:", i, "\n")
        indiv_data <- trk %>% filter(id == i) %>% arrange(t_) 
        
        if (nrow(indiv_data) == 0) {
          cat("No data for ID", i, "\n")
          next
        }
        
        start_date <- as.Date(min(indiv_data$t_))
        end_date <- as.Date(max(indiv_data$t_))
        windows <- create_3day_windows(start_date, end_date)
        
        for (window in windows) {    # looping by non-overlapping window
          processed_count <- processed_count + 1
          window_start <- window$start
          window_end <- window$end
          center_date <- window$center_date
          
          cat("Processing window", processed_count, "/", total_combinations, "for ID", i, "\n")
          cat("Window:", as.character(window_start), "to", as.character(window_end), "\n")
          cat("Center date:", as.character(center_date), "\n")
          
          window_data <- indiv_data %>%
            filter(t_ >= window_start, t_ <= window_end)
          
          cat("Points in window:", nrow(window_data), "\n")
          
          if (nrow(window_data) >= 5) {   # only compute metrics if >5 points within the window
            
            # keep clean track data for KDE (no mutations that break track class)
            clean_window_data <- window_data %>% arrange(t_)
            
            # compute total distance traveled on separate data (this breaks track class but we don't need it for KDE)
            distance_data <- clean_window_data %>%
              mutate(
                lead_x = lead(x_),
                lead_y = lead(y_),
                step_length = sqrt((lead_x - x_)^2 + (lead_y - y_)^2)
              )
            
            total_distance_m <- sum(distance_data$step_length, na.rm = TRUE)
            total_distance_km <- total_distance_m / 1000
            
            # compute proportions of points per state (use clean data)
            state_counts <- table(clean_window_data$state)
            state_props <- state_counts / sum(state_counts)
            
            # compute KDE areas using clean track data that maintains track class
            kde_results <- compute_kde_areas(clean_window_data, i, center_date, window_start, window_end)
            
            if (kde_results$success) {
              # combine all results
              combined_results <- list(
                ID = i,
                center_date = center_date,
                window_start = window_start,
                window_end = window_end,
                n_points = nrow(clean_window_data),
                window_days = 3,
                total_distance_km = total_distance_km,
                
                # overall UD areas
                area_50_overall_km2 = kde_results$area_50_overall_km2,
                area_95_overall_km2 = kde_results$area_95_overall_km2,
                
                # state-specific UD areas for all 3 states
                area_50_state_1_km2 = kde_results$area_50_state_1_km2,
                area_95_state_1_km2 = kde_results$area_95_state_1_km2,
                n_points_state_1 = kde_results$n_points_state_1,
                
                area_50_state_2_km2 = kde_results$area_50_state_2_km2,
                area_95_state_2_km2 = kde_results$area_95_state_2_km2,
                n_points_state_2 = kde_results$n_points_state_2,
                
                area_50_state_3_km2 = kde_results$area_50_state_3_km2,
                area_95_state_3_km2 = kde_results$area_95_state_3_km2,
                n_points_state_3 = kde_results$n_points_state_3,
                
                # state proportions
                prop_state_1 = if("1" %in% names(state_props)) as.numeric(state_props["1"]) else 0,
                prop_state_2 = if("2" %in% names(state_props)) as.numeric(state_props["2"]) else 0,
                prop_state_3 = if("3" %in% names(state_props)) as.numeric(state_props["3"]) else 0
              )
              
              MW_results <- append(MW_results, list(combined_results))
              cat("Successfully computed all UDs for ID", i, "window ending", as.character(center_date), "\n")
            }
            
          } else {
            cat("Insufficient data points (", nrow(window_data), ") for ID", i, "window ending", as.character(center_date), "\n")
          }
          cat("---\n")
        }
        cat("Completed ID:", i, "\n")
        cat("========================================\n")
      }
      cat("Processing complete! Generated", length(MW_results), "results.\n")
      
      # convert results to dataframe
      if (length(MW_results) > 0) {
        cat("Converting results to dataframe...\n")
        
        tryCatch({
          # inspect the structure of results to identify issues
          cat("Checking result structures...\n")
          lengths <- sapply(MW_results, length)
          cat("Result lengths:", unique(lengths), "\n")
          
          # check if any results failed
          success_flags <- sapply(MW_results, function(x) {
            if("success" %in% names(x)) return(x$success) else return(TRUE)
          })
          cat("Failed results:", sum(!success_flags), "out of", length(MW_results), "\n")
          
          # filter out failed results
          successful_results <- MW_results[success_flags]
          cat("Processing", length(successful_results), "successful results.\n")
          
          if(length(successful_results) == 0) {
            cat("No successful results to process.\n")
            MW_df <- NULL
          } else {
            # define all expected columns with default values
            expected_cols <- list(
              ID = NA_character_,
              center_date = as.Date(NA),
              window_start = as.POSIXct(NA),
              window_end = as.POSIXct(NA),
              n_points = NA_real_,
              window_days = 3,
              total_distance_km = NA_real_,
              area_50_overall_km2 = NA_real_,
              area_95_overall_km2 = NA_real_,
              area_50_state_1_km2 = NA_real_,
              area_95_state_1_km2 = NA_real_,
              n_points_state_1 = NA_real_,
              area_50_state_2_km2 = NA_real_,
              area_95_state_2_km2 = NA_real_,
              n_points_state_2 = NA_real_,
              area_50_state_3_km2 = NA_real_,
              area_95_state_3_km2 = NA_real_,
              n_points_state_3 = NA_real_,
              prop_state_1 = NA_real_,
              prop_state_2 = NA_real_,
              prop_state_3 = NA_real_
            )
            
            # standardize each result to have all expected columns
            standardized_results <- lapply(successful_results, function(result) {
              # start with expected structure
              std_result <- expected_cols
              
              # fill in actual values where they exist
              for(col_name in names(expected_cols)) {
                if(col_name %in% names(result) && !is.null(result[[col_name]])) {
                  std_result[[col_name]] <- result[[col_name]]
                }
              }
              
              return(std_result)
            })
            
            # now convert to dataframe
            MW_df <- do.call(rbind.data.frame, standardized_results)
            
            # ensure proper data types
            MW_df <- MW_df %>%
              mutate(
                ID = as.character(ID),
                center_date = as.Date(center_date),
                window_start = as.POSIXct(window_start),
                window_end = as.POSIXct(window_end),
                n_points = as.numeric(n_points),
                window_days = as.numeric(window_days),
                total_distance_km = as.numeric(total_distance_km),
                
                # overall areas
                area_50_overall_km2 = as.numeric(area_50_overall_km2),
                area_95_overall_km2 = as.numeric(area_95_overall_km2),
                
                # state-specific areas for all 3 states
                area_50_state_1_km2 = as.numeric(area_50_state_1_km2),
                area_95_state_1_km2 = as.numeric(area_95_state_1_km2),
                n_points_state_1 = as.numeric(n_points_state_1),
                
                area_50_state_2_km2 = as.numeric(area_50_state_2_km2),
                area_95_state_2_km2 = as.numeric(area_95_state_2_km2),
                n_points_state_2 = as.numeric(n_points_state_2),
                
                area_50_state_3_km2 = as.numeric(area_50_state_3_km2),
                area_95_state_3_km2 = as.numeric(area_95_state_3_km2),
                n_points_state_3 = as.numeric(n_points_state_3),
                
                # proportions
                prop_state_1 = as.numeric(prop_state_1),
                prop_state_2 = as.numeric(prop_state_2),
                prop_state_3 = as.numeric(prop_state_3)
              )
            
            cat("Successfully created dataframe with", nrow(MW_df), "rows and", ncol(MW_df), "columns.\n")
            print(head(MW_df, 3))
            
            # print summary statistics for all states (only for non-NA values)
            cat("\nSummary of UD areas:\n")
            cat("Overall 50% UD: mean =", round(mean(MW_df$area_50_overall_km2, na.rm = TRUE), 3), "km²\n")
            cat("Overall 95% UD: mean =", round(mean(MW_df$area_95_overall_km2, na.rm = TRUE), 3), "km²\n")
            
            cat("State 1 (Resting) 50% UD: mean =", round(mean(MW_df$area_50_state_1_km2, na.rm = TRUE), 3), "km²\n")
            cat("State 1 (Resting) 95% UD: mean =", round(mean(MW_df$area_95_state_1_km2, na.rm = TRUE), 3), "km²\n")
            
            cat("State 2 (Foraging) 50% UD: mean =", round(mean(MW_df$area_50_state_2_km2, na.rm = TRUE), 3), "km²\n")
            cat("State 2 (Foraging) 95% UD: mean =", round(mean(MW_df$area_95_state_2_km2, na.rm = TRUE), 3), "km²\n")
            
            cat("State 3 (Traveling) 50% UD: mean =", round(mean(MW_df$area_50_state_3_km2, na.rm = TRUE), 3), "km²\n")
            cat("State 3 (Traveling) 95% UD: mean =", round(mean(MW_df$area_95_state_3_km2, na.rm = TRUE), 3), "km²\n")
          }
          
        }, error = function(e) {
          cat("Error creating dataframe:", e$message, "\n")
          # debug information
          cat("Debugging information:\n")
          if(length(MW_results) > 0) {
            cat("First result structure:\n")
            print(str(MW_results[[1]]))
            if(length(MW_results) > 1) {
              cat("Second result structure:\n")
              print(str(MW_results[[2]]))
            }
          }
          MW_df <- NULL
        })
      } else {
        cat("No results to convert - list is empty.\n")
        MW_df <- NULL
      }
      
      # save final result as .csv
      if (!is.null(MW_df)) {
        write.csv(MW_df, "merged_3day_windows_multi_UD.csv", row.names = FALSE)
        cat("Results saved to merged_3day_windows_multi_UD.csv\n")
      } else {
        cat("No data to save.\n")
      }
      