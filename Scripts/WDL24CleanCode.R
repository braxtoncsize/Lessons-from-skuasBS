##############################################################################
# SKUA CASE STUDY -- CONSIDERATIONS FOR DISEASE/MOVEMENT ECOLOGY STUDY #
# Weddell Sea Cleaning and Bayesian Multilevel Regression Analyses #
# Last updated: 22Mar26 BS #
# R version 4.5.2 (2025-10-31) #

##############################################################################

# Packages
library(readr)
library(viridis)
library(diveMove)
library(tidyr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(sf)
library(brms)
library(bayesplot)

# Load data
skua_wdl24_gps <- read_csv("Raw Data/WDL24_Raw.csv", col_types = cols(Date = col_date(format = "%d/%m/%Y"), Time = col_time(format = "%H:%M:%S"), Time_complete = col_datetime(format = "%m/%d/%Y %H:%M")))

# Store Argos data
skua_wdl24_argos <- skua_wdl24_gps[skua_wdl24_gps$Data_type == "Argos",]

# Then select GPS data frame
skua_wdl24_gps <- skua_wdl24_gps[skua_wdl24_gps$Data_type == "GPS",]


## Clean up ARGOS data ---------------

    # Replace GPS reference with ARGOS
    skua_wdl24_argos$Latitude <- as.vector(skua_wdl24_argos$Latitude)
    skua_wdl24_argos$Longitude <- as.vector(skua_wdl24_argos$Longitude)
    
    # Remove invalid latitude/longitude
    skua_wdl24_argos <- skua_wdl24_argos %>% 
      filter(between(Latitude, -90, 90), between(Longitude, -180, 180))
    
    # Remove duplicates/near-duplicates within 10 seconds
    skua_wdl24_argos$Time_complete_local <- as.POSIXct(skua_wdl24_argos$Time_complete_local, format = "%m/%d/%Y %H:%M")
    
    nodups <- skua_wdl24_argos %>%
      group_by(ID) %>%
      do(distinct(., Time_complete_local, .keep_all = TRUE)) %>%
      do(mutate(
        .,
        dup = difftime(Time_complete_local, lag(Time_complete_local), units = "secs") < 10
      )) %>%
      do(arrange(., order(Time_complete_local))) %>%
      dplyr::filter(!dup) %>%
      dplyr::select(-dup) %>%
      arrange(Time_complete_local)
    
    dup.new <- nrow(skua_wdl24_argos) - nrow(nodups)
    cat(sprintf("%d duplicate time &/or near duplicate location/time records removed\n", dup.new)) 
    
    # Filter out records exceeding maximum travel speed
    vmax <- 30  # m/s
    
    speedfiltered <- nodups %>%
      do(mutate(., keep = grpSpeedFilter(cbind(Time_complete_local, Longitude, Latitude), speed.thr = vmax)))
    
    d5.rep.new <- nrow(nodups) - sum(speedfiltered$keep)
    cat(sprintf("\n%d records with travel rates > %d m/s will be ignored by SSM filter\n", d5.rep.new, vmax))
    
    # Remove flagged records
    skua_wdl24_argos <- speedfiltered %>%
      filter(keep == "TRUE") %>%
      mutate(lc = "G") %>%
      dplyr::select(Individual, ID, Date, Time, Latitude, Longitude, Fix, CRC, Time_complete, Time_complete_local, 
             Data_type, Species, Metal, Capture_ID, Island, Time_complete_local_deployment, 
             Latitude_deployment, Longitude_deployment, Comment) %>%
      ungroup()
    
    summary(skua_wdl24_argos)
    
    # Count points per individual
    skua_wdl24_argos %>%
      dplyr::group_by(ID) %>%
      dplyr::summarise(count = n())
    
    # Build ARGOS heatmap summarizing fixes
    ARGOSdata_recording <- as.data.frame(table(skua_wdl24_argos$ID, skua_wdl24_argos$Date))
    colnames(ARGOSdata_recording) <- c("ID", "date", "points")
    ARGOSdata_recording$date <- as.Date(ARGOSdata_recording$date)
    
    ARGOSfull_data <- expand_grid(
      date = seq(min(ARGOSdata_recording$date), max(ARGOSdata_recording$date), by = "day"),
      ID = unique(ARGOSdata_recording$ID)
    ) %>%
      left_join(ARGOSdata_recording, by = c("date", "ID")) %>%
      mutate(points = replace_na(points, 0))  # Fill NAs with 0
    
    heatmapARGOS <- ggplot(data = ARGOSfull_data, aes(x = date, y = ID, fill = points)) +
      geom_tile() +
      scale_fill_viridis(limits = c(0, max(ARGOSdata_recording$points)), option = "D") +
      scale_x_date(date_breaks = "5 days", date_labels = "%b %d") +
      xlab("Date") + ylab("ID") +
      labs(fill = "Argos fixes per day") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
## Clean up GPS data ---------------
    
    # Remove any long > 180 or <-180, any lat > 90 or <-90
    skua_wdl24_gps$Latitude <- as.vector(skua_wdl24_gps$Latitude)
    skua_wdl24_gps$Longitude <- as.vector(skua_wdl24_gps$Longitude)
    
    skua_wdl24_gps <- skua_wdl24_gps %>% filter(Latitude < 90)
    skua_wdl24_gps <- skua_wdl24_gps %>% filter(Latitude > -90)
    
    skua_wdl24_gps <- skua_wdl24_gps %>% filter(Longitude < 180)
    skua_wdl24_gps <- skua_wdl24_gps %>% filter(Longitude > -180)
    
    
    # Remove any duplicates or near duplicates that occurred within 10 s
    
    skua_wdl24_gps$Time_complete_local <- as.POSIXct(skua_wdl24_gps$Time_complete_local, format = "%m/%d/%Y %H:%M")
    
    nodups <- skua_wdl24_gps %>%
      group_by(ID) %>%
      do(distinct(., Time_complete_local, .keep_all = TRUE)) %>%
      do(mutate(
        .,
        dup = difftime(Time_complete_local, lag(Time_complete_local), units = "secs") < 10)) %>%
      do(arrange(., order(Time_complete_local))) %>%
      dplyr:: filter(.,!dup) %>%
      dplyr:: select(-dup) %>%
      arrange(Time_complete_local)
    
    dup.new <- nrow(skua_wdl24_gps) - nrow(nodups)
    cat(sprintf("%d duplicate time &/or near duplicate location/time records removed\n", dup.new)) 
    
    # Remove fixes with extreme travel rates
    
    vmax <- 30  ## meters/s
    
    speedfiltered <- nodups %>%
      do(mutate(., keep = grpSpeedFilter(cbind(Time_complete_local, Longitude, Latitude), speed.thr = vmax)))
    
    d5.rep.new <- nrow(nodups) - sum(speedfiltered$keep)
    cat(sprintf(
      paste(
        "\n%d records with travel rates >",
        vmax,
        "m/s will be ignored by SSM filter\n"
      ),
      d5.rep.new
    )) 
    
    # Remove all locations flagged to discard during the above prefiltering stage before running the SSM
    skua_wdl24_gps <- speedfiltered %>%
      filter(keep == "TRUE") %>%
      mutate(lc = "G") %>%
      dplyr :: select(Individual, ID, Date, Time, Latitude, Longitude, Fix, CRC, Time_complete, Time_complete_local, Data_type, Species, Metal, Capture_ID, Island, Time_complete_local_deployment, Latitude_deployment, Longitude_deployment, Comment) %>%
      ungroup()
    
    # Last screening, remove unrealistic points that remain...
    skua_wdl24_gps <- skua_wdl24_gps[skua_wdl24_gps$Latitude <=0, ]
    skua_wdl24_gps <- skua_wdl24_gps[skua_wdl24_gps$Longitude <=0, ]
    skua_wdl24_gps <- skua_wdl24_gps[skua_wdl24_gps$Longitude >= -80, ]
    
    # Count points per individual
    skua_wdl24_gps %>%
      dplyr::group_by(ID) %>%
      dplyr::summarise(count = n())
    
    # Build GPS heatmap summarizing fixes
    data_recording <- as.data.frame(table(skua_wdl24_gps$ID,skua_wdl24_gps$Date))
    colnames(data_recording)<-c("ID","date","points")
    data_recording$date <- as.Date(data_recording$date)
    
    full_data <- expand_grid(
      date = seq(min(data_recording$date), max(data_recording$date), by = "day"),
      ID = unique(data_recording$ID)
    ) %>%
      left_join(data_recording, by = c("date", "ID")) %>%
      mutate(points = replace_na(points, 0))  # Replace missing values with 0
    
    heatmapGPS <- ggplot(data = full_data, aes(x = date, y = ID, fill = points)) +
      geom_tile()+
      scale_fill_viridis(limits = c(0, max(data_recording$points)),option="D")+
      scale_x_date(
        date_breaks = "5 days",
        date_labels = "%b %d")+
      xlab("Date") +
      ylab("ID") +
      labs(fill = "GPS fixes per day") + 
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
      
## Compute composite homeranges -----------------------
      
      # Convert data to sf object with lon/lat, project as needed
      gps_sf <- st_as_sf(skua_wdl24_gps, coords = c("Longitude", "Latitude"), crs = 4326)
      arg_sf <- st_as_sf(skua_wdl24_argos, coords = c("Longitude", "Latitude"), crs = 4326)
      
      gps_proj_sf <- st_transform(gps_sf, 3031)
      arg_proj_sf <- st_transform(arg_sf, 3031)
      
      gps_proj_df <- gps_proj_sf %>%
        mutate(
          x = st_coordinates(.)[,1],
          y = st_coordinates(.)[,2]
        ) %>%
        st_drop_geometry()
      
      arg_proj_df <- arg_proj_sf %>%
        mutate(
          x = st_coordinates(.)[,1],
          y = st_coordinates(.)[,2]
        ) %>%
        st_drop_geometry()
      
      gps_proj_df <- gps_proj_df %>% filter(Individual != "Tyler")
      arg_proj_df <- arg_proj_df %>% filter(Individual != "Tyler")
      
      # Make GPS, Argos, and combined track... use this last one for creating shared trast
      skua_track_GPS <- make_track(gps_proj_df, x, y, Time_complete_local, crs = 3031)
      
      skua_track_arg <- make_track(arg_proj_df, x, y, Time_complete_local, crs = 3031)
      
      all_tracks <- bind_rows(skua_track_arg, skua_track_GPS)
      
      # Make shared trast
      shared_trast <- make_trast(all_tracks, res = 10000)
      
      # Make KDEs
      
      GPS_kde <- hr_kde(skua_track_GPS, trast = shared_trast, levels = c(0.5, 0.95))
      
      arg_kde <- hr_kde(skua_track_arg, trast = shared_trast, levels = c(0.5, 0.95))
      
      # Compare KDEs from Argos data and GPS data
      
      kdeOVERLAP <- hr_overlap(GPS_kde, arg_kde, type = "vi", conditional = FALSE)
      
      # Extract isopleths and reproject
      GPS_iso <- hr_isopleths(GPS_kde)
      arg_iso <- hr_isopleths(arg_kde)
      
      GPS_iso <- st_transform(GPS_iso, 4326)
      arg_iso <- st_transform(arg_iso, 4326)
      
      bbox <- st_bbox(GPS_iso)
      buffer_x <- (bbox["xmax"] - bbox["xmin"]) * 0.2
      buffer_y <- (bbox["ymax"] - bbox["ymin"]) * 0.2
      
      # Load and reproject coastline
      land <- ne_countries(scale = "medium", returnclass = "sf")
      land_antarctica <- st_transform(land, 4326)
      
      GPS_kde_plot <- ggplot() +
        geom_sf(data = land_antarctica, fill = "gray85", color = "gray60") +
        geom_sf(data = GPS_iso, aes(fill = as.factor(level)), color = "black", alpha = 0.5) +
        scale_fill_manual(values = c("0.5" = "#332288", "0.95" = "#88CCEE"), name = "Isopleth") +
        coord_sf(
          crs = st_crs(4326),
          xlim = c(bbox["xmin"] - buffer_x, bbox["xmax"] + buffer_x),
          ylim = c(bbox["ymin"] - buffer_y, bbox["ymax"] + buffer_y),
          expand = FALSE
        ) +
        theme_bw()
      
      ARG_kde_plot <- ggplot() +
        geom_sf(data = land_antarctica, fill = "gray85", color = "gray60") +
        geom_sf(data = arg_iso, aes(fill = as.factor(level)), color = "black", alpha = 0.5) +
        scale_fill_manual(values = c("0.5" = "#332288", "0.95" = "#88CCEE"), name = "Isopleth") +
        coord_sf(
          crs = st_crs(4326),
          xlim = c(bbox["xmin"] - buffer_x, bbox["xmax"] + buffer_x),
          ylim = c(bbox["ymin"] - buffer_y, bbox["ymax"] + buffer_y),
          expand = FALSE
        ) +
        theme_bw() 
      
      KDE_plot <- ggarrange(
        ARG_kde_plot, 
        GPS_kde_plot,
        labels = c("A", "B"),
        ncol = 2,
        common.legend = TRUE,
        legend = "right")
      
      ggsave(filename = "CompositeWDLkde.tiff", plot = KDE_plot, width = 8, height = 6, dpi = 300)
      
      # Calculate kde areas
      KDEareaGPS <- hr_area(GPS_kde)
      KDEareaARG <- hr_area(arg_kde)
      
      KDEareaGPS <- KDEareaGPS %>% select(-level, -what)
      KDEareaARG <- KDEareaARG %>% select(-level, -what)
      allareas <- cbind(KDEareaGPS, KDEareaARG)
      
# Compute individual homeranges --------------------
      
      # Prep necessary for MW analysis
      ##convert data to sf object with lon/lat
      gps_sf <- st_as_sf(skua_wdl24_gps, coords = c("Longitude", "Latitude"), crs = 4326)
      
      arg_sf <- st_as_sf(skua_wdl24_argos, coords = c("Longitude", "Latitude"), crs = 4326)
      
      gps_proj_sf <- st_transform(gps_sf, 3031)
      arg_proj_sf <- st_transform(arg_sf, 3031)
      
      gps_proj_df <- gps_proj_sf %>%
        mutate(
          X = st_coordinates(.)[,1],
          Y = st_coordinates(.)[,2]
        ) %>%
        st_drop_geometry()
      
      arg_proj_df <- arg_proj_sf %>%
        mutate(
          X = st_coordinates(.)[,1],
          Y = st_coordinates(.)[,2]
        ) %>%
        st_drop_geometry()
      
      gps_proj_df <- gps_proj_df %>% filter(Individual != "Tyler")
      arg_proj_df <- arg_proj_df %>% filter(Individual != "Tyler")

###########
      ##make tracks
      
      ##combined track... will use this for KDEs 
      skua_track_GPS <- make_track(gps_proj_df, X, Y, Time_complete_local, id = Individual, crs = 3031)
      
      skua_track_arg <- make_track(arg_proj_df, X, Y, Time_complete_local, id = Individual, crs = 3031)
      
      all_tracks <- bind_rows(skua_track_arg, skua_track_GPS)
      
      ##make shared trast
      shared_trast <- make_trast(all_tracks, res = 10000)
      
      ##split data by individual
      gps_individual_tracks <- skua_track_GPS %>%
        group_by(id) %>%
        group_split()
      
      arg_individual_tracks <- skua_track_arg %>%
        group_by(id) %>%
        group_split()
      
      
      ##make KDEs
      
      GPS_kdes <- map(gps_individual_tracks, ~ hr_kde(.x, trast = shared_trast, levels = (0.5)))
      
      arg_kdes <- map(arg_individual_tracks, ~ hr_kde(.x, trast = shared_trast, levels = c(0.5, 0.95)))
      
      ##compare KDEs from Argos data and GPS data
      kde1 <- hr_area(GPS_kdes[[1]], levels = 0.95)
      kde2 <- hr_area(GPS_kdes[[2]], levels = 0.95)
      
      kde1 <- hr_overlap(GPS_kdes[[1]], arg_kdes[[1]], type = "vi", conditional = FALSE)
      kde2 <- hr_overlap(GPS_kdes[[2]], arg_kdes[[2]], type = "vi", conditional = FALSE)
      kde3 <- hr_overlap(GPS_kdes[[3]], arg_kdes[[3]], type = "vi", conditional = FALSE)
      kde4 <- hr_overlap(GPS_kdes[[4]], arg_kdes[[4]], type = "vi", conditional = FALSE)
      kde5 <- hr_overlap(GPS_kdes[[5]], arg_kdes[[5]], type = "vi", conditional = FALSE)
      kde6 <- hr_overlap(GPS_kdes[[6]], arg_kdes[[6]], type = "vi", conditional = FALSE)
      kde7 <- hr_overlap(GPS_kdes[[7]], arg_kdes[[7]], type = "vi", conditional = FALSE)
      
      KDEoverlaps <- rbind(kde1, kde2, kde3, kde4, kde5, kde6, kde7)
    
######################################################################
      ######################################################################
      #ARGOS MOVING WINDOW - NON-OVERLAPPING 3-DAY WINDOWS -----------------
      # Load required libraries
      library(dplyr)
      library(amt)           # For mk_track, hr_kde, hr_isopleths
      library(sf)            # For spatial data handling
      
      # Ensure 'time' is POSIXct BEFORE creating the track
      data_hmmNA_joinedsf = arg_proj_df
      colnames(data_hmmNA_joinedsf)[8] <- "time"
      
      cat("Original time class:", class(data_hmmNA_joinedsf$time), "\n")
      cat("Sample time values:", head(data_hmmNA_joinedsf$time), "\n")
      
      if (is.numeric(data_hmmNA_joinedsf$time)) {
        cat("Converting from UNIX timestamp to POSIXct...\n")
        data_hmmNA_joinedsf$time <- as.POSIXct(data_hmmNA_joinedsf$time, origin = "1970-01-01", tz = "UTC")
      } else {
        cat("Converting existing time to POSIXct...\n")
        data_hmmNA_joinedsf$time <- as.POSIXct(data_hmmNA_joinedsf$time, tz = "UTC")
      }
      
      str(data_hmmNA_joinedsf$time)
      cat("Converted time range:", as.character(range(data_hmmNA_joinedsf$time, na.rm = TRUE)), "\n")
      
      # Create track
      trk <- mk_track(data_hmmNA_joinedsf, .x = x, .y = y, .t = time, id = ID, type = Data_type)
      str(trk$t_)
      cat("Track time range:", as.character(range(trk$t_, na.rm = TRUE)), "\n")
      
      # Function to create non-overlapping 3-day windows
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
              center_date = current_date + 1  # Middle day of the 3-day window
            )))
          }
          
          current_date <- current_date + 3  # Move by 3 days for non-overlapping windows
        }
        
        return(windows)
      }
      
      # Non-overlapping window computation
      MW_kderesults <- list()
      processed_count <- 0
      total_combinations <- 0
      
      # Calculate total combinations
      for (i in unique(trk$id)) {
        indiv_data <- trk %>% filter(id == i) %>% as_tibble() %>% arrange(t_)
        if (nrow(indiv_data) > 0) {
          start_date <- as.Date(min(indiv_data$t_))
          end_date <- as.Date(max(indiv_data$t_))
          windows <- create_3day_windows(start_date, end_date)
          total_combinations <- total_combinations + length(windows)
        }
      }
      cat("Total ID-window combinations to process:", total_combinations, "\n")
      
      for (i in unique(trk$id)) {
        cat("Processing ID:", i, "\n")
        indiv_data <- trk %>% filter(id == i) %>% as_tibble() %>% arrange(t_)
        
        if (nrow(indiv_data) == 0) {
          cat("No data for ID", i, "\n")
          next
        }
        
        start_date <- as.Date(min(indiv_data$t_))
        end_date <- as.Date(max(indiv_data$t_))
        windows <- create_3day_windows(start_date, end_date)
        
        for (window in windows) {
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
          
          if (nrow(window_data) >= 5) {
            # Compute total distance traveled (Euclidean distance)
            window_data <- window_data %>%
              arrange(t_) %>%
              mutate(
                lead_x = lead(x_),
                lead_y = lead(y_),
                step_length = sqrt((lead_x - x_)^2 + (lead_y - y_)^2)
              )
            
            total_distance_m <- sum(window_data$step_length, na.rm = TRUE)
            total_distance_km <- total_distance_m / 1000
            
            tryCatch({
              window_track <- mk_track(window_data, .x = x_, .y = y_, .t = t_, id = id)
              kde <- window_track %>% hr_kde()
              
              # Calculate both 50% and 95% isopleths
              hr_50 <- kde %>% hr_isopleths(levels = 0.50)
              hr_95 <- kde %>% hr_isopleths(levels = 0.95)
              
              # Combine both levels into one sf object
              hr_combined <- rbind(hr_50, hr_95)
              
              # Add metadata to combined object
              hr_combined$ID <- i
              hr_combined$center_date <- center_date
              hr_combined$window_start <- window_start
              hr_combined$window_end <- window_end
              hr_combined$n_points <- nrow(window_data)
              hr_combined$window_days <- 3  # Always 3 days
              hr_combined$total_distance_km <- total_distance_km
              
              MW_kderesults <- append(MW_kderesults, list(hr_combined))
              cat("Successfully computed KDE for ID", i, "window ending", as.character(center_date), "\n")
              
            }, error = function(e) {
              cat("Error computing KDE for ID", i, "window ending", as.character(center_date), ":", e$message, "\n")
            })
          } else {
            cat("Insufficient data points (", nrow(window_data), ") for ID", i, "window ending", as.character(center_date), "\n")
          }
          cat("---\n")
        }
        cat("Completed ID:", i, "\n")
        cat("========================================\n")
      }
      cat("Processing complete! Generated", length(MW_kderesults), "KDE results.\n")
      
      # Convert results to dataframe
      if (length(MW_kderesults) > 0) {
        cat("Converting results to dataframe...\n")
        
        tryCatch({
          MW_kderesults_df <- do.call(rbind, MW_kderesults)
          
          cat("Successfully created sf dataframe with", nrow(MW_kderesults_df), "home range polygons.\n")
          
          MW_summary_df <- MW_kderesults_df %>%
            dplyr::mutate(
              center_date = as.Date(center_date),
              window_start = as.POSIXct(window_start),
              window_end = as.POSIXct(window_end),
              area_km2 = as.numeric(st_area(MW_kderesults_df)) / 1000000,
              total_distance_km = as.numeric(total_distance_km)
            ) %>%
            st_drop_geometry() %>%
            dplyr::select(ID, center_date, window_start, window_end, n_points, window_days,
                          area_km2, total_distance_km, level)
          
          cat("Created summary dataframe with", nrow(MW_summary_df), "rows.\n")
          print(head(MW_summary_df, 6))  # Show more rows to see both levels
          
        }, error = function(e) {
          cat("Error creating dataframe:", e$message, "\n")
          cat("Attempting manual extraction...\n")
          
          MW_basic_df <- data.frame(
            ID = sapply(MW_kderesults, function(x) x$ID),
            center_date = sapply(MW_kderesults, function(x) as.character(x$center_date)),
            window_start = sapply(MW_kderesults, function(x) as.character(x$window_start)),
            window_end = sapply(MW_kderesults, function(x) as.character(x$window_end)),
            n_points = sapply(MW_kderesults, function(x) x$n_points),
            window_days = sapply(MW_kderesults, function(x) x$window_days),
            total_distance_km = sapply(MW_kderesults, function(x) x$total_distance_km),
            level = sapply(MW_kderesults, function(x) x$level)
          )
          
          MW_basic_df$window_start <- as.POSIXct(MW_basic_df$window_start)
          MW_basic_df$window_end <- as.POSIXct(MW_basic_df$window_end)
          
          cat("Created basic dataframe with", nrow(MW_basic_df), "rows.\n")
          print(head(MW_basic_df, 6))
        })
      } else {
        cat("No results to convert - list is empty.\n")
      }
      
      MW_summary_df$center_date <- as.Date(MW_summary_df$center_date)
      
      # Join with metadata
      MW_summary_df$type <- "Argos"
      ARGOSMW <- MW_summary_df
      
      # GPS Moving window - NON-OVERLAPPING 3-DAY WINDOWS ------------------
      # Ensure 'time' is POSIXct BEFORE creating the track
      data_hmmNA_joinedsf = gps_proj_df
      colnames(data_hmmNA_joinedsf)[8] <- "time"
      
      cat("Original time class:", class(data_hmmNA_joinedsf$time), "\n")
      cat("Sample time values:", head(data_hmmNA_joinedsf$time), "\n")
      
      if (is.numeric(data_hmmNA_joinedsf$time)) {
        cat("Converting from UNIX timestamp to POSIXct...\n")
        data_hmmNA_joinedsf$time <- as.POSIXct(data_hmmNA_joinedsf$time, origin = "1970-01-01", tz = "UTC")
      } else {
        cat("Converting existing time to POSIXct...\n")
        data_hmmNA_joinedsf$time <- as.POSIXct(data_hmmNA_joinedsf$time, tz = "UTC")
      }
      
      str(data_hmmNA_joinedsf$time)
      cat("Converted time range:", as.character(range(data_hmmNA_joinedsf$time, na.rm = TRUE)), "\n")
      
      # Create track
      trk <- mk_track(data_hmmNA_joinedsf, .x = x, .y = y, .t = time, id = ID, type = Data_type)
      str(trk$t_)
      cat("Track time range:", as.character(range(trk$t_, na.rm = TRUE)), "\n")
      
      # Non-overlapping window computation for GPS
      MW_kderesults <- list()
      processed_count <- 0
      total_combinations <- 0
      
      # Calculate total combinations
      for (i in unique(trk$id)) {
        indiv_data <- trk %>% filter(id == i) %>% as_tibble() %>% arrange(t_)
        if (nrow(indiv_data) > 0) {
          start_date <- as.Date(min(indiv_data$t_))
          end_date <- as.Date(max(indiv_data$t_))
          windows <- create_3day_windows(start_date, end_date)
          total_combinations <- total_combinations + length(windows)
        }
      }
      cat("Total ID-window combinations to process:", total_combinations, "\n")
      
      for (i in unique(trk$id)) {
        cat("Processing ID:", i, "\n")
        indiv_data <- trk %>% filter(id == i) %>% as_tibble() %>% arrange(t_)
        
        if (nrow(indiv_data) == 0) {
          cat("No data for ID", i, "\n")
          next
        }
        
        start_date <- as.Date(min(indiv_data$t_))
        end_date <- as.Date(max(indiv_data$t_))
        windows <- create_3day_windows(start_date, end_date)
        
        for (window in windows) {
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
          
          if (nrow(window_data) >= 5) {
            # Compute total distance traveled (Euclidean distance)
            window_data <- window_data %>%
              arrange(t_) %>%
              mutate(
                lead_x = lead(x_),
                lead_y = lead(y_),
                step_length = sqrt((lead_x - x_)^2 + (lead_y - y_)^2)
              )
            
            total_distance_m <- sum(window_data$step_length, na.rm = TRUE)
            total_distance_km <- total_distance_m / 1000
            
            tryCatch({
              window_track <- mk_track(window_data, .x = x_, .y = y_, .t = t_, id = id)
              kde <- window_track %>% hr_kde()
              
              # Calculate both 50% and 95% isopleths
              hr_50 <- kde %>% hr_isopleths(levels = 0.50)
              hr_95 <- kde %>% hr_isopleths(levels = 0.95)
              
              # Combine both levels into one sf object
              hr_combined <- rbind(hr_50, hr_95)
              
              # Add metadata to combined object
              hr_combined$ID <- i
              hr_combined$center_date <- center_date
              hr_combined$window_start <- window_start
              hr_combined$window_end <- window_end
              hr_combined$n_points <- nrow(window_data)
              hr_combined$window_days <- 3  # Always 3 days
              hr_combined$total_distance_km <- total_distance_km
              
              MW_kderesults <- append(MW_kderesults, list(hr_combined))
              cat("Successfully computed KDE for ID", i, "window ending", as.character(center_date), "\n")
              
            }, error = function(e) {
              cat("Error computing KDE for ID", i, "window ending", as.character(center_date), ":", e$message, "\n")
            })
          } else {
            cat("Insufficient data points (", nrow(window_data), ") for ID", i, "window ending", as.character(center_date), "\n")
          }
          cat("---\n")
        }
        cat("Completed ID:", i, "\n")
        cat("========================================\n")
      }
      cat("Processing complete! Generated", length(MW_kderesults), "KDE results.\n")
      
      # Convert results to dataframe
      if (length(MW_kderesults) > 0) {
        cat("Converting results to dataframe...\n")
        
        tryCatch({
          MW_kderesults_df <- do.call(rbind, MW_kderesults)
          
          cat("Successfully created sf dataframe with", nrow(MW_kderesults_df), "home range polygons.\n")
          
          MW_summary_df <- MW_kderesults_df %>%
            dplyr::mutate(
              center_date = as.Date(center_date),
              window_start = as.POSIXct(window_start),
              window_end = as.POSIXct(window_end),
              area_km2 = as.numeric(st_area(MW_kderesults_df)) / 1000000,
              total_distance_km = as.numeric(total_distance_km)
            ) %>%
            st_drop_geometry() %>%
            dplyr::select(ID, center_date, window_start, window_end, n_points, window_days,
                          area_km2, total_distance_km, level)
          
          cat("Created summary dataframe with", nrow(MW_summary_df), "rows.\n")
          print(head(MW_summary_df, 6))  # Show more rows to see both levels
          
        }, error = function(e) {
          cat("Error creating dataframe:", e$message, "\n")
          cat("Attempting manual extraction...\n")
          
          MW_basic_df <- data.frame(
            ID = sapply(MW_kderesults, function(x) x$ID),
            center_date = sapply(MW_kderesults, function(x) as.character(x$center_date)),
            window_start = sapply(MW_kderesults, function(x) as.character(x$window_start)),
            window_end = sapply(MW_kderesults, function(x) as.character(x$window_end)),
            n_points = sapply(MW_kderesults, function(x) x$n_points),
            window_days = sapply(MW_kderesults, function(x) x$window_days),
            total_distance_km = sapply(MW_kderesults, function(x) x$total_distance_km),
            level = sapply(MW_kderesults, function(x) x$level)
          )
          
          MW_basic_df$window_start <- as.POSIXct(MW_basic_df$window_start)
          MW_basic_df$window_end <- as.POSIXct(MW_basic_df$window_end)
          
          cat("Created basic dataframe with", nrow(MW_basic_df), "rows.\n")
          print(head(MW_basic_df, 6))
        })
      } else {
        cat("No results to convert - list is empty.\n")
      }
      
      MW_summary_df$center_date <- as.Date(MW_summary_df$center_date)
      
      # Join with metadata
      MW_summary_df$type <- "GPS"
      GPSMW <- MW_summary_df
      
      # Merge Argos and GPS moving windows
      fullwindow <- rbind(GPSMW, ARGOSMW)
      
      # Optional: Create separate dataframes for each level if needed
      fullwindow_50 <- fullwindow %>% filter(level == 0.5)
      fullwindow_95 <- fullwindow %>% filter(level == 0.95)
      
      cat("Final results:\n")
      cat("Combined data:", nrow(fullwindow), "rows\n")
      cat("50% UD data:", nrow(fullwindow_50), "rows\n")
      cat("95% UD data:", nrow(fullwindow_95), "rows\n")
      
      # Show level distribution
      cat("Level distribution:\n")
      print(table(fullwindow$level, fullwindow$type))
      
      # Strip "ID " prefix for plotting
      fullwindow_50 <- fullwindow_50 %>%
        mutate(ID_stripped = gsub("^ID\\s*", "", ID))
      
      ggplot(fullwindow_50, aes(x = ID_stripped, y = area_km2,
                                fill = type)) +
        geom_boxplot(position = position_dodge(width = 0.8)) +
        scale_y_log10(
          name = "50% UD Area (km²)",
          breaks = c(0.01, 0.1, 1, 10, 100, 1000, 10000),
          labels = c("0.01", "0.1", "1", "10", "100", "1000", "10000")
        ) +
        scale_fill_manual(values = c("#30123B", "#1AE4B6")) +
        labs(
          x = "ID",
          y = "Log 50% UD Area",
          fill = "Logger Type"
        ) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right"
        )
      
      ggsave("WDLidArea.png", dpi = 600, bg = "white", width = 5, height = 5, units = "in")

###########################################
    
##Stats! ----------------------
    
    fullwindow_50$log_area_50 <- log(fullwindow_50$area_km2)
    fullwindow_95$log_area_95 <- log(fullwindow_95$area_km2)
    fullwindow_50$log_distance <- log(fullwindow_50$total_distance_km)
    season_start <- as.Date("2024-12-20")
    
    fullwindow_50 <- fullwindow_50 %>%
      mutate(
        window_start = as.Date(window_start),
        day_of_season = as.integer(window_start - season_start) + 1
      )
    
    fullwindow_95 <- fullwindow_95 %>%
      mutate(
        window_start = as.Date(window_start),
        day_of_season = as.integer(window_start - season_start) + 1
      )
    
# area 50% UD km2
    bayes_area50WDL <- brm(bf(log_area_50 ~ type * day_of_season + (1 | ID),
                           sigma ~ type),
                           data = fullwindow_50,
                           family = gaussian(),
                           prior = c(
                             prior(normal(0, 2), class = b),
                             prior(normal(0, 1), class = b, dpar = sigma)),
                           chains = 4,iter = 2000, cores = 4, seed = 123,
                           save_pars = save_pars(all = TRUE))
    
    shinystan::launch_shinystan(bayes_area50WDL)
    
# distance km
    bayes_distanceWDL <- brm(bf(log_distance ~ type * day_of_season + (1 | ID),
                                    sigma ~ type),
                                 data = fullwindow_50,
                                 family = gaussian(),
                             prior = c(
                               prior(normal(0, 2), class = b),
                               prior(normal(0, 2), class = b, dpar = sigma)),
                            chains = 4, iter = 2000, cores = 4, seed = 123,
                            save_pars = save_pars(all = TRUE))
    
    shinystan::launch_shinystan(bayes_distanceWDL)
    
# area 95% UD km2
    bayes_area95WDL <- brm(bf(log_area_95 ~ type * day_of_season + (1 | ID),
                              sigma ~ type),
                           data = fullwindow_95,
                           family = gaussian(),
                           prior = c(
                             prior(normal(0, 2), class = b),
                             prior(normal(0, 1), class = b, dpar = sigma)),
                           chains = 4, iter = 2000, cores = 4, seed = 123,
                           save_pars = save_pars(all = TRUE))
    
    shinystan::launch_shinystan(bayes_area95WDL)
   
####### Ridge plots
    
    library(tidyverse)
    library(ggridges)
    library(viridis)
    library(patchwork)
    
    # Function to extract posterior samples - INCLUDE INTERCEPT
    extract_posterior_samples_WDL <- function(model, model_name) {
      posterior_samples <- as_draws_df(model)
      
      # Get ALL fixed effect columns (including Intercept)
      fe_cols <- posterior_samples %>%
        select(starts_with("b_")) %>%
        names()
      
      posterior_samples %>%
        select(all_of(fe_cols), .draw) %>%
        pivot_longer(cols = -.draw, names_to = "term", values_to = "value") %>%
        mutate(
          model = model_name,
          term = str_remove(term, "^b_"),
          submodel = if_else(str_detect(term, "^sigma_"), "Sigma sub-model", "Mean sub-model"),
          term_clean = case_when(
            str_detect(term, "^sigma_") ~ str_remove(term, "^sigma_"),
            TRUE ~ term
          ),
          term_clean = case_when(
            term_clean == "Intercept" ~ "Intercept",
            str_detect(term_clean, "typeGPS$") ~ "GPS",
            term_clean == "day_of_season" ~ "Day of Season", 
            str_detect(term_clean, "typeGPS:day_of_season") ~ "GPS × Day of Season",
            TRUE ~ term_clean
          ),
          term_group = case_when(
            term_clean %in% c("Day of Season", "GPS × Day of Season") ~ "Season terms",
            TRUE ~ "Main terms"
          )
        )
    }
    
    # Extract posterior samples for all models
    posterior_WDL <- suppressWarnings({
      map_dfr(names(models_WDL), function(model_name) {
        extract_posterior_samples_WDL(models_WDL[[model_name]], model_name)
      })
    })
    
    # Create ordered factors - separate orders for each group
    main_term_order <- c("GPS", "Intercept")
    season_term_order <- c("GPS × Day of Season", "Day of Season")
    
    posterior_WDL <- posterior_WDL %>%
      mutate(
        submodel_ordered = factor(submodel, levels = c("Mean sub-model", "Sigma sub-model")),
        model = factor(model, levels = c("Distance", "50% UD Area", "95% UD Area"))
      )
    
    # Split into main terms and season terms
    main_posterior <- posterior_WDL %>% 
      filter(term_group == "Main terms") %>%
      mutate(term_clean = factor(term_clean, levels = main_term_order))
    
    season_posterior <- posterior_WDL %>% 
      filter(term_group == "Season terms") %>%
      mutate(term_clean = factor(term_clean, levels = season_term_order))
    
    # Main terms ridge plot (including Intercept)
    p_main_ridges <- ggplot(main_posterior, 
                            aes(x = value, y = term_clean, 
                                fill = submodel_ordered, color = submodel_ordered)) +
      geom_density_ridges(alpha = 0.7, scale = 0.9, size = 0.5) +
      geom_vline(xintercept = 0, color = "darkgrey", linewidth = 0.8, alpha = 0.8) +
      facet_grid(rows = vars(model), scales = "free_y") +
      scale_fill_manual(values = c("Mean sub-model" = "#30123B", "Sigma sub-model" = "#1AE4B6")) +
      scale_color_manual(values = c("Mean sub-model" = "#30123B", "Sigma sub-model" = "#1AE4B6")) +
      theme_bw() +
      theme(
        strip.text = element_text(size = 11),
        strip.background = element_rect(fill = "gray95", color = "white"),
        legend.position = "top",
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.title = element_text()
      ) +
      labs(
        x = "Effect Size (log scale)",
        y = NULL,
        fill = "Sub-model",
        color = "Sub-model"
      )
    
    # Season terms ridge plot  
    p_season_ridges <- ggplot(season_posterior, 
                              aes(x = value, y = term_clean, 
                                  fill = submodel_ordered, color = submodel_ordered)) +
      geom_density_ridges(alpha = 0.7, scale = 0.9, size = 0.5) +
      geom_vline(xintercept = 0, color = "darkgrey", linewidth = 0.8, alpha = 0.8) +
      facet_grid(rows = vars(model), scales = "free_y") +
      scale_fill_manual(values = c("Mean sub-model" = "#30123B", "Sigma sub-model" = "#1AE4B6")) +
      scale_color_manual(values = c("Mean sub-model" = "#30123B", "Sigma sub-model" = "#1AE4B6")) +
      theme_bw() +
      theme(
        strip.text = element_text(size = 11),
        strip.background = element_rect(fill = "gray95", color = "white"),
        legend.position = "none",  # Remove legend from second plot
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12)
      ) +
      labs(
        x = "Effect Size (log scale)",
        y = NULL
      )
    
    # Combine plots
    combined_ridges <- p_main_ridges / p_season_ridges
    ggsave("WDL_ridge_effects.png", combined_ridges, bg = "white", dpi = 600, width = 10, height = 10, units = "in")
    
    # Individual plots if you prefer
    ggsave("WDL_main_ridges.png", p_main_ridges, bg = "white", dpi = 600, width = 6, height = 6, units = "in")
    ggsave("WDL_season_ridges.png", p_season_ridges, bg = "white", dpi = 600, width = 6, height = 6, units = "in")
    
######### Visualize effects

    library(dplyr)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(tidyr)
    library(broom.mixed)
    
    # Function to extract both mean and sigma effects
    extract_all_effects <- function(model, model_name) {
      suppressWarnings(tidy(model, effects = "fixed", conf.int = TRUE)) %>%
        filter(!str_detect(term, "Intercept")) %>%
        mutate(
          model = model_name,
          effect_type = case_when(
            str_detect(term, "sigma") ~ "Sigma Effects",
            TRUE ~ "Mean Effects"
          ),
          exp_estimate = exp(estimate),
          exp_conf.low = exp(conf.low),
          exp_conf.high = exp(conf.high),
          term_clean = case_when(
            str_detect(term, "sigma") ~ str_remove(term, "b_sigma_"),
            TRUE ~ str_remove(term, "b_")
          ),
          term_clean = str_replace_all(term_clean, c(
            "type" = "",
            ":" = " × "
          )),
          term_clean = str_trim(term_clean)
        ) %>%
        select(model, effect_type, term_clean, estimate, exp_estimate, exp_conf.low, exp_conf.high)
    }
    
    # Your three models
    models <- list(
      "Distance" = bayes_distanceWDL,
      "Overall 50% UD" = bayes_area50WDL,
      "Overall 95% UD" = bayes_area95WDL
    )
    
    # Extract all effects
    all_effects_WDL <- map_dfr(names(models), ~extract_all_effects(models[[.]], .))
    
    # Clean labels
    all_effects_WDL <- all_effects_WDL %>%
      mutate(
        term_clean = case_when(
          term_clean %in% c("sigma_GPS", "GPS") ~ "GPS",
          term_clean == "day_of_season" ~ "Day of Season",
          term_clean == "GPS × day_of_season" ~ "GPS × Day of Season",
          term_clean == "sigma_GPS × day_of_season" ~ "GPS × Day of Season",
          TRUE ~ term_clean
        ),
        term_group = case_when(
          term_clean %in% c("Day of Season", "GPS × Day of Season") ~ "Season terms",
          TRUE ~ "Main terms"
        )
      )
    
    plot_main <- all_effects_WDL %>%
      filter(term_group == "Main terms") %>%
      ggplot(aes(x = exp_estimate, y = term_clean)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7, linewidth = 0.8) +
      geom_errorbarh(aes(xmin = exp_conf.low, xmax = exp_conf.high, color = effect_type),
                     height = 0.25, alpha = 0.8, linewidth = 1) +
      geom_point(aes(color = effect_type), size = 3.5, alpha = 0.9) +
      facet_grid(model ~ ., scales = "free_y") +
      scale_colour_viridis_d(name = "Sub-model", option = "turbo") +
      scale_x_continuous(
        name = "Multiplicative Effect",
        breaks = c(0.3, 0.7, 1.0, 1.4, 2.0, 3.0),
        labels = c("0.3×", "0.7×", "1×", "1.4×", "2×", "3×")
      ) +
      labs(y = "Comparison") +
      theme_minimal() +
      theme(
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "gray95", color = "white"),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(alpha = 0.3),
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 11)
      )
    
    plot_main
    
    plot_season <- all_effects_WDL %>%
      filter(term_group == "Season terms") %>%
      ggplot(aes(x = exp_estimate, y = term_clean)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7, linewidth = 0.8) +
      geom_errorbarh(aes(xmin = exp_conf.low, xmax = exp_conf.high, color = effect_type),
                     height = 0.25, alpha = 0.8, linewidth = 1) +
      geom_point(aes(color = effect_type), size = 3.5, alpha = 0.9) +
      facet_grid(model ~ ., scales = "free_y") +
      scale_colour_viridis_d(name = "Sub-model", option = "turbo") +
      scale_x_continuous(
        name = "Multiplicative Effect per 1-day increase",
        breaks = c(0.95, 1.00, 1.05, 1.10, 1.15),
        labels = c("0.95×", "1×", "1.05×", "1.10×", "1.15×")
      ) +
      labs(y = "Comparison") +
      theme_minimal() +
      theme(
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "gray95", color = "white"),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(alpha = 0.3),
        legend.position = "none",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 11)
      )
    
    plot_season
    
    library(patchwork)
    together <- plot_main / plot_season
    ggsave("WDLeffects.png", together, bg = "white", dpi = 600, width = 10, height = 10, units = "in")
    
    ce <- conditional_effects(bayes_area50WDL)
    
    plotBAR <- plot(ce)[[1]] +
      scale_color_discrete(name = "Logger Type", labels = c("Argos", "GPS"), palette = c("#30123B", "#1AE4B6" )) +
      scale_fill_discrete(name = "Logger Type", labels = c("Argos", "GPS"), palette = c("#30123B", "#1AE4B6" )) + 
      xlab("Logger Type") +
      ylab("Log 50% UD Area") +
      theme_bw()
    
    plotLINE <- plot(ce)[[3]] +
      scale_color_discrete(name = "Logger Type", labels = c("Argos", "GPS"), palette = c("#30123B", "#1AE4B6" )) +
      scale_fill_discrete(name = "Logger Type", labels = c("Argos", "GPS"), palette = c("#30123B", "#1AE4B6" )) + 
      xlab("Day of Season") +
      ylab("Log 50% UD Area") +
      theme_bw()

    ggsave("WDLlinepred.png", dpi = 600, bg = "white", width = 5, height = 5, units = "in")
    
    
################ MAP
    
    library(sf)
    library(dplyr)
    library(ggplot2)
    library(ggspatial)
    library(patchwork)  # for combining plots
    
    # -----------------------------
    # 1. GPS Track data (high resolution)
    # -----------------------------
    data_indiv_gps <- gps_proj_df %>%
      filter(ID == "ID 283452", !is.na(x), !is.na(y)) %>%
      arrange(Time_complete)
    
    points_sf_gps <- st_as_sf(data_indiv_gps, coords = c("x", "y"), crs = 3031)
    
    track_line_gps <- points_sf_gps %>%
      summarise(do_union = FALSE) %>%
      st_cast("LINESTRING")
    
    # -----------------------------
    # 2. ARG Track data (larger scale, messier) - SAME ID
    # -----------------------------
    data_indiv_arg <- arg_proj_df %>%
      filter(ID == "ID 283452", !is.na(x), !is.na(y)) %>%  # Same ID filter as GPS
      arrange(Time_complete)  # adjust timestamp column name if different
    
    points_sf_arg <- st_as_sf(data_indiv_arg, coords = c("x", "y"), crs = 3031)
    
    track_line_arg <- points_sf_arg %>%
      summarise(do_union = FALSE) %>%
      st_cast("LINESTRING")
    
    # -----------------------------
    # 3. Calculate bounding boxes
    # -----------------------------
    # GPS bbox (smaller scale)
    bbox_gps <- st_bbox(points_sf_gps)
    xmin_gps <- as.numeric(bbox_gps["xmin"])
    xmax_gps <- as.numeric(bbox_gps["xmax"])
    ymin_gps <- as.numeric(bbox_gps["ymin"])
    ymax_gps <- as.numeric(bbox_gps["ymax"])
    
    x_center_gps <- (xmin_gps + xmax_gps) / 2
    y_center_gps <- (ymin_gps + ymax_gps) / 2
    width_gps <- xmax_gps - xmin_gps
    height_gps <- ymax_gps - ymin_gps
    half_size_gps <- max(width_gps, height_gps) / 2 * 1.2
    
    xlim_gps <- c(x_center_gps - half_size_gps, x_center_gps + half_size_gps)
    ylim_gps <- c(y_center_gps - half_size_gps, y_center_gps + half_size_gps)
    
    # Create GPS bbox polygon for showing on ARG map
    plot_bbox_gps <- st_as_sfc(
      st_bbox(c(
        xmin = xlim_gps[1],
        ymin = ylim_gps[1],
        xmax = xlim_gps[2],
        ymax = ylim_gps[2]
      ), crs = st_crs(points_sf_gps))
    )
    
    # ARG bbox (larger scale)
    bbox_arg <- st_bbox(points_sf_arg)
    xmin_arg <- as.numeric(bbox_arg["xmin"])
    xmax_arg <- as.numeric(bbox_arg["xmax"])
    ymin_arg <- as.numeric(bbox_arg["ymin"])
    ymax_arg <- as.numeric(bbox_arg["ymax"])
    
    x_center_arg <- (xmin_arg + xmax_arg) / 2
    y_center_arg <- (ymin_arg + ymax_arg) / 2
    width_arg <- xmax_arg - xmin_arg
    height_arg <- ymax_arg - ymin_arg
    half_size_arg <- max(width_arg, height_arg) / 2 * 1.2
    
    xlim_arg <- c(x_center_arg - half_size_arg, x_center_arg + half_size_arg)
    ylim_arg <- c(y_center_arg - half_size_arg, y_center_arg + half_size_arg)
    
    plot_bbox_arg <- st_as_sfc(
      st_bbox(c(
        xmin = xlim_arg[1],
        ymin = ylim_arg[1],
        xmax = xlim_arg[2],
        ymax = ylim_arg[2]
      ), crs = st_crs(points_sf_arg))
    )
    
    # -----------------------------
    # 4. Read and prepare SCAR ADD polygon
    # -----------------------------
    add_poly <- st_read("add_coastline_high_res_polygon_v7_12.shp")
    add_poly <- st_transform(add_poly, st_crs(points_sf_gps))
    add_poly <- st_make_valid(add_poly)
    
    # Crop coastline for each map extent
    add_local_gps <- st_intersection(add_poly, plot_bbox_gps)
    add_local_arg <- st_intersection(add_poly, plot_bbox_arg)
    
    # -----------------------------
    # 5. Create plots
    # -----------------------------
    # ARG plot (larger scale) with GPS bbox overlay
    plot_arg <- ggplot() +
      geom_sf(data = add_local_arg, fill = "grey80", colour = "black", linewidth = 0.2) +
      geom_sf(data = track_line_arg, colour = "#30123B", linewidth = 0.5, alpha = 0.7) +  # ARG tracks
      geom_sf(data = plot_bbox_gps, fill = NA, colour = "black", linewidth = 0.5) +  # GPS extent box
      coord_sf(xlim = xlim_arg, ylim = ylim_arg, expand = FALSE) +
      theme_classic() +
      theme(
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 10, face = "bold")
      ) +
      annotation_scale(location = "bl", width_hint = 0.3) +
      ggtitle("Argos-PTT")
    
    # GPS plot (smaller scale, detailed)
    plot_gps <- ggplot() +
      geom_sf(data = add_local_gps, fill = "grey80", colour = "black", linewidth = 0.2) +
      geom_sf(data = track_line_gps, colour = "#1AE4B6", linewidth = 0.8) +
      coord_sf(xlim = xlim_gps, ylim = ylim_gps, expand = FALSE) +
      theme_classic() +
      theme(
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 10, face = "bold")
      ) +
      annotation_scale(location = "bl", width_hint = 0.3) +
      ggtitle("GPS")
    
    # -----------------------------
    # 6. Combine plots
    # -----------------------------
    combined_plot <- plot_arg + plot_gps
    
    # Save the combined figure
    ggsave("combined_tracking_maps_ID283452.png", combined_plot, 
           dpi = 600, bg = "white", width = 8, height = 5)
    
    # Optional: Save individual plots too
    ggsave("ARG_map_ID283452.png", plot_arg, dpi = 600, bg = "white", width = 5, height = 5)
    ggsave("GPS_map_ID283452.png", plot_gps, dpi = 600, bg = "white", width = 5, height = 5)
    
    