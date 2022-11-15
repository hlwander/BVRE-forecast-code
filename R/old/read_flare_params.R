#' Read FLARE parameters
#' @param files vector; of files
#' @param type character; can be either "forecast" or "all". Load just the
#'  forecast parameters or all the parameters
#' @param summary logical; return summary statistics (mean, sd, median)
#' @importFrom ncdf4 nc_open ncvar_get ncatt_get nc_close
#' @importFrom reshape2 melt
#' @importFrom plyr ddply
#'
read_flare_params <- function(files, type = "forecast", summary = FALSE) {
  out <- lapply(files, function(f) {
    if(!file.exists(f)) {
      stop("File ", f, " does not exist.")
    }

    nc <- ncdf4::nc_open(f)
    t <- ncdf4::ncvar_get(nc,'time')
    local_tzone <- ncdf4::ncatt_get(nc, 0)$local_time_zone_of_simulation
    full_time <- as.POSIXct(t,
                            origin = '1970-01-01 00:00.00 UTC',
                            tz = "UTC")
    full_time_day <- lubridate::as_date(full_time)
    nsteps <- length(full_time_day)
    forecast <- ncdf4::ncvar_get(nc, 'forecast')
    # depths <- round(ncdf4::ncvar_get(nc, 'depth'),2)

    var_names <- names(nc$var)
    output_type <- rep(NA, length(var_names))
    target_name <- rep(NA, length(var_names))
    for(i in 1:length(var_names)){
      tmp <- ncdf4::ncatt_get(nc, varid = var_names[i],attname = "long_name")$value
      output_type[i] <- stringr::str_split(tmp, ":")[[1]][1]
      target_name[i] <- stringr::str_split(tmp, ":")[[1]][2]
    }

    par_names <- var_names[output_type  == "parameter"]

    if(length(which(forecast == 1)) > 0 & type == "forecast"){
      extract_index <- which(forecast == 1)[1]
    }else{
      extract_index <- 1:length(forecast)
    }


    par_list <- list()
    if(length(par_names) > 0){
      par_list <- lapply(par_names, function(x) {
        pars <- ncdf4::ncvar_get(nc, x)[extract_index, ]
        df <- data.frame(datetime = as.Date(full_time_day[extract_index]),
                   # ens = paste0("mem_", formatC(1:length(pars), width = 3, format = "d", flag = "0")),
                   value = pars)
        if(dim(df)[2] > 2) {
          colnames(df)[-1] <- paste0("mem_", formatC(1:(ncol(df)-1), width = 3, format = "d", flag = "0"))
          df <- reshape2::melt(df, id.vars = "datetime")
          colnames(df)[2] <- "member"
        } else {
          df$member <- paste0("mem_", formatC(1:length(pars), width = 3, format = "d", flag = "0"))
        }
        return(df)
      })
      names(par_list) <- par_names

    } else {
      ncdf4::nc_close(nc)
      stop("No parameters in ", f)
    }
    ncdf4::nc_close(nc)

    mlt <- reshape2::melt(par_list, id.vars = c("datetime", "member"))
    mlt <- mlt[, c("datetime", "member", "L1", "value")]
    colnames(mlt)[3] <- "parameter"
    return(mlt)
  })
  dat <- do.call(rbind, out)

  if(summary) {
    dat <- plyr::ddply(dat, c("datetime", "parameter"), function(x) {
      data.frame(mean = mean(x$value), median = median(x$value), sd = sd(x$value))
    })
  }
  return(dat)
}
