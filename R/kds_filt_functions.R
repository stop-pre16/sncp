#' Function to do KDS-Filt using random draws for number of points
#'
#' Wrapper function that does KDS-Filt, then gets cluster assignments
#'
#' @param sim_pattern Data frame with standard long format CT info for just diseased tissue
#' @param xwin Range of x values
#' @param ywin Range of y values
#' @param bw_method_c Either numeric smoothing bandwidth, or character string for methods (bw_ppl, bw_diggle, bw_CvL)
#' @param thresh Quantile to use for filtering noise
#' @param n_sim_data Number of simulation datasets to use for constructing null distribution
#'
#' @return
#' Returns a data frame with x and y coordinates, and the feature or noise classification for each point
#' @export
#'

KDS_filt <- function(sim_pattern,
                     xwin,
                     ywin,
                     bw_method_c = 'bw.ppl',
                     thresh = .5,
                     n_sim_data = 20){
  sim_pattern_full <- sim_pattern
  # print("denoising")
  if(is.numeric(bw_method_c)){
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    sigma2 <- bw_method_c
    sim_pattern$dens_d <- as.numeric(dentools::fastDMVNorm_diag_norm_sum_drop_safe(A = as.matrix(sim_pattern[, 1:2]),
                                                                              B = as.matrix(sim_pattern[, 1:2]),
                                                                              sigma2 = sigma2))

    sim_thresh <- quantile(dentools::para_surface_dist_r(n_points = n_pts_orig,
                                                         xwin = xwin,
                                                         ywin = ywin,
                                                         sigma2 = sigma2,
                                                         n_it = n_sim_data),
                           probs = thresh)

    sim_pattern_sub <- sim_pattern %>%
      dplyr::filter(dens_d > sim_thresh) %>%
      dplyr::select(x, y)
    n_pts_sub <- nrow(sim_pattern_sub)

    while(n_pts != n_pts_sub & n_pts_sub > 0){
      # print('denoising again')
      sim_pattern <- sim_pattern_sub
      n_pts <- nrow(sim_pattern)
      xwin_sub <- range(sim_pattern$x)
      ywin_sub <- range(sim_pattern$y)
      LM_sim_sub <- (xwin_sub[2] - xwin_sub[1]) * (ywin_sub[2] - ywin_sub[1])
      # dat_ppp_sub <- ppp(sim_pattern$x, sim_pattern$y, owin(xwin, ywin))
      # dens_ppp <- density(x = dat_ppp, sigma = bw_method, edge = T, diggle = T, leaveoneout = T)
      # sigma2 <- attributes(dens_ppp)$sigma^2
      # if(bw_method_c == "bw.CvL"){
      #   sigma2 <- attributes(dens_ppp)$sigma
      # }

      sim_pattern$dens_d <- as.numeric(dentools::fastDMVNorm_diag_norm_sum_drop_safe(A = as.matrix(sim_pattern[, 1:2]),
                                                                                B = as.matrix(sim_pattern[, 1:2]),
                                                                                sigma2 = sigma2))


      sim_thresh <- quantile(dentools::para_surface_dist_r(n_points = n_pts,
                                                           xwin = xwin,
                                                           ywin = ywin,
                                                           sigma2 = sigma2,
                                                           n_it = n_sim_data),
                             probs = thresh)
      sim_pattern_sub <- sim_pattern %>%
        dplyr::filter(dens_d > sim_thresh) %>%
        dplyr::select(x, y)
      n_pts_sub <- nrow(sim_pattern_sub)
    }
    if(nrow(sim_pattern_sub) == 0){
      ret <- sim_pattern_full %>% dplyr::mutate(type = 'noise')
    }
    else{
      sim_pattern_sub <- sim_pattern_sub %>% dplyr::mutate(type = 'feature')
      ret <- sim_pattern_full %>% dplyr::left_join(sim_pattern_sub, by = c('x', 'y'))
      ret$type <- ifelse(is.na(ret$type), yes = 'noise', no = ret$type)
    }
  }
  else{
    if(bw_method_c == 'bw.ppl'){
      bw_method <- bw.ppl
    }
    else if(bw_method_c == 'bw.diggle'){
      bw_method <- bw.diggle
    }
    else{
      bw_method <- bw.CvL
    }
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    dat_ppp <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin(xwin, ywin))
    dens_ppp <- tryCatch(density(x = dat_ppp, sigma = bw_method, edge = T, diggle = T, leaveoneout = T),
                         error = function(e){return(NA)})
    if(any(is.na(dens_ppp))){
      return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
    }
    sigma2 <- attributes(dens_ppp)$sigma^2
    if(bw_method_c == "bw.CvL"){
      sigma2 <- attributes(dens_ppp)$sigma
    }
    sim_pattern$dens_d <- as.numeric(dentools::fastDMVNorm_diag_norm_sum_drop_safe(A = as.matrix(sim_pattern[, 1:2]),
                                                                              B = as.matrix(sim_pattern[, 1:2]),
                                                                              sigma2 = sigma2))

    sim_thresh <- quantile(dentools::para_surface_dist_r(n_points = n_pts_orig,
                                                         xwin = xwin,
                                                         ywin = ywin,
                                                         sigma2 = sigma2,
                                                         n_it = n_sim_data),
                           probs = thresh)

    sim_pattern_sub <- sim_pattern %>%
      dplyr::filter(dens_d > sim_thresh) %>%
      dplyr::select(x, y)
    n_pts_sub <- nrow(sim_pattern_sub)

    while(n_pts != n_pts_sub & n_pts_sub > 0){
      # print('denoising again')
      sim_pattern <- sim_pattern_sub
      n_pts <- nrow(sim_pattern)
      xwin_sub <- range(sim_pattern$x)
      ywin_sub <- range(sim_pattern$y)
      LM_sim_sub <- (xwin_sub[2] - xwin_sub[1]) * (ywin_sub[2] - ywin_sub[1])
      dat_ppp_sub <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin(xwin, ywin))
      dens_ppp <- tryCatch(density(x = dat_ppp_sub, sigma = bw_method, edge = T, diggle = T, leaveoneout = T),
                           error = function(e){return(NA)})
      if(any(is.na(dens_ppp))){
        return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
      }
      sigma2 <- attributes(dens_ppp)$sigma^2
      if(bw_method_c == "bw.CvL"){
        sigma2 <- attributes(dens_ppp)$sigma
      }

      sim_pattern$dens_d <- as.numeric(dentools::fastDMVNorm_diag_norm_sum_drop_safe(A = as.matrix(sim_pattern[, 1:2]),
                                                                                B = as.matrix(sim_pattern[, 1:2]),
                                                                                sigma2 = sigma2))


      sim_thresh <- quantile(dentools::para_surface_dist_r(n_points = n_pts,
                                                           xwin = xwin,
                                                           ywin = ywin,
                                                           sigma2 = sigma2,
                                                           n_it = n_sim_data),
                             probs = thresh)
      sim_pattern_sub <- sim_pattern %>%
        dplyr::filter(dens_d > sim_thresh) %>%
        dplyr::select(x, y)
      n_pts_sub <- nrow(sim_pattern_sub)
    }
    if(nrow(sim_pattern_sub) == 0){
      ret <- sim_pattern_full %>% dplyr::mutate(type = 'noise')
    }
    else{
      sim_pattern_sub <- sim_pattern_sub %>% dplyr::mutate(type = 'feature')
      ret <- sim_pattern_full %>% dplyr::left_join(sim_pattern_sub, by = c('x', 'y'))
      ret$type <- ifelse(is.na(ret$type), yes = 'noise', no = ret$type)
    }
  }

  return(ret)
}


#' Function to do KDS-Filt using random draws for number of points
#'
#' Wrapper function that does KDS-Filt. This version is for pixelated data
#'
#' @param sim_pattern Data frame with standard long format CT info for just diseased tissue
#' @param all_points Data frame with standard long format that has all points in a given slice
#' @param bw_method_c Either numeric smoothing bandwidth, or character string for methods (bw_ppl, bw_diggle, bw_CvL)
#' @param thresh Quantile to use for filtering noise
#' @param n_sim_data Number of simulation datasets to use for constructing null distribution
#'
#' @return
#' Returns a data frame with x and y coordinates, and the feature or noise classification for each point
#'
#' @export
#'
#'
#'
#'
KDS_filt_pix <- function(sim_pattern,
                         all_points,
                         bw_method_c = 'bw.ppl',
                         thresh = .5,
                         n_sim_data = 20){
  sim_pattern_full <- sim_pattern
  # print("denoising")
  if(is.numeric(bw_method_c)){
    LM_sim <- nrow(all_points)
    xwin <- range(all_points[, 1])
    ywin <- range(all_points[, 2])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    sigma2 <- bw_method_c
    sim_pattern$dens_d <- as.numeric(dentools::fastDMVNorm_diag_norm_sum_drop_safe(A = as.matrix(sim_pattern[, 1:2]),
                                                                                   B = as.matrix(sim_pattern[, 1:2]),
                                                                                   sigma2 = sigma2))

    sim_thresh <- quantile(dentools::para_surface_dist_r_pix(n_points = n_pts_orig,
                                                             all_points = all_points,
                                                             sigma2 = sigma2,
                                                             n_it = n_sim_data),
                           probs = thresh)

    sim_pattern_sub <- sim_pattern %>%
      dplyr::filter(dens_d > sim_thresh) %>%
      dplyr::select(x, y)
    n_pts_sub <- nrow(sim_pattern_sub)

    while(n_pts != n_pts_sub & n_pts_sub > 0){
      # print('denoising again')
      sim_pattern <- sim_pattern_sub
      n_pts <- nrow(sim_pattern)
      xwin_sub <- range(sim_pattern$x)
      ywin_sub <- range(sim_pattern$y)
      LM_sim_sub <- (xwin_sub[2] - xwin_sub[1]) * (ywin_sub[2] - ywin_sub[1])
      # dat_ppp_sub <- ppp(sim_pattern$x, sim_pattern$y, owin(xwin, ywin))
      # dens_ppp <- density(x = dat_ppp, sigma = bw_method, edge = T, diggle = T, leaveoneout = T)
      # sigma2 <- attributes(dens_ppp)$sigma^2
      # if(bw_method_c == "bw.CvL"){
      #   sigma2 <- attributes(dens_ppp)$sigma
      # }

      sim_pattern$dens_d <- as.numeric(dentools::fastDMVNorm_diag_norm_sum_drop_safe(A = as.matrix(sim_pattern[, 1:2]),
                                                                                     B = as.matrix(sim_pattern[, 1:2]),
                                                                                     sigma2 = sigma2))


      sim_thresh <- quantile(dentools::para_surface_dist_r_pix(n_points = n_pts,
                                                               all_points = all_points,
                                                               sigma2 = sigma2,
                                                               n_it = n_sim_data),
                             probs = thresh)
      sim_pattern_sub <- sim_pattern %>%
        dplyr::filter(dens_d > sim_thresh) %>%
        dplyr::select(x, y)
      n_pts_sub <- nrow(sim_pattern_sub)
    }
    if(nrow(sim_pattern_sub) == 0){
      ret <- sim_pattern_full %>% dplyr::mutate(type = 'noise')
    }
    else{
      sim_pattern_sub <- sim_pattern_sub %>% dplyr::mutate(type = 'feature')
      ret <- sim_pattern_full %>% dplyr::left_join(sim_pattern_sub, by = c("x", "y"))
      ret$type <- ifelse(is.na(ret$type), yes = 'noise', no = ret$type)
    }
  }
  else{
    if(bw_method_c == 'bw.ppl'){
      bw_method <- bw.ppl
    }
    else if(bw_method_c == 'bw.diggle'){
      bw_method <- bw.diggle
    }
    else{
      bw_method <- bw.CvL
    }
    LM_sim <- nrow(all_points)
    xwin <- range(all_points[, 1])
    ywin <- range(all_points[, 2])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    dat_ppp <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin(xwin, ywin))
    dens_ppp <- tryCatch(density(x = dat_ppp, sigma = bw_method, edge = T, diggle = T, leaveoneout = T),
                         error = function(e){return(NA)})
    if(any(is.na(dens_ppp))){
      return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
    }
    sigma2 <- attributes(dens_ppp)$sigma^2
    if(bw_method_c == "bw.CvL"){
      sigma2 <- attributes(dens_ppp)$sigma
    }
    sim_pattern$dens_d <- as.numeric(dentools::fastDMVNorm_diag_norm_sum_drop_safe(A = as.matrix(sim_pattern[, 1:2]),
                                                                              B = as.matrix(sim_pattern[, 1:2]),
                                                                              sigma2 = sigma2))

    sim_thresh <- quantile(dentools::para_surface_dist_r_pix(n_points = n_pts_orig,
                                                             xwin = xwin,
                                                             ywin = ywin,
                                                             sigma2 = sigma2,
                                                             n_it = n_sim_data),
                           probs = thresh)

    sim_pattern_sub <- sim_pattern %>%
      dplyr::filter(dens_d > sim_thresh) %>%
      dplyr::select(x, y)
    n_pts_sub <- nrow(sim_pattern_sub)

    while(n_pts != n_pts_sub & n_pts_sub > 0){
      # print('denoising again')
      sim_pattern <- sim_pattern_sub
      n_pts <- nrow(sim_pattern)
      xwin_sub <- range(sim_pattern$x)
      ywin_sub <- range(sim_pattern$y)
      LM_sim_sub <- (xwin_sub[2] - xwin_sub[1]) * (ywin_sub[2] - ywin_sub[1])
      dat_ppp_sub <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin(xwin, ywin))
      dens_ppp <- tryCatch(density(x = dat_ppp_sub, sigma = bw_method, edge = T, diggle = T, leaveoneout = T),
                           error = function(e){return(NA)})
      if(any(is.na(dens_ppp))){
        return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
      }
      sigma2 <- attributes(dens_ppp)$sigma^2
      if(bw_method_c == "bw.CvL"){
        sigma2 <- attributes(dens_ppp)$sigma
      }

      sim_pattern$dens_d <- as.numeric(dentools::fastDMVNorm_diag_norm_sum_drop_safe(A = as.matrix(sim_pattern[, 1:2]),
                                                                                B = as.matrix(sim_pattern[, 1:2]),
                                                                                sigma2 = sigma2))


      sim_thresh <- quantile(dentools::para_surface_dist_r_pix(n_points = n_pts,
                                                               all_points = all_points,
                                                               sigma2 = sigma2,
                                                               n_it = n_sim_data),
                             probs = thresh)
      sim_pattern_sub <- sim_pattern %>%
        dplyr::filter(dens_d > sim_thresh) %>%
        dplyr::select(x, y)
      n_pts_sub <- nrow(sim_pattern_sub)
    }
    if(nrow(sim_pattern_sub) == 0){
      ret <- sim_pattern_full %>% dplyr::mutate(type = 'noise')
    }
    else{
      sim_pattern_sub <- sim_pattern_sub %>% dplyr::mutate(type = 'feature')
      ret <- sim_pattern_full %>% dplyr::left_join(sim_pattern_sub, by = c("x", "y"))
      ret$type <- ifelse(is.na(ret$type), yes = 'noise', no = ret$type)
    }
  }

  return(ret)
}


