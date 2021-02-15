rk_bound = function(x = c(0, 0),
                    vmat = diag(5, nrow = 2, ncol = 2),
                    xwin = c(-100, 100),
                    ywin = c(-100, 100),
                    npix = 1000){
  sx = seq(xwin[1], xwin[2], length.out = npix + 1)[2:(npix + 1)]
  sy = seq(ywin[1], ywin[2], length.out = npix + 1)[2:(npix + 1)]
  rx = (sx[2] - sx[1])
  ry = (sy[2] - sy[1])
  sx = sx - (rx / 2)
  sy = sy - (ry / 2)
  gpts = as.matrix(expand.grid(x = sx, y = sy))
  # dens = mvnfast::dmvn(X = gpts, mu = x, sigma = vmat)
  dens = mvtnorm::dmvnorm(x = gpts, mean = x, sigma = vmat)
  rkb = sum(dens^2 * (rx * ry))
  return(rkb)
}

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
#' Wrapper function that does KDS-Filt, then gets cluster assignments
#'
#' @param sim_pattern Data frame with standard long format CT info for just diseased tissue
#' @param xwin Range of x values
#' @param ywin Range of y values
#' @param bw_method_c Either numeric smoothing bandwidth, or character string for methods (bw_ppl, bw_diggle, bw_CvL)
#' @param thresh Quantile to use for filtering noise
#' @param n_sim_data Number of simulation datasets to use for constructing null distribution
#' @param edge_correct Logical variable denoting if kernel intensity estimates should be edge corrected
#'
#' @return
#' Returns a data frame with x and y coordinates, and the feature or noise classification for each point
#' @export
#'

KDS_filt2 <- function(sim_pattern,
                      xwin,
                      ywin,
                      bw_method_c = 'bw.ppl',
                      thresh = .5,
                      n_sim_data = 20,
                      edge_correct = T){
  sim_pattern_full <- sim_pattern
  # print("denoising")
  if(is.numeric(bw_method_c)){
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    sigma2 <- bw_method_c
    dat_ppp = spatstat::ppp(x = sim_pattern[, 1], y = sim_pattern[, 2], window = owin(xwin, ywin))
    sim_pattern$dens_d <- as.numeric(spatstat::density.ppp(x = dat_ppp,
                                                           sigma = sigma2,
                                                           edge = edge_correct,
                                                           diggle = edge_correct,
                                                           leaveoneout = T,
                                                           at = 'points'))

    res_sim = do.call(c, lapply(1:n_sim_data, FUN = function(i){
      sim_pts_sub = data.frame(x = runif(n = n_pts, min = xwin[1], max = xwin[2]),
                               y = runif(n = n_pts, min = ywin[1], max = ywin[2]))
      dat_ppp_sub = spatstat::ppp(sim_pts_sub$x, sim_pts_sub$y, owin(xwin, ywin))
      dens_ppp_sub = spatstat::density.ppp(x = dat_ppp_sub,
                                           sigma = sigma2,
                                           edge = edge_correct,
                                           diggle = edge_correct,
                                           leaveoneout = T,
                                           at = 'points')
      return(as.numeric(dens_ppp_sub))
    }))

    sim_thresh <- quantile(res_sim, probs = thresh)

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


      sim_pattern$dens_d <- as.numeric(spatstat::density.ppp(x = dat_ppp_sub,
                                                             sigma = sigma2,
                                                             edge = edge_correct,
                                                             diggle = edge_correct,
                                                             leaveoneout = T,
                                                             at = 'points'))


      res_sim = do.call(c, lapply(1:n_sim_data, FUN = function(i){
        sim_pts_sim = data.frame(x = runif(n = n_pts, min = xwin[1], max = xwin[2]),
                                 y = runif(n = n_pts, min = ywin[1], max = ywin[2]))
        dat_ppp_sim = spatstat::ppp(sim_pts_sim$x, sim_pts_sim$y, owin(xwin, ywin))
        dens_ppp_sim = spatstat::density.ppp(x = dat_ppp_sim,
                                             sigma = sigma2,
                                             edge = edge_correct,
                                             diggle = edge_correct,
                                             leaveoneout = T,
                                             at = 'points')
        return(as.numeric(dens_ppp_sim))
      }))

      sim_thresh <- quantile(res_sim, probs = thresh)
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

  else if(bw_method_c == 'Hpi'){
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    dat_ppp <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin(xwin, ywin))
    vcm_hpi = ks::Hpi.diag(sim_pattern[, c('x', 'y')])
    dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp,
                                               varcov = vcm_hpi,
                                               edge = edge_correct,
                                               diggle = edge_correct,
                                               leaveoneout = T,
                                               at = 'points'),
                         error = function(e){return(NA)})
    if(any(is.na(dens_ppp))){
      return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
    }

    sim_pattern$dens_d <- as.numeric(dens_ppp)

    res_sim = do.call(c, lapply(1:n_sim_data, FUN = function(i){
      sim_pts_sim = data.frame(x = runif(n = n_pts, min = xwin[1], max = xwin[2]),
                               y = runif(n = n_pts, min = ywin[1], max = ywin[2]))
      dat_ppp_sim = spatstat::ppp(sim_pts_sim$x, sim_pts_sim$y, owin(xwin, ywin))

      dens_ppp_sim = spatstat::density.ppp(x = dat_ppp_sim,
                                           varcov = vcm_hpi,
                                           edge = edge_correct,
                                           diggle = edge_correct,
                                           leaveoneout = T,
                                           at = 'points')

      return(as.numeric(dens_ppp_sim))
    }))

    sim_thresh <- quantile(res_sim, probs = thresh)

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
      vcm_hpi = ks::Hpi.diag(sim_pattern[, c('x', 'y')])
      dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp_sub,
                                                 varcov = vcm_hpi,
                                                 edge = edge_correct,
                                                 diggle = edge_correct,
                                                 leaveoneout = T,
                                                 at = 'points'),
                           error = function(e){return(NA)})
      if(any(is.na(dens_ppp))){
        return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
      }

      sigma2 <- attributes(dens_ppp)$varcov
      sim_pattern$dens_d <- as.numeric(dens_ppp)

      res_sim = do.call(c, lapply(1:n_sim_data, FUN = function(i){
        sim_pts_sim = data.frame(x = runif(n = n_pts, min = xwin[1], max = xwin[2]),
                                 y = runif(n = n_pts, min = ywin[1], max = ywin[2]))
        dat_ppp_sim = spatstat::ppp(sim_pts_sim$x, sim_pts_sim$y, owin(xwin, ywin))

        dens_ppp_sim = spatstat::density.ppp(x = dat_ppp_sim,
                                             varcov = vcm_hpi,
                                             edge = edge_correct,
                                             diggle = edge_correct,
                                             leaveoneout = T,
                                             at = 'points')


        return(as.numeric(dens_ppp_sim))
      }))

      sim_thresh <- quantile(res_sim, probs = thresh)
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
    else if(bw_method_c == 'bw.scott'){
      bw_method <- bw.scott
    }
    else if(bw_method_c == 'bw.scott.iso'){
      bw_method <- bw.scott.iso
    }
    else{
      bw_method <- bw.CvL
    }
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    dat_ppp <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin(xwin, ywin))
    dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp,
                                               sigma = bw_method,
                                               edge = edge_correct,
                                               diggle = edge_correct,
                                               leaveoneout = T,
                                               at = 'points'),
                         error = function(e){return(NA)})
    if(any(is.na(dens_ppp))){
      return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
    }
    if(bw_method_c == 'bw.scott'){
      sigma2 <- attributes(dens_ppp)$varcov
    }
    else{
      sigma2 <- attributes(dens_ppp)$sigma
    }
    sim_pattern$dens_d <- as.numeric(dens_ppp)

    res_sim = do.call(c, lapply(1:n_sim_data, FUN = function(i){
      sim_pts_sim = data.frame(x = runif(n = n_pts, min = xwin[1], max = xwin[2]),
                               y = runif(n = n_pts, min = ywin[1], max = ywin[2]))
      dat_ppp_sim = spatstat::ppp(sim_pts_sim$x, sim_pts_sim$y, owin(xwin, ywin))
      if(bw_method_c == 'bw.scott'){
        dens_ppp_sim = spatstat::density.ppp(x = dat_ppp_sim,
                                             varcov = sigma2,
                                             edge = edge_correct,
                                             diggle = edge_correct,
                                             leaveoneout = T,
                                             at = 'points')
      }
      else{
        dens_ppp_sim = spatstat::density.ppp(x = dat_ppp_sim,
                                             sigma = sigma2,
                                             edge = edge_correct,
                                             diggle = edge_correct,
                                             leaveoneout = T,
                                             at = 'points')
      }
      return(as.numeric(dens_ppp_sim))
    }))

    sim_thresh <- quantile(res_sim, probs = thresh)

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
      dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp_sub,
                                                 sigma = bw_method,
                                                 edge = edge_correct,
                                                 diggle = edge_correct,
                                                 leaveoneout = T,
                                                 at = 'points'),
                           error = function(e){return(NA)})
      if(any(is.na(dens_ppp))){
        return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
      }
      if(bw_method_c == 'bw.scott'){
        sigma2 <- attributes(dens_ppp)$varcov
      }
      else{
        sigma2 <- attributes(dens_ppp)$sigma
      }

      sim_pattern$dens_d <- as.numeric(dens_ppp)

      res_sim = do.call(c, lapply(1:n_sim_data, FUN = function(i){
        sim_pts_sim = data.frame(x = runif(n = n_pts, min = xwin[1], max = xwin[2]),
                                 y = runif(n = n_pts, min = ywin[1], max = ywin[2]))
        dat_ppp_sim = spatstat::ppp(sim_pts_sim$x, sim_pts_sim$y, owin(xwin, ywin))
        if(bw_method_c == 'bw.scott'){
          dens_ppp_sim = spatstat::density.ppp(x = dat_ppp_sim,
                                               varcov = sigma2,
                                               edge = edge_correct,
                                               diggle = edge_correct,
                                               leaveoneout = T,
                                               at = 'points')
        }
        else{
          dens_ppp_sim = spatstat::density.ppp(x = dat_ppp_sim,
                                               sigma = sigma2,
                                               edge = edge_correct,
                                               diggle = edge_correct,
                                               leaveoneout = T,
                                               at = 'points')
        }
        return(as.numeric(dens_ppp_sim))
      }))

      sim_thresh <- quantile(res_sim, probs = thresh)
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
#' Wrapper function that does KDS-Filt, then gets cluster assignments
#'
#' @param sim_pattern Data frame with standard long format CT info for just diseased tissue
#' @param xwin Range of x values
#' @param ywin Range of y values
#' @param bw_method_c Either numeric smoothing bandwidth, or character string for methods (bw_ppl, bw_diggle, bw_CvL)
#' @param thresh Quantile to use for filtering noise
#' @param n_sim_data Number of simulation datasets to use for constructing null distribution
#' @param edge_correct Logical variable denoting if kernel intensity estimates should be edge corrected
#'
#' @return
#' Returns a data frame with x and y coordinates, and the feature or noise classification for each point
#' @export
#'

KDS_filt2_asym <- function(sim_pattern,
                           xwin,
                           ywin,
                           bw_method_c = 'bw.ppl',
                           thresh = .5,
                           n_sim_data = 20,
                           edge_correct = T){
  sim_pattern_full <- sim_pattern
  # print("denoising")
  if(is.numeric(bw_method_c)){
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    sigma2 <- bw_method_c
    dat_ppp = spatstat::ppp(x = sim_pattern[, 1], y = sim_pattern[, 2], window = owin(xwin, ywin))
    sim_pattern$dens_d <- as.numeric(spatstat::density.ppp(x = dat_ppp,
                                                           sigma = sigma2,
                                                           edge = edge_correct,
                                                           diggle = edge_correct,
                                                           leaveoneout = T,
                                                           at = 'points')) / (n_pts - 1)

    hpi = 1 / LM_sim
    R_k = 1/(2*sqrt(pi))
    hpi_sd = sqrt(hpi * R_k^2 / (sigma2^2 * (n_pts - 1)))

    sim_thresh <- qnorm(p = thresh, mean = hpi, sd = hpi_sd)

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


      sim_pattern$dens_d <- as.numeric(spatstat::density.ppp(x = dat_ppp_sub,
                                                             sigma = sigma2,
                                                             edge = edge_correct,
                                                             diggle = edge_correct,
                                                             leaveoneout = T,
                                                             at = 'points')) / (n_pts - 1)


      hpi = 1 / LM_sim
      R_k = 1/(2*sqrt(pi))
      hpi_sd = sqrt(hpi * R_k^2 / (sigma2^2 * (n_pts - 1)))

      sim_thresh <- qnorm(p = thresh, mean = hpi, sd = hpi_sd)

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
    else if(bw_method_c == 'bw.scott'){
      bw_method <- bw.scott
    }
    else if(bw_method_c == 'bw.scott.iso'){
      bw_method <- bw.scott.iso
    }
    else{
      bw_method <- bw.CvL
    }
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    dat_ppp <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin(xwin, ywin))
    dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp,
                                               sigma = bw_method,
                                               edge = edge_correct,
                                               diggle = edge_correct,
                                               leaveoneout = T,
                                               at = 'points'),
                         error = function(e){return(NA)})
    if(any(is.na(dens_ppp))){
      return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
    }
    sim_pattern$dens_d <- as.numeric(dens_ppp) / (n_pts - 1)

    hpi = 1 / LM_sim
    R_k = 1/(2*sqrt(pi))
    if(bw_method_c == 'bw.scott'){
      sigma2 <- attributes(dens_ppp)$varcov
      hpi_sd = sqrt(hpi * R_k^2 / (sqrt(det(sigma2)) * (n_pts - 1)))
    }
    else{
      sigma2 <- attributes(dens_ppp)$sigma
      hpi_sd = sqrt(hpi * R_k^2 / (sigma2^2 * (n_pts - 1)))
    }

    sim_thresh <- qnorm(p = thresh, mean = hpi, sd = hpi_sd)

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
      dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp_sub,
                                                 sigma = bw_method,
                                                 edge = edge_correct,
                                                 diggle = edge_correct,
                                                 leaveoneout = T,
                                                 at = 'points'),
                           error = function(e){return(NA)})
      if(any(is.na(dens_ppp))){
        return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
      }

      sim_pattern$dens_d <- as.numeric(dens_ppp) / (n_pts - 1)

      hpi = 1 / LM_sim
      R_k = 1/(2*sqrt(pi))
      if(bw_method_c == 'bw.scott'){
        sigma2 <- attributes(dens_ppp)$varcov
        hpi_sd = sqrt(hpi * R_k^2 / (sqrt(det(sigma2)) * (n_pts - 1)))
      }
      else{
        sigma2 <- attributes(dens_ppp)$sigma
        hpi_sd = sqrt(hpi * R_k^2 / (sigma2^2 * (n_pts - 1)))
      }

      sim_thresh <- qnorm(p = thresh, mean = hpi, sd = hpi_sd)
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
#' Wrapper function that does KDS-Filt, then gets cluster assignments
#'
#' @param sim_pattern Data frame with standard long format CT info for just diseased tissue
#' @param xwin Range of x values
#' @param ywin Range of y values
#' @param bw_method_c Either numeric smoothing bandwidth, or character string for methods (bw_ppl, bw_diggle, bw_CvL)
#' @param thresh Quantile to use for filtering noise
#' @param n_sim_data Number of simulation datasets to use for constructing null distribution
#' @param edge_correct Logical variable denoting if kernel intensity estimates should be edge corrected
#'
#' @return
#' Returns a data frame with x and y coordinates, and the feature or noise classification for each point
#' @export
#'

KDS_filt2_asym_bound <- function(sim_pattern,
                           xwin,
                           ywin,
                           bw_method_c = 'bw.ppl',
                           thresh = .5,
                           p_edge = .7,
                           edge_correct = T){
  R_k = 1/(2*sqrt(pi))
  sim_pattern_full <- sim_pattern
  # print("denoising")
  if(is.numeric(bw_method_c)){
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    sigma2 <- bw_method_c
    dat_ppp = spatstat::ppp(x = sim_pattern[, 1], y = sim_pattern[, 2], window = owin(xwin, ywin))
    sim_pattern$dens_d <- as.numeric(spatstat::density.ppp(x = dat_ppp,
                                                           sigma = sigma2,
                                                           edge = edge_correct,
                                                           diggle = edge_correct,
                                                           leaveoneout = T,
                                                           at = 'points')) / (n_pts - 1)

    hpi = 1 / LM_sim

    # sk = do.call(c, lapply(1:npts, function(i){
    #   ret = as.numeric(pmvnorm(lower = c(xwin[1], ywin[1]),
    #                            upper = c(xwin[2], ywin[2]),
    #                            sigma = diag(sigma2, nrow = 2, ncol = 2),
    #                            mean = sim_pattern[i, 1:2]))
    #   return(ret)
    # }))

    sim_thresh = do.call(c, lapply(1:n_pts, function(i){
      sk = as.numeric(pmvnorm(lower = c(xwin[1], ywin[1]),
                              upper = c(xwin[2], ywin[2]),
                              sigma = diag(sigma2^2, nrow = 2, ncol = 2),
                              mean = as.numeric(sim_pattern[i, 1:2])))
      if(sk < p_edge){
        vb = hpi * rk_bound(x = as.numeric(sim_pattern[i, 1:2]),
                            vmat = diag(sigma2^2, nrow = 2, ncol = 2)) / (n_pts - 1) / (sk^2)
      }
      else{
        vb = hpi * R_k^2 / (sigma2^2 * (n_pts - 1))
      }
      ret = qnorm(p = thresh, mean = hpi, sd = sqrt(vb))
    }))

    # R_k = 1/(2*sqrt(pi))
    # hpi_sd = sqrt(hpi * R_k^2 / (sigma2^2 * (n_pts - 1)))
    #
    # sim_thresh <- qnorm(p = thresh, mean = hpi, sd = hpi_sd)

    sim_pattern_sub <- sim_pattern %>%
      dplyr::mutate(sim_thresh = sim_thresh) %>%
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


      sim_pattern$dens_d <- as.numeric(spatstat::density.ppp(x = dat_ppp_sub,
                                                             sigma = sigma2,
                                                             edge = edge_correct,
                                                             diggle = edge_correct,
                                                             leaveoneout = T,
                                                             at = 'points')) / (n_pts - 1)


      hpi = 1 / LM_sim
      sim_thresh = do.call(c, lapply(1:n_pts, function(i){
        sk = as.numeric(pmvnorm(lower = c(xwin[1], ywin[1]),
                                upper = c(xwin[2], ywin[2]),
                                sigma = diag(sigma2^2, nrow = 2, ncol = 2),
                                mean = as.numeric(sim_pattern[i, 1:2])))
        if(sk < p_edge){
          vb = hpi * rk_bound(x = as.numeric(sim_pattern[i, 1:2]),
                              vmat = diag(sigma2^2, nrow = 2, ncol = 2)) / (n_pts - 1) / (sk^2)
        }
        else{
          vb = hpi * R_k^2 / (sigma2^2 * (n_pts - 1))
        }
        ret = qnorm(p = thresh, mean = hpi, sd = sqrt(vb))
      }))

      sim_pattern_sub <- sim_pattern %>%
        dplyr::mutate(sim_thresh = sim_thresh) %>%
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
    else if(bw_method_c == 'bw.scott'){
      bw_method <- bw.scott
    }
    else if(bw_method_c == 'bw.scott.iso'){
      bw_method <- bw.scott.iso
    }
    else{
      bw_method <- bw.CvL
    }
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    dat_ppp <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin(xwin, ywin))
    dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp,
                                               sigma = bw_method,
                                               edge = edge_correct,
                                               diggle = edge_correct,
                                               leaveoneout = T,
                                               at = 'points'),
                         error = function(e){return(NA)})
    if(any(is.na(dens_ppp))){
      return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
    }
    sim_pattern$dens_d <- as.numeric(dens_ppp) / (n_pts - 1)

    hpi = 1 / LM_sim
    # R_k = 1/(2*sqrt(pi))
    if(bw_method_c == 'bw.scott'){
      sigma2 <- attributes(dens_ppp)$varcov
      sim_thresh = do.call(c, lapply(1:n_pts, function(i){
        sk = as.numeric(pmvnorm(lower = c(xwin[1], ywin[1]),
                                upper = c(xwin[2], ywin[2]),
                                sigma = sigma2,
                                mean = as.numeric(sim_pattern[i, 1:2])))
        vb = hpi * rk_bound(x = as.numeric(sim_pattern[i, 1:2]),
                            vmat = sigma2) / (n_pts - 1) / (sk^2)
        ret = qnorm(p = thresh, mean = hpi, sd = sqrt(vb))
      }))
    }
    else{
      sigma2 <- attributes(dens_ppp)$sigma
      sim_thresh = do.call(c, lapply(1:n_pts, function(i){
        sk = as.numeric(pmvnorm(lower = c(xwin[1], ywin[1]),
                                upper = c(xwin[2], ywin[2]),
                                sigma = diag(sigma2^2, nrow = 2, ncol = 2),
                                mean = as.numeric(sim_pattern[i, 1:2])))
        if(sk < p_edge){
          vb = hpi * rk_bound(x = as.numeric(sim_pattern[i, 1:2]),
                              vmat = diag(sigma2^2, nrow = 2, ncol = 2)) / (n_pts - 1) / (sk^2)
        }
        else{
          vb = hpi * R_k^2 / (sigma2^2 * (n_pts - 1))
        }
        ret = qnorm(p = thresh, mean = hpi, sd = sqrt(vb))
      }))
      # hpi_sd = sqrt(hpi * R_k^2 / (sigma2^2 * (n_pts - 1)))
    }

    # sim_thresh <- qnorm(p = thresh, mean = hpi, sd = hpi_sd)

    sim_pattern_sub <- sim_pattern %>%
      dplyr::mutate(sim_thresh = sim_thresh) %>%
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
      dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp_sub,
                                                 sigma = bw_method,
                                                 edge = edge_correct,
                                                 diggle = edge_correct,
                                                 leaveoneout = T,
                                                 at = 'points'),
                           error = function(e){return(NA)})
      if(any(is.na(dens_ppp))){
        return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
      }

      sim_pattern$dens_d <- as.numeric(dens_ppp) / (n_pts - 1)

      hpi = 1 / LM_sim
      # R_k = 1/(2*sqrt(pi))
      if(bw_method_c == 'bw.scott'){
        sigma2 <- attributes(dens_ppp)$varcov
        sim_thresh = do.call(c, lapply(1:n_pts, function(i){
          sk = as.numeric(pmvnorm(lower = c(xwin[1], ywin[1]),
                                  upper = c(xwin[2], ywin[2]),
                                  sigma = sigma2,
                                  mean = as.numeric(sim_pattern[i, 1:2])))
          vb = hpi * rk_bound(x = as.numeric(sim_pattern[i, 1:2]),
                              vmat = sigma2) / (n_pts - 1) / (sk^2)
          ret = qnorm(p = thresh, mean = hpi, sd = sqrt(vb))
        }))
      }
      else{
        sigma2 <- attributes(dens_ppp)$sigma
        sim_thresh = do.call(c, lapply(1:n_pts, function(i){
          sk = as.numeric(pmvnorm(lower = c(xwin[1], ywin[1]),
                                  upper = c(xwin[2], ywin[2]),
                                  sigma = diag(sigma2^2, nrow = 2, ncol = 2),
                                  mean = as.numeric(sim_pattern[i, 1:2])))
          if(sk < p_edge){
            vb = hpi * rk_bound(x = as.numeric(sim_pattern[i, 1:2]),
                                vmat = diag(sigma2^2, nrow = 2, ncol = 2)) / (n_pts - 1) / (sk^2)
          }
          else{
            vb = hpi * R_k^2 / (sigma2^2 * (n_pts - 1))
          }
          ret = qnorm(p = thresh, mean = hpi, sd = sqrt(vb))
        }))
      }


      sim_pattern_sub <- sim_pattern %>%
        dplyr::mutate(sim_thresh = sim_thresh) %>%
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
#' Wrapper function that does KDS-Filt, then gets cluster assignments
#'
#' @param sim_pattern Data frame with standard long format CT info for just diseased tissue
#' @param xwin Range of x values
#' @param ywin Range of y values
#' @param bw_method_c Either numeric smoothing bandwidth, or character string for methods (bw_ppl, bw_diggle, bw_CvL)
#' @param thresh Quantile to use for filtering noise
#' @param n_sim_data Number of simulation datasets to use for constructing null distribution
#' @param edge_correct Logical variable denoting if kernel intensity estimates should be edge corrected
#'
#' @return
#' Returns a data frame with x and y coordinates, and the feature or noise classification for each point
#' @export
#'

KDS_filt2_asym_bound2 <- function(sim_pattern,
                                 xwin,
                                 ywin,
                                 bw_method_c = 'bw.ppl',
                                 thresh = .5,
                                 p_edge = .7,
                                 edge_correct = T){
  R_k = 1/(2*sqrt(pi))
  sim_pattern_full <- sim_pattern
  # print("denoising")
  if(is.numeric(bw_method_c)){
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    sigma2 <- bw_method_c
    dat_ppp = spatstat::ppp(x = sim_pattern[, 1], y = sim_pattern[, 2], window = owin(xwin, ywin))
    dens_ppp = spatstat::density.ppp(x = dat_ppp,
                                     sigma = sigma2,
                                     edge = edge_correct,
                                     diggle = edge_correct,
                                     leaveoneout = T,
                                     at = 'points',
                                     spill = 1)
    sim_pattern$dens_d <- dens_ppp$result / (n_pts - 1)
    sim_pattern$sk = dens_ppp$edg

    hpi = 1 / LM_sim

    sim_thresh = do.call(c, lapply(1:n_pts, function(i){
      sk = sim_pattern$sk[i]
      if(sk < p_edge){
        vb = hpi * rk_bound(x = as.numeric(sim_pattern[i, 1:2]),
                            vmat = diag(sigma2^2, nrow = 2, ncol = 2)) / (n_pts - 1) / (sk^2)
      }
      else{
        vb = hpi * R_k^2 / (sigma2^2 * (n_pts - 1))
      }
      ret = qnorm(p = thresh, mean = hpi, sd = sqrt(vb))
    }))
    # sk = do.call(c, lapply(1:npts, function(i){
    #   ret = as.numeric(pmvnorm(lower = c(xwin[1], ywin[1]),
    #                            upper = c(xwin[2], ywin[2]),
    #                            sigma = diag(sigma2, nrow = 2, ncol = 2),
    #                            mean = sim_pattern[i, 1:2]))
    #   return(ret)
    # }))

    # sim_thresh = do.call(c, lapply(1:n_pts, function(i){
    #   sk = as.numeric(pmvnorm(lower = c(xwin[1], ywin[1]),
    #                           upper = c(xwin[2], ywin[2]),
    #                           sigma = diag(sigma2^2, nrow = 2, ncol = 2),
    #                           mean = as.numeric(sim_pattern[i, 1:2])))
    #   if(sk < .90){
    #     vb = hpi * rk_bound(x = as.numeric(sim_pattern[i, 1:2]),
    #                         vmat = diag(sigma2^2, nrow = 2, ncol = 2)) / (n_pts - 1) / (sk^2)
    #   }
    #   else{
    #     vb = hpi * R_k^2 / (sigma2^2 * (n_pts - 1))
    #   }
    #   ret = qnorm(p = thresh, mean = hpi, sd = sqrt(vb))
    # }))

    # R_k = 1/(2*sqrt(pi))
    # hpi_sd = sqrt(hpi * R_k^2 / (sigma2^2 * (n_pts - 1)))
    #
    # sim_thresh <- qnorm(p = thresh, mean = hpi, sd = hpi_sd)

    sim_pattern_sub <- sim_pattern %>%
      dplyr::mutate(sim_thresh = sim_thresh) %>%
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
      dens_ppp = spatstat::density.ppp(x = dat_ppp_sub,
                                       sigma = sigma2,
                                       edge = edge_correct,
                                       diggle = edge_correct,
                                       leaveoneout = T,
                                       at = 'points',
                                       spill = 1)
      sim_pattern$dens_d <- dens_ppp$result / (n_pts - 1)
      sim_pattern$sk = dens_ppp$edg

      hpi = 1 / LM_sim

      sim_thresh = do.call(c, lapply(1:n_pts, function(i){
        sk = sim_pattern$sk[i]
        if(sk < p_edge){
          vb = hpi * rk_bound(x = as.numeric(sim_pattern[i, 1:2]),
                              vmat = diag(sigma2^2, nrow = 2, ncol = 2)) / (n_pts - 1) / (sk^2)
        }
        else{
          vb = hpi * R_k^2 / (sigma2^2 * (n_pts - 1))
        }
        ret = qnorm(p = thresh, mean = hpi, sd = sqrt(vb))
      }))
      # sim_thresh = do.call(c, lapply(1:n_pts, function(i){
      #   sk = as.numeric(pmvnorm(lower = c(xwin[1], ywin[1]),
      #                           upper = c(xwin[2], ywin[2]),
      #                           sigma = diag(sigma2^2, nrow = 2, ncol = 2),
      #                           mean = as.numeric(sim_pattern[i, 1:2])))
      #   if(sk < .90){
      #     vb = hpi * rk_bound(x = as.numeric(sim_pattern[i, 1:2]),
      #                         vmat = diag(sigma2^2, nrow = 2, ncol = 2)) / (n_pts - 1) / (sk^2)
      #   }
      #   else{
      #     vb = hpi * R_k^2 / (sigma2^2 * (n_pts - 1))
      #   }
      #   ret = qnorm(p = thresh, mean = hpi, sd = sqrt(vb))
      # }))

      sim_pattern_sub <- sim_pattern %>%
        dplyr::mutate(sim_thresh = sim_thresh) %>%
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
    else if(bw_method_c == 'bw.scott'){
      bw_method <- bw.scott
    }
    else if(bw_method_c == 'bw.scott.iso'){
      bw_method <- bw.scott.iso
    }
    else{
      bw_method <- bw.CvL
    }
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    dat_ppp <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin(xwin, ywin))
    dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp,
                                               sigma = bw_method,
                                               edge = edge_correct,
                                               diggle = edge_correct,
                                               leaveoneout = T,
                                               at = 'points',
                                               spill = 1),
                         error = function(e){return(NA)})
    if(any(is.na(dens_ppp))){
      return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
    }
    sim_pattern$dens_d <- as.numeric(dens_ppp$result) / (n_pts - 1)
    sim_pattern$sk = dens_ppp$edg
    hpi = 1 / LM_sim
    # R_k = 1/(2*sqrt(pi))
    if(bw_method_c == 'bw.scott'){
      sigma2 <- dens_ppp$varcov
      sim_thresh = do.call(c, lapply(1:n_pts, function(i){
        sk = sim_pattern$sk[i]
        if(sk < p_edge){
          vb = hpi * rk_bound(x = as.numeric(sim_pattern[i, 1:2]),
                              vmat = sigma2) / (n_pts - 1) / (sk^2)
        }
        else{
          vb = hpi * R_k^2 / (sqrt(det(sigma2)) * (n_pts - 1))
        }
        ret = qnorm(p = thresh, mean = hpi, sd = sqrt(vb))
      }))
    }
    else{
      sigma2 <- dens_ppp$sigma
      sim_thresh = do.call(c, lapply(1:n_pts, function(i){
        sk = sim_pattern$sk[i]
        if(sk < p_edge){
          vb = hpi * rk_bound(x = as.numeric(sim_pattern[i, 1:2]),
                              vmat = diag(sigma2^2, nrow = 2, ncol = 2)) / (n_pts - 1) / (sk^2)
        }
        else{
          vb = hpi * R_k^2 / (sigma2^2 * (n_pts - 1))
        }
        ret = qnorm(p = thresh, mean = hpi, sd = sqrt(vb))
      }))
      # hpi_sd = sqrt(hpi * R_k^2 / (sigma2^2 * (n_pts - 1)))
    }

    # sim_thresh <- qnorm(p = thresh, mean = hpi, sd = hpi_sd)

    sim_pattern_sub <- sim_pattern %>%
      dplyr::mutate(sim_thresh = sim_thresh) %>%
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
      dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp_sub,
                                                 sigma = bw_method,
                                                 edge = edge_correct,
                                                 diggle = edge_correct,
                                                 leaveoneout = T,
                                                 at = 'points',
                                                 spill = 1),
                           error = function(e){return(NA)})
      if(any(is.na(dens_ppp))){
        return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
      }
      sim_pattern$dens_d <- as.numeric(dens_ppp$result) / (n_pts - 1)
      sim_pattern$sk = dens_ppp$edg

      hpi = 1 / LM_sim
      # R_k = 1/(2*sqrt(pi))
      if(bw_method_c == 'bw.scott'){
        sigma2 <- dens_ppp$varcov
        sim_thresh = do.call(c, lapply(1:n_pts, function(i){
          sk = sim_pattern$sk[i]
          if(sk < p_edge){
            vb = hpi * rk_bound(x = as.numeric(sim_pattern[i, 1:2]),
                                vmat = sigma2) / (n_pts - 1) / (sk^2)
          }
          else{
            vb = hpi * R_k^2 / (sqrt(det(sigma2)) * (n_pts - 1))
          }
          ret = qnorm(p = thresh, mean = hpi, sd = sqrt(vb))
        }))
      }
      else{
        sigma2 <- dens_ppp$sigma
        sim_thresh = do.call(c, lapply(1:n_pts, function(i){
          sk = sim_pattern$sk[i]
          if(sk < p_edge){
            vb = hpi * rk_bound(x = as.numeric(sim_pattern[i, 1:2]),
                                vmat = diag(sigma2^2, nrow = 2, ncol = 2)) / (n_pts - 1) / (sk^2)
          }
          else{
            vb = hpi * R_k^2 / (sigma2^2 * (n_pts - 1))
          }
          ret = qnorm(p = thresh, mean = hpi, sd = sqrt(vb))
        }))
      }


      sim_pattern_sub <- sim_pattern %>%
        dplyr::mutate(sim_thresh = sim_thresh) %>%
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
#' Wrapper function that does KDS-Filt, then gets cluster assignments
#'
#' @param sim_pattern Data frame with standard long format CT info for just diseased tissue
#' @param xwin Range of x values
#' @param ywin Range of y values
#' @param bw_method_c Either numeric smoothing bandwidth, or character string for methods (bw_ppl, bw_diggle, bw_CvL)
#' @param thresh Quantile to use for filtering noise
#' @param n_sim_data Number of simulation datasets to use for constructing null distribution
#' @param edge_correct Logical variable denoting if kernel intensity estimates should be edge corrected
#'
#' @return
#' Returns a data frame with x and y coordinates, and the feature or noise classification for each point
#' @export
#'

KDS_filt2_ci <- function(sim_pattern,
                         xwin,
                         ywin,
                         bw_method_c = 'bw.ppl',
                         level = .05,
                         n_sim_data = 20,
                         edge_correct = T){
  sim_pattern = as.data.frame(sim_pattern)
  sim_pattern_full <- sim_pattern
  # print("denoising")
  if(is.numeric(bw_method_c)){
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    sigma2 <- bw_method_c
    bw = sigma2^2
    dat_ppp = spatstat::ppp(x = sim_pattern[, 1],
                            y = sim_pattern[, 2],
                            window = spatstat::owin(xwin, ywin))
    sim_pattern$dens_d <- as.numeric(spatstat::density.ppp(x = dat_ppp,
                                                           sigma = sigma2,
                                                           edge = edge_correct,
                                                           diggle = edge_correct,
                                                           leaveoneout = T,
                                                           at = 'points')) / (n_pts - 1)

    hpi = (n_pts - 1) / n_pts / LM_sim

    ciq = level
    R_k = 1/(2*sqrt(pi))
    z_t = abs(qnorm(level))

    sim_pattern$lci = sim_pattern$dens_d - z_t * sqrt(sim_pattern$dens_d * R_k^2 / (bw * (n_pts - 1)))

    sim_pattern_sub <- sim_pattern %>%
      dplyr::filter(lci > hpi) %>%
      dplyr::select(x, y)

    n_pts_sub <- nrow(sim_pattern_sub)

    while(n_pts != n_pts_sub & n_pts_sub > 0){
      # print('denoising again')
      sim_pattern <- sim_pattern_sub
      n_pts <- nrow(sim_pattern)
      xwin_sub <- range(sim_pattern$x)
      ywin_sub <- range(sim_pattern$y)
      LM_sim_sub <- (xwin_sub[2] - xwin_sub[1]) * (ywin_sub[2] - ywin_sub[1])
      dat_ppp_sub <- spatstat::ppp(sim_pattern$x,
                                   sim_pattern$y,
                                   spatstat::owin(xwin, ywin))


      sim_pattern$dens_d <- as.numeric(spatstat::density.ppp(x = dat_ppp_sub,
                                                             sigma = sigma2,
                                                             edge = edge_correct,
                                                             diggle = edge_correct,
                                                             leaveoneout = T,
                                                             at = 'points')) / (n_pts - 1)



      hpi = (n_pts - 1) / n_pts / LM_sim

      ciq = level
      R_k = (1/(2*sqrt(pi)))
      z_t = abs(qnorm(level))

      sim_pattern$lci = sim_pattern$dens_d - z_t * sqrt(sim_pattern$dens_d * R_k^2 / (bw * (n_pts - 1)))

      sim_pattern_sub <- sim_pattern %>%
        dplyr::filter(lci > hpi) %>%
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
    else if(bw_method_c == 'bw.scott'){
      bw_method <- bw.scott
    }
    else if(bw_method_c == 'bw.scott.iso'){
      bw_method <- bw.scott.iso
    }
    else{
      bw_method <- bw.CvL
    }
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    dat_ppp <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin(xwin, ywin))
    dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp,
                                               sigma = bw_method,
                                               edge = edge_correct,
                                               diggle = edge_correct,
                                               leaveoneout = T,
                                               at = 'points'),
                         error = function(e){return(NA)})
    if(any(is.na(dens_ppp))){
      return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
    }
    if(bw_method_c == 'bw.scott'){
      sigma2 <- attributes(dens_ppp)$varcov
    }
    else{
      sigma2 <- attributes(dens_ppp)$sigma
    }
    sim_pattern$dens_d <- as.numeric(dens_ppp) / (n_pts - 1)
    bw = sigma2^2
    hpi = (n_pts - 1) / n_pts / LM_sim

    ciq = level
    R_k = (1/(2*sqrt(pi)))
    z_t = abs(qnorm(level))

    sim_pattern$lci = sim_pattern$dens_d - z_t * sqrt(sim_pattern$dens_d * R_k^2 / (bw * (n_pts - 1)))

    sim_pattern_sub <- sim_pattern %>%
      dplyr::filter(lci > hpi) %>%
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
      dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp_sub,
                                                 sigma = bw_method,
                                                 edge = edge_correct,
                                                 diggle = edge_correct,
                                                 leaveoneout = T,
                                                 at = 'points'),
                           error = function(e){return(NA)})
      if(any(is.na(dens_ppp))){
        return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
      }
      if(bw_method_c == 'bw.scott'){
        sigma2 <- attributes(dens_ppp)$varcov
      }
      else{
        sigma2 <- attributes(dens_ppp)$sigma
      }
      bw = sigma2^2
      sim_pattern$dens_d <- as.numeric(dens_ppp) / (n_pts - 1)

      hpi = (n_pts - 1) / n_pts / LM_sim

      ciq = level
      R_k = (1/(2*sqrt(pi)))
      z_t = abs(qnorm(level))

      sim_pattern$lci = sim_pattern$dens_d - z_t * sqrt(sim_pattern$dens_d * R_k^2 / (bw * (n_pts - 1)))

      sim_pattern_sub <- sim_pattern %>%
        dplyr::filter(lci > hpi) %>%
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
#' Wrapper function that does KDS-Filt, then gets cluster assignments
#'
#' @param sim_pattern Data frame with standard long format CT info for just diseased tissue
#' @param xwin Range of x values
#' @param ywin Range of y values
#' @param bw_method_c Either numeric smoothing bandwidth, or character string for methods (bw_ppl, bw_diggle, bw_CvL)
#' @param thresh Quantile to use for filtering noise
#' @param n_sim_data Number of simulation datasets to use for constructing null distribution
#' @param edge_correct Logical variable denoting if kernel intensity estimates should be edge corrected
#'
#' @return
#' Returns a data frame with x and y coordinates, and the feature or noise classification for each point
#' @export
#'

KDS_filt2_ci2 <- function(sim_pattern,
                          xwin,
                          ywin,
                          bw_method_c = 'bw.ppl',
                          level = .05,
                          n_sim_data = 20,
                          edge_correct = T){
  sim_pattern = as.data.frame(sim_pattern)
  sim_pattern_full <- sim_pattern
  # print("denoising")
  if(is.numeric(bw_method_c)){
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    sigma2 <- bw_method_c
    bw = sigma2^2
    dat_ppp = spatstat::ppp(x = sim_pattern[, 1],
                            y = sim_pattern[, 2],
                            window = spatstat::owin(xwin, ywin))
    sim_pattern$dens_d <- as.numeric(spatstat::density.ppp(x = dat_ppp,
                                                           sigma = sigma2,
                                                           edge = edge_correct,
                                                           diggle = edge_correct,
                                                           leaveoneout = T,
                                                           at = 'points')) / (n_pts - 1)

    hpi = (n_pts - 1) / n_pts / LM_sim

    ciq = level
    R_k = 1/(2*sqrt(pi))
    z_t = abs(qnorm(level))

    sim_pattern$zs = sim_pattern$dens_d - hpi / sqrt(hpi * R_k^2 / (bw * (n_pts - 1)))

    sim_pattern_sub <- sim_pattern %>%
      dplyr::filter(zs > z_t) %>%
      dplyr::select(x, y)

    n_pts_sub <- nrow(sim_pattern_sub)

    while(n_pts != n_pts_sub & n_pts_sub > 0){
      # print('denoising again')
      sim_pattern <- sim_pattern_sub
      n_pts <- nrow(sim_pattern)
      xwin_sub <- range(sim_pattern$x)
      ywin_sub <- range(sim_pattern$y)
      LM_sim_sub <- (xwin_sub[2] - xwin_sub[1]) * (ywin_sub[2] - ywin_sub[1])
      dat_ppp_sub <- spatstat::ppp(sim_pattern$x,
                                   sim_pattern$y,
                                   spatstat::owin(xwin, ywin))


      sim_pattern$dens_d <- as.numeric(spatstat::density.ppp(x = dat_ppp_sub,
                                                             sigma = sigma2,
                                                             edge = edge_correct,
                                                             diggle = edge_correct,
                                                             leaveoneout = T,
                                                             at = 'points')) / (n_pts - 1)



      hpi = (n_pts - 1) / n_pts / LM_sim

      ciq = level
      R_k = (1/(2*sqrt(pi)))
      z_t = abs(qnorm(level))

      sim_pattern$zs = sim_pattern$dens_d - hpi / sqrt(hpi * R_k^2 / (bw * (n_pts - 1)))

      sim_pattern_sub <- sim_pattern %>%
        dplyr::filter(zs > z_t) %>%
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
    else if(bw_method_c == 'bw.scott'){
      bw_method <- bw.scott
    }
    else if(bw_method_c == 'bw.scott.iso'){
      bw_method <- bw.scott.iso
    }
    else{
      bw_method <- bw.CvL
    }
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    dat_ppp <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin(xwin, ywin))
    dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp,
                                               sigma = bw_method,
                                               edge = edge_correct,
                                               diggle = edge_correct,
                                               leaveoneout = T,
                                               at = 'points'),
                         error = function(e){return(NA)})
    if(any(is.na(dens_ppp))){
      return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
    }
    if(bw_method_c == 'bw.scott'){
      sigma2 <- attributes(dens_ppp)$varcov
    }
    else{
      sigma2 <- attributes(dens_ppp)$sigma
    }
    sim_pattern$dens_d <- as.numeric(dens_ppp) / (n_pts - 1)
    bw = sigma2^2
    hpi = (n_pts - 1) / n_pts / LM_sim

    ciq = level
    R_k = (1/(2*sqrt(pi)))
    z_t = abs(qnorm(level))

    sim_pattern$zs = sim_pattern$dens_d - hpi / sqrt(hpi * R_k^2 / (bw * (n_pts - 1)))

    sim_pattern_sub <- sim_pattern %>%
      dplyr::filter(zs > z_t) %>%
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
      dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp_sub,
                                                 sigma = bw_method,
                                                 edge = edge_correct,
                                                 diggle = edge_correct,
                                                 leaveoneout = T,
                                                 at = 'points'),
                           error = function(e){return(NA)})
      if(any(is.na(dens_ppp))){
        return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
      }
      if(bw_method_c == 'bw.scott'){
        sigma2 <- attributes(dens_ppp)$varcov
      }
      else{
        sigma2 <- attributes(dens_ppp)$sigma
      }
      bw = sigma2^2
      sim_pattern$dens_d <- as.numeric(dens_ppp) / (n_pts - 1)

      hpi = (n_pts - 1) / n_pts / LM_sim

      ciq = level
      R_k = (1/(2*sqrt(pi)))
      z_t = abs(qnorm(level))

      sim_pattern$zs = sim_pattern$dens_d - hpi / sqrt(hpi * R_k^2 / (bw * (n_pts - 1)))

      sim_pattern_sub <- sim_pattern %>%
        dplyr::filter(zs > z_t) %>%
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
#' Wrapper function that does KDS-Filt, then gets cluster assignments
#'
#' @param sim_pattern Data frame with standard long format CT info for just diseased tissue
#' @param xwin Range of x values
#' @param ywin Range of y values
#' @param bw_method_c Either numeric smoothing bandwidth, or character string for methods (bw_ppl, bw_diggle, bw_CvL)
#' @param thresh Quantile to use for filtering noise
#' @param n_sim_data Number of simulation datasets to use for constructing null distribution
#' @param edge_correct Logical variable denoting if kernel intensity estimates should be edge corrected
#'
#' @return
#' Returns a data frame with x and y coordinates, and the feature or noise classification for each point
#' @export
#'

KDS_filt2_r <- function(sim_pattern,
                        xwin,
                        ywin,
                        bw_method_c = 'bw.ppl',
                        thresh = .5,
                        n_sim_data = 20,
                        edge_correct = T){
  sim_pattern_full <- sim_pattern
  # print("denoising")
  if(is.numeric(bw_method_c)){
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    sigma2 <- bw_method_c
    dat_ppp = spatstat::ppp(x = sim_pattern[, 1], y = sim_pattern[, 2], window = owin(xwin, ywin))
    sim_pattern$dens_d <- as.numeric(spatstat::density.ppp(x = dat_ppp,
                                                           sigma = sigma2,
                                                           edge = edge_correct,
                                                           diggle = edge_correct,
                                                           leaveoneout = T,
                                                           at = 'points'))

    res_sim = do.call(c, lapply(1:n_sim_data, FUN = function(i){
      n_pts_rand = rpois(n = 1, lambda = n_pts)
      sim_pts_sim = data.frame(x = runif(n = n_pts_rand, min = xwin[1], max = xwin[2]),
                               y = runif(n = n_pts_rand, min = ywin[1], max = ywin[2]))
      dat_ppp_sub = spatstat::ppp(sim_pts_sim$x, sim_pts_sim$y, owin(xwin, ywin))
      dens_ppp_sub = spatstat::density.ppp(x = dat_ppp_sub,
                                           sigma = sigma2,
                                           edge = edge_correct,
                                           diggle = edge_correct,
                                           leaveoneout = T,
                                           at = 'points')
      return(as.numeric(dens_ppp_sub))
    }))

    sim_thresh <- quantile(res_sim, probs = thresh)

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


      sim_pattern$dens_d <- as.numeric(spatstat::density.ppp(x = dat_ppp_sub,
                                                             sigma = sigma2,
                                                             edge = edge_correct,
                                                             diggle = edge_correct,
                                                             leaveoneout = T,
                                                             at = 'points'))


      res_sim = do.call(c, lapply(1:n_sim_data, FUN = function(i){
        n_pts_rand = rpois(n = 1, lambda = n_pts)
        sim_pts_sim = data.frame(x = runif(n = n_pts_rand, min = xwin[1], max = xwin[2]),
                                 y = runif(n = n_pts_rand, min = ywin[1], max = ywin[2]))
        dat_ppp_sim = spatstat::ppp(sim_pts_sim$x, sim_pts_sim$y, owin(xwin, ywin))
        dens_ppp_sim = spatstat::density.ppp(x = dat_ppp_sim,
                                             sigma = sigma2,
                                             edge = edge_correct,
                                             diggle = edge_correct,
                                             leaveoneout = T,
                                             at = 'points')
        return(as.numeric(dens_ppp_sim))
      }))

      sim_thresh <- quantile(res_sim, probs = thresh)
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
    else if(bw_method_c == 'bw.scott'){
      bw_method <- bw.scott
    }
    else if(bw_method_c == 'bw.scott.iso'){
      bw_method <- bw.scott.iso
    }
    else if(bw_method_c == 'null'){
      bw_method = NULL
    }
    else{
      bw_method <- bw.CvL
    }
    LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    dat_ppp <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin(xwin, ywin))
    dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp,
                                               sigma = bw_method,
                                               edge = edge_correct,
                                               diggle = edge_correct,
                                               leaveoneout = T,
                                               at = 'points'),
                         error = function(e){return(NA)})
    if(any(is.na(dens_ppp))){
      return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
    }
    if(bw_method_c == 'bw.scott'){
      sigma2 <- attributes(dens_ppp)$varcov
    }
    else{
      sigma2 <- attributes(dens_ppp)$sigma
    }
    sim_pattern$dens_d <- as.numeric(dens_ppp)

    res_sim = do.call(c, lapply(1:n_sim_data, FUN = function(i){
      n_pts_rand = rpois(n = 1, lambda = n_pts)
      sim_pts_sim = data.frame(x = runif(n = n_pts_rand, min = xwin[1], max = xwin[2]),
                               y = runif(n = n_pts_rand, min = ywin[1], max = ywin[2]))
      dat_ppp_sim = spatstat::ppp(sim_pts_sim$x, sim_pts_sim$y, owin(xwin, ywin))
      if(bw_method_c == 'bw.scott'){
        dens_ppp_sim = spatstat::density.ppp(x = dat_ppp_sim,
                                             varcov = sigma2,
                                             edge = edge_correct,
                                             diggle = edge_correct,
                                             leaveoneout = T,
                                             at = 'points')
      }
      else{
        dens_ppp_sim = spatstat::density.ppp(x = dat_ppp_sim,
                                             sigma = sigma2,
                                             edge = edge_correct,
                                             diggle = edge_correct,
                                             leaveoneout = T,
                                             at = 'points')
      }
      return(as.numeric(dens_ppp_sim))
    }))

    sim_thresh <- quantile(res_sim, probs = thresh)

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
      dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp_sub,
                                                 sigma = bw_method,
                                                 edge = edge_correct,
                                                 diggle = edge_correct,
                                                 leaveoneout = T,
                                                 at = 'points'),
                           error = function(e){return(NA)})
      if(any(is.na(dens_ppp))){
        return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
      }
      if(bw_method_c == 'bw.scott'){
        sigma2 <- attributes(dens_ppp)$varcov
      }
      else{
        sigma2 <- attributes(dens_ppp)$sigma
      }

      sim_pattern$dens_d <- as.numeric(dens_ppp)

      res_sim = do.call(c, lapply(1:n_sim_data, FUN = function(i){
        n_pts_rand = rpois(n = 1, lambda = n_pts)
        sim_pts_sim = data.frame(x = runif(n = n_pts_rand, min = xwin[1], max = xwin[2]),
                                 y = runif(n = n_pts_rand, min = ywin[1], max = ywin[2]))
        dat_ppp_sim = spatstat::ppp(sim_pts_sim$x, sim_pts_sim$y, owin(xwin, ywin))
        if(bw_method_c == 'bw.scott'){
          dens_ppp_sim = spatstat::density.ppp(x = dat_ppp_sim,
                                               varcov = sigma2,
                                               edge = edge_correct,
                                               diggle = edge_correct,
                                               leaveoneout = T,
                                               at = 'points')
        }
        else{
          dens_ppp_sim = spatstat::density.ppp(x = dat_ppp_sim,
                                               sigma = sigma2,
                                               edge = edge_correct,
                                               diggle = edge_correct,
                                               leaveoneout = T,
                                               at = 'points')
        }
        return(as.numeric(dens_ppp_sim))
      }))

      sim_thresh <- quantile(res_sim, probs = thresh)
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
    else if(bw_method_c == 'bw.scott'){
      bw_method <- bw.scott
    }
    else if(bw_method_c == 'bw.scott.iso'){
      bw_method <- bw.scott.iso
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
KDS_filt2_pix <- function(sim_pattern,
                          all_points,
                          bw_method_c = 'bw.ppl',
                          thresh = .5,
                          n_sim_data = 20,
                          con_hull = F,
                          edge_correct = T){
  sim_pattern_full <- sim_pattern
  if(con_hull == T){
    owin_obs = spatstat::convexhull.xy(x = all_points$x, y = all_points$y)
  }
  else{
    xwin <- range(all_points[, 1])
    ywin <- range(all_points[, 2])
    owin_obs = owin(xwin, ywin)
  }
  # print("denoising")
  if(is.numeric(bw_method_c)){
    LM_sim <- nrow(all_points)

    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    sigma2 <- bw_method_c
    dat_ppp <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin_obs)
    sim_pattern$dens_d <- as.numeric(spatstat::density.ppp(x = dat_ppp,
                                                           sigma = sigma2,
                                                           edge = edge_correct,
                                                           diggle = edge_correct,
                                                           leaveoneout = T,
                                                           at = 'points'))

    res_sim = do.call(c, lapply(1:n_sim_data, FUN = function(i){
      sim_pts_sub = all_points[sample(1:nrow(all_points), size = n_pts, replace = F), ]
      dat_ppp_sub = spatstat::ppp(sim_pts_sub[, 1], sim_pts_sub[, 2], owin_obs)
      dens_ppp_sub = spatstat::density.ppp(x = dat_ppp_sub,
                                           sigma = sigma2,
                                           edge = edge_correct,
                                           diggle = edge_correct,
                                           leaveoneout = T,
                                           at = 'points')
      return(as.numeric(dens_ppp_sub))
    }))

    sim_thresh <- quantile(res_sim, probs = thresh)

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
      dat_ppp_sub <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin_obs)


      sim_pattern$dens_d <- as.numeric(spatstat::density.ppp(x = dat_ppp_sub,
                                                             sigma = sigma2,
                                                             edge = edge_correct,
                                                             diggle = edge_correct,
                                                             leaveoneout = T,
                                                             at = 'points'))


      res_sim = do.call(c, lapply(1:n_sim_data, FUN = function(i){
        sim_pts_sub = data.frame(x = runif(n = n_pts, min = xwin[1], max = xwin[2]),
                                 y = runif(n = n_pts, min = ywin[1], max = ywin[2]))
        dat_ppp_sub = spatstat::ppp(sim_pts_sub$x, sim_pts_sub$y, owin_obs)
        dens_ppp_sub = spatstat::density.ppp(x = dat_ppp_sub,
                                             sigma = sigma2,
                                             edge = edge_correct,
                                             diggle = edge_correct,
                                             leaveoneout = T,
                                             at = 'points')
        return(as.numeric(dens_ppp_sub))
      }))

      sim_thresh <- quantile(res_sim, probs = thresh)

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
  else if(bw_method_c == 'Hpi'){
    # LM_sim <- (xwin[2] - xwin[1]) * (ywin[2] - ywin[1])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    dat_ppp <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin_obs)
    vcm_hpi = ks::Hpi.diag(sim_pattern[, c('x', 'y')])
    dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp,
                                               varcov = vcm_hpi,
                                               edge = edge_correct,
                                               diggle = edge_correct,
                                               leaveoneout = T,
                                               at = 'points'),
                         error = function(e){return(NA)})
    if(any(is.na(dens_ppp))){
      return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
    }

    sim_pattern$dens_d <- as.numeric(dens_ppp)

    res_sim = do.call(c, lapply(1:n_sim_data, FUN = function(i){
      sim_pts_sim = all_points[sample(1:nrow(all_points), size = n_pts, replace = F), ]
      dat_ppp_sim = spatstat::ppp(sim_pts_sim$x, sim_pts_sim$y, owin_obs)

      dens_ppp_sim = spatstat::density.ppp(x = dat_ppp_sim,
                                           varcov = vcm_hpi,
                                           edge = edge_correct,
                                           diggle = edge_correct,
                                           leaveoneout = T,
                                           at = 'points')

      return(as.numeric(dens_ppp_sim))
    }))

    sim_thresh <- quantile(res_sim, probs = thresh)

    sim_pattern_sub <- sim_pattern %>%
      dplyr::filter(dens_d > sim_thresh) %>%
      dplyr::select(x, y)
    n_pts_sub <- nrow(sim_pattern_sub)

    while(n_pts != n_pts_sub & n_pts_sub > 1){
      # print('denoising again')
      sim_pattern <- sim_pattern_sub
      n_pts <- nrow(sim_pattern)
      # xwin_sub <- range(sim_pattern$x)
      # ywin_sub <- range(sim_pattern$y)
      # LM_sim_sub <- (xwin_sub[2] - xwin_sub[1]) * (ywin_sub[2] - ywin_sub[1])
      dat_ppp_sub <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin_obs)
      vcm_hpi = ks::Hpi.diag(sim_pattern[, c('x', 'y')])
      dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp_sub,
                                                 varcov = vcm_hpi,
                                                 edge = edge_correct,
                                                 diggle = edge_correct,
                                                 leaveoneout = T,
                                                 at = 'points'),
                           error = function(e){return(NA)})
      if(any(is.na(dens_ppp))){
        return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
      }

      sigma2 <- attributes(dens_ppp)$varcov
      sim_pattern$dens_d <- as.numeric(dens_ppp)

      res_sim = do.call(c, lapply(1:n_sim_data, FUN = function(i){
        sim_pts_sim = all_points[sample(1:nrow(all_points), size = n_pts, replace = F), ]
        dat_ppp_sim = spatstat::ppp(sim_pts_sim$x, sim_pts_sim$y, owin_obs)

        dens_ppp_sim = spatstat::density.ppp(x = dat_ppp_sim,
                                             varcov = vcm_hpi,
                                             edge = edge_correct,
                                             diggle = edge_correct,
                                             leaveoneout = T,
                                             at = 'points')


        return(as.numeric(dens_ppp_sim))
      }))

      sim_thresh <- quantile(res_sim, probs = thresh)
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
    else if(bw_method_c == 'bw.scott'){
      bw_method <- bw.scott
    }
    else if(bw_method_c == 'bw.scott.iso'){
      bw_method <- bw.scott.iso
    }
    else{
      bw_method <- bw.CvL
    }
    LM_sim <- nrow(all_points)
    xwin <- range(all_points[, 1])
    ywin <- range(all_points[, 2])
    n_pts <- nrow(sim_pattern)
    n_pts_orig <- n_pts
    dat_ppp <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin_obs)
    dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp,
                                               sigma = bw_method,
                                               edge = edge_correct,
                                               diggle = edge_correct,
                                               leaveoneout = T,
                                               at = 'points'),
                         error = function(e){return(NA)})
    if(any(is.na(dens_ppp))){
      return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
    }
    if(bw_method_c == 'bw.scott'){
      sigma2 <- attributes(dens_ppp)$varcov
    }
    else{
      sigma2 <- attributes(dens_ppp)$sigma
    }
    sim_pattern$dens_d <- as.numeric(dens_ppp)

    res_sim = do.call(c, lapply(1:n_sim_data, FUN = function(i){
      sim_pts_sub = all_points[sample(1:nrow(all_points), size = n_pts, replace = F), ]
      dat_ppp_sub = spatstat::ppp(sim_pts_sub[, 1], sim_pts_sub[, 2], owin_obs)
      if(bw_method_c == 'bw.scott'){
        dens_ppp_sub = spatstat::density.ppp(x = dat_ppp_sub,
                                             varcov = sigma2,
                                             edge = edge_correct,
                                             diggle = edge_correct,
                                             leaveoneout = T,
                                             at = 'points')
      }
      else{
        dens_ppp_sub = spatstat::density.ppp(x = dat_ppp_sub,
                                             sigma = sigma2,
                                             edge = edge_correct,
                                             diggle = edge_correct,
                                             leaveoneout = T,
                                             at = 'points')
      }
      return(as.numeric(dens_ppp_sub))
    }))

    sim_thresh <- quantile(res_sim, probs = thresh)

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
      dat_ppp_sub <- spatstat::ppp(sim_pattern$x, sim_pattern$y, owin_obs)
      dens_ppp <- tryCatch(spatstat::density.ppp(x = dat_ppp_sub,
                                                 sigma = bw_method,
                                                 edge = edge_correct,
                                                 diggle = edge_correct,
                                                 leaveoneout = T,
                                                 at = 'points'),
                           error = function(e){return(NA)})
      if(any(is.na(dens_ppp))){
        return(sim_pattern_full %>% dplyr::mutate(type = 'noise'))
      }
      if(bw_method_c == 'bw.scott'){
        sigma2 <- attributes(dens_ppp)$varcov
      }
      else{
        sigma2 <- attributes(dens_ppp)$sigma
      }

      sim_pattern$dens_d <- as.numeric(dens_ppp)

      res_sim = do.call(c, lapply(1:n_sim_data, FUN = function(i){
        sim_pts_sub = all_points[sample(1:nrow(all_points), size = n_pts, replace = F), ]
        dat_ppp_sub = spatstat::ppp(sim_pts_sub[, 1], sim_pts_sub[, 2], owin_obs)
        if(bw_method_c == 'bw.scott'){
          dens_ppp_sub = spatstat::density.ppp(x = dat_ppp_sub,
                                               varcov = sigma2,
                                               edge = edge_correct,
                                               diggle = edge_correct,
                                               leaveoneout = T,
                                               at = 'points')
        }
        else{
          dens_ppp_sub = spatstat::density.ppp(x = dat_ppp_sub,
                                               sigma = sigma2,
                                               edge = edge_correct,
                                               diggle = edge_correct,
                                               leaveoneout = T,
                                               at = 'points')
        }
        return(as.numeric(dens_ppp_sub))
      }))

      sim_thresh <- quantile(res_sim, probs = thresh)
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



