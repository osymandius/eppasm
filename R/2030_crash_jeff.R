library(eppasm)
library(dplyr)
library(magrittr)
library(broom)

############### FUNCTIONS TO RUN ##############

tidy_test <- function(bw, ancsite=TRUE) {
  
  idvars <- data.frame(country = "country",
                       eppregion = "eppregion",
                       modlab = "modlab")
  
  ss <- fp$ss
  
  param_list <- list(fnCreateParam(theta_ur, fp))
  fp_list <- lapply(param_list, function(par) update(fp, list=par))
  mod_list <- lapply(fp_list, simmod)
  
  b_site <- Map(sample_b_site, mod_list, fp_list, list(likdat$ancsite.dat), resid = FALSE)
  
  b_site_sigma <- sapply(b_site, anclik::sample.sigma2)
  
  b_site_df <- estci2(do.call(cbind, b_site))
  ancsite_b <- data.frame(idvars, site = rownames(b_site_df), b_site_df)
  
  newdata <- expand.grid(site = unique(likdat$ancsite.dat$df$site),
                         year = 1985:2020,
                         type = "ancss",
                         age = 15,
                         agspan = 35,
                         n = 300)
  new_df <- ancsite_pred_df(newdata, fp)
  
  ancsite_pred <- mapply(sample_ancsite_pred, mod_list, fp_list,
                         b_site = b_site,
                         MoreArgs = list(newdata = new_df))
  
  ancsite_pred <- pnorm(ancsite_pred)
  
  
  ancsite_pred <- data.frame(newdata, estci2(ancsite_pred))
  #
  ancsite_pred <- merge(ancsite_pred,
                        likdat$ancsite.dat$df[c("site", "year", "type", "age", "agspan", "n", "prev", "pstar", "W", "v")],
                        by = c("site", "year", "type", "age", "agspan"),
                        suffixes = c("_sim", "_obs"), all.x=TRUE)
  #
  ancsite_pred <- data.frame(idvars, ancsite_pred)
  
  return(ancsite_pred)
  
} 

min_year <- function(bw) {
  
  get_min_year <- tidy_test(bw) %>%
    group_by(site) %>%
    filter(!is.na(n_obs)) %>%
    summarise(min = min(year))
  
  years <- unlist(get_min_year[2])
  sites <- as.character(get_min_year$site)
  
  
  nest_sites <- tidy_test(bw) %>%
    group_by(site) %>%
    nest()
  
  filter_years <- function(a, b) {
    a%<>%
      filter(year>=b)
  }
  
  filtered_min_year <- map2(nest_sites$data, years, filter_years)
  
  for (i in 1:length(filtered_min_year)) {
    filtered_min_year[[i]]$site <- sites[i]
  }
  
  filtered_min_year <- data.frame(do.call(rbind, filtered_min_year))
  
  filtered_min_year %<>%
    mutate(prev = ifelse(is.na(prev) & year>2011, mean, prev)) %>%
    mutate(n_obs = ifelse(year>2011, 300, n_obs))
  
  return(filtered_min_year)
  
}

###############################################

files <- c("~/Documents/Data/Spectrum files/2018 final/SSA/Botswana_ 2018 updated ART.PJNZ", "~/Documents/2018-12 2020 targets/Botswana_ 2030.PJNZ")

pjnz <- files[1]
bw <- prepare_spec_fit(pjnz, proj.end=2021.5)

fp <- prepare_directincid(pjnz)

theta_ur <- c(-0.63758, -2.76655, -1.26204, 1996.65945, 0.00778, 0.05195,
              0.05103, 0.032, 0.01765, 0.01154, -0.00028, 0.01627, -0.00051,
              0.01439, -0.00937, -0.01135, 0.03692, 0.14959, 0.00803, 0.02424,
              -0.03548, 3.65223, -0.02515, -4.74563, 0.26259, -6.90124, 0.01583)

fp <- attr(bw$Urban, "specfp")

fp <- prepare_rhybrid(fp, rw_start = 2005, rw_dk = 1)

## Set some flags that are set in fitmod(), (later improve this code...)
fp$ancsitedata <- TRUE
fp$ancrt <- "both"
fp$logitiota <- TRUE
fp$rw_start <- 2005
fp$incidmod <- "eppspectrum"

param <- fnCreateParam(theta_ur, fp)
fp_par <- update(fp, list = param)


## Simulate the model once.

mod <- simmod(fp_par)
## Prepare likelihood and calculate the likelihood once

## Prior
lprior(theta_ur, fp)

## Prepare likelihood data
likdat <- prepare_likdat(attr(bw$Urban, "eppd"), fp)

#########################################

plot4 <- ggplot()+
  geom_line(data=min_year(bw) %>% filter(!is.na(prev)) %>% filter(year<2012), aes(x=year, y=prev, color=site)) +
  geom_line(data=min_year(bw) %>% filter(!is.na(prev)) %>% filter(year>2010), aes(x=year, y=prev, color=site), linetype = 2)

grid.arrange(plot1, plot2, plot3, plot4)

