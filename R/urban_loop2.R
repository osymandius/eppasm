setwd("~/Documents/GitHub/EPP")

library(dplyr)
library(magrittr)

devtools::install_github("mrc-ide/eppasm@new-master")
library(eppasm)
install.packages("testthat")
devtools::load_all()

files <- c("~/Documents/2018-12 2020 targets/Uganda  2030.pjnz")
# prop_75 <- vector("list", length(files))

inc_75 <- rep(list(vector("list", 100)), length(files))
mort_75 <- rep(list(vector("list", 100)), length(files))

survey_years <- c("2020", "2022")
inc_readout <- vector("list", length(survey_years))
names(inc_readout)[1:length(survey_years)] <- survey_years

num_surv_iter <- 2

for (j in 1:length(files)) {

  ## Preparing EPP-ASM inputs
  # pjnz <- "~/Documents/2018-12 2020 targets/Botswana_ 2030.PJNZ"
  # pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="eppasm")
  pjnz <- files[1]
  bw <- prepare_spec_fit(pjnz, proj.end=2029.5)
  clean_hhs <- attr(bw$Urban, "eppd")$hhs
  
  fp <- prepare_directincid(pjnz)
  
  ####### Direct incidence #########
  
  fp_new <- prepare_directincid(pjnz)
  
  x <- c(2010, 2020)
  y <- c(fp[["incidinput"]][41], fp[["incidinput"]][41]*0.25)
  
  test2 <- data.frame(approx(x, y, xout=2010:2020))
  
  fp_new[["incidinput"]][41:51] <- test2$y
  fp_new[["incidinput"]][52:length(fp_new[["incidinput"]])] <- fp_new[["incidinput"]][51]
  
  plot(fp_new[["incidinput"]])
  
  mod_new <- simmod(fp_new)
  
  ########## Simulate HHS #############
  
  for (k in 1:length(survey_years)) {
  
    
    # sim_2020 <- data.frame(rbinom(100, 3000, prev(mod_new)[51])/3000) # 51st item in prev(mod_new) is the 2020 prevalence generated from 2020 incidence at 25% of 2010 incidence.
    sim_survey <- data.frame(matrix(ncol=length(survey_years)))
    colnames(sim_survey) <- survey_years
    sim_survey[1:100,1] <- data.frame(rbinom(100, 3000, prev(mod_new)[as.numeric(survey_years[1])-1969])/3000) 
    sim_survey[1:100,2] <- data.frame(rbinom(100, 3000, prev(mod_new)[as.numeric(survey_years[2])-1969])/3000) 
    sim_survey[1:100,3] <- data.frame(rbinom(100, 3000, prev(mod_new)[as.numeric(survey_years[3])-1969])/3000) 
    sim_survey[1:100,4] <- data.frame(rbinom(100, 3000, prev(mod_new)[as.numeric(survey_years[4])-1969])/3000) 
    
    # prop_75 <- data.frame("prop_above_75" = 0)
    sim_hhs <- vector("list", 100)
    
    for (i in 1:num_surv_iter) {

      attr(bw$Urban, "eppd")$hhs <- clean_hhs %>%
          rbind(data.frame("year" = as.numeric(survey_years[1]), "sex" = "both", "agegr"="15-49", "n" = 3000,"prev" = sim_survey[i,1], "se"=sqrt((sim_survey[i,1]*(1-sim_survey[i,1]))/3000),"deff"= 2, "deff_approx"= 2, "used"=TRUE)) %>%
          rbind(data.frame("year" = as.numeric(survey_years[2]), "sex" = "both", "agegr"="15-49", "n" = 3000,"prev" = sim_survey[i,2], "se"=sqrt((sim_survey[i,2]*(1-sim_survey[i,2]))/3000),"deff"= 2, "deff_approx"= 2, "used"=TRUE)) %>%
          rbind(data.frame("year" = as.numeric(survey_years[3]), "sex" = "both", "agegr"="15-49", "n" = 3000,"prev" = sim_survey[i,3], "se"=sqrt((sim_survey[i,3]*(1-sim_survey[i,3]))/3000),"deff"= 2, "deff_approx"= 2, "used"=TRUE)) %>%
          rbind(data.frame("year" = as.numeric(survey_years[4]), "sex" = "both", "agegr"="15-49", "n" = 3000,"prev" = sim_survey[i,4], "se"=sqrt((sim_survey[i,4]*(1-sim_survey[i,4]))/3000),"deff"= 2, "deff_approx"= 2, "used"=TRUE)) %>%
          filter(year<=as.numeric(survey_years[k]))

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

      ## Calculate likelihood
      ll(theta_ur, fp, likdat)

      ## Components of likelihood calculation
      ll_hhsage(mod, likdat$hhs.dat)

      ll_ancsite(mod, fp_par,
                 coef = c(fp_par$ancbias, fp_par$ancrtsite.beta),
                 vinfl = fp_par$v.infl,
                 dat = likdat$ancsite.dat)

      ll_ancrtcens(mod, likdat$ancrtcens.dat, fp_par)

      ## Fitting the EPP-ASM model
      bwfit <- list()

      bwfit$Urban <- fitmod(bw$Urban, eppmod = "rhybrid", rw_start = 2005,
                              B0=1e4, B=1e2, opt_iter = 1, number_k=50)
      # bwfit$Rural <- fitmod(bw$Rural, eppmod = "rhybrid", rw_start = 2005,
      #                        B0=1e3, B=1e2, opt_iter = 1, number_k=50)

        #' When fitting, the random-walk based models only simulate through the end of the
        #  data period. The `extend_projection()` function extends the random walk for $r(t)$
        #  through the end of the projection period.

      bwfit <- lapply(bwfit, extend_projection, proj_years = 60)


        ## Simulating model outptus

      bwout <- Map(tidy_output, bwfit, "r-hybrid", attr(bw$Urban, "eppd")$country, names(bwfit))

      bwaggr <- aggr_specfit(bwfit)

      reduction <- data.frame(
        "incidence_red" = 1-sapply(1:3000, function(x) {
            attr(bwaggr[[x]], "incid15to49")[as.numeric(survey_years[k])-1969]/attr(bwaggr[[x]], "incid15to49")[41]
        }),
        "mort_red" = sapply(1:3000, function(x) {

          foo <- attr(bwaggr[[x]], "hivdeaths") %>%
            aperm(c(3,2,1)) %>%
            rowSums(dim=1)

          mort_red <- 1-(foo[as.numeric(survey_years[k])-1969]/foo[41])

        }))

      inc_75[[j]][i] <- sum(reduction$incidence_red >= 0.75)/3000
      mort_75[[j]][i] <- sum(reduction$mort_red >= 0.75)/3000
    } # end of i loop
  
    names(inc_75)[j] <- attr(bw$Urban, "eppd")$country
    names(mort_75)[j] <- attr(bw$Urban, "eppd")$country

   inc_test <- unlist(inc_75)
   # mort_test <- unlist(mort_75)

   if (sum(inc_test/num_surv_iter)<0.5) {
    inc_readout[k] <- "No"
   } else {
    inc_readout[k] <- "Yes"
    break
   }
    
  } #end of k loop

} # end of j loop

test <- Map(tidy_edit, bwfit, "r-hybrid", attr(bw$Urban, "eppd")$country, names(bwfit))

head(test2)
  ggplot(data=attr(bw$Urban, "eppd")$hhs, aes(x=year, y=prev)) +
    geom_point() +
    geom_line(data=bwout$Urban$core %>% filter(indicator=="prev"), aes(x=year, y=mean))
  
  ggplot()+
    geom_line(data=bwout$Urban$core %>% filter(indicator=="incid" & year>=2010), aes(x=year, y=mean))

  # data=data.frame("prev"=attr(bwaggr[[350]], "prev15to49"), "year" = 1970:2029)

  tidy_edit <- function(fit, modlab, country=NA, eppregion=NA, ancsite = TRUE){
    
    idvars <- data.frame(country = country,
                         eppregion = eppregion,
                         modlab = modlab)
    print(paste(country, eppregion))
    
    ss <- fit$fp$ss
    
    ## simulate model projections
    param_list <- lapply(seq_len(nrow(fit$resample)), function(ii) fnCreateParam(fit$resample[ii,], fit$fp))
    
    
    fp_list <- lapply(param_list, function(par) update(fit$fp, list=par))
    mod_list <- lapply(fp_list, simmod)
    
    ## Site-level ANC outputs

      b_site <- Map(sample_b_site, mod_list, fp_list,
                    list(fit$likdat$ancsite.dat), resid = FALSE)
      
      b_site_sigma <- sapply(b_site, anclik::sample.sigma2)
      
      b_site_df <- estci2(do.call(cbind, b_site))
      ancsite_b <- data.frame(idvars, site = rownames(b_site_df), b_site_df)
      
      newdata <- expand.grid(site = unique(fit$likdat$ancsite.dat$df$site),
                             year = 1985:2029,
                             type = "ancss",
                             age = 15,
                             agspan = 35,
                             n = 300)
      new_df <- ancsite_pred_df(newdata, fit$fp)
      
      ancsite_pred <- mapply(sample_ancsite_pred, mod_list, fp_list,
                             b_site = b_site,
                             MoreArgs = list(newdata = new_df))
      
      ancsite_pred <- pnorm(ancsite_pred)
       
       
      ancsite_pred <- data.frame(newdata, estci2(ancsite_pred))
      # 
      ancsite_pred <- merge(ancsite_pred,
                            fit$likdat$ancsite.dat$df[c("site", "year", "type", "age", "agspan", "n", "prev", "pstar", "W", "v")],
                            by = c("site", "year", "type", "age", "agspan"),
                            suffixes = c("_sim", "_obs"), all.x=TRUE)
      # 
      ancsite_pred <- data.frame(idvars, ancsite_pred)
    
    # out <- list(core = core,
    #             ageprevdat=ageprevdat,
    #             pregprev = pregprev,
    #             agegr3prev = agegr3prev,
    #             ancsite_pred = ancsite_pred,
    #             ancsite_b = ancsite_b)
    ## ageincid = ageincid,
    ## ageinfections = ageinfections,
    ## relincid = relincid)
      
      return(ancsite_pred)
  }
