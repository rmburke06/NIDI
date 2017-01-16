# My functions

## Function to get categorical variables by visit number

prev_byvis <- function(var, data){
  newvar <- as.vector(unlist(subset(data, select = var)), mode = "numeric")
  data$newvar <- newvar
  newdat <- filter(data, visit_number != 3 & visit_number != 3.2 &
                     visit_number != 9 & visit_number != 9.2 & visit_number!= 9.3 &
                     visit_number != 9.9)
  table <- table(newdat$newvar, newdat$visit_number)
  props <- txtRound(100*prop.table(table, 2), digits = 1)
  prop_vec <- props[2, ]
  return(prop_vec)
}

summary_cat <- function(var, visit, data) {
  table <- table(var, 
                 filter(data, visit_number == visit)$visit_number)
  props <- txtRound(100*prop.table(table, 2), digits = 1)
  n <- table[2]
  summary <- cbind(n, props[2])
  return(summary)
}

summary_cat2 <- function(var, visit, data) {
  newvar <- as.vector(unlist(subset(data, select = var)), mode = "numeric")
  data$newvar <- newvar
  table <- table(filter(data, visit_number == visit)$newvar, 
                 filter(data, visit_number == visit)$visit_number)
  props <- txtRound(100*prop.table(table, 2), digits = 1)
  n <- table[2]
  summary <- cbind(n, props[2])
  return(summary)
}


## Function to get continuous variables for tables
summary_cont2 <- function(var, digits1, digits2) {
  mean <- txtRound(data.frame(mean(var, na.rm = TRUE)), digits = digits1)
  sd <- txtRound(data.frame(sd(var, na.rm = TRUE)), digits = digits2)
  summary <- paste(mean, sd, sep = " &plusmn; ")
  summary2 <- cbind(summary, "&ndash;")
  return(summary2)
}
## Function to get medians for skewed variables
summary_median <- function(var, digits1, digits2) {
  median <- txtRound(data.frame(median(var, na.rm = TRUE)), digits = digits1)
  iqr0 <- txtRound(data.frame(quantile(var, probs = c(0.25, 0.75), na.rm = TRUE)), digits = digits2)
  iqr <- paste(iqr0[1, ], iqr0[2, ], sep = ", ")
  iqr2 <- paste("(", iqr, sep = "")
  iqr3 <- paste(iqr2, ")", sep = "")
  summary <- paste(median, iqr3, sep = " ")
  summary2 <- cbind(summary, "&ndash;")
  return(summary2)
}


summary_median2 <- function(var, digits1, digits2) {
  median <- txtRound(data.frame(median(var, na.rm = TRUE)), digits = digits1)
  iqr0 <- txtRound(data.frame(quantile(var, probs = c(0.25, 0.75), na.rm = TRUE)), digits = digits2)
  iqr <- paste(iqr0[1, ], iqr0[2, ], sep = ", ")
  iqr2 <- paste("(", iqr, sep = "")
  iqr3 <- paste(iqr2, ")", sep = "")
  summary <- cbind(median, iqr3)
  colnames(summary) <- c("Median", "IQR (25%, 75%)")
  return(summary)
}

#Function to get prevalences
prevfun <- function(x){
  tab <- prop.table(table(x))
  percent <- 100*prop.table(tab)
  prev <- percent[2]
  return(prev)
}

# Function to get and save odds ratios
get_ORs <- function(model){
  CIs <- txtRound(as.data.frame(exp(confint(model))), digits = 2)
  smallCIs <- paste("(", CIs[,1], ",&nbsp;", CIs[,2], ")", sep = "")
  ORs <- txtRound(as.data.frame(exp(coef(model))), digits = 2)
  vars <- dim(ORs)[1]
  Walds <- (summary(model)$coefficients[ , 1])^2 / 
    (summary(model)$coefficients[ , 2])^2
  Pval <- txtPval((1 - pchisq(Walds, 1)), lim.2dec = 10^-2, 
                  lim.sig = 10^-4, html = TRUE)
  combo <- cbind(ORs, smallCIs, Pval)[2:vars, ]
  colnames(combo) <- c("OR", "95% CI", "P value")
  rownames(combo) <- rownames(confint(model))[2:vars]
  combo[] <- lapply(combo, as.character)
  return(combo)
}

# Create a function to save betas

get_betas <- function(model){
  vsub <- Vectorize(sub)
  CIs <- txtRound(as.data.frame(confint(model)), digits = 2)
  negCIs <- vsub("-", "&#8209;", CIs)
  smallCIs <- paste("(", negCIs[,1], ",&nbsp;",negCIs[,2], ")", sep = "")
  betas <- txtRound(as.data.frame(coef(model)), digits = 2)
  negbetas <- vsub("-", "&#x2011;", betas)
  vars <- dim(negbetas)[1]
  Pvals <- txtPval(coef(summary(model))[, "Pr(>|t|)"], lim.2dec = 10^-2, 
                   lim.sig = 10^-4, html = TRUE)
  combo <- cbind(negbetas, smallCIs, Pvals)[2:vars, ]
  colnames(combo) <- c("Beta", "95% CI", "P value")
  rownames(combo) <- rownames(confint(model))[2:vars]
  return(combo)
}

get_betas_transform <- function(model, whichbase){
  vsub <- Vectorize(sub)
  SEs <- as.data.frame(summary(model)$coefficients[ , 2])
  betas <- as.data.frame(summary(model)$coefficients[ , 1])
  expbetas <- txtRound(whichbase^betas, digits = 2)
  LCLs <- txtRound(whichbase^(betas - 1.96*SEs), digits = 2)
  UCLs <- txtRound(whichbase^(betas + 1.96*SEs), digits = 2)
  CIcomb <- cbind(LCLs, UCLs)
  CIs <- paste("(", CIcomb[,1], ",&nbsp;",CIcomb[,2], ")", sep = "")
  vars <- dim(betas)[1]
  Pvals <- txtPval(coef(summary(model))[, "Pr(>|t|)"], lim.2dec = 10^-2, 
                   lim.sig = 10^-4, html = TRUE)
  combo <- cbind(expbetas[2:vars, ], CIs[2:vars], Pvals[2:vars])
  colnames(combo) <- c("Ratio", "CI", "P value")
  rownames(combo) <- rownames(summary(model)$coefficients)[2:vars]
  return(combo)
}

get_betas_transform_percent <- function(model, whichbase){
  vsub <- Vectorize(sub)
  SEs <- as.data.frame(summary(model)$coefficients[ , 2])
  betas <- as.data.frame(summary(model)$coefficients[ , 1])
  expbetas <- txtRound(100*((whichbase^betas) - 1), digits = 1)
  LCLs <- txtRound(100*(whichbase^(betas - 1.96*SEs)-1), digits = 1)
  UCLs <- txtRound(100*(whichbase^(betas + 1.96*SEs) - 1), digits = 1)
  CIcomb <- cbind(LCLs, UCLs)
  CIs <- paste("(", CIcomb[,1], ",&nbsp;",CIcomb[,2], ")", sep = "")
  vars <- dim(betas)[1]
  Pvals <- txtPval(coef(summary(model))[, "Pr(>|t|)"], lim.2dec = 10^-2, 
                   lim.sig = 10^-4, html = TRUE)
  combo <- cbind(expbetas[2:vars, ], CIs[2:vars], Pvals[2:vars])
  colnames(combo) <- c("Percent Change", "CI", "P value")
  rownames(combo) <- rownames(summary(model)$coefficients)[2:vars]
  return(combo)
}



# Summary to get betas for random effects models
get_betas_random <- function(model){
  vsub <- Vectorize(sub)
  SEs <- txtRound(as.data.frame(summary(model)$coefficients[ , 2]), digits = 3)
  betas <- txtRound(as.data.frame(summary(model)$coefficients[ , 1]), digits = 3)
  negbetas <- vsub("-", "&#x2011;", betas)
  vars <- dim(negbetas)[1]
  ftests <- anova(model)
  Pvals <- txtPval(ftests[ , 6], lim.2dec = 10^-2, 
                   lim.sig = 10^-4, html = TRUE)
  combo <- cbind(negbetas[2:vars, ], SEs[2:vars, ], Pvals)
  colnames(combo) <- c("Beta", "SE", "P value")
  rownames(combo) <- rownames(summary(model)$coefficients)[2:vars]
  return(combo)
}

get_betas_random_transform <- function(model, whichbase){
  vsub <- Vectorize(sub)
  SEs <- as.data.frame(summary(model)$coefficients[ , 2])
  betas <- as.data.frame(summary(model)$coefficients[ , 1])
  expbetas <- txtRound(whichbase^betas, digits = 2)
  LCLs <- txtRound(whichbase^(betas - 1.96*SEs), digits = 2)
  UCLs <- txtRound(whichbase^(betas + 1.96*SEs), digits = 2)
  CIcomb <- cbind(LCLs, UCLs)
  CIs <- paste("(", CIcomb[,1], ",&nbsp;",CIcomb[,2], ")", sep = "")
  vars <- dim(betas)[1]
  ftests <- anova(model)
  Pvals <- txtPval(ftests[ , 6], lim.2dec = 10^-2, 
                   lim.sig = 10^-4, html = TRUE)
  combo <- cbind(expbetas[2:vars, ], CIs[2:vars], Pvals)
  colnames(combo) <- c("Ratio", "CI", "P value")
  rownames(combo) <- rownames(summary(model)$coefficients)[2:vars]
  return(combo)
}

## Function for ORs from GEEs
get_ORs_gee <- function(model){
  ORs1 <- as.data.frame(summary(model)$coefficients[ , 1])
  ORs <- txtRound(exp(ORs1), digits = 2)
  SEs <- as.data.frame(summary(model)$coefficients[ , 4])
  UCL <- txtRound(exp(ORs1 + (1.96*SEs)), digits = 2)
  LCL <- txtRound(exp(ORs1 - (1.96*SEs)), digits = 2)
  CIs <- cbind(LCL, UCL)
  smallCIs <- paste("(", CIs[ , 1], ",&nbsp;", CIs[ , 2], ")", sep = "")
  vars <- dim(ORs)[1]
  Walds <- (summary(model)$coefficients[ , 5])^2
  Pval <- txtPval((1 - pchisq(Walds, 1)), lim.2dec = 10^-2, 
                  lim.sig = 10^-4, html = TRUE)
  combo <- cbind(ORs, smallCIs, Pval)[2:vars, ]
  colnames(combo) <- c("OR", "95% CI", "P value")
  rownames(combo) <- rownames(ORs1)[2:vars]
  return(combo)
}

# ORs from GLMER
get_ORs_glmer <- function(model){
  ORs1 <- as.data.frame(summary(model)$coefficients[ , 1])
  ORs <- txtRound(exp(ORs1), digits = 2)
  SEs <- as.data.frame(summary(model)$coefficients[ , 2])
  UCL <- txtRound(exp(ORs1 + (1.96*SEs)), digits = 2)
  LCL <- txtRound(exp(ORs1 - (1.96*SEs)), digits = 2)
  CIs <- cbind(LCL, UCL)
  smallCIs <- paste("(", CIs[ , 1], ",&nbsp;", CIs[ , 2], ")", sep = "")
  vars <- dim(ORs)[1]
  Walds <- (summary(model)$coefficients[ , 1])^2 / 
    (summary(model)$coefficients[ , 2])^2
  Pval <- txtPval((1 - pchisq(Walds, 1)), lim.2dec = 10^-2, 
                  lim.sig = 10^-4, html = TRUE)
  combo <- cbind(ORs, smallCIs, Pval)[2:vars, ]
  colnames(combo) <- c("OR", "95% CI", "P value")
  rownames(combo) <- rownames(ORs1)[2:vars]
  return(combo)
}


# Generalized WALD for GEE
generalized_wald <- function(parameters, model){s
  VI <- solve(model$robust.variance[parameters, parameters])
  res <- as.numeric(coef(model)[parameters] %*% VI %*% coef(model)[parameters])
  pval <- txtPval(1 - pchisq(res, df = length(parameters)), lim.2dec = 10^-2,
                  lim.sig = 10^-4, html = TRUE)
  return(pval)
}

# Function to get ORs in data frame for model
get_OR_biv <- function(outcome, predictor, dat1){
  dat2 <- subset(dat1, select = c(outcome, predictor))
  model <- glm(dat2[ ,1] ~ dat2[ ,2], data = dat2)
  vars <- dim(as.data.frame(exp(coef(model))))[1]
  ortab <- get_ORs(model)
  rownames(ortab) <- rownames(confint(model))[2:vars]
  return(ortab)
}


# Function to get LM info to plot

lm_info = function(x, y, df){
  m <- lm(y ~ x, df)
  pval <- txtPval(coef(summary(m))[2, "Pr(>|t|)"], lim.2dec = 10^-2, 
                  lim.sig = 10^-4, html = FALSE)
  beta <- signif(coef(m)[2], digits = 2)
  r2 <- round(summary(m)$adj.r.squared, digits = 2)
  info <- paste("Beta = ", beta, "; p = ", pval, "; Adj. R2 = ", r2, sep = "")
  return(info)
}


## NIDI Functions


# Quantile function
dec1fn <- function(x){
  dec1 <- quantile(x, probs = (0.1), na.rm = TRUE)
  return(dec1)
}

newcorrfn <- function(x, visit, data, mom){
  # First just filter dataset to make smaller and more manageable
  subdat <- filter(data, visit_number == visit)
  if(mom == "b"){
    subdat2 <- subset(subdat, select = c("twinid", "visit_number", x, "crp_b",
                                         "agp_b", "inflcat_b", "ln_crp_b",
                                         "ln_agp_b"))
    subdat3 <- mutate(subdat2, ln_x = log(subdat2[ ,3]), x = subdat2[ , 3],
                      CRP = crp_b, AGP = agp_b, ln_CRP = ln_crp_b, 
                      ln_AGP = ln_agp_b, inflcat = inflcat_b)
    extdec <- c(-2.102639, -0.468499)
  }
  if(mom == "m"){
    subdat2 <- subset(subdat, select = c("twinid", "visit_number", x, "crp_m",
                                         "agp_m", "inflcat_m", "ln_crp_m",
                                         "ln_agp_m"))
    subdat3 <- mutate(subdat2, ln_x = log(subdat2[ ,3]), x = subdat2[ , 3],
                      CRP = crp_m, AGP = agp_m, ln_CRP = ln_crp_m, 
                      ln_AGP = ln_agp_m, inflcat = inflcat_m)
    extdec <- c(-1.75152, -0.625036)
  }
  grouped <- group_by(subdat3, inflcat) %>%
    filter(!is.na(inflcat))
  # Next grab CRP and AGP references
  ## CRP - by inflammation group
  crps <- summarise_each(grouped, funs(mean(., na.rm = TRUE), 
                                       median(., na.rm = TRUE), 
                                       dec1fn, 
                                       min(., na.rm = TRUE)), CRP)
  type <- c(rep("crp", dim(crps)[1]))
  crp <- cbind(type, crps)
  colnames(crp) <- c("var", "inflammation", "mean", "median", "decile.1", "min")
  crpln <- summarise_each(grouped, funs(mean(., na.rm = TRUE), 
                                        median(., na.rm = TRUE), 
                                        dec1fn, 
                                        min(., na.rm = TRUE)),
                          ln_CRP)
  crpln.2 <- cbind(crpln, extdec[1])
  type <- c(rep("crp", dim(crpln.2)[1]))
  ln.crp <- cbind(type, crpln.2)
  colnames(ln.crp) <- c("var", "inflammation", "mean", "median", "decile.1", "min", "ext.decile")
  ## CRP - overall
  crps.all <- summarise_each(subdat3, funs(mean(., na.rm = TRUE), 
                                           median(., na.rm = TRUE), 
                                           dec1fn, 
                                           min(., na.rm = TRUE)), CRP)
  type <- c(rep("crp_all", dim(crps.all)[1]))
  crp_all <- cbind(type, crps.all)
  colnames(crp_all) <- c("var", "mean", "median", "decile.1", "min")
  crpln_all <- summarise_each(subdat3, funs(mean(., na.rm = TRUE), 
                                            median(., na.rm = TRUE), 
                                            dec1fn, 
                                            min(., na.rm = TRUE)), ln_CRP)
  crpln.2_all <- cbind(crpln_all, extdec[1])
  type <- c(rep("crp_all", dim(crpln.2_all)[1]))
  ln.crp_all <- cbind(type, crpln.2_all)
  colnames(ln.crp_all) <- c("var", "mean", "median", "decile.1", "min", "ext.decile")
  
  ## AGP - by inflammation group
  agps <- summarise_each(grouped, funs(mean(., na.rm = TRUE), 
                                       median(., na.rm = TRUE), 
                                       dec1fn, 
                                       min(., na.rm = TRUE)), 
                         AGP)
  type <- c(rep("agp", dim(agps)[1]))
  agp <- cbind(type, agps)
  colnames(agp) <- c("var", "inflammation", "mean", "median", "decile.1", "min")
  agpln <- summarise_each(grouped, funs(mean(., na.rm = TRUE), 
                                        median(., na.rm = TRUE), 
                                        dec1fn, 
                                        min(., na.rm = TRUE)), 
                          ln_AGP)
  agpln.2 <- cbind(agpln, extdec[2])
  type <- c(rep("agp", dim(agpln.2)[1]))
  ln.agp <- cbind(type, agpln.2)
  colnames(ln.agp) <- c("var", "inflammation", "mean", "median", "decile.1", "min", "ext.decile")  
  infl.refs <- rbind(crp, agp)
  infl.refs.ln <- rbind(ln.crp, ln.agp)
  ## AGP - overall
  agps.all <- summarise_each(subdat3, funs(mean(., na.rm = TRUE), 
                                           median(., na.rm = TRUE), 
                                           dec1fn, 
                                           min(., na.rm = TRUE)), AGP)
  type <- c(rep("agp_all", dim(agps.all)[1]))
  agp_all <- cbind(type, agps.all)
  colnames(agp_all) <- c("var", "mean", "median", "decile.1", "min")
  agpln_all <- summarise_each(subdat3, funs(mean(., na.rm = TRUE), 
                                            median(., na.rm = TRUE), 
                                            dec1fn, 
                                            min(., na.rm = TRUE)), ln_AGP)
  agpln.2_all <- cbind(agpln_all, extdec[1])
  type <- c(rep("agp_all", dim(agpln.2_all)[1]))
  ln.agp_all <- cbind(type, agpln.2_all)
  colnames(ln.agp_all) <- c("var",  "mean", "median", "decile.1", "min", "ext.decile")
  
  
  # Next get means for saving and making CF
  ## Get means and save
  vals <- summarise_each(grouped, funs(mean(., na.rm = TRUE), 
                                       median(., na.rm = TRUE), 
                                       dec1fn), x)
  colnames(vals) <- c("inflammation", "mean", "median", "decile.1")
  cf.ref0 <- cbind(vals[1, 2],vals[1, 3], vals[1, 4])
  median.ratios <- mutate(vals, cf.median.0 = median / cf.ref0[,2] )
  ln.vals <- summarise_each(grouped, funs(mean(., na.rm = TRUE),
                                          median(., na.rm = TRUE),
                                          dec1fn), ln_x)
  colnames(ln.vals) <- c("inflammation", "mean", "median", "decile.1")
  ### Create CF
  cf.refs <- cbind(ln.vals[1, 2], ln.vals[1, 3], ln.vals[1, 4])
  cf <- mutate(ln.vals, cf.mean = mean / cf.refs[,1],
               cf.median = median / cf.refs[ ,2],
               cf.decile.1 = decile.1 / cf.refs[ ,3],
               cf2.mean = mean - cf.refs[,1])
  cf$cf.median.0 <- median.ratios$cf.median.0
  cf <- mutate(cf, cf2.mean.exp = exp(cf2.mean))
  # Difference of the means of the lns is the ratio of the geometric means
  
  # Run linear regressions and save values
  linreg <- lm(ln_x ~ ln_CRP + ln_AGP, subdat3)
  coef <- coef(linreg)
  linreg.int <- lm(ln_x ~ ln_CRP + ln_AGP + ln_CRP*ln_AGP, subdat3)
  coef.int <- coef(linreg.int)
  
  # Do corrections
  ## Thurnham correction
  subdat3$x_th <- ifelse(subdat3$inflcat == 0, subdat3$x, 
                         ifelse(subdat3$inflcat == 1, subdat3$x / cf[2,9],
                                ifelse(subdat3$inflcat == 2, subdat3$x / cf[3,9],
                                       ifelse(subdat3$inflcat == 3, 
                                              subdat3$x / cf[4,9],
                                              NA))))
  subdat3$x_th <- unlist(subdat3$x_th)
  subdat3$x_th1 <- ifelse(subdat3$inflcat == 0, subdat3$x, 
                          ifelse(subdat3$inflcat == 1, subdat3$x / cf[2,10],
                                 ifelse(subdat3$inflcat == 2, subdat3$x / cf[3,10],
                                        ifelse(subdat3$inflcat == 3, 
                                               subdat3$x / cf[4,10],
                                               NA))))
  subdat3$x_th1 <- unlist(subdat3$x_th1)
  
  ## Linear regressions
  ### Linear regression correction with mean reference - no-interaction
  subdat3 <- mutate(subdat3, 
                    x_lc.1 = ifelse(ln_CRP > ln.crp[1,3] & ln_AGP > ln.agp[1,3],
                                    exp(ln_x - coef[2]*(ln_CRP - ln.crp[1,3]) -
                                          coef[3]*(ln_AGP - ln.agp[1,3])),
                                    ifelse(ln_CRP > ln.crp[1,3] & ln_AGP <= ln.agp[1,3],
                                           exp(ln_x - coef[2]*(ln_CRP - ln.crp[1,3])),
                                           ifelse(ln_CRP <= ln.crp[1,3] & 
                                                    ln_AGP > ln.agp[1,3],
                                                  exp(ln_x - coef[3]*(ln_AGP - 
                                                                        ln.agp[1,3])),
                                                  exp(ln_x)))))
  ### Linear regression correction with mean reference - with interaction
  subdat3 <- mutate(subdat3, 
                    x_lc.1.1 = ifelse(ln_CRP > ln.crp[1,3] & ln_AGP > ln.agp[1,3],
                                      exp(ln_x - coef.int[2]*(ln_CRP - ln.crp[1,3]) -
                                            coef.int[3]*(ln_AGP - ln.agp[1,3]) -
                                            coef.int[4]*(ln_AGP - ln.agp[1,3])*
                                            (ln_CRP - ln.crp[1, 3])),
                                      ifelse(ln_CRP <= ln.crp[1,3] & 
                                               ln_AGP > ln.agp[1,3], 
                                             exp(ln_x - 
                                                   coef.int[3]*(ln_AGP - ln.agp[1,3]) -
                                                   coef.int[4]*(ln_AGP - ln.agp[1,3])*
                                                   (ln_CRP)),
                                             ifelse(ln_CRP > ln.crp[1,3] & 
                                                      ln_AGP <= ln.agp[1,3],
                                                    exp(ln_x - coef.int[2]*(ln_CRP - 
                                                                              ln.crp[1,3]) -
                                                          coef.int[4]*(ln_AGP)*
                                                          (ln_CRP - ln.crp[1, 3])),
                                                    exp(ln_x - 
                                                          coef.int[4]*(ln_AGP)*(ln_CRP))))))
  
  ### Linear regression correction with median reference - no interaction
  subdat3 <- mutate(subdat3, 
                    x_lc.2 = ifelse(ln_CRP > ln.crp[1,4] & ln_AGP > ln.agp[1,4],
                                    exp(ln_x - coef[2]*(ln_CRP - ln.crp[1,4]) -
                                          coef[3]*(ln_AGP - ln.agp[1,4])),
                                    ifelse(ln_CRP > ln.crp[1,4] & ln_AGP <= ln.agp[1,4],
                                           exp(ln_x - coef[2]*(ln_CRP - ln.crp[1,4])),
                                           ifelse(ln_CRP <= ln.crp[1,4] & 
                                                    ln_AGP > ln.agp[1,4],
                                                  exp(ln_x - coef[3]*(ln_AGP - 
                                                                        ln.agp[1,4])),
                                                  exp(ln_x)))))
  ### Linear regression correction with median reference - with interaction
  subdat3 <- mutate(subdat3, 
                    x_lc.2.1 = ifelse(ln_CRP > ln.crp[1,4] & ln_AGP > ln.agp[1,4],
                                      exp(ln_x - coef.int[2]*(ln_CRP - ln.crp[1,4]) -
                                            coef.int[3]*(ln_AGP - ln.agp[1,4]) -
                                            coef.int[4]*(ln_AGP - ln.agp[1,4])*
                                            (ln_CRP - ln.crp[1,4])),
                                      ifelse(ln_CRP <= ln.crp[1,4] & 
                                               ln_AGP > ln.agp[1,4], 
                                             exp(ln_x - 
                                                   coef.int[3]*(ln_AGP - ln.agp[1,4]) -
                                                   coef.int[4]*(ln_AGP - ln.agp[1,4])*
                                                   (ln_CRP)),
                                             ifelse(ln_CRP > ln.crp[1,4] & 
                                                      ln_AGP <= ln.agp[1,4],
                                                    exp(ln_x - coef.int[2]*(ln_CRP - 
                                                                              ln.crp[1,4]) -
                                                          coef.int[4]*(ln_AGP)*
                                                          (ln_CRP - ln.crp[1,4])),
                                                    exp(ln_x - 
                                                          coef.int[4]*(ln_AGP)*(ln_CRP))))))
  
  
  ### Linear regression correction with decile reference - no interaction
  subdat3 <- mutate(subdat3, 
                    x_lc.3 = ifelse(ln_CRP > ln.crp[1,5] & ln_AGP > ln.agp[1,5],
                                    exp(ln_x - coef[2]*(ln_CRP - ln.crp[1,5]) -
                                          coef[3]*(ln_AGP - ln.agp[1,5])),
                                    ifelse(ln_CRP > ln.crp[1,5] & ln_AGP <= ln.agp[1,5],
                                           exp(ln_x - coef[2]*(ln_CRP - ln.crp[1,5])),
                                           ifelse(ln_CRP <= ln.crp[1,5] & 
                                                    ln_AGP > ln.agp[1,5],
                                                  exp(ln_x - coef[3]*(ln_AGP - 
                                                                        ln.agp[1,5])),
                                                  exp(ln_x)))))
  
  ### Linear regression correction with decile reference - CORRECT ALL
  subdat3 <- mutate(subdat3, 
                    x_lc.3_allcorr = exp(ln_x - coef[2]*(ln_CRP - ln.crp[1,5]) -
                                           coef[3]*(ln_AGP - ln.agp[1,5])))
  
  ### Linear regression correction with decile reference - DECILE OF ALL vs. UNINFLAMED
  subdat3 <- mutate(subdat3, 
                    x_lc.3.5 = ifelse(ln_CRP > ln.crp_all[1,5] & ln_AGP > ln.agp_all[1,5],
                                      exp(ln_x - coef[2]*(ln_CRP - ln.crp_all[1,5]) -
                                            coef[3]*(ln_AGP - ln.agp_all[1,5])),
                                      ifelse(ln_CRP > ln.crp_all[1,5] & ln_AGP <= ln.agp_all[1,5],
                                             exp(ln_x - coef[2]*(ln_CRP - ln.crp_all[1,5])),
                                             ifelse(ln_CRP <= ln.crp_all[1,5] & 
                                                      ln_AGP > ln.agp_all[1,5],
                                                    exp(ln_x - coef[3]*(ln_AGP - 
                                                                          ln.agp_all[1,5])),
                                                    exp(ln_x)))))
  
  ### Linear regression correction with decile reference - with interaction
  subdat3 <- mutate(subdat3, 
                    x_lc.3.1 = ifelse(ln_CRP > ln.crp[1,5] & ln_AGP > ln.agp[1,5],
                                      exp(ln_x - coef.int[2]*(ln_CRP - ln.crp[1,5]) -
                                            coef.int[3]*(ln_AGP - ln.agp[1,5]) -
                                            coef.int[4]*(ln_AGP - ln.agp[1,5])*
                                            (ln_CRP - ln.crp[1,5])),
                                      ifelse(ln_CRP <= ln.crp[1,5] & 
                                               ln_AGP > ln.agp[1,5], 
                                             exp(ln_x - 
                                                   coef.int[3]*(ln_AGP - ln.agp[1,5]) -
                                                   coef.int[4]*(ln_AGP - ln.agp[1,5])*
                                                   (ln_CRP)),
                                             ifelse(ln_CRP > ln.crp[1,5] & 
                                                      ln_AGP <= ln.agp[1,5],
                                                    exp(ln_x - coef.int[2]*(ln_CRP - 
                                                                              ln.crp[1,5]) -
                                                          coef.int[4]*(ln_AGP)*
                                                          (ln_CRP - ln.crp[1,5])),
                                                    exp(ln_x - 
                                                          coef.int[4]*(ln_AGP)*(ln_CRP))))))
  
  ### Linear regression correction with minimum reference - no interaction
  subdat3 <- mutate(subdat3, 
                    x_lc.4 = ifelse(ln_CRP > ln.crp[1,6] & ln_AGP > ln.agp[1,6],
                                    exp(ln_x - coef[2]*(ln_CRP - ln.crp[1,6]) -
                                          coef[3]*(ln_AGP - ln.agp[1,6])),
                                    ifelse(ln_CRP > ln.crp[1,6] & ln_AGP <= ln.agp[1,6],
                                           exp(ln_x - coef[2]*(ln_CRP - ln.crp[1,6])),
                                           ifelse(ln_CRP <= ln.crp[1,6] & 
                                                    ln_AGP > ln.agp[1,6],
                                                  exp(ln_x - coef[3]*(ln_AGP - 
                                                                        ln.agp[1,6])),
                                                  exp(ln_x)))))
  
  ### Linear regression correction with minimum reference - with interaction
  subdat3 <- mutate(subdat3, 
                    x_lc.4.1 = ifelse(ln_CRP > ln.crp[1,6] & ln_AGP > ln.agp[1,6],
                                      exp(ln_x - coef.int[2]*(ln_CRP - ln.crp[1,6]) -
                                            coef.int[3]*(ln_AGP - ln.agp[1,6]) -
                                            coef.int[4]*(ln_AGP - ln.agp[1,6])*
                                            (ln_CRP - ln.crp[1,6])),
                                      ifelse(ln_CRP <= ln.crp[1,6] & 
                                               ln_AGP > ln.agp[1,6], 
                                             exp(ln_x - 
                                                   coef.int[3]*(ln_AGP - ln.agp[1,6]) -
                                                   coef.int[4]*(ln_AGP - ln.agp[1,6])*
                                                   (ln_CRP)),
                                             ifelse(ln_CRP > ln.crp[1,6] & 
                                                      ln_AGP <= ln.agp[1,6],
                                                    exp(ln_x - coef.int[2]*(ln_CRP - 
                                                                              ln.crp[1,6]) -
                                                          coef.int[4]*(ln_AGP)*
                                                          (ln_CRP - ln.crp[1,6])),
                                                    exp(ln_x - 
                                                          coef.int[4]*(ln_AGP)*(ln_CRP))))))
  
  ### Linear regression correction with external decile reference - no interaction
  subdat3 <- mutate(subdat3, 
                    x_lc.5 = ifelse(ln_CRP > ln.crp[1,7] & ln_AGP > ln.agp[1,7],
                                    exp(ln_x - coef[2]*(ln_CRP - ln.crp[1,7]) -
                                          coef[3]*(ln_AGP - ln.agp[1,7])),
                                    ifelse(ln_CRP > ln.crp[1,7] & ln_AGP <= ln.agp[1,7],
                                           exp(ln_x - coef[2]*(ln_CRP - ln.crp[1,7])),
                                           ifelse(ln_CRP <= ln.crp[1,7] & 
                                                    ln_AGP > ln.agp[1,7],
                                                  exp(ln_x - coef[3]*(ln_AGP - 
                                                                        ln.agp[1,7])),
                                                  exp(ln_x)))))
  
  ### Linear regression correction with external decile reference - with interaction
  subdat3 <- mutate(subdat3, 
                    x_lc.5.1 = ifelse(ln_CRP > ln.crp[1,7] & ln_AGP > ln.agp[1,7],
                                      exp(ln_x - coef.int[2]*(ln_CRP - ln.crp[1,7]) -
                                            coef.int[3]*(ln_AGP - ln.agp[1,7]) -
                                            coef.int[4]*(ln_AGP - ln.agp[1,7])*
                                            (ln_CRP - ln.crp[1,7])),
                                      ifelse(ln_CRP <= ln.crp[1,7] & 
                                               ln_AGP > ln.agp[1,7], 
                                             exp(ln_x - 
                                                   coef.int[3]*(ln_AGP - ln.agp[1,7]) -
                                                   coef.int[4]*(ln_AGP - ln.agp[1,7])*
                                                   (ln_CRP)),
                                             ifelse(ln_CRP > ln.crp[1,7] & 
                                                      ln_AGP <= ln.agp[1,7],
                                                    exp(ln_x - coef.int[2]*(ln_CRP - 
                                                                              ln.crp[1,7]) -
                                                          coef.int[4]*(ln_AGP)*
                                                          (ln_CRP - ln.crp[1,7])),
                                                    exp(ln_x - 
                                                          coef.int[4]*(ln_AGP)*(ln_CRP))))))
  
  return(subdat3)
  
}

# Function to rename interim variables 
renamefn <- function(data, var, mom){
  if(var == "ferritin" & mom == "m"){
    newdat <- rename(data, ln_fer_m = ln_x, fer_m_th.1 = x_th, 
                     fer_m_th.1.1 = x_th1,
                     fer_m_lc.1 = x_lc.1, fer_m_lc1.1 = x_lc.1.1,
                     fer_m_lc.2 = x_lc.2, fer_m_lc2.1 = x_lc.2.1,
                     fer_m_lc.3 = x_lc.3, fer_m_lc3.1 = x_lc.3.1,
                     fer_m_lc.4 = x_lc.4, fer_m_lc4.1 = x_lc.4.1,
                     fer_m_lc.5 = x_lc.5, fer_m_lc5.1 = x_lc.5.1,
                     fer_m_lc.3_allcorr = x_lc.3_allcorr,
                     fer_m_lc.3.5 = x_lc.3.5)}
  if(var == "ferritin" & mom == "b"){
    newdat <- rename(data, ln_fer_b = ln_x, fer_b_th.1 = x_th,
                     fer_b_th.1.1 = x_th1,
                     fer_b_lc.1 = x_lc.1, fer_b_lc1.1 = x_lc.1.1,
                     fer_b_lc.2 = x_lc.2, fer_b_lc2.1 = x_lc.2.1,
                     fer_b_lc.3 = x_lc.3, fer_b_lc3.1 = x_lc.3.1,
                     fer_b_lc.4 = x_lc.4, fer_b_lc4.1 = x_lc.4.1,
                     fer_b_lc.5 = x_lc.5, fer_b_lc5.1 = x_lc.5.1,
                     fer_b_lc.3_allcorr = x_lc.3_allcorr,
                     fer_b_lc.3.5 = x_lc.3.5)}
  if(var == "stfr" & mom == "m"){
    newdat <- rename(data, ln_stfr_m = ln_x, stfr_m_th.1 = x_th,
                     stfr_m_th.1.1 = x_th1,
                     stfr_m_lc.1 = x_lc.1, stfr_m_lc1.1 = x_lc.1.1,
                     stfr_m_lc.2 = x_lc.2, stfr_m_lc2.1 = x_lc.2.1,
                     stfr_m_lc.3 = x_lc.3, stfr_m_lc3.1 = x_lc.3.1,
                     stfr_m_lc.4 = x_lc.4, stfr_m_lc4.1 = x_lc.4.1,
                     stfr_m_lc.5 = x_lc.5, stfr_m_lc5.1 = x_lc.5.1,
                     stfr_m_lc.3_allcorr = x_lc.3_allcorr,
                     stfr_m_lc.3.5 = x_lc.3.5)}
  if(var == "stfr" & mom == "b"){
    newdat <- rename(data, ln_stfr_b = ln_x, stfr_b_th.1 = x_th,
                     stfr_b_th.1.1 = x_th1,
                     stfr_b_lc.1 = x_lc.1, stfr_b_lc1.1 = x_lc.1.1,
                     stfr_b_lc.2 = x_lc.2, stfr_b_lc2.1 = x_lc.2.1,
                     stfr_b_lc.3 = x_lc.3, stfr_b_lc3.1 = x_lc.3.1,
                     stfr_b_lc.4 = x_lc.4, stfr_b_lc4.1 = x_lc.4.1,
                     stfr_b_lc.5 = x_lc.5, stfr_b_lc5.1 = x_lc.5.1,
                     stfr_b_lc.3_allcorr = x_lc.3_allcorr,
                     stfr_b_lc.3.5 = x_lc.3.5)}
  if(var == "rbp" & mom == "m"){
    newdat <- rename(data, ln_rbp_m = ln_x, rbp_m_th.1 = x_th,
                     rbp_m_th.1.1 = x_th1,
                     rbp_m_lc.1 = x_lc.1, rbp_m_lc1.1 = x_lc.1.1,
                     rbp_m_lc.2 = x_lc.2, rbp_m_lc2.1 = x_lc.2.1,
                     rbp_m_lc.3 = x_lc.3, rbp_m_lc3.1 = x_lc.3.1,
                     rbp_m_lc.4 = x_lc.4, rbp_m_lc4.1 = x_lc.4.1,
                     rbp_m_lc.5 = x_lc.5, rbp_m_lc5.1 = x_lc.5.1,
                     rbp_m_lc.3_allcorr = x_lc.3_allcorr,
                     rbp_m_lc.3.5 = x_lc.3.5)}
  if(var == "rbp" & mom == "b"){
    newdat <- rename(data, ln_rbp_b = ln_x, rbp_b_th.1 = x_th,
                     rbp_b_th.1.1 = x_th1,
                     rbp_b_lc.1 = x_lc.1, rbp_b_lc1.1 = x_lc.1.1,
                     rbp_b_lc.2 = x_lc.2, rbp_b_lc2.1 = x_lc.2.1,
                     rbp_b_lc.3 = x_lc.3, rbp_b_lc3.1 = x_lc.3.1,
                     rbp_b_lc.4 = x_lc.4, rbp_b_lc4.1 = x_lc.4.1,
                     rbp_b_lc.5 = x_lc.5, rbp_b_lc5.1 = x_lc.5.1,
                     rbp_b_lc.3_allcorr = x_lc.3_allcorr,
                     rbp_b_lc.3.5 = x_lc.3.5)}
  if(var == "hb" & mom == "m"){
    newdat <- rename(data, ln_hb_m = ln_x, hb_m_th.1 = x_th,
                     hb_m_th.1.1 = x_th1,
                     hb_m_lc.1 = x_lc.1, hb_m_lc1.1 = x_lc.1.1,
                     hb_m_lc.2 = x_lc.2, hb_m_lc2.1 = x_lc.2.1,
                     hb_m_lc.3 = x_lc.3, hb_m_lc3.1 = x_lc.3.1,
                     hb_m_lc.4 = x_lc.4, hb_m_lc4.1 = x_lc.4.1,
                     hb_m_lc.5 = x_lc.5, hb_m_lc5.1 = x_lc.5.1,
                     hb_m_lc.3_allcorr = x_lc.3_allcorr,
                     hb_m_lc.3.5 = x_lc.3.5)}
  if(var == "hb" & mom == "b"){
    newdat <- rename(data, ln_hb_b = ln_x, hb_b_th.1 = x_th,
                     hb_b_th.1.1 = x_th1,
                     hb_b_lc.1 = x_lc.1, hb_b_lc1.1 = x_lc.1.1,
                     hb_b_lc.2 = x_lc.2, hb_b_lc2.1 = x_lc.2.1,
                     hb_b_lc.3 = x_lc.3, hb_b_lc3.1 = x_lc.3.1,
                     hb_b_lc.4 = x_lc.4, hb_b_lc4.1 = x_lc.4.1,
                     hb_b_lc.5 = x_lc.5, hb_b_lc5.1 = x_lc.5.1,
                     hb_b_lc.3_allcorr = x_lc.3_allcorr,
                     hb_b_lc.3.5 = x_lc.3.5)}
  if(var == "zn" & mom == "m"){
    newdat <- rename(data, ln_zn_m = ln_x, zn_m_th.1 = x_th,
                     zn_m_th.1.1 = x_th1,
                     zn_m_lc.1 = x_lc.1, zn_m_lc1.1 = x_lc.1.1,
                     zn_m_lc.2 = x_lc.2, zn_m_lc2.1 = x_lc.2.1,
                     zn_m_lc.3 = x_lc.3, zn_m_lc3.1 = x_lc.3.1,
                     zn_m_lc.4 = x_lc.4, zn_m_lc4.1 = x_lc.4.1,
                     zn_m_lc.5 = x_lc.5, zn_m_lc5.1 = x_lc.5.1,
                     zn_m_lc.3_allcorr = x_lc.3_allcorr,
                     zn_m_lc.3.5 = x_lc.3.5)}
  if(var == "zn" & mom == "b"){
    newdat <- rename(data, ln_zn_b = ln_x, zn_b_th.1 = x_th,
                     zn_b_th.1.1 = x_th1,
                     zn_b_lc.1 = x_lc.1, zn_b_lc1.1 = x_lc.1.1,
                     zn_b_lc.2 = x_lc.2, zn_b_lc2.1 = x_lc.2.1,
                     zn_b_lc.3 = x_lc.3, zn_b_lc3.1 = x_lc.3.1,
                     zn_b_lc.4 = x_lc.4, zn_b_lc4.1 = x_lc.4.1,
                     zn_b_lc.5 = x_lc.5, zn_b_lc5.1 = x_lc.5.1,
                     zn_b_lc.3_allcorr = x_lc.3_allcorr,
                     zn_b_lc.3.5 = x_lc.3.5)}
  if(var == "retinol" & mom == "m"){
    newdat <- rename(data, ln_ret_m = ln_x, ret_m_th.1 = x_th,
                     ret_m_th.1.1 = x_th1,
                     ret_m_lc.1 = x_lc.1, ret_m_lc1.1 = x_lc.1.1,
                     ret_m_lc.2 = x_lc.2, ret_m_lc2.1 = x_lc.2.1,
                     ret_m_lc.3 = x_lc.3, ret_m_lc3.1 = x_lc.3.1,
                     ret_m_lc.4 = x_lc.4, ret_m_lc4.1 = x_lc.4.1,
                     ret_m_lc.5 = x_lc.5, ret_m_lc5.1 = x_lc.5.1,
                     ret_m_lc.3_allcorr = x_lc.3_allcorr,
                     ret_m_lc.3.5 = x_lc.3.5)}
  if(var == "retinol" & mom == "b"){
    newdat <- rename(data, ln_ret_b = ln_x, ret_b_th.1 = x_th,
                     ret_b_th.1.1 = x_th1,
                     ret_b_lc.1 = x_lc.1, ret_b_lc1.1 = x_lc.1.1,
                     ret_b_lc.2 = x_lc.2, ret_b_lc2.1 = x_lc.2.1,
                     ret_b_lc.3 = x_lc.3, ret_b_lc3.1 = x_lc.3.1,
                     ret_b_lc.4 = x_lc.4, ret_b_lc4.1 = x_lc.4.1,
                     ret_b_lc.5 = x_lc.5, ret_b_lc5.1 = x_lc.5.1,
                     ret_b_lc.3_allcorr = x_lc.3_allcorr,
                     ret_b_lc.3.5 = x_lc.3.5)}
  
  newdat <- mutate(newdat, x = NULL, CRP = NULL, AGP = NULL, ln_CRP = NULL,
                   ln_AGP = NULL, inflcat = NULL)
  return(newdat)
}




## NEW FUNCTION FROM R COOKBOOK
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }}}
