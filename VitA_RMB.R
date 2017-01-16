# RMB Vitamin A

setwd("C:/Users/LXX8/Documents/Dissertation Followup/vitamin a")


# Get necessary packages
library(dplyr)
library(ggplot2)
library(htmlTable)
library(reshape2)

# Read data
small <- read.csv("NIDI_streamlined.csv")


### FIRST EXPLORE DEAD INFANTS AND ANY RELATIONSHIP
# Figure out which infants died and create a marker for them

table(small$b_dead)

deadones <- subset(small, b_dead == 1, select = c("twinid", "visit_number", "mellizo_trillizo", "male"))

# Add variable for dead ones

small <- mutate(small, dead_baby = ifelse(twinid %in% c("134-B", "183-B", "213-B", "235-B", "764-G"), 1, 0))

table(small$dead_baby)

deadones2 <- subset(small, dead_baby == 1, select = c("twinid", "visit_number", "mellizo_trillizo", "male", "rbp_b", "vad_b_un",
                                                      "vad_rbp_b_lc.3.5"))

baby764g <- subset(small, twinid == "764-G", select = c("twinid", "visit_number", "mellizo_trillizo", "male"))

vis2 <- filter(small, visit_number == 2)
tab1 <- table(vis2$vad_b_un, vis2$dead_baby)
tab2 <- table(vis2$vad_rbp_b_lc.3.5, vis2$dead_baby)

fisher.test(tab2)

ggplot(aes(y = rbp_b, x = agemo2, colour = dead_baby), data = small) + geom_point()

# Only three infants died after the first measurement; all were VAD but hard to assess

small <- mutate(small, visfact = factor(visit_number))

smaller <- filter(small, visit_number == 1 | visit_number == 2 | visit_number == 6 | visit_number == 8)
  
### NOW LETs PLOT AGP and CRP vs RETINOL
## CRP and retinol - neg rel for 6,8 mo, pos for 2 mo
vis1 <- filter(smaller, visit_number == 1)
vis2 <- filter(smaller, visit_number == 2)
vis6 <- filter(smaller, visit_number == 6)
vis8 <- filter(smaller, visit_number == 8)


## Moms
ggplot(data= smaller, aes(x = ln_crp_m, y = retinol_umolL_m, colour = visfact)) + geom_point(size  = 2) +
  scale_color_manual(values = c("red", "blue", "green", "orange")) + geom_smooth(method = "lm") + 
  annotate("text",  label = paste("1 month postpartum: ", lm_info(y = vis1$retinol_umolL_m, x = vis1$ln_crp_m, 
                                    df = vis1)), x = 2.5, y = 3) +
  annotate("text",  label = paste("6 months postpartum: ", lm_info(y = vis6$retinol_umolL_m, x = vis6$ln_crp_m, 
                                    df = vis6)),  x = 2.5, y = 2.75) +
  ggtitle("Retinol by CRP by Visit, mothers")


ggplot(data= smaller, aes(x = ln_agp_m, y = retinol_umolL_m, colour = visfact)) + geom_point(size  = 2) +
  scale_color_manual(values = c("red", "blue", "green", "orange")) + geom_smooth(method = "lm") + 
  annotate("text",  label = paste("1 month postpartum: ", lm_info(y = vis1$retinol_umolL_m, x = vis1$ln_agp_m, 
                                    df = vis1)),           x = 0.5, y = 2.75) +
  annotate("text",  label = paste("6 months postpartum: ", lm_info(y = vis6$retinol_umolL_m, x = vis6$ln_agp_m, 
                                    df = vis6)),            x = 0.5, y = 2.5) +
  ggtitle("Retinol by AGP by Visit, mothers")

ggplot(data= smaller, aes(x = ln_crp_m, y = rbp_m, colour = visfact)) + geom_point(size  = 2) +
  scale_color_manual(values = c("red", "blue", "green", "orange")) + geom_smooth(method = "lm") + 
  annotate("text",  label = paste("1 month postpartum: ", lm_info(y = vis1$rbp_m, x = vis1$ln_crp_m, 
                                    df = vis1)), 
           x =3, y = 4) +
  annotate("text",  label = paste("6 months postpartum: ", lm_info(y = vis6$rbp_m, x = vis6$ln_crp_m, 
                                    df = vis6)), 
           x = 3, y = 3.5) +
  ggtitle("RBP by CRP by Visit, mothers")


ggplot(data= smaller, aes(x = ln_agp_m, y = rbp_m, colour = visfact)) + geom_point(size  = 2) +
  scale_color_manual(values = c("red", "blue", "green", "orange")) + geom_smooth(method = "lm") + 
  annotate("text",  label = paste("1 month postpartum: ", lm_info(y = vis1$rbp_m, x = vis1$ln_agp_m, 
                                    df = vis1)),      x = 0.75, y = 4) +
  annotate("text",  label = paste("6 months postpartum: ", lm_info(y = vis6$rbp_m, x = vis6$ln_agp_m, 
                                    df = vis6)),          x = 0.75, y = 3.5) +
  ggtitle("RBP by AGP by Visit, mothers")

## infants



ggplot(data= smaller, aes(x = ln_crp_b, y = retinol_umolL_b, colour = visfact)) + geom_point(size  = 2) +
  scale_color_manual(values = c("red", "green", "blue", "orange")) + geom_smooth(method = "lm") + 
  annotate("text",  label = paste("2 mo: ", lm_info(y = vis2$retinol_umolL_b, x = vis2$ln_crp_b, 
                                    df = vis2)), 
                                    x = 2.5, y = 2) +
  annotate("text",  label = paste("6 mo:"  , lm_info(y = vis6$retinol_umolL_b, x = vis6$ln_crp_b, 
                                    df = vis6)), 
                                    x = 2.5, y = 1.9) +
  annotate("text",  label = paste("12-18 mo:"  , lm_info(y = filter(smaller, visit_number == 8)$retinol_umolL_b, 
                                    x = filter(smaller, visit_number == 8)$ln_crp_b, 
                                    df = filter(smaller, visit_number == 8))), 
                                    x = 2.5, y = 1.8 )
  
             


ggplot(data= smaller, aes(x = ln_crp_b, y = retinol_umolL_b)) + geom_point(size  = 2) +
   geom_smooth(method = "lm") +
  annotate("text",  label = lm_info(y =smaller$retinol_umolL_b, 
                                    x = smaller$ln_crp_b, 
                                    df = smaller), 
           x = 2.5, y = 2)




## CRP and RBP -- neg for 6,8, null for 2mo
ggplot(data= smaller, aes(x = ln_crp_b, y = rbp_b, colour = visfact)) + geom_point(size  = 2) +
  scale_color_manual(values = c("orange", "red", "green", "blue")) + geom_smooth(method = "lm") + 
  annotate("text",  label = paste("2 mo: ", lm_info(y = vis2$rbp_b, x = vis2$ln_crp_b, 
                                    df = vis2)), 
           x = 2.5, y = 2) +
  annotate("text",  label = paste("6 mo: ", lm_info(y = vis6$rbp_b, x = vis6$ln_crp_b, 
                                    df = vis6)), 
           x = 2.5, y = 1.9) +
  annotate("text",  label = paste("12-18 mo: ", lm_info(y = filter(smaller, visit_number == 8)$rbp_b, 
                                    x = filter(smaller, visit_number == 8)$ln_crp_b, 
                                    df = filter(smaller, visit_number == 8))), 
           x = 2.5, y = 1.8)



## AGP and reinol -- neg for 6, 8, mull for 2
ggplot(data= smaller, aes(x = ln_agp_b, y = retinol_umolL_b, colour = visfact)) + geom_point(size  = 2) +
  scale_color_manual(values = c("orange", "red", "green", "blue")) + geom_smooth(method = "lm") + 
  annotate("text",  label = paste("2 mo: ", lm_info(y = vis2$retinol_umolL_b, x = vis2$ln_agp_b, 
                                    df = vis2)), 
           x = 0.5, y = 2) +
  annotate("text",  label =paste("6 mo: ",  lm_info(y = vis6$retinol_umolL_b, x = vis6$ln_agp_b, 
                                    df = vis6)), 
           x = 0.5, y = 1.8) +
  annotate("text",  label = paste("12 - 18 mo: ", 
                                  lm_info(y = filter(smaller, visit_number == 8)$retinol_umolL_b, 
                                    x = filter(smaller, visit_number == 8)$ln_agp_b, 
                                    df = filter(smaller, visit_number == 8))), 
           x = 0.5, y = 1.6 )



ggplot(data= smaller, aes(x = ln_agp_b, y = retinol_umolL_b)) + geom_point(size  = 2) +
 geom_smooth(method = "lm") +
  annotate("text",  label = lm_info(y = smaller$retinol_umolL_b, x = smaller$ln_agp_b, 
                                    df = smaller), 
           x = -2, y = 0.5) 
  
  
  
   

ggplot(data= smaller, aes(x = ln_agp_b, y = rbp_b, colour = visfact)) + geom_point(size  = 2) +
  scale_color_manual(values = c("orange", "red", "green", "blue")) + geom_smooth(method = "lm") +
  annotate("text",  label = paste("2 mo: ", lm_info(y = vis2$rbp_b, x = vis2$ln_agp_b, 
                                    df = vis2)), 
           x = 0.5, y = 2) +
  annotate("text",  label = paste("6 mo: ", lm_info(y = vis6$rbp_b, x = vis6$ln_agp_b, 
                                    df = vis6)), 
           x =0.5, y = 1.9) +
  annotate("text",  label =paste("12-18 mo: ",  lm_info(y = filter(smaller, visit_number == 8)$rbp_b, 
                                    x = filter(smaller, visit_number == 8)$ln_agp_b, 
                                    df = filter(smaller, visit_number == 8))),
           x =0.5, y = 1.8)
           
           
           
           

ggplot(data= smaller, aes(x = ln_agp_b, y = retinol_umolL_b)) + geom_point(size  = 2) +
   geom_smooth(method = "lm") + 
   annotate("text",  label = lm_info(smaller$retinol_umolL_b, smaller$ln_agp_b, df = smaller), x = 0.5, y = 2)


## Plot overlay - CRP
smalldat_b <- subset(smaller, !is.na(retinol_umolL_b), select = c("twinid", "rbp_b", "retinol_umolL_b", "crp_b", "ln_crp_b"))
meltedb <- melt(smalldat_b, id = c("twinid", "crp_b", "ln_crp_b"))
ggplot(aes(x = ln_crp_b, y = value, colour = variable), data = meltedb) + geom_point(size = 2) + geom_smooth(method = "lm") +
  ggtitle("Retinol and RBP by CRP, Infants") +
  annotate("text",  label = lm_info(y = smaller$retinol_umolL_b, 
                                    x = smaller$ln_crp_b , 
                                    df = smaller),
           x =2.5, y = 1.3)+
  annotate("text",  label = lm_info(y = smaller$rbp_b, 
                                    x = smaller$ln_crp_b , 
                                    df = smaller),
           x =2.5, y = 0.5)

smalldat_m <- subset(smaller, !is.na(retinol_umolL_m), select = c("twinid", "rbp_m", "retinol_umolL_m", "crp_m", "ln_crp_m"))
meltedm <- melt(smalldat_m, id = c("twinid", "crp_m", "ln_crp_m"))
ggplot(aes(x = ln_crp_m, y = value, colour = variable), data = meltedm) + geom_point(size = 2) + geom_smooth(method = "lm") +
  ggtitle("Retinol and RBP by CRP, Mothers") +
  annotate("text",  label = lm_info(y = filter(meltedm, variable == "rbp_m")$value, 
                                    x = filter(meltedm, variable == "rbp_m")$ln_crp_m , 
                                    df = filter(meltedm, variable == "rbp_m")),
           x =0, y = 1)+
  annotate("text",  label = lm_info(y = filter(meltedm, variable == "retinol_umolL_m")$value, 
                                    x = filter(meltedm, variable == "retinol_umolL_m")$ln_crp_m , 
                                    df = filter(meltedm, variable == "retinol_umolL_m")),
           x =2, y = 2.5)


## Plot overlay - AGP
smalldat_b <- subset(smaller, !is.na(retinol_umolL_b), select = c("twinid", "rbp_b", "retinol_umolL_b", "agp_b", "ln_agp_b"))
meltedb <- melt(smalldat_b, id = c("twinid", "agp_b", "ln_agp_b"))
ggplot(aes(x = ln_agp_b, y = value, colour = variable), data = meltedb) + geom_point(size = 2) + geom_smooth(method = "lm")+
  ggtitle("Retinol and RBP by AGP, Infants") +
  annotate("text",  label = lm_info(y = smaller$retinol_umolL_b, 
                                    x = smaller$ln_crp_b , 
                                    df = smaller),
           x =0.5, y = 1.3)+
  annotate("text",  label = lm_info(y = smaller$rbp_b, 
                                    x = smaller$ln_crp_b , 
                                    df = smaller),
           x =0.5, y = 0.5)


smalldat_m <- subset(smaller, !is.na(retinol_umolL_m), select = c("twinid", "rbp_m", "retinol_umolL_m", "agp_m", "ln_agp_m"))
meltedm <- melt(smalldat_m, id = c("twinid", "agp_m", "ln_agp_m"))
ggplot(aes(x = ln_agp_m, y = value, colour = variable), data = meltedm) + geom_point(size = 2) + geom_smooth(method = "lm")+
  ggtitle("Retinol and RBP by AGP, Mothers") +
  annotate("text",  label = lm_info(y = filter(meltedm, variable == "rbp_m")$value, 
                                    x = filter(meltedm, variable == "rbp_m")$ln_agp_m , 
                                    df = filter(meltedm, variable == "rbp_m")),
           x =-1, y = 1)+
  annotate("text",  label = lm_info(y = filter(meltedm, variable == "retinol_umolL_m")$value, 
                                    x = filter(meltedm, variable == "retinol_umolL_m")$ln_agp_m , 
                                    df = filter(meltedm, variable == "retinol_umolL_m")),
           x =0, y = 2.5)

## AGP and RBP - neg for 8, null for 6, pos for 2
ggplot(data= smaller, aes(x = ln_agp_b, y = rbp_b, colour = visfact)) + geom_point(size  = 2) +
  scale_color_manual(values = c("red", "green", "blue")) + geom_smooth(method = "lm")
  
## RBP and Retinol -- null for 2, pos for 6 and 8
ggplot(data= smaller, aes(x = rbp_b, y = retinol_umolL_b, colour = visfact)) + geom_point(size  = 2) +
  scale_color_manual(values = c("red", "green", "blue")) + geom_smooth(method = "lm")

ggplot(data= smaller, aes(x = rbp_b, y = retinol_umolL_b)) + geom_point(size  = 2) +
   geom_smooth(method = "lm")

rbp_ret <- lm_info(x =smaller$rbp_b, y = smaller$retinol_umolL_b, df = smaller)


## CORRECT RETINOL FOR INFLAMMATION? Does not seem to be an effect unless stratified -- but small data...
visit2.b.ret <- renamefn(newcorrfn("retinol_umolL_b", 2, small, "b"), "retinol", "b")
visit6.b.ret <- renamefn(newcorrfn("retinol_umolL_b", 6, small, "b"), "retinol", "b")
visit8.b.ret <- renamefn(newcorrfn("retinol_umolL_b", 8, small, "b"), "retinol", "b")

## short correlation function
shortcorrfn <- function(x, data, mom){
  # First just filter dataset to make smaller and more manageable
  subdat <-data
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
  
  # Run linear regressions and save values
  linreg <- lm(ln_x ~ ln_CRP + ln_AGP, subdat3)
  coef <- coef(linreg)
  linreg.int <- lm(ln_x ~ ln_CRP + ln_AGP + ln_CRP*ln_AGP, subdat3)
  coef.int <- coef(linreg.int)
  
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
  subdat3$x_th1 <- ifelse(subdat3$inflcat == 0, subdat3$x, 
                          ifelse(subdat3$inflcat == 1, subdat3$x / cf[2,10],
                                 ifelse(subdat3$inflcat == 2, subdat3$x / cf[3,10],
                                        ifelse(subdat3$inflcat == 3, 
                                               subdat3$x / cf[4,10],
                                               NA))))
  subdat3$x_th1 <- unlist(subdat3$x_th1)
  
  ## Linear regressions
  
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
  
    return(subdat3)
  
}

shortcorrfn_vis <- function(x, visit, data, mom){
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
  
  # Run linear regressions and save values
  linreg <- lm(ln_x ~ ln_CRP + ln_AGP, subdat3)
  coef <- coef(linreg)
  linreg.int <- lm(ln_x ~ ln_CRP + ln_AGP + ln_CRP*ln_AGP, subdat3)
  coef.int <- coef(linreg.int)
  
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
  subdat3$x_th1 <- ifelse(subdat3$inflcat == 0, subdat3$x, 
                          ifelse(subdat3$inflcat == 1, subdat3$x / cf[2,10],
                                 ifelse(subdat3$inflcat == 2, subdat3$x / cf[3,10],
                                        ifelse(subdat3$inflcat == 3, 
                                               subdat3$x / cf[4,10],
                                               NA))))
  subdat3$x_th1 <- unlist(subdat3$x_th1)
  
  ## Linear regressions
  
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
  
  return(subdat3)
  
}

# Function to rename interim variables 
retnamefn <- function(data, var, mom){
  if(var == "retinol" & mom == "m"){
    newdat <- rename(data, ln_ret_m = ln_x, 
                     ret_m_th.1.1 = x_th1,
                     ret_m_lc.5 = x_lc.5,
                     ret_m_lc.3.5 = x_lc.3.5)}
  if(var == "retinol" & mom == "b"){
    newdat <- rename(data, ln_ret_b = ln_x, 
                     ret_b_th.1.1 = x_th1,
                      ret_b_lc.5 = x_lc.5,
                     ret_b_lc.3.5 = x_lc.3.5)}
  
  newdat <- mutate(newdat, x = NULL, CRP = NULL, AGP = NULL, ln_CRP = NULL,
                   ln_AGP = NULL, inflcat = NULL)
  return(newdat)
}

## LOOK AT 
ret.b <- retnamefn(shortcorrfn("retinol_umolL_b", small, "b"), "retinol", "b")
ret.m <- retnamefn(shortcorrfn("retinol_umolL_m", small, "m"), "retinol", "m")

corr_ret <- merge(ret.b, ret.m, by = c("twinid", "visit_number")) %>%
  mutate(agp_b = NULL, inflcat_b = NULL, crp_b = NULL, crp_m = NULL, agp_m = NULL, inflcat_m = NULL,
         ln_ret_b = NULL, ln_ret_m = NULL, ln_crp_b = NULL, ln_agp_b = NULL, ln_agp_m = NULL, ln_crp_m = NULL)

small_wcorret <- merge(small, corr_ret, by = c("twinid", "visit_number")) %>%
  mutate(retinol_umolL_m.y = NULL, retinol_umolL_b.y = NULL) %>%
  rename(retinol_umolL_m = retinol_umolL_m.x, retinol_umolL_b = retinol_umolL_b.x)

## Now correct retinol by visit number
ret.b.vis2 <- retnamefn(shortcorrfn_vis("retinol_umolL_b", 2, small_wcorret, "b"), "retinol", "b")
ret.b.vis6 <- retnamefn(shortcorrfn_vis("retinol_umolL_b", 6, small_wcorret, "b"), "retinol", "b")
ret.b.vis8 <- retnamefn(shortcorrfn_vis("retinol_umolL_b", 8, small_wcorret, "b"), "retinol", "b")
ret_b_vis <- rbind(ret.b.vis2, ret.b.vis6, ret.b.vis8)  %>%
  mutate(agp_b = NULL, inflcat_b = NULL, crp_b = NULL, crp_m = NULL, agp_m = NULL, inflcat_m = NULL,
         ln_ret_b = NULL, ln_ret_m = NULL, ln_crp_b = NULL, ln_agp_b = NULL, ln_agp_m = NULL, ln_crp_m = NULL)
## Now correct retinol by visit number-- MOMS
ret.m.vis1 <- retnamefn(shortcorrfn_vis("retinol_umolL_m", 1, small_wcorret, "m"), "retinol", "m")
ret.m.vis6 <- retnamefn(shortcorrfn_vis("retinol_umolL_m", 6, small_wcorret, "m"), "retinol", "m")
ret_m_vis <- rbind(ret.m.vis1, ret.m.vis6)  %>%
  mutate(agp_b = NULL, inflcat_b = NULL, crp_b = NULL, crp_m = NULL, agp_m = NULL, inflcat_m = NULL,
         ln_ret_b = NULL, ln_ret_m = NULL, ln_crp_b = NULL, ln_agp_b = NULL, ln_agp_m = NULL, ln_crp_m = NULL)
# Merge them together
ret_vis <- merge(ret_b_vis, ret_m_vis, by = c("twinid", "visit_number"), all = TRUE)


## Merge again
small_w2corret <- merge(ret_vis, small_wcorret, by = c("twinid", "visit_number"), all = TRUE) %>%
  mutate(retinol_umolL_b.x = NULL, retinol_umolL_m.x = NULL) %>R%
  rename(retinol_umolL_b = retinol_umolL_b.y, retinol_umolL_m = retinol_umolL_m.y,
         ret_b_lc.3.5_vis = ret_b_lc.3.5.x,  ret_b_lc.3.5_all = ret_b_lc.3.5.y, 
         ret_m_lc.3.5_vis = ret_m_lc.3.5.x,  ret_m_lc.3.5_all = ret_m_lc.3.5.y, 
         ret_b_th.1.1_vis = ret_b_th.1.1.x,  ret_b_th.1.1_all = ret_b_th.1.1.y, 
         ret_m_th.1.1_vis = ret_m_th.1.1.x,  ret_m_th.1.1_all = ret_m_th.1.1.y, 
         ret_b_lc.5_vis = ret_b_lc.5.x,  ret_b_lc.5_all = ret_b_lc.5.y, 
         ret_m_lc.5_vis = ret_m_lc.5.x,  ret_m_lc.5_all = ret_m_lc.5.y) %>%
  mutate(vad_ret_b_un = ifelse(retinol_umolL_b < 0.7, 1, ifelse(!is.na(retinol_umolL_b), 0, NA)),
         vad_ret_b_3.5_vis = ifelse(ret_b_lc.3.5_vis < 0.7, 1, ifelse(!is.na(ret_b_lc.3.5_vis), 0, NA)),
         vad_ret_m_3.5_vis = ifelse(ret_m_lc.3.5_vis < 0.7, 1, ifelse(!is.na(ret_m_lc.3.5_vis), 0, NA)),
         vad_ret_b_3.5_all = ifelse(ret_b_lc.3.5_all < 0.7, 1, ifelse(!is.na(ret_b_lc.3.5_all), 0, NA)),
         vad_ret_m_3.5_all = ifelse(ret_m_lc.3.5_all < 0.7, 1, ifelse(!is.na(ret_m_lc.3.5_all), 0, NA)),
         vad_rbp_b_un_33 = ifelse(rbp_b < 0.33, 1, ifelse(!is.na(rbp_b), 0, NA)))


## Add 




## Now get numbers for table
summary(baby2comp$ret_b_lc.3.5_vis)
sd(baby2comp$ret_b_lc.3.5_vis, na.rm = TRUE)

summary(baby_two.six$ret_b_lc.3.5_vis)
sd(baby_two.six$ret_b_lc.3.5_vis, na.rm = TRUE)

summary(baby_all3$ret_b_lc.3.5_vis)
sd(baby_all3$ret_b_lc.3.5_vis, na.rm = TRUE)


table(baby2comp$vad_rbp_b_un_33)
table(baby2comp$vad_ret_b_3.5_vis)



## Look at relationship of supplementation to VAD - not sig in models, but means are supposedly different...
###NS for logistic mod
vasmod <- glm(vad_rbp_b_lc.3.5 ~ vita_pre6, data = filter(small_w2corret, visit_number == 6), family = binomial)
### try linmod
vasmodlin <- lm(log(rbp_b_lc.3.5)~ vita_pre6, data =  filter(small_w2corret, visit_number == 6))
summary(vasmodlin)
t.test(filter(small_w2corret, visit_number == 6)$rbp_b_lc.3.5, filter(small_w2corret, visit_number == 6)$vita_pre6)
wilcox.test(filter(small_w2corret, visit_number == 6)$rbp_b_lc.3.5, filter(small_w2corret, visit_number == 6)$vita_pre6)


### Look at anemia in infatns VAD at 2 mo
vis2 <- filter(small_w2corret, visit_number == 2)
chisq.test(table(vis2$anemia_b_alt2_v6, vis2$vad_rbp_b_lc.3.5))

chisq.test(table(vis2$anemia_b_alt2_v6, vis2$vad_rbp_b_un_33))
fisher.test(table(vis2$anemia_b_alt2_v6, vis2$vad_ret_b_3.5_vis))


#### Look at vitamin A status and later effects
ggplot(data= vis2, aes(x = rbp_b, y = num.anymorb2to6)) + geom_point(size  = 2) +
  geom_vline(xintercept = 0.7, linetype = 2) + geom_vline(xintercept = 0.33, linetype = 2) + 
  geom_smooth(method = "lm")

ggplot(data= vis2, aes(x = rbp_b, y = num.morb.2to6)) + geom_point(size  = 2) +
  geom_vline(xintercept = 0.7, linetype = 2) + geom_vline(xintercept = 0.33, linetype = 2) + 
  geom_smooth(method = "lm")

ggplot(data= vis2, aes(x = rbp_b, y = num.diarr2to6)) + geom_point(size  = 2) +
  geom_vline(xintercept = 0.7, linetype = 2) + geom_vline(xintercept = 0.33, linetype = 2) + 
  geom_smooth(method = "lm")

ggplot(data= vis2, aes(x = rbp_b, y = num.cough2to6)) + geom_point(size  = 2) +
  geom_vline(xintercept = 0.7, linetype = 2) + geom_vline(xintercept = 0.33, linetype = 2) + 
  geom_smooth(method = "lm")

ggplot(data= vis2, aes(x = rbp_b, y = num.fever2to6)) + geom_point(size  = 2) +
  geom_vline(xintercept = 0.7, linetype = 2) + geom_vline(xintercept = 0.33, linetype = 2) + 
  geom_smooth(method = "lm")

# Retinol

ggplot(data= vis2, aes(x = retinol_umolL_b, y = num.anymorb2to6)) + geom_point(size  = 2) +
  geom_vline(xintercept = 0.7, linetype = 2) + geom_vline(xintercept = 0.33, linetype = 2) + 
  geom_smooth(method = "lm")

ggplot(data= vis2, aes(x = retinol_umolL_b, y = num.morb.2to6)) + geom_point(size  = 2) +
  geom_vline(xintercept = 0.7, linetype = 2) + geom_vline(xintercept = 0.33, linetype = 2) + 
  geom_smooth(method = "lm")

ggplot(data= vis2, aes(x = retinol_umolL_b, y = num.diarr2to6)) + geom_point(size  = 2) +
  geom_vline(xintercept = 0.7, linetype = 2) + geom_vline(xintercept = 0.33, linetype = 2) + 
  geom_smooth(method = "lm")

ggplot(data= vis2, aes(x = retinol_umolL_b, y = num.cough2to6)) + geom_point(size  = 2) +
  geom_vline(xintercept = 0.7, linetype = 2) + geom_vline(xintercept = 0.33, linetype = 2) + 
  geom_smooth(method = "lm")

ggplot(data= vis2, aes(x = retinol_umolL_b, y = num.fever2to6)) + geom_point(size  = 2) +
  geom_vline(xintercept = 0.7, linetype = 2) + geom_vline(xintercept = 0.33, linetype = 2) + 
  geom_smooth(method = "lm")


boxplot( rbp_b~ num.anymorb2to6, vis2)

## Attempt poisson model
vis2 <- mutate(vis2, monthsv26 = agemo.v6 - agemo.v2) %>%
  mutate(logmon2_6 = log(monthsv26))


poissmod_crude <- glm(formula = num.anymorb2to6 ~ vad_b_un + offset(logmon2_6),
                       family = poisson, data = vis2)



poissmod_adj <- glm(formula = num.anymorb2to6 ~ vad_b_un + m_job + offset(logmon2_6),
                      family = poisson, data = vis2)


poissmod_crude_rbp <- glm(formula = num.anymorb2to6 ~ rbp_b + offset(logmon2_6),
                      family = poisson, data = vis2)

poissmod_crude_ret <- glm(formula = num.anymorb2to6 ~ retinol_umolL_b + offset(logmon2_6),
                          family = poisson, data = vis2)

poissmod_adj_rbp <- glm(formula = num.anymorb2to6 ~ rbp_b + crp_b + agp_b + m_job + offset(logmon2_6),
                        family = poisson, data = vis2)


poissmod_adj_ret <- glm(formula = num.anymorb2to6 ~  retinol_umolL_b + crp_b + agp_b + m_job + offset(logmon2_6),
                        family = poisson, data = vis2)

## Try with excluding morbidities at visit 2
poissmod_crude <- glm(formula = num.anymorbpost2to6 ~ vad_b_un + offset(logmon2_6),
                      family = poisson, data = vis2)


poissmod_crude_rbp <- glm(formula = num.anymorb2to6 ~ rbp_b + offset(logmon2_6),
                          family = poisson, data = vis2)

poissmod_crude_ret <- glm(formula = num.anymorb2to6 ~ retinol_umolL_b + offset(logmon2_6),
                          family = poisson, data = vis2)

poissmod_adj_rbp <- glm(formula = num.anymorb2to6 ~ rbp_b + crp_b + agp_b + m_job + offset(logmon2_6),
                        family = poisson, data = vis2)


poissmod_adj_ret <- glm(formula = num.anymorb2to6 ~  retinol_umolL_b + crp_b + agp_b + m_job + offset(logmon2_6),
                        family = poisson, data = vis2)



poissmod_crude_lc <- glm(formula = num.anymorbpost2to6 ~ vad_rbp_b_lc.3.5 + offset(logmon2_6),
                      family = poisson, data = vis2)

poissmod_crude_rbp_lc <- glm(formula = num.anymorbpost2to6 ~ rbp_b_lc.3.5 + offset(logmon2_6),
                          family = poisson, data = vis2)

poissmod_adj_rbp_post <- glm(formula = num.anymorbpost2to6 ~ rbp_b + crp_b + agp_b + m_job + offset(logmon2_6),
                        family = poisson, data = vis2)

poissmod_adj_rbp_post_lc <- glm(formula = num.anymorbpost2to6 ~ rbp_b_lc.3.5 + m_job + offset(logmon2_6),
                             family = poisson, data = vis2)


## Try with total morbidities






## Export small dataset for donnie

ret_rbp_dat <- subset(small_w2corret, visit_number %in% c(1,2,6,8),
                      select = c("twinid", "visit_number", "male", "date_visita_mmddyyyy", "rbp_b", "rbp_m", "crp_b", "crp_m",
                                 "agp_b", "agp_m", "ret_m_lc.3.5_all", "ret_b_lc.3.5_all", "ret_m_lc.3.5_vis", "ret_b_lc.3.5_vis", 
                                 "rbp_b_lc.3.5", "rbp_m_lc.3.5", "retinol_umolL_b", "retinol_umolL_m", "retinol_ugL_b",
                                 "retinol_ugL_m"))

write.csv(ret_rbp_dat , "ret_rbp_dat.csv")
