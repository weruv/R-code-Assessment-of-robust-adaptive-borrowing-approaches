library(RBesT)
library(parallel)
library(ggplot2)
library(latex2exp)
library(ggh4x)
library(cowplot)

## ===========================================================================================================
## average TIE rate and power- hybrid control
## ===========================================================================================================
source("function_avgOC.R")

nt <- 20                                             # sample size in treatment arm
nc <- 20                                             # sample size in control arm
n_ext <- 15                                         # historical data sample size
mu_ext <- 0                       ## external data mean


design_prior_type <-  c("info", "mix", "unit","fixed_pp","vague")          ## different kinds of design priors

design_prior_mean<-seq(-2.5, 2.5, by = 0.01)

grid <- expand.grid(nt = nt, nc = nc, n_ext = n_ext, mu_ext = mu_ext, design_prior_type = design_prior_type,design_prior_mean=design_prior_mean)

ncores <- 80

avgOCs <- do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))),
                                  function(params) {
                                    with(params, {
                                      sigma <- 1
                                      sig2 <- sigma^2                  # variance of one observation
                                      
                                      sigmapc <- sqrt(sig2 / n_ext)    # SD control arm prior
                                      sigmapt <- sqrt(sig2 / 0)        # SD treat arm prior
                                      muc <- 0                         # prior mean control
                                      
                                      mut <- 0.4                       # prior mean treatment
                                      priorPars <- rbind(c(muc, mut), c(sigmapc, sigmapt))
                                      
                                      sigmapc_downweighted_external <- sqrt(sig2 / (n_ext*0.5))    # downweight external data
                                      alpha.b <- 0.025              # frequentist TIE
                                      th0 <- 0
                                      
                                      bound <- FALSE #
                                      nmc <- 1e6                        # number of Monte Carlo samples
                                      
                                      th1 <- 0.8330874                    # point for Power evaluation
                                      delta <- c(th0, th1)
                                      
                                      p_parameter <- c(0.2,0.3,0.5,0.8) # power parameter of the fixed power analysis prior
                                      unit_sd <- 1  ## 1 for the case considered here, not always one
                                      
                                      grid_mixture <- expand.grid(robust_component_sd = c(1, 100), robust_comp_location = c("ext", "current.mean"), prior_weight = c(0.5))
                                      
                                      
                                      set.seed(123)
                                      Avg.OC <- sapply(delta, function(x) {
                                        if (design_prior_type == "info") {
                                          thc <- rnorm(nmc, priorPars[1, 1]+design_prior_mean, priorPars[2, 1])
                                        } else if (design_prior_type == "mix") {
                                          thc <- rmix(mixnorm(informative = c(0.5, c(priorPars[1, 1]+design_prior_mean, priorPars[2, 1])), vague = c(1 - 0.5, c(priorPars[1, 1]+design_prior_mean, unit_sd)), sigma = 1), nmc)
                                        } else if (design_prior_type == "unit") {
                                          thc <- rnorm(nmc, priorPars[1, 1]+design_prior_mean, sqrt(sig2))
                                        }else if (design_prior_type == "fixed_pp") {
                                          thc <- rnorm(nmc, priorPars[1, 1]+design_prior_mean, sigmapc_downweighted_external)
                                        } else if (design_prior_type == "vague") {
                                          thc <- rnorm(nmc, priorPars[1, 1]+design_prior_mean, 10^4)
                                        }
                                        
                                        # sample data
                                        yc <- rnorm(nmc, thc, sigma / sqrt(nc))
                                        yt <- rnorm(nmc, thc + x, sigma / sqrt(nt))
                                        
                                        
                                        # Weight vector, 1 for monte carlo
                                        wmatT <- rep(1 / nmc, nmc)
                                        
                                        priorPars.v <- priorPars
                                        priorPars.v[2, ] <- c(Inf, Inf)
                                        p.h1.unif <- pH1(
                                          yc = yc, yt = yt, nc = nc, nt = nt,
                                          priorPars = priorPars.v, w.mix = w.mix, th0 = th0, robust_comp_location = robust_comp_location, dE = mu_ext, power_parameter = power_param
                                        )
                                        
                                        
                                        p.h1.info <- pH1(
                                          yc = yc, yt = yt, nc = nc, nt = nt,
                                          priorPars = priorPars, w.mix = w.mix, th0 = th0, robust_comp_location = robust_comp_location, dE = mu_ext, power_parameter = power_param
                                        )
                                        
                                        p.h1.mix <- mapply(
                                          function(prior_weight, robust_comp_location, robust_component_sd) {
                                            pH1(
                                              prior = "mix", yc = yc, yt = yt, nc = nc, nt = nt,
                                              priorPars = priorPars, w.mix = prior_weight, th0 = th0, robust_comp_location = robust_comp_location, dE = mu_ext, robust_comp_sd = robust_component_sd, power_parameter = power_param
                                            )
                                          },
                                          grid_mixture$prior_weight, grid_mixture$robust_comp_location, grid_mixture$robust_component_sd
                                        )
                                        
                                        
                                        
                                        p.h1.fixed_pp <- sapply(p_parameter, function(power_parameter) {
                                          pH1(
                                            prior = "fixed_pp", yc = yc, yt = yt, nc = nc, nt = nt,
                                            priorPars = priorPars, w.mix = w.mix, th0 = th0, robust_comp_location = robust_comp_location, dE = mu_ext, power_parameter = power_parameter
                                          )
                                        })
                                        
                                        
                                        p.h1.EB <- pH1(
                                          prior = "EB", yc = yc, yt = yt, nc = nc, nt = nt,
                                          priorPars = priorPars, w.mix = w.mix, bound = bound, th0 = th0, robust_comp_location = robust_comp_location, dE = mu_ext, power_parameter = power_param
                                        )
                                        
                                        # test decision
                                        p.h1.list <- cbind(p.h1.info, p.h1.unif, p.h1.mix, p.h1.fixed_pp, p.h1.EB)
                                       
                                        
                                        dec.list <- cbind(
                                          p.h1.info > (1 - alpha.b),
                                          p.h1.unif > (1 - alpha.b),
                                          p.h1.mix > (1 - alpha.b),
                                          p.h1.fixed_pp > (1 - alpha.b),
                                          p.h1.EB > (1 - alpha.b)
                                        )
                                        
                                        
                                        AvgOC <- wmatT %*% dec.list
                                        AvgOC
                                      })
                                      names_fixedpp <- paste0("fixed_pp", "_", p_parameter)
                                      names_mix <- do.call(paste, c(list("mix"), grid_mixture, sep = "_"))
                                      rownames(Avg.OC) <- c("info", "freq", names_mix, names_fixedpp, "EB")
                                      colnames(Avg.OC) <- c("average_TIE", "average_Power")
                                      
                                      Avg.OC <- as.data.frame(Avg.OC)
                                      Avg.OC <- cbind(Avg.OC, mu_ext = mu_ext, nt, nc, n_ext, design_prior_type = design_prior_type,design_prior_mean)
                                      Avg.OC$analysis_prior_type <- rownames(Avg.OC)
                                      
                                      Avg.OC$robust_location <- ifelse(grepl("mix", Avg.OC$analysis_prior_type), sub(".*_([^_]+)_[^_]+$", "\\1", Avg.OC$analysis_prior_type), NA)
                                      Avg.OC$robust_sd <- as.numeric(ifelse(grepl("mix", Avg.OC$analysis_prior_type), sub("^[^_]+_([^_]+)_.*", "\\1", Avg.OC$analysis_prior_type), NA))
                                      Avg.OC$w.mix <- as.numeric(ifelse(grepl("mix", Avg.OC$analysis_prior_type), sub(".*_", "", Avg.OC$analysis_prior_type), NA))
                                      Avg.OC$power_parameter <- as.numeric(ifelse(grepl("fixed", Avg.OC$analysis_prior_type), sub(".*_", "", Avg.OC$analysis_prior_type), NA))
                                      return(Avg.OC)
                                    })
                                  },
                                  mc.cores = ncores
))


## ================================================================================================================ ##
## average TIE rate and power for different design priors but informative analysis prior
## ================================================================================================================ ##

avgOCs$analysis_prior_type <- as.factor(avgOCs$analysis_prior_type)

avgOCs_info_analysisonly <- avgOCs[avgOCs$analysis_prior_type %in% c("info"), ]

avgOCs_info_analysisonly <- avgOCs_info_analysisonly[avgOCs_info_analysisonly$design_prior_type!="vague", ]

avgOCs_info_analysisonly$design_prior_type <- factor(avgOCs_info_analysisonly$design_prior_type)

levels(avgOCs_info_analysisonly$design_prior_type) <- c(
  fixed_pp = TeX("$\\Fixed \\, power \\, prior"),
  info = TeX("$Informative"),
  mix = TeX("$\\RMP$"),
  unit = TeX("$\\Unit-information$")
)

avg_tie_n_power_labels <- colnames(avgOCs_info_analysisonly)[grepl("average_TIE|average_Power", colnames(avgOCs_info_analysisonly))]
avg_tie_n_power_labels2 <- c("average TIE rate", "average power")
avg_tie_n_power_labels3 <- c(0.025, 0.75)
avg_tie_n_power_labels4 <- list(c(0.4, 0.8), c(0.69, 0.25))


avgOCs_plot<-lapply(1:2, function(index){
  ggplot() +
    geom_line(data = avgOCs_info_analysisonly, aes(x =design_prior_mean , y = .data[[avg_tie_n_power_labels[index]]], color = as.factor(design_prior_type)), size = 1.2)+
    xlab("")+
    scale_color_discrete(labels = scales::label_parse())+
    guides(
      color = guide_legend(label.hjust = 0, title = "Design prior used")
    )+
    xlab(TeX("$\\Design \\, prior\\, mean\\, (\\mu^d)$"))+
    ylab(avg_tie_n_power_labels2[index])+
    theme_bw()+
    theme(axis.text = element_text(color = c("black", "black", "black", "black", "black")),
          axis.ticks = element_line(color = c("black", "black", "black", "black", "black")))+
    theme(axis.text=element_text(size=32,family = "serif"),
          axis.title.x=element_text(size=32,family = "serif"),
          axis.title.y=element_text(size=32,family = "serif"),
          legend.text = element_text(size=26,family = "serif"),
          legend.title=element_text(size=26,family = "serif"),
          legend.position =avg_tie_n_power_labels4[[index]],
          legend.box="vertical",
           legend.key.width =  unit(3, "line"),
          legend.spacing.x = unit(0.25, 'cm'),
          legend.spacing.y = unit(0.3, 'cm'),
          strip.text.x = element_text(size=26,family = "serif"),
          strip.text.y = element_text(size=26,family = "serif")
    )+
    theme(panel.background = element_rect(fill = 'white', colour =NA))+
    theme(axis.line = element_line(colour = "black",size=1),
          panel.grid.minor = element_blank(), 
          panel.border = element_blank()) 
})


## ================================================================================================================ ##
## average TIE rate and power for different design priors but fixed power analysis prior
## ================================================================================================================ ##
plot_names <- c("fixed_pp_0.2", "fixed_pp_0.5", "fixed_pp_0.8")

y_axis_variable <- colnames(avgOCs)[grepl("average_TIE|average_Power", colnames(avgOCs))]
y_axis_variable_labels <- c("average TIE rate", "average power")
y_axis_hlines <- c(0.025, 0.75)

avgOCs_fixedpp_as_analysisprior<-lapply(1:2, function(index){
  lapply(plot_names, function(plot){
    
    title_name <- sub(".*_", "", plot)
    avgOCs_powerprioronly <- avgOCs[grepl(plot, avgOCs$analysis_prior_type), ]
    avgOCs_powerprioronly <- avgOCs_powerprioronly[avgOCs_powerprioronly$design_prior_type!="vague", ]
    avgOCs_powerprioronly$design_prior_type <- factor(avgOCs_powerprioronly$design_prior_type)

    levels(avgOCs_powerprioronly$design_prior_type) <- c(
      fixed_pp = TeX("$\\Fixed \\, power \\, prior"),
      info = TeX("$Informative"),
      mix = TeX("$\\RMP$"),
      unit = TeX("$\\Unit-information$")
    )
    
    ggplot() +
     geom_line(data = avgOCs_powerprioronly, aes(x =design_prior_mean , y =  .data[[y_axis_variable[index]]], color = as.factor(design_prior_type)), size = 1.2)+
      xlab("")+
      scale_color_discrete(labels = scales::label_parse())+
      guides(
        color = guide_legend(label.hjust = 0, title = "Design prior used")
      )+
      xlab(TeX("$\\Design\\, prior\\, mean\\, (\\mu^d)$"))+
      ylab(y_axis_variable_labels[index])+
      scale_y_continuous(limits = c(0,1),breaks = c(0, 0.25, 0.5, .75, 1))+
      ggtitle(TeX(paste0("$\\alpha_0 = ", title_name, "$")))+
      theme_bw()+
      theme(axis.text = element_text(color = c("black", "black", "black", "black", "black")),
            axis.ticks = element_line(color = c("black", "black", "black", "black", "black")))+
      theme(axis.text=element_text(size=30,family = "serif"),
            axis.title.x=element_text(size=30,family = "serif",margin = margin(b = 25)),
            axis.title.y=element_text(size=30,family = "serif"),
            legend.text = element_text(size=30,family = "serif"),
            legend.title=element_text(size=30,family = "serif"),
            legend.position ="bottom",
            legend.box="vertical",
            legend.key.width =  unit(3, "line"),
            legend.spacing.x = unit(0.25, 'cm'),
            legend.spacing.y = unit(0.3, 'cm'),
            strip.text.x = element_text(size=26,family = "serif"),
            strip.text.y = element_text(size=26,family = "serif"),
            plot.title = element_text(
              size = 26,
              hjust = 0.5
            ),
            plot.margin = margin(25,25, 35, 25),
            legend.box.margin = margin(t = 25,b=25),
            legend.margin = margin(b = 25)
      )+
      theme(panel.background = element_rect(fill = 'white', colour =NA))+
      theme(axis.line = element_line(colour = "black",size=1),
            panel.grid.minor = element_blank(), 
            panel.border = element_blank())
    
  })
})

## ================================================================================================================ ##
## average TIE rate and power for different design priors but dynamic analysis priors
## ================================================================================================================ ##
avgOCs_EB_n_mixture_analysisonly <- avgOCs[avgOCs$analysis_prior_type %in% c("EB", "mix_1_ext_0.5", "mix_1_current.mean_0.5"), ]
avgOCs_EB_n_mixture_analysisonly$analysisprior_title <- "Analysis prior used"
avgOCs_EB_n_mixture_analysisonly$analysis_prior_type <- factor(avgOCs_EB_n_mixture_analysisonly$analysis_prior_type)

avgOCs_EB_n_mixture_analysisonly <- avgOCs_EB_n_mixture_analysisonly[avgOCs_EB_n_mixture_analysisonly$design_prior_type!="vague", ]

levels(avgOCs_EB_n_mixture_analysisonly$analysis_prior_type) <- c(
  EB = TeX("$\\EB \\, power \\, prior"),
  mix_1_current.mean_0.5 = expression(atop("RMP", (mu[robust] == bar(y)[c] ~ ", unit-info, w=0.5"))),
  mix_1_ext_0.5 = expression(atop("RMP", (mu[robust] == bar(y)[ext] ~ ", unit-info, w=0.5")))
)
avgOCs_EB_n_mixture_analysisonly$design_prior_type <- factor(avgOCs_EB_n_mixture_analysisonly$design_prior_type)

levels(avgOCs_EB_n_mixture_analysisonly$design_prior_type) <- c(
  fixed_pp = TeX("$\\Fixed \\, power \\, prior"),
  info = TeX("$Informative"),
  mix = TeX("$\\RMP$"),
  unit = TeX("$\\Unit-information$")
)

avg_tie_n_power_labels <- colnames(avgOCs_EB_n_mixture_analysisonly)[grepl("average_TIE|average_Power", colnames(avgOCs_EB_n_mixture_analysisonly))]
avg_tie_n_power_labels2 <- c("average TIE rate", "average power")
avg_tie_n_power_labels3 <- c(0.025, 0.75)
avg_tie_n_power_labels4 <- list(c("bottom"), c("bottom"))

avgOCs_plot2<-lapply(1:2, function(index){
  ggplot() +
    geom_hline(yintercept = avg_tie_n_power_labels3[index], size = 1.2) +
    geom_line(data = avgOCs_EB_n_mixture_analysisonly, aes(x =design_prior_mean , y = .data[[avg_tie_n_power_labels[index]]], color = as.factor(design_prior_type)), size = 1.2)+
    facet_nested(.~analysisprior_title*analysis_prior_type,labeller = labeller(analysis_prior_type = label_parsed),
                 strip = strip_nested(text_x = elem_list_text(size = c(24,24,24,24))))+
    xlab("")+
    scale_color_discrete(labels = scales::label_parse())+
    guides(
      color = guide_legend(label.hjust = 0, title = "Design prior used",nrow = 2,title.position = "top",title.hjust = 0.5)
    )+
    xlab(TeX("$\\Design\\, prior\\, mean\\, (\\mu^d)$"))+
    ylab(avg_tie_n_power_labels2[index])+
    theme_bw()+
    theme(axis.text = element_text(color = c("black", "black", "black", "black", "black")),
          axis.ticks = element_line(color = c("black", "black", "black", "black", "black")))+
    theme(axis.text=element_text(size=32,family = "serif"),
          axis.title.x=element_text(size=32,family = "serif",margin = margin(b = 25)),
          axis.title.y=element_text(size=32,family = "serif"),
          legend.text = element_text(size=32,family = "serif"),
          legend.title=element_text(size=32,family = "serif"),
          legend.position =avg_tie_n_power_labels4[[index]],
          legend.box="vertical",
          legend.key.width =  unit(3, "line"),
          legend.spacing.x = unit(1, 'cm'),
          legend.spacing.y = unit(0.3, 'cm')
    )+
    theme(panel.background = element_rect(fill = 'white', colour =NA))+
    theme(axis.line = element_line(colour = "black",size=1),
          panel.grid.minor = element_blank(), 
          panel.border = element_blank())
})

## ================================================================================================================ ##
## Recalibrating the Bayesian test decision
## ================================================================================================================ ##

threshold_value <- seq(0, 1, by = 0.025)

grid <- expand.grid(nt = nt, nc = nc, n_ext = n_ext, mu_ext = mu_ext, design_prior_type = design_prior_type, threshold_value = threshold_value, design_prior_mean = design_prior_mean)

ncores <- 80

average_OC_recalibrateddecision <- do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))),
                                                           function(params) {
                                                             with(params, {
                                                               sigma <- 1
                                                               sig2 <- sigma^2                  # variance of one observation
                                                               
                                                               sigmapc <- sqrt(sig2 / n_ext)    # SD control arm prior
                                                               sigmapt <- sqrt(sig2 / 0)        # SD treat arm prior
                                                               muc <- 0                         # prior mean control
                                                               mut <- 0.4                       # prior mean treatment
                                                               priorPars <- rbind(c(muc, mut), c(sigmapc, sigmapt))
                                                               
                                                               sigmapc_downweighted <- sqrt(sig2 / (n_ext*0.5))    # # downweight external data
                                                               alpha.b <- threshold_value              # frequentist TIE
                                                               th0 <- 0
                                                               
                                                               bound <- FALSE #
                                                               nmc <- 1e6                        # number of Monte Carlo samples
                                                               
                                                               th1 <- 0.8330874                  # point for Power evaluation
                                                               delta <- c(th0, th1)
                                                               
                                                               p_parameter <- c(0.2,0.3,0.5,0.8) # power parameter of the fixed power analysis prior
                                                               unit_sd <- 1 ## 1 for the case considered here, not always one
                                                               
                                                               grid_mixture <- expand.grid(robust_component_sd = c(1, 100), robust_comp_location = c("ext", "current.mean"), prior_weight = c(0.25, 0.5, 0.75))
                                                               
                                                               
                                                               set.seed(123)
                                                               Avg.OC <- sapply(delta, function(x) {
                                                                 if (design_prior_type == "info") {
                                                                   thc <- rnorm(nmc, priorPars[1, 1]+design_prior_mean, priorPars[2, 1])
                                                                 } else if (design_prior_type == "mix") {
                                                                   thc <- rmix(mixnorm(informative = c(0.5, c(priorPars[1, 1]+design_prior_mean, priorPars[2, 1])), vague = c(1 - 0.5, c(priorPars[1, 1]+design_prior_mean, unit_sd)), sigma = 1), nmc)
                                                                 } else if (design_prior_type == "unit") {
                                                                   thc <- rnorm(nmc, priorPars[1, 1]+design_prior_mean, sqrt(sig2))
                                                                 }else if (design_prior_type == "fixed_pp") {
                                                                   thc <- rnorm(nmc, priorPars[1, 1]+design_prior_mean, sigmapc_downweighted_external)
                                                                 }else if (design_prior_type == "vague") {
                                                                   thc <- rnorm(nmc, priorPars[1, 1]+design_prior_mean, 10^4)
                                                                 }
                                                                 
                                                                 # sample data
                                                                 yc <- rnorm(nmc, thc, sigma / sqrt(nc))
                                                                 yt <- rnorm(nmc, thc + x, sigma / sqrt(nt))
                                                                 
                                                                 # Weight vector, 1 for monte carlo
                                                                 wmatT <- rep(1 / nmc, nmc)
                                                                 
                                                                 priorPars.v <- priorPars
                                                                 priorPars.v[2, ] <- c(Inf, Inf)
                                                                 p.h1.unif <- pH1(
                                                                   yc = yc, yt = yt, nc = nc, nt = nt,
                                                                   priorPars = priorPars.v, w.mix = w.mix, th0 = th0, robust_comp_location = robust_comp_location, dE = mu_ext, power_parameter = power_param
                                                                 )
                                                                 
                                                                 
                                                                 p.h1.info <- pH1(
                                                                   yc = yc, yt = yt, nc = nc, nt = nt,
                                                                   priorPars = priorPars, w.mix = w.mix, th0 = th0, robust_comp_location = robust_comp_location, dE = mu_ext, power_parameter = power_param
                                                                 )
                                                                 
                                                                 p.h1.mix <- mapply(
                                                                   function(prior_weight, robust_comp_location, robust_component_sd) {
                                                                     pH1(
                                                                       prior = "mix", yc = yc, yt = yt, nc = nc, nt = nt,
                                                                       priorPars = priorPars, w.mix = prior_weight, th0 = th0, robust_comp_location = robust_comp_location, dE = mu_ext, robust_comp_sd = robust_component_sd, power_parameter = power_param
                                                                     )
                                                                   },
                                                                   grid_mixture$prior_weight, grid_mixture$robust_comp_location, grid_mixture$robust_component_sd
                                                                 )
                                                                 
                                                                 
                                                                 
                                                                 p.h1.fixed_pp <- sapply(p_parameter, function(power_parameter) {
                                                                   pH1(
                                                                     prior = "fixed_pp", yc = yc, yt = yt, nc = nc, nt = nt,
                                                                     priorPars = priorPars, w.mix = w.mix, th0 = th0, robust_comp_location = robust_comp_location, dE = mu_ext, power_parameter = power_parameter
                                                                   )
                                                                 })
                                                                 
                                                                 
                                                                 p.h1.EB <- pH1(
                                                                   prior = "EB", yc = yc, yt = yt, nc = nc, nt = nt,
                                                                   priorPars = priorPars, w.mix = w.mix, bound = bound, th0 = th0, robust_comp_location = robust_comp_location, dE = mu_ext, power_parameter = power_param
                                                                 )
                                                                 
                                                                 # test decision
                                                                 p.h1.list <- cbind(p.h1.info, p.h1.unif, p.h1.mix, p.h1.fixed_pp, p.h1.EB)
                                                                 
                                                                 dec.list <- cbind(
                                                                   p.h1.info > (1 - alpha.b),
                                                                   p.h1.unif > (1 - alpha.b),
                                                                   p.h1.mix > (1 - alpha.b),
                                                                   p.h1.fixed_pp > (1 - alpha.b),
                                                                   p.h1.EB > (1 - alpha.b)
                                                                 )
                                                                 
                                                                 
                                                                 AvgOC <- wmatT %*% dec.list
                                                                 AvgOC
                                                               })
                                                               names_fixedpp <- paste0("fixed_pp", "_", p_parameter)
                                                               names_mix <- do.call(paste, c(list("mix"), grid_mixture, sep = "_"))
                                                               rownames(Avg.OC) <- c("info", "freq", names_mix, names_fixedpp, "EB")
                                                               colnames(Avg.OC) <- c("average_TIE", "average_Power")
                                                               
                                                               Avg.OC <- as.data.frame(Avg.OC)
                                                               Avg.OC <- cbind(Avg.OC, mu_ext = mu_ext, nt, nc, n_ext, design_prior_type = design_prior_type,threshold_value,design_prior_mean)
                                                               Avg.OC$analysis_prior_type <- rownames(Avg.OC)
                                                               
                                                               Avg.OC$robust_location <- ifelse(grepl("mix", Avg.OC$analysis_prior_type), sub(".*_([^_]+)_[^_]+$", "\\1", Avg.OC$analysis_prior_type), NA)
                                                               Avg.OC$robust_sd <- as.numeric(ifelse(grepl("mix", Avg.OC$analysis_prior_type), sub("^[^_]+_([^_]+)_.*", "\\1", Avg.OC$analysis_prior_type), NA))
                                                               Avg.OC$w.mix <- as.numeric(ifelse(grepl("mix", Avg.OC$analysis_prior_type), sub(".*_", "", Avg.OC$analysis_prior_type), NA))
                                                               Avg.OC$power_parameter <- as.numeric(ifelse(grepl("fixed", Avg.OC$analysis_prior_type), sub(".*_", "", Avg.OC$analysis_prior_type), NA))
                                                               return(Avg.OC)
                                                             })
                                                           },
                                                           mc.cores = ncores
))

average_OC_recalibrateddecision_subset <- subset(average_OC_recalibrateddecision, analysis_prior_type %in% c("mix_1_ext_0.5", "mix_1_current.mean_0.5", "EB", "info", "freq", "fixed_pp_0.2", "fixed_pp_0.8") & design_prior_mean %in% c(-3, -0.5, 0, 0.5, 3))

average_OC_recalibrateddecision_subset <- subset(average_OC_recalibrateddecision_subset, design_prior_type!="vague")

average_OC_recalibrateddecision_subset$analysis_prior_type <- as.factor(average_OC_recalibrateddecision_subset$analysis_prior_type)
average_OC_recalibrateddecision_subset$line_type <- ifelse(average_OC_recalibrateddecision_subset$analysis_prior_type == "freq", "dashed", "solid")

old_levels <- levels(average_OC_recalibrateddecision_subset$analysis_prior_type)

my_colors <- setNames(
  c(
    "darkorchid", "#377EB8", "#4DAF4A", "darkgoldenrod",
    "red", "azure4", "deeppink4"
  ),
  old_levels
)

levels(average_OC_recalibrateddecision_subset$analysis_prior_type) <- c(
  EB = TeX("$EB \\, power \\, prior"),
  fixed_pp_0.2 = TeX("$\\Fixed \\, power \\, prior \\, (\\alpha=0.2)"),
  fixed_pp_0.8 = TeX("$\\Fixed \\, power \\, prior \\, (\\alpha=0.8)"),
  freq = TeX("$\\No \\, borrowing"),
  info = TeX("$\\Full \\, borrowing"),
  mix_1_current.mean_0.5 = TeX("$\\RMP \\, (\\mu_{robust}=\\bar{y}_c, \\, unit-info, w=0.5)$"),
  mix_1_ext_0.5 = TeX("$\\RMP \\, (\\mu_{robust}=\\bar{y}_{ext}, \\, unit-info, w=0.5)$")
)

# Rename the color vector to match the new levels
names(my_colors) <- levels(average_OC_recalibrateddecision_subset$analysis_prior_type)

average_OC_recalibrateddecision_subset$designprior_title <- as.factor("Design_prior_used")
levels(average_OC_recalibrateddecision_subset$designprior_title) <- c(Design_prior_used = TeX("$\\Design \\, prior\\, used$"))
average_OC_recalibrateddecision_subset$conflict_level_title <- as.factor("conflict_level")
levels(average_OC_recalibrateddecision_subset$conflict_level_title) <- c(conflict_level = TeX("$\\Design \\, prior\\, mean\\, (\\mu^d)$"))

average_OC_recalibrateddecision_subset$design_prior_type <- factor(average_OC_recalibrateddecision_subset$design_prior_type)

levels(average_OC_recalibrateddecision_subset$design_prior_type) <- c(
  fixed_pp = TeX("$\\Fixed \\, power \\, prior"),
  info = TeX("$Informative"),
  mix = TeX("$\\RMP$"),
  unit = TeX("$\\Unit-information$")
)

ggplot() +
  geom_line(data = average_OC_recalibrateddecision_subset, aes(x =average_TIE , y = average_Power, color = as.factor(analysis_prior_type),linetype = line_type), size = 0.8)+
  facet_nested(designprior_title*design_prior_type~conflict_level_title*design_prior_mean,labeller =labeller(.cols = label_parsed,.rows = label_parsed))+
  scale_color_manual(
    values = my_colors,
    labels = scales::label_parse()
  )+
  scale_linetype_manual(
    values = c("solid" = "solid", "dashed" = "dashed"),
    guide = "none"
  )+
  guides(
    color = guide_legend(label.hjust = 0,title.position = "top",title.hjust = 0.5, title = "Analysis prior used",nrow = 4,override.aes = list(
      linetype = c(rep("solid", 3),"dashed", rep("solid", 3))
    ))
  )+
  xlab("average TIE rate")+
  ylab("average Power")+
  theme_bw()+
  theme(axis.text = element_text(color = c("black", "black", "black", "black", "black")),
        axis.ticks = element_line(color = c("black", "black", "black", "black", "black")))+
  theme(axis.text=element_text(size=20,family = "serif"),
        axis.title.x=element_text(size=22,family = "serif"),
        axis.title.y=element_text(size=22,family = "serif"),
        legend.text = element_text(size=25,family = "serif"),
        legend.title=element_text(size=25,family = "serif"),
        legend.position ="bottom",
        legend.box="vertical",
        #legend.key.width =  unit(3.5, "line"),
        legend.key.width =  unit(3, "line"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.spacing.y = unit(0.3, 'cm'),
        strip.text.x = element_text(size=19,family = "serif"),
        strip.text.y = element_text(size=19,family = "serif"),
        plot.title = element_text(
          size = 26,
          hjust = 0.5
        )
  )+
  theme(panel.background = element_rect(fill = 'white', colour =NA))+
  theme(axis.line = element_line(colour = "black",size=1),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())



## ================================================================================================================ ##
## compromise decision: CD-discard and CD-constraint
## ================================================================================================================ ##
#This script computes average OCs with compromise decision

source("Functions_Subm_N_compromise.R")
source("function_avgOC_compromise.R")

#####parameters

outcome <- "normal"
sigma <- 1        # variance of one observation
mupc <- 0         # prior mean control

alpha.up <- 0.075             # max TIE
alpha.low <- 0.01             # min TIE
alpha.b <- 0.025              # frequentist TIE
th0 <- 0                      # null hypothesis border
nmc <- 1e6                     # number of Monte Carlo simulations
seed <- 1234                  # random seed

th1 <- 0.8330874                 # point for Power evaluation
delta <- c(th0, th1)

nc <- 20                         ## current control sample size
nt <- 20                         ## current treatment sample size
n_ext <- 15                      ## historical data sample size
mu_ext <- 0                      ## external data mean


design_prior_type <- c("info", "mix", "unit", "fixed_pp","vague") ## different kinds of design priors

design_prior_mean <- seq(-2.5, 2.5, by = 0.01)

t_value <- c(0.5, 1, 4, 10)


grid <- expand.grid(nt = nt, nc = nc, n_ext = n_ext, mu_ext = mu_ext, design_prior_type = design_prior_type,design_prior_mean=design_prior_mean,t_value=t_value)


ncores <- 80
dat.AvgOC=do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))),
                                  function(params) {
                                    with(params, {
                                      #borrowing parameters
                                      cp=t_value  #steepness
                                      ct=t_value
                                      cf=ct*sqrt(1/nc+1/n_ext) #point of full discard
                                      
                                      
                                      alpha.low.mix.05=0.012433     #min TIE with mix mu_robust=bar(y)
                                      alpha.up.mix.05=0.049854         #max TIE with mix mu_robust=bar(y)
                                      
                                      alpha.low.eb=0.015785   #min TIE with EB pp
                                      alpha.up.eb=0.072952    #max TIE with EB pp
                                      
                                      alpha.low.mix.05_RMPext=0.010082   #min TIE with mix mu_robust=bar(y_ext)
                                      alpha.up.mix.05_RMPext=0.055762  #max TIE with mix mu_robust=bar(y_ext)
                                      
                                      
                                      #prior parameters
                                      sigma <- 1
                                      sig2 <- sigma^2                  # variance of one observation
                                      
                                      sigmapc <- sqrt(sig2 / n_ext)    # SD control arm prior
                                      sigmapt <- sqrt(sig2 / 0)        # SD treat arm prior
                                      priorPars = rbind(c(mupc,mupc),c(1/sqrt(n_ext),Inf))
                                      
                                      sigmapc_downweighted_external <- sqrt(sig2 / (n_ext*0.5))    # downweight external data
                                      p_parameter <- c(0.2,0.5,0.8)
                                      unit_sd <- 1  ## 1 for the case considered here, not always one
                                      
                                      
                                      Avg.OC=sapply(delta, function(x)
                                      {
                                        set.seed(seed)
                                        
                                        #sample control param.
                                        if (design_prior_type == "info") {
                                          thc <- rnorm(nmc, priorPars[1, 1]+design_prior_mean, priorPars[2, 1])
                                        } else if (design_prior_type == "mix") {
                                          thc <- rmix(mixnorm(informative = c(0.5, c(priorPars[1, 1]+design_prior_mean, priorPars[2, 1])), vague = c(1 - 0.5, c(priorPars[1, 1]+design_prior_mean, unit_sd)), sigma = 1), nmc)
                                        } else if (design_prior_type == "unit") {
                                          thc <- rnorm(nmc, priorPars[1, 1]+design_prior_mean, sqrt(sig2))
                                        } else if (design_prior_type == "fixed_pp") {
                                          thc <- rnorm(nmc, priorPars[1, 1]+design_prior_mean, sigmapc_downweighted_external)
                                        }else if (design_prior_type == "vague") {
                                          thc <- rnorm(nmc, priorPars[1, 1]+design_prior_mean, 10^4)
                                        }
                                        
                                        
                                        #sample data
                                        yc=rnorm(nmc,thc,sigma/sqrt(nc))
                                        yt=rnorm(nmc,thc+x,sigma/sqrt(nt))
                                        
                                        #Weight vector, 1 for monte carlo
                                        wmatT=rep(1/nmc,nmc)
                                        
                                        priorPars.v=priorPars
                                        priorPars.v[2,]=c(Inf,Inf)
                                        
                                        p.h1.unif <- pH1(
                                          yc = yc, yt = yt, nc = nc, nt = nt,
                                          priorPars = priorPars.v, w.mix = w.mix, th0 = th0, robust_comp_location = robust_comp_location, dE = mu_ext, power_parameter = power_param
                                        )
                                        
                                        
                                        kappa.f.mix.05=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,
                                                                  priorPars=priorPars,
                                                                  outcome=outcome, mc.results=mc.results,
                                                                  alpha.b=alpha.b, alpha.up=alpha.up.mix.05, alpha.low=alpha.low.mix.05, th0=th0,
                                                                  cf=cf,cp=cp,CD='Discard')
                                        
                                        kappa.f.c.mix.05=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,
                                                                    priorPars=priorPars,
                                                                    outcome=outcome, mc.results=mc.results,
                                                                    alpha.b=alpha.b, alpha.up=alpha.up.mix.05, alpha.low=alpha.low.mix.05, th0=th0,
                                                                    cf=cf,cp=cp,CD='Constrain')
                                        
                                        
                                        kappa.f.mix.05_RMPext=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,
                                                                         priorPars=priorPars,
                                                                         outcome=outcome, mc.results=mc.results,
                                                                         alpha.b=alpha.b, alpha.up=alpha.up.mix.05_RMPext, alpha.low=alpha.low.mix.05_RMPext, th0=th0,
                                                                         cf=cf,cp=cp,CD='Discard')
                                        
                                        kappa.f.c.mix.05_RMPext=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,
                                                                           priorPars=priorPars,
                                                                           outcome=outcome, mc.results=mc.results,
                                                                           alpha.b=alpha.b, alpha.up=alpha.up.mix.05_RMPext, alpha.low=alpha.low.mix.05_RMPext, th0=th0,
                                                                           cf=cf,cp=cp,CD='Constrain')
                                        
                                        kappa.f.eb=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,
                                                              priorPars=priorPars,
                                                              outcome=outcome, mc.results=mc.results,
                                                              alpha.b=alpha.b, alpha.up=alpha.up.eb, alpha.low=alpha.low.eb, th0=th0,
                                                              cf=cf,cp=cp,CD='Discard')
                                        
                                        kappa.f.c.eb=cThreshold(yc=yc,yt=yt,nc=nc,nt=nt,
                                                                priorPars=priorPars,
                                                                outcome=outcome, mc.results=mc.results,
                                                                alpha.b=alpha.b, alpha.up=alpha.up.eb, alpha.low=alpha.low.eb, th0=th0,
                                                                cf=cf,cp=cp,CD='Constrain')
                                        
                                        #test decision
                                        p.h1.list=cbind(p.h1.unif)
                                        
                                        dec.list = cbind(
                                          p.h1.unif>(1-alpha.b),
                                          p.h1.unif>(1-kappa.f.mix.05),
                                          p.h1.unif>(1-kappa.f.c.mix.05),
                                          p.h1.unif>(1-kappa.f.mix.05_RMPext),
                                          p.h1.unif>(1-kappa.f.c.mix.05_RMPext),
                                          p.h1.unif>(1-kappa.f.eb),
                                          p.h1.unif>(1-kappa.f.c.eb)
                                        )
                                        
                                        
                                        AvgOC= wmatT %*% dec.list
                                        AvgOC
                                      })
                                      
                                      rownames(Avg.OC) <- c("freq","CD-Discard.05","CD-Constraint.05","CD-Discard.05_RMPext","CD-Constraint.05_RMPext","CD-Discard.EB","CD-Constraint.EB")
                                      
                                      colnames(Avg.OC) <- c("average_TIE", "average_Power")
                                      
                                      Avg.OC <- as.data.frame(Avg.OC)
                                      Avg.OC <- cbind(Avg.OC, mu_ext = mu_ext, nt, nc, n_ext, design_prior_type = design_prior_type,design_prior_mean,cp,ct,alpha.low.mix.05_RMPext,alpha.up.mix.05_RMPext,alpha.up.eb,alpha.low.eb)
                                      Avg.OC$approach_type <- rownames(Avg.OC)
                                      
                                      return(Avg.OC)
                                      
                                    })
                                  },mc.cores=ncores))

## CD discard
dat.AvgOC <- subset(dat.AvgOC, approach_type != "freq")

dat.AvgOC <- subset(dat.AvgOC, design_prior_type != "vague")

dat.AvgOC$design_prior_type <- factor(dat.AvgOC$design_prior_type)

levels(dat.AvgOC$design_prior_type) <- c(
  fixed_pp = TeX("$\\Fixed \\, power \\, prior"),
  info = TeX("$Informative"),
  mix = TeX("$\\RMP$"),
  unit = TeX("$\\Unit-information$")
)
dat.AvgOC$c_n_p_title <- as.factor("t_and_p")
levels(dat.AvgOC$c_n_p_title) <- c(t_and_p = TeX("$Values \\, of \\, t \\, and \\, p \\, (t=p)"))

dat.AvgOC_discard <- dat.AvgOC[dat.AvgOC$approach_type %in% c("CD-Discard.05", "CD-Discard.05_RMPext", "CD-Discard.EB"), ]
dat.AvgOC_discard$approach_type <- factor(dat.AvgOC_discard$approach_type)

levels(dat.AvgOC_discard$approach_type) <- c(
  "CD-Discard.05" = expression(atop("CD-Discard", "(Calibrated to RMP with " ~ mu[robust] * " = " * bar(y)[c] * ")")),
  "CD-Discard.05_RMPext" = expression(atop("CD-Discard", "(Calibrated to RMP with " ~ mu[robust] * " = " * bar(y)[ext] * ")")),
  "CD-Discard.EB" = expression(atop("CD-Discard", "(Calibrated to EB power prior)"))
)

ggplot() +
  geom_hline(yintercept = 0.025, size = 1.2) +
  geom_line(data = dat.AvgOC_discard, aes(x =design_prior_mean , y = average_TIE, color = as.factor(design_prior_type)), size = 1.2)+
  xlab("")+
  facet_nested(approach_type~c_n_p_title*ct,labeller = labeller(.cols = label_parsed,.rows = label_parsed))+
  scale_color_discrete(labels = scales::label_parse())+
  guides(
    color = guide_legend(label.hjust = 0, title = "Design prior used",nrow = 3)
  )+
  xlab(TeX("$\\Design \\, prior\\, mean\\, (\\mu^d)$"))+
  ylab("average TIE rate")+
  theme_bw()+
  theme(axis.text = element_text(color = c("black", "black", "black", "black", "black")),
        axis.ticks = element_line(color = c("black", "black", "black", "black", "black")))+
  theme(axis.text=element_text(size=28,family = "serif"),
        axis.title.x=element_text(size=28,family = "serif"),
        axis.title.y=element_text(size=30,family = "serif"),
        legend.text = element_text(size=25,family = "serif"),
        legend.title=element_text(size=25,family = "serif"),
        legend.position ="bottom",
        legend.box="vertical",
        legend.key.width =  unit(3, "line"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.spacing.y = unit(0.3, 'cm'),
        strip.text.x = element_text(size=24,family = "serif"),
        strip.text.y = element_text(size=24,family = "serif")
  )+
  theme(panel.background = element_rect(fill = 'white', colour =NA))+
  theme(axis.line = element_line(colour = "black",size=1),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

### CD-Constraint
dat.AvgOC_compromise_constraint <- dat.AvgOC[dat.AvgOC$approach_type %in% c("CD-Constraint.05", "CD-Constraint.05_RMPext", "CD-Constraint.EB"), ]
dat.AvgOC_compromise_constraint$approach_type <- factor(dat.AvgOC_compromise_constraint$approach_type)

levels(dat.AvgOC_compromise_constraint$approach_type) <- c(
  "CD-Constraint.05" = expression(atop("CD-Constraint", "(Calibrated to RMP with " ~ mu[robust] * " = " * bar(y)[c] * ")")),
  "CD-Constraint.05_RMPext" = expression(atop("CD-Constraint", "(Calibrated to RMP with " ~ mu[robust] * " = " * bar(y)[ext] * ")")),
  "CD-Constraint.EB" = expression(atop("CD-Constraint", "(Calibrated to EB power prior)"))
)

ggplot() +
  geom_hline(yintercept = 0.025, size = 1.2) +
  geom_line(data = dat.AvgOC_compromise_constraint[dat.AvgOC_compromise_constraint$cp==4,], aes(x =design_prior_mean , y = average_TIE, color = as.factor(design_prior_type)), size = 1.2)+
  xlab("")+
  facet_nested(.~approach_type,labeller = labeller(.cols = label_parsed))+
  scale_color_discrete(labels = scales::label_parse())+
  guides(
    color = guide_legend(label.hjust = 0, title = "Design prior used",nrow = 3)
  )+
  xlab(TeX("$\\Design \\, prior\\, mean\\, (\\mu^d)$"))+
  ylab("average TIE rate")+
  theme_bw()+
  theme(axis.text = element_text(color = c("black", "black", "black", "black", "black")),
        axis.ticks = element_line(color = c("black", "black", "black", "black", "black")))+
  theme(axis.text=element_text(size=28,family = "serif"),
        axis.title.x=element_text(size=28,family = "serif"),
        axis.title.y=element_text(size=30,family = "serif"),
        legend.text = element_text(size=25,family = "serif"),
        legend.title=element_text(size=25,family = "serif"),
        legend.position ="bottom",
        legend.box="vertical",
        legend.key.width =  unit(3, "line"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.spacing.y = unit(0.3, 'cm'),
        strip.text.x = element_text(size=20,family = "serif"),
        strip.text.y = element_text(size=20,family = "serif")
  )+
  theme(panel.background = element_rect(fill = 'white', colour =NA))+
  theme(axis.line = element_line(colour = "black",size=1),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())



## ================================================================================================================ ##
## average MSE of the treatment effect
## ================================================================================================================ ##
#devtools::install_github("DKFZ-biostats/ESS")
library(ESS)
library(RBesT)


## ===========================================================================================================
## average MSE- hybrid control prior but evaluating MSE of the treatment effect
## ===========================================================================================================

nt <- 20                     # sample size in treatment arm
nc <- 20                     # sample size in control arm
n_ext <- 15                  # historical data sample size
mu_ext <- 0                   ## external data mean

design_prior_type <-  c("info", "mix", "unit","fixed_pp","vague")          ## different kinds of design priors
design_prior_mean <- seq(-2.5, 2.5, by = 0.01)

grid <- expand.grid(nt = nt, nc = nc, n_ext = n_ext, mu_ext = mu_ext, design_prior_type = design_prior_type, design_prior_mean = design_prior_mean)

ncores <- 85
nsims <- 1e6

average_MSE.trteffectE <- do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))),
                                                  function(params) {
                                                    with(params, {
                                                      sigma <- 1
                                                      sig2 <- sigma^2               # variance of one observation
                                                      
                                                      sigmapc <- sqrt(sig2 / n_ext) # SD control arm prior
                                                      sigmapt <- sqrt(sig2 / 0)     # SD treat arm prior
                                                      muc <- 0                  # prior mean control
                                                      mut <- 0.4                # prior mean treatment
                                                      priorPars <- rbind(c(muc, mut), c(sigmapc, sigmapt))
                                                      sigmapc_downweighted <- sqrt(sig2 / (n_ext*0.5))    # downweight external data
                                                      
                                                      alpha.b <- 0.025 # frequentist TIE
                                                      th0 <- 0
                                                      
                                                      bound <- FALSE #
                                                      nmc <- nsims # number of Monte Carlo samples
                                                      
                                                      th1 <- 0.8330874 # point for Power evaluation
                                                      delta <- c(th0, th1)
                                                      
                                                      #p_parameter <- c(0.2, 0.3, 0.5, 0.8) 
                                                      unit_sd <- 1  ## 1 for the case considered here, not always one
                                                      
                                                      grid_mixture <- expand.grid(n_robust = c(1,100), robust_component_location = c("ext","current_observed_mean"), w = c(0.25, 0.5, 0.75))
                                                      
                                                      
                                                      set.seed(123)
                                                      
                                                      if (design_prior_type == "info") {
                                                        thc <- rnorm(nmc, priorPars[1, 1]+design_prior_mean, priorPars[2, 1])
                                                      } else if (design_prior_type == "mix") {
                                                        thc <- rmix(mixnorm(informative = c(0.5, c(priorPars[1, 1]+design_prior_mean, priorPars[2, 1])), vague = c(1 - 0.5, c(priorPars[1, 1]+design_prior_mean, unit_sd)), sigma = 1), nmc)
                                                      } else if (design_prior_type == "unit") {
                                                        thc <- rnorm(nmc, priorPars[1, 1]+design_prior_mean, sqrt(sig2))
                                                      }else if (design_prior_type == "fixed_pp") {
                                                        thc <- rnorm(nmc, priorPars[1, 1]+design_prior_mean, sigmapc_downweighted_external)
                                                      } else if (design_prior_type == "vague") {
                                                        thc <- rnorm(nmc, priorPars[1, 1]+design_prior_mean, 10^4)
                                                      }
                                                      
                                                      # sample data
                                                      yc <- rnorm(nmc, thc, sigma / sqrt(nc))
                                                      yt <- rnorm(nmc, thc, sigma / sqrt(nt))
                                                      
                                                      
                                                      MSE_results_mix <- Map(
                                                        function(w, robust_component_location, n_robust) {
                                                          data.frame(w, robust_component_location, n_robust,
                                                                     avg_MSE =
                                                                       
                                                                       if (robust_component_location == "current_observed_mean") {
                                                                         mean(sapply(1:length(yc), function(z) {
                                                                           mixprior_ctrl <- mixcombine(mixnorm(informative = c(1, c(mu_ext, sigma / sqrt(n_ext))), sigma = sigma), mixnorm(vague = c(1, c(yc[z], n_robust)), sigma = sigma, param = "mn"), weight = c(w, (1 - w)))
                                                                           suppressMessages({
                                                                             posterior.sum.ctrl <- postmix(mixprior_ctrl, m = yc[z], n = nc)
                                                                           })
                                                                           
                                                                           post_mean.ctrl <- (posterior.sum.ctrl["m", 1] * posterior.sum.ctrl["w", 1]) + (posterior.sum.ctrl["m", 2] * posterior.sum.ctrl["w", 2])
                                                                           mixprior_trt <- mixnorm(informative = c(1, c(mut, 10^3)), sigma = sigma)
                                                                           suppressMessages({
                                                                             posterior.sum.trt <- postmix(mixprior_trt, m = yt[z], n = nt)
                                                                           })
                                                                           post_mean.trt <- (posterior.sum.trt["m", 1] * posterior.sum.trt["w", 1])
                                                                           return(((post_mean.trt-post_mean.ctrl) - 0)^2)
                                                                         }))
                                                                       } else if (robust_component_location == "ext") {
                                                                         mean(sapply(1:length(yc), function(z) {
                                                                           mixprior_ctrl <- mixcombine(mixnorm(informative = c(1, c(mu_ext, sigma / sqrt(n_ext))), sigma = sigma), mixnorm(vague = c(1, c(mu_ext, n_robust)), sigma = sigma, param = "mn"), weight = c(w, (1 - w)))
                                                                           suppressMessages({
                                                                             posterior.sum.ctrl <- postmix(mixprior_ctrl, m = yc[z], n = nc)
                                                                           })
                                                                           
                                                                           post_mean.ctrl <- (posterior.sum.ctrl["m", 1] * posterior.sum.ctrl["w", 1]) + (posterior.sum.ctrl["m", 2] * posterior.sum.ctrl["w", 2])
                                                                           mixprior_trt <- mixnorm(informative = c(1, c(mut, 10^3)), sigma = sigma)
                                                                           suppressMessages({
                                                                             posterior.sum.trt <- postmix(mixprior_trt, m = yt[z], n = nt)
                                                                           })
                                                                           post_mean.trt <- (posterior.sum.trt["m", 1] * posterior.sum.trt["w", 1])
                                                                           return(((post_mean.trt-post_mean.ctrl) - 0)^2)
                                                                         }))
                                                                       }, analysis_prior_type = "mix"
                                                          )
                                                        },
                                                        grid_mixture$w, grid_mixture$robust_component_location, grid_mixture$n_robust
                                                      )
                                                      MSE_results_mix <- do.call(rbind, MSE_results_mix)
                                                      
                                                      
                                                      MSE_results_info <- data.frame(avg_MSE = mean(sapply(1:length(yc), function(z) {
                                                        mixprior_ctrl <- mixnorm(informative = c(1, c(mu_ext, sigma / sqrt(n_ext))), sigma = sigma)
                                                        suppressMessages({
                                                          posterior.sum.ctrl <- postmix(mixprior_ctrl, m = yc[z], n = nc)
                                                        })
                                                        
                                                        post_mean.ctrl <- (posterior.sum.ctrl["m", 1] * posterior.sum.ctrl["w", 1]) 
                                                        mixprior_trt <- mixnorm(informative = c(1, c(mut, 10^3)), sigma = sigma)
                                                        suppressMessages({
                                                          posterior.sum.trt <- postmix(mixprior_trt, m = yt[z], n = nt)
                                                        })
                                                        post_mean.trt <- (posterior.sum.trt["m", 1] * posterior.sum.trt["w", 1])
                                                        return(((post_mean.trt-post_mean.ctrl) - 0)^2)
                                                      })), analysis_prior_type = "info")
                                                      
                                                      
                                                      MSE_results_freq <- data.frame(avg_MSE = mean(sapply(1:length(yc), function(z) {
                                                        mixprior_ctrl <- mixnorm(informative = c(1, c(mu_ext, 10^4)), sigma = sigma)
                                                        suppressMessages({
                                                          posterior.sum.ctrl <- postmix(mixprior_ctrl, m = yc[z], n = nc)
                                                        })
                                                        
                                                        post_mean.ctrl <- (posterior.sum.ctrl["m", 1] * posterior.sum.ctrl["w", 1]) 
                                                        mixprior_trt <- mixnorm(informative = c(1, c(mut, 10^3)), sigma = sigma)
                                                        suppressMessages({
                                                          posterior.sum.trt <- postmix(mixprior_trt, m = yt[z], n = nt)
                                                        })
                                                        post_mean.trt <- (posterior.sum.trt["m", 1] * posterior.sum.trt["w", 1])
                                                        return(((post_mean.trt-post_mean.ctrl) - 0)^2)
                                                      })), analysis_prior_type = "freq")
                                                      
                                                      MSE_results_EB <- data.frame(avg_MSE = mean(sapply(1:length(yc), function(z) {
                                                        mixprior_ctrl <- as.powerprior(mixnorm(informative = c(1, c(mu_ext, sigma / sqrt(n_ext))), sigma = sigma))
                                                        suppressMessages({
                                                          posterior.sum.ctrl <- postmix(mixprior_ctrl, m = yc[z], n = nc)
                                                        })
                                                        
                                                        post_mean.ctrl <- (posterior.sum.ctrl["m", 1] * posterior.sum.ctrl["w", 1]) 
                                                        mixprior_trt <- mixnorm(informative = c(1, c(mut, 10^3)), sigma = sigma)
                                                        suppressMessages({
                                                          posterior.sum.trt <- postmix(mixprior_trt, m = yt[z], n = nt)
                                                        })
                                                        post_mean.trt <- (posterior.sum.trt["m", 1] * posterior.sum.trt["w", 1])
                                                        return(((post_mean.trt-post_mean.ctrl) - 0)^2)
                                                      })), analysis_prior_type = "EB")
                                                      
                                                      MSE_results <- dplyr::bind_rows(MSE_results_mix, MSE_results_info, MSE_results_freq,MSE_results_EB)
                                                      
                                                      return(cbind(MSE_results, mu_ext, design_prior_type, nt, nc, n_ext,design_prior_mean))
                                                      
                                                    })
                                                  },
                                                  mc.cores = ncores
))

average_MSE.trteffectE$design_prior_analysis_prior <- as.factor(ifelse(
  grepl("mix", average_MSE.trteffectE$analysis_prior_type), #
  paste0(average_MSE.trteffectE$design_prior_type, "_", average_MSE.trteffectE$analysis_prior_type, "_", average_MSE.trteffectE$w, "_", average_MSE.trteffectE$robust_component_location, "_", average_MSE.trteffectE$n_robust), #
  paste0(average_MSE.trteffectE$design_prior_type, "_", average_MSE.trteffectE$analysis_prior_type) #
))


average_MSE.trteffectE$analysis_prior <- as.factor(ifelse(
  grepl("mix", average_MSE.trteffectE$analysis_prior_type), #
  paste0(average_MSE.trteffectE$analysis_prior_type, "_", average_MSE.trteffectE$w, "_", average_MSE.trteffectE$robust_component_location, "_", average_MSE.trteffectE$n_robust), #
  paste0(average_MSE.trteffectE$analysis_prior_type) #
))

average_MSE.trteffectE_subset <- subset(average_MSE.trteffectE, analysis_prior %in% c("info", "EB", "mix_0.5_current_observed_mean_1", "mix_0.5_ext_1"))

average_MSE.trteffectE_subset$analysis_prior <- factor(average_MSE.trteffectE_subset$analysis_prior)

levels(average_MSE.trteffectE_subset$analysis_prior) <- c(
  EB = TeX("$\\EB \\, power \\, prior"),
  info = TeX("$\\Informative \\, prior"),
  mix_0.5_current_observed_mean_1 = expression(atop("RMP", (mu[robust] == bar(y)[c] ~ ", unit-info, w=0.5"))),
  mix_0.5_ext_1 = expression(atop("RMP", (mu[robust] == bar(y)[ext] ~ ", unit-info, w=0.5")))
)


average_MSE.trteffectE_subset$design_prior_type <- factor(average_MSE.trteffectE_subset$design_prior_type)

levels(average_MSE.trteffectE_subset$design_prior_type) <- c(
  fixed_pp = TeX("$\\Fixed \\, power \\, prior"),
  info = TeX("$Informative"),
  mix = TeX("$\\RMP$"),
  unit = TeX("$\\Unit-information$"),
  vague = TeX("$\\Vague$")
)
average_MSE.trteffectE_subset$designprior <- "Design prior"
average_MSE.trteffectE_subset$analysisprior <- "Analysis prior used"

ggplot() +
  geom_line(data = average_MSE.trteffectE_subset, aes(x = design_prior_mean, y = avg_MSE,color=as.factor(design_prior_type)),size=1.4)+
  facet_nested(.~analysisprior*analysis_prior,scales = "free_y",independent = "y",labeller = labeller(analysis_prior=label_parsed))+
  scale_color_discrete(labels = scales::label_parse())+
  guides(
    color = guide_legend(label.hjust = 0, title = "Design prior used",nrow = 3)
  )+
  xlab(TeX("$\\Design \\, prior\\, mean\\, (\\mu^d)$"))+
  ylab("average MSE")+
  theme_bw()+
  theme(axis.text = element_text(color = c("black", "black", "black", "black", "black")),
        axis.ticks = element_line(color = c("black", "black", "black", "black", "black")))+
  theme(axis.text=element_text(size=28,family = "serif"),
        axis.title.x=element_text(size=28,family = "serif"),
        axis.title.y=element_text(size=30,family = "serif"),
        legend.text = element_text(size=30,family = "serif"),
        legend.title=element_text(size=30,family = "serif"),
        legend.position ="bottom",
        legend.box="vertical",
        legend.key.width =  unit(3, "line"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.spacing.y = unit(0.3, 'cm'),
        strip.text.x = element_text(size=18,family = "serif"),
        strip.text.y = element_text(size=18,family = "serif")
  )+
  theme(panel.background = element_rect(fill = 'white', colour =NA))+
  theme(axis.line = element_line(colour = "black",size=1),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())