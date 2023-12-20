which_model <- commandArgs(trailingOnly = TRUE)
which_model <- as.numeric(which_model[1])

library(tidyverse)
library(brms)
library(cmdstanr)
library(loo)
options(mc.cores = 4)

load("/ff_ern_alldata.RData")
set_cmdstan_path("/cmdstan-2.29.2")
save_path <- "/FF_ERN_brms"

n_chains <- 4
n_cores <- 4
n_iter <- 5000
n_warmup <- 15000
n_seed <- 102823
n_threads <- 6


pr_ls <- c(
  set_prior("normal(0,3)", class = "b"),
  set_prior("lkj(3)", class = "L"),
  set_prior("student_t(10, 0, 2)", class = "sd")
)

erp_resp_onlysngtrl$site <- as.factor(ifelse(erp_resp_onlysngtrl$subjid < 2000, "usf", "byu"))
erp_resp_onlysngtrl <- erp_resp_onlysngtrl %>%
  separate(task_event, into = c("task", "event"), sep = "_") %>%
  mutate(task = as.factor(task),
         event = as.factor(event))

df_ern <- erp_resp_onlysngtrl[erp_resp_onlysngtrl$erp_comp == "ern",]
df_pe <- erp_resp_onlysngtrl[erp_resp_onlysngtrl$erp_comp == "pe",]


model_formulas <- list(
  #model 1
  bf(erp ~ 0 + event + (0 + event|p|subjid),
     sigma ~ 0 + event + (0 + event|p|subjid)),
  
  #model 2
  bf(erp ~ 0 + event + event:task + (0 + event|p|subjid),
     sigma ~ 0 + event + event:task + (0 + event|p|subjid)),
  
  #model 3
  bf(erp ~ 0 + event + event:site + (0 + event|p|subjid),
     sigma ~ 0 + event + event:site + (0 + event|p|subjid)),
  
  #model 4
  bf(erp ~ 0 + event + event:task + event:site + event:task:site + (0 + event|p|subjid),
     sigma ~ 0 + event + event:task + event:site + event:task:site +(0 + event|p|subjid))
)

# Function to fit, summarize, and add criterion to a model
ff_brm_fitmodel <- function(formula, model_name, df_in,
                                   pr_ls, n_chains, n_cores,
                                   n_iter, n_warmup, n_seed,
                                   n_threads, save_path) {
  model <- brm(formula,
               data = df_in,
               family = gaussian(),
               prior = pr_ls,
               chains = n_chains,
               cores = n_cores,
               iter = n_iter + n_warmup,
               warmup = n_warmup,
               seed = n_seed,
               sample_prior = "yes",
               backend = "cmdstanr",
               threads = threading(n_threads),
               file = file.path(save_path, model_name))
  
  print(summary(model))
  
  add_criterion(model,
                c("loo"),
                pointwise = TRUE,
                file = file.path(save_path, model_name),
                force_save = TRUE)
  
  return(model)
}

# Use the function based on the value of which_model
if (1 <= which_model && which_model <= 4) {
  model_name <- paste0("brm_ern_m", which_model)
  assign(model_name, 
         ff_brm_fitmodel(model_formulas[[which_model]], model_name, df_ern,
                         pr_ls, n_chains, n_cores,
                         n_iter, n_warmup, n_seed,
                         n_threads, save_path))
} else if (5 <= which_model && which_model <= 8) {
  which_model <- which_model - 4
  model_name <- paste0("brm_pe_m", which_model)
  assign(model_name, 
         ff_brm_fitmodel(model_formulas[[which_model]], model_name, df_pe,
                         pr_ls, n_chains, n_cores,
                         n_iter, n_warmup, n_seed,
                         n_threads, save_path))
} else {
  stop("Invalid value for which_model")
}