inp_args <- commandArgs(trailingOnly = TRUE)
which_model <- as.numeric(inp_args[1])
which_chain <- as.numeric(inp_args[2])

library(tidyverse)
library(brms)
library(cmdstanr)


load("/ff_ern_alldata.RData")
set_cmdstan_path("/cmdstan-2.29.2")
save_path <- "/FF_ERN_brms"

n_chains <- 1
n_cores <- 1
n_iter <- 5000
n_warmup <- 15000
n_seed <- 102823
n_threads <- 24


pr_ls <- c(
  set_prior("normal(0,3)", class = "b"),
  set_prior("lkj(3)", class = "L"),
  set_prior("student_t(10, 0, 2)", class = "sd")
)

erp_resp_onlysngtrl$site <- as.factor(ifelse(erp_resp_onlysngtrl$subjid < 2000, "usf", "byu"))

erp_resp_onlysngtrl$task_event <- as.factor(erp_resp_onlysngtrl$task_event)


df_ern <- erp_resp_onlysngtrl[erp_resp_onlysngtrl$erp_comp == "ern",] %>%
  arrange(subjid,task_event) %>%
  group_by(subjid,task_event) %>%
  mutate(trial = row_number())

df_pe <- erp_resp_onlysngtrl[erp_resp_onlysngtrl$erp_comp == "pe",] %>%
  arrange(subjid,task_event) %>%
  group_by(subjid,task_event) %>%
  mutate(trial = row_number())



model_formulas <- list(
  bf(erp ~ 0 + task_event + (0 + task_event|p|subjid) + (0 + task_event|q|trial),
     sigma ~ 0 + task_event + (0 + task_event|p|subjid) + (0 + task_event|q|trial)),
  
  bf(erp ~ 0 + task_event + task_event:site + (0 + task_event|p|subjid) + (0 + task_event|q|trial),
     sigma ~ 0 + task_event + task_event:site + (0 + task_event|p|subjid) + (0 + task_event|q|trial))
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
  
  return(model)
}

# Use the function based on the value of which_model
if (1 <= which_model && which_model <= 2) {
  model_name <- paste0("brm_ern_icc_m", which_model,"_",which_chain)
  assign(model_name, 
         ff_brm_fitmodel(model_formulas[[which_model]], model_name, df_ern,
                         pr_ls, n_chains, n_cores,
                         n_iter, n_warmup, n_seed + which_chain,
                         n_threads, save_path))
} else if (3 <= which_model && which_model <= 4) {
  which_model <- which_model - 2
  model_name <- paste0("brm_pe_icc_m", which_model,"_",which_chain)
  assign(model_name, 
         ff_brm_fitmodel(model_formulas[[which_model]], model_name, df_pe,
                         pr_ls, n_chains, n_cores,
                         n_iter, n_warmup, n_seed + which_chain,
                         n_threads, save_path))
} else {
  stop("Invalid value for which_model")
}