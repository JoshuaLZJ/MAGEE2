#functions for discoal simulations
library(popgen.tools)
library(foreach)

discoal_path <- "/home/zhulaw/Documents/Honours/bin/discoal"
popsize_CHB <- readRDS('popsize_CHB.rds')
rand_n <- readRDS('random_number.rds')

discoal_neutral <- function(nsims=100, rand_n=rand_n){
  all = list()
  foreach (i=1:nsims) %do% {
    all[[i]] <- discoal_sim(
      mu=0.73*(1.1*10^(-8)),
      recomb_rate=0.873*10^(-8),
      Ne=40600,
      genome_length=50000,
      samplesize=100,
      s = 0,
      discoal_path=discoal_path,
      fix_time = NA,
      seed=c(rand_n[2*i-1], rand_n[2*i]),
      sweep="neutral",
      start_freq = 0,
      popsize_changes = popsize_CHB,
      demes = NA,
      sample_dist = NA,
      deme_join = NULL)
  }
  return(all)
}

discoal_ssv <- function(nsims=100, rand_n=rand_n){
  all = list()
  foreach (i=1:nsims) %do% {
    all[[i]] <- discoal_sim(
      mu=0.73*(1.1*10^(-8)),
      recomb_rate=0.873*10^(-8),
      Ne=40600,
      genome_length=50000,
      samplesize=100,
      s = 0.01,
      discoal_path=discoal_path,
      fix_time = 0,
      seed=c(rand_n[2*i-1], rand_n[2*i]),
      sweep='soft',
      start_freq = 0.1,
      popsize_changes = popsize_CHB,
      demes = NA,
      sample_dist = NA,
      deme_join = NULL
    )
  }
  return(all)
}

discoal_sdn <- function(nsims=100, rand_n=rand_n){
  all = list()
  foreach (i=1:nsims) %do% {
    all[[i]] <- discoal_sim(
      mu=0.73*(1.1*10^(-8)),
      recomb_rate=0.873*10^(-8),
      Ne=40600,
      genome_length=50000,
      samplesize=100,
      s = 0.01,
      discoal_path=discoal_path,
      fix_time = 0,
      seed=c(rand_n[2*i-1], rand_n[2*i]),
      sweep='hard',
      popsize_changes = popsize_CHB,
      demes = NA,
      sample_dist = NA,
      deme_join = NULL
    )
  }
  return(all)
}

sum_stats_all_sims <- function(all, nsims=100){
  all_sum_stats = list()
  for (i in c(1:nsims)){
    print(i)
    all_sum_stats[[i]] <- sum_stats(all[[i]], nwins=2, ID=0, snp=1000)
  }
  return (all_sum_stats)
}

sdn <- discoal_sdn(nsims=100, rand_n=rand_n)
ssv <- discoal_ssv(nsims=100, rand_n=rand_n)
neutral <- discoal_neutral(nsims=100, rand_n=rand_n)

sum_stats_sdn <- sum_stats_all_sims(sdn, nsims=100)
sum_stats_ssv <- sum_stats_all_sims(ssv, nsims=100)
sum_stats_neutral <- sum_stats_all_sims(neutral, nsims=100)