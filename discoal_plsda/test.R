test_neutral <- discoal_sim(
  mu=0.73*(1.1*10^(-8)),
  recomb_rate=0.873*10^(-8),
  Ne=40600,
  genome_length=50000,
  samplesize=100,
  s = 0,
  discoal_path=discoal_path,
  fix_time = NA,
  seed=c(rand_n[2*11-1], rand_n[2*11]),
  sweep="neutral",
  start_freq = 0,
  popsize_changes = popsize_CHB,
  demes = NA,
  sample_dist = NA,
  deme_join = NULL)

test_sdn <- discoal_sim(
  mu=0.73*(1.1*10^(-8)),
  recomb_rate=0.873*10^(-8),
  Ne=40600,
  genome_length=50000,
  samplesize=100,
  s = 0.001,
  discoal_path=discoal_path,
  fix_time = 0,
  seed=c(rand_n[2*11-1], rand_n[2*11]),
  sweep='hard',
  popsize_changes = popsize_CHB,
  demes = NA,
  sample_dist = NA,
  deme_join = NULL
)

test_ssv <-discoal_sim(
  mu=0.73*(1.1*10^(-8)),
  recomb_rate=0.873*10^(-8),
  Ne=40600,
  genome_length=50000,
  samplesize=100,
  s = 0.001,
  discoal_path=discoal_path,
  fix_time = 0,
  seed=c(rand_n[2*11-1], rand_n[2*11]),
  sweep='soft',
  start_freq = 0.1,
  popsize_changes = popsize_CHB,
  demes = NA,
  sample_dist = NA,
  deme_join = NULL
)

test_sdn$num_seg
test_ssv$num_seg
test_neutral$num_seg