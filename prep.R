dat = readRDS('../analysis/data/r_obj/pmb-rif_all-clean.rds')
dat[, effect := round(effect, 3)]

# incomplete data interferes with surface analysis methods
dat = dat[!(cond=='intra' & strain=='ATCC27853' &
  !(c1 %in% c(0, 0.0016, 0.016, 0.048, 0.16, 0.48, 
      4.8, 16, 48, 160, 480, 1600) &
   c2 %in% c(0, 0.0004, 0.003, 0.004, 0.012, 0.03, 
     0.12, 0.3, 1, 1.2, 3, 4, 30, 40, 60, 200)))
]
# arrive at informative and uniform concentration range
dat = dat[(d1 %between% c(0.003, 100) | d1==0) & (d2 %between% c(0.003, 100) | d2==0)]

dat[is.na(effect), effect := -4.5]


# # Response surface analysis
# rs = getRS(dat, 'rs')
# 
# 
# # Time-kill ------------------------------
# 
# dtk = fread('input/dat/time-kill.tsv') %>% 
#   remNonCharcoal()
# 
# .myLabels <- c("neg", "RIF 0.3x MIC", "PMB 1x MIC", "combination", 
#   "RIF 3x MIC", "PMB 10x MIC", "PMB 0.3x MIC", "PMB 3x MIC")
# 
# dtk[, new_time_h := ifelse(time_h == 24, 10, time_h)] 
# xaxs = c(0, 2, 4, 6, 10)
# xlabs = c(0, 2, 4, 6, 24)
# dtk[, group := factor(group, levels=1:8, labels=.myLabels)]
