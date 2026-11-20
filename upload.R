# read data
rm(list=ls(all=TRUE))
# setwd("~/TWT")
mycsvfile=list.files(pattern = '*.csv')
list2env(lapply(setNames(mycsvfile,make.names(gsub('*.csv$',"",mycsvfile))),
                read.csv,header=T),envir = .GlobalEnv)


P_Pemb <- function(t) {
  return(exp(-Hllogis(
    t * 42/30,
    PFS_Pemb[,4],
    PFS_Pemb[,5]
  )))
}

P_Pacl <- function(t) {
  return(exp(-Hsurvspline2(
    t * 42/30,
    PFS_Pacl[,4],
    PFS_Pacl[,5],
    PFS_Pacl[,6],
    PFS_Pacl[,7],
    knots = knots1
  )))
}

S_Pemb <- function(t) {
  return(exp(-Hllogis(
    t * 42/30,
    OS_Pemb[,4],
    OS_Pemb[,5]
  )))
}

S_Pacl <- function(t) {
  return(exp(-Hsurvspline2(
    t * 42/30,
    OS_Pacl[,4],
    OS_Pacl[,5],
    OS_Pacl[,6],
    OS_Pacl[,7],
    knots = knots2
  )))
}


# setwd("E:\\TTFCEA\\CEA")

# 转换概率
mod_par <- define_parameters(
  natrual_death = get_who_mr(
    age = 61 + (model_time * 42/ 365),
    # country = "USA",
    local = T
  ) %>%
    rescale_prob(p = ., to = 42 / 365),
  #Pemb
  pNTD_Pemb = natrual_death,
  pNTP_Pemb = (P_Pemb(model_time) - P_Pemb(model_time + 1))/(P_Pemb(model_time)),
  pNTN_Pemb = 1 - pNTP_Pemb - pNTD_Pemb,
  pPTD_Pemb = ifelse(((S_Pemb(model_time) - S_Pemb(model_time + 1))/
                      (S_Pemb(model_time) - P_Pemb(model_time)) - natrual_death) > 1 |
                      (S_Pemb(model_time) - P_Pemb(model_time)) < 0, 1,
    (S_Pemb(model_time) - S_Pemb(model_time + 1))/
    (S_Pemb(model_time) - P_Pemb(model_time)) - natrual_death),
  pPTP_Pemb = 1 - pPTD_Pemb,
  pDTD_Pemb = 1,
  #Placebo
  pNTD_Pacl = natrual_death,
  pNTP_Pacl = (P_Pacl(model_time) - P_Pacl(model_time + 1))/P_Pacl(model_time),
  pNTN_Pacl = 1 - pNTP_Pacl - pNTD_Pacl,
  pPTD_Pacl = ifelse(((S_Pacl(model_time) - S_Pacl(model_time + 1))/
                      (S_Pacl(model_time) - P_Pacl(model_time)) - natrual_death) > 1 |
                      (S_Pacl(model_time) - P_Pacl(model_time) < 0), 1,
    (S_Pacl(model_time) - S_Pacl(model_time + 1))/
    (S_Pacl(model_time) - P_Pacl(model_time)) - natrual_death),
  pPTP_Pacl = 1 - pPTD_Pacl,
  pDTD_Pacl = 1
)

#停药概率 不良反应概率
mod_par <- modify(
  mod_par,
  nDiscontinue_Pemb = 136/320,
  nDiscontinue_Pacl = 109/318,

  pDiscontinue_Pemb = rescale_prob(nDiscontinue_Pemb, to = 42/(30*48)),
  pDiscontinue_Pacl = rescale_prob(nDiscontinue_Pacl, to = 42/(30*48)),
  
  nAE_1to2_Pemb = 55/320,
  nAE_3to5_Pemb = 265/320,
  nAE_1to2_Pacl = 93/318,
  nAE_3to5_Pacl = 225/318,
  
  pAE_1to2_Pemb = rescale_prob(nAE_1to2_Pemb, to = 42/(30*48)),
  pAE_3to5_Pemb = rescale_prob(nAE_3to5_Pemb, to = 42/(30*48)),
  pAE_1to2_Pacl = rescale_prob(nAE_1to2_Pacl, to = 42/(30*48)),
  pAE_3to5_Pacl = rescale_prob(nAE_3to5_Pacl, to = 42/(30*48))
)

#其他
mod_par <- modify(
  mod_par,
  
  Discount_rate = 0.03,
  Discount = rescale_discount_rate(Discount_rate, from = 365, to = 42),
  BSA =	1.67,
  Body_weight = 70,
  Ccr = 70,
  
  nAnemia_Pemb = 38/320,
  nPeripheral_neuropathy_Pemb = 17/320,
  nNeutropenia_Pemb = 45/320,
  nNausea_Pemb = 4/320,
  nFatigue_Pemb = 17/320,
  nDiarrhoea_Pemb = 10/320,
  # nNeutrophil_count_decreased_Pemb = 59/320,
  nPneumonitis_Pemb = 2/320,
  nAdrenal_insufficiency_Pemb = 7/320,
  nInfusion_reaction_Pemb = 6/320,

  nAnemia_Pacl = 25/318,
  nPeripheral_neuropathy_Pacl = 19/318,
  nNeutropenia_Pacl = 44/318,
  nNausea_Pacl = 4/318,
  nFatigue_Pacl = 19/318,
  nDiarrhoea_Pacl = 3/318,
  # nNeutrophil_count_decreased_Pacl = 43/318,
  nPneumonitis_Pacl = 3/318,
  # nAdrenal_insufficiency_Pacl = 0,
  nInfusion_reaction_Pacl = 2/318,
  
  pAnemia_Pemb = rescale_prob(nAnemia_Pemb, to = 42/(30*48)),
  pPeripheral_neuropathy_Pemb = rescale_prob(nPeripheral_neuropathy_Pemb, to = 42/(30*48)),
  pNeutropenia_Pemb = rescale_prob(nNeutropenia_Pemb, to = 42/(30*48)),
  pNausea_Pemb = rescale_prob(nNausea_Pemb, to = 42/(30*48)),
  pFatigue_Pemb = rescale_prob(nFatigue_Pemb, to = 42/(30*48)),
  pDiarrhoea_Pemb = rescale_prob(nDiarrhoea_Pemb, to = 42/(30*48)),
  # pNeutrophil_count_decreased_Pemb = rescale_prob(nNeutrophil_count_decreased_Pemb, to = 42/(30*48)),
  pPneumonitis_Pemb = rescale_prob(nPneumonitis_Pemb, to = 42/(30*48)),
  pAdrenal_insufficiency_Pemb = rescale_prob(nAdrenal_insufficiency_Pemb, to = 42/(30*48)),
  pInfusion_reaction_Pemb = rescale_prob(nInfusion_reaction_Pemb, to = 42/(30*48)),
  
  pAnemia_Pacl = rescale_prob(nAnemia_Pacl, to = 42/(30*48)),
  pPeripheral_neuropathy_Pacl = rescale_prob(nPeripheral_neuropathy_Pacl, to = 42/(30*48)),
  pNeutropenia_Pacl = rescale_prob(nNeutropenia_Pacl, to = 42/(30*48)),
  pNausea_Pacl = rescale_prob(nNausea_Pacl, to = 42/(30*48)),
  pFatigue_Pacl = rescale_prob(nFatigue_Pacl, to = 42/(30*48)),
  pDiarrhoea_Pacl = rescale_prob(nDiarrhoea_Pacl, to = 42/(30*48)),
  # pNeutrophil_count_decreased_Pacl = rescale_prob(nNeutrophil_count_decreased_Pacl, to = 42/(30*48)),
  pPneumonitis_Pacl = rescale_prob(nPneumonitis_Pacl, to = 42/(30*48)),
  pInfusion_reaction_Pacl = rescale_prob(nInfusion_reaction_Pacl, to = 42/(30*48)),
  
  cAnemia = 2199.69,
  cPeripheral_neuropathy = 1971.13,
  cNeutropenia = 1387.11,
  cNausea = 4542.63,
  cFatigue = 922.96,
  cDiarrhoea = 4530.92,
  # cNeutrophil_count_decreased = 1387.11,
  cPneumonitis = 7705.55,
  cAdrenal_insufficiency = 4127.04,
  cInfusion_reaction = 15962.51,
  
  cAE_Pemb = pAnemia_Pemb * cAnemia + 
    pPeripheral_neuropathy_Pemb * cPeripheral_neuropathy + pNeutropenia_Pemb * cNeutropenia +
    pNausea_Pemb * cNausea + pFatigue_Pemb * cFatigue + pDiarrhoea_Pemb * cDiarrhoea + 
    pPneumonitis_Pemb * cPneumonitis +
    pAdrenal_insufficiency_Pemb * cAdrenal_insufficiency + pInfusion_reaction_Pemb * cInfusion_reaction,
  cAE_Pacl = pAnemia_Pacl * cAnemia +
    pPeripheral_neuropathy_Pacl * cPeripheral_neuropathy + pNeutropenia_Pacl * cNeutropenia +
    pNausea_Pacl * cNausea + pFatigue_Pacl * cFatigue + pDiarrhoea_Pacl * cDiarrhoea + 
    pPneumonitis_Pacl * cPneumonitis + pInfusion_reaction_Pacl * cInfusion_reaction
)


#后续治疗
mod_par <- modify(
  mod_par,

  nBevacizumab_Pemb = 235/322,

  nBevacizumab_Pacl = 236/321,

  cBevacizumab = 7.3201,
  
  cpembrolizumab = 60,#60.291
  cPaclitaxel  = 0.135,
  
  cAdministration_first = 133.27,
  cAdministration_seq = 66.47,
  cEnd_of_life = 19066.1,
  cLaboratory_tests_and_scans_per_cycle = 285.31,
  cFollow_up_per_cycle = 39.24,
  cBest_supportive_care_per_cycle = 2771.47,

)

#效用值
mod_par <- modify(
  mod_par,
  uPF	= 0.84, 
  uPD = 0.79,
  uAEs_1to2 =	0.01,
  uAEs_3to5 = 0.16,
  
  uPF_cycle	= uPF * 42/365, 
  uPD_cycle = uPD * 42/365,
  uAEs_1to2_cycle =	uAEs_1to2 * 42/365,
  uAEs_3to5_cycle	= uAEs_3to5 * 42/365,
  
  uPF_Pemb = uPF_cycle - uAEs_1to2_cycle * pAE_1to2_Pemb -
    uAEs_3to5_cycle * pAE_3to5_Pemb,
  uPF_Pacl = uPF_cycle - uAEs_1to2_cycle * pAE_1to2_Pacl -
    uAEs_3to5_cycle * pAE_3to5_Pacl,
  uP = uPD_cycle,
  uD = 0
)

#成本

mod_par <- modify(
  mod_par,
  cpembrolizumab_cycle = ifelse(model_time <= 18, cpembrolizumab * 400, 0),
  cPaclitaxel_cycle = (6 * 80 * BSA * cPaclitaxel),
  cBevacizumab_cycle = (3 * 10 * Body_weight * cBevacizumab),
  
  cAdministration_Pemb_cycle =
    6 * cAdministration_first +
    (ifelse(model_time <= 18, 1, 0) + 3 * nBevacizumab_Pemb) *
    cAdministration_seq,
  
  cAdministration_Pacl_cycle =
    6 * cAdministration_first +
    3 * nBevacizumab_Pacl * cAdministration_seq,
  
  cPFS_Pemb = (cpembrolizumab_cycle + cPaclitaxel_cycle + (cBevacizumab_cycle * nBevacizumab_Pemb)) +
    cAE_Pemb + cAdministration_Pemb_cycle +
    cFollow_up_per_cycle +
    cLaboratory_tests_and_scans_per_cycle,
  
  cPFS_Pacl = (cPaclitaxel_cycle + (cBevacizumab_cycle * nBevacizumab_Pacl)) +
    cAE_Pacl + cAdministration_Pacl_cycle +
    cFollow_up_per_cycle +
    cLaboratory_tests_and_scans_per_cycle,
  
  cPD_Pemb = 
    cFollow_up_per_cycle +
    cLaboratory_tests_and_scans_per_cycle +
    cBest_supportive_care_per_cycle,
  
  cPD_Pacl = 
    cFollow_up_per_cycle +
    cLaboratory_tests_and_scans_per_cycle +
    cBest_supportive_care_per_cycle,
    
)



# Markov Model

mat_trans_Pemb <- define_transition(
  state_names=c("PFS", "PD", "Death"),
  pNTN_Pemb, pNTP_Pemb,  pNTD_Pemb,
  0,        pPTP_Pemb,  pPTD_Pemb,
  0,        0,         1
)
# plot(mat_trans_Pemb)
mat_trans_Pacl <- define_transition(
  state_names=c("PFS", "PD", "Death"),
  pNTN_Pacl, pNTP_Pacl,  pNTD_Pacl,
  0,        pPTP_Pacl,  pPTD_Pacl,
  0,        0,         1
)
# plot(mat_trans_Pacl)

### 定义各个状态下的参数值

state_PFS <- define_state(
  cost = dispatch_strategy(
    Pemb = discount(cPFS_Pemb, r = Discount),
    Pacl = discount(cPFS_Pacl, r = Discount)),
  qaly = dispatch_strategy(
    Pemb = discount(uPF_Pemb, r = Discount),
    Pacl = discount(uPF_Pacl, r = Discount)),
  ly = 42/365
)

state_PD <- define_state(
  cost = dispatch_strategy(
    Pemb = discount(cPD_Pemb, r = Discount),
    Pacl = discount(cPD_Pacl, r = Discount)
  ),
  qaly = discount(uP, r = Discount),
  ly = 42/365
)



###定义治疗方案
strat_Pemb <- define_strategy(
  transition = mat_trans_Pemb,
  PFS = state_PFS,
  PD = state_PD,
  Death = state_D)

strat_Pacl <- define_strategy(
  transition = mat_trans_Pacl,
  PFS = state_PFS,
  PD = state_PD,
  Death = state_D)

#load("res_mod.RData")


##运行模型
start <- Sys.time()

res_mod <- run_model(
   parameters = mod_par,
   init = c(1000, 0, 0),
   Pemb = strat_Pemb,
   Pacl = strat_Pacl,
   cost = cost,
   effect = qaly,
   method = "life-table")

end <- Sys.time()
end - start



#成本效果图
# plot(res_mod,type="ce")
### 单因素敏感性分析
ae_ci <- function(events, total) {
  unname(stats::binom.test(events, total, conf.level = 0.95)$conf.int)
}

ci_AE_1to2_Pemb     <- ae_ci(55, 320)
ci_AE_1to2_Pacl     <- ae_ci(93, 318)
ci_AE_3to5_Pemb     <- ae_ci(265, 320)
ci_AE_3to5_Pacl     <- ae_ci(225, 318)

ci_Anemia_Pemb <- ae_ci(38, 320)
ci_Peripheral_neuropathy_Pemb     <- ae_ci(17, 320)
ci_Neutropenia_Pemb     <- ae_ci(45, 320)
ci_Nausea_Pemb     <- ae_ci(4, 320)
ci_Fatigue_Pemb     <- ae_ci(17, 320)
ci_Diarrhoea_Pemb <- ae_ci(10, 320)
ci_Pneumonitis_Pemb     <- ae_ci(2, 320)
ci_Adrenal_insufficiency_Pemb     <- ae_ci(7, 320)
ci_Infusion_reaction_Pemb     <- ae_ci(6, 320)

ci_Anemia_Pacl <- ae_ci(25, 318)
ci_Peripheral_neuropathy_Pacl     <- ae_ci(19, 318)
ci_Neutropenia_Pacl     <- ae_ci(44, 318)
ci_Nausea_Pacl     <- ae_ci(4, 318)
ci_Fatigue_Pacl     <- ae_ci(19, 318)
ci_Diarrhoea_Pacl <- ae_ci(3, 318)
ci_Pneumonitis_Pacl     <- ae_ci(3, 318)
ci_Infusion_reaction_Pacl     <- ae_ci(2, 318)

ci_Bevacizumab_Pemb <- ae_ci(235, 322)
ci_Bevacizumab_Pacl <- ae_ci(236, 321)



def_dsa <- define_dsa(
  # nDiscontinue_Pemb, Low25("nDiscontinue_Pemb"), High25("nDiscontinue_Pemb"),
  # nDiscontinue_Pacl, Low25("nDiscontinue_Pacl"), High25("nDiscontinue_Pacl"),

  # Discount_rate, 0, 0.08,
  BSA, Low25("BSA"), High25("BSA"),
  Body_weight, Low25("Body_weight"), High25("Body_weight"),
  # Ccr, Low25("Ccr"), High25("Ccr"),

  nAE_1to2_Pemb, ci_AE_1to2_Pemb[1], ci_AE_1to2_Pemb[2],
  nAE_3to5_Pemb, ci_AE_3to5_Pemb[1], ci_AE_3to5_Pemb[2],
  nAE_1to2_Pacl, ci_AE_1to2_Pacl[1], ci_AE_1to2_Pacl[2],
  nAE_3to5_Pacl, ci_AE_3to5_Pacl[1], ci_AE_3to5_Pacl[2],

  nPneumonitis_Pemb, ci_Pneumonitis_Pemb[1], ci_Pneumonitis_Pemb[2],
  nAnemia_Pemb, ci_Anemia_Pemb[1], ci_Anemia_Pemb[2],
  nPeripheral_neuropathy_Pemb, ci_Peripheral_neuropathy_Pemb[1], ci_Peripheral_neuropathy_Pemb[2],
  nNeutropenia_Pemb, ci_Neutropenia_Pemb[1], ci_Neutropenia_Pemb[2],
  nNausea_Pemb, ci_Nausea_Pemb[1], ci_Nausea_Pemb[2],
  nFatigue_Pemb, ci_Fatigue_Pemb[1], ci_Fatigue_Pemb[2],
  nDiarrhoea_Pemb, ci_Diarrhoea_Pemb[1], ci_Diarrhoea_Pemb[2],
  nAdrenal_insufficiency_Pemb, ci_Adrenal_insufficiency_Pemb[1], ci_Adrenal_insufficiency_Pemb[2],
  nInfusion_reaction_Pemb, ci_Infusion_reaction_Pemb[1], ci_Infusion_reaction_Pemb[2],
  
  nPneumonitis_Pacl, ci_Pneumonitis_Pacl[1], ci_Pneumonitis_Pacl[2],
  nAnemia_Pacl, ci_Anemia_Pacl[1], ci_Anemia_Pacl[2],
  nPeripheral_neuropathy_Pacl, ci_Peripheral_neuropathy_Pacl[1], ci_Peripheral_neuropathy_Pacl[2],
  nNeutropenia_Pacl, ci_Neutropenia_Pacl[1], ci_Neutropenia_Pacl[2],
  nNausea_Pacl, ci_Nausea_Pacl[1], ci_Nausea_Pacl[2],
  nFatigue_Pacl, ci_Fatigue_Pacl[1], ci_Fatigue_Pacl[2],
  nDiarrhoea_Pacl, ci_Diarrhoea_Pacl[1], ci_Diarrhoea_Pacl[2],
  nInfusion_reaction_Pacl, ci_Infusion_reaction_Pacl[1], ci_Infusion_reaction_Pacl[2],

  cAnemia, Low25("cAnemia"), High25("cAnemia"),
  cPeripheral_neuropathy, Low25("cPeripheral_neuropathy"), High25("cPeripheral_neuropathy"),
  cNeutropenia, Low25("cNeutropenia"), High25("cNeutropenia"),
  cNausea, Low25("cNausea"), High25("cNausea"),
  cFatigue, Low25("cFatigue"), High25("cFatigue"),
  cDiarrhoea, Low25("cDiarrhoea"), High25("cDiarrhoea"),
  cAdrenal_insufficiency, Low25("cAdrenal_insufficiency"), High25("cAdrenal_insufficiency"),
  cPneumonitis, Low25("cPneumonitis"), High25("cPneumonitis"),
  cInfusion_reaction, Low25("cInfusion_reaction"), High25("cInfusion_reaction"),

  nBevacizumab_Pemb, ci_Bevacizumab_Pemb[1], ci_Bevacizumab_Pemb[2],

  nBevacizumab_Pacl, ci_Bevacizumab_Pacl[1], ci_Bevacizumab_Pacl[2],


  cpembrolizumab, Low25("cpembrolizumab"), High25("cpembrolizumab"),
  cPaclitaxel, Low25("cPaclitaxel"), High25("cPaclitaxel"),
  cBevacizumab, Low25("cBevacizumab"), High25("cBevacizumab"),
  # cPembrolizumab, Low25("cPembrolizumab"), High25("cPembrolizumab"),
  # cAtezolizumab, Low25("cAtezolizumab"), High25("cAtezolizumab"),

  cAdministration_first, Low25("cAdministration_first"), High25("cAdministration_first"),
  cAdministration_seq, Low25("cAdministration_seq"), High25("cAdministration_seq"),
  cEnd_of_life, Low25("cEnd_of_life"), High25("cEnd_of_life"),
  cLaboratory_tests_and_scans_per_cycle, Low25("cLaboratory_tests_and_scans_per_cycle"), High25("cLaboratory_tests_and_scans_per_cycle"),
  cFollow_up_per_cycle, Low25("cFollow_up_per_cycle"), High25("cFollow_up_per_cycle"),
  cBest_supportive_care_per_cycle, Low25("cBest_supportive_care_per_cycle"), High25("cBest_supportive_care_per_cycle"),
  # cHospitalization_per_cycle, Low25("cHospitalization_per_cycle"), High25("cHospitalization_per_cycle"),

  uPF, Low10("uPF"), High10("uPF"),
  uPD, Low10("uPD"), High10("uPD"),
  uAEs_1to2, Low10("uAEs_1to2"), High10("uAEs_1to2"),
  uAEs_3to5, Low10("uAEs_3to5"), High10("uAEs_3to5")
)

gc()



# use_cluster(4, cluster = clust)
# use_cluster(1)
res_dsa <- run_dsa(res_mod, dsa = def_dsa)

close_cluster()

end <- Sys.time()
end - start

save(res_dsa, file = "res_dsa.RData")



### 概率敏感性分析

def_psa<-define_psa(
  
  # nDiscontinue_Pemb~beta(48, 320-48),
  # nDiscontinue_Pacl~beta(27, 318-27),

  BSA~normal(mean = get_param_value("BSA"),
             sd = get_normal_sd("BSA")),
  Body_weight~normal(mean = get_param_value("Body_weight"),
                     sd = get_normal_sd("Body_weight")),
  # Ccr~normal(mean = get_param_value("Ccr"),
  #            sd = get_normal_sd("Ccr")),
  
  nAE_1to2_Pemb~beta(55, 320 - 55),
  nAE_3to5_Pemb~beta(265, 320 - 265),
  nAE_1to2_Pacl~beta(93, 318 - 93),
  nAE_3to5_Pacl~beta(225, 318 - 225),

  nPneumonitis_Pemb~beta(2, 320 - 2),
  nAnemia_Pemb~beta(38, 320 - 38),
  nPeripheral_neuropathy_Pemb~beta(17, 320 - 17),
  nNeutropenia_Pemb~beta(45, 320 - 45),
  nNausea_Pemb~beta(4, 320 - 4),
  nFatigue_Pemb~beta(17, 320 - 17),
  nDiarrhoea_Pemb~beta(10, 320 - 10),
  nAdrenal_insufficiency_Pemb~beta(7, 320 - 7),
  nInfusion_reaction_Pemb~beta(6, 320 - 6),
  
  nPneumonitis_Pacl~beta(3, 318 - 3),
  nAnemia_Pacl~beta(25, 318 - 25),
  nPeripheral_neuropathy_Pacl~beta(19, 318 - 19),
  nNeutropenia_Pacl~beta(44, 318 - 44),
  nNausea_Pacl~beta(4, 318 - 4),
  nFatigue_Pacl~beta(19, 318 - 19),
  nDiarrhoea_Pacl~beta(3, 318 - 3),
  nInfusion_reaction_Pacl~beta(2, 318 - 2),
  
  cPneumonitis~gamma(mean = get_param_value("cPneumonitis"),
                   sd = get_gamma_sd("cPneumonitis")),
  cAnemia~gamma(mean = get_param_value("cAnemia"),
                sd = get_gamma_sd("cAnemia")),
  cPeripheral_neuropathy~gamma(mean = get_param_value("cPeripheral_neuropathy"),
                    sd = get_gamma_sd("cPeripheral_neuropathy")),
  cNeutropenia~gamma(mean = get_param_value("cNeutropenia"),
                 sd = get_gamma_sd("cNeutropenia")),
  cNausea~gamma(mean = get_param_value("cNausea"),
                  sd = get_gamma_sd("cNausea")),
  cFatigue~gamma(mean = get_param_value("cFatigue"),
                          sd = get_gamma_sd("cFatigue")),
  cDiarrhoea~gamma(mean = get_param_value("cDiarrhoea"),
                              sd = get_gamma_sd("cDiarrhoea")),
  cAdrenal_insufficiency~gamma(mean = get_param_value("cAdrenal_insufficiency"),
                sd = get_gamma_sd("cAdrenal_insufficiency")),
  cInfusion_reaction~gamma(mean = get_param_value("cInfusion_reaction"),
                           sd = get_gamma_sd("cInfusion_reaction")),
  

  nBevacizumab_Pemb~beta(235, 322 - 235),
  # 
  nBevacizumab_Pacl~beta(236, 321 - 236),
  
  cpembrolizumab~gamma(mean = get_param_value("cpembrolizumab"),
             sd = get_gamma_sd("cpembrolizumab")),
  cPaclitaxel~gamma(mean = get_param_value("cPaclitaxel"),
                      sd = get_gamma_sd("cPaclitaxel")),
  cBevacizumab~gamma(mean = get_param_value("cBevacizumab"),
                   sd = get_gamma_sd("cBevacizumab")),
  
  cAdministration_first~gamma(mean = get_param_value("cAdministration_first"),
                        sd = get_gamma_sd("cAdministration_first")),
  cAdministration_seq~gamma(mean = get_param_value("cAdministration_seq"),
                              sd = get_gamma_sd("cAdministration_seq")),
  cEnd_of_life~gamma(mean = get_param_value("cEnd_of_life"),
                     sd = get_gamma_sd("cEnd_of_life")),
  cLaboratory_tests_and_scans_per_cycle~gamma(mean = get_param_value("cLaboratory_tests_and_scans_per_cycle"),
                                              sd = get_gamma_sd("cLaboratory_tests_and_scans_per_cycle")),
  cFollow_up_per_cycle~gamma(mean = get_param_value("cFollow_up_per_cycle"),
                             sd = get_gamma_sd("cFollow_up_per_cycle")),
  cBest_supportive_care_per_cycle~gamma(mean = get_param_value("cBest_supportive_care_per_cycle"),
                                        sd = get_gamma_sd("cBest_supportive_care_per_cycle")),
  # cHospitalization_per_cycle~gamma(mean = get_param_value("cHospitalization_per_cycle"),
  #                                sd = get_gamma_sd("cHospitalization_per_cycle")),


  uPF~beta(get_alpha("uPF"), get_beta("uPF")),
  uPD~beta(get_alpha("uPD"), get_beta("uPD")),
  uAEs_1to2~beta(get_alpha("uAEs_1to2"), get_beta("uAEs_1to2")),
  uAEs_3to5~beta(get_alpha("uAEs_3to5"), get_beta("uAEs_3to5"))
)
