library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


act_dattest_base = readRDS("data_systematic_review_2018/act_dattest_base")
## Modelling the impact 
stan_Acte <- rstan::stan(file="stat_models/probability_estimates_lq and kq_random_effect_mean.stan", 
                         data=act_dattest_base, ## Data from Sherrard-Smith et al (2018) (1)
                         warmup=1000,
                         control = list(adapt_delta = 0.8,
                                        max_treedepth = 20),
                         iter=2000, chains=4)

sum_dattest_base = readRDS("data_systematic_review_2018/sum_dattest_base")
## Modelling the impact 
stan_Sumi <- rstan::stan(file="stat_models/probability_estimates_lq and kq_random_effect_mean.stan", 
                         data=sum_dattest_base, ## Data from Sherrard-Smith et al (2018) (1)
                         warmup=1000,
                         control = list(adapt_delta = 0.8,
                                        max_treedepth = 20),
                         iter=2000, chains=4)


#library(shinystan)  ## can use this to check the model diagnostics
#launch_shinystan(stan_base)

base_Acte <- extract(stan_Acte)
base_Sumi <- extract(stan_Sumi)

time=seq(1,365,length=365)
mean_valssp_Acte  = 1 / (1 + exp(-mean(base_Acte$alpha1) - mean(base_Acte$alpha2)*time))
mean_valssp_ActeU = 1 / (1 + exp(-quantile(base_Acte$alpha1,0.9) - quantile(base_Acte$alpha2,0.9)*time))
mean_valssp_ActeL = 1 / (1 + exp(-quantile(base_Acte$alpha1,0.1) - quantile(base_Acte$alpha2,0.1)*time))

mean_valsfp_Acte  = (1 / (1 + exp(-mean(base_Acte$beta1) - mean(base_Acte$beta2)*time)))
mean_valsfp_ActeU = (1 / (1 + exp(-quantile(base_Acte$beta1,0.9) - quantile(base_Acte$beta2,0.9)*time)))
mean_valsfp_ActeL = (1 / (1 + exp(-quantile(base_Acte$beta1,0.1) - quantile(base_Acte$beta2,0.1)*time)))

##setting same depreciation as mortality for deterrence
mean_valsdet_Acte  = 1 / (1 + exp(-mean(base_Acte$omega1) - mean(base_Acte$alpha2)*time))
mean_valsdet_ActeU = 1 / (1 + exp(-quantile(base_Acte$omega1,0.9) - quantile(base_Acte$alpha2,0.9)*time))
mean_valsdet_ActeL = 1 / (1 + exp(-quantile(base_Acte$omega1,0.1) - quantile(base_Acte$alpha2,0.1)*time))


mean_valssp_Sumi  = 1 / (1 + exp(-mean(base_Sumi$alpha1) - mean(base_Sumi$alpha2)*time))
mean_valssp_SumiU = 1 / (1 + exp(-quantile(base_Sumi$alpha1,0.9) - quantile(base_Sumi$alpha2,0.9)*time))
mean_valssp_SumiL = 1 / (1 + exp(-quantile(base_Sumi$alpha1,0.1) - quantile(base_Sumi$alpha2,0.1)*time))

mean_valsfp_Sumi  = (1 / (1 + exp(-mean(base_Sumi$beta1) - mean(base_Sumi$beta2)*time)))
mean_valsfp_SumiU = (1 / (1 + exp(-quantile(base_Sumi$beta1,0.9) - quantile(base_Sumi$beta2,0.9)*time)))
mean_valsfp_SumiL = (1 / (1 + exp(-quantile(base_Sumi$beta1,0.1) - quantile(base_Sumi$beta2,0.1)*time)))

##setting same depreciation as mortality for deterrence
mean_valsdet_Sumi  = 1 / (1 + exp(-mean(base_Sumi$omega1) - mean(base_Sumi$alpha2)*time))
mean_valsdet_SumiU = 1 / (1 + exp(-quantile(base_Sumi$omega1,0.9) - quantile(base_Sumi$alpha2,0.9)*time))
mean_valsdet_SumiL = 1 / (1 + exp(-quantile(base_Sumi$omega1,0.1) - quantile(base_Sumi$alpha2,0.1)*time))


feed1 = (1 - mean_valssp_Acte) * mean_valsfp_Acte * (1 - mean_valsdet_Acte)
death1 = mean_valssp_Acte * (1 - mean_valsdet_Acte)
rep1 = (1 - (death1 + feed1)) * (1 - mean_valsdet_Acte)
deter1 = mean_valsdet_Acte  

TOTS = feed1 + rep1 + death1 + deter1

feed1a = (1 - mean_valssp_Sumi) * mean_valsfp_Sumi * (1 - mean_valsdet_Sumi)
death1a = mean_valssp_Sumi * (1 - mean_valsdet_Sumi)
rep1a = (1 - (death1a + feed1a)) * (1 - mean_valsdet_Sumi)
deter1a = mean_valsdet_Sumi  

TOTSa = feed1a + rep1a + death1a + deter1a

first_line = feed1 / TOTS
second_line = (feed1 + rep1) / TOTS
third_line = (feed1 + rep1 + deter1 ) / TOTS

Time = 1:365
Time2 = rev(Time)
minimal = rep(0,length(Time))
maximal = rep(1,length(Time))

plot(rev(first_line) ~ Time2,ylim=c(0,1),yaxt="n",
     ylab="Probable outcome (%)",xlab="Time in days",xaxt="n",
     main="Systematic review",
     cex.axis=1.4,cex.lab=1.4,bty="n",pch="")
axis(2,las=2,at=seq(0,1,0.2),labels = seq(0,100,20),cex.axis = 1.4)
axis(1,at=seq(0,365,120),labels = seq(0,365,120),cex.axis = 1.4)
polygon(c(Time2,rev(Time2)),c(rev(first_line),rev(minimal)),col=adegenet::transp("red",0.4),border = NA)
polygon(c(Time2,rev(Time2)),c(rev(second_line),first_line),col=adegenet::transp("orange",0.6),border = NA)
polygon(c(Time2,rev(Time2)),c(rev(third_line),second_line),col=adegenet::transp("darkgreen",0.4),border = NA)
polygon(c(Time2,rev(Time2)),c(maximal,third_line),col=adegenet::transp("royalblue",0.6),border = NA)


first_line = feed1a / TOTSa
second_line = (feed1a + rep1a) / TOTSa
third_line = (feed1a + rep1a + deter1a ) / TOTSa

Time = 1:365
Time2 = rev(Time)
minimal = rep(0,length(Time))
maximal = rep(1,length(Time))

plot(rev(first_line) ~ Time2,ylim=c(0,1),yaxt="n",
     ylab="Probable outcome (%)",xlab="Time in days",xaxt="n",
     main="Systematic review",
     cex.axis=1.4,cex.lab=1.4,bty="n",pch="")
axis(2,las=2,at=seq(0,1,0.2),labels = seq(0,100,20),cex.axis = 1.4)
axis(1,at=seq(0,365,120),labels = seq(0,365,120),cex.axis = 1.4)
polygon(c(Time2,rev(Time2)),c(rev(first_line),rev(minimal)),col=adegenet::transp("red",0.4),border = NA)
polygon(c(Time2,rev(Time2)),c(rev(second_line),first_line),col=adegenet::transp("orange",0.6),border = NA)
polygon(c(Time2,rev(Time2)),c(rev(third_line),second_line),col=adegenet::transp("darkgreen",0.4),border = NA)
polygon(c(Time2,rev(Time2)),c(maximal,third_line),col=adegenet::transp("royalblue",0.6),border = NA)


NEW_DAT = read.csv("data_cone_bioassay\\raw_data_2016_2017.csv",header=TRUE)
NEW_DAT$total_exposed = NEW_DAT$exposed_mosquito_1 + NEW_DAT$exposed_mosquito_2 + 
  NEW_DAT$exposed_mosquito_3 + NEW_DAT$exposed_mosquito_4 
NEW_DAT$total_mortality_24 = NEW_DAT$mortalidade_24h_1 + NEW_DAT$mortalidade_24h_2 + 
  NEW_DAT$mortalidade_24h_3 + NEW_DAT$mortalidade_24h_4
NEW_DAT$total_mortality_72 = NEW_DAT$mortalidade_72h_1 + NEW_DAT$mortalidade_72h_2 + 
  NEW_DAT$mortalidade_72h_3 + NEW_DAT$mortalidade_72h_4


NEW_DAT2 = data.frame(local = rep("Actellic_2016-2017-location",nrow(NEW_DAT)),
                      test_date = NEW_DAT$test_date,
                      month_official = NEW_DAT$Month,
                      house_type = NEW_DAT$house_type,
                      participation = NEW_DAT$participation,
                      mosquito_species = NEW_DAT$mosquito_species,
                      total_exposed = NEW_DAT$total_exposed,
                      total_mortality_24 = NEW_DAT$total_mortality_24,
                      total_mortality_72 = NEW_DAT$total_mortality_72)



Sumi_new_dat = read.csv("data_cone_bioassay\\raw_data_2018_2019_SumiShield (1).csv",header=TRUE)
Sumi_new_dat$total_exposed = Sumi_new_dat$exposed_mosquito_1 + 
  Sumi_new_dat$exposed_mosquito_2 + Sumi_new_dat$exposed_mosquito_3  
Sumi_new_dat$total_mortality_24 = Sumi_new_dat$mortalidade_24h_1 + 
  Sumi_new_dat$mortalidade_24h_2 + Sumi_new_dat$mortalidade_24h_3
Sumi_new_dat$total_mortality_72 = Sumi_new_dat$mortalidade_72h_1 + 
  Sumi_new_dat$mortalidade_72h_2 + Sumi_new_dat$mortalidade_72h_3
Sumi_new_dat$date_of_test = Sumi_new_dat$test_date


Sumi_new_dat = tidyr::separate(Sumi_new_dat, "test_date", c("Year", "Month", "Day"), sep = "-")

Sumi_new_dat2 = data.frame(local = rep("Sumi_2017-2018-location",nrow(Sumi_new_dat)),
                           test_date = Sumi_new_dat$date_of_test,
                           month_official = Sumi_new_dat$month_test,
                           house_type = Sumi_new_dat$house_type,
                           participation = Sumi_new_dat$participation,
                           mosquito_species = Sumi_new_dat$mosquito_species,
                           total_exposed = Sumi_new_dat$total_exposed,
                           total_mortality_24 = Sumi_new_dat$total_mortality_24,
                           total_mortality_72 = Sumi_new_dat$total_mortality_72)


## Decision is to use the 2016-2017 data for Actellic and then the 2017-2018 for Sumi 72 hours
cone_bios = rbind(Sumi_new_dat2, NEW_DAT2)

Con_bio_Acte_d_t_mud = 
  Con_bio_Acte_n_t_mud = 
  Con_bio_Sumi_d_t_mud = 
  Con_bio_Sumi_n_t_mud = 
  Con_bio_Acte_d_t_cem = 
  Con_bio_Acte_n_t_cem = 
  Con_bio_Sumi_d_t_cem = 
  Con_bio_Sumi_n_t_cem = numeric(length(unique(cone_bios$month_official)))

for(m in 1:length(unique(cone_bios$month_official))){
  Con_bio_Acte_d_t_mud[m] =  sum(cone_bios$total_mortality_24[cone_bios$house_type == "Mud" &
                                                                cone_bios$local == "Actellic_2016-2017-location" & cone_bios$month_official == unique(cone_bios$month_official)[m]])
  Con_bio_Acte_n_t_mud[m] =  sum(cone_bios$total_exposed[cone_bios$house_type == "Mud" &
                                                           cone_bios$local == "Actellic_2016-2017-location" &
                                                           cone_bios$month_official == unique(cone_bios$month_official)[m]])
  Con_bio_Sumi_d_t_mud[m] =  sum(cone_bios$total_mortality_72[cone_bios$house_type == "Mud" &
                                                                cone_bios$local == "Sumi_2017-2018-location" &
                                                                cone_bios$month_official == unique(cone_bios$month_official)[m]])
  Con_bio_Sumi_n_t_mud[m] =  sum(cone_bios$total_exposed[cone_bios$house_type == "Mud" &
                                                           cone_bios$local == "Sumi_2017-2018-location" &
                                                           cone_bios$month_official == unique(cone_bios$month_official)[m]])
  
  Con_bio_Acte_d_t_cem[m] =  sum(cone_bios$total_mortality_24[cone_bios$house_type == "Cement" & cone_bios$local == "Actellic_2016-2017-location" & cone_bios$month_official == unique(cone_bios$month_official)[m]])
  Con_bio_Acte_n_t_cem[m] =  sum(cone_bios$total_exposed[cone_bios$house_type == "Cement" & cone_bios$local == "Actellic_2016-2017-location" & cone_bios$month_official == unique(cone_bios$month_official)[m]])
  Con_bio_Sumi_d_t_cem[m] =  sum(cone_bios$total_mortality_72[cone_bios$house_type == "Cement" & cone_bios$local == "Sumi_2017-2018-location" & cone_bios$month_official == unique(cone_bios$month_official)[m]])
  Con_bio_Sumi_n_t_cem[m] =  sum(cone_bios$total_exposed[cone_bios$house_type == "Cement" &
                                                           cone_bios$local == "Sumi_2017-2018-location" &
                                                           cone_bios$month_official == unique(cone_bios$month_official)[m]])
}


## Add in a line to demonstrate the residual efficacy estimated by Mercy in MOZAMBIQUE

N_data = 12
Con_bio_d_t_mud = Con_bio_Acte_d_t_mud[1:12]
Con_bio_n_t_mud = Con_bio_Acte_n_t_mud[1:12]
Con_bio_d_t_cem = Con_bio_Acte_d_t_cem[1:12]
Con_bio_n_t_cem= Con_bio_Acte_n_t_cem[1:12]
time_sequence = c(1:12)*30

data_list_mud = list(N = N_data, ## number
                     d_t = Con_bio_d_t_mud,## repeat for each chemistry
                     n_t = Con_bio_n_t_mud,
                     time = time_sequence,
                     N_eff = 1, ## eg '2' for 2 wall types
                     eff = rep(1,N_data))##[the number of reps for each group in your data]

data_list_cem = list(N = N_data, ## number
                     d_t = Con_bio_d_t_cem,
                     n_t = Con_bio_n_t_cem,
                     time = time_sequence,
                     N_eff = 1, ## eg '2' for 2 wall types
                     eff = rep(1,N_data))##[the number of reps for each group in your data]


stan_model_mud <- stan(file="stat_models/log_mod.stan", 
                       data=data_list_mud, 
                       warmup=500,
                       control = list(adapt_delta = 0.9,
                                      max_treedepth = 20),
                       iter=1000, chains=4)

stan_model_cem <- stan(file="stat_models/log_mod.stan", 
                       data=data_list_cem, 
                       warmup=500,
                       control = list(adapt_delta = 0.9,
                                      max_treedepth = 20),
                       iter=1000, chains=4)

base_moz1 <- extract(stan_model_mud) ## can use this to extract the model parameter estimates
base_moz2 <- extract(stan_model_cem) ## can use this to extract the model parameter estimates

## plot it against your data!
d_t1 = Con_bio_d_t_mud
n_t1 = Con_bio_n_t_mud
DEAD1 = d_t1/n_t1

d_t2 = Con_bio_d_t_cem
n_t2 = Con_bio_n_t_cem
DEAD2 = d_t2/n_t2

time = seq(1,365,by=1)

mean_prediction_mud = 1 / (1 + exp(-mean(base_moz1$alpha1[,1]) - 
                                     mean(base_moz1$alpha2[,1])*time))
max_prediction_mud = 1 / (1 + exp(-quantile(base_moz1$alpha1[,1],0.9) - 
                                    quantile(base_moz1$alpha2[,1],0.9)*time))
min_prediction_mud = 1 / (1 + exp(-quantile(base_moz1$alpha1[,1],0.1) - 
                                    quantile(base_moz1$alpha2[,1],0.1)*time))

mean_prediction_cem = 1 / (1 + exp(-mean(base_moz2$alpha1[,1]) - 
                                     mean(base_moz2$alpha2[,1])*time))
max_prediction_cem = 1 / (1 + exp(-quantile(base_moz2$alpha1[,1],0.9) - 
                                    quantile(base_moz2$alpha2[,1],0.9)*time))
min_prediction_cem = 1 / (1 + exp(-quantile(base_moz2$alpha1[,1],0.1) - 
                                    quantile(base_moz2$alpha2[,1],0.1)*time))

## Plotting these produces part of Figure 2 main manuscript


## The mean prediction is weighted by the proportion of households
## with mud or cement walls in each village

## percent_mud: 40% for Matutuine, 97% for Boane
## percent_cem: 60% for Matutuine,  3% for Boane

##e.g for Matuttuine
PERCENT_MUD = 0.40
PERCENT_CEM = 0.60

percent_mud = PERCENT_MUD
percent_cem = PERCENT_CEM

mean_prediction = (mean_prediction_mud*percent_mud) + (mean_prediction_cem*percent_cem)
max_prediction = (max_prediction_mud*percent_mud) + (max_prediction_cem*percent_cem)
min_prediction = (min_prediction_mud*percent_mud) + (min_prediction_cem*percent_cem)

feed2 = (1 - mean_prediction) * mean_valsfp_Acte * (1 - mean_valsdet_Acte)
death2 = mean_prediction  * (1 - mean_valsdet_Acte)


feed2l = (1 - min_prediction) * mean_valsfp_ActeL * (1 - mean_valsdet_ActeL)
death2l = min_prediction  * (1 - mean_valsdet_ActeL)
feed2u = (1 - max_prediction) * mean_valsfp_ActeU * (1 - mean_valsdet_ActeU)
death2u = max_prediction  * (1 - mean_valsdet_ActeU)

## This is repeated for Boane and Sumishield data and combined into a data file

N_data = 12
Con_bio_d_t_mud = Con_bio_Sumi_d_t_mud[1:12]
Con_bio_n_t_mud = Con_bio_Sumi_n_t_mud[1:12]
Con_bio_d_t_cem = Con_bio_Sumi_d_t_cem[1:12]
Con_bio_n_t_cem= Con_bio_Sumi_n_t_cem[1:12]
time_sequence = c(1:12)*30

data_list_mud = list(N = N_data, ## number
                     d_t = Con_bio_d_t_mud,## repeat for each chemistry
                     n_t = Con_bio_n_t_mud,
                     time = time_sequence,
                     N_eff = 1, ## eg '2' for 2 wall types
                     eff = rep(1,N_data))##[the number of reps for each group in your data]

data_list_cem = list(N = N_data, ## number
                     d_t = Con_bio_d_t_cem,
                     n_t = Con_bio_n_t_cem,
                     time = time_sequence,
                     N_eff = 1, ## eg '2' for 2 wall types
                     eff = rep(1,N_data))##[the number of reps for each group in your data]


stan_model_mud <- stan(file="stat_models/log_mod.stan", 
                       data=data_list_mud, 
                       warmup=500,
                       control = list(adapt_delta = 0.9,
                                      max_treedepth = 20),
                       iter=1000, chains=4)

stan_model_cem <- stan(file="stat_models/log_mod.stan", 
                       data=data_list_cem, 
                       warmup=500,
                       control = list(adapt_delta = 0.9,
                                      max_treedepth = 20),
                       iter=1000, chains=4)

base_moz1 <- extract(stan_model_mud) ## can use this to extract the model parameter estimates
base_moz2 <- extract(stan_model_cem) ## can use this to extract the model parameter estimates

## plot it against your data!
d_t1 = Con_bio_d_t_mud
n_t1 = Con_bio_n_t_mud
DEAD1 = d_t1/n_t1

d_t2 = Con_bio_d_t_cem
n_t2 = Con_bio_n_t_cem
DEAD2 = d_t2/n_t2

time = seq(1,365,by=1)

mean_prediction_mud = 1 / (1 + exp(-mean(base_moz1$alpha1[,1]) - 
                                     mean(base_moz1$alpha2[,1])*time))
max_prediction_mud = 1 / (1 + exp(-quantile(base_moz1$alpha1[,1],0.9) - 
                                    quantile(base_moz1$alpha2[,1],0.9)*time))
min_prediction_mud = 1 / (1 + exp(-quantile(base_moz1$alpha1[,1],0.1) - 
                                    quantile(base_moz1$alpha2[,1],0.1)*time))

mean_prediction_cem = 1 / (1 + exp(-mean(base_moz2$alpha1[,1]) - 
                                     mean(base_moz2$alpha2[,1])*time))
max_prediction_cem = 1 / (1 + exp(-quantile(base_moz2$alpha1[,1],0.9) - 
                                    quantile(base_moz2$alpha2[,1],0.9)*time))
min_prediction_cem = 1 / (1 + exp(-quantile(base_moz2$alpha1[,1],0.1) - 
                                    quantile(base_moz2$alpha2[,1],0.1)*time))

## Plotting these produces part of Figure 2 main manuscript


## The mean prediction is weighted by the proportion of households
## with mud or cement walls in each village

## percent_mud: 40% for Matutuine, 97% for Boane
## percent_cem: 60% for Matutuine,  3% for Boane

##e.g for Matuttuine
PERCENT_MUD = 0.40
PERCENT_CEM = 0.60

percent_mud = PERCENT_MUD
percent_cem = PERCENT_CEM

mean_prediction = (mean_prediction_mud*percent_mud) + (mean_prediction_cem*percent_cem)
max_prediction = (max_prediction_mud*percent_mud) + (max_prediction_cem*percent_cem)
min_prediction = (min_prediction_mud*percent_mud) + (min_prediction_cem*percent_cem)

feed2a = (1 - mean_prediction) * mean_valsfp_Sumi * (1 - mean_valsdet_Sumi)
death2a = mean_prediction  * (1 - mean_valsdet_Sumi)


feed2al = (1 - min_prediction) * mean_valsfp_SumiL * (1 - mean_valsdet_SumiL)
death2al = min_prediction  * (1 - mean_valsdet_SumiL)
feed2au = (1 - max_prediction) * mean_valsfp_SumiU * (1 - mean_valsdet_SumiU)
death2au = max_prediction  * (1 - mean_valsdet_ActeU)

dta1 = data.frame(
  time = time,
  acte_death2 = death2,
  acte_death2u = death2u,
  acte_death2l = death2l,
  acte_fed2 = feed2,
  acte_fed2l = feed2l,
  acte_fed2u = feed2u,
  
  sumi_death2 = death2a,
  sumi_death2u = death2au,
  sumi_death2l = death2al,
  sumi_fed2 = feed2a,
  sumi_fed2l = feed2al,
  sumi_fed2u = feed2au
)


actellic_details = dta1[,1:7]
sumishield_details = dta1[,c(1,8:13)]

write.csv(actellic_details,"data_generated_analysis_3/actellic_details.csv")
write.csv(sumishield_details,"data_generated_analysis_3/sumishield_details.csv")


##############################
##
## 3 Estimated impact of ITNs, from Churcher et al. 2016
##
###############################
is.pbo = 0 #says whether pbo net (0 = standard, 1= PBO)
species =  1 ##  species parameters are generic as we do not yet have enough info!
metric = 1 #1 = best guess, 2= lower 95% confidence interval 3 upper

#Assay to hut mortality conversion		
alpha1=	array(c(rep(0.63445,3),rep(0.012,3),rep(1.294,3)),c(3,3))
alpha2=	array(c(rep(3.997,3),rep(3.171,3),rep(5.119,3)),c(3,3))

#Benefit of PBO in assay		
beta1=	array(c(rep(3.407,2),2.527,rep(2.666,2),1.528,rep(4.331,2),3.547),c(3,3))
beta2=	array(c(rep(5.88,2),0.891,rep(4.754,2),(0.128),rep(6.956,2),1.882),c(3,3))
beta3=	array(c(rep(0.783,2),0,rep(1.038,2),0,rep(0.543,2),0),c(3,3))

#Deterency from mortality		
delta1=	array(c(rep(0.071,3),rep(0.17,3),rep(0.255,3)),c(3,3))
delta2=	array(c(rep(1.257,3),rep(0.627,3),rep(2.073,3)),c(3,3))
delta3=	array(c(rep(-1.517,3),rep(4.03,3),rep(0.646,3)),c(3,3))

#Success from mortality		
theta1=	array(c(rep(0.025,3),rep(0.007,3),rep(0.034,3)),c(3,3))
theta2=	array(c(rep(3.317,3),rep(2.919,3),rep(4.899,3)),c(3,3))

#Decay in insecticide non-PBO net		
mup=	array(c(rep(-2.36,3),rep(2.948,3),rep(1.821,3)),c(3,3))
rhop=	array(c(rep(-3.05,3),rep(3.762,3),rep(2.322,3)),c(3,3))


kp0=0.699
net_halflife=2.64
#1-0.796 Bioko bradley study
##1-0.11 Kagera West study
#c(0.922,0.455)#

surv_bioassay=0		#measure of resistance 0=no resistance 1=100% survival in discriminating dose bioassay}

#Benefit of PBO in assay		
mort_assay=1-surv_bioassay 
mort_hut_a = alpha1[species,metric] + alpha2[species,metric]*(mort_assay-0.5)			              	#relationship mortality in bioassay -> hut trial, logit scale}
mort_hut   = exp(mort_hut_a)/(1+exp(mort_hut_a))

det_hut_a = delta1[species,metric]+delta2[species,metric]*(mort_hut-0.5)+delta3[species,metric]*(mort_hut-0.5)^2	#relationship hut trial mortality -> deterrence}
det_hut   = ifelse(det_hut_a<0,0,det_hut_a)			                  #censored to stop becoming negative}
suc_hut   = theta1[species,metric] *exp(theta2[species,metric] *(1-mort_hut))				              #relationship hut trial mortality -> success}
rep_hut   = 1-suc_hut-mort_hut

n1n0 = 1-det_hut
kp1  = n1n0*suc_hut
jp1  = n1n0*rep_hut+(1-n1n0)
lp1  = n1n0*mort_hut

r_ITN0  = (1-kp1/kp0)*(jp1/(lp1+jp1))		          	#probability of dying with an encounter with ITN (max)}
d_ITN0  = (1-kp1/kp0)*(lp1/(lp1+jp1))		          	#probability of repeating behaviour (max)}
s_ITN0  = 1-d_ITN0-r_ITN0   

mort_max_a = alpha1[species,metric] + alpha2[species,metric]*(1-0.5)				          #maximum mortality seen in huts, used to adjust}
mort_max   = exp(mort_max_a)/(1+exp(mort_max_a))

mort_min_a = alpha1[species,metric] + alpha2[species,metric]*(0-0.5)				          #maximum mortality seen in huts, used to adjust}
mort_min   = exp(mort_min_a)/(1+exp(mort_min_a))

det_max_a = delta1[species,metric]+delta2[species,metric]*(mort_max-0.5)+delta3[species,metric]*(mort_max-0.5)^2	#relationship hut trial mortality -> deterrence}
det_max   = ifelse(det_max_a<0,0,det_max_a)			                  #censored to stop becoming negative}
suc_max   = theta1[species,metric] *exp(theta2[species,metric] *(1-mort_max))				              #relationship hut trial mortality -> success}
rep_max   = 1-suc_max-mort_max

n1n0_max = 1-det_max
kp1_max  = n1n0_max*suc_max
jp1_max  = n1n0_max*rep_max+(1-n1n0_max)
lp1_max  = n1n0_max*mort_max

r_ITN0_max  = (1-kp1_max/kp0)*(jp1_max/(lp1_max+jp1_max))		          	#probability of dying with an encounter with ITN (max)}
d_ITN0_max  = (1-kp1_max/kp0)*(lp1_max/(lp1_max+jp1_max))		          	#probability of repeating behaviour (max)}
s_ITN0_max  = 1-d_ITN0_max-r_ITN0_max


det_min_a = delta1[species,metric]+delta2[species,metric]*(mort_min-0.5)+delta3[species,metric]*(mort_min-0.5)^2	#relationship hut trial mortality -> deterrence}
det_min   = ifelse(det_min_a<0,0,det_min_a)			                  #censored to stop becoming negative}
suc_min   = theta1[species,metric] *exp(theta2[species,metric] *(1-mort_min))				              #relationship hut trial mortality -> success}
rep_min   = 1-suc_min-mort_min

n1n0_min= 1-det_min
kp1_min  = n1n0_min*suc_min
jp1_min  = n1n0_min*rep_min+(1-n1n0_min)
lp1_min  = n1n0_min*mort_min

r_ITN0_min  = (1-kp1_min/kp0)*(jp1_min/(lp1_min+jp1_min))		          	#probability of dying with an encounter with ITN (max)}
d_ITN0_min  = (1-kp1_min/kp0)*(lp1_min/(lp1_min+jp1_min))		          	#probability of repeating behaviour (max)}
s_ITN0_min  = 1-d_ITN0_min-r_ITN0_min

#{halflife}
my_max_washes_a = mup[species,metric] +rhop[species,metric]*(mort_max-0.5)		
my_max_washes   = log(2)/(exp(my_max_washes_a)/(1+exp(my_max_washes_a)))

wash_decay_rate_a = mup[species,metric] +rhop[species,metric]*(mort_hut-0.5)
wash_decay_rate   = log(2)/(exp(wash_decay_rate_a)/(1+exp(wash_decay_rate_a)))
itn_half_life     = wash_decay_rate/my_max_washes*net_halflife

## adjusted to match Griffin et al 2015 Natt Comms
Griff_d_ITN0<-0.51
Griff_r_ITN0<-0.31  ###THINK THIS NEEDS TO BE CHECKED
Griff_s_ITN0<-1-Griff_d_ITN0-Griff_r_ITN0

##ERG parameterisations
##mortality parameters modified to match Jamies paper 
##success paramater scaled to start at jamies paper values and go to elife parameters
ERG_d_ITN0 <- d_ITN0/d_ITN0_max*Griff_d_ITN0
ERG_s_ITN0 <- (Griff_s_ITN0)+(s_ITN0-s_ITN0_max)/(s_ITN0_min-s_ITN0_max)*(s_ITN0_min-Griff_s_ITN0)
ERG_r_ITN0 <- 1-ERG_d_ITN0-ERG_s_ITN0

ERG_r_ITN0;ERG_d_ITN0;itn_half_life

itn_loss = log(2)/itn_half_life
ITN_interval=3*365
time=1:2000
## decay in efficacy of net over time
## **** this is wrong need to work this out
ITN_decay = exp(-(time/ITN_interval)*itn_loss)

r_ITN_min=0.24 
d_ITN = ERG_d_ITN0 * ITN_decay 	 		## insecticide mortality rate 
r_ITN = r_ITN_min + (ERG_r_ITN0 - r_ITN_min)*ITN_decay 
s_ITN = 1 - d_ITN - r_ITN			## successful protected human biting 

# d_ITN;r_ITN;s_ITN

par(mfrow = c(2,3))


w_Acte = yy_Acte = z_Acte = w_Sumi = yy_Sumi = z_Sumi = array(dim=c(180,4,9)) 
## row = time series
## col 1 = no intervention
## col 2 = nets only
## col 3 = spray only
## col 4 = both nets and spray
## dimension 3 stores the uncertainty from phi parameter and efficacy estimates

#############################
## 

##sensitivity analysis
##Species are different in each location 

# use PNAS paper to estimate some range for
# the values of phi B and I (ref 6 below)

PHI_B_mut = c(0.85,0.6,0.95) ## probability of bites in bed (mean, lower and upper)
PHI_I_mut = c(0.90,0.68,0.99) ## probability of bites indoors

PHI_B_boa = c(0.85,0.6,0.95)  ## probability of bites in bed
PHI_I_boa = c(0.90,0.68,0.99) ## probability of bites indoors

## From Griffin et al 2010 (4)
k0 = 0.699 ## expected feeding in the absence of interventions
ksA = actellic_details[,5]
lsA = actellic_details[,2]
jsA = 1 - actellic_details[,2] - actellic_details[,5]

s_IRS_Acte = ksA/k0 ##feed2
r_IRS_Acte = (1 - ksA/k0)*(jsA/(lsA+jsA)) ##rep2

ksS = sumishield_details[,5]
lsS = sumishield_details[,2]
jsS = 1 - sumishield_details[,2] - sumishield_details[,5]

s_IRS_Sumi = ksS/k0 ##feed2
r_IRS_Sumi = (1 - ksS/k0)*(jsS/(lsS+jsS)) ##rep2

## lower
ksAl = actellic_details[,6]
lsAl = actellic_details[,3]
jsAl = 1 - actellic_details[,6] - actellic_details[,3]

s_IRS_Actel = ksAl/k0 ##feed2
r_IRS_Actel = (1 - ksAl/k0)*(jsAl/(lsAl+jsAl)) ##rep2

ksSl = sumishield_details[,6]
lsSl = sumishield_details[,3]
jsSl = 1 - sumishield_details[,6] - sumishield_details[,3]

s_IRS_Sumil = ksSl/k0 ##feed2
r_IRS_Sumil = (1 - ksSl/k0)*(jsSl/(lsSl+jsSl)) ##rep2

## upper
ksAu = actellic_details[,7]
lsAu = actellic_details[,4]
jsAu = 1 - actellic_details[,7] - actellic_details[,4]

s_IRS_Acteu = ksAu/k0 ##feed2
r_IRS_Acteu = (1 - ksAu/k0)*(jsAu/(lsAu+jsAu)) ##rep2

ksSu = sumishield_details[,7]
lsSu = sumishield_details[,4]
jsSu = 1 - sumishield_details[,7] - sumishield_details[,4]

s_IRS_Sumiu = ksSu/k0 ##feed2
r_IRS_Sumiu = (1 - ksSu/k0)*(jsSu/(lsSu+jsSu)) ##rep2

##
w_Acte[,1,] = w_Sumi[,1,] = rep(1,180) 
## Probability that a mosquito bites and survives in the presence of indoor vector control

cl = c(0,3,6)
for(j in 1:3){
  for(i in 1:180){
    PHI_B = PHI_B_mut[j]
    PHI_I = PHI_I_mut[j]
    w_Acte[i,2,1+cl[j]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 
    ## probability of surviving biting given that there is ITN
    w_Acte[i,3,1+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acte[i])*s_IRS_Acte[i]	
    ##			probability of surviving biting given that there is IRS
    w_Acte[i,4,1+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acte[i])*s_ITN[i+547]*s_IRS_Acte[i] + (PHI_I - PHI_B)*(1-r_IRS_Acte[i])*s_IRS_Acte[i] ## probability of surviving biting given that there is ITN & IRS
    
    w_Acte[i,2,2+cl[j]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 
    ## probability of surviving biting given that there is ITN
    w_Acte[i,3,2+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Actel[i])*s_IRS_Actel[i]	
    ##			probability of surviving biting given that there is IRS
    w_Acte[i,4,2+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Actel[i])*s_ITN[i+547]*s_IRS_Actel[i] + (PHI_I - PHI_B)*(1-r_IRS_Actel[i])*s_IRS_Actel[i] ## probability of surviving biting given that there is ITN & IRS
    
    w_Acte[i,2,3+cl[j]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 
    ## probability of surviving biting given that there is ITN
    w_Acte[i,3,3+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acteu[i])*s_IRS_Acteu[i]	
    ##			probability of surviving biting given that there is IRS
    w_Acte[i,4,3+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acteu[i])*s_ITN[i+547]*s_IRS_Acteu[i] + (PHI_I - PHI_B)*(1-r_IRS_Acteu[i])*s_IRS_Acteu[i] ## probability of surviving biting given that there is ITN & IRS
  }  
}

cl = c(0,3,6) ## to capture the ranging uncertainty estimates
for(j in 1:3){
  for(i in 1:180){
    PHI_B = PHI_B_mut[j]
    PHI_I = PHI_I_mut[j]
    
    w_Sumi[i,2,1+cl[j]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 
    ## probability of surviving biting given that there is ITN
    w_Sumi[i,3,1+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumi[i])*s_IRS_Sumi[i]	
    ##			probability of surviving biting given that there is IRS
    w_Sumi[i,4,1+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumi[i])*s_ITN[i+547]*s_IRS_Sumi[i] + 
      (PHI_I - PHI_B)*(1-r_IRS_Sumi[i])*s_IRS_Sumi[i] 
    ## probability of surviving biting given that there is ITN & IRS
    
    w_Sumi[i,2,2+cl[j]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 
    ## probability of surviving biting given that there is ITN
    w_Sumi[i,3,2+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumil[i])*s_IRS_Sumil[i]	
    ##			probability of surviving biting given that there is IRS
    w_Sumi[i,4,2+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumil[i])*s_ITN[i+547]*s_IRS_Sumil[i] + 
      (PHI_I - PHI_B)*(1-r_IRS_Sumil[i])*s_IRS_Sumil[i] 
    ## probability of surviving biting given that there is ITN & IRS
    
    w_Sumi[i,2,3+cl[j]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 
    ## probability of surviving biting given that there is ITN
    w_Sumi[i,3,3+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumiu[i])*s_IRS_Sumiu[i]	
    ##			probability of surviving biting given that there is IRS
    w_Sumi[i,4,3+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumiu[i])*s_ITN[i+547]*s_IRS_Sumiu[i] + 
      (PHI_I - PHI_B)*(1-r_IRS_Sumiu[i])*s_IRS_Sumiu[i] 
    ## probability of surviving biting given that there is ITN & IRS
  }
}


## Probability of any bite (if there is IRS, a mosquito may bite and then die immediately afterwards)
yy_Acte[,1,] = w_Acte[,1,] 
yy_Acte[,2,] = w_Acte[,2,]

yy_Sumi[,1,] = w_Sumi[,1,] 
yy_Sumi[,2,] = w_Sumi[,2,]

for(j in 1:3){
  for(i in 1:180){
    PHI_B = PHI_B_mut[j]
    PHI_I = PHI_I_mut[j]
    yy_Acte[i,3,1+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acte[i])
    yy_Acte[i,4,1+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acte[i])*s_ITN[i+547] + 
      (PHI_I - PHI_B)*(1-r_IRS_Acte[i])
    
    yy_Acte[i,3,2+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Actel[i])
    yy_Acte[i,4,2+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Actel[i])*s_ITN[i+547] + 
      (PHI_I - PHI_B)*(1-r_IRS_Actel[i])
    
    yy_Acte[i,3,3+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acteu[i])
    yy_Acte[i,4,3+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acteu[i])*s_ITN[i+547] + 
      (PHI_I - PHI_B)*(1-r_IRS_Acteu[i])
  }
}

for(j in 1:3){
  for(i in 1:180){
    PHI_B = PHI_B_boa[j]
    PHI_I = PHI_I_boa[j]
    yy_Sumi[i,3,1+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumi[i])
    yy_Sumi[i,4,1+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumi[i])*s_ITN[i+547] + 
      (PHI_I - PHI_B)*(1-r_IRS_Sumi[i])
    
    yy_Sumi[i,3,2+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumil[i])
    yy_Sumi[i,4,2+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumil[i])*s_ITN[i+547] + 
      (PHI_I - PHI_B)*(1-r_IRS_Sumil[i])
    
    yy_Sumi[i,3,3+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumiu[i])
    yy_Sumi[i,4,3+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumiu[i])*s_ITN[i+547] + 
      (PHI_I - PHI_B)*(1-r_IRS_Sumiu[i])
  }  
}

## 
plot(yy_Acte[,1,1] ~ time[1:180],ylim=c(0,1),pch="",
     ylab = "Probability mosquito bites",
     col="black",
     main = "",cex.main=1.2,xlim=c(-60,180),xaxt="n",
     xlab="Time in months",yaxt="n",cex.lab=1.4,cex.axis=1.4,cex=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,1,0.2),cex.lab=1.4,cex.axis=1.4)
axis(1,at=seq(-60,150,30)+15,labels = c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),
     cex.axis = 1.4)

colsd=c("darkred","red","orange","blue")
for(i in 3){
  lines(yy_Acte[,i,1] ~ time[1:180],col="darkblue",lwd=2)
  lines(yy_Sumi[,i,1] ~ time[1:180],col="aquamarine3",lwd=2)
}


## Including max uncertainty for phiI, phiB and efficacy
polygon(c(time[1:180],rev(time[1:180])),
        c(yy_Acte[,3,5],rev(yy_Acte[,3,9])),
        col=adegenet::transp("darkblue",0.3),border=NA)
polygon(c(time[1:180],rev(time[1:180])),
        c(yy_Sumi[,3,5],rev(yy_Sumi[,3,9])),
        col=adegenet::transp("aquamarine3",0.3),border=NA)

legend("topleft",legend = c("Actellic 300CS","SumiShield"),
       col = c("darkblue","aquamarine3"),lwd = 2, lty=c(1,1),cex=1.2,bty="n")


## 
plot(w_Acte[,1,1] ~ time[1:180],ylim=c(0,1),pch="",
     ylab = "Probability mosquito bites and survives",
     col="black",
     main = "",cex.main=1.2,xlim=c(-60,180),xaxt="n",
     xlab="Time in months",yaxt="n",cex.lab=1.4,cex.axis=1.4,cex=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,1,0.2),cex.lab=1.4,cex.axis=1.4)
axis(1,at=seq(-60,150,30)+15,labels = c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),
     cex.axis = 1.4)

colsd=c("darkred","red","orange","blue")
for(i in 3){
  lines(w_Acte[,i,1] ~ time[1:180],col="darkblue",lwd = 2)
  lines(w_Sumi[,i,1] ~ time[1:180],col="aquamarine3",lwd = 2)
}
polygon(c(time[1:180],rev(time[1:180])),
        c(w_Acte[,3,5],rev(w_Acte[,3,9])),
        col=adegenet::transp("darkblue",0.3),border=NA)
polygon(c(time[1:180],rev(time[1:180])),
        c(w_Sumi[,3,5],rev(w_Sumi[,3,9])),
        col=adegenet::transp("aquamarine3",0.3),border=NA)


## Probability repelled
z_Acte[,1,1] = 0
z_Sumi[,1,1] = 0

for(j in 1:3){
  for(i in 1:180){
    PHI_B = PHI_B_boa[j]
    PHI_I = PHI_I_boa[j]
    
    z_Acte[i,2,1+cl[j]] = PHI_B*r_ITN[i+547]
    z_Acte[i,3,1+cl[j]] = PHI_I*r_IRS_Acte[i]
    z_Acte[i,4,1+cl[j]] = PHI_B*(r_IRS_Acte[i] + (1-r_IRS_Acte[i])*r_ITN[i+547]) + (PHI_I - PHI_B)*r_IRS_Acte[i]
    
    z_Acte[i,2,2+cl[j]] = PHI_B*r_ITN[i+547]
    z_Acte[i,3,2+cl[j]] = PHI_I*r_IRS_Actel[i]
    z_Acte[i,4,2+cl[j]] = PHI_B*(r_IRS_Actel[i] + (1-r_IRS_Actel[i])*r_ITN[i+547]) + 
      (PHI_I - PHI_B)*r_IRS_Actel[i]
    
    z_Acte[i,2,3+cl[j]] = PHI_B*r_ITN[i+547]
    z_Acte[i,3,3+cl[j]] = PHI_I*r_IRS_Acteu[i]
    z_Acte[i,4,3+cl[j]] = PHI_B*(r_IRS_Acteu[i] + (1-r_IRS_Acteu[i])*r_ITN[i+547]) + 
      (PHI_I - PHI_B)*r_IRS_Acteu[i]
    
    z_Sumi[i,2,1+cl[j]] = PHI_B*r_ITN[i+547]
    z_Sumi[i,3,1+cl[j]] = PHI_I*r_IRS_Sumi[i]
    z_Sumi[i,4,1+cl[j]] = PHI_B*(r_IRS_Sumi[i] + (1-r_IRS_Sumi[i])*r_ITN[i+547]) + 
      (PHI_I - PHI_B)*r_IRS_Sumi[i]
    
    z_Sumi[i,2,2+cl[j]] = PHI_B*r_ITN[i+547]
    z_Sumi[i,3,2+cl[j]] = PHI_I*r_IRS_Sumil[i]
    z_Sumi[i,4,2+cl[j]] = PHI_B*(r_IRS_Sumil[i] + (1-r_IRS_Sumil[i])*r_ITN[i+547]) + 
      (PHI_I - PHI_B)*r_IRS_Sumil[i]
    
    z_Sumi[i,2,3+cl[j]] = PHI_B*r_ITN[i+547]
    z_Sumi[i,3,3+cl[j]] = PHI_I*r_IRS_Sumiu[i]
    z_Sumi[i,4,3+cl[j]] = PHI_B*(r_IRS_Sumiu[i] + (1-r_IRS_Sumiu[i])*r_ITN[i+547]) + 
      (PHI_I - PHI_B)*r_IRS_Sumiu[i]
  }
}
## 
plot(z_Acte[,1,1] ~ time[1:180],ylim=c(0,1),pch="",
     ylab = "Probability mosquito is repelled",
     col="black",
     main = "",cex.main=1.2,xlim=c(-60,180),xaxt="n",
     xlab="Time in months",yaxt="n",cex.lab=1.4,cex.axis=1.4,cex=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,1,0.2),cex.lab=1.4,cex.axis=1.4)
axis(1,at=seq(-60,150,30)+15,labels = c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),
     cex.axis = 1.4)

colsd=c("darkred","red","orange","blue")
for(i in 3){
  lines(z_Acte[,i,1] ~ time[1:180],col="darkblue",lwd=2)
  lines(z_Sumi[,i,1] ~ time[1:180],col="aquamarine3",lwd=2)
}

polygon(c(time[1:180],rev(time[1:180])),
        c(z_Acte[,3,5],rev(z_Acte[,3,9])),
        col=adegenet::transp("darkblue",0.3),border=NA)
polygon(c(time[1:180],rev(time[1:180])),
        c(z_Sumi[,3,5],rev(z_Sumi[,3,9])),
        col=adegenet::transp("aquamarine3",0.3),border=NA)

## Campaigns tend to take up to a few months to complete
## We assume the ratio of houses monitored per start month
## reflects the proportion of houses covered by the spray campaign
## in that month

## Houses were tracked from November, December and January (1 in Feb)
## 129/(129+88+27+1) is 52.65% of houses in the community have max protection in Nov
## 88/(129+88+27+1) is 35.92% of houses max protection Dec (52% - any modifications have 1 mont old protection)
## 27/(129+88+27+1) is 11.02% max protection Jan (35.92% - modif 1 month old, 52% - modif 2 month old)


## derived ITN/IRS quantities
## prob bites and survives

w_Acte1 = yy_Acte1 = z_Acte1 = w_Sumi1 = yy_Sumi1 = z_Sumi1 = array(dim=c(365,18,4,9))
## row = time series
## cols indicate the ranging initiation of the spray campaign to be combined
## dimension 3: 1 will be the effect if there is no intervention
## dimension 3: is with ITNs only
## dimension 3: is with IRS only no loss in coverage
## dimension 3: is with IRS only loss in coverage
## dimensions 4: the uncertainty for efficacy and phi I, B

#############################
## 

## Data from spray campaign
prop_houses_sprayed_WeeklyB = 0.97*c(0,	0.027219794,	## August
                                     0.077014558,	0.136261919,	0.196901742,	0.250817066, ## Sept
                                     0.301047746,	0.347687015,	0.395348206,	0.464541108, ## Oct
                                     0.521818166,	0.581372602,	0.643327061,	0.710948828, ## Nov
                                     0.777620537,	0.847130239,	0.911498024,	0.930365931, ## Dec
                                     0.947042313,	0.96389906,	0.983400068,	0.991357264,	 ## Jan
                                     0.99238783,	1) ## Feb
prop_houses_sprayed_WeeklyM = 0.96*c(0.065358837,	0.193716785,## August
                                     0.334444596,	0.432141444,	0.533708885,	0.614060416,##sep
                                     0.667531537,	0.711462198,	0.769981823,	0.880751007,##oct
                                     0.919045204,	0.944027976,	0.97443187,	0.990922736,  ##nov
                                     0.991873307,	0.99744169,	 1,             1,          ##dec
                                     1,        1) ## jan

ksA = lsA = jsA = array(dim=c(365,18,9))
for(w in 1:17){
  ksA[,1,1] =  actellic_details[,5][1:365]
  ksA[,w+1,1] = c(rep(k0,w*7),actellic_details[,5][1:(365-7*w)])
  
  ksA[,1,2] =  actellic_details[,6][1:365]
  ksA[,w+1,2] = c(rep(k0,w*7),actellic_details[,6][1:(365-7*w)])
  
  ksA[,1,3] =  actellic_details[,7][1:365]
  ksA[,w+1,3] = c(rep(k0,w*7),actellic_details[,7][1:(365-7*w)])
  
  
  lsA[,1,1] = actellic_details[,2]
  lsA[,w+1,1] = c(rep(0,w*7),actellic_details[,2][1:(365-7*w)])
  
  lsA[,1,2] = actellic_details[,3]
  lsA[,w+1,2] = c(rep(0,w*7),actellic_details[,3][1:(365-7*w)])
  
  lsA[,1,3] = actellic_details[,4]
  lsA[,w+1,3] = c(rep(0,w*7),actellic_details[,4][1:(365-7*w)])
  
  
  jsA[,w,1] = 1 - ksA[,w,1] - lsA[,w,1]
  jsA[,w,2] = 1 - ksA[,w,2] - lsA[,w,2]
  jsA[,w,3] = 1 - ksA[,w,3] - lsA[,w,3]
  
}
ksA[,,4:6] = ksA[,,1:3]
ksA[,,7:9] = ksA[,,1:3]

lsA[,,4:6] = lsA[,,1:3]
lsA[,,7:9] = lsA[,,1:3]

jsA[,,4:6] = jsA[,,1:3]
jsA[,,7:9] = jsA[,,1:3]

jsA[,18,] = 1 - ksA[,18,] - lsA[,18,]



s_IRS_Acte1 = r_IRS_Acte1 = array(dim=c(365,18,9))
for(j in 1:9){
  for(w in 1:18){
    s_IRS_Acte1[,w,j] = ksA[,w,j]/k0 ##feed2 = when IRS is implemented in month 1 (Nov)
    r_IRS_Acte1[,w,j] = (1 - ksA[,w,j]/k0)*(jsA[,w,j]/(lsA[,w,j]+jsA[,w,j])) ##rep2 
    
  }
  
}



ksS = lsS = jsS = array(dim=c(365,18,9))
for(w in 1:17){
  ksS[,1,1] =  sumishield_details[,5][1:365]
  ksS[,w+1,1] = c(rep(k0,w*7),sumishield_details[,5][1:(365-7*w)])
  
  ksS[,1,2] =  sumishield_details[,6][1:365]
  ksS[,w+1,2] = c(rep(k0,w*7),sumishield_details[,6][1:(365-7*w)])
  
  ksS[,1,3] =  sumishield_details[,7][1:365]
  ksS[,w+1,3] = c(rep(k0,w*7),sumishield_details[,7][1:(365-7*w)])
  
  
  lsS[,1,1] = sumishield_details[,2]
  lsS[,w+1,1] = c(rep(0,w*7),sumishield_details[,2][1:(365-7*w)])
  
  lsS[,1,2] = sumishield_details[,3]
  lsS[,w+1,2] = c(rep(0,w*7),sumishield_details[,3][1:(365-7*w)])
  
  lsS[,1,3] = sumishield_details[,4]
  lsS[,w+1,3] = c(rep(0,w*7),sumishield_details[,4][1:(365-7*w)])
  
  
  jsS[,w,1] = 1 - ksS[,w,1] - lsS[,w,1]
  jsS[,w,2] = 1 - ksS[,w,2] - lsS[,w,2]
  jsS[,w,3] = 1 - ksS[,w,3] - lsS[,w,3]
  
}
ksS[,,4:6] = ksS[,,1:3]
ksS[,,7:9] = ksS[,,1:3]

lsS[,,4:6] = lsS[,,1:3]
lsS[,,7:9] = lsS[,,1:3]

jsS[,,4:6] = jsS[,,1:3]
jsS[,,7:9] = jsS[,,1:3]

jsS[,18,] = 1 - ksS[,18,] - lsS[,18,]


s_IRS_Sumi1 = r_IRS_Sumi1 = array(dim=c(365,18,9))
for(j in 1:9){
  for(w in 1:18){
    s_IRS_Sumi1[,w,j] = ksS[,w,j]/k0 ##feed2 = when IRS is implemented in month 1 (Nov)
    r_IRS_Sumi1[,w,j] = (1 - ksS[,w,j]/k0)*(jsS[,w,j]/(lsS[,w,j]+jsS[,w,j])) ##rep2 
  }
}


w_Acte1[,,1,] = w_Sumi1[,,1,] = rep(1,365) 
## Probability that a mosquito bites and survives in the presence of indoor vector control
for(k in 1:3){
  for(j in 1:18){
    for(i in 1:365){
      PHI_B = PHI_B_mut[k]
      PHI_I = PHI_I_mut[k]
      w_Acte1[i,j,2,1+cl[k]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 
      ## probability of surviving biting given that there is ITN
      w_Acte1[i,j,3,1+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acte1[i,j,1])*s_IRS_Acte1[i,j,1]	
      ##			probability of surviving biting given that there is IRS
      w_Acte1[i,j,4,1+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acte1[i,j,1])*s_ITN[i+547]*s_IRS_Acte1[i,j,1] + 
        (PHI_I - PHI_B)*(1-r_IRS_Acte1[i,j,1])*s_IRS_Acte1[i,j,1] 
      ## probability of surviving biting given that there is ITN & IRS
      
      w_Sumi1[i,j,2,1+cl[k]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 
      ## probability of surviving biting given that there is ITN
      w_Sumi1[i,j,3,1+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumi1[i,j,1])*s_IRS_Sumi1[i,j,1]	
      ##			probability of surviving biting given that there is IRS
      w_Sumi1[i,j,4,1+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumi1[i,j,1])*s_ITN[i+547]*s_IRS_Sumi1[i,j,1] + 
        (PHI_I - PHI_B)*(1-r_IRS_Sumi1[i,j,1])*s_IRS_Sumi1[i,j,1] 
      ## probability of surviving biting given that there is ITN & IRS
      
      
      w_Acte1[i,j,2,2+cl[k]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 
      ## probability of surviving biting given that there is ITN
      w_Acte1[i,j,3,2+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acte1[i,j,2])*s_IRS_Acte1[i,j,2]	
      ##			probability of surviving biting given that there is IRS
      w_Acte1[i,j,4,2+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acte1[i,j,2])*s_ITN[i+547]*s_IRS_Acte1[i,j,2] + 
        (PHI_I - PHI_B)*(1-r_IRS_Acte1[i,j,2])*s_IRS_Acte1[i,j,2] 
      ## probability of surviving biting given that there is ITN & IRS
      
      w_Sumi1[i,j,2,2+cl[k]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 
      ## probability of surviving biting given that there is ITN
      w_Sumi1[i,j,3,2+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumi1[i,j,2])*s_IRS_Sumi1[i,j,2]	
      ##			probability of surviving biting given that there is IRS
      w_Sumi1[i,j,4,2+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumi1[i,j,2])*s_ITN[i+547]*s_IRS_Sumi1[i,j,2] + 
        (PHI_I - PHI_B)*(1-r_IRS_Sumi1[i,j,2])*s_IRS_Sumi1[i,j,2] 
      ## probability of surviving biting given that there is ITN & IRS
      
      
      w_Acte1[i,j,2,3+cl[k]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 
      ## probability of surviving biting given that there is ITN
      w_Acte1[i,j,3,3+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acte1[i,j,3])*s_IRS_Acte1[i,j,3]	
      ##			probability of surviving biting given that there is IRS
      w_Acte1[i,j,4,3+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acte1[i,j,3])*s_ITN[i+547]*s_IRS_Acte1[i,j,3] + 
        (PHI_I - PHI_B)*(1-r_IRS_Acte1[i,j,3])*s_IRS_Acte1[i,j,3] 
      ## probability of surviving biting given that there is ITN & IRS
      
      w_Sumi1[i,j,2,3+cl[k]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 
      ## probability of surviving biting given that there is ITN
      w_Sumi1[i,j,3,3+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumi1[i,j,3])*s_IRS_Sumi1[i,j,3]	
      ##			probability of surviving biting given that there is IRS
      w_Sumi1[i,j,4,3+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumi1[i,j,3])*s_ITN[i+547]*s_IRS_Sumi1[i,j,3] + 
        (PHI_I - PHI_B)*(1-r_IRS_Sumi1[i,j,3])*s_IRS_Sumi1[i,j,3] 
      ## probability of surviving biting given that there is ITN & IRS
      
    }
    
  }
  
}


prop_this_weekM = c(prop_houses_sprayed_WeeklyM[1],diff(prop_houses_sprayed_WeeklyM)[1:17])
prop_this_weekB = c(prop_houses_sprayed_WeeklyB[1],diff(prop_houses_sprayed_WeeklyB)[1:17])

w_Acte = yy_Acte = z_Acte = w_Sumi = yy_Sumi = z_Sumi = array(dim=c(365,4,9) )

for(j in 1:9){
  for(i in 1:365){
    w_Acte[i,1,j] = sum(w_Acte1[i,,1,j] * prop_this_weekM)
    w_Acte[i,2,j] = sum(w_Acte1[i,,2,j] * prop_this_weekM)
    w_Acte[i,3,j] = sum(w_Acte1[i,,3,j] * prop_this_weekM)
    w_Acte[i,4,j] = sum(w_Acte1[i,,4,j] * prop_this_weekM)
    
    w_Sumi[i,1,j] = sum(w_Sumi1[i,,1,j] * prop_this_weekB)
    w_Sumi[i,2,j] = sum(w_Sumi1[i,,2,j] * prop_this_weekB)
    w_Sumi[i,3,j] = sum(w_Sumi1[i,,3,j] * prop_this_weekB)
    w_Sumi[i,4,j] = sum(w_Sumi1[i,,4,j] * prop_this_weekB)
  }
  
}


## Probability of any bite (if there is IRS, a mosquito may bite and then die immediately afterwards)
yy_Acte[,1,] = w_Acte[,1,] 
yy_Acte[,2,] = w_Acte[,2,]

for(k in 1:3){
  for(j in 1:18){
    for(i in 1:365){
      PHI_B = PHI_B_mut[k]
      PHI_I = PHI_I_mut[k]
      
      yy_Acte1[i,j,3,1+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acte1[i,j,1])
      yy_Acte1[i,j,4,1+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acte1[i,j,1])*s_ITN[i+547] + 
        (PHI_I - PHI_B)*(1-r_IRS_Acte1[i,j,1])
      
      yy_Sumi1[i,j,3,1+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumi1[i,j,1])
      yy_Sumi1[i,j,4,1+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumi1[i,j,1])*s_ITN[i+547] + 
        (PHI_I - PHI_B)*(1-r_IRS_Sumi1[i,j,1])
      
      yy_Acte1[i,j,3,2+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acte1[i,j,2])
      yy_Acte1[i,j,4,2+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acte1[i,j,2])*s_ITN[i+547] + 
        (PHI_I - PHI_B)*(1-r_IRS_Acte1[i,j,2])
      
      yy_Sumi1[i,j,3,2+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumi1[i,j,2])
      yy_Sumi1[i,j,4,2+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumi1[i,j,2])*s_ITN[i+547] + 
        (PHI_I - PHI_B)*(1-r_IRS_Sumi1[i,j,2])
      
      yy_Acte1[i,j,3,3+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acte1[i,j,3])
      yy_Acte1[i,j,4,3+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acte1[i,j,3])*s_ITN[i+547] + 
        (PHI_I - PHI_B)*(1-r_IRS_Acte1[i,j,3])
      
      yy_Sumi1[i,j,3,3+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumi1[i,j,3])
      yy_Sumi1[i,j,4,3+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumi1[i,j,3])*s_ITN[i+547] + 
        (PHI_I - PHI_B)*(1-r_IRS_Sumi1[i,j,3])
      
    }
  }  
}

for(j in 1:9){
  for(i in 1:365){
    yy_Acte[i,3,j] = sum(yy_Acte1[i,,3,j] * prop_this_weekM)
    yy_Acte[i,4,j] = sum(yy_Acte1[i,,4,j] * prop_this_weekM)
    
    yy_Sumi[i,3,j] = sum(yy_Sumi1[i,,3,j] * prop_this_weekB)
    yy_Sumi[i,4,j] = sum(yy_Sumi1[i,,4,j] * prop_this_weekB)
  }  
}


## supplementary figure 1a
plot(yy_Acte[,1,1] ~ time[1:365],ylim=c(0,1),pch="",
     ylab = "Probability mosquito bites",
     col="black",
     main = "",cex.main=1.2,xlim=c(1,240),xaxt="n",
     xlab="Time in months",yaxt="n",cex.lab=1.4,cex.axis=1.4,cex=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,1,0.2),cex.lab=1.4,cex.axis=1.4)
axis(1,at=seq(15,240,30),labels = c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),
     cex.axis = 1.4)

colsd=c("darkred","red","orange","blue")
for(i in 3){
  lines(yy_Acte[1:240,i,1] ~ time[1:240],col="darkblue",lwd=2)
  lines(yy_Sumi[1:240,i,1] ~ time[1:240],col="aquamarine3",lwd=2)
}

polygon(c(time[1:240],rev(time[1:240])),
        c(yy_Acte[1:240,3,5],rev(yy_Acte[1:240,3,9])),
        col=adegenet::transp("darkblue",0.3),border=NA)
polygon(c(time[1:240],rev(time[1:240])),
        c(yy_Sumi[1:240,3,5],rev(yy_Sumi[1:240,3,9])),
        col=adegenet::transp("aquamarine3",0.3),border=NA)


legend("bottomleft",legend = c("Actellic 300CS","SumiShield"),
       col = c("darkblue","aquamarine3"),lwd = 2, lty=c(1,1),cex=1.2,bty="n")

## 
plot(w_Acte[,1,1] ~ time[1:365],ylim=c(0,1),pch="",
     ylab = "Probability mosquito bites and survives",
     col="black",
     main = "",cex.main=1.2,xlim=c(1,240),xaxt="n",
     xlab="Time in months",yaxt="n",cex.lab=1.4,cex.axis=1.4,cex=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,1,0.2),cex.lab=1.4,cex.axis=1.4)
axis(1,at=seq(15,240,30),labels = c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),
     cex.axis = 1.4)

colsd=c("darkred","red","orange","blue")
for(i in 3){
  lines(w_Acte[1:240,i,1] ~ time[1:240],col="darkblue",lwd = 2)
  lines(w_Sumi[1:240,i,1] ~ time[1:240],col="aquamarine3",lwd = 2)
}

polygon(c(time[1:240],rev(time[1:240])),
        c(w_Acte[1:240,3,5],rev(w_Acte[1:240,3,9])),
        col=adegenet::transp("darkblue",0.3),border=NA)
polygon(c(time[1:240],rev(time[1:240])),
        c(w_Sumi[1:240,3,5],rev(w_Sumi[1:240,3,9])),
        col=adegenet::transp("aquamarine3",0.3),border=NA)


## Probability repelled
z_Acte[,1,] = 0
z_Sumi[,1,] = 0

for(k in 1:3){
  for(j in 1:18){
    for(i in 1:365){
      PHI_B = PHI_B_mut[k]
      PHI_I = PHI_I_mut[k]
      z_Acte1[i,j,2,1+cl[k]] = PHI_B*r_ITN[i+547]
      z_Acte1[i,j,3,1+cl[k]] = PHI_I*r_IRS_Acte1[i,j,1]
      z_Acte1[i,j,4,1+cl[k]] = PHI_B*(r_IRS_Acte1[i,j,1] + (1-r_IRS_Acte1[i,j,1])*r_ITN[i+547]) + 
        (PHI_I - PHI_B)*r_IRS_Acte1[i,j,1]
      
      z_Sumi1[i,j,2,1+cl[k]] = PHI_B*r_ITN[i+547]
      z_Sumi1[i,j,3,1+cl[k]] = PHI_I*r_IRS_Sumi1[i,j,1]
      z_Sumi1[i,j,4,1+cl[k]] = PHI_B*(r_IRS_Sumi1[i,j,1] + (1-r_IRS_Sumi1[i,j,1])*r_ITN[i+547]) + 
        (PHI_I - PHI_B)*r_IRS_Sumi1[i,j,1]
      
      
      z_Acte1[i,j,2,2+cl[k]] = PHI_B*r_ITN[i+547]
      z_Acte1[i,j,3,2+cl[k]] = PHI_I*r_IRS_Acte1[i,j,2]
      z_Acte1[i,j,4,2+cl[k]] = PHI_B*(r_IRS_Acte1[i,j,2] + (1-r_IRS_Acte1[i,j,2])*r_ITN[i+547]) + 
        (PHI_I - PHI_B)*r_IRS_Acte1[i,j,2]
      
      z_Sumi1[i,j,2,2+cl[k]] = PHI_B*r_ITN[i+547]
      z_Sumi1[i,j,3,2+cl[k]] = PHI_I*r_IRS_Sumi1[i,j,2]
      z_Sumi1[i,j,4,2+cl[k]] = PHI_B*(r_IRS_Sumi1[i,j,2] + (1-r_IRS_Sumi1[i,j,2])*r_ITN[i+547]) + 
        (PHI_I - PHI_B)*r_IRS_Sumi1[i,j,2]
      
      z_Acte1[i,j,2,3+cl[k]] = PHI_B*r_ITN[i+547]
      z_Acte1[i,j,3,3+cl[k]] = PHI_I*r_IRS_Acte1[i,j,3]
      z_Acte1[i,j,4,3+cl[k]] = PHI_B*(r_IRS_Acte1[i,j,3] + (1-r_IRS_Acte1[i,j,3])*r_ITN[i+547]) + 
        (PHI_I - PHI_B)*r_IRS_Acte1[i,j,3]
      
      z_Sumi1[i,j,2,3+cl[k]] = PHI_B*r_ITN[i+547]
      z_Sumi1[i,j,3,3+cl[k]] = PHI_I*r_IRS_Sumi1[i,j,3]
      z_Sumi1[i,j,4,3+cl[k]] = PHI_B*(r_IRS_Sumi1[i,j,3] + (1-r_IRS_Sumi1[i,j,3])*r_ITN[i+547]) + 
        (PHI_I - PHI_B)*r_IRS_Sumi1[i,j,3]
      
    }  
  }
  
}

for(j in 1:9){
  for(i in 1:365){
    z_Acte[i,2,j] = sum(z_Acte1[i,,3,j] * prop_this_weekM)
    z_Acte[i,3,j] = sum(z_Acte1[i,,3,j] * prop_this_weekM)
    z_Acte[i,4,j] = sum(z_Acte1[i,,4,j] * prop_this_weekM)
    
    z_Sumi[i,2,j] = sum(z_Sumi1[i,,3,j] * prop_this_weekB)
    z_Sumi[i,3,j] = sum(z_Sumi1[i,,3,j] * prop_this_weekB)
    z_Sumi[i,4,j] = sum(z_Sumi1[i,,4,j] * prop_this_weekB)
  }
  
}


## supplementary figure 1c
plot(z_Acte[,1,1] ~ time[1:365],ylim=c(0,1),pch="",
     ylab = "Probability mosquito is repelled",
     col="black",
     main = "",cex.main=1.2,xlim=c(1,240),xaxt="n",
     xlab="Time in months",yaxt="n",cex.lab=1.4,cex.axis=1.4,cex=1.4)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,1,0.2),cex.lab=1.4,cex.axis=1.4)
axis(1,at=seq(15,240,30),labels = c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),
     cex.axis = 1.4)

colsd=c("darkred","red","orange","blue")
for(i in 3){
  lines(z_Acte[1:240,i,1] ~ time[1:240],col="darkblue",lwd=2)
  lines(z_Sumi[1:240,i,1] ~ time[1:240],col="aquamarine3",lwd=2)
}


polygon(c(time[1:240],rev(time[1:240])),
        c(z_Acte[1:240,3,5],rev(z_Acte[1:240,3,9])),
        col=adegenet::transp("darkblue",0.3),border=NA)
polygon(c(time[1:240],rev(time[1:240])),
        c(z_Sumi[1:240,3,5],rev(z_Sumi[1:240,3,9])),
        col=adegenet::transp("aquamarine3",0.3),border=NA)


par(xpd=NA,cex = 1.11)

text(x = -650, y = 2.8,"(A)")
text(x = -340, y = 2.8,"(B)")
text(x = -38, y = 2.8,"(C)")

text(x = -650, y = 1.2,"(D)")
text(x = -340, y = 1.2,"(E)")
text(x = -38, y = 1.2,"(F)")









actellic_details = dta1[,1:7]#readRDS("data/actellic_details_v2.Rdata")
sumishield_details = dta1[,c(1,8:13)]#readRDS("data/sumishield_details_v2.Rdata")

time = 1:365
par(mfrow=c(2,2))

## derived ITN/IRS quantities
## prob bites and survives

w_Acte = yy_Acte = z_Acte = 
  w_Sumi = yy_Sumi = z_Sumi = array(dim=c(240,4,9)) 

#############################
## 

##Species are different in each location 
PHI_B_mut = c(0.85,0.6,0.95) ## probability of bites in bed
PHI_I_mut = c(0.90,0.68,0.99) ## probability of bites indoors

PHI_B_boa = c(0.85,0.6,0.95)  ## probability of bites in bed
PHI_I_boa = c(0.90,0.68,0.99) ## probability of bites indoors


k0 = 0.699
ksA = actellic_details[,5]
lsA = actellic_details[,2]
jsA = 1 - actellic_details[,2] - actellic_details[,5]

s_IRS_Acte = ksA/k0 ##feed2
r_IRS_Acte = (1 - ksA/k0)*(jsA/(lsA+jsA)) ##rep2

ksS = sumishield_details[,5]
lsS = sumishield_details[,2]
jsS = 1 - sumishield_details[,2] - sumishield_details[,5]

s_IRS_Sumi = ksS/k0 ##feed2
r_IRS_Sumi = (1 - ksS/k0)*(jsS/(lsS+jsS)) ##rep2

## lower
ksAl = actellic_details[,6]
lsAl = actellic_details[,3]
jsAl = 1 - actellic_details[,6] - actellic_details[,3]

s_IRS_Actel = ksAl/k0 ##feed2
r_IRS_Actel = (1 - ksAl/k0)*(jsAl/(lsAl+jsAl)) ##rep2

ksSl = sumishield_details[,6]
lsSl = sumishield_details[,3]
jsSl = 1 - sumishield_details[,6] - sumishield_details[,3]

s_IRS_Sumil = ksSl/k0 ##feed2
r_IRS_Sumil = (1 - ksSl/k0)*(jsSl/(lsSl+jsSl)) ##rep2

## upper
ksAu = actellic_details[,7]
lsAu = actellic_details[,4]
jsAu = 1 - actellic_details[,7] - actellic_details[,4]

s_IRS_Acteu = ksAu/k0 ##feed2
r_IRS_Acteu = (1 - ksAu/k0)*(jsAu/(lsAu+jsAu)) ##rep2

ksSu = sumishield_details[,7]
lsSu = sumishield_details[,4]
jsSu = 1 - sumishield_details[,7] - sumishield_details[,4]

s_IRS_Sumiu = ksSu/k0 ##feed2
r_IRS_Sumiu = (1 - ksSu/k0)*(jsSu/(lsSu+jsSu)) ##rep2

w_Acte[,1,] = w_Sumi[,1,] = rep(1,240) ## Probability that a mosquito bites and survives in the presence of indoor vector control

cl = c(0,3,6)
for(j in 1:3){
  for(i in 1:240){
    PHI_B = PHI_B_mut[j]
    PHI_I = PHI_I_mut[j]
    w_Acte[i,2,1+cl[j]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 ## probability of surviving biting given that there is ITN
    w_Acte[i,3,1+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acte[i])*s_IRS_Acte[i]	##			probability of surviving biting given that there is IRS
    w_Acte[i,4,1+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acte[i])*s_ITN[i+547]*s_IRS_Acte[i] + (PHI_I - PHI_B)*(1-r_IRS_Acte[i])*s_IRS_Acte[i] ## probability of surviving biting given that there is ITN & IRS
    
    w_Acte[i,2,2+cl[j]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 ## probability of surviving biting given that there is ITN
    w_Acte[i,3,2+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Actel[i])*s_IRS_Actel[i]	##			probability of surviving biting given that there is IRS
    w_Acte[i,4,2+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Actel[i])*s_ITN[i+547]*s_IRS_Actel[i] + (PHI_I - PHI_B)*(1-r_IRS_Actel[i])*s_IRS_Actel[i] ## probability of surviving biting given that there is ITN & IRS
    
    w_Acte[i,2,3+cl[j]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 ## probability of surviving biting given that there is ITN
    w_Acte[i,3,3+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acteu[i])*s_IRS_Acteu[i]	##			probability of surviving biting given that there is IRS
    w_Acte[i,4,3+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acteu[i])*s_ITN[i+547]*s_IRS_Acteu[i] + (PHI_I - PHI_B)*(1-r_IRS_Acteu[i])*s_IRS_Acteu[i] ## probability of surviving biting given that there is ITN & IRS
  }  
}

cl = c(0,3,6)
for(j in 1:3){
  for(i in 1:240){
    PHI_B = PHI_B_mut[j]
    PHI_I = PHI_I_mut[j]
    
    w_Sumi[i,2,1+cl[j]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 ## probability of surviving biting given that there is ITN
    w_Sumi[i,3,1+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumi[i])*s_IRS_Sumi[i]	##			probability of surviving biting given that there is IRS
    w_Sumi[i,4,1+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumi[i])*s_ITN[i+547]*s_IRS_Sumi[i] + (PHI_I - PHI_B)*(1-r_IRS_Sumi[i])*s_IRS_Sumi[i] ## probability of surviving biting given that there is ITN & IRS
    
    w_Sumi[i,2,2+cl[j]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 ## probability of surviving biting given that there is ITN
    w_Sumi[i,3,2+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumil[i])*s_IRS_Sumil[i]	##			probability of surviving biting given that there is IRS
    w_Sumi[i,4,2+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumil[i])*s_ITN[i+547]*s_IRS_Sumil[i] + (PHI_I - PHI_B)*(1-r_IRS_Sumil[i])*s_IRS_Sumil[i] ## probability of surviving biting given that there is ITN & IRS
    
    w_Sumi[i,2,3+cl[j]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 ## probability of surviving biting given that there is ITN
    w_Sumi[i,3,3+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumiu[i])*s_IRS_Sumiu[i]	##			probability of surviving biting given that there is IRS
    w_Sumi[i,4,3+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumiu[i])*s_ITN[i+547]*s_IRS_Sumiu[i] + (PHI_I - PHI_B)*(1-r_IRS_Sumiu[i])*s_IRS_Sumiu[i] ## probability of surviving biting given that there is ITN & IRS
  }
}
yy_Acte[,1,] = w_Acte[,1,] 
yy_Acte[,2,] = w_Acte[,2,]

yy_Sumi[,1,] = w_Sumi[,1,] 
yy_Sumi[,2,] = w_Sumi[,2,]

for(j in 1:3){
  for(i in 1:240){
    PHI_B = PHI_B_mut[j]
    PHI_I = PHI_I_mut[j]
    yy_Acte[i,3,1+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acte[i])
    yy_Acte[i,4,1+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acte[i])*s_ITN[i+547] + (PHI_I - PHI_B)*(1-r_IRS_Acte[i])
    
    yy_Acte[i,3,2+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Actel[i])
    yy_Acte[i,4,2+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Actel[i])*s_ITN[i+547] + (PHI_I - PHI_B)*(1-r_IRS_Actel[i])
    
    yy_Acte[i,3,3+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acteu[i])
    yy_Acte[i,4,3+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acteu[i])*s_ITN[i+547] + (PHI_I - PHI_B)*(1-r_IRS_Acteu[i])
  }
}

for(j in 1:3){
  for(i in 1:240){
    PHI_B = PHI_B_boa[j]
    PHI_I = PHI_I_boa[j]
    yy_Sumi[i,3,1+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumi[i])
    yy_Sumi[i,4,1+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumi[i])*s_ITN[i+547] + (PHI_I - PHI_B)*(1-r_IRS_Sumi[i])
    
    yy_Sumi[i,3,2+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumil[i])
    yy_Sumi[i,4,2+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumil[i])*s_ITN[i+547] + (PHI_I - PHI_B)*(1-r_IRS_Sumil[i])
    
    yy_Sumi[i,3,3+cl[j]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumiu[i])
    yy_Sumi[i,4,3+cl[j]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumiu[i])*s_ITN[i+547] + (PHI_I - PHI_B)*(1-r_IRS_Sumiu[i])
  }  
}

## Probability repelled
z_Acte[,1,1] = 0
z_Sumi[,1,1] = 0

for(j in 1:3){
  for(i in 1:240){
    PHI_B = PHI_B_boa[j]
    PHI_I = PHI_I_boa[j]
    
    z_Acte[i,2,1+cl[j]] = PHI_B*r_ITN[i+547]
    z_Acte[i,3,1+cl[j]] = PHI_I*r_IRS_Acte[i]
    z_Acte[i,4,1+cl[j]] = PHI_B*(r_IRS_Acte[i] + (1-r_IRS_Acte[i])*r_ITN[i+547]) + (PHI_I - PHI_B)*r_IRS_Acte[i]
    
    z_Acte[i,2,2+cl[j]] = PHI_B*r_ITN[i+547]
    z_Acte[i,3,2+cl[j]] = PHI_I*r_IRS_Actel[i]
    z_Acte[i,4,2+cl[j]] = PHI_B*(r_IRS_Actel[i] + (1-r_IRS_Actel[i])*r_ITN[i+547]) + (PHI_I - PHI_B)*r_IRS_Actel[i]
    
    z_Acte[i,2,3+cl[j]] = PHI_B*r_ITN[i+547]
    z_Acte[i,3,3+cl[j]] = PHI_I*r_IRS_Acteu[i]
    z_Acte[i,4,3+cl[j]] = PHI_B*(r_IRS_Acteu[i] + (1-r_IRS_Acteu[i])*r_ITN[i+547]) + (PHI_I - PHI_B)*r_IRS_Acteu[i]
    
    z_Sumi[i,2,1+cl[j]] = PHI_B*r_ITN[i+547]
    z_Sumi[i,3,1+cl[j]] = PHI_I*r_IRS_Sumi[i]
    z_Sumi[i,4,1+cl[j]] = PHI_B*(r_IRS_Sumi[i] + (1-r_IRS_Sumi[i])*r_ITN[i+547]) + (PHI_I - PHI_B)*r_IRS_Sumi[i]
    
    z_Sumi[i,2,2+cl[j]] = PHI_B*r_ITN[i+547]
    z_Sumi[i,3,2+cl[j]] = PHI_I*r_IRS_Sumil[i]
    z_Sumi[i,4,2+cl[j]] = PHI_B*(r_IRS_Sumil[i] + (1-r_IRS_Sumil[i])*r_ITN[i+547]) + (PHI_I - PHI_B)*r_IRS_Sumil[i]
    
    z_Sumi[i,2,3+cl[j]] = PHI_B*r_ITN[i+547]
    z_Sumi[i,3,3+cl[j]] = PHI_I*r_IRS_Sumiu[i]
    z_Sumi[i,4,3+cl[j]] = PHI_B*(r_IRS_Sumiu[i] + (1-r_IRS_Sumiu[i])*r_ITN[i+547]) + (PHI_I - PHI_B)*r_IRS_Sumiu[i]
  }
}
## waning usage of IRS with time


prop_mod_Acte_temp = c(0.891472868,0.81212081,0.639282555,0.481047262,0.386204738,0.304873387,0.166894353)

##And sorted for weeks:
prop_mod_Acte = c(rep(1,10),rep(prop_mod_Acte_temp,each=4))## for time series of spray campaign

true_cover_irs_Acte =  rep(prop_mod_Acte*0.96,each=7)

#House coverage: Matutuine district 96 %
irs_cov_no_loss_Acte = rep(0.96,30*8)
irs_cov_Acte = true_cover_irs_Acte


prop_mod_Sumi_temp = c(0.96460177,0.966330299,0.945153088,0.923357485,0.917585479,0.904508888,0.895619142)

##And sorted for weeks:
prop_mod_Sumi = c(rep(1,10),rep(prop_mod_Sumi_temp,each=4))## for time series of spray campaign

true_cover_irs_Sumi = rep(prop_mod_Sumi*prop_houses_sprayed_WeeklyB,each=7)

#House coverage: Boane (sumi) district 97 %, Manhica district (Palmeira) 98 % 
irs_cov_no_loss_Sumi = rep(0.97,30*8)
irs_cov_Sumi = rep(prop_mod_Sumi*0.97,each=7)

time = 1:365
plot(irs_cov_no_loss_Acte[61:240] ~ time[61:240],ylab = "Community IRS cover (%)",
     ylim=c(0,1),col="black",pch="",
     main = "",cex.main=1.2,xlim=c(1,240),xaxt="n",
     xlab="Time in months",yaxt="n",cex.lab=1.2,cex.axis=1.2,cex=1.2)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.2,cex.axis=1.2)
axis(1,at=seq(0,230,30)+15,labels = c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),cex.axis = 1.4)


lines(irs_cov_no_loss_Acte[61:240] ~ time[61:240],lty=1,lwd=2,col = "darkblue") ## IRS only no loss 
lines(irs_cov_Acte[61:240] ~ time[61:240],lty=3,lwd=2,col = "darkblue") ## IRS with loss

lines(irs_cov_no_loss_Sumi[61:240] ~ c(time[61:240]+1),lty=1,lwd=2,col = "aquamarine3")
lines(irs_cov_Sumi[61:240] ~ time[61:240],lty=3,lwd=2,col = "aquamarine3")

legend("bottomleft",legend = c("Matutuine","Boane","IRS cover, no loss", "IRS cover, observed loss"),
       col = c("darkblue","aquamarine3","black","black"),lwd = 2, lty=c(NA,NA,1,3),
       pch=c(15,15,NA,NA),cex=1.2,bty="n")

cov1A = cov1S = cov2A = cov2S = array(dim=c(240,4))

# itn_cov_Acte = 0.52
# itn_cov_Sumi = 0.73

## Table 1 data
matu_net_cov = c(27.9,  ## nov
                 mean(c(32.6,42.1)), ##dec
                 mean(c(43.4,43.2,33.3)),##jan
                 mean(c(55.0,44.3,48.2)),##feb
                 mean(c(61.2,56.8,40.7)),##mar
                 mean(c(51.9,61.4,66.7)),##apr
                 mean(c(48.9,55.6)),##may
                 29.6)##june

itn_cov_Acte_temp = c(rep(mean(c(matu_net_cov)),10),
                      rep(mean(c(32.6,42.1)),4), ##dec
                      rep(mean(c(43.4,43.2,33.3)),4),##jan
                      rep(mean(c(55.0,44.3,48.2)),4),##feb
                      rep(mean(c(61.2,56.8,40.7)),4),##mar
                      rep(mean(c(51.9,61.4,66.7)),4),##apr
                      rep(mean(c(48.9,55.6)),4),##may
                      rep(29.6,4))##june



boan_net_cov = c(57.5,  ## nov
                 mean(c(64.6,74.5)), ##dec
                 mean(c(68.1,71.9,67.1)),##jan
                 mean(c(62.8,71.2,79)),##feb
                 mean(c(70.8,82.4,86.8)),##mar
                 mean(c(65.5,81.8,81.6)),##apr
                 mean(c(78.4,84.2)),##may
                 10.5)##june

itn_cov_Boan_temp = c(rep(mean(c(boan_net_cov)),10),
                      rep(mean(c(64.6,74.5)),4), ##dec
                      rep(mean(c(68.1,71.9,67.1)),4),##jan
                      rep(mean(c(62.8,71.2,79)),4),##feb
                      rep(mean(c(70.8,82.4,86.8)),4),##mar
                      rep(mean(c(65.5,81.8,81.6)),4),##apr
                      rep(mean(c(78.4,84.2)),4),##may
                      rep(10.5,4))

itn_cov_Acte = rep(itn_cov_Acte_temp,each=7)/100
itn_cov_Sumi = rep(itn_cov_Boan_temp,each=7)/100

## Here we are creating a matrix
## with the coverage or use of nets waning with time
## and the coverage of IRS either staying fixed, or also waning 
## when walls are washed
cov1A[,1] = 1
cov1A[,2] = itn_cov_Acte[1:240] ## ITN only
cov1A[,3] = irs_cov_no_loss_Acte ## IRS only
cov1A[,4] = itn_cov_Acte[1:240]*irs_cov_no_loss_Acte ## both interventions

cov2A[,1] = 1
cov2A[,2] = itn_cov_Acte[1:240] ## ITN only
cov2A[,3] = irs_cov_Acte[1:240] ## IRS only
cov2A[,4] = itn_cov_Acte[1:240]*irs_cov_Acte[1:240]## both interventions

cov1S[,1] = 1
cov1S[,2] = itn_cov_Sumi[1:240] ## ITN only
cov1S[,3] = irs_cov_no_loss_Sumi ## IRS only
cov1S[,4] = itn_cov_Sumi[1:240]*irs_cov_no_loss_Sumi ## both interventions

cov2S[,1] = 1
cov2S[,2] = itn_cov_Sumi[1:240] ## ITN only
cov2S[,3] = irs_cov_Sumi[1:240] ## IRS only
cov2S[,4] = itn_cov_Sumi[1:240]*irs_cov_Sumi[1:240] ## both interventions



## Entomological model parameters to estimate 

Q0 = 0.98  ## this is anthropophagy - we can use human blood index
chi = 0.86 ## this is endophily (pi_i)

fv0 = 0.333 ## biting rate 1 bite every 3 days
tau1 = 0.69 ## duration of host seeking, assumed to be constant between species (delta_10, altered to delta_1 when interventions used)
tau2 = 1/fv0-tau1  ## indoor feeding endophagy (delta_2)
av0 = Q0*fv0
mu0 = 0.132 ## background mortality from external sources
p10 = exp(-mu0*tau1) ## 
p2 = exp(-mu0*tau2)  ## probability of surviving resting period in absence of intervntion


## These are the adjusted w, z, when coverage is changing
## so these are the intervention coverages
zhi1A = whi1A = zhi2A = whi2A = array(dim=c(240,4,9))
zhi1S = whi1S = zhi2S = whi2S = array(dim=c(240,4,9))
for(i in 1:9){
  zhi1A[,,i] = cov1A*z_Acte[1:240,,i]
  whi1A[,,i] = cov1A*w_Acte[1:240,,i]
  
  zhi1S[,,i] =cov1S*z_Sumi[1:240,,i]
  whi1S[,,i] =cov1S*w_Sumi[1:240,,i]
  
  
  zhi2A[,,i] =cov2A*z_Acte[1:240,,i]
  whi2A[,,i] =cov2A*w_Acte[1:240,,i]
  
  zhi2S[,,i] =cov2S*z_Sumi[1:240,,i]
  whi2S[,,i] =cov2S*w_Sumi[1:240,,i]
  
}


zbar1A = wbar1A = zbar2A = wbar2A = array(dim=c(240,4,9)) 
zbar1S = wbar1S = zbar2S = wbar2S = array(dim=c(240,4,9)) 
for(j in 1:9){
  for(i in 1:4){
    zbar1A[,i,j] = Q0*zhi1A[,i,j]
    wbar1A[,i,j] = (1 - Q0) + Q0*whi1A[,i,j]
    zbar2A[,i,j] = Q0*zhi2A[,i,j]
    wbar2A[,i,j] = (1 - Q0) + Q0*whi2A[,i,j]
    
    zbar1S[,i,j] = Q0*zhi1S[,i,j]
    wbar1S[,i,j] = (1 - Q0) + Q0*whi1S[,i,j]
    zbar2S[,i,j] = Q0*zhi2S[,i,j]
    wbar2S[,i,j] = (1 - Q0) + Q0*whi2S[,i,j]
  }
  
}


## From Walker et al 2016
## Mosquito feeding rate (tau1 is delta10, tau2 is delta2 in the methods)
fR1A = 1 / ((tau1/(1 - zbar1A)) + tau2)
mu1A = -fR1A*log((wbar1A*p10/(1 - zbar1A*p10))*p2) 
Q1A = 1 - (1-Q0)/wbar1A

fR2A = 1 / ((tau1/(1 - zbar2A)) + tau2)
mu2A = -fR2A*log((wbar2A*p10/(1 - zbar2A*p10))*p2) 
Q2A = 1 - (1-Q0)/wbar2A


fR1S = 1 / ((tau1/(1 - zbar1S)) + tau2)
mu1S = -fR1S*log((wbar1S*p10/(1 - zbar1S*p10))*p2) 
Q1S = 1 - (1-Q0)/wbar1S

fR2S = 1 / ((tau1/(1 - zbar2S)) + tau2)
mu2S = -fR2S*log((wbar2S*p10/(1 - zbar2S*p10))*p2) 
Q2S = 1 - (1-Q0)/wbar2S


## Rate at which a person in the popn is bitten by mosquitoes is
lambda1A = lambda2A = array(dim = c(240,4,9))
lambda1S = lambda2S = array(dim = c(240,4,9))

for(j in 1:9){
  for(i in 1:4){
    lambda1A[,,j] = (Q1A[,,j]*fR1A[,,j]*yy_Acte[,i,j])/whi1A[,i,j]
    lambda2A[,,j] = (Q2A[,,j]*fR2A[,,j]*yy_Acte[,i,j])/whi2A[,i,j]
    
    lambda1S[,,j] = (Q1S[,,j]*fR1S[,,j]*yy_Sumi[,i,j])/whi1S[,i,j]
    lambda2S[,,j] = (Q2S[,,j]*fR2S[,,j]*yy_Sumi[,i,j])/whi2S[,i,j]
    
  }
}


 


## Additional infectious bites per person per year 
Estimated_added_EIR = array(dim=c(240,2,9))
Estimated_added_EIR[,1,] = (lambda2A[,4,] - lambda1A[,4,])
Estimated_added_EIR[,2,] = (lambda2S[,4,] - lambda1S[,4,])

c(sum(Estimated_added_EIR[1:30,1,])/30,sum(Estimated_added_EIR[31:60,1,])/30,
  sum(Estimated_added_EIR[61:90,1,])/30,sum(Estimated_added_EIR[91:120,1,])/30,
  sum(Estimated_added_EIR[121:150,1,])/30,sum(Estimated_added_EIR[151:180,1,])/30)

c(sum(Estimated_added_EIR[1:30,2,])/30,sum(Estimated_added_EIR[31:60,2,])/30,
  sum(Estimated_added_EIR[61:90,2,])/30,sum(Estimated_added_EIR[91:120,2,])/30,
  sum(Estimated_added_EIR[121:150,2,])/30,sum(Estimated_added_EIR[151:180,2,])/30)


## Additional infectious bites per person per year 
Estimated_propn_increase_EIR = array(dim=c(240,2,9))
Estimated_propn_increase_EIR[,1,] = (lambda2A[,4,] - lambda1A[,4,])/lambda2A[,4,]
Estimated_propn_increase_EIR[,2,] = (lambda2S[,4,] - lambda1S[,4,])/lambda2S[,4,]

c(sum(Estimated_propn_increase_EIR[1:30,1,])/30,sum(Estimated_propn_increase_EIR[31:60,1,])/30,
  sum(Estimated_propn_increase_EIR[61:90,1,])/30,sum(Estimated_propn_increase_EIR[91:120,1,])/30,
  sum(Estimated_propn_increase_EIR[121:150,1,])/30,sum(Estimated_propn_increase_EIR[151:180,1,])/30)

c(sum(Estimated_propn_increase_EIR[1:30,2,])/30,sum(Estimated_propn_increase_EIR[31:60,2,])/30,
  sum(Estimated_propn_increase_EIR[61:90,2,])/30,sum(Estimated_propn_increase_EIR[91:120,2,])/30,
  sum(Estimated_propn_increase_EIR[121:150,2,])/30,sum(Estimated_propn_increase_EIR[151:180,2,])/30)

mean(Estimated_propn_increase_EIR[15:45,1,])##sep
mean(Estimated_propn_increase_EIR[46:75,1,])##oct
mean(Estimated_propn_increase_EIR[76:105,1,])##nov
mean(Estimated_propn_increase_EIR[106:135,1,])##dec
mean(Estimated_propn_increase_EIR[136:165,1,])##jan
mean(Estimated_propn_increase_EIR[166:195,1,])##feb
mean(Estimated_propn_increase_EIR[196:225,1,])##mar
mean(Estimated_propn_increase_EIR[226:240,1,])##part april

mean(Estimated_propn_increase_EIR[15:45,2,])##sep
mean(Estimated_propn_increase_EIR[46:75,2,])##oct
mean(Estimated_propn_increase_EIR[76:105,2,])##nov
mean(Estimated_propn_increase_EIR[106:135,2,])##dec
mean(Estimated_propn_increase_EIR[136:165,2,])##jan
mean(Estimated_propn_increase_EIR[166:195,2,])##feb
mean(Estimated_propn_increase_EIR[196:225,2,])##mar
mean(Estimated_propn_increase_EIR[226:240,2,])##part april

## Figure 5b

plot(Estimated_propn_increase_EIR[61:240,1,1] ~ time[61:240],ylim=c(0,1),pch="",
     ylab = "",
     col="black",
     main = "",cex.main=1.2,xlim=c(1,240),xaxt="n",
     xlab="Time in months",yaxt="n",cex.lab=1.2,cex.axis=1.2,cex=1.2)
mtext(side=2, line =4.3,
      "Relative increase in daily bites ")
mtext(side=2, line =3, "due to modifications (%)")
axis(2,las=2,at=seq(0,1,0.2),label=seq(0,100,20),cex.lab=1.2,cex.axis=1.2)

axis(1,at=seq(0,230,30)+15,labels = c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),cex.axis = 1.4)

colsd = c("darkblue","aquamarine3")
for(i in 1:2){
  lines(Estimated_propn_increase_EIR[61:240,i,1] ~ time[61:240],col=colsd[i],lty=2,lwd=2)
}

min_EIR = max_EIR = array(dim=c(length(61:240),2))
for(j in 1:2){
  for(i in 1:length(61:240)){
    min_EIR[i,j] = min(Estimated_propn_increase_EIR[60+i,j,])
    max_EIR[i,j] = max(Estimated_propn_increase_EIR[60+i,j,])
    
  }  
}

polygon(c(time[61:240],rev(time[61:240])),
        c(min_EIR[,1],rev(max_EIR[,1])),
        col=adegenet::transp("darkblue",0.3),border=NA)
polygon(c(time[61:240],rev(time[61:240])),
        c(min_EIR[,2],rev(max_EIR[,2])),
        col=adegenet::transp("aquamarine",0.3),border=NA)



Estimated_added_EIR = array(dim=c(240,2,9))
Estimated_added_EIR[,1,] = (lambda2A[,3,] - lambda1A[,3,])
Estimated_added_EIR[,2,] = (lambda2S[,3,] - lambda1S[,3,])

## Additional infectious bites per person per year 
Estimated_propn_increase_EIR = array(dim=c(240,2,9))
Estimated_propn_increase_EIR[,1,] = (lambda2A[,3,] - lambda1A[,3,])/lambda2A[,3,]
Estimated_propn_increase_EIR[,2,] = (lambda2S[,3,] - lambda1S[,3,])/lambda2S[,3,]

for(i in 1:2){
  lines(Estimated_propn_increase_EIR[61:240,i,1] ~ time[61:240],col=colsd[i],lty=1,lwd=1)
}

min_EIR = max_EIR = array(dim=c(length(61:240),2))
for(j in 1:2){
  for(i in 1:length(61:240)){
    min_EIR[i,j] = min(Estimated_propn_increase_EIR[60+i,j,])
    max_EIR[i,j] = max(Estimated_propn_increase_EIR[60+i,j,])
    
  }  
}

polygon(c(time[61:240],rev(time[61:240])),
        c(min_EIR[,1],rev(max_EIR[,1])),
        col=adegenet::transp("darkblue",0.3),border=NA)
polygon(c(time[61:240],rev(time[61:240])),
        c(min_EIR[,2],rev(max_EIR[,2])),
        col=adegenet::transp("aquamarine",0.3),border=NA)


mean(Estimated_propn_increase_EIR[1:30,1,])
mean(Estimated_propn_increase_EIR[31:60,1,])
mean(Estimated_propn_increase_EIR[61:92,1,])
mean(Estimated_propn_increase_EIR[93:122,1,])
mean(Estimated_propn_increase_EIR[123:153,1,])
mean(Estimated_propn_increase_EIR[154:180,1,])


legend("topleft",legend = c("Matutuine (assuming no LLIN)","Boane (assuming no LLIN)",
                            "Matutuine (LLIN use)","Boane (LLIN use)"),
       col = c("darkblue","aquamarine3","darkblue","aquamarine3"),lwd = c(2,2,1,1), lty=c(1,1,2,2),cex=1.2,bty="n")



####################################
##
##   Work out the probability that a feeding attempt by mosquito from species i ends in blood feeding on a person 
##   whilst considering the prolonged nature of application campaigns
##
####################################################

## Campaigns tend to take up to a few months to complete
## We assume the ratio of houses monitored per start month
## reflects the proportion of houses covered by the spray campaign
## in that month

## Houses were tracked from November, December and January (1 in Feb)
## 129/(129+88+27+1) is 52.65% of houses in the community have max protection in Nov
## 88/(129+88+27+1) is 35.92% of houses max protection Dec (52% - any modifications have 1 mont old protection)
## 27/(129+88+27+1) is 11.02% max protection Jan (35.92% - modif 1 month old, 52% - modif 2 month old)


## derived LLIN/IRS quantities
## prob bites and survives

w_Acte = yy_Acte = z_Acte = w_Sumi = yy_Sumi = z_Sumi = array(dim=c(365,4,9)) 
w_Acte1 = yy_Acte1 = z_Acte1 = w_Sumi1 = yy_Sumi1 = z_Sumi1 = array(dim=c(365,18,4,9))
# 
# 
# w_Acte1 = yy_Acte1 = z_Acte1 = w_Sumi1 = yy_Sumi1 = z_Sumi1 = array(dim=c(180,4)) 
# w_Acte2 = yy_Acte2 = z_Acte2 = w_Sumi2 = yy_Sumi2 = z_Sumi2 = array(dim=c(180,4)) 
# w_Acte3 = yy_Acte3 = z_Acte3 = w_Sumi3 = yy_Sumi3 = z_Sumi3 = array(dim=c(180,4)) 
# 

#############################
## 

##Species are different in each location 
PHI_B_mut = c(0.85,0.6,0.95) ## probability of bites in bed
PHI_I_mut = c(0.90,0.68,0.99) ## probability of bites indoors

PHI_B_boa = c(0.85,0.6,0.95)  ## probability of bites in bed
PHI_I_boa = c(0.90,0.68,0.99) ## probability of bites indoors

## Data for the proportion of houses sprayed
## over the course of the trial
prop_houses_sprayed_WeeklyB = 0.97*c(0,	0.027219794,	## August
                                     0.077014558,	0.136261919,	0.196901742,	0.250817066, ## Sept
                                     0.301047746,	0.347687015,	0.395348206,	0.464541108, ## Oct
                                     0.521818166,	0.581372602,	0.643327061,	0.710948828, ## Nov
                                     0.777620537,	0.847130239,	0.911498024,	0.930365931, ## Dec
                                     0.947042313,	0.96389906,	0.983400068,	0.991357264,	 ## Jan
                                     0.99238783,	1) ## Feb
prop_houses_sprayed_WeeklyM = 0.96*c(0.065358837,	0.193716785,## August
                                     0.334444596,	0.432141444,	0.533708885,	0.614060416,##sep
                                     0.667531537,	0.711462198,	0.769981823,	0.880751007,##oct
                                     0.919045204,	0.944027976,	0.97443187,	0.990922736,  ##nov
                                     0.991873307,	0.99744169,	 1,             1,          ##dec
                                     1,        1) ## jan

ksA = lsA = jsA = array(dim=c(365,18,9))
for(w in 1:17){
  ksA[,1,1] =  actellic_details[,5][1:365]
  ksA[,w+1,1] = c(rep(k0,w*7),actellic_details[,5][1:(365-7*w)])
  
  ksA[,1,2] =  actellic_details[,6][1:365]
  ksA[,w+1,2] = c(rep(k0,w*7),actellic_details[,6][1:(365-7*w)])
  
  ksA[,1,3] =  actellic_details[,7][1:365]
  ksA[,w+1,3] = c(rep(k0,w*7),actellic_details[,7][1:(365-7*w)])
  
  
  lsA[,1,1] = actellic_details[,2]
  lsA[,w+1,1] = c(rep(0,w*7),actellic_details[,2][1:(365-7*w)])
  
  lsA[,1,2] = actellic_details[,3]
  lsA[,w+1,2] = c(rep(0,w*7),actellic_details[,3][1:(365-7*w)])
  
  lsA[,1,3] = actellic_details[,4]
  lsA[,w+1,3] = c(rep(0,w*7),actellic_details[,4][1:(365-7*w)])
  
  
  jsA[,w,1] = 1 - ksA[,w,1] - lsA[,w,1]
  jsA[,w,2] = 1 - ksA[,w,2] - lsA[,w,2]
  jsA[,w,3] = 1 - ksA[,w,3] - lsA[,w,3]
  
}
ksA[,,4:6] = ksA[,,1:3]
ksA[,,7:9] = ksA[,,1:3]

lsA[,,4:6] = lsA[,,1:3]
lsA[,,7:9] = lsA[,,1:3]

jsA[,,4:6] = jsA[,,1:3]
jsA[,,7:9] = jsA[,,1:3]

jsA[,18,] = 1 - ksA[,18,] - lsA[,18,]



s_IRS_Acte1 = r_IRS_Acte1 = array(dim=c(365,18,9))
for(j in 1:9){
  for(w in 1:18){
    s_IRS_Acte1[,w,j] = ksA[,w,j]/k0 ##feed2 = when IRS is implemented in month 1 (Nov)
    r_IRS_Acte1[,w,j] = (1 - ksA[,w,j]/k0)*(jsA[,w,j]/(lsA[,w,j]+jsA[,w,j])) ##rep2 
    
  }
  
}



ksS = lsS = jsS = array(dim=c(365,18,9))
for(w in 1:17){
  ksS[,1,1] =  sumishield_details[,5][1:365]
  ksS[,w+1,1] = c(rep(k0,w*7),sumishield_details[,5][1:(365-7*w)])
  
  ksS[,1,2] =  sumishield_details[,6][1:365]
  ksS[,w+1,2] = c(rep(k0,w*7),sumishield_details[,6][1:(365-7*w)])
  
  ksS[,1,3] =  sumishield_details[,7][1:365]
  ksS[,w+1,3] = c(rep(k0,w*7),sumishield_details[,7][1:(365-7*w)])
  
  
  lsS[,1,1] = sumishield_details[,2]
  lsS[,w+1,1] = c(rep(0,w*7),sumishield_details[,2][1:(365-7*w)])
  
  lsS[,1,2] = sumishield_details[,3]
  lsS[,w+1,2] = c(rep(0,w*7),sumishield_details[,3][1:(365-7*w)])
  
  lsS[,1,3] = sumishield_details[,4]
  lsS[,w+1,3] = c(rep(0,w*7),sumishield_details[,4][1:(365-7*w)])
  
  
  jsS[,w,1] = 1 - ksS[,w,1] - lsS[,w,1]
  jsS[,w,2] = 1 - ksS[,w,2] - lsS[,w,2]
  jsS[,w,3] = 1 - ksS[,w,3] - lsS[,w,3]
  
}
ksS[,,4:6] = ksS[,,1:3]
ksS[,,7:9] = ksS[,,1:3]

lsS[,,4:6] = lsS[,,1:3]
lsS[,,7:9] = lsS[,,1:3]

jsS[,,4:6] = jsS[,,1:3]
jsS[,,7:9] = jsS[,,1:3]

jsS[,18,] = 1 - ksS[,18,] - lsS[,18,]


s_IRS_Sumi1 = r_IRS_Sumi1 = array(dim=c(365,18,9))
for(j in 1:9){
  for(w in 1:18){
    s_IRS_Sumi1[,w,j] = ksS[,w,j]/k0 ##feed2 = when IRS is implemented in month 1 (Nov)
    r_IRS_Sumi1[,w,j] = (1 - ksS[,w,j]/k0)*(jsS[,w,j]/(lsS[,w,j]+jsS[,w,j])) ##rep2 
  }
}

w_Acte1[,,1,] = w_Sumi1[,,1,] = rep(1,365) 
## Probability that a mosquito bites and survives in the presence of indoor vector control
for(k in 1:3){
  for(j in 1:18){
    for(i in 1:365){
      PHI_B = PHI_B_mut[k]
      PHI_I = PHI_I_mut[k]
      w_Acte1[i,j,2,1+cl[k]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 ## probability of surviving biting given that there is ITN
      w_Acte1[i,j,3,1+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acte1[i,j,1])*s_IRS_Acte1[i,j,1]	##			probability of surviving biting given that there is IRS
      w_Acte1[i,j,4,1+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acte1[i,j,1])*s_ITN[i+547]*s_IRS_Acte1[i,j,1] + (PHI_I - PHI_B)*(1-r_IRS_Acte1[i,j,1])*s_IRS_Acte1[i,j,1] ## probability of surviving biting given that there is ITN & IRS
      
      w_Sumi1[i,j,2,1+cl[k]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 ## probability of surviving biting given that there is ITN
      w_Sumi1[i,j,3,1+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumi1[i,j,1])*s_IRS_Sumi1[i,j,1]	##			probability of surviving biting given that there is IRS
      w_Sumi1[i,j,4,1+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumi1[i,j,1])*s_ITN[i+547]*s_IRS_Sumi1[i,j,1] + (PHI_I - PHI_B)*(1-r_IRS_Sumi1[i,j,1])*s_IRS_Sumi1[i,j,1] ## probability of surviving biting given that there is ITN & IRS
      
      
      w_Acte1[i,j,2,2+cl[k]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 ## probability of surviving biting given that there is ITN
      w_Acte1[i,j,3,2+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acte1[i,j,2])*s_IRS_Acte1[i,j,2]	##			probability of surviving biting given that there is IRS
      w_Acte1[i,j,4,2+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acte1[i,j,2])*s_ITN[i+547]*s_IRS_Acte1[i,j,2] + (PHI_I - PHI_B)*(1-r_IRS_Acte1[i,j,2])*s_IRS_Acte1[i,j,2] ## probability of surviving biting given that there is ITN & IRS
      
      w_Sumi1[i,j,2,2+cl[k]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 ## probability of surviving biting given that there is ITN
      w_Sumi1[i,j,3,2+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumi1[i,j,2])*s_IRS_Sumi1[i,j,2]	##			probability of surviving biting given that there is IRS
      w_Sumi1[i,j,4,2+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumi1[i,j,2])*s_ITN[i+547]*s_IRS_Sumi1[i,j,2] + (PHI_I - PHI_B)*(1-r_IRS_Sumi1[i,j,2])*s_IRS_Sumi1[i,j,2] ## probability of surviving biting given that there is ITN & IRS
      
      
      w_Acte1[i,j,2,3+cl[k]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 ## probability of surviving biting given that there is ITN
      w_Acte1[i,j,3,3+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acte1[i,j,3])*s_IRS_Acte1[i,j,3]	##			probability of surviving biting given that there is IRS
      w_Acte1[i,j,4,3+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acte1[i,j,3])*s_ITN[i+547]*s_IRS_Acte1[i,j,3] + (PHI_I - PHI_B)*(1-r_IRS_Acte1[i,j,3])*s_IRS_Acte1[i,j,3] ## probability of surviving biting given that there is ITN & IRS
      
      w_Sumi1[i,j,2,3+cl[k]] = 1 - PHI_B + PHI_B*s_ITN[i+547]				 ## probability of surviving biting given that there is ITN
      w_Sumi1[i,j,3,3+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumi1[i,j,3])*s_IRS_Sumi1[i,j,3]	##			probability of surviving biting given that there is IRS
      w_Sumi1[i,j,4,3+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumi1[i,j,3])*s_ITN[i+547]*s_IRS_Sumi1[i,j,3] + (PHI_I - PHI_B)*(1-r_IRS_Sumi1[i,j,3])*s_IRS_Sumi1[i,j,3] ## probability of surviving biting given that there is ITN & IRS
      
    }
    
  }
  
}

prop_this_weekM = c(prop_houses_sprayed_WeeklyM[1],diff(prop_houses_sprayed_WeeklyM)[1:17])
prop_this_weekB = c(prop_houses_sprayed_WeeklyB[1],diff(prop_houses_sprayed_WeeklyB)[1:17])

w_Acte = yy_Acte = z_Acte = w_Sumi = yy_Sumi = z_Sumi = array(dim=c(365,4,9) )

for(j in 1:9){
  for(i in 1:365){
    w_Acte[i,1,j] = sum(w_Acte1[i,,1,j] * prop_this_weekM)
    w_Acte[i,2,j] = sum(w_Acte1[i,,2,j] * prop_this_weekM)
    w_Acte[i,3,j] = sum(w_Acte1[i,,3,j] * prop_this_weekM)
    w_Acte[i,4,j] = sum(w_Acte1[i,,4,j] * prop_this_weekM)
    
    w_Sumi[i,1,j] = sum(w_Sumi1[i,,1,j] * prop_this_weekB)
    w_Sumi[i,2,j] = sum(w_Sumi1[i,,2,j] * prop_this_weekB)
    w_Sumi[i,3,j] = sum(w_Sumi1[i,,3,j] * prop_this_weekB)
    w_Sumi[i,4,j] = sum(w_Sumi1[i,,4,j] * prop_this_weekB)
  }  
}






## Probability of any bite (if there is IRS, a mosquito may bite and then die immediately afterwards)
yy_Acte[,1,] = w_Acte[,1,] 
yy_Acte[,2,] = w_Acte[,2,]

# yy_Sumi[,1] = w_Sumi[,1] 
# yy_Sumi[,2] = w_Sumi[,2]

for(k in 1:3){
  for(j in 1:18){
    for(i in 1:365){
      PHI_B = PHI_B_mut[k]
      PHI_I = PHI_I_mut[k]
      
      yy_Acte1[i,j,3,1+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acte1[i,j,1])
      yy_Acte1[i,j,4,1+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acte1[i,j,1])*s_ITN[i+547] + (PHI_I - PHI_B)*(1-r_IRS_Acte1[i,j,1])
      
      yy_Sumi1[i,j,3,1+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumi1[i,j,1])
      yy_Sumi1[i,j,4,1+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumi1[i,j,1])*s_ITN[i+547] + (PHI_I - PHI_B)*(1-r_IRS_Sumi1[i,j,1])
      
      yy_Acte1[i,j,3,2+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acte1[i,j,2])
      yy_Acte1[i,j,4,2+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acte1[i,j,2])*s_ITN[i+547] + (PHI_I - PHI_B)*(1-r_IRS_Acte1[i,j,2])
      
      yy_Sumi1[i,j,3,2+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumi1[i,j,2])
      yy_Sumi1[i,j,4,2+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumi1[i,j,2])*s_ITN[i+547] + (PHI_I - PHI_B)*(1-r_IRS_Sumi1[i,j,2])
      
      yy_Acte1[i,j,3,3+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Acte1[i,j,3])
      yy_Acte1[i,j,4,3+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Acte1[i,j,3])*s_ITN[i+547] + (PHI_I - PHI_B)*(1-r_IRS_Acte1[i,j,3])
      
      yy_Sumi1[i,j,3,3+cl[k]] = 1 - PHI_I + PHI_I*(1-r_IRS_Sumi1[i,j,3])
      yy_Sumi1[i,j,4,3+cl[k]] = 1 - PHI_I + PHI_B*(1-r_IRS_Sumi1[i,j,3])*s_ITN[i+547] + (PHI_I - PHI_B)*(1-r_IRS_Sumi1[i,j,3])
      
    }
  }  
}

for(j in 1:9){
  for(i in 1:365){
    yy_Acte[i,3,j] = sum(yy_Acte1[i,,3,j] * prop_this_weekM)
    yy_Acte[i,4,j] = sum(yy_Acte1[i,,4,j] * prop_this_weekM)
    
    yy_Sumi[i,3,j] = sum(yy_Sumi1[i,,3,j] * prop_this_weekB)
    yy_Sumi[i,4,j] = sum(yy_Sumi1[i,,4,j] * prop_this_weekB)
  }  
}


z_Acte[,1,] = 0
z_Sumi[,1,] = 0

for(k in 1:3){
  for(j in 1:18){
    for(i in 1:365){
      PHI_B = PHI_B_mut[k]
      PHI_I = PHI_I_mut[k]
      z_Acte1[i,j,2,1+cl[k]] = PHI_B*r_ITN[i+547]
      z_Acte1[i,j,3,1+cl[k]] = PHI_I*r_IRS_Acte1[i,j,1]
      z_Acte1[i,j,4,1+cl[k]] = PHI_B*(r_IRS_Acte1[i,j,1] + (1-r_IRS_Acte1[i,j,1])*r_ITN[i+547]) + (PHI_I - PHI_B)*r_IRS_Acte1[i,j,1]
      
      z_Sumi1[i,j,2,1+cl[k]] = PHI_B*r_ITN[i+547]
      z_Sumi1[i,j,3,1+cl[k]] = PHI_I*r_IRS_Sumi1[i,j,1]
      z_Sumi1[i,j,4,1+cl[k]] = PHI_B*(r_IRS_Sumi1[i,j,1] + (1-r_IRS_Sumi1[i,j,1])*r_ITN[i+547]) + (PHI_I - PHI_B)*r_IRS_Sumi1[i,j,1]
      
      
      z_Acte1[i,j,2,2+cl[k]] = PHI_B*r_ITN[i+547]
      z_Acte1[i,j,3,2+cl[k]] = PHI_I*r_IRS_Acte1[i,j,2]
      z_Acte1[i,j,4,2+cl[k]] = PHI_B*(r_IRS_Acte1[i,j,2] + (1-r_IRS_Acte1[i,j,2])*r_ITN[i+547]) + (PHI_I - PHI_B)*r_IRS_Acte1[i,j,2]
      
      z_Sumi1[i,j,2,2+cl[k]] = PHI_B*r_ITN[i+547]
      z_Sumi1[i,j,3,2+cl[k]] = PHI_I*r_IRS_Sumi1[i,j,2]
      z_Sumi1[i,j,4,2+cl[k]] = PHI_B*(r_IRS_Sumi1[i,j,2] + (1-r_IRS_Sumi1[i,j,2])*r_ITN[i+547]) + (PHI_I - PHI_B)*r_IRS_Sumi1[i,j,2]
      
      z_Acte1[i,j,2,3+cl[k]] = PHI_B*r_ITN[i+547]
      z_Acte1[i,j,3,3+cl[k]] = PHI_I*r_IRS_Acte1[i,j,3]
      z_Acte1[i,j,4,3+cl[k]] = PHI_B*(r_IRS_Acte1[i,j,3] + (1-r_IRS_Acte1[i,j,3])*r_ITN[i+547]) + (PHI_I - PHI_B)*r_IRS_Acte1[i,j,3]
      
      z_Sumi1[i,j,2,3+cl[k]] = PHI_B*r_ITN[i+547]
      z_Sumi1[i,j,3,3+cl[k]] = PHI_I*r_IRS_Sumi1[i,j,3]
      z_Sumi1[i,j,4,3+cl[k]] = PHI_B*(r_IRS_Sumi1[i,j,3] + (1-r_IRS_Sumi1[i,j,3])*r_ITN[i+547]) + (PHI_I - PHI_B)*r_IRS_Sumi1[i,j,3]
      
    }  
  }
  
}

for(j in 1:9){
  for(i in 1:365){
    z_Acte[i,2,j] = sum(z_Acte1[i,,3,j] * prop_this_weekM)
    z_Acte[i,3,j] = sum(z_Acte1[i,,3,j] * prop_this_weekM)
    z_Acte[i,4,j] = sum(z_Acte1[i,,4,j] * prop_this_weekM)
    
    z_Sumi[i,2,j] = sum(z_Sumi1[i,,3,j] * prop_this_weekB)
    z_Sumi[i,3,j] = sum(z_Sumi1[i,,3,j] * prop_this_weekB)
    z_Sumi[i,4,j] = sum(z_Sumi1[i,,4,j] * prop_this_weekB)
  }
  
}

## waning usage of IRS with time
## as well as altered cover over time from 3 month time line of campaign
prop_houses_sprayed_WeeklyB = 0.97*c(0,	0.027219794,	## August
                                     0.077014558,	0.136261919,	0.196901742,	0.250817066, ## Sept
                                     0.301047746,	0.347687015,	0.395348206,	0.464541108, ## Oct
                                     0.521818166,	0.581372602,	0.643327061,	0.710948828, ## Nov
                                     0.777620537,	0.847130239,	0.911498024,	0.930365931, ## Dec
                                     0.947042313,	0.96389906,	0.983400068,	0.991357264,	 ## Jan
                                     0.99238783,	1,1,1,## Feb
                                     rep(1,12)## mar-may
) 
prop_houses_sprayed_WeeklyM = 0.96*c(0.065358837,	0.193716785,## August
                                     0.334444596,	0.432141444,	0.533708885,	0.614060416,##sep
                                     0.667531537,	0.711462198,	0.769981823,	0.880751007,##oct
                                     0.919045204,	0.944027976,	0.97443187,	0.990922736,  ##nov
                                     0.991873307,	0.99744169,	 1,             1,          ##dec
                                     1,        1,1,1,## jan
                                     rep(1,16)) ## feb-may



## from Table 1 main manuscript 

## Accounting for modifications and added rooms in OCT cohort we have
cov_a_oct = c(0.891472868,0.776978417,0.556962025, 0.393063584,0.302702703, 0.218274112)
## Accounting for modifications and added rooms in NOV cohort we have
cov_a_nov = c(0.863636364,0.694736842,0.538461538,0.426086957,0.341269841,0.288888889)
## Accounting for modifications and added rooms in DEC cohort we have
cov_a_dec = c(0.851851852,0.714285714,0.655172414,0.6,0.566666667,0.548387097)

## Accounting for modifications and added rooms in OCT cohort we have
cov_s_oct = c(0.96460177,0.956140351,0.930434783,0.922413793,0.922413793,0.914529915)
## Accounting for modifications and added rooms in NOV cohort we have
cov_s_nov = c(0.973856209,0.948387097,0.918238994,0.918238994,0.900621118,0.895061728)
## Accounting for modifications and added rooms in DEC cohort we have
cov_s_dec = c(0.960526316,0.935064935,0.909090909,0.897435897,0.897435897,0.897435897)

## And weighting for the proportion of people represented in each survey 129:88:27 for Mut and 113:153:76 for Boane

# prop_mod_Acte = 1 - c(rep(0,10), ## august &  sep & oct
#                   rep(14/129,4), ## nov
#                   rep(14/129,4)+rep((12+7)/(117+88),4), ## dec
#                   rep(14/129,4)+rep((12+7)/(117+88),4)+rep((20+10+4)/(117+86+27),4), ## jan
#                   rep(14/129,4)+rep((12+7)/(117+88),4)+rep((20+10+4)/(117+86+27),4)+rep((20+10+3)/(116+85+25),4), ## feb
#                   rep(14/129,4)+rep((12+7)/(117+88),4)+rep((20+10+4)/(117+86+27),4)+rep((20+10+3)/(116+85+25),4)+rep((12+7+1)/(115+85+25),4), ## mar
#                   rep(14/129,4)+rep((12+7)/(117+88),4)+rep((20+10+4)/(117+86+27),4)+rep((20+10+3)/(116+85+25),4)+rep((12+7+1)/(115+85+25),4)+rep((16+6+1)/(84+25),4), ## apr
#                   rep(14/129,4)+rep((12+7)/(117+88),4)+rep((20+10+4)/(117+86+27),4)+rep((20+10+3)/(116+85+25),4)+rep((12+7+1)/(115+85+25),4)+rep((16+6+1)/(84+25),4)+rep((4+1)/(81+25),4)) ## may
prop_mod_Acte_temp = c(0.891472868,0.81212081,0.639282555,0.481047262,0.386204738,0.304873387,0.166894353)

##And sorted for weeks:
prop_mod_Acte = c(rep(1,10),rep(prop_mod_Acte_temp,each=4))## for time series of spray campaign

true_cover_irs_Acte =  rep(prop_mod_Acte*prop_houses_sprayed_WeeklyM,each=7)

#House coverage: Matutuine district 96 %
irs_cov_no_loss_Acte = rep(0.96,30*8)
irs_cov_Acte = true_cover_irs_Acte

prop_mod_Sumi_temp = c(0.96460177,0.966330299,0.945153088,0.923357485,0.917585479,0.904508888,0.895619142)

##And sorted for weeks:
prop_mod_Sumi = c(rep(1,10),rep(prop_mod_Sumi_temp,each=4))## for time series of spray campaign

true_cover_irs_Sumi = rep(prop_mod_Sumi*prop_houses_sprayed_WeeklyB,each=7)

#House coverage: Boane (sumi) district 97 %, Manhica district (Palmeira) 98 % 
irs_cov_no_loss_Sumi = rep(0.97,30*8)
irs_cov_Sumi = true_cover_irs_Sumi

plot(irs_cov_no_loss_Acte[1:240] ~ time[1:240],ylab = "Community IRS cover (%)",
     ylim=c(0,1),col="black",pch="",
     main = "",cex.main=1.2,xlim=c(1,240),xaxt="n",
     xlab="Time in months",yaxt="n",cex.lab=1.2,cex.axis=1.2,cex=1.2)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.2,cex.axis=1.2)
axis(1,at=seq(0,230,30)+15,labels = c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),cex.axis = 1.4)


lines(irs_cov_no_loss_Acte[61:240] ~ time[61:240],lty=1,lwd=2,col = "darkblue")
lines(irs_cov_Acte[61:240] ~ time[61:240],lty=3,lwd=2,col = "darkblue")

lines(irs_cov_no_loss_Sumi[61:240] ~ c(time[61:240]+1),lty=1,lwd=2,col = "aquamarine3")
lines(irs_cov_Sumi[61:240] ~ time[61:240],lty=3,lwd=2,col = "aquamarine3")

cov1A = cov1S = cov2A = cov2S = array(dim=c(240,4))

# itn_cov_Acte = 0.52
# itn_cov_Sumi = 0.73

## Table 1 data
matu_net_cov = c(27.9,  ## nov
                 mean(c(32.6,42.1)), ##dec
                 mean(c(43.4,43.2,33.3)),##jan
                 mean(c(55.0,44.3,48.2)),##feb
                 mean(c(61.2,56.8,40.7)),##mar
                 mean(c(51.9,61.4,66.7)),##apr
                 mean(c(48.9,55.6)),##may
                 29.6)##june

itn_cov_Acte_temp = c(rep(mean(c(matu_net_cov)),10),
                      rep(mean(c(32.6,42.1)),4), ##dec
                      rep(mean(c(43.4,43.2,33.3)),4),##jan
                      rep(mean(c(55.0,44.3,48.2)),4),##feb
                      rep(mean(c(61.2,56.8,40.7)),4),##mar
                      rep(mean(c(51.9,61.4,66.7)),4),##apr
                      rep(mean(c(48.9,55.6)),4),##may
                      rep(29.6,4))##june



boan_net_cov = c(57.5,  ## nov
                 mean(c(64.6,74.5)), ##dec
                 mean(c(68.1,71.9,67.1)),##jan
                 mean(c(62.8,71.2,79)),##feb
                 mean(c(70.8,82.4,86.8)),##mar
                 mean(c(65.5,81.8,81.6)),##apr
                 mean(c(78.4,84.2)),##may
                 10.5)##june

itn_cov_Boan_temp = c(rep(mean(c(boan_net_cov)),10),
                      rep(mean(c(64.6,74.5)),4), ##dec
                      rep(mean(c(68.1,71.9,67.1)),4),##jan
                      rep(mean(c(62.8,71.2,79)),4),##feb
                      rep(mean(c(70.8,82.4,86.8)),4),##mar
                      rep(mean(c(65.5,81.8,81.6)),4),##apr
                      rep(mean(c(78.4,84.2)),4),##may
                      rep(10.5,4))

itn_cov_Acte = rep(itn_cov_Acte_temp,each=7)/100
itn_cov_Sumi = rep(itn_cov_Boan_temp,each=7)/100

## Here we are creating a matrix
## with the coverage or use of nets waning with time
## and the coverage of IRS either staying fixed, or also waning 
## when walls are washed
cov1A[,1] = 1
cov1A[,2] = itn_cov_Acte[1:240] ## ITN only
cov1A[,3] = irs_cov_no_loss_Acte ## IRS only
cov1A[,4] = itn_cov_Acte[1:240]*irs_cov_no_loss_Acte ## both interventions

cov2A[,1] = 1
cov2A[,2] = itn_cov_Acte[1:240] ## ITN only
cov2A[,3] = irs_cov_Acte[1:240] ## IRS only
cov2A[,4] = itn_cov_Acte[1:240]*irs_cov_Acte[1:240]## both interventions

cov1S[,1] = 1
cov1S[,2] = itn_cov_Sumi[1:240] ## ITN only
cov1S[,3] = irs_cov_no_loss_Sumi ## IRS only
cov1S[,4] = itn_cov_Sumi[1:240]*irs_cov_no_loss_Sumi ## both interventions

cov2S[,1] = 1
cov2S[,2] = itn_cov_Sumi[1:240] ## ITN only
cov2S[,3] = irs_cov_Sumi[1:240] ## IRS only
cov2S[,4] = itn_cov_Sumi[1:240]*irs_cov_Sumi[1:240] ## both interventions

## Entomological model parameters to estimate 

Q0 = 0.98  ## this is anthropophagy - we can use human blood index
chi = 0.86 ## this is endophily (pi_i)

fv0 = 0.333 ## biting rate 1 bite every 3 days
tau1 = 0.69 ## duration of host seeking, assumed to be constant between species (delta_10, altered to delta_1 when interventions used)
tau2 = 1/fv0-tau1  ## indoor feeding endophagy (delta_2)
av0 = Q0*fv0
mu0 = 0.132 ## background mortality from external sources
p10 = exp(-mu0*tau1) ## 
p2 = exp(-mu0*tau2)  ## probability of surviving resting period in absence of intervntion


## These are the adjusted w, z, when coverage is changing
## so these are the intervention coverages
zhi1A = whi1A = zhi2A = whi2A = array(dim=c(240,4,9))
zhi1S = whi1S = zhi2S = whi2S = array(dim=c(240,4,9))
for(i in 1:9){
  zhi1A[,,i]=cov1A*z_Acte[1:240,,i]
  whi1A[,,i]=cov1A*w_Acte[1:240,,i]
  
  zhi1S[,,i]=cov1S*z_Sumi[1:240,,i]
  whi1S[,,i]=cov1S*w_Sumi[1:240,,i]
  
  
  zhi2A[,,i]=cov2A*z_Acte[1:240,,i]
  whi2A[,,i]=cov2A*w_Acte[1:240,,i]
  
  zhi2S[,,i]=cov2S*z_Sumi[1:240,,i]
  whi2S[,,i]=cov2S*w_Sumi[1:240,,i]
  
}



zbar1A = wbar1A = zbar2A = wbar2A = array(dim=c(240,4,9)) 
zbar1S = wbar1S = zbar2S = wbar2S = array(dim=c(240,4,9))
for(j in 1:9){
  for(i in 1:4){
    zbar1A[,i,j] = Q0*zhi1A[,i,j]
    wbar1A[,i,j] = (1 - Q0) + Q0*whi1A[,i,j]
    zbar2A[,i,j] = Q0*zhi2A[,i,j]
    wbar2A[,i,j] = (1 - Q0) + Q0*whi2A[,i,j]
    
    zbar1S[,i,j] = Q0*zhi1S[,i,j]
    wbar1S[,i,j] = (1 - Q0) + Q0*whi1S[,i,j]
    zbar2S[,i,j] = Q0*zhi2S[,i,j]
    wbar2S[,i,j] = (1 - Q0) + Q0*whi2S[,i,j]
  }  
}


## From Walker et al 2016
## Mosquito feeding rate (tau1 is delta10, tau2 is delta2 in the methods)
fR1A = 1 / ((tau1/(1 - zbar1A)) + tau2)
mu1A = -fR1A*log((wbar1A*p10/(1 - zbar1A*p10))*p2) 
Q1A = 1 - (1-Q0)/wbar1A

fR2A = 1 / ((tau1/(1 - zbar2A)) + tau2)
mu2A = -fR2A*log((wbar2A*p10/(1 - zbar2A*p10))*p2) 
Q2A = 1 - (1-Q0)/wbar2A


fR1S = 1 / ((tau1/(1 - zbar1S)) + tau2)
mu1S = -fR1S*log((wbar1S*p10/(1 - zbar1S*p10))*p2) 
Q1S = 1 - (1-Q0)/wbar1S

fR2S = 1 / ((tau1/(1 - zbar2S)) + tau2)
mu2S = -fR2S*log((wbar2S*p10/(1 - zbar2S*p10))*p2) 
Q2S = 1 - (1-Q0)/wbar2S


## Rate at which a person in the popn is bitten by mosquitoes is
lambda1A = lambda2A = array(dim = c(240,4,9))
lambda1S = lambda2S = array(dim = c(240,4,9))
for(j in 1:9){
  for(i in 1:4){
    lambda1A[,,j] = (Q1A[,,j]*fR1A[,,j]*yy_Acte[1:240,i,j])/whi1A[1:240,i,j]
    lambda2A[,,j] = (Q2A[,,j]*fR2A[,,j]*yy_Acte[1:240,i,j])/whi2A[1:240,i,j]
    
    lambda1S[,,j] = (Q1S[,,j]*fR1S[,,j]*yy_Sumi[1:240,i,j])/whi1S[1:240,i,j]
    lambda2S[,,j] = (Q2S[,,j]*fR2S[,,j]*yy_Sumi[1:240,i,j])/whi2S[1:240,i,j]
  }  
}




## Additional infectious bites per person per year 
Estimated_propn_increase_EIR = array(dim=c(240,2,9))
Estimated_propn_increase_EIR[,1,] = (lambda2A[,4,] - lambda1A[,4,])/lambda2A[,4,]
Estimated_propn_increase_EIR[,2,] = (lambda2S[,4,] - lambda1S[,4,])/lambda2S[,4,]

c(sum(Estimated_propn_increase_EIR[1:30,1,])/30,sum(Estimated_propn_increase_EIR[31:60,1,])/30,
  sum(Estimated_propn_increase_EIR[61:90,1,])/30,sum(Estimated_propn_increase_EIR[91:120,1,])/30,
  sum(Estimated_propn_increase_EIR[121:150,1,])/30,sum(Estimated_propn_increase_EIR[151:180,1,])/30)

c(sum(Estimated_propn_increase_EIR[1:30,2,])/30,sum(Estimated_propn_increase_EIR[31:60,2,])/30,
  sum(Estimated_propn_increase_EIR[61:90,2,])/30,sum(Estimated_propn_increase_EIR[91:120,2,])/30,
  sum(Estimated_propn_increase_EIR[121:150,2,])/30,sum(Estimated_propn_increase_EIR[151:180,2,])/30)

mean(Estimated_propn_increase_EIR[15:45,1,])##sep
mean(Estimated_propn_increase_EIR[46:75,1,])##oct
mean(Estimated_propn_increase_EIR[76:105,1,])##nov
mean(Estimated_propn_increase_EIR[106:135,1,])##dec
mean(Estimated_propn_increase_EIR[136:165,1,])##jan
mean(Estimated_propn_increase_EIR[166:195,1,])##feb
mean(Estimated_propn_increase_EIR[196:225,1,])##mar
mean(Estimated_propn_increase_EIR[226:240,1,])##part april

mean(Estimated_propn_increase_EIR[15:45,2,])##sep
mean(Estimated_propn_increase_EIR[46:75,2,])##oct
mean(Estimated_propn_increase_EIR[76:105,2,])##nov
mean(Estimated_propn_increase_EIR[106:135,2,])##dec
mean(Estimated_propn_increase_EIR[136:165,2,])##jan
mean(Estimated_propn_increase_EIR[166:195,2,])##feb
mean(Estimated_propn_increase_EIR[196:225,2,])##mar
mean(Estimated_propn_increase_EIR[226:240,2,])##part april


plot(Estimated_propn_increase_EIR[,1,1] ~ time[1:240],ylim=c(0,1),pch="",
     ylab = "",
     col="black",
     main = "",cex.main=1.2,xlim=c(1,240),xaxt="n",
     xlab="Time in months",yaxt="n",cex.lab=1.2,cex.axis=1.2,cex=1.2)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.lab=1.2,cex.axis=1.2)
axis(1,at=seq(0,230,30)+15,labels = c("Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr"),cex.axis = 1.4)
mtext(side=2, line =4.3,
      "Relative increase in daily bites")
mtext(side=2, line =3, "due to spray campaign & modifications (%)")

colsd = c("darkblue","aquamarine3")
for(i in 1:2){
  lines(Estimated_propn_increase_EIR[61:240,i,1] ~ time[61:240],col=colsd[i],lty=2,lwd=2)
}

min_EIR = max_EIR = array(dim=c(length(61:240),2))
for(j in 1:2){
  for(i in 1:length(61:240)){
    min_EIR[i,j] = min(Estimated_propn_increase_EIR[60+i,j,])
    max_EIR[i,j] = max(Estimated_propn_increase_EIR[60+i,j,])
    
  }  
}

polygon(c(time[61:240],rev(time[61:240])),
        c(min_EIR[,1],rev(max_EIR[,1])),
        col=adegenet::transp("darkblue",0.3),border=NA)
polygon(c(time[61:240],rev(time[61:240])),
        c(min_EIR[,2],rev(max_EIR[,2])),
        col=adegenet::transp("aquamarine",0.3),border=NA)



Estimated_added_EIR = array(dim=c(240,2,9))
Estimated_added_EIR[,1,] = (lambda2A[,3,] - lambda1A[,3,])
Estimated_added_EIR[,2,] = (lambda2S[,3,] - lambda1S[,3,])

## Additional infectious bites per person per year 
Estimated_propn_increase_EIR = array(dim=c(240,2,9))
Estimated_propn_increase_EIR[,1,] = (lambda2A[,3,] - lambda1A[,3,])/lambda2A[,3,]
Estimated_propn_increase_EIR[,2,] = (lambda2S[,3,] - lambda1S[,3,])/lambda2S[,3,]

for(i in 1:2){
  lines(Estimated_propn_increase_EIR[61:240,i,1] ~ time[61:240],col=colsd[i],lty=1,lwd=1)
}

min_EIR = max_EIR = array(dim=c(length(61:240),2))
for(j in 1:2){
  for(i in 1:length(61:240)){
    min_EIR[i,j] = min(Estimated_propn_increase_EIR[60+i,j,])
    max_EIR[i,j] = max(Estimated_propn_increase_EIR[60+i,j,])
    
  }  
}

polygon(c(time[61:240],rev(time[61:240])),
        c(min_EIR[,1],rev(max_EIR[,1])),
        col=adegenet::transp("darkblue",0.3),border=NA)
polygon(c(time[61:240],rev(time[61:240])),
        c(min_EIR[,2],rev(max_EIR[,2])),
        col=adegenet::transp("aquamarine",0.3),border=NA)


