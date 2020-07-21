setwd('/Users/czhao/Documents/CYZ GITHUB/Weitz Group/COVID-19/Lopman COVID collab/Serological_Shielding/')

source('CYZ_reorg/R/input_SD_Simulations.R')

# Get Pars
pars_doNothing = pars_default

pars_doNothing[['socialDistancing_other']]=0.25
pars_doNothing[['socialDistancing_other_c']]=0.25
pars_doNothing[['p_reduced']]=0.1
pars_doNothing[['p_reduced_c']]=0.1
pars_doNothing[['daily_tests']]=0

pars_doNothing[['tStart_distancing']]=500
pars_doNothing[['tStart_test']]=500
pars_doNothing[['tStart_target']]=500
pars_doNothing[['tStart_school']]=500
pars_doNothing[['tStart_reopen']]=500

# Run Model
model_out = ode(y = X0
                , times = pars_doNothing$times
                , fun = seir_model_shields_rcfc_nolatent
                , parms = pars_doNothing
                , method='ode45')
model_out = as.data.frame(model_out)[,-1]

Infected=pars_doNothing$N - sum(model_out[366, c(pars_doNothing$S_ids, pars_doNothing$S_pos_ids)])
Pinfect=Infected/pars_doNothing$N
Infected; Pinfect

Deaths=sum(model_out[366, pars_doNothing$D_ids])
Deaths
