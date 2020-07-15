#General setup
rm(list=ls()) #Clear workspace
require(deSolve)
require(ggplot2)
require(reshape2)
require(cowplot)
require(here)
require(plyr)
require(beepr)

source('reorg/utils.R')


# (1) Model Parameters ----------------------------------------------------
model_pars=list()

# Population
model_pars[['N']]=323*10^6
model_pars[['agefrac.0']]=c(0.12,0.13,0.13,0.13,0.13,0.13,0.11,0.06,0.04,0.02) # from Weitz model
model_pars[['agestruc']]=c(sum(model_pars$agefrac.0[1:2])
                         , sum(model_pars$agefrac.0[3:6], 0.5*model_pars$agefrac.0[7])
                         , sum(0.5*model_pars$agefrac.0[7], model_pars$agefrac.0[8:10]))

# Timeline
model_pars[['nDays']]=366 # days
model_pars[['times']]=1:model_pars$nDays

# Indices
model_pars[['nSubgroups']]=5 #ALWAYS: c, a, rc, fc, e
N=model_pars$nSubgroups

model_pars[['nCompartments']]=13 #ALWAYS: S, E, Isym, Iasym, Hsub, Hcri, D, R, S_pos, E_pos, Isym_pos, Iasym_pos, R_pos
Nc=model_pars$nCompartments

model_pars[['varNames']] = c('S_c', 'S_a', 'S_rc', 'S_fc', 'S_e'
                             , 'E_c', 'E_a', 'E_rc', 'E_fc', 'E_e'
                             , 'Isym_c', 'Isym_a', 'Isym_rc', 'Isym_fc', 'Isym_e'
                             , 'Iasym_c', 'Iasym_a', 'Iasym_rc', 'Iasym_fc', 'Iasym_e'
                             , 'Hsub_c', 'Hsub_a', 'Hsub_rc', 'Hsub_fc', 'Hsub_e'
                             , 'Hcri_c', 'Hcri_a', 'Hcri_rc', 'Hcri_fc', 'Hcri_e'
                             , 'D_c', 'D_a', 'D_rc', 'D_fc', 'D_e'
                             , 'R_c', 'R_a', 'R_rc', 'R_fc', 'R_e'
                             , 'S_c_pos', 'S_a_pos', 'S_rc_pos', 'S_fc_pos', 'S_e_pos'
                             , 'E_c_pos', 'E_a_pos', 'E_rc_pos', 'E_fc_pos', 'E_e_pos'
                             , 'Isym_c_pos', 'Isym_a_pos', 'Isym_rc_pos', 'Isym_fc_pos', 'Isym_e_pos'
                             , 'Iasym_c_pos', 'Iasym_a_pos', 'Iasym_rc_pos', 'Iasym_fc_pos', 'Iasym_e_pos'
                             , 'R_c_pos', 'R_a_pos', 'R_rc_pos', 'R_fc_pos', 'R_e_pos')
temp_idxNames = 1:(N*Nc)
names(temp_idxNames) = model_pars$varNames
model_pars[['idxNames']] = temp_idxNames

# Base Compartment Indices
model_pars[['S_ids']]=0*N+(1:N)
model_pars[['E_ids']]=1*N+(1:N)
model_pars[['Isym_ids']]=2*N+(1:N)
model_pars[['Iasym_ids']]=3*N+(1:N)
model_pars[['Hsub_ids']]=4*N+(1:N)
model_pars[['Hcri_ids']]=5*N+(1:N)
model_pars[['D_ids']]=6*N+(1:N)
model_pars[['R_ids']]=7*N+(1:N)

# Positive Test Compartment Indices
model_pars[['S_pos_ids']]=8*N+(1:N)
model_pars[['E_pos_ids']]=9*N+(1:N)
model_pars[['Isym_pos_ids']]=10*N+(1:N)
model_pars[['Iasym_pos_ids']]=11*N+(1:N)
model_pars[['R_pos_ids']]=12*N+(1:N)


# IDs by SubGroup
model_pars[['c_ids']]=(1:Nc-1)*N+1
model_pars[['a_ids']]=(1:Nc-1)*N+2
model_pars[['rc_ids']]=(1:Nc-1)*N+3
model_pars[['fc_ids']]=(1:Nc-1)*N+4
model_pars[['e_ids']]=(1:Nc-1)*N+5

model_pars[['nTotSubComp']]=N*Nc


# (2) Contact Parameters --------------------------------------------------
contact_pars=list()

contact_pars[['frac_home']]=0.316
contact_pars[['frac_full']]=0.0565
  #contact_pars[['prob.home']]=0.4
  #contact_pars[['prob.reduced']]==0.678
contact_pars[['frac_reduced']]=1-contact_pars$frac_full-contact_pars$frac_home

#3x3 data from Prem et al
contact_pars[['AllContacts']]=matrix(c(9.75, 2.57, 0.82, 5.97, 10.32, 2.25, 0.39, 0.46, 1.20), nrow=3)
contact_pars[['WorkContacts']]=matrix(c(0.20, 0.28, 0, 0.64, 4.73, 0, 0, 0, 0), nrow=3)
contact_pars[['SchoolContacts']]=matrix(c(4.32, 0.47, 0.02, 1.10, 0.32, 0.04, 0.01, 0.01, 0.03), nrow=3)
contact_pars[['HomeContacts']]=matrix(c(2.03, 1.02, 0.50, 2.37, 1.82, 0.68, 0.24, 0.14, 0.62), nrow=3)
contact_pars[['OtherContacts']]=matrix(c(3.20, 0.80, 0.30, 1.86, 3.45, 1.53, 0.14, 0.32, 0.55), nrow=3)

# 5x5 Expansions
contact_pars[['WorkContacts_5x5']]=Expand_5x5(contactMatrix=contact_pars$WorkContacts
                                              , p_Home=contact_pars$frac_home
                                              , p_Reduced=contact_pars$frac_reduced
                                              , p_Full=contact_pars$frac_full)
contact_pars[['HomeContacts_5x5']]=Expand_5x5(contactMatrix=contact_pars$HomeContacts
                                              , p_Home=contact_pars$frac_home
                                              , p_Reduced=contact_pars$frac_reduced
                                              , p_Full=contact_pars$frac_full)
contact_pars[['SchoolContacts_5x5']]=Expand_5x5(contactMatrix=contact_pars$SchoolContacts
                                                , p_Home=contact_pars$frac_home
                                                , p_Reduced=contact_pars$frac_reduced
                                                , p_Full=contact_pars$frac_full)
contact_pars[['OtherContacts_5x5']]=Expand_5x5(contactMatrix=contact_pars$OtherContacts
                                               , p_Home=contact_pars$frac_home
                                               , p_Reduced=contact_pars$frac_reduced
                                               , p_Full=contact_pars$frac_full)


# (3) Intervention Pars ---------------------------------------------------
intervention_pars=list()

# Test Data for Cellex (from April 2020)
intervention_pars[['sensitivity']]=0.94
intervention_pars[['specificity']]=0.96
intervention_pars[['daily_tests']]=10**3

# Other intervention parameters
intervention_pars[['tStart_distancing']]=70
intervention_pars[['tStart_target']]=115
intervention_pars[['tStart_school']]=230

intervention_pars[['socialDistancing_other']]=.25
intervention_pars[['p_reduced']]=0.1 # proportion of contacts reduced 
intervention_pars[['p_full']]=1 # proportion of contacts reduced for full contact adults

intervention_pars[['alpha']]=4 # shielding. Note this is not alpha_JSW, but (alpha_JSW+1)
intervention_pars[['c']]=1

intervention_pars[['socialDistancing_other_c']]=1-(0.75*intervention_pars$c)
intervention_pars[['p_reduced_c']]=1-(0.9*intervention_pars$c)

intervention_pars[['tStart_test']]=107 # can change when in the outbreak testing becomes available
intervention_pars[['tStart_target']]=115


# (4) Epi Pars ------------------------------------------------------------
epi_pars=list()

epi_pars[['R0']]=2.9                # Note on R0: with base structure 63.28q
epi_pars[['q']]=epi_pars$R0/63.28   # Probability of transmission from children
epi_pars[['asymp.red']]=0.55        # Relative infectiousness of asymptomatic vs symptomatic case

epi_pars[['gamma_e']]=1/3           # Latent period (He et al)
epi_pars[['gamma_a']]=1/7           # Recovery rate, undocumented (Kissler et al)
epi_pars[['gamma_s']]=1/7           # Recovery rate, undocumented (Kissler et al)
epi_pars[['gamma_h']]=1/15          # Recovery rate, hospitalized cases (Zhou et al--China study)
epi_pars[['p']]=0.3#0.14            # Fraction 'Symptomatic'documented' (Shaman's paper)

epi_pars[['hosp_frac']]=c(0.061, 0.182, 0.417)  # From MMWR
epi_pars[['hosp_crit']]=c(0, 0.063, 0.173)      # From CDC, MMWR
epi_pars[['crit_die']]=c(0, 0.5, 0.5)           # Obtained from initial fitting

epi_pars[['hosp_frac_5']]=epi_pars$hosp_frac[c(1,2,2,2,3)]
epi_pars[['hosp_crit_5']]=epi_pars$hosp_crit[c(1,2,2,2,3)]
epi_pars[['crit_die_5']]=epi_pars$crit_die[c(1,2,2,2,3)]  


# (5) Inits ---------------------------------------------------------------

# Initial conditions
inits=list()

inits[['I_c0']]=60
inits[['I_a0']]=20
inits[['I_rc0']]=50
inits[['I_fc0']]=1
inits[['I_e0']]=40

X0=rep(0,model_pars$nTotSubComp)

# Remove from susceptible pool ...
X0[1]=model_pars$agestruc[1]*model_pars$N-inits$I_c0 
X0[2]=model_pars$agestruc[2]*contact_pars$frac_home*model_pars$N-inits$I_a0 
X0[3]=model_pars$agestruc[2]*contact_pars$frac_reduced*model_pars$N-inits$I_rc0 
X0[4]=model_pars$agestruc[2]*contact_pars$frac_full*model_pars$N-inits$I_fc0 
X0[5]=model_pars$agestruc[3]*model_pars$N-inits$I_e0 

# ... and add to infected
X0[11]=inits$I_c0
X0[12]=inits$I_a0
X0[13]=inits$I_rc0
X0[14]=inits$I_fc0
X0[15]=inits$I_e0

names(X0) = model_pars$varNames


# (6) Remainder -----------------------------------------------------------

#Might want a ramp up period for these, ignore for now
tswitch1.dat=data.frame(times=model_pars$times, test.switch1=c(rep(1, 366)))
tswitch2.dat=data.frame(times=model_pars$times, test.switch2=c(rep(0, 366)))

sw1fxn=approxfun(tswitch1.dat$times, tswitch1.dat$test.switch1, rule=2)
sw2fxn=approxfun(tswitch2.dat$times, tswitch2.dat$test.switch2, rule=2)

# Combine
pars = c(model_pars, contact_pars, intervention_pars, epi_pars)
