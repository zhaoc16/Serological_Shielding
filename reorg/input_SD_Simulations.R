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
model_pars[['nDays']]=365 # days
model_pars[['times']]=0:model_pars$nDays

# Indices
model_pars[['subgroups']] = c('c', 'a', 'rc', 'fc', 'e')
model_pars[['nSubgroups']] = length(model_pars$subgroups) #ALWAYS: c, a, rc, fc, e
N=model_pars$nSubgroups

model_pars[['compartments']] = c('S', 'E', 'Isym', 'Iasym', 'Hsub', 'Hcri', 'D', 'R', 'S_pos', 'E_pos', 'Isym_pos', 'Iasym_pos', 'R_pos')
model_pars[['nCompartments']]=length(model_pars$compartments) #ALWAYS: S, E, Isym, Iasym, Hsub, Hcri, D, R, S_pos, E_pos, Isym_pos, Iasym_pos, R_pos
Nc=model_pars$nCompartments

# Variable Names
model_pars[['varNames']] = Get_Var_Names(model_pars$subgroups, model_pars$compartments)

temp_idxNames = 1:(N*Nc)
names(temp_idxNames) = model_pars$varNames
model_pars[['idxNames']] = temp_idxNames

# Variables in Matrix Format
temp.varMat = matrix(model_pars$varNames, ncol = Nc)
colnames(temp.varMat) = model_pars$compartments
rownames(temp.varMat) = model_pars$subgroups
model_pars[['varMat']] = temp.varMat #columns are compartments, rows are population subclasses

temp.idxMat = matrix(model_pars$idxNames, ncol = Nc)
colnames(temp.idxMat) = model_pars$compartments
rownames(temp.idxMat) = model_pars$subgroups
model_pars[['idxMat']] = temp.idxMat #columns are compartments, rows are population subclasses

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
intervention_pars[['sensitivity']]=1.00
intervention_pars[['specificity']]=0.998
intervention_pars[['daily_tests']]=10**3

# Other intervention parameters
intervention_pars[['tStart_distancing']]=70
intervention_pars[['tStart_test']]=107    # can change when in the outbreak testing becomes available
intervention_pars[['tStart_target']]=115
intervention_pars[['tStart_school']]=230
intervention_pars[['tStart_reopen']]=500

intervention_pars[['socialDistancing_other']]=0.25
intervention_pars[['p_reduced']]=0.5      # proportion of contacts reduced 
intervention_pars[['p_full']]=1           # proportion of contacts reduced for full contact adults

intervention_pars[['alpha']]=1.2          # shielding. Note this is not alpha_JSW, but (alpha_JSW+1)
intervention_pars[['c']]=1

intervention_pars[['socialDistancing_other_c']]=0.25
intervention_pars[['p_reduced_c']]=0.1


# (4) Epi Pars ------------------------------------------------------------
epi_pars=list()

epi_pars[['R0']]=3.1                # Note on R0: with base structure 63.28q
epi_pars[['q']]=epi_pars$R0/79.27   # Probability of transmission from children
epi_pars[['asymp_red']]=0.55        # Relative infectiousness of asymptomatic vs symptomatic case

epi_pars[['gamma_e']]=1/3           # Latent period (He et al)
epi_pars[['gamma_a']]=1/7           # Recovery rate, undocumented (Kissler et al)
epi_pars[['gamma_s']]=1/7           # Recovery rate, undocumented (Kissler et al)
epi_pars[['gamma_hs']]=1/5          # LOS for subcritical cases (medrxiv paper)
epi_pars[['gamma_hc']]=1/7          # LOS for critical cases (medrxiv paper)
epi_pars[['p']]=0.5#0.14            # Fraction 'Symptomatic'documented' (Shaman's paper)

epi_pars[['hosp_frac']]=c(0.002, 0.056, 0.224)  # From MMWR
epi_pars[['hosp_crit']]=c(0.001, 0.0048, 0.099) # From CDC, MMWR
epi_pars[['crit_die']]=c(0, 0.5, 0.5)           # Obtained from initial fitting

epi_pars[['hosp_frac_5']]=epi_pars$hosp_frac[c(1,2,2,2,3)]
epi_pars[['hosp_crit_5']]=epi_pars$hosp_crit[c(1,2,2,2,3)]
epi_pars[['crit_die_5']]=epi_pars$crit_die[c(1,2,2,2,3)]  

# Overwrite Using Ferguson Parms


# (5) Inits ---------------------------------------------------------------

# Initial conditions
inits=list()

inits[['I_c0']]=60
inits[['I_a0']]=20
inits[['I_rc0']]=50
inits[['I_fc0']]=1
inits[['I_e0']]=40


# (6) Remainder -----------------------------------------------------------

#Might want a ramp up period for these, ignore for now
tswitch1.dat=data.frame(times=model_pars$times, test.switch1=c(rep(1, 366)))
tswitch2.dat=data.frame(times=model_pars$times, test.switch2=c(rep(0, 366)))

sw1fxn=approxfun(tswitch1.dat$times, tswitch1.dat$test.switch1, rule=2)
sw2fxn=approxfun(tswitch2.dat$times, tswitch2.dat$test.switch2, rule=2)

# Combine
# pars_default = list('model_pars' = model_pars
#                     , 'contact_pars' = contact_pars
#                     , 'intervention_pars'= intervention_pars
#                     , 'epi_pars' = epi_pars)

pars_default = c(model_pars, contact_pars, intervention_pars, epi_pars, inits)


# (7) Test Cases ----------------------------------------------------------

X0=rep(0,pars_default$nTotSubComp)

# Remove from susceptible pool ...
X0[1]=pars_default$agestruc[1]*pars_default$N-inits$I_c0 
X0[2]=pars_default$agestruc[2]*pars_default$frac_home*pars_default$N-inits$I_a0 
X0[3]=pars_default$agestruc[2]*pars_default$frac_reduced*pars_default$N-inits$I_rc0 
X0[4]=pars_default$agestruc[2]*pars_default$frac_full*pars_default$N-inits$I_fc0 
X0[5]=pars_default$agestruc[3]*pars_default$N-inits$I_e0 

# ... and add to infected
X0[11]=inits$I_c0
X0[12]=inits$I_a0
X0[13]=inits$I_rc0
X0[14]=inits$I_fc0
X0[15]=inits$I_e0

names(X0) = pars_default$varNames


# 
# # Test
# set.seed(1234)
# Xtest = round(runif(length(X0))*pars_default$N)
# names(Xtest) = names(X0)
# 
# seir_model_shields_rcfc_nolatent(0, X0, pars_default)
