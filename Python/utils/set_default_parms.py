import numpy as np
from utils.helpers import Expand_Contact_Matrices, Expand_10x10

# pars ----------------------------------------------------------------------------------------
pars = {}

# Timeline
pars['nDays'] = 366 # days

# Population Structure
pars['N'] = 323*10**6
pars['agefrac_0'] = np.array([0.12,0.13,0.13,0.13,0.13,0.13,0.11,0.06,0.04,0.02])
pars['agestruc'] = np.array((np.sum(pars['agefrac_0'][0:2])
                            , np.sum(pars['agefrac_0'][2:6]) + 0.5*pars['agefrac_0'][6]
                            , 0.5*pars['agefrac_0'][6] + np.sum(pars['agefrac_0'][7:10])))

# IDs
pars['nSubgroups'] = 5 #ALWAYS: c, a, rc, fc, e
N = pars['nSubgroups']

pars['nCompartments'] = 13 #ALWAYS: S, E, Isym, Iasym, Hsub, Hcri, D, R, S_pos, E_pos, Isym_pos, Iasym_pos, R_pos
Nc = pars['nCompartments']

# Base Compartment Indices
pars['S_ids'] = np.arange(start = (0*N), stop = (1*N))
pars['E_ids'] = np.arange(start = (1*N), stop = (2*N))
pars['Isym_ids'] = np.arange(start = (2*N), stop = (3*N))
pars['Iasym_ids'] = np.arange(start = (3*N), stop = (4*N))
pars['Hsub_ids'] = np.arange(start = (4*N), stop = (5*N))
pars['Hcri_ids'] = np.arange(start = (5*N), stop = (6*N))
pars['D_ids'] = np.arange(start = (6*N), stop = (7*N))
pars['R_ids'] = np.arange(start = (7*N), stop = (8*N))

# Positive Test Compartment Indices
pars['S_pos_ids'] = np.arange(start = (8*N), stop = (9*N))
pars['E_pos_ids'] = np.arange(start = (9*N), stop = (10*N))
pars['Isym_pos_ids'] = np.arange(start = (10*N), stop = (11*N))
pars['Iasym_pos_ids'] = np.arange(start = (11*N), stop = (12*N))
pars['R_pos_ids'] = np.arange(start = (12*N), stop = (13*N))


# IDs by SubGroup
pars['c_ids'] = np.arange(Nc)*N + 0
pars['a_ids'] = np.arange(Nc)*N + 1
pars['rc_ids'] = np.arange(Nc)*N + 2
pars['fc_ids'] = np.arange(Nc)*N + 3
pars['e_ids'] = np.arange(Nc)*N + 4

pars['nTotSubComp'] = N * Nc


# Epi Pars -----------------------------------------------------------------------------------------
epiPars = {}

epiPars['R0'] = 2.9
epiPars['q'] = epiPars['R0']/63.28 # Probability of transmission from children
epiPars['asymp_red'] = 0.55        # Relative infectiousness of asymptomatic vs symptomatic case

epiPars['gamma_e'] = 1/3           # Latent period (He et al)
epiPars['gamma_a']= 1/7            # Recovery rate, undocumented (Kissler et al)
epiPars['gamma_s']= 1/7            # Recovery rate, undocumented (Kissler et al)
epiPars['gamma_hs']= 1/5           # LOS for subcritical cases (medrxiv paper)
epiPars['gamma_hc']= 1/7           # LOS for critical cases (medrxiv paper)
epiPars['p'] = 0.3#0.14            # Fraction 'Symptomatic'documented' (Shaman's paper)

epiPars['hosp_frac'] = np.array((0.061, 0.182, 0.417))  #From MMWR
epiPars['hosp_crit'] = np.array((0, 0.063, 0.173))      #From CDC, MMWR
epiPars['crit_die'] = np.array((0, 0.5, 0.5))           #Obtained from initial fitting

# Expand epiPars
epiPars['hosp_frac_5'] = epiPars['hosp_frac'][np.array((0,1,1,1,2))]
epiPars['hosp_crit_5'] = epiPars['hosp_crit'][np.array((0,1,1,1,2))]
epiPars['crit_die_5'] = epiPars['crit_die'][np.array((0,1,1,1,2))]


# Contact Interaction Pars ---------------------------------------------------------------------------
contactPars = {}

# Contact Matrices: 3x3 data from Prem et al
contactPars['AllContacts'] = np.array((9.75, 2.57, 0.82, 5.97, 10.32, 2.25, 0.39, 0.46, 1.20)).reshape((3,3)).T
contactPars['WorkContacts'] = np.array((0.20, 0.28, 0, 0.64, 4.73, 0, 0, 0, 0)).reshape((3,3)).T
contactPars['SchoolContacts'] = np.array((4.32, 0.47, 0.02, 1.10, 0.32, 0.04, 0.01, 0.01, 0.03)).reshape((3,3)).T
contactPars['HomeContacts'] = np.array((2.03, 1.02, 0.50, 2.37, 1.82, 0.68, 0.24, 0.14, 0.62)).reshape((3,3)).T
contactPars['OtherContacts'] = np.array((3.20, 0.80, 0.30, 1.86, 3.45, 1.53, 0.14, 0.32, 0.55)).reshape((3,3)).T

contactPars['frac_home'] = 0.316
contactPars['frac_reduced'] = 0.6275
contactPars['frac_full'] = 0.0565

# Include Expansions
contactPars['WorkContacts_5x5'] = Expand_Contact_Matrices(contactPars['WorkContacts']
                                                          , contactPars['frac_home']
                                                          , contactPars['frac_reduced']
                                                          , contactPars['frac_full'])
contactPars['SchoolContacts_5x5'] = Expand_Contact_Matrices(contactPars['SchoolContacts']
                                                            , contactPars['frac_home']
                                                            , contactPars['frac_reduced']
                                                            , contactPars['frac_full'])
contactPars['HomeContacts_5x5'] = Expand_Contact_Matrices(contactPars['HomeContacts']
                                                          , contactPars['frac_home']
                                                          , contactPars['frac_reduced']
                                                          , contactPars['frac_full'])
contactPars['OtherContacts_5x5'] = Expand_Contact_Matrices(contactPars['OtherContacts']
                                                           , contactPars['frac_home']
                                                           , contactPars['frac_reduced']
                                                           , contactPars['frac_full'])


# Intervention Pars ---------------------------------------------------------------------------
intvPars = {}

# Test Data for Cellex
intvPars['sensitivity'] = 0.94
intvPars['specificity'] = 0.96
intvPars['daily_tests'] = 10**3

# Other intervention parameters
intvPars['tStart_distancing'] = 70
intvPars['tStart_target'] = 115
intvPars['tStart_school'] = 230

intvPars['socialDistancing_other'] = .25
intvPars['p_reduced'] = 0.1 # proportion of contacts reduced 
intvPars['p_full'] = 1 # proportion of contacts reduced for full contact adults

intvPars['alpha'] = 4 # shielding. Note this is not alpha_JSW, but (alpha_JSW + 1)
intvPars['c'] = 1

intvPars['socialDistancing_other_c'] = 1-(0.75*intvPars['c'])
intvPars['p_reduced_c'] = 1-(0.9*intvPars['c'])

intvPars['tStart_test'] = 107 # can change when in the outbreak testing becomes available
intvPars['tStart_target'] = 115


# Initial conditions ---------------------------------------------------------------------------
inits = {}

inits['I_c0'] = 60
inits['I_a0'] = 20
inits['I_rc0'] = 50
inits['I_fc0'] = 1
inits['I_e0'] = 40

X0 = np.zeros(pars['nTotSubComp'])

# Remove from susceptible pool ...
X0[0] = pars['agestruc'][0] * pars['N'] - inits['I_c0'] 
X0[1] = pars['agestruc'][1] * contactPars['frac_home'] * pars['N'] - inits['I_a0'] 
X0[2] = pars['agestruc'][1] * contactPars['frac_reduced'] * pars['N'] - inits['I_rc0'] 
X0[3] = pars['agestruc'][1] * contactPars['frac_full'] * pars['N'] - inits['I_fc0'] 
X0[4] = pars['agestruc'][2] * pars['N'] - inits['I_e0'] 

# ... and add to infected
X0[10] = inits['I_c0']
X0[11] = inits['I_a0']
X0[12] = inits['I_rc0']
X0[13] = inits['I_fc0']
X0[14] = inits['I_e0']