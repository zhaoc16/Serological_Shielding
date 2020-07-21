import numpy as np
from utils.helpers import Expand_Contact_Matrices, Expand_10x10

def seir_model_shields_rcfc_nolatent_time(X,t,parsList):
    '''
    - Subscripts are defined as follows: c = children, a = non-essential adults, rc = reduced contact adults, fc = full contact adults, e = elderly
    - For infection equations: The I compartments correspond to infectious periods, not periods in which people are symptomatic
    - I's are cases that are severe enough to be documented eventually and Ia are undocumented cases
    - Since we are not fitting to data we do not need to model specifically when the cases are documented 
    - pos indicates people who have tested positive for COVID-19 by antibody test
    '''
    
    # Load All Parameter Sets
    pars = parsList[0]
    epiPars = parsList[1]
    contactPars = parsList[2]
    intvPars = parsList[3]
    
    def time_switch1(t):
        return(0+(t<=intvPars['tStart_test']))

    def time_switch2(t):
        return(0+(t>intvPars['tStart_test']))
    
    # Update Pars
    intvPars['socialDistancing_other_c'] = 1-(0.75*intvPars['c'])
    intvPars['p_reduced_c'] = 1-(0.9*intvPars['c'])
    
    # Load all variables
    S_c, S_a, S_rc, S_fc, S_e = X[pars['S_ids']]
    E_c, E_a, E_rc, E_fc, E_e = X[pars['E_ids']]
    Isym_c, Isym_a, Isym_rc, Isym_fc, Isym_e = X[pars['Isym_ids']]
    Iasym_c, Iasym_a, Iasym_rc, Iasym_fc, Iasym_e = X[pars['Iasym_ids']]
    Hsub_c, Hsub_a, Hsub_rc, Hsub_fc, Hsub_e = X[pars['Hsub_ids']]
    Hcri_c, Hcri_a, Hcri_rc, Hcri_fc, Hcri_e = X[pars['Hcri_ids']]
    D_c, D_a, D_rc, D_fc, D_e = X[pars['D_ids']]
    R_c, R_a, R_rc, R_fc, R_e = X[pars['R_ids']]

    S_c_pos, S_a_pos, S_rc_pos, S_fc_pos, S_e_pos = X[pars['S_pos_ids']]
    E_c_pos, E_a_pos, E_rc_pos, E_fc_pos, E_e_pos = X[pars['E_pos_ids']]
    Isym_c_pos, Isym_a_pos, Isym_rc_pos, Isym_fc_pos, Isym_e_pos = X[pars['Isym_pos_ids']]
    Iasym_c_pos, Iasym_a_pos, Iasym_rc_pos, Iasym_fc_pos, Iasym_e_pos = X[pars['Iasym_pos_ids']]
    R_c_pos, R_a_pos, R_rc_pos, R_fc_pos, R_e_pos = X[pars['R_pos_ids']]

    # Infectiousness of each subgroup
    ifc_c_gen, ifc_a_gen, ifc_rc_gen, ifc_fc_gen, ifc_e_gen = epiPars['asymp_red'] * X[pars['Iasym_ids']] + X[pars['Isym_ids']]
    ifc_c_pos, ifc_a_pos, ifc_rc_pos, ifc_fc_pos, ifc_e_pos = epiPars['asymp_red'] * X[pars['Iasym_pos_ids']] + X[pars['Isym_pos_ids']]
    ifc_comb = np.array([ifc_c_pos, ifc_c_gen, ifc_a_pos, ifc_a_gen, ifc_rc_pos, ifc_rc_gen, ifc_fc_pos, ifc_fc_gen, ifc_e_pos, ifc_e_gen])

    # Totals
    tot_c, tot_a, tot_rc, tot_fc, tot_e = [np.sum(X[pars[sub + '_ids']]) for sub in ('c', 'a', 'rc', 'fc', 'e')] # These should be constant
    tot_c_gen, tot_a_gen, tot_rc_gen, tot_fc_gen, tot_e_gen = [np.sum(X[pars[sub + '_ids'][np.arange(8)]]) for sub in ('c', 'a', 'rc', 'fc', 'e')]
    tot_c_pos, tot_a_pos, tot_rc_pos, tot_fc_pos, tot_e_pos = [np.sum(X[pars[sub + '_ids'][np.arange(start=8,stop=13)]]) for sub in ('c', 'a', 'rc', 'fc', 'e')]
    tot_comb = np.array([tot_c_pos, tot_c_gen, tot_a_pos, tot_a_gen, tot_rc_pos, tot_rc_gen, tot_fc_pos, tot_fc_gen, tot_e_pos, tot_e_gen])
    
    # Released Fraction
    frac_released = np.array((tot_c_pos/tot_c, tot_a_pos/tot_a, tot_rc_pos/tot_rc, tot_fc_pos/tot_fc, tot_e_pos/tot_e))
    frac_distanced = 1-frac_released
    frac_comb = np.stack((frac_released, frac_distanced), axis=-1).flatten() # interleaved
    
    # Derive contact matrices:
    # - Matrix entry i,j is the number of contacts/day that an individual in group i has with group j
    # - These are normalized for population size and corrected for symmetry using standard methods
    
    # These will be useful later
    temp_even_10x10_Indices = np.meshgrid(2*np.arange(5), 2*np.arange(5))
    temp_rep_0011223344 = np.repeat(np.arange(5), 2)
    
    # Work Contacts
    CM_Work_High = contactPars['WorkContacts_5x5'] * frac_released
    CM_Work_Low = contactPars['WorkContacts_5x5'] * frac_distanced
    CM_Work_Comb = Expand_10x10(CM_Work_High, CM_Work_Low)
    
    # Broad Distancing Matrix
    # - Children have no work contacts
    # - Non-essential workers don't work
    # - Reduced contact workers reduce contacts by preduced and those contacts are proportional to 
    #   prevalence in the population
    # - Full contact workers keep all of their baseline contacts
    # - The distribution of those contacts is based on prevalence in the population
    CM_Work_Comb_Distancing = np.zeros((10,10))
    CM_Work_Comb_Distancing[4:6,:] = contactPars['WorkContacts_5x5'][2,temp_rep_0011223344]*frac_comb*intvPars['p_reduced']
    CM_Work_Comb_Distancing[6:8,:] = contactPars['WorkContacts_5x5'][3,temp_rep_0011223344]*frac_comb*intvPars['p_full']

    # Targeted distancing work
    # - DISTANCING CONDITIONAL ON TESTING: People working from home regain their workplace contacts if they test positive 
    # - People in reduced and full contact occupations have the distribution of their workplace contacts changed by alpha    
    # - We use fixed shielding (except we multiply by alpha and not alpha+1), 
    #   so we correct for the situation where alpha times prevalence is greater than 1
    alpha_Scale = frac_released*intvPars['alpha']
    alpha_Scale = np.array([1.0 if scale>1 else scale for scale in alpha_Scale])
    
    temp_rc = CM_Work_Comb[4,:]*intvPars['p_reduced_c']
    
    CM_Work_Comb_TargetedDistancing = np.zeros((10,10))
    CM_Work_Comb_TargetedDistancing[2,:] = np.stack((contactPars['WorkContacts_5x5'][1,:], np.zeros(5)), axis=-1).flatten()
    CM_Work_Comb_TargetedDistancing[4:6,0::2] = np.array([temp_rc[2*i] + temp_rc[2*i+1]*alpha_Scale[i] for i in np.arange(5)])
    CM_Work_Comb_TargetedDistancing[4:6,1::2] = np.array([temp_rc[2*i] + temp_rc[2*i+1]*(1-alpha_Scale[i]) for i in np.arange(5)])
    CM_Work_Comb_TargetedDistancing[6:8,:] = np.copy(CM_Work_Comb_TargetedDistancing[4:6,:])/intvPars['p_reduced_c']
    
    # School Contacts
    CM_School_High = contactPars['SchoolContacts_5x5'] * frac_released
    CM_School_Low = contactPars['SchoolContacts_5x5'] * frac_distanced
    CM_School_Comb = Expand_10x10(CM_School_High, CM_School_Low)

    CM_School_Comb_reopen = np.zeros((10,10))
    CM_School_Comb_reopen[temp_even_10x10_Indices] = contactPars['SchoolContacts_5x5'].T # meshgrid transposed
        
    # Home Contacts
    CM_Home_High = contactPars['HomeContacts_5x5'] * frac_released
    CM_Home_Low = contactPars['HomeContacts_5x5'] * frac_distanced
    CM_Home_Comb = Expand_10x10(CM_Home_High, CM_Home_Low)
    
    # Other Contacts
    CM_Other_High = contactPars['OtherContacts_5x5'] * frac_released
    CM_Other_Low = contactPars['OtherContacts_5x5'] * frac_distanced
    CM_Other_Comb = Expand_10x10(CM_Other_High, CM_Other_Low)
    
    CM_Other_Comb_Distancing = np.copy(CM_Other_Comb)*intvPars['socialDistancing_other']
    
    CM_Other_Comb_TargetedDistancing = Expand_10x10(contactPars['OtherContacts_5x5']*intvPars['socialDistancing_other_c']*alpha_Scale
                                                    ,contactPars['OtherContacts_5x5']*intvPars['socialDistancing_other_c']*(1-alpha_Scale))
    CM_Other_Comb_TargetedDistancing2 = np.copy(CM_Other_Comb_TargetedDistancing)
    CM_Other_Comb_TargetedDistancing2[2*np.arange(5),:] = CM_Other_Comb_TargetedDistancing2[2*np.arange(5),:] / intvPars['socialDistancing_other_c']
    
    # Matrices Change with respect to Control Strategy
    if (t<intvPars['tStart_distancing']):
        CM = CM_Work_Comb + CM_School_Comb + CM_Home_Comb + CM_Other_Comb
    elif (t < intvPars['tStart_target']):
        CM = CM_Work_Comb_Distancing + CM_Home_Comb + CM_Other_Comb_Distancing
    elif (t < intvPars['tStart_school']):
        CM = CM_Other_Comb_TargetedDistancing2 + CM_Home_Comb + CM_Other_Comb_TargetedDistancing
    else:
        CM = CM_Other_Comb_TargetedDistancing2 + CM_Home_Comb + CM_Other_Comb_TargetedDistancing + CM_School_Comb_reopen #XXX ori code has baseline

    # Force of Infections
    tot_comb[0::2] = np.array([1.0 if ifc_in==0 else ifc_in for ifc_in in tot_comb[0::2] ]) # avoid division by 0
    foi_comb = np.array([epiPars['q']*np.sum(CM[i,:]*ifc_comb/tot_comb) for i in np.arange(10)])
    foi_c_pos, foi_c_gen, foi_a_pos, foi_a_gen, foi_rc_pos, foi_rc_gen, foi_fc_pos, foi_fc_gen, foi_e_pos, foi_e_gen = foi_comb
    
    # Testing
    child_tests = intvPars['daily_tests'] * pars['agestruc'][0] # Get nTests available by group
    adult_tests_h = intvPars['daily_tests'] * pars['agestruc'][1] * contactPars['frac_home']
    adult_tests_rc = intvPars['daily_tests'] * pars['agestruc'][1] * contactPars['frac_reduced']
    adult_tests_fc = intvPars['daily_tests'] * pars['agestruc'][1] * contactPars['frac_full']
    el_tests = intvPars['daily_tests'] * pars['agestruc'][2]
            
    # Divide by people eligible to be tested to get proportion tested per day
    test_c = child_tests / np.sum((S_c, E_c, Iasym_c, Hsub_c, Hcri_c, R_c))
    test_a = adult_tests_h / np.sum((S_a, E_a, Iasym_a, Hsub_a, Hcri_a, R_a))
    test_rc = adult_tests_rc / np.sum((S_rc, E_rc, Iasym_rc, Hsub_rc, Hcri_rc, R_rc))
    test_fc = adult_tests_fc / np.sum((S_fc, E_fc, Iasym_fc, Hsub_fc, Hcri_fc, R_fc))
    test_e = el_tests / np.sum((S_e, E_e, Iasym_e, Hsub_e, Hcri_e, R_e))
    
    test_comb = np.array((test_c, test_a, test_rc, test_fc, test_e))
    
    test_switch1 = time_switch1(t)
    test_switch2 = time_switch2(t)

    # Model Equations
    dS_c, dS_a, dS_rc, dS_fc, dS_e = - foi_comb[1::2]*X[pars['S_ids']] \
                                     - (1-intvPars['specificity'])*test_comb*test_switch2*X[pars['S_ids']]
    dE_c, dE_a, dE_rc, dE_fc, dE_e = foi_comb[1::2]*X[pars['S_ids']] \
                                     - epiPars['gamma_e']*X[pars['E_ids']] \
                                     - (1-intvPars['specificity'])*test_comb*test_switch2*X[pars['E_ids']]
    dIsym_c, dIsym_a, dIsym_rc, dIsym_fc, dIsym_e = epiPars['gamma_e']*X[pars['E_ids']]*epiPars['p'] \
                                                    - epiPars['gamma_s']*X[pars['Isym_ids']]
    dIasym_c, dIasym_a, dIasym_rc, dIasym_fc, dIasym_e = epiPars['gamma_e']*X[pars['E_ids']]*(1-epiPars['p']) \
                                                         - epiPars['gamma_a']*X[pars['Iasym_ids']] \
                                                         - (1-intvPars['specificity'])*test_comb*test_switch2*X[pars['Iasym_ids']]
    dHsub_c, dHsub_a, dHsub_rc, dHsub_fc, dHsub_e = epiPars['gamma_s']*X[pars['Isym_ids']]*(epiPars['hosp_frac_5']-epiPars['hosp_crit_5']) \
                                                    + epiPars['gamma_s']*X[pars['Isym_pos_ids']]*(epiPars['hosp_frac_5']-epiPars['hosp_crit_5']) \
                                                    - epiPars['gamma_hs']*X[pars['Hsub_ids']]
    dHcri_c, dHcri_a, dHcri_rc, dHcri_fc, dHcri_e = epiPars['gamma_s']*X[pars['Isym_ids']]*epiPars['hosp_crit_5'] \
                                                    + epiPars['gamma_s']*X[pars['Isym_pos_ids']]*epiPars['hosp_crit_5'] \
                                                    - epiPars['gamma_hc']*X[pars['Hcri_ids']]
    dD_c, dD_a, dD_rc, dD_fc, dD_e = epiPars['gamma_hc']*X[pars['Hcri_ids']]*epiPars['crit_die_5']
    dR_c, dR_a, dR_rc, dR_fc, dR_e = (1-epiPars['hosp_frac_5'])*epiPars['gamma_s']*X[pars['Isym_ids']] \
                                     + epiPars['gamma_a']*X[pars['Iasym_ids']] \
                                     + (1-intvPars['sensitivity'])*test_switch2*epiPars['gamma_hs']*X[pars['Hsub_ids']] \
                                     + (1-intvPars['sensitivity'])*test_switch2*epiPars['gamma_hc']*X[pars['Hcri_ids']]*(1-epiPars['crit_die_5']) \
                                     - intvPars['sensitivity']*test_comb*test_switch2*X[pars['R_ids']] \
                                     + epiPars['gamma_hc']*test_switch1*X[pars['Hcri_ids']]*(1-epiPars['crit_die_5']) \
                                     + test_switch1*epiPars['gamma_hs']*X[pars['Hsub_ids']]

    dS_c_pos, dS_a_pos, dS_rc_pos, dS_fc_pos, dS_e_pos = (1-intvPars['specificity'])*test_comb*test_switch2*X[pars['S_ids']] \
                                                         - foi_comb[0::2]*X[pars['S_pos_ids']]
    dE_c_pos, dE_a_pos, dE_rc_pos, dE_fc_pos, dE_e_pos = (1-intvPars['specificity'])*test_comb*test_switch2*X[pars['E_ids']] \
                                                         + foi_comb[0::2]*X[pars['S_pos_ids']] \
                                                         - epiPars['gamma_e']*X[pars['E_pos_ids']]
    dIsym_c_pos, dIsym_a_pos, dIsym_rc_pos, dIsym_fc_pos, dIsym_e_pos = epiPars['gamma_e']*X[pars['E_pos_ids']]*epiPars['p'] \
                                                                        - epiPars['gamma_s']*X[pars['Isym_pos_ids']]
    dIasym_c_pos, dIasym_a_pos, dIasym_rc_pos, dIasym_fc_pos, dIasym_e_pos = (1-intvPars['specificity'])*test_comb*test_switch2*X[pars['Iasym_ids']] \
                                                                             + epiPars['gamma_e']*X[pars['E_pos_ids']]*(1-epiPars['p']) \
                                                                             - epiPars['gamma_a']*X[pars['Iasym_pos_ids']]
    dR_c_pos, dR_a_pos, dR_rc_pos, dR_fc_pos, dR_e_pos = intvPars['sensitivity']*test_comb*test_switch2*X[pars['R_ids']] \
                                                         + intvPars['sensitivity']*epiPars['gamma_hc']*test_switch2*X[pars['Hcri_ids']]*(1-epiPars['crit_die_5']) \
                                                         + intvPars['sensitivity']*epiPars['gamma_hs']*test_switch2*X[pars['Hsub_ids']] \
                                                         + (1-epiPars['hosp_frac_5'])*(epiPars['gamma_s']*X[pars['Isym_pos_ids']]) \
                                                         + epiPars['gamma_a']*X[pars['Iasym_pos_ids']]
    
    dXdt = np.array((dS_c, dS_a, dS_rc, dS_fc, dS_e
                   , dE_c, dE_a, dE_rc, dE_fc, dE_e
                   , dIsym_c, dIsym_a, dIsym_rc, dIsym_fc, dIsym_e
                   , dIasym_c, dIasym_a, dIasym_rc, dIasym_fc, dIasym_e
                   , dHsub_c, dHsub_a, dHsub_rc, dHsub_fc, dHsub_e
                   , dHcri_c, dHcri_a, dHcri_rc, dHcri_fc, dHcri_e
                   , dD_c, dD_a, dD_rc, dD_fc, dD_e
                   , dR_c, dR_a, dR_rc, dR_fc, dR_e
                   , dS_c_pos, dS_a_pos, dS_rc_pos, dS_fc_pos, dS_e_pos
                   , dE_c_pos, dE_a_pos, dE_rc_pos, dE_fc_pos, dE_e_pos
                   , dIsym_c_pos, dIsym_a_pos, dIsym_rc_pos, dIsym_fc_pos, dIsym_e_pos
                   , dIasym_c_pos, dIasym_a_pos, dIasym_rc_pos, dIasym_fc_pos, dIasym_e_pos
                   , dR_c_pos, dR_a_pos, dR_rc_pos, dR_fc_pos, dR_e_pos))
  
    return dXdt