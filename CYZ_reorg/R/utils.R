Interleave = function(MatA, MatB, byCol = TRUE){
  # This function interleaves two matrices byCol or byRow
  ncolA = ncol(MatA)
  ncolB = ncol(MatB)
  nrowA = nrow(MatA)
  nrowB = nrow(MatB)
  
  if(byCol){
    ret = matrix(0, ncol = ncolA+ncolB, nrow = nrowA)
    ret[,(1:ncolA)*2-1] = MatA
    ret[,(1:ncolB)*2] = MatB
  }else{
    ret = matrix(0, ncol = ncolA, nrow = nrowA+nrowB)
    ret[(1:nrowA)*2-1,] = MatA
    ret[(1:nrowB)*2,] = MatB
  }
  
  return(ret)
}

Expand_5x5 = function(contactMatrix, p_Home, p_Reduced, p_Full){
  # This function expands the 3x3 contact matrices into 5x5 ones by separating
  # adult contacts based on home, reduced, or full contact probabilities.
  temp_rowexpand = contactMatrix[c(1,2,2,2,3),]
  temp_expand = temp_rowexpand[,c(1,2,2,2,3)]
  temp_expand[,c(2,3,4)] = t(t(temp_expand[,c(2,3,4)])* c(p_Home, p_Reduced, p_Full)) # R handles element-wise oddly.
  
  return(temp_expand)
}

Expand_10x10 = function(HighContact, LowContact, ZInflate=FALSE){
  # This function takes a high and low contact matrix and combines them
  # This function can also expand a single matrix by appending zeros.
  if(ZInflate){
    LowContact = matrix(0, ncol = ncol(HighContact), nrow = nrow(HighContact))
  }
  
  temp = Interleave(HighContact, LowContact, byCol=TRUE)
  ret = temp[rep(1:5, each=2),]
  
  if(ZInflate){
    ret[(1:5)*2,] = 0
  }
  
  return(ret)
}

Get_Var_Names = function(Subgroups=model_pars$subgroups, Compartments=model_pars$compartments){
  # This function is to keep the input file cleaner.
  v.out = apply(expand.grid(Subgroups, Compartments), 1, function(x){
    v.split_in = strsplit(x[2], '_')[[1]]
    if(length(v.split_in)-1){
      ret = paste(v.split_in[1], '_', x[1], '_', v.split_in[2], sep='', collapse='')
    }else{
      ret = paste(v.split_in[1], '_', x[1], sep='', collapse='')
    }
    return(ret)
  })
  
  return(v.out)
}

Vec_to_Mat = function(v.Vars, Subgroups=model_pars$subgroups, Compartments=model_pars$compartments){
  # This function takes the variable vector and converts to a matrix
  ret = matrix(v.Vars, ncol=length(Compartments))
  colnames(ret) = Compartments
  rownames(ret) = Subgroups
  
  return(ret)
}

Mat_to_Vec = function(mat.Vars, Pars = model_pars){
  # This function takes the variable matrix and flattens it
  ret = as.vector(mat.Vars)
  names(ret) = Pars$varNames
  
  return(ret)
}

seir_model_shields_rcfc_nolatent = function(t, X, Pars) {
  # Core transmission model
  #   - Subscripts are defined as follows: c=children, a=non-essential adults, 
  #     rc=reduced contact adults, fc=full contact adults, e=elderly
  #   - For infection equations: The I compartments correspond to infectious 
  #     periods, not periods in which people are symptomatic
  #   - Is are cases that are severe enough to be documented eventually and Ia are undocumented cases
  #   - Since we are not fitting to data we do not need to model specifically when the cases are documented 
  #   - _pos indicates people who have tested positive for COVID-19 by antibody test
  
  
  test.switch1<-sw1fxn(t)
  test.switch2<-sw2fxn(t)
  
  
  # (0) Setup ---------------------------------------------------------------

  # Load all pars
  for(i_par in 1:length(Pars)){
    assign(names(Pars)[i_par], Pars[[i_par]])
  }
  
  # Update Pars
  socialDistancing_other = 1-(0.75*c)
  p_reduced_c = 1-(0.9*c)
  
  # Load all Vars
  for(i_var in 1:length(X)){
    assign(names(X)[i_var], X[[i_var]])
  }
  
  # Load in Matrix Form
  mat_X = Vec_to_Mat(X, Subgroups=subgroups, Compartments=compartments)
  mat_X_pos = mat_X[,grepl('_pos', colnames(mat_X))]
  mat_X_gen = mat_X[,!grepl('_pos', colnames(mat_X))]
  
  # Number of infections at time t for each subgroup
  infec_gen = asymp_red*mat_X[,'Iasym'] + mat_X[,'Isym']
  infec_pos = asymp_red*mat_X[,'Iasym_pos'] + mat_X[,'Isym_pos']
  
  # Population sizes/testing status at time t    
  tot = rowSums(mat_X[,colnames(mat_X)!='D'])
  tot_pos = rowSums(mat_X_pos)
  tot_gen = rowSums(mat_X_gen[,colnames(mat_X_gen)!='D'])

  #Fraction of population in each group who has tested positive by time
  frac_released = tot_pos/tot
  frac_distanced = c(rep(1,5))-frac_released
  frac_comb = as.vector(rbind(frac_released, frac_distanced)) # interleaved

  tot_pos[tot_pos==0] = tot[tot_pos==0] # avoid division by 0

  
  # (1) Derive Contact Matrices ---------------------------------------------
  
  # Home contacts
  HighContact_Home = HomeContacts_5x5 %*% diag(frac_released)
  LowContact_Home = HomeContacts_5x5 %*% diag(frac_distanced)
  HomeContacts_10x10 = Expand_10x10(HighContact_Home, LowContact_Home)
  
  # School contacts 
  HighContact_School = SchoolContacts_5x5 %*% diag(frac_released)
  LowContact_School = SchoolContacts_5x5 %*% diag(frac_distanced)
  SchoolContacts_10x10_Baseline = Expand_10x10(HighContact_School, LowContact_School)
  SchoolContacts_10x10_reopen=Expand_10x10(SchoolContacts_5x5, 0, ZInflate=T)  #With general distancing, schools closed
  
  # Work contacts
  HighContact_Work = WorkContacts_5x5 %*% diag(frac_released)
  LowContact_Work = WorkContacts_5x5 %*% diag(frac_distanced)
  WorkContacts_10x10_Baseline = Expand_10x10(HighContact_Work, LowContact_Work) #Baseline matrix
  
  # Broad distancing matrix
  #   - Children have no work contacts
  #   - Non-essential workers don't work
  #   - Reduced contact workers reduce contacts by p_reduced and those contacts are 
  #     proportional to prevalence in the population
  WorkContacts_10x10_Distancing = matrix(0, ncol=10, nrow=10)
  WorkContacts_10x10_Distancing[5:6,] = WorkContacts_10x10_Baseline[5:6,]*p_reduced
  WorkContacts_10x10_Distancing[7:8,] = WorkContacts_10x10_Baseline[7:8,]*p_full
  
  # Targeted distancing work
  #   - DISTANCING CONDITIONAL ON TESTING: People working from home regain their 
  #     workplace contacts if they test positive 
  #   - People in reduced and full contact occupations have the distribution of 
  #     their workplace contacts changed by alpha
  #   - We use fixed shielding (except we multiply by alpha and not alpha +1), 
  #     so we correct for the situation where alpha times prevalence is greater than 1
  Scale=frac_released*alpha
  Scale[Scale>1]=1  
  
  WorkContacts_10x10_TargetedDistancing = matrix(0, nrow=10, ncol=10)
  DistancedWork_RedContact = colSums(matrix(WorkContacts_10x10_Baseline[5,], nrow=2))
  WorkContacts_10x10_TargetedDistancing[3,] = as.vector(rbind(WorkContacts_5x5[2,], rep(0,5)))
  WorkContacts_10x10_TargetedDistancing[5,] = as.vector(rbind(Scale*DistancedWork_RedContact
                                                              , (1-Scale)*DistancedWork_RedContact))*p_reduced_c
  
  WorkContacts_10x10_TargetedDistancing[6,]=WorkContacts_10x10_TargetedDistancing[5,]
  WorkContacts_10x10_TargetedDistancing[7,]=WorkContacts_10x10_TargetedDistancing[5,]/p_reduced_c
  WorkContacts_10x10_TargetedDistancing[8,]=WorkContacts_10x10_TargetedDistancing[5,]/p_reduced_c
  
  # Other contacts
  #   - BASELINE OTHER CONTACTS MATRIX: People mix freely, contacts are based 
  #     on proportion in the population
  HighContact_Other = OtherContacts_5x5 %*% diag(frac_released)
  LowContact_Other = OtherContacts_5x5 %*% diag(frac_distanced)
  OtherContacts_10x10_Baseline = Expand_10x10(HighContact_Other, LowContact_Other)
  
  #   - UNTARGETED SOCIAL DISTANCING OTHER CONTACTS: All contacts for social 
  #     distancing are reduced by a fraction
  OtherContacts_10x10_Distancing = OtherContacts_10x10_Baseline*socialDistancing_other
  
  #   - TARGETED SOCIAL DISTANCING OTHER CONTACTS: Shielding
  #   - Contact probabilities proportional to prevalence in the population times
  #     a scaling factor to account for shielding
  #   - Keep the total contacts the same by subtracting from the total number 
  #     of expected contacts between age classes
  #   - Once the scaling factor*prevalence exceeds 1, set to 1 (so that all 
  #     contacts are with a positive individual)
  temp.other = Expand_10x10(OtherContacts_5x5, OtherContacts_5x5)
  temp.other = t(t(temp.other)*as.vector(rbind(Scale, (1-Scale))))*socialDistancing_other_c
  OtherContacts_10x10_TargetedDistancing = temp.other
  
  #   - For test positive individuals, keep the distribution of other contacts 
  #     shielded but change the number to pre-pandemic levels
  OtherContacts_10x10_TargetedDistancing_c=OtherContacts_10x10_TargetedDistancing
  OtherContacts_10x10_TargetedDistancing_c[(1:5)*2-1,]=OtherContacts_10x10_TargetedDistancing[(1:5)*2-1,]/socialDistancing_other_c

  
  # (2) Evolve Contact Matrices over time -----------------------------------
  
  #Change matrices used over time
  if(t<tStart_distancing | t>=tStart_reopen){ #Use baseline matrices until social distancing starts
    CM = HomeContacts_10x10 + 
      SchoolContacts_10x10_Baseline + 
      WorkContacts_10x10_Baseline +
      OtherContacts_10x10_Baseline 
  }
  
  if(t>=tStart_distancing & t<tStart_reopen){
    # if(t>=tStart_distancing & t<tStart_target){ #Use these matrices under general social distancing without testing
    CM = HomeContacts_10x10 +
      WorkContacts_10x10_Distancing +
      OtherContacts_10x10_Distancing
  }
  
  if(t>=tStart_target & t<tStart_school){ #Use these matrices between the start of targeting and schools reopening
    CM = HomeContacts_10x10 +
      WorkContacts_10x10_TargetedDistancing + 
      OtherContacts_10x10_TargetedDistancing_c
  }
  
  if(t>=tStart_school){ #Keep most matrices the same but add schools back in
    CM = HomeContacts_10x10 + 
      WorkContacts_10x10_TargetedDistancing + 
      OtherContacts_10x10_TargetedDistancing_c + 
      SchoolContacts_10x10_Baseline
  }
  
  
  # (3) Calculate v.fois ------------------------------------------------------
  
  # Force of infection by group 
  temp.I = as.vector(rbind(infec_pos/tot_pos, infec_gen/tot_gen))
  temp.I = t(t(CM)*temp.I)*q
  v.fois = rowSums(temp.I)
  names(v.fois) = c('foi_c_pos', 'foi_c_gen', 'foi_a_pos', 'foi_a_gen'
                  , 'foi_rc_pos', 'foi_rc_gen', 'foi_fc_pos', 'foi_fc_gen'
                  , 'foi_e_pos', 'foi_e_gen')
  
  # Fois by pos or gen
  v.fois_gen = v.fois[(1:5)*2]
  v.fois_pos = v.fois[(1:5)*2-1]
  
  
  # (4) Sero Testing --------------------------------------------------------

  # Get number of tests available by group
  v.avail_tests = daily_tests*agestruc[c(1,2,2,2,3)]*c(1,frac_home,frac_reduced,frac_full,1)
  v.eligible_pop = rowSums(mat_X_gen[,c('S', 'E', 'Iasym', 'Hsub', 'Hcri', 'R')])
  
  # Divide by number of test-eligible people for daily test fraction
  v.daily_tests = v.avail_tests/v.eligible_pop
  
  # When prevalence and testing are high, might not be enough eligible:
  v.daily_tests[v.daily_tests>1] = 1 
  
  
  # (5) Model Equations -----------------------------------------------------
  v.hosp_frac = hosp_frac[c(1,2,2,2,3)]
  v.hosp_crit = hosp_crit[c(1,2,2,2,3)]
  v.crit_die = crit_die[c(1,2,2,2,3)]
  
  dS = -v.fois_gen*mat_X[,'S'] - 
    (1-specificity)*test.switch2*v.daily_tests*mat_X[,'S']
  
  dE = v.fois_gen*mat_X[,'S'] - 
    gamma_e*mat_X[,'E'] - 
    (1-specificity)*test.switch2*v.daily_tests*mat_X[,'E']
  
  dIsym = p*gamma_e*mat_X[,'E'] - 
    gamma_s*mat_X[,'Isym']
  
  dIasym = (1-p)*gamma_e*mat_X[,'E'] - 
    gamma_a*mat_X[,'Iasym'] - 
    (1-specificity)*test.switch2*v.daily_tests*mat_X[,'Iasym']
  
  dHsub = gamma_s*(v.hosp_frac-v.hosp_crit)*mat_X[,'Isym'] + 
    gamma_s*(v.hosp_frac-v.hosp_crit)*mat_X[,'Isym_pos'] - 
    gamma_hs*mat_X[,'Hsub']
  
  dHcri = gamma_s*v.hosp_crit*mat_X[,'Isym'] + 
    gamma_s*v.hosp_crit*mat_X[,'Isym_pos'] - 
    gamma_hc*mat_X[,'Hcri']
  
  dD = gamma_hc*v.crit_die*mat_X[,'Hcri']
  
  dR = gamma_s*(1-v.hosp_frac)*mat_X[,'Isym'] + 
    gamma_a*mat_X[,'Iasym'] + 
    (1-sensitivity)*test.switch2*gamma_hs*mat_X[,'Hsub'] + 
    (1-sensitivity)*test.switch2*gamma_hc*(1-v.crit_die)*mat_X[,'Hcri'] - 
    sensitivity*test.switch2*v.daily_tests*mat_X[,'R'] + 
    gamma_hc*test.switch1*(1-v.crit_die)*mat_X[,'Hcri'] + 
    gamma_hs*test.switch1*mat_X[,'Hsub']
  
  dS_pos = (1-specificity)*test.switch2*v.daily_tests*mat_X[,'S'] - 
    v.fois_pos*mat_X[,'S_pos']
  
  dE_pos = (1-specificity)*test.switch2*v.daily_tests*mat_X[,'E'] + 
    v.fois_pos*mat_X[,'S_pos'] - 
    gamma_e*mat_X[,'E_pos']
  
  dIsym_pos = p*gamma_e*mat_X[,'E_pos'] -
    gamma_s*mat_X[,'Isym_pos']
  
  dIasym_pos = (1-specificity)*test.switch2*v.daily_tests*mat_X[,'Iasym'] + 
    (1-p)*gamma_e*mat_X[,'E_pos'] - 
    gamma_a*mat_X[,'Iasym_pos']
  
  dR_pos = sensitivity*test.switch2*v.daily_tests*mat_X[,'R'] + 
    sensitivity*test.switch2*gamma_hc*(1-v.crit_die)*mat_X[,'Hcri'] + 
    sensitivity*test.switch2*gamma_hs*mat_X[,'Hsub'] + 
    gamma_s*(1-v.hosp_frac)*mat_X[,'Isym_pos'] +
    gamma_a*mat_X[,'Iasym_pos']
  
  res = c(dS, dE, dIsym, dIasym, dHsub, dHcri, dD, dR, dS_pos, dE_pos, dIsym_pos, dIasym_pos, dR_pos)
  names(res) = varNames
  
  return(list(res))
}

  #seir_model_shields_rcfc_nolatent(0, X0, pars)




