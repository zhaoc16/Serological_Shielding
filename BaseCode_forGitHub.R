#General setup
rm(list=ls()) #Clear workspace
require(deSolve)
require(ggplot2)
require(reshape2)
require(cowplot)
require(here)
require(plyr)
require(beepr)

#Subscripts are defined as follows: c=children, a=non-essential adults, rc=reduced contact adults, fc=full contact adults, e=elderly
#For infection equations: The I compartments correspond to infectious periods, not periods in which people are symptomatic
#Is are cases that are severe enough to be documented eventually and Ia are undocumented cases
#Since we are not fitting to data we do not need to model specifically when the cases are documented 
#.pos indicates people who have tested positive for COVID-19 by antibody test

#Define transmission model ##############################
seir_model_shields_rcfc_nolatent <- function(t, x, params) {
  test.switch1<-sw1fxn(t)
  test.switch2<-sw2fxn(t)
  
  # children
  S.c = x[1]   #Susceptible children
  E.c = x[2]  #Latently infected children
  Is.c = x[3]   #Infectious children, severe enough to eventually become documented
  Ia.c = x[4]   #Infectious children, mild enough that does not ever become a documented case
  Hs.c = x[5]   #Hospitalized children, subcritical case
  Hc.c = x[6]   #Hospitalized children, critical case
  D.c = x[7]    #Child deaths
  R.c = x[8]    #Recovered children
  S.c.pos = x[9] #.pos indicates test positive, so this is susceptible children who have tested as positive
  E.c.pos = x[10]
  Is.c.pos = x[11]
  Ia.c.pos = x[12]
  R.c.pos = x[13]  
  
  # non-essential adults (.a )
  S.a = x[14]   
  E.a = x[15]
  Is.a = x[16]   
  Ia.a = x[17]   
  Hs.a = x[18]
  Hc.a = x[19]
  D.a = x[20]
  R.a = x[21]
  S.a.pos = x[22] #test positive
  E.a.pos = x[23]
  Is.a.pos = x[24]
  Ia.a.pos = x[25]
  R.a.pos = x[26]  
  
  #reduced contact adults (.rc )
  S.rc = x[27]   
  E.rc = x[28]
  Is.rc = x[29]   
  Ia.rc = x[30]   
  Hs.rc = x[31]
  Hc.rc = x[32]
  D.rc = x[33]
  R.rc = x[34]
  S.rc.pos = x[35]
  E.rc.pos = x[36]
  Is.rc.pos = x[37]
  Ia.rc.pos = x[38]
  R.rc.pos = x[39]  
  
  
  # full contact working adults (.fc)
  S.fc = x[40]   
  E.fc = x[41]
  Is.fc = x[42]   
  Ia.fc = x[43]   
  Hs.fc = x[44]
  Hc.fc = x[45]
  D.fc = x[46]
  R.fc = x[47]
  S.fc.pos = x[48]
  E.fc.pos = x[49]
  Is.fc.pos = x[50]
  Ia.fc.pos = x[51]
  R.fc.pos = x[52]
  
  # elderly (.e)
  S.e = x[53]   
  E.e = x[54]
  Is.e = x[55]   
  Ia.e = x[56]   
  Hs.e = x[57]
  Hc.e = x[58]
  D.e = x[59]
  R.e = x[60]
  S.e.pos = x[61]
  E.e.pos = x[62]
  Is.e.pos = x[63]
  Ia.e.pos = x[64]
  R.e.pos = x[65]  
  
  
  #Number of infections at time t for each subgroup
  infec.c.gen = asymp.red*Ia.c + Is.c
  infec.c.pos = asymp.red*Ia.c.pos + Is.c.pos
  infec.a.gen = asymp.red*Ia.a + Is.a 
  infec.a.pos = asymp.red*Ia.a.pos + Is.a.pos
  infec.rc.gen= asymp.red*Ia.rc + Is.rc 
  infec.rc.pos= asymp.red*Ia.rc.pos + Is.rc.pos
  infec.fc.gen= asymp.red*Ia.fc + Is.fc
  infec.fc.pos= asymp.red*Ia.fc.pos + Is.fc.pos
  infec.e.gen = asymp.red*Ia.e + Is.e  
  infec.e.pos = asymp.red*Ia.e.pos+ Is.e.pos 
  
  #Population sizes/testing status at time t    
  tot.c = sum(S.c, E.c, Is.c, Ia.c, Hs.c, Hc.c, D.c, R.c, 
              S.c.pos, E.c.pos, Is.c.pos, Ia.c.pos, R.c.pos)
  tot.a = sum(S.a, E.a, Is.a, Ia.a, Hs.a, Hc.a, D.a, R.a,
              S.a.pos, E.a.pos, Is.a.pos, Ia.c.pos, R.a.pos)
  tot.rc = sum(S.rc, E.rc, Is.rc, Ia.rc, Hs.rc, Hc.rc, D.rc, R.rc,
               S.rc.pos, E.rc.pos, Is.rc.pos, Ia.rc.pos, R.rc.pos)
  tot.fc = sum(S.fc, E.fc, Is.fc, Ia.fc, Hs.fc, Hc.fc, D.fc, R.fc, 
               S.fc.pos, E.fc.pos, Is.fc.pos, Ia.fc.pos, R.fc.pos)
  tot.e = sum(S.e, E.e, Is.e, Ia.e, Hs.e, Hc.e, D.e, R.e,
              S.e.pos, E.e.pos, Is.e.pos, Ia.e.pos, R.e.pos)
  
  tot.c.pos = sum(S.c.pos, E.c.pos, Is.c.pos, Ia.c.pos, R.c.pos)
  tot.a.pos = sum(S.a.pos, E.a.pos, Is.a.pos, Ia.a.pos, R.a.pos) 
  tot.rc.pos = sum(S.rc.pos, E.rc.pos, Is.rc.pos, Ia.rc.pos, R.rc.pos)
  tot.fc.pos = sum(S.fc.pos, E.fc.pos, Is.fc.pos, Ia.fc.pos, R.fc.pos)
  tot.e.pos = sum(S.e.pos, E.e.pos, Is.e.pos, Ia.e.pos, R.e.pos)
  
  tot.c.gen = tot.c - tot.c.pos
  tot.a.gen = tot.a - tot.a.pos
  tot.rc.gen = tot.rc - tot.rc.pos
  tot.fc.gen = tot.fc - tot.fc.pos
  tot.e.gen = tot.e - tot.e.pos
  
  #Fraction of population in each group who has tested positive by time
  FracReleased<-c(tot.c.pos/tot.c, tot.a.pos/tot.a, tot.rc.pos/tot.rc, 
                  tot.fc.pos/tot.fc, tot.e.pos/tot.e) 
  FracDistanced<-c(rep(1,5))-FracReleased
  
  #Correct for zeros, CANNOT MOVE THIS PART UP
  if(tot.c.pos==0){tot.c.pos <- tot.c}
  
  if(tot.a.pos==0){tot.a.pos <-tot.a}
  
  if(tot.rc.pos==0){tot.rc.pos <-tot.rc}
  
  if(tot.fc.pos==0){tot.fc.pos <-tot.fc}
  
  if(tot.e.pos==0){tot.e.pos <-tot.e}
  
  #Derive contact matrices
  
  #Matrix entry i,j is the number of contacts/day that an individual in group i has with group j
  #These are normalized for population size and corrected for symmetry using standard methods
  
  #Uncomment lines below for troubleshooting matrix changes
  #FracReleased<-c(0.01, 0.1, 0.2, 0.3, 0.05)
  #FracDistanced<-c(rep(1,5))-FracReleased
  
  #Home contacts
  HighContact_Home<-matrix(0, ncol=5, nrow=5)
  HighContact_Home[,1]<-HomeContacts_5x5[,1]*FracReleased[1]
  HighContact_Home[,2]<-HomeContacts_5x5[,2]*FracReleased[2]
  HighContact_Home[,3]<-HomeContacts_5x5[,3]*FracReleased[3]
  HighContact_Home[,4]<-HomeContacts_5x5[,4]*FracReleased[4]
  HighContact_Home[,5]<-HomeContacts_5x5[,5]*FracReleased[5]
  
  LowContact_Home<-matrix(0, ncol=5, nrow=5)
  LowContact_Home[,1]<-HomeContacts_5x5[,1]*FracDistanced[1]
  LowContact_Home[,2]<-HomeContacts_5x5[,2]*FracDistanced[2]
  LowContact_Home[,3]<-HomeContacts_5x5[,3]*FracDistanced[3]
  LowContact_Home[,4]<-HomeContacts_5x5[,4]*FracDistanced[4]
  LowContact_Home[,5]<-HomeContacts_5x5[,5]*FracDistanced[5]
  
  HomeContacts_10x10<-matrix(0, ncol=10, nrow=10)
  for(i in 1:10){
    rowfind<-c(rep(1:5, each=2))
    HomeVec<-c(HighContact_Home[rowfind[i],1], LowContact_Home[rowfind[i], 1], 
               HighContact_Home[rowfind[i],2], LowContact_Home[rowfind[i],2], 
               HighContact_Home[rowfind[i],3], LowContact_Home[rowfind[i],3],
               HighContact_Home[rowfind[i],4], LowContact_Home[rowfind[i], 4], 
               HighContact_Home[rowfind[i],5], LowContact_Home[rowfind[i], 5]) 
    HomeContacts_10x10[i,]<-HomeVec 
  }
  
  #School contacts
  
  #Baseline
  HighContact_School<-matrix(0, ncol=5, nrow=5)
  HighContact_School[,1]<-SchoolContacts_5x5[,1]*FracReleased[1]
  HighContact_School[,2]<-SchoolContacts_5x5[,2]*FracReleased[2]
  HighContact_School[,3]<-SchoolContacts_5x5[,3]*FracReleased[3]
  HighContact_School[,4]<-SchoolContacts_5x5[,4]*FracReleased[4]
  HighContact_School[,5]<-SchoolContacts_5x5[,5]*FracReleased[5]
  
  LowContact_School<-matrix(0, ncol=5, nrow=5)
  LowContact_School[,1]<-SchoolContacts_5x5[,1]*FracDistanced[1]
  LowContact_School[,2]<-SchoolContacts_5x5[,2]*FracDistanced[2]
  LowContact_School[,3]<-SchoolContacts_5x5[,3]*FracDistanced[3]
  LowContact_School[,4]<-SchoolContacts_5x5[,4]*FracDistanced[4]
  LowContact_School[,5]<-SchoolContacts_5x5[,5]*FracDistanced[5]
  
  SchoolContacts_10x10_Baseline<-matrix(0, ncol=10, nrow=10)
  for(i in 1:10){
    rowfind<-c(rep(1:5, each=2))
    SchoolVec<-c(HighContact_School[rowfind[i],1], LowContact_School[rowfind[i], 1], 
                 HighContact_School[rowfind[i],2], LowContact_School[rowfind[i],2], 
                 HighContact_School[rowfind[i],3], LowContact_School[rowfind[i],3],
                 HighContact_School[rowfind[i],4], LowContact_School[rowfind[i],4], 
                 HighContact_School[rowfind[i],5], LowContact_School[rowfind[i],5]) 
    SchoolContacts_10x10_Baseline[i,]<-SchoolVec 
  }
  
  #With general distancing, schools closed
  
  #Targeted distancing
  SchoolContacts_10x10_reopen<-matrix(0, nrow=10, ncol=10)
  for(i in 1:5){
    for(j in 1:5){
      SchoolContacts_10x10_reopen[i*2-1, j*2-1]<-SchoolContacts_5x5[i,j]
    }
  }
  
  #Work contacts
  
  #Baseline matrix
  HighContact_Work<-matrix(0, ncol=5, nrow=5)
  HighContact_Work[,1]<-WorkContacts_5x5[,1]*FracReleased[1]
  HighContact_Work[,2]<-WorkContacts_5x5[,2]*FracReleased[2]
  HighContact_Work[,3]<-WorkContacts_5x5[,3]*FracReleased[3]
  HighContact_Work[,4]<-WorkContacts_5x5[,4]*FracReleased[4]
  HighContact_Work[,5]<-WorkContacts_5x5[,5]*FracReleased[5]
  
  LowContact_Work<-matrix(0, ncol=5, nrow=5)
  LowContact_Work[,1]<-WorkContacts_5x5[,1]*FracDistanced[1]
  LowContact_Work[,2]<-WorkContacts_5x5[,2]*FracDistanced[2]
  LowContact_Work[,3]<-WorkContacts_5x5[,3]*FracDistanced[3]
  LowContact_Work[,4]<-WorkContacts_5x5[,4]*FracDistanced[4]
  LowContact_Work[,5]<-WorkContacts_5x5[,5]*FracDistanced[5]
  
  WorkContacts_10x10_Baseline<-matrix(0, ncol=10, nrow=10)
  for(i in 1:10){
    rowfind<-c(rep(1:5, each=2))
    WorkVec<-c(HighContact_Work[rowfind[i],1], LowContact_Work[rowfind[i],1], 
               HighContact_Work[rowfind[i],2], LowContact_Work[rowfind[i],2], 
               HighContact_Work[rowfind[i],3], LowContact_Work[rowfind[i],3],
               HighContact_Work[rowfind[i],4], LowContact_Work[rowfind[i],4], 
               HighContact_Work[rowfind[i],5], LowContact_Work[rowfind[i],5]) 
    WorkContacts_10x10_Baseline[i,]<-WorkVec 
  }
  
  #Broad distancing matrix
  WorkContacts_10x10_Distancing<-matrix(0, ncol=10, nrow=10)
  #Children have no work contacts
  #Non-essential workers don't work
  #Reduced contact workers reduce contacts by preduced and those contacts are proportional to 
  #prevalence in the population
  WorkContacts_10x10_Distancing[5,]<-c(WorkContacts_5x5[3,1]*FracReleased[1]*preduced, 
                                       WorkContacts_5x5[3,1]*FracDistanced[1]*preduced, 
                                       WorkContacts_5x5[3,2]*FracReleased[2]*preduced, WorkContacts_5x5[3,2]*FracDistanced[2]*preduced,
                                       WorkContacts_5x5[3,3]*FracReleased[3]*preduced, WorkContacts_5x5[3,3]*FracDistanced[3]*preduced,
                                       WorkContacts_5x5[3,4]*FracReleased[4]*preduced, WorkContacts_5x5[3,4]*FracDistanced[4]*preduced,
                                       WorkContacts_5x5[3,5]*FracReleased[5]*preduced, WorkContacts_5x5[3,5]*FracDistanced[5]*preduced) 
  WorkContacts_10x10_Distancing[6,]<-WorkContacts_10x10_Distancing[5,]
  #Full contact workers keep all of their baseline contacts
  #The distribution of those contacts is based on prevalence in the population
  WorkContacts_10x10_Distancing[7,]<-c(WorkContacts_5x5[4,1]*FracReleased[1]*pfull, WorkContacts_5x5[4,1]*FracDistanced[1]*pfull, 
                                       WorkContacts_5x5[4,2]*FracReleased[2]*pfull, WorkContacts_5x5[3,2]*FracDistanced[2]*pfull,
                                       WorkContacts_5x5[4,3]*FracReleased[3]*pfull, WorkContacts_5x5[3,3]*FracDistanced[3]*pfull,
                                       WorkContacts_5x5[4,4]*FracReleased[4]*pfull, WorkContacts_5x5[3,4]*FracDistanced[4]*pfull,
                                       WorkContacts_5x5[4,5]*FracReleased[5]*pfull, WorkContacts_5x5[3,5]*FracDistanced[5]*pfull)
  WorkContacts_10x10_Distancing[8,]<-WorkContacts_10x10_Distancing[7,]
  
  #Targeted distancing work
  #DISTANCING CONDITIONAL ON TESTING: People working from home regain their workplace contacts if they test positive 
  #People in reduced and full contact occupations have the distribution of their workplace contacts changed by alpha
  WorkContacts_10x10_TargetedDistancing<-matrix(0, nrow=10, ncol=10)
  WorkContacts_10x10_TargetedDistancing[3,]<-c(WorkContacts_5x5[2,1], 0, WorkContacts_5x5[2,2], 0, WorkContacts_5x5[2,3], 0,
                                               WorkContacts_5x5[2,4], 0, WorkContacts_5x5[2,5], 0)
  
  #We use fixed shielding (except we multiply by alpha and not alpha +1), 
  #so we correct for the situation where alpha times prevalence is greater than 1
  Scale<-numeric(5)
  Scale<-FracReleased*alpha
  for(i in 1:5){
    if(Scale[i]>1){Scale[i]<-1} 
  }
  
  WorkContacts_10x10_TargetedDistancing[1,]<-c(rep(0, 10))
  DistancedWork_RedContact<-WorkContacts_10x10_Baseline[5,]*preduced2
  WorkContacts_10x10_TargetedDistancing[5,1]<-sum(DistancedWork_RedContact[1], DistancedWork_RedContact[2])*Scale[1]
  WorkContacts_10x10_TargetedDistancing[5,2]<-sum(DistancedWork_RedContact[1], DistancedWork_RedContact[2])*(1-Scale[1])
  WorkContacts_10x10_TargetedDistancing[5,3]<-sum(DistancedWork_RedContact[3], DistancedWork_RedContact[4])*Scale[2]
  WorkContacts_10x10_TargetedDistancing[5,4]<-sum(DistancedWork_RedContact[3], DistancedWork_RedContact[4])*(1-Scale[2])
  WorkContacts_10x10_TargetedDistancing[5,5]<-sum(DistancedWork_RedContact[5], DistancedWork_RedContact[6])*Scale[3]
  WorkContacts_10x10_TargetedDistancing[5,6]<-sum(DistancedWork_RedContact[5], DistancedWork_RedContact[6])*(1-Scale[3])
  WorkContacts_10x10_TargetedDistancing[5,7]<-sum(DistancedWork_RedContact[7], DistancedWork_RedContact[8])*Scale[4]
  WorkContacts_10x10_TargetedDistancing[5,8]<-sum(DistancedWork_RedContact[7], DistancedWork_RedContact[8])*(1-Scale[4])
  WorkContacts_10x10_TargetedDistancing[5,9]<-sum(DistancedWork_RedContact[9], DistancedWork_RedContact[10])*Scale[5]
  WorkContacts_10x10_TargetedDistancing[5,10]<-sum(DistancedWork_RedContact[9], DistancedWork_RedContact[10])*(1-Scale[5])
  WorkContacts_10x10_TargetedDistancing[6,]<-WorkContacts_10x10_TargetedDistancing[5,]
  WorkContacts_10x10_TargetedDistancing[7,]<-WorkContacts_10x10_TargetedDistancing[5,]/preduced2
  WorkContacts_10x10_TargetedDistancing[8,]<-WorkContacts_10x10_TargetedDistancing[5,]/preduced2
  ##Other contacts
  
  #BASELINE OTHER CONTACTS MATRIX: People mix freely, contacts are based on proportion in the 
  #population
  HighContact_Other<-matrix(0, ncol=5, nrow=5)
  HighContact_Other[,1]<-OtherContacts_5x5[,1]*FracReleased[1]
  HighContact_Other[,2]<-OtherContacts_5x5[,2]*FracReleased[2]
  HighContact_Other[,3]<-OtherContacts_5x5[,3]*FracReleased[3]
  HighContact_Other[,4]<-OtherContacts_5x5[,4]*FracReleased[4]
  HighContact_Other[,5]<-OtherContacts_5x5[,5]*FracReleased[5]
  
  LowContact_Other<-matrix(0, ncol=5, nrow=5)
  LowContact_Other[,1]<-OtherContacts_5x5[,1]*FracDistanced[1]
  LowContact_Other[,2]<-OtherContacts_5x5[,2]*FracDistanced[2]
  LowContact_Other[,3]<-OtherContacts_5x5[,3]*FracDistanced[3]
  LowContact_Other[,4]<-OtherContacts_5x5[,4]*FracDistanced[4]
  LowContact_Other[,5]<-OtherContacts_5x5[,5]*FracDistanced[5]
  
  OtherContacts_10x10_Baseline<-matrix(0, ncol=10, nrow=10)
  for(i in 1:10){
    rowfind<-c(rep(1:5, each=2))
    OtherVec<-c(HighContact_Other[rowfind[i],1], LowContact_Other[rowfind[i], 1], 
                HighContact_Other[rowfind[i],2], LowContact_Other[rowfind[i],2], 
                HighContact_Other[rowfind[i],3], LowContact_Other[rowfind[i],3],
                HighContact_Other[rowfind[i],4], LowContact_Other[rowfind[i],4], 
                HighContact_Other[rowfind[i],5], LowContact_Other[rowfind[i],5]) 
    OtherContacts_10x10_Baseline[i,]<-OtherVec 
  }
  
  #UNTARGETED SOCIAL DISTANCING OTHER CONTACTS: All contacts for social distancing are reduced by 
  #a fraction
  OtherContacts_10x10_Distancing<-OtherContacts_10x10_Baseline*sd.other
  
  #TARGETED SOCIAL DISTANCING OTHER CONTACTS: Shielding
  #Contact probabilities proportional to prevalence in the population times a scaling factor to account for shielding
  #Keep the total contacts the same by subtracting from the total number of expected contacts between age classes
  #Once the scaling factor*prevalence exceeds 1, set to 1 (so that all contacts are with a positive individual)
  OtherContacts_10x10_TargetedDistancing<-matrix(0, nrow=10, ncol=10)
  
  for(i in 1:5){
    for(j in 1:5){
      EvenVec<-c(2, 4, 6, 8, 10)
      OddVec<-c(1,3,5,7,9)
      #Odd rows
      OtherContacts_10x10_TargetedDistancing[2*i-1,OddVec[j]]<-OtherContacts_5x5[i,j]*sd.other2*Scale[j]
      OtherContacts_10x10_TargetedDistancing[2*i-1,EvenVec[j]]<-OtherContacts_5x5[i,j]*sd.other2*(1-Scale[j])
      #Even rows
      OtherContacts_10x10_TargetedDistancing[2*i,OddVec[j]]<-OtherContacts_5x5[i,j]*sd.other2*Scale[j]
      OtherContacts_10x10_TargetedDistancing[2*i,EvenVec[j]]<-OtherContacts_5x5[i,j]*sd.other2*(1-Scale[j])
    }
  }
  
  
  #For test positive individuals, keep the distribution of other contacts shielded but change the number to pre-pandemic levels
  OtherContacts_10x10_TargetedDistancing2<-OtherContacts_10x10_TargetedDistancing
  OtherContacts_10x10_TargetedDistancing2[1,]<-OtherContacts_10x10_TargetedDistancing[1,]*(1/sd.other2)
  OtherContacts_10x10_TargetedDistancing2[3,]<-OtherContacts_10x10_TargetedDistancing[3,]*(1/sd.other2)
  OtherContacts_10x10_TargetedDistancing2[5,]<-OtherContacts_10x10_TargetedDistancing[5,]*(1/sd.other2)
  OtherContacts_10x10_TargetedDistancing2[7,]<-OtherContacts_10x10_TargetedDistancing[7,]*(1/sd.other2)
  OtherContacts_10x10_TargetedDistancing2[9,]<-OtherContacts_10x10_TargetedDistancing[9,]*(1/sd.other2)
  
  
  #Change matrices used over time
  if(t<start.time.distance){ #Use baseline matrices until social distancing starts
    CM<-HomeContacts_10x10+SchoolContacts_10x10_Baseline+WorkContacts_10x10_Baseline+
      OtherContacts_10x10_Baseline 
  }
  #if(t>=start.time.distance & t<time.reopen){ #Use this line if changing between distancing scenarios
  if(t>=start.time.distance & t<start.time.target){ #Use these matrices under general social distancing without testing
    CM<-HomeContacts_10x10 + WorkContacts_10x10_Distancing + OtherContacts_10x10_Distancing
  }
  if(t>=start.time.target & t<school.start.time){ #Use these matrices between the start of targeting and schools reopening
    CM<-HomeContacts_10x10 + WorkContacts_10x10_TargetedDistancing + 
      OtherContacts_10x10_TargetedDistancing2}
  if(t>=school.start.time){ #Keep most matrices the same but add schools back in
    CM<-HomeContacts_10x10 + WorkContacts_10x10_TargetedDistancing + 
      OtherContacts_10x10_TargetedDistancing2+SchoolContacts_10x10_Baseline
  }
  
  #Force of infection by group 
  #Children
  foi.c.pos = q*((CM[1,1]*infec.c.pos/(tot.c.pos)) + (CM[1,2]*infec.c.gen/(tot.c.gen))+
                   (CM[1,3]*infec.a.pos/(tot.a.pos)) + (CM[1,4]*infec.a.gen/(tot.a.gen))+
                   (CM[1,5]*infec.rc.pos/(tot.rc.pos)) + (CM[1,6]*infec.rc.gen/(tot.rc.gen))+
                   (CM[1,7]*infec.fc.pos/(tot.fc.pos)) + (CM[1,8]*infec.fc.gen/(tot.fc.gen))+
                   (CM[1,9]*infec.e.pos/(tot.e.pos)) + (CM[1,10]*infec.e.gen/(tot.e.gen)))
  
  foi.c.gen = q*((CM[2,1]*infec.c.pos/(tot.c.pos)) + (CM[2,2]*infec.c.gen/(tot.c.gen))+
                   (CM[2,3]*infec.a.pos/(tot.a.pos)) + (CM[2,4]*infec.a.gen/(tot.a.gen))+
                   (CM[2,5]*infec.rc.pos/(tot.rc.pos)) + (CM[2,6]*infec.rc.gen/(tot.rc.gen))+
                   (CM[2,7]*infec.fc.pos/(tot.fc.pos)) + (CM[2,8]*infec.fc.gen/(tot.fc.gen))+
                   (CM[2,9]*infec.e.pos/(tot.e.pos)) + (CM[2,10]*infec.e.gen/(tot.e.gen)))
  
  #Homebound adults
  foi.a.pos = q*((CM[3,1]*infec.c.pos/(tot.c.pos)) + (CM[3,2]*infec.c.gen/(tot.c.gen))+
                   (CM[3,3]*infec.a.pos/(tot.a.pos)) + (CM[3,4]*infec.a.gen/(tot.a.gen))+
                   (CM[3,5]*infec.rc.pos/(tot.rc.pos)) + (CM[3,6]*infec.rc.gen/(tot.rc.gen))+
                   (CM[3,7]*infec.fc.pos/(tot.fc.pos)) + (CM[3,8]*infec.fc.gen/(tot.fc.gen))+
                   (CM[3,9]*infec.e.pos/(tot.e.pos)) + (CM[3,10]*infec.e.gen/(tot.e.gen)))
  
  foi.a.gen = q*((CM[4,1]*infec.c.pos/(tot.c.pos)) + (CM[4,2]*infec.c.gen/(tot.c.gen))+
                   (CM[4,3]*infec.a.pos/(tot.a.pos)) + (CM[4,4]*infec.a.gen/(tot.a.gen))+
                   (CM[4,5]*infec.rc.pos/(tot.rc.pos)) + (CM[4,6]*infec.rc.gen/(tot.rc.gen))+
                   (CM[4,7]*infec.fc.pos/(tot.fc.pos)) + (CM[4,8]*infec.fc.gen/(tot.fc.gen))+
                   (CM[4,9]*infec.e.pos/(tot.e.pos)) + (CM[4,10]*infec.e.gen/(tot.e.gen)))
  
  #Reduced contact adults
  foi.rc.pos = q*((CM[5,1]*infec.c.pos/(tot.c.pos)) + (CM[5,2]*infec.c.gen/(tot.c.gen))+
                    (CM[5,3]*infec.a.pos/(tot.a.pos)) + (CM[5,4]*infec.a.gen/(tot.a.gen))+
                    (CM[5,5]*infec.rc.pos/(tot.rc.pos)) + (CM[5,6]*infec.rc.gen/(tot.rc.gen))+
                    (CM[5,7]*infec.fc.pos/(tot.fc.pos)) + (CM[5,8]*infec.fc.gen/(tot.fc.gen))+
                    (CM[5,9]*infec.e.pos/(tot.e.pos)) + (CM[5,10]*infec.e.gen/(tot.e.gen)))
  
  foi.rc.gen = q*((CM[6,1]*infec.c.pos/(tot.c.pos)) + (CM[6,2]*infec.c.gen/(tot.c.gen))+
                    (CM[6,3]*infec.a.pos/(tot.a.pos)) + (CM[6,4]*infec.a.gen/(tot.a.gen))+
                    (CM[6,5]*infec.rc.pos/(tot.rc.pos)) + (CM[6,6]*infec.rc.gen/(tot.rc.gen))+
                    (CM[6,7]*infec.fc.pos/(tot.fc.pos)) + (CM[6,8]*infec.fc.gen/(tot.fc.gen))+
                    (CM[6,9]*infec.e.pos/(tot.e.pos)) + (CM[6,10]*infec.e.gen/(tot.e.gen)))
  
  #Full contact adults
  foi.fc.pos = q*((CM[7,1]*infec.c.pos/(tot.c.pos)) + (CM[7,2]*infec.c.gen/(tot.c.gen))+
                    (CM[7,3]*infec.a.pos/(tot.a.pos)) + (CM[7,4]*infec.a.gen/(tot.a.gen))+
                    (CM[7,5]*infec.rc.pos/(tot.rc.pos)) + (CM[7,6]*infec.rc.gen/(tot.rc.gen))+
                    (CM[7,7]*infec.fc.pos/(tot.fc.pos)) + (CM[7,8]*infec.fc.gen/(tot.fc.gen))+
                    (CM[7,9]*infec.e.pos/(tot.e.pos)) + (CM[7,10]*infec.e.gen/(tot.e.gen)))
  
  foi.fc.gen = q*((CM[8,1]*infec.c.pos/(tot.c.pos)) + (CM[8,2]*infec.c.gen/(tot.c.gen))+
                    (CM[8,3]*infec.a.pos/(tot.a.pos)) + (CM[8,4]*infec.a.gen/(tot.a.gen))+
                    (CM[8,5]*infec.rc.pos/(tot.rc.pos)) + (CM[8,6]*infec.rc.gen/(tot.rc.gen))+
                    (CM[8,7]*infec.fc.pos/(tot.fc.pos)) + (CM[8,8]*infec.fc.gen/(tot.fc.gen))+
                    (CM[8,9]*infec.e.pos/(tot.e.pos)) + (CM[8,10]*infec.e.gen/(tot.e.gen)))
  
  #Elderly
  foi.e.pos = q*((CM[9,1]*infec.c.pos/(tot.c.pos)) + (CM[9,2]*infec.c.gen/(tot.c.gen))+
                   (CM[9,3]*infec.a.pos/(tot.a.pos)) + (CM[9,4]*infec.a.gen/(tot.a.gen))+
                   (CM[9,5]*infec.rc.pos/(tot.rc.pos)) + (CM[9,6]*infec.rc.gen/(tot.rc.gen))+
                   (CM[9,7]*infec.fc.pos/(tot.fc.pos)) + (CM[9,8]*infec.fc.gen/(tot.fc.gen))+
                   (CM[9,9]*infec.e.pos/(tot.e.pos)) + (CM[9,10]*infec.e.gen/(tot.e.gen)))
  
  foi.e.gen = q*((CM[10,1]*infec.c.pos/(tot.c.pos)) + (CM[10,2]*infec.c.gen/(tot.c.gen))+
                   (CM[10,3]*infec.a.pos/(tot.a.pos)) + (CM[10,4]*infec.a.gen/(tot.a.gen))+
                   (CM[10,5]*infec.rc.pos/(tot.rc.pos)) + (CM[10,6]*infec.rc.gen/(tot.rc.gen))+
                   (CM[10,7]*infec.fc.pos/(tot.fc.pos)) + (CM[10,8]*infec.fc.gen/(tot.fc.gen))+
                   (CM[10,9]*infec.e.pos/(tot.e.pos)) + (CM[10,10]*infec.e.gen/(tot.e.gen)))
  
  #Testing
  
  #Get number of tests available by group
  child.tests<-daily.tests*agestruc[1]
  adult.tests.h<-daily.tests*agestruc[2]*prob.home
  adult.tests.rc<-daily.tests*agestruc[2]*prob.reduced
  adult.tests.fc<-daily.tests*agestruc[2]*prob.full
  el.tests<-daily.tests*agestruc[3]
  
  #Divide by number of people eligible to be tested to get proportion tested per day
  test.c<-child.tests/sum(S.c, E.c, Ia.c, Hs.c, Hc.c, R.c)
  if(test.c>1){test.c<-1} #When prevalence and testing are high, might not be enough people eligible to be tested so correct this here
  test.a<-adult.tests.h/sum(S.a, E.a, Ia.a, Hs.a, Hc.a, R.a)
  if(test.a>1){test.a<-1}
  test.rc<-adult.tests.rc/sum(S.rc, E.rc, Ia.rc, Hs.rc, Hc.rc, R.rc)
  if(test.rc>1){test.rc<-1}
  test.fc<-adult.tests.fc/sum(S.fc, E.fc, Ia.fc, Hs.fc, Hc.fc, R.fc)
  if(test.fc>1){test.fc<-1}
  test.e<-el.tests/sum(S.e, E.e, Ia.e, Hs.e, Hc.e, R.e)
  if(test.e>1){test.e<-1}
  
  #Model equations
  dS.c = - foi.c.gen*S.c - (1-specificity)*test.c*test.switch2*S.c 
  dE.c = foi.c.gen*S.c - gamma_e*E.c - (1-specificity)*test.c*test.switch2*E.c
  dIs.c = gamma_e*E.c*p - gamma_s*Is.c 
  dIa.c = gamma_e*E.c*(1-p) - gamma_a*Ia.c - (1-specificity)*test.c*test.switch2*Ia.c
  dHs.c = gamma_s*Is.c*(hosp_frac[1]-hosp_crit[1]) + gamma_s*Is.c.pos*(hosp_frac[1]-hosp_crit[1]) - gamma_h*Hs.c
  dHc.c = gamma_s*Is.c*hosp_crit[1] + gamma_s*Is.c.pos*hosp_crit[1]- gamma_h*Hc.c
  dD.c = gamma_h*Hc.c*crit_die[1]
  dR.c = (1-hosp_frac[1])*gamma_s*Is.c + gamma_a*Ia.c + (1-sensitivity)*test.switch2*gamma_h*Hs.c + 
    (1-sensitivity)*gamma_h*test.switch2*Hc.c*(1-crit_die[1])  - sensitivity*test.c*test.switch2*R.c+
    gamma_h*test.switch1*Hc.c*(1-crit_die[1])+test.switch1*gamma_h*Hs.c 
  dS.c.pos = (1-specificity)*test.c*test.switch2*S.c - foi.c.pos*S.c.pos 
  dE.c.pos = (1-specificity)*test.c*test.switch2*E.c + foi.c.pos*S.c.pos - gamma_e*E.c.pos 
  dIs.c.pos = gamma_e*E.c.pos*p - gamma_s*Is.c.pos
  dIa.c.pos = (1-specificity)*test.c*test.switch2*Ia.c + gamma_e*E.c.pos*(1-p) - gamma_a*Ia.c.pos
  dR.c.pos = sensitivity*test.c*test.switch2*R.c + sensitivity*gamma_h*test.switch2*Hc.c*(1-crit_die[1]) + sensitivity*gamma_h*test.switch2*Hs.c +
    (1-hosp_frac[1])*gamma_s*Is.c.pos + gamma_a*Ia.c.pos
  
  dS.a = - foi.a.gen*S.a - (1-specificity)*test.a*test.switch2*S.a 
  dE.a = foi.a.gen*S.a - gamma_e*E.a - (1-specificity)*test.a*test.switch2*E.a
  dIs.a = gamma_e*E.a*p - gamma_s*Is.a 
  dIa.a = gamma_e*E.a*(1-p) - gamma_a*Ia.a - (1-specificity)*test.a*test.switch2*Ia.a
  dHs.a = gamma_s*Is.a*(hosp_frac[2]-hosp_crit[2]) + gamma_s*Is.a.pos*(hosp_frac[2]-hosp_crit[2]) - gamma_h*Hs.a
  dHc.a = gamma_s*Is.a*hosp_crit[2] + gamma_s*Is.a.pos*hosp_crit[2] - gamma_h*Hc.a
  dD.a = gamma_h*Hc.a*crit_die[2]
  dR.a = (1-hosp_frac[2])*gamma_s*Is.a + gamma_a*Ia.a + (1-sensitivity)*test.switch2*gamma_h*Hs.a + 
    (1-sensitivity)*test.switch2*gamma_h*Hc.a*(1-crit_die[2])  - sensitivity*test.a*test.switch2*R.a+
    gamma_h*test.switch1*Hc.a*(1-crit_die[2])+test.switch1*gamma_h*Hs.a 
  dS.a.pos = (1-specificity)*test.a*test.switch2*S.a - foi.a.pos*S.a.pos 
  dE.a.pos = (1-specificity)*test.a*test.switch2*E.a + foi.a.pos*S.a.pos - gamma_e*E.a.pos
  dIs.a.pos = gamma_e*E.a.pos*p - gamma_s*Is.a.pos
  dIa.a.pos = (1-specificity)*test.a*test.switch2*Ia.a + gamma_e*E.a.pos*(1-p) - gamma_a*Ia.a.pos
  dR.a.pos = sensitivity*test.a*test.switch2*R.a + sensitivity*test.switch2*gamma_h*Hc.a*(1-crit_die[2]) + sensitivity*test.switch2*gamma_h*Hs.a +
    (1-hosp_frac[2])*gamma_s*Is.a.pos + gamma_a*Ia.a.pos
  
  dS.rc = - foi.rc.gen*S.rc - (1-specificity)*test.rc*test.switch2*S.rc 
  dE.rc = foi.rc.gen*S.rc - gamma_e*E.rc - (1-specificity)*test.rc*test.switch2*E.rc
  dIs.rc = gamma_e*E.rc*p - gamma_s*Is.rc 
  dIa.rc = gamma_e*E.rc*(1-p) - gamma_a*Ia.rc - (1-specificity)*test.rc*test.switch2*Ia.rc
  dHs.rc = gamma_s*Is.rc*(hosp_frac[2]-hosp_crit[2]) + gamma_s*Is.rc.pos*(hosp_frac[2]-hosp_crit[2]) - gamma_h*Hs.rc
  dHc.rc = gamma_s*Is.rc*hosp_crit[2] + gamma_s*Is.rc.pos*hosp_crit[2]- gamma_h*Hc.rc
  dD.rc = gamma_h*Hc.rc*crit_die[2]
  dR.rc = (1-hosp_frac[2])*gamma_s*Is.rc + gamma_a*Ia.rc + (1-sensitivity)*test.switch2*gamma_h*Hs.rc + 
    (1-sensitivity)*test.switch2*gamma_h*Hc.rc*(1-crit_die[2])  - sensitivity*test.rc*test.switch2*R.rc+
    gamma_h*test.switch1*Hc.rc*(1-crit_die[2])+test.switch1*gamma_h*Hs.rc 
  dS.rc.pos = (1-specificity)*test.rc*test.switch2*S.rc - foi.rc.pos*S.rc.pos 
  dE.rc.pos = (1-specificity)*test.rc*test.switch2*E.rc + foi.rc.pos*S.rc.pos - gamma_e*E.rc.pos 
  dIs.rc.pos = gamma_e*E.rc.pos*p - gamma_s*Is.rc.pos
  dIa.rc.pos = (1-specificity)*test.rc*test.switch2*Ia.rc + gamma_e*E.rc.pos*(1-p) - gamma_a*Ia.rc.pos
  dR.rc.pos = sensitivity*test.rc*test.switch2*R.rc + sensitivity*test.switch2*gamma_h*Hc.rc*(1-crit_die[2]) + sensitivity*test.switch2*gamma_h*Hs.rc +
    (1-hosp_frac[2])*gamma_s*Is.rc.pos + gamma_a*Ia.rc.pos 
  
  
  dS.fc = - foi.fc.gen*S.fc - (1-specificity)*test.fc*test.switch2*S.fc 
  dE.fc = foi.fc.gen*S.fc - gamma_e*E.fc - (1-specificity)*test.fc*test.switch2*E.fc
  dIs.fc = gamma_e*E.fc*p - gamma_s*Is.fc 
  dIa.fc = gamma_e*E.fc*(1-p) - gamma_a*Ia.fc - (1-specificity)*test.fc*test.switch2*Ia.fc
  dHs.fc = gamma_s*Is.fc*(hosp_frac[2]-hosp_crit[2]) + gamma_s*Is.fc.pos*(hosp_frac[2]-hosp_crit[2]) - gamma_h*Hs.fc
  dHc.fc = gamma_s*Is.fc*hosp_crit[2] + gamma_s*Is.fc.pos*hosp_crit[2]- gamma_h*Hc.fc
  dD.fc = gamma_h*Hc.fc*crit_die[2]
  dR.fc = (1-hosp_frac[2])*gamma_s*Is.fc + gamma_a*Ia.fc + (1-sensitivity)*test.switch2*gamma_h*Hs.fc + 
    (1-sensitivity)*test.switch2*gamma_h*Hc.fc*(1-crit_die[2])  - sensitivity*test.fc*test.switch2*R.fc +
    gamma_h*test.switch1*Hc.fc*(1-crit_die[2])+test.switch1*gamma_h*Hs.fc 
  dS.fc.pos = (1-specificity)*test.fc*test.switch2*S.fc - foi.fc.pos*S.fc.pos 
  dE.fc.pos = (1-specificity)*test.fc*test.switch2*E.fc + foi.fc.pos*S.fc.pos - gamma_e*E.fc.pos 
  dIs.fc.pos = gamma_e*E.fc.pos*p - gamma_s*Is.fc.pos
  dIa.fc.pos = (1-specificity)*test.fc*test.switch2*Ia.fc + gamma_e*E.fc.pos*(1-p) - gamma_a*Ia.fc.pos
  dR.fc.pos = sensitivity*test.fc*test.switch2*R.fc + sensitivity*test.switch2*gamma_h*Hc.fc*(1-crit_die[2]) + sensitivity*test.switch2*gamma_h*Hs.fc +
    (1-hosp_frac[2])*gamma_s*Is.fc.pos + gamma_a*Ia.fc.pos
  
  dS.e = - foi.e.gen*S.e - (1-specificity)*test.e*test.switch2*S.e 
  dE.e = foi.e.gen*S.e - gamma_e*E.e - (1-specificity)*test.e*test.switch2*E.e
  dIs.e = gamma_e*E.e*p - gamma_s*Is.e 
  dIa.e = gamma_e*E.e*(1-p) - gamma_a*Ia.e - (1-specificity)*test.e*test.switch2*Ia.e
  dHs.e = gamma_s*Is.e*(hosp_frac[3]-hosp_crit[3]) + gamma_s*Is.e.pos*(hosp_frac[3]-hosp_crit[3]) - gamma_h*Hs.e
  dHc.e = gamma_s*Is.e*hosp_crit[3] + gamma_s*Is.e.pos*hosp_crit[3]- gamma_h*Hc.e
  dD.e = gamma_h*Hc.e*crit_die[3]
  dR.e = (1-hosp_frac[3])*gamma_s*Is.a + gamma_a*Ia.e + (1-sensitivity)*test.switch2*gamma_h*Hs.e + 
    (1-sensitivity)*test.switch2*gamma_h*Hc.e*(1-crit_die[3])  - sensitivity*test.e*test.switch2*R.e+
    gamma_h*test.switch1*Hc.e*(1-crit_die[3])+test.switch1*gamma_h*Hs.e 
  dS.e.pos = (1-specificity)*test.e*test.switch2*S.e - foi.e.pos*S.e.pos 
  dE.e.pos = (1-specificity)*test.e*test.switch2*E.e + foi.e.pos*S.e.pos - gamma_e*E.e.pos 
  dIs.e.pos = gamma_e*E.e.pos*p - gamma_s*Is.e.pos
  dIa.e.pos = (1-specificity)*test.e*test.switch2*Ia.e + gamma_e*E.e.pos*(1-p) - gamma_a*Ia.e.pos
  dR.e.pos = sensitivity*test.e*test.switch2*R.e + sensitivity*gamma_h*test.switch2*Hc.e*(1-crit_die[3]) + sensitivity*gamma_h*test.switch2*Hs.e +
    (1-hosp_frac[3])*gamma_s*Is.e.pos + gamma_a*Ia.e.pos
  
  
  
  res = c(dS.c, dE.c, dIs.c, dIa.c, dHs.c, dHc.c, dD.c, dR.c, dS.c.pos, dE.c.pos, dIs.c.pos, dIa.c.pos, dR.c.pos,
          dS.a, dE.a, dIs.a, dIa.a, dHs.a, dHc.a, dD.a, dR.a, dS.a.pos, dE.a.pos, dIs.a.pos, dIa.a.pos, dR.a.pos,
          dS.rc, dE.rc, dIs.rc, dIa.rc, dHs.rc, dHc.rc, dD.rc, dR.rc, 
          dS.rc.pos, dE.rc.pos, dIs.rc.pos, dIa.rc.pos, dR.rc.pos,
          dS.fc, dE.fc, dIs.fc, dIa.fc, dHs.fc, dHc.fc, dD.fc, dR.fc, 
          dS.fc.pos, dE.fc.pos, dIs.fc.pos, dIa.fc.pos, dR.fc.pos,
          dS.e, dE.e, dIs.e, dIa.e, dHs.e, dHc.e, dD.e, dR.e, dS.c.pos, dE.e.pos, dIs.e.pos, dIa.e.pos, dR.e.pos)
  
  list(res)
}

#Model parameters

# Population
N=323*10^6
agefrac.0=c(0.12,0.13,0.13,0.13,0.13,0.13,0.11,0.06,0.04,0.02) # from Weitz model
agestruc=c(sum(agefrac.0[1:2]), sum(agefrac.0[3:6], 0.5*agefrac.0[7]), 
           sum(0.5*agefrac.0[7], agefrac.0[8:10]))

#Matric parameters
#3x3 data from Prem et al
AllContacts<-matrix(data=c(9.75, 2.57, 0.82, 5.97, 10.32, 2.25, 0.39, 0.46, 1.20), nrow=3)
WorkContacts<-matrix(data=c(0.20, 0.28, 0, 0.64, 4.73, 0, 0, 0, 0), nrow=3)
SchoolContacts<-matrix(data=c(4.32, 0.47, 0.02, 1.10, 0.32, 0.04, 0.01, 0.01, 0.03), nrow=3)
HomeContacts<-matrix(data=c(2.03, 1.02, 0.50, 2.37, 1.82, 0.68, 0.24, 0.14, 0.62), nrow=3)
OtherContacts<-matrix(data=c(3.20, 0.80, 0.30, 1.86, 3.45, 1.53, 0.14, 0.32, 0.55), nrow=3)
#Expand the 3x3 matrix to a 5x5 based on the fraction of the population in each of the worker
#subgroups
expand_5x5<-function(oldmatrix, phome, preduced, pfull){
  newmatrix<-matrix(NA, ncol=5, nrow=5)
  newmatrix[1,1]<-oldmatrix[1,1]
  newmatrix[2:4,1]<-oldmatrix[2,1]
  newmatrix[5,1]<-oldmatrix[3,1]
  newmatrix[1,5]<-oldmatrix[1,3]
  newmatrix[2:4,5]<-oldmatrix[2,3]
  newmatrix[5,5]<-oldmatrix[3,3]
  vec1<-c(oldmatrix[1,2]*phome, oldmatrix[1,2]*preduced, oldmatrix[1,2]*pfull)
  newmatrix[1,2:4]<-vec1
  vec2<-c(oldmatrix[2,2]*phome, oldmatrix[2,2]*preduced, oldmatrix[2,2]*pfull)
  newmatrix[2,2:4]<-vec2
  newmatrix[3,2:4]<-vec2
  newmatrix[4,2:4]<-vec2
  vec3<-c(oldmatrix[3,2]*phome, oldmatrix[3,2]*preduced, oldmatrix[3,2]*pfull)
  newmatrix[5,2:4]<-vec3
  newmatrix
}

prob.home<-0.316
#prob.home<-0.4
prob.full<-0.0565
prob.reduced<-1-prob.full-prob.home
#prob.reduced<-0.678
preduced<-0.5
pfull<-1
sd.other<-0.25
alpha<-1.2

WorkContacts_5x5<-expand_5x5(oldmatrix=WorkContacts, 
                             phome=prob.home, preduced=prob.reduced, pfull=prob.full)
HomeContacts_5x5<-expand_5x5(oldmatrix=HomeContacts, 
                             phome=prob.home, preduced=prob.reduced, pfull=prob.full)
SchoolContacts_5x5<-expand_5x5(oldmatrix=SchoolContacts, 
                               phome=prob.home, preduced=prob.reduced, pfull=prob.full)
OtherContacts_5x5<-expand_5x5(oldmatrix=OtherContacts, 
                              phome=prob.home, preduced=prob.reduced, pfull=prob.full)

#Natural history/transmission parameters
#Note on R0: with base structure 63.28q
R0=2.9
q=R0/63.28   #probability of transmission from children
asymp.red=0.55 #relative infectiousness of asymptomatic infections compared to 
#symptomatic infections


gamma_e=1/3     # Latent period (He et al)
gamma_a=1/7     # Recovery rate, undocumented (Kissler et al)
gamma_s=1/7    # Recovery rate, undocumented (Kissler et al)
gamma_h=1/15    # Recovery rate, hospitalized cases (Zhou et al--China study)
p=0.14           # Fraction 'Symptomatic'documented' (Shaman's paper)

hosp_frac=c(0.061, 0.182, 0.417) #From MMWR
hosp_crit=c(0, 0.063, 0.173) #From CDC, MMWR
crit_die=c(0, 0.5, 0.5) #Obtained from initial fitting

##TESTING
#For Cellex
sensitivity=0.94
specificity=0.96

testvec<-rep(0, 366)
times <- seq(0, 365, length = 366)

#Initial conditions
initcase.c=60
initcase.a=20
initcase.e=40
initcase.rc=50
initcase.fc=1

# base model
start = c(S.c = agestruc[1]*N - initcase.c, E.c = 0, Is.c = initcase.c, Ia.c = 0, Hs.c = 0, Hc.c = 0, D.c = 0, R.c = 0, 
          S.c.pos = 0, E.c.pos = 0, Is.c.pos = 0, Ia.c.pos = 0, R.c.pos =0,
          
          S.a = agestruc[2]*N*prob.home - initcase.a, E.a = 0, Is.a = initcase.a, Ia.a = 0, Hs.a = 0, Hc.a = 0, D.a = 0, R.a = 0, 
          S.a.pos = 0, E.a.pos = 0, Is.a.pos = 0, Ia.a.pos = 0, R.a.pos =0,
          
          S.rc = agestruc[2]*N*prob.reduced-initcase.rc, E.rc = 0, Is.rc = initcase.rc, Ia.rc = 0, Hs.rc = 0, Hc.rc = 0, D.rc = 0, R.rc = 0, 
          S.rc.pos = 0, E.rc.pos = 0, Is.rc.pos = 0, Ia.rc.pos = 0, R.rc.pos =0,
          
          S.fc = agestruc[2]*N*prob.full-initcase.fc, E.fc = 0, Is.fc = initcase.fc, Ia.fc = 0, Hs.fc = 0, Hc.fc = 0, D.fc = 0, R.fc = 0, 
          S.fc.pos = 0, E.fc.pos = 0, Is.fc.pos = 0, Ia.fc.pos = 0, R.fc.pos =0,
          S.e = agestruc[3]*N - initcase.e, E.e = 0, Is.e = initcase.e, Ia.e = 0, Hs.e = 0, Hc.e= 0, D.e = 0, R.e = 0, 
          S.e.pos = 0, E.e.pos = 0, Is.e.pos = 0, Ia.e.pos = 0, R.e.pos =0)

t = (0:365)

#Example simulation: do nothing
start.time.distance<-500
sd.other<-0.25
sd.other2<-sd.other
preduced<-0.1
preduced2<-preduced
start.time.test<-364 #can change when in the outbreak testing becomes available
start.time.target<-500
time.reopen<-500
times <- seq(0, 365, length = 366)
daily.tests<-0
school.start.time<-500

#Might want a ramp up period for these, ignore for now
tswitch1.dat<-data.frame(times=times, test.switch1=c(rep(1, 366)))
tswitch2.dat<-data.frame(times=times, test.switch2=c(rep(0, 366)))

#Use these lines if testing is used
#tswitch1.dat<-data.frame(times=times, test.switch1=c(rep(1, start.time.test+1),
#                                                     rep(0, length(times)-(start.time.test+1))))
#tswitch2.dat<-data.frame(times=times, test.switch2=c(rep(0, start.time.test+1),
#                                                     rep(1, length(times)-(start.time.test+1))))

sw1fxn<-approxfun(tswitch1.dat$times, tswitch1.dat$test.switch1, rule=2)
sw2fxn<-approxfun(tswitch2.dat$times, tswitch2.dat$test.switch2, rule=2)

params = c('agestruc'=agestruc, 
           'HomeContacts_5x5'=HomeContacts_5x5, 'WorkContacts_5x5'=WorkContacts_5x5, 'SchoolContacts_5x5'=SchoolContacts_5x5, 
           'OtherContacts_5x5'=OtherContacts_5x5, 
           'preduced'=preduced, 'pfull'=pfull, 'sd.other'=sd.other,
           'preduced2'=preduced, 'sd.other2'=sd.other,
           'alpha'=alpha, 'daily.tests'=daily.tests,
           'start.time.distance'=start.time.distance, 
           'time.reopen'=time.reopen,
           #'start.time.target'=start.time.target, 
           #'school.start.time'=school.start.time,
           'gamma_e'=gamma_e, 'gamma_a'=gamma_a, 'gamma_s'=gamma_s, 'gamma_h'=gamma_h, 
           'p'=p, 'hosp_frac'=hosp_frac, 'hosp_crit'=hosp_crit, 'crit_die'=crit_die,
           'sensitivity'=sensitivity, 'specificity'=specificity) 

model_out = ode(y = start, times = t, fun = seir_model_shields_rcfc_nolatent, parms = par, 
                method='ode45')

model_out<-as.data.frame(model_out)
Infected<-N-sum(model_out$S.c[366], model_out$S.c.pos[366], model_out$S.a[366], model_out$S.a.pos[366], 
                model_out$S.rc[366], model_out$S.rc.pos[366], model_out$S.fc[366], model_out$S.fc.pos[366],
                model_out$S.e[366], model_out$S.e.pos[366])
Pinfect<-Infected/N
Infected; Pinfect
Deaths<-sum(model_out$D.a[366], model_out$D.c[366], model_out$D.e[366], model_out$D.rc[366], 
            model_out$D.fc[366])
Deaths

null.model<-ddply(model_out, .(time), summarize, 
                  CriticalCare=sum(Hc.fc, Hc.rc, Hc.c, Hc.a, Hc.e, Hc.c),
                  Deaths=sum(D.c, D.a, D.e, D.rc, D.fc),
                  CI=(N-sum(S.c, S.c.pos, S.a, S.a.pos, S.rc, S.rc.pos, S.fc, S.fc.pos, S.e, S.e.pos))/N)
plot(null.model$CriticalCare)

null.model$Deaths[60]; null.model$Deaths[90]