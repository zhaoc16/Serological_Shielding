#Two tier social mixing model with code only
#Setup######
require(deSolve)
require(ggplot2)
require(reshape2)
require(cowplot)
require(here)
require(plyr)

#Define transmission model ##############################
seir_model_shields_rcfc <- function(t, x, RM_p) {
  with(as.list(c(RM_p, x)), {
    #Time-varying testing functions
    test.c<-cfxn(t)
    test.a<-afxn(t)
    test.rc<-rcfxn(t)
    test.fc <- fcfxn(t)
    test.e<-efxn(t)
    test.switch1<-sw1fxn(t)
    test.switch2<-sw2fxn(t)
    
    # children
    S.c = x[1]   
    E.c = x[2]
    Is.c = x[3]   
    Ia.c = x[4]   
    Hs.c = x[5]
    Hc.c = x[6]
    D.c = x[7]
    R.c = x[8]
    S.c.pos = x[9]
    E.c.pos = x[10]
    Is.c.pos = x[11]
    Ia.c.pos = x[12]
    R.c.pos = x[13]  
    
    # adults
    S.a = x[14]   
    E.a = x[15]
    Is.a = x[16]   
    Ia.a = x[17]   
    Hs.a = x[18]
    Hc.a = x[19]
    D.a = x[20]
    R.a = x[21]
    S.a.pos = x[22]
    E.a.pos = x[23]
    Is.a.pos = x[24]
    Ia.a.pos = x[25]
    R.a.pos = x[26]  
    
    #reduced contact adults
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
    
    
    # full contact working adults
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
    
    # elderly
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
    infec.c.gen = asymp.red*(E.c + Ia.c) + Is.c
    infec.c.pos = asymp.red*(E.c.pos + Ia.c.pos) + Is.c.pos
    infec.a.gen = asymp.red*(E.a + Ia.a) + Is.a 
    infec.a.pos = asymp.red*(E.a.pos + Ia.a.pos) + Is.a.pos
    infec.rc.gen= asymp.red*(E.rc + Ia.rc) + Is.rc 
    infec.rc.pos= asymp.red*(E.rc.pos + Ia.rc.pos) + Is.rc.pos
    infec.fc.gen= asymp.red*(E.fc + Ia.fc) + Is.fc
    infec.fc.pos= asymp.red*(E.fc.pos + Ia.fc.pos) + Is.fc.pos
    infec.e.gen = asymp.red*(E.e + Ia.e) + Is.e  
    infec.e.pos = asymp.red*(E.e.pos + Ia.e.pos)+ Is.e.pos 
    
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
    tot.a.pos = sum(S.a.pos, E.a.pos, Is.a.pos, Ia.c.pos, R.a.pos)
    tot.rc.pos = sum(S.rc.pos, E.rc.pos, Is.rc.pos, Ia.rc.pos, R.rc.pos)
    tot.fc.pos = sum(S.fc.pos, E.fc.pos, Is.fc.pos, Ia.fc.pos, R.fc.pos)
    tot.e.pos = sum(S.e.pos, E.e.pos, Is.e.pos, Ia.e.pos, R.e.pos)
    
    tot.c.gen = tot.c - tot.c.pos
    tot.a.gen = tot.a - tot.a.pos
    tot.rc.gen = tot.rc - tot.rc.pos
    tot.fc.gen = tot.fc - tot.fc.pos
    tot.e.gen = tot.e - tot.e.pos
    
    #Fractions positive by time
    FracReleased<-c(tot.c.pos/tot.c, tot.a.pos/tot.a, tot.rc.pos/tot.rc, 
                    tot.fc.pos/tot.fc, tot.e.pos/tot.e) 
    FracDistanced<-c(rep(1,5))-FracReleased
    
    #Correct for zeros, CANNOT MOVE THIS PART UP
    if(sum(S.c.pos, E.c.pos, Is.c.pos, Ia.c.pos, R.c.pos)==0){
      tot.c.pos <- tot.c}
    
    if(sum(S.a.pos, E.a.pos, Is.a.pos, Ia.c.pos, R.a.pos)==0){
      tot.a.pos <-tot.a}
    
    if(sum(S.rc.pos, E.rc.pos, Is.rc.pos, Ia.rc.pos, R.rc.pos)==0){
      tot.rc.pos <-tot.rc}
    
    if(sum(S.fc.pos, E.fc.pos, Is.fc.pos, Ia.fc.pos, R.fc.pos)==0){
      tot.fc.pos <-tot.fc}
    
    if(sum(S.e.pos, E.e.pos, Is.e.pos, Ia.e.pos, R.e.pos)==0){
      tot.e.pos <-tot.e}
    
    #Derive matrices 
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
    #Reduced contact workers reduce contacts by 50% and those contacts proportional to 
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
    WorkContacts_10x10_Distancing[8,]<-WorkContacts_10x10_Distancing[8,]
    
    #Targeted distancing work
    #DISTANCING CONDITIONAL ON TESTING: People working from home regain their workplace contacts if they test positive 
    #Only test positive individuals are allowed to work in reduced or full contact positions, and those contacts are proportional
    #to population prevalence because they are 'forward-facing' occupations.
    #Could consider other scenarios (like the need to have more reduced/full contact workers active than 
    #simply those that test positive)
    #We could also consider having people who go back to work have fewer workplace contacts based on 
    #the number of people allowed to go back to work and also activities within workplace
    #venues
    WorkContacts_10x10_TargetedDistancing<-matrix(0, nrow=10, ncol=10)
    for(i in 1:5){
      for(j in 1:5){
        WorkContacts_10x10_TargetedDistancing[i*2-1, j*2-1]<-WorkContacts_5x5[i,j]
      }
    }
    
    WorkContacts_10x10_TargetedDistancing[5,]<-c(WorkContacts_5x5[3,1]*FracReleased[1], WorkContacts_5x5[3,1]*FracDistanced[1],
                                                 WorkContacts_5x5[3,2]*FracReleased[2], WorkContacts_5x5[3,2]*FracDistanced[2],
                                                 WorkContacts_5x5[3,3]*FracReleased[3], WorkContacts_5x5[3,3]*FracDistanced[3],
                                                 WorkContacts_5x5[3,4]*FracReleased[4], WorkContacts_5x5[3,4]*FracDistanced[4],
                                                 WorkContacts_5x5[3,5]*FracReleased[5], WorkContacts_5x5[3,5]*FracDistanced[5])
    WorkContacts_10x10_TargetedDistancing[7,]<-c(WorkContacts_5x5[4,1]*FracReleased[1], WorkContacts_5x5[4,1]*FracDistanced[1],
                                                 WorkContacts_5x5[4,2]*FracReleased[2], WorkContacts_5x5[4,2]*FracDistanced[2],
                                                 WorkContacts_5x5[4,3]*FracReleased[3], WorkContacts_5x5[4,3]*FracDistanced[3],
                                                 WorkContacts_5x5[4,4]*FracReleased[4], WorkContacts_5x5[4,4]*FracDistanced[4],
                                                 WorkContacts_5x5[4,5]*FracReleased[5], WorkContacts_5x5[4,5]*FracDistanced[5])
    
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
    
    #TARGETED SOCIAL DISTANCING OTHER CONTACTS: Flexible shielding
    OtherContacts_10x10_TargetedDistancing<-matrix(0, nrow=10, ncol=10)
    DistancedOther_5x5<-OtherContacts_5x5*sd.other
    LowContact_OtherTargeted<-DistancedOther_5x5/(alpha+1)
    HighContact_OtherTargeted<-LowContact_OtherTargeted*alpha
    
    for(i in 1:10){
      rowfind<-c(rep(1:5, each=2))
      OtherVecTarget<-c(HighContact_OtherTargeted[rowfind[i],1], LowContact_OtherTargeted[rowfind[i], 1], 
                        HighContact_OtherTargeted[rowfind[i],2], LowContact_OtherTargeted[rowfind[i],2], 
                        HighContact_OtherTargeted[rowfind[i],3], LowContact_OtherTargeted[rowfind[i],3],
                        HighContact_OtherTargeted[rowfind[i],4], LowContact_OtherTargeted[rowfind[i],4], 
                        HighContact_OtherTargeted[rowfind[i],5], LowContact_OtherTargeted[rowfind[i],5]) 
      OtherContacts_10x10_TargetedDistancing[i,]<-OtherVecTarget 
    }
    
    OtherContacts_10x10_TargetedDistancing2<-OtherContacts_10x10_TargetedDistancing
    OtherContacts_10x10_TargetedDistancing2[1,]<-OtherContacts_10x10_TargetedDistancing[1,]*(1/sd.other)
    OtherContacts_10x10_TargetedDistancing2[3,]<-OtherContacts_10x10_TargetedDistancing[3,]*(1/sd.other)
    OtherContacts_10x10_TargetedDistancing2[5,]<-OtherContacts_10x10_TargetedDistancing[5,]*(1/sd.other)
    OtherContacts_10x10_TargetedDistancing2[7,]<-OtherContacts_10x10_TargetedDistancing[7,]*(1/sd.other)
    OtherContacts_10x10_TargetedDistancing2[9,]<-OtherContacts_10x10_TargetedDistancing[9,]*(1/sd.other)
    
    
    if(t<start.time.distance){
      CM<-HomeContacts_10x10+SchoolContacts_10x10_Baseline+WorkContacts_10x10_Baseline+
        OtherContacts_10x10_Baseline 
    }
    if(t>=start.time.distance & t<start.time.target){
      CM<-HomeContacts_10x10 + WorkContacts_10x10_Distancing + OtherContacts_10x10_Distancing
    }
    if(t>=start.time.target){
      CM<-HomeContacts_10x10 + SchoolContacts_10x10_reopen + WorkContacts_10x10_TargetedDistancing + OtherContacts_10x10_TargetedDistancing2}
    
    
    #Force of infection by group 
    #Children
    foi.c.pos = q*((CM[1,1]*infec.c.pos/(tot.c.pos)) + (CM[1,2]*infec.c.gen/(tot.c.gen))+
                     (CM[1,3]*infec.a.pos/(tot.a.pos)) + (CM[1,4]*infec.a.gen/(tot.a.gen))+
                     (CM[1,5]*infec.rc.pos/(tot.rc.pos)) + (CM[1,6]*infec.rc.gen/(tot.rc.gen))+
                     (CM[1,7]*infec.fc.pos/(tot.fc.pos)) + (CM[1,8]*infec.fc.gen/(tot.fc.gen))+
                     (CM[1,9]*infec.a.pos/(tot.a.pos)) + (CM[1,10]*infec.a.gen/(tot.a.gen)))
    
    foi.c.gen = q*((CM[2,1]*infec.c.pos/(tot.c.pos)) + (CM[2,2]*infec.c.gen/(tot.c.gen))+
                     (CM[2,3]*infec.a.pos/(tot.a.pos)) + (CM[2,4]*infec.a.gen/(tot.a.gen))+
                     (CM[2,5]*infec.rc.pos/(tot.rc.pos)) + (CM[2,6]*infec.rc.gen/(tot.rc.gen))+
                     (CM[2,7]*infec.fc.pos/(tot.fc.pos)) + (CM[2,8]*infec.fc.gen/(tot.fc.gen))+
                     (CM[2,9]*infec.a.pos/(tot.a.pos)) + (CM[2,10]*infec.a.gen/(tot.a.gen)))
    
    #Homebound adults
    foi.a.pos = q*((CM[3,1]*infec.c.pos/(tot.c.pos)) + (CM[3,2]*infec.c.gen/(tot.c.gen))+
                     (CM[3,3]*infec.a.pos/(tot.a.pos)) + (CM[3,4]*infec.a.gen/(tot.a.gen))+
                     (CM[3,5]*infec.rc.pos/(tot.rc.pos)) + (CM[3,6]*infec.rc.gen/(tot.rc.gen))+
                     (CM[3,7]*infec.fc.pos/(tot.fc.pos)) + (CM[3,8]*infec.fc.gen/(tot.fc.gen))+
                     (CM[3,9]*infec.a.pos/(tot.a.pos)) + (CM[3,10]*infec.a.gen/(tot.a.gen)))
    
    foi.a.gen = q*((CM[4,1]*infec.c.pos/(tot.c.pos)) + (CM[4,2]*infec.c.gen/(tot.c.gen))+
                     (CM[4,3]*infec.a.pos/(tot.a.pos)) + (CM[4,4]*infec.a.gen/(tot.a.gen))+
                     (CM[4,5]*infec.rc.pos/(tot.rc.pos)) + (CM[4,6]*infec.rc.gen/(tot.rc.gen))+
                     (CM[4,7]*infec.fc.pos/(tot.fc.pos)) + (CM[4,8]*infec.fc.gen/(tot.fc.gen))+
                     (CM[4,9]*infec.a.pos/(tot.a.pos)) + (CM[4,10]*infec.a.gen/(tot.a.gen)))
    
    #Reduced contact adults
    foi.rc.pos = q*((CM[5,1]*infec.c.pos/(tot.c.pos)) + (CM[5,2]*infec.c.gen/(tot.c.gen))+
                      (CM[5,3]*infec.a.pos/(tot.a.pos)) + (CM[5,4]*infec.a.gen/(tot.a.gen))+
                      (CM[5,5]*infec.rc.pos/(tot.rc.pos)) + (CM[5,6]*infec.rc.gen/(tot.rc.gen))+
                      (CM[5,7]*infec.fc.pos/(tot.fc.pos)) + (CM[5,8]*infec.fc.gen/(tot.fc.gen))+
                      (CM[5,9]*infec.a.pos/(tot.a.pos)) + (CM[5,10]*infec.a.gen/(tot.a.gen)))
    
    foi.rc.gen = q*((CM[6,1]*infec.c.pos/(tot.c.pos)) + (CM[6,2]*infec.c.gen/(tot.c.gen))+
                      (CM[6,3]*infec.a.pos/(tot.a.pos)) + (CM[6,4]*infec.a.gen/(tot.a.gen))+
                      (CM[6,5]*infec.rc.pos/(tot.rc.pos)) + (CM[6,6]*infec.rc.gen/(tot.rc.gen))+
                      (CM[6,7]*infec.fc.pos/(tot.fc.pos)) + (CM[6,8]*infec.fc.gen/(tot.fc.gen))+
                      (CM[6,9]*infec.a.pos/(tot.a.pos)) + (CM[6,10]*infec.a.gen/(tot.a.gen)))
    
    #Full contact adults
    foi.fc.pos = q*((CM[7,1]*infec.c.pos/(tot.c.pos)) + (CM[7,2]*infec.c.gen/(tot.c.gen))+
                      (CM[7,3]*infec.a.pos/(tot.a.pos)) + (CM[7,4]*infec.a.gen/(tot.a.gen))+
                      (CM[7,5]*infec.rc.pos/(tot.rc.pos)) + (CM[7,6]*infec.rc.gen/(tot.rc.gen))+
                      (CM[7,7]*infec.fc.pos/(tot.fc.pos)) + (CM[7,8]*infec.fc.gen/(tot.fc.gen))+
                      (CM[7,9]*infec.a.pos/(tot.a.pos)) + (CM[7,10]*infec.a.gen/(tot.a.gen)))
    
    foi.fc.gen = q*((CM[8,1]*infec.c.pos/(tot.c.pos)) + (CM[8,2]*infec.c.gen/(tot.c.gen))+
                      (CM[8,3]*infec.a.pos/(tot.a.pos)) + (CM[8,4]*infec.a.gen/(tot.a.gen))+
                      (CM[8,5]*infec.rc.pos/(tot.rc.pos)) + (CM[8,6]*infec.rc.gen/(tot.rc.gen))+
                      (CM[8,7]*infec.fc.pos/(tot.fc.pos)) + (CM[8,8]*infec.fc.gen/(tot.fc.gen))+
                      (CM[8,9]*infec.a.pos/(tot.a.pos)) + (CM[8,10]*infec.a.gen/(tot.a.gen)))
    
    #Elderly
    foi.e.pos = q*((CM[9,1]*infec.c.pos/(tot.c.pos)) + (CM[9,2]*infec.c.gen/(tot.c.gen))+
                     (CM[9,3]*infec.a.pos/(tot.a.pos)) + (CM[9,4]*infec.a.gen/(tot.a.gen))+
                     (CM[9,5]*infec.rc.pos/(tot.rc.pos)) + (CM[9,6]*infec.rc.gen/(tot.rc.gen))+
                     (CM[9,7]*infec.fc.pos/(tot.fc.pos)) + (CM[9,8]*infec.fc.gen/(tot.fc.gen))+
                     (CM[9,9]*infec.a.pos/(tot.a.pos)) + (CM[9,10]*infec.a.gen/(tot.a.gen)))
    
    foi.e.gen = q*((CM[10,1]*infec.c.pos/(tot.c.pos)) + (CM[10,2]*infec.c.gen/(tot.c.gen))+
                     (CM[10,3]*infec.a.pos/(tot.a.pos)) + (CM[10,4]*infec.a.gen/(tot.a.gen))+
                     (CM[10,5]*infec.rc.pos/(tot.rc.pos)) + (CM[10,6]*infec.rc.gen/(tot.rc.gen))+
                     (CM[10,7]*infec.fc.pos/(tot.fc.pos)) + (CM[10,8]*infec.fc.gen/(tot.fc.gen))+
                     (CM[10,9]*infec.a.pos/(tot.a.pos)) + (CM[10,10]*infec.a.gen/(tot.a.gen)))
    
    #Model equations
    dS.c = - foi.c.gen*S.c - (1-specificity)*test.c*S.c 
    dE.c = foi.c.gen*S.c - gamma_e*E.c - (1-specificity)*test.c*E.c
    dIs.c = gamma_e*E.c*p - gamma_s*Is.c 
    dIa.c = gamma_e*E.c*(1-p) - gamma_a*Ia.c - (1-specificity)*test.c*Ia.c
    dHs.c = gamma_s*Is.c*hosp_frac[1]*(1-hosp_crit[1]) + gamma_s*Is.c.pos*hosp_frac[1]*(1-hosp_crit[1]) - gamma_h*Hs.c
    dHc.c = gamma_s*Is.c*hosp_frac[1]*hosp_crit[1] + gamma_s*Is.c.pos*hosp_frac[1]*hosp_crit[1]- gamma_h*Hc.c
    dD.c = gamma_h*Hc.c*crit_die
    dR.c = (1-hosp_frac[1])*gamma_s*Is.c + gamma_a*Ia.c + (1-sensitivity)*test.switch2*gamma_h*Hs.c + 
      (1-sensitivity)*gamma_h*test.switch2*Hc.c*(1-crit_die)  - sensitivity*test.c*R.c+
      gamma_h*test.switch1*Hc.c*(1-crit_die)+test.switch1*gamma_h*Hs.c 
    dS.c.pos = (1-specificity)*test.c*S.c - foi.c.pos*S.c.pos 
    dE.c.pos = (1-specificity)*test.c*E.c + foi.c.pos*S.c.pos - gamma_e*E.c.pos 
    dIs.c.pos = gamma_e*E.c.pos*p - gamma_s*Is.c.pos
    dIa.c.pos = (1-specificity)*test.c*Ia.c + gamma_e*E.c.pos*(1-p) - gamma_a*Ia.c.pos
    dR.c.pos = sensitivity*test.c*R.c + sensitivity*gamma_h*test.switch2*Hc.c*(1-crit_die) + sensitivity*gamma_h*test.switch2*Hs.c +
      (1-hosp_frac[1])*gamma_s*Is.c.pos + gamma_a*Ia.c.pos
    
    dS.a = - foi.a.gen*S.a - (1-specificity)*test.a*S.a 
    dE.a = foi.a.gen*S.a - gamma_e*E.a - (1-specificity)*test.a*E.a
    dIs.a = gamma_e*E.a*p - gamma_s*Is.a 
    dIa.a = gamma_e*E.a*(1-p) - gamma_a*Ia.a - (1-specificity)*test.a*Ia.a
    dHs.a = gamma_s*Is.a*hosp_frac[2]*(1-hosp_crit[2]) + gamma_s*Is.a.pos*hosp_frac[2]*(1-hosp_crit[2]) - gamma_h*Hs.a
    dHc.a = gamma_s*Is.a*hosp_frac[2]*hosp_crit[2] + gamma_s*Is.a.pos*hosp_frac[2]*hosp_crit[2] - gamma_h*Hc.a
    dD.a = gamma_h*Hc.a*crit_die
    dR.a = (1-hosp_frac[2])*gamma_s*Is.a + gamma_a*Ia.a + (1-sensitivity)*test.switch2*gamma_h*Hs.a + 
      (1-sensitivity)*test.switch2*gamma_h*Hc.a*(1-crit_die)  - sensitivity*test.a*R.a+
      gamma_h*test.switch1*Hc.a*(1-crit_die)+test.switch1*gamma_h*Hs.a 
    dS.a.pos = (1-specificity)*test.a*S.a - foi.a.pos*S.a.pos 
    dE.a.pos = (1-specificity)*test.a*E.a + foi.a.pos*S.a.pos - gamma_e*E.a.pos
    dIs.a.pos = gamma_e*E.a.pos*p - gamma_s*Is.a.pos
    dIa.a.pos = (1-specificity)*test.a*Ia.a + gamma_e*E.a.pos*(1-p) - gamma_a*Ia.a.pos
    dR.a.pos = sensitivity*test.a*R.a + sensitivity*test.switch2*gamma_h*Hc.a*(1-crit_die) + sensitivity*test.switch2*gamma_h*Hs.a +
      (1-hosp_frac[2])*gamma_s*Is.a.pos + gamma_a*Ia.a.pos
    
    dS.rc = - foi.rc.gen*S.rc - (1-specificity)*test.rc*S.rc 
    dE.rc = foi.rc.gen*S.rc - gamma_e*E.rc - (1-specificity)*test.rc*E.rc
    dIs.rc = gamma_e*E.rc*p - gamma_s*Is.rc 
    dIa.rc = gamma_e*E.rc*(1-p) - gamma_a*Ia.rc - (1-specificity)*test.rc*Ia.rc
    dHs.rc = gamma_s*Is.rc*hosp_frac[2]*(1-hosp_crit[2]) + gamma_s*Is.rc.pos*hosp_frac[2]*(1-hosp_crit[2]) - gamma_h*Hs.rc
    dHc.rc = gamma_s*Is.rc*hosp_frac[2]*hosp_crit[2] + gamma_s*Is.rc.pos*hosp_frac[2]*hosp_crit[2]- gamma_h*Hc.rc
    dD.rc = gamma_h*Hc.rc*crit_die
    dR.rc = (1-hosp_frac[2])*gamma_s*Is.rc + gamma_a*Ia.rc + (1-sensitivity)*test.switch2*gamma_h*Hs.rc + 
      (1-sensitivity)*test.switch2*gamma_h*Hc.rc*(1-crit_die)  - sensitivity*test.rc*R.rc+
      gamma_h*test.switch1*Hc.rc*(1-crit_die)+test.switch1*gamma_h*Hs.rc 
    dS.rc.pos = (1-specificity)*test.rc*S.rc - foi.rc.pos*S.rc.pos 
    dE.rc.pos = (1-specificity)*test.rc*E.rc + foi.rc.pos*S.rc.pos - gamma_e*E.rc.pos 
    dIs.rc.pos = gamma_e*E.rc.pos*p - gamma_s*Is.rc.pos
    dIa.rc.pos = (1-specificity)*test.rc*Ia.rc + gamma_e*E.rc.pos*(1-p) - gamma_a*Ia.rc.pos
    dR.rc.pos = sensitivity*test.rc*R.rc + sensitivity*test.switch2*gamma_h*Hc.rc*(1-crit_die) + sensitivity*test.switch2*gamma_h*Hs.rc +
      (1-hosp_frac[2])*gamma_s*Is.rc.pos + gamma_a*Ia.rc.pos 
    
    
    dS.fc = - foi.fc.gen*S.fc - (1-specificity)*test.fc*S.fc 
    dE.fc = foi.fc.gen*S.fc - gamma_e*E.fc - (1-specificity)*test.fc*E.fc
    dIs.fc = gamma_e*E.fc*p - gamma_s*Is.fc 
    dIa.fc = gamma_e*E.fc*(1-p) - gamma_a*Ia.fc - (1-specificity)*test.fc*Ia.fc
    dHs.fc = gamma_s*Is.fc*hosp_frac[2]*(1-hosp_crit[2]) + gamma_s*Is.fc.pos*hosp_frac[2]*(1-hosp_crit[2]) - gamma_h*Hs.fc
    dHc.fc = gamma_s*Is.fc*hosp_frac[2]*hosp_crit[2] + gamma_s*Is.fc.pos*hosp_frac[2]*hosp_crit[2]- gamma_h*Hc.fc
    dD.fc = gamma_h*Hc.fc*crit_die
    dR.fc = (1-hosp_frac[2])*gamma_s*Is.fc + gamma_a*Ia.fc + (1-sensitivity)*test.switch2*gamma_h*Hs.fc + 
      (1-sensitivity)*test.switch2*gamma_h*Hc.fc*(1-crit_die)  - sensitivity*test.c*R.fc +
      gamma_h*test.switch1*Hc.fc*(1-crit_die)+test.switch1*gamma_h*Hs.fc 
    dS.fc.pos = (1-specificity)*test.fc*S.fc - foi.fc.pos*S.fc.pos 
    dE.fc.pos = (1-specificity)*test.fc*E.fc + foi.fc.pos*S.fc.pos - gamma_e*E.fc.pos 
    dIs.fc.pos = gamma_e*E.fc.pos*p - gamma_s*Is.fc.pos
    dIa.fc.pos = (1-specificity)*test.fc*Ia.fc + gamma_e*E.fc.pos*(1-p) - gamma_a*Ia.fc.pos
    dR.fc.pos = sensitivity*test.fc*R.fc + sensitivity*test.switch2*gamma_h*Hc.fc*(1-crit_die) + sensitivity*test.switch2*gamma_h*Hs.fc +
      (1-hosp_frac[2])*gamma_s*Is.fc.pos + gamma_a*Ia.fc.pos
    
    dS.e = - foi.e.gen*S.e - (1-specificity)*test.e*S.e 
    dE.e = foi.e.gen*S.e - gamma_e*E.e - (1-specificity)*test.e*E.e
    dIs.e = gamma_e*E.e*p - gamma_s*Is.e 
    dIa.e = gamma_e*E.e*(1-p) - gamma_a*Ia.e - (1-specificity)*test.e*Ia.e
    dHs.e = gamma_s*Is.e*hosp_frac[3]*(1-hosp_crit[3]) + gamma_s*Is.e.pos*hosp_frac[3]*(1-hosp_crit[3]) - gamma_h*Hs.e
    dHc.e = gamma_s*Is.e*hosp_frac[3]*hosp_crit[3] + gamma_s*Is.e.pos*hosp_frac[3]*hosp_crit[3]- gamma_h*Hc.e
    dD.e = gamma_h*Hc.e*crit_die
    dR.e = (1-hosp_frac[3])*gamma_s*Is.a + gamma_a*Ia.e + (1-sensitivity)*test.switch2*gamma_h*Hs.e + 
      (1-sensitivity)*test.switch2*gamma_h*Hc.e*(1-crit_die)  - sensitivity*test.e*R.e+
      gamma_h*test.switch1*Hc.e*(1-crit_die)+test.switch1*gamma_h*Hs.e 
    dS.e.pos = (1-specificity)*test.e*S.e - foi.e.pos*S.e.pos 
    dE.e.pos = (1-specificity)*test.e*E.e + foi.e.pos*S.e.pos - gamma_e*E.e.pos 
    dIs.e.pos = gamma_e*E.e.pos*p - gamma_s*Is.e.pos
    dIa.e.pos = (1-specificity)*test.e*Ia.e + gamma_e*E.e.pos*(1-p) - gamma_a*Ia.e.pos
    dR.e.pos = sensitivity*test.e*R.e + sensitivity*gamma_h*Hc.e*(1-crit_die) + sensitivity*gamma_h*Hs.e +
      (1-hosp_frac[3])*gamma_s*Is.e.pos + gamma_a*Ia.e.pos
    
    
    
    res = c(dS.c, dE.c, dIs.c, dIa.c, dHs.c, dHc.c, dD.c, dR.c, dS.c.pos, dE.c.pos, dIs.c.pos, dIa.c.pos, dR.c.pos,
            dS.a, dE.a, dIs.a, dIa.a, dHs.a, dHc.a, dD.a, dR.a, dS.a.pos, dE.a.pos, dIs.a.pos, dIa.a.pos, dR.a.pos,
            dS.rc, dE.rc, dIs.rc, dIa.rc, dHs.rc, dHc.rc, dD.rc, dR.rc, 
            dS.rc.pos, dE.rc.pos, dIs.rc.pos, dIa.rc.pos, dR.rc.pos,
            dS.fc, dE.fc, dIs.fc, dIa.fc, dHs.fc, dHc.fc, dD.fc, dR.fc, 
            dS.fc.pos, dE.fc.pos, dIs.fc.pos, dIa.fc.pos, dR.fc.pos,
            dS.e, dE.e, dIs.e, dIa.e, dHs.e, dHc.e, dD.e, dR.e, dS.c.pos, dE.e.pos, dIs.e.pos, dIa.e.pos, dR.e.pos)
    
    list(res)
  })
}

# Population
N=331*10^6
agefrac.0=c(0.12,0.13,0.13,0.13,0.13,0.13,0.11,0.06,0.04,0.02) # from Weitz model
agestruc=c(sum(agefrac.0[1:2]), sum(agefrac.0[3:6]), sum(agefrac.0[7:10]))

###PART 1: Input 3x3 matrices###

AllContacts<-matrix(data=c(9.75, 2.57, 0.82, 5.97, 10.32, 2.25, 0.39, 0.46, 1.20), nrow=3)
WorkContacts<-matrix(data=c(0.20, 0.28, 0, 0.64, 4.73, 0, 0, 0, 0), nrow=3)
SchoolContacts<-matrix(data=c(4.32, 0.47, 0.02, 1.10, 0.32, 0.04, 0.01, 0.01, 0.03), nrow=3)
HomeContacts<-matrix(data=c(2.03, 1.02, 0.50, 2.37, 1.82, 0.68, 0.24, 0.14, 0.62), nrow=3)
OtherContacts<-matrix(data=c(3.20, 0.80, 0.30, 1.86, 3.45, 1.53, 0.14, 0.32, 0.55), nrow=3)


###PART 2: Expand 3x3 matrices to 5x5 matrices######

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

prob.home<-0.3
prob.reduced<-0.64
prob.full<-0.06
preduced<-0.5
pfull<-1
sd.other<-0.5
alpha<-1.2

WorkContacts_5x5<-expand_5x5(oldmatrix=WorkContacts, 
                             phome=prob.home, preduced=prob.reduced, pfull=prob.full)
HomeContacts_5x5<-expand_5x5(oldmatrix=HomeContacts, 
                             phome=prob.home, preduced=prob.reduced, pfull=prob.full)
SchoolContacts_5x5<-expand_5x5(oldmatrix=SchoolContacts, 
                               phome=prob.home, preduced=prob.reduced, pfull=prob.full)
OtherContacts_5x5<-expand_5x5(oldmatrix=OtherContacts, 
                              phome=prob.home, preduced=prob.reduced, pfull=prob.full)

#Note on R0: with base structure 119.525q=R0
q = 0.025    #probability of transmission from children
asymp.red=0.5 #relative infectiousness of asymptomatic infections compared to 
#symptomatic infections


# Natural History Parameters
gamma_e=1/4     # Latent period
gamma_a=1/6     # Recovery rate, asymptomatic
gamma_s=1/10    # Recovery rate, symptomatic
gamma_h=1/10    # Recovery rate, hostpitalized cases
p=0.5           # Fraction Symptomatic
#pars$p=c(0.99, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.5, 0.5) # Structured, from Weitz model

#beta_a=2.5/10   # Transmission for asymptomatic # from Weitz model
#beta_s=5/10     # Transmission for symptomatic # from Weitz model

#hosp_frac.0=c(0.1, 0.3, 1.2, 3.2, 4.9, 10.2, 16.6, 24.3, 27.3, 27.3)/100 # from Weitz model, collapse for wider age groups
hosp_frac=c(0.002, 0.049, 0.21) 
#hosp_frac=c(0.061, 0.182, 0.417) #From CDC and New York, probably biased high

#hosp_crit.0=c(5, 5, 5, 5, 6.3, 12.2, 27.4, 43.2, 70.9, 70.9)/100 # from Weitz model, collapse for wider age groups
hosp_crit=c(0.05, 0.07, 0.43)

crit_die=0.5

#For Wantai total antibody ELISA, make specificity slightly less than 1
sensitivity=0.90
specificity=0.999

#Pulse functions
testvec<-rep(0, 366)
times <- seq(0, 365, length = 366)

start.time.test<-45 #can change when in the outbreak testing becomes available
start.time.distance<-30
start.time.target<-50

#adults
testvec.a<-testvec
freq.a<-7 #once started, frequency of testing
index.a<-seq(from=start.time.test, to=366, by=freq.a) #make a vector with time indices for testing
frac.a<-1 #what fraction to test at each time point


#children
testvec.c<-testvec
freq.c<-7 #once started, frequency of testing
index.c<-seq(from=start.time.test, to=366, by=freq.c) #make a vector with time indices for testing
frac.c<-1 

#elderly
testvec.e<-testvec
freq.e<-7 #once started, frequency of testing
index.e<-seq(from=start.time.test, to=366, by=freq.e) #make a vector with time indices for testing
frac.e<-1 #what fraction/rate to test at each time point...technically a rate, so can be over 1

#adult, reduced contact
testvec.rc<-testvec
freq.rc<-7 #once started, frequency of testing
index.rc<-seq(from=start.time.test, to=366, by=freq.rc) #make a vector with time indices for testing
frac.rc<-1 #what fraction/rate to test at each time point...technically a rate so can be over 1

#adult, full contact
testvec.fc<-testvec
freq.fc<-7 #once started, frequency of testing
index.fc<-seq(from=start.time.test, to=366, by=freq.fc) #make a vector with time indices for testing
frac.fc<-1 #what fraction/rate to test at each time point

for(i in 1:366){
  if (i %in% unique(index.a)){testvec.a[i]<-frac.a}
  if (i %in% unique(index.c)){testvec.c[i]<-frac.c}
  if (i %in% unique(index.e)){testvec.e[i]<-frac.e}
  if (i %in% unique(index.rc)){testvec.rc[i]<-frac.rc}
  if (i %in% unique(index.fc)){testvec.fc[i]<-frac.fc}
}

c.testdat<-data.frame(times=times, test.c=testvec.c)
a.testdat<-data.frame(times=times, test.a=testvec.a)
rc.testdat<-data.frame(times=times, test.rc=testvec.rc)
fc.testdat<-data.frame(times=times, test.fc=testvec.fc)
e.testdat<-data.frame(times=times, test.e=testvec.e)
#Might want a ramp up period for these, ignore for now
tswitch1.dat<-data.frame(times=times, test.switch1=c(rep(1, start.time.test+1),
                                                     rep(0, length(times)-(start.time.test+1))))
tswitch2.dat<-data.frame(times=times, test.switch1=c(rep(0, start.time.test+1),
                                                     rep(1, length(times)-(start.time.test+1))))

cfxn<- approxfun(c.testdat$times, c.testdat$test.c, rule = 2)
afxn<- approxfun(a.testdat$times, a.testdat$test.a, rule = 2)
rcfxn<- approxfun(rc.testdat$times, rc.testdat$test.rc, rule = 2)
fcfxn<- approxfun(fc.testdat$times, fc.testdat$test.fc, rule = 2)
efxn<- approxfun(e.testdat$times, e.testdat$test.e, rule = 2)
sw1fxn<-approxfun(tswitch1.dat$times, tswitch1.dat$test.switch1, rule=2)
sw2fxn<-approxfun(tswitch2.dat$times, tswitch2.dat$test.switch2, rule=2)

#Example 1: do nothing
start.time.distance<-500
start.time.target<-500

initcase.c=10
initcase.a=100
initcase.e=10
initcase.rc=5
initcase.fc=5

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
#t = (0:364)


par = c(agestruc, 
        HomeContacts_5x5, WorkContacts_5x5, SchoolContacts_5x5, OtherContacts_5x5, 
        preduced, pfull, sd.other, alpha,
        start.time.distance, start.time.target, 
        gamma_e, gamma_a, gamma_s, gamma_h, p, 
        hosp_frac, hosp_crit, crit_die,
        sensitivity, specificity) 

model_out = ode(y = start, times = t, fun = seir_model_shields_rcfc, parms = par, 
                rtol=1e-8, hmax=0.1, method='ode45')


model_out <- as.data.frame(model_out)

null.model<-ddply(model_out, .(time), summarize, 
                  CriticalCare=sum(Hc.fc, Hc.rc, Hc.c, Hc.a, Hc.e, Hc.c),
                  Deaths=sum(D.c, D.a, D.e, D.rc, D.fc))
ggplot(data=null.model, aes(x=time, y=CriticalCare))+geom_line()+theme_classic()+
  labs(x='time', y='Critical Care Cases')

##Example 2: Untargeted social distancing#####
start.time.distance<-75
start.time.target<-500
par = c(agestruc, 
        HomeContacts_5x5, WorkContacts_5x5, SchoolContacts_5x5, OtherContacts_5x5, 
        preduced, pfull, sd.other, alpha,
        start.time.distance, start.time.target, 
        gamma_e, gamma_a, gamma_s, gamma_h, p, 
        hosp_frac, hosp_crit, crit_die,
        sensitivity, specificity) 

model_out = as.data.frame(ode(y = start, times = t, fun = seir_model_shields_rcfc, parms = par, 
                              rtol=1e-8, hmax=0.1, method='ode45'))

sd.model<-ddply(model_out, .(time), summarize, 
                CriticalCare=sum(Hc.fc, Hc.rc, Hc.c, Hc.a, Hc.e, Hc.c),
                Deaths=sum(D.c, D.a, D.e, D.rc, D.fc))

##Example 3: Targeted social distancing##########
start.time.distance<-75
start.time.target<-120
par = c(agestruc, 
        HomeContacts_5x5, WorkContacts_5x5, SchoolContacts_5x5, OtherContacts_5x5, 
        preduced, pfull, sd.other, alpha,
        start.time.distance, start.time.target, 
        gamma_e, gamma_a, gamma_s, gamma_h, p, 
        hosp_frac, hosp_crit, crit_die,
        sensitivity, specificity) 

model_out = as.data.frame(ode(y = start, times = t, fun = seir_model_shields_rcfc, parms = par, 
                              rtol=1e-8, hmax=0.1, method='ode45'))

targetsd.model<-ddply(model_out, .(time), summarize, 
                      CriticalCare=sum(Hc.fc, Hc.rc, Hc.c, Hc.a, Hc.e, Hc.c),
                      Deaths=sum(D.c, D.a, D.e, D.rc, D.fc))

ggplot()+#geom_line(data=null.model, aes(x=time, y=CriticalCare))+
  geom_line(data=sd.model, aes(x=time, y=CriticalCare))+
  geom_line(data=targetsd.model, aes(x=time, y=CriticalCare))+
  theme_classic()+
  labs(x='time', y='Critical Care Cases')

#These lines are on top of each other but start to diverge when targeting begins earlier
