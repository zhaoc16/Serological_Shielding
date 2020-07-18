source('reorg/input_SD_Simulations.R')

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

# null.model=ddply(model_out, .(time), summarize,
#                  CriticalCare=sum(Hc.fc, Hc.rc, Hc.c, Hc.a, Hc.e, Hc.c),
#                  Deaths=sum(D.c, D.a, D.e, D.rc, D.fc),
#                  CI=(N-sum(S.c, S.c.pos, S.a, S.a.pos, S.rc, S.rc.pos, S.fc, S.fc.pos, S.e, S.e.pos))/N)
# plot(null.model$CriticalCare)
# 
# null.model$Deaths[60]; null.model$Deaths[90]

# model_out$S.c[366], model_out$S.c.pos[366],
# model_out$S.a[366], model_out$S.a.pos[366], 
# model_out$S.rc[366], model_out$S.rc.pos[366],
# model_out$S.fc[366], model_out$S.fc.pos[366],
# model_out$S.e[366], model_out$S.e.pos[366]