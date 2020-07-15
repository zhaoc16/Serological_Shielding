model_out = ode(y = start, times = t, fun = seir_model_shields_rcfc_nolatent, parms = par, 
                method='ode45')

model_out=as.data.frame(model_out)
Infected=N-sum(model_out$S.c[366], model_out$S.c.pos[366], model_out$S.a[366], model_out$S.a.pos[366], 
               model_out$S.rc[366], model_out$S.rc.pos[366], model_out$S.fc[366], model_out$S.fc.pos[366],
               model_out$S.e[366], model_out$S.e.pos[366])
Pinfect=Infected/N
Infected; Pinfect
Deaths=sum(model_out$D.a[366], model_out$D.c[366], model_out$D.e[366], model_out$D.rc[366], 
           model_out$D.fc[366])
Deaths

null.model=ddply(model_out, .(time), summarize, 
                 CriticalCare=sum(Hc.fc, Hc.rc, Hc.c, Hc.a, Hc.e, Hc.c),
                 Deaths=sum(D.c, D.a, D.e, D.rc, D.fc),
                 CI=(N-sum(S.c, S.c.pos, S.a, S.a.pos, S.rc, S.rc.pos, S.fc, S.fc.pos, S.e, S.e.pos))/N)
plot(null.model$CriticalCare)

null.model$Deaths[60]; null.model$Deaths[90]