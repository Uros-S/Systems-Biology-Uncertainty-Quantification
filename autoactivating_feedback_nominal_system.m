function dXdt = autoactivating_feedback_nominal_system(t,X,PAR,vars)

 x = X;

 dXdt = PAR.K + PAR.a*(x^PAR.n)/(1+x^PAR.n) - x;
end