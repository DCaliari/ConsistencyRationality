function LL = Log_lik_EB(pi_par,F)

LL = - sum(log(F*pi_par),'all');

end