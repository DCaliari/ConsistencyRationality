%% Understanding

Understanding = xlsread('DATI.xlsx','Variables','AT1:AU146');
Understanding = sum(Understanding,2);
Understanding = 10-Understanding;


[r1, p1] = corr(Raven, Understanding, 'Type', 'Pearson')


[r1, p1] = corr(WARP_D, Understanding, 'Type', 'Spearman')
[r1, p1] = corr(WARP_D, Understanding, 'Type', 'Pearson')

[r1, p1] = corr(WARP_L, Understanding, 'Type', 'Spearman')
[r1, p1] = corr(WARP_L, Understanding, 'Type', 'Pearson')

[r1, p1] = corr(Monotonicity, Understanding, 'Type', 'Pearson')
[r1, p1] = corr(Impatience, Understanding, 'Type', 'Pearson')
[r1, p1] = corr(FOSD, Understanding, 'Type', 'Pearson')
[r1, p1] = corr(SOSD, Understanding, 'Type', 'Pearson')



WARP_T = WARP_D + WARP_L;



%% ALL TOGETHER


id_heu = time_ind < 0.91 | risk_ind >2.2 | deliberate_time' ==1 |  deliberate' ==1 | risk_ind< 0.2 | time_ind>0.99;
id_rum = time_ind >= 0.91 & time_ind <= 0.99 & risk_ind <= 2.2 & risk_ind >=0.2 & deliberate' ==0 & deliberate_time' ==0;


[r1, p1] = corr(WARP_T(id_heu), Understanding(id_heu), 'Type', 'Pearson')
[r1, p1] = corr(WARP_T(id_rum), Understanding(id_rum), 'Type', 'Pearson')




id_heu = time_ind < 0.91 | deliberate_time' ==1 |  time_ind>0.99;
id_rum = time_ind >= 0.91 & time_ind <= 0.99 & deliberate_time' ==0;


[r1, p1] = corr(WARP_D(id_heu), Understanding(id_heu), 'Type', 'Pearson')
[r1, p1] = corr(WARP_D(id_rum), Understanding(id_rum), 'Type', 'Pearson')



id_heu =  risk_ind >2.2 |  deliberate' ==1 | risk_ind< 0.2;
id_rum =  risk_ind <= 2.2 & risk_ind >=0.2 & deliberate' ==0;

[r1, p1] = corr(WARP_L(id_heu), Understanding(id_heu), 'Type', 'Pearson')
[r1, p1] = corr(WARP_L(id_rum), Understanding(id_rum), 'Type', 'Pearson')



%% Not weighted raven

nraven = xlsread('DATI.xlsx','Variables','P1:P146');



[r1, p1] = corr(Raven, nraven, 'Type', 'Pearson')


[r1, p1] = corr(WARP_D, nraven, 'Type', 'Spearman')
[r1, p1] = corr(WARP_D, nraven, 'Type', 'Pearson')

[r1, p1] = corr(WARP_L, nraven, 'Type', 'Spearman')
[r1, p1] = corr(WARP_L, nraven, 'Type', 'Pearson')

[r1, p1] = corr(Monotonicity, nraven, 'Type', 'Pearson')
[r1, p1] = corr(Impatience, nraven, 'Type', 'Pearson')
[r1, p1] = corr(FOSD, nraven, 'Type', 'Pearson')
[r1, p1] = corr(SOSD, nraven, 'Type', 'Pearson')



WARP_T = WARP_D + WARP_L;



%% ALL TOGETHER


id_heu = time_ind < 0.91 | risk_ind >2.2 | deliberate_time' ==1 |  deliberate' ==1 | risk_ind< 0.2 | time_ind>0.99;
id_rum = time_ind >= 0.91 & time_ind <= 0.99 & risk_ind <= 2.2 & risk_ind >=0.2 & deliberate' ==0 & deliberate_time' ==0;


[r1, p1] = corr(WARP_T(id_heu), nraven(id_heu), 'Type', 'Pearson')
[r1, p1] = corr(WARP_T(id_rum), nraven(id_rum), 'Type', 'Pearson')




id_heu = time_ind < 0.91 | deliberate_time' ==1 |  time_ind>0.99;
id_rum = time_ind >= 0.91 & time_ind <= 0.99 & deliberate_time' ==0;


[r1, p1] = corr(WARP_D(id_heu), nraven(id_heu), 'Type', 'Pearson')
[r1, p1] = corr(WARP_D(id_rum), nraven(id_rum), 'Type', 'Pearson')



id_heu =  risk_ind >2.2 |  deliberate' ==1 | risk_ind< 0.2;
id_rum =  risk_ind <= 2.2 & risk_ind >=0.2 & deliberate' ==0;

[r1, p1] = corr(WARP_L(id_heu), nraven(id_heu), 'Type', 'Pearson')
[r1, p1] = corr(WARP_L(id_rum), nraven(id_rum), 'Type', 'Pearson')




%% CRT

CRT = xlsread('DATI.xlsx','Variables','AJ1:AJ146');



[r1, p1] = corr(Raven, CRT, 'Type', 'Pearson')


[r1, p1] = corr(WARP_D, CRT, 'Type', 'Spearman')
[r1, p1] = corr(WARP_D, CRT, 'Type', 'Pearson')

[r1, p1] = corr(WARP_L, CRT, 'Type', 'Spearman')
[r1, p1] = corr(WARP_L, CRT, 'Type', 'Pearson')

[r1, p1] = corr(Monotonicity, CRT, 'Type', 'Pearson')
[r1, p1] = corr(Impatience, CRT, 'Type', 'Pearson')
[r1, p1] = corr(FOSD, CRT, 'Type', 'Pearson')
[r1, p1] = corr(SOSD, CRT, 'Type', 'Pearson')



WARP_T = WARP_D + WARP_L;



%% ALL TOGETHER


id_heu = time_ind < 0.91 | risk_ind >2.2 | deliberate_time' ==1 |  deliberate' ==1 | risk_ind< 0.2 | time_ind>0.99;
id_rum = time_ind >= 0.91 & time_ind <= 0.99 & risk_ind <= 2.2 & risk_ind >=0.2 & deliberate' ==0 & deliberate_time' ==0;


[r1, p1] = corr(WARP_T(id_heu), CRT(id_heu), 'Type', 'Pearson')
[r1, p1] = corr(WARP_T(id_rum), CRT(id_rum), 'Type', 'Pearson')




id_heu = time_ind < 0.91 | deliberate_time' ==1 |  time_ind>0.99;
id_rum = time_ind >= 0.91 & time_ind <= 0.99 & deliberate_time' ==0;


[r1, p1] = corr(WARP_D(id_heu), CRT(id_heu), 'Type', 'Pearson')
[r1, p1] = corr(WARP_D(id_rum), CRT(id_rum), 'Type', 'Pearson')



id_heu =  risk_ind >2.2 |  deliberate' ==1 | risk_ind< 0.2;
id_rum =  risk_ind <= 2.2 & risk_ind >=0.2 & deliberate' ==0;

[r1, p1] = corr(WARP_L(id_heu), CRT(id_heu), 'Type', 'Pearson')
[r1, p1] = corr(WARP_L(id_rum), CRT(id_rum), 'Type', 'Pearson')


