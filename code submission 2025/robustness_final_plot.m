
WARP_T = WARP_D + WARP_L;


ra_all = [ra_crra', ra_cara', ra_crrace', ra_carace'];
rl_all = [rl_crra', rl_cara', rl_crrace', rl_carace'];


ti_all = [ti_exp', ti_hyp', ti_betadelta'];
tp_all = [tp_exp', tp_hyp', tp_betadelta'];


%% ALL TOGETHER

for z=1:3
    for j=1:4
        ti = ti_all(:,z);
        tp = tp_all(:,z);
        ra = ra_all(:,j);
        rl = rl_all(:,j);

id_heu = ti | rl | deliberate_time ==1 |  deliberate ==1 | ra | tp;
id_rum = ~ti & ~tp & ~rl & ~ra & deliberate ==0 & deliberate_time ==0;

% Correlation coefficients for the main panel
[r1(z,j), p1(z,j)] = corr(WARP_T(id_heu), Raven(id_heu));
[r2(z,j), p2(z,j)] = corr(WARP_T(id_rum), Raven(id_rum));

    end
end

for z=1:3
        ti = ti_all(:,z);
        tp = tp_all(:,z);

id_heu = ti | deliberate_time ==1 |  tp;
id_rum = ~ti & ~tp & deliberate_time ==0;

% Correlation coefficients for "Time" panel
[r1_time(z,1), p1_time(z,1)] = corr(WARP_D(id_heu), Raven(id_heu));
[r2_time(z,1), p2_time(z,1)] = corr(WARP_D(id_rum), Raven(id_rum));

end



for j=1:4
        ra = ra_all(:,j);
        rl = rl_all(:,j);

id_heu =  ra |  deliberate ==1 | rl;
id_rum =  ~ra & ~rl & deliberate ==0;

% Correlation coefficients for "Risk" panel
[r1_risk(j,1), p1_risk(j,1)] = corr(WARP_L(id_heu), Raven(id_heu));
[r2_risk(j,1), p2_risk(j,1)] = corr(WARP_L(id_rum), Raven(id_rum));

end
