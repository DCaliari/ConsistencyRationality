clear
clc

structural_estimation_3rdSubmission_CRRACARAEXPBETA  % estimates under the different specifications

load warpD.mat
load warpL.mat
WARP_D = WARPxyzw;
WARP_L = WARPabcd;
Raven = xlsread('DATI.xlsx','Variables','S1:S146');


deliberate_time =  xlsread('DATI.xlsx','Variables','FN1:FN146');
deliberate = xlsread('DATI.xlsx','Variables','EE1:EE146');

WARP_T = WARPabcd + WARPxyzw;

%% Figure 10 - heterogeneity in correlation WARP violations and Raven scores %%

id_heu = ti_exp' | ra_crra' | deliberate_time ==1 |  deliberate ==1 | rl_crra' | tp_exp';
id_rum = ~ti_exp' & ~tp_exp' & ~ra_crra' & ~rl_crra' & deliberate ==0 & deliberate_time ==0;

n1 = sum(id_heu);
n2 = sum(id_rum);


% Define a tiled layout with 2 rows: one for the main panel and one for the smaller panels
tiledlayout(2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');  % 2 rows, 2 columns

% Main panel: Aggregate (Top row spanning both columns)
nexttile([1 2]) % Span the main panel across both columns in the top row
scatter(WARP_T(id_heu), Raven(id_heu));  % Heuristic group scatter
hold on;
scatter(WARP_T(id_rum), Raven(id_rum));  % Rumination group scatter

% Regression lines for heuristic and rumination groups
% Heuristic group regression line
p_heu = polyfit(WARP_T(id_heu), Raven(id_heu), 1);
yfit_heu = polyval(p_heu, WARP_T(id_heu));
plot(WARP_T(id_heu), yfit_heu, 'b-', 'LineWidth', 1.5); % Blue line for heuristic

% Rumination group regression line
p_rum = polyfit(WARP_T(id_rum), Raven(id_rum), 1);
yfit_rum = polyval(p_rum, WARP_T(id_rum));
plot(WARP_T(id_rum), yfit_rum, 'r-', 'LineWidth', 1.5); % Red line for rumination

% Correlation coefficients for the main panel
[r1, p1] = corr(WARP_T(id_heu), Raven(id_heu));
[r2, p2] = corr(WARP_T(id_rum), Raven(id_rum));

% Display correlation coefficients on main panel
text(1, max(Raven)*0.9, ['Ext & Rand: r = ', num2str(r1, '%.2f'), ', p = ', num2str(p1, '%.3f')], 'Color', 'blue', 'FontName', 'Times New Roman');
text(1, max(Raven)*0.8, ['Others: r = ', num2str(r2, '%.2f'), ', p = ', num2str(p2, '%.3f')], 'Color', 'red', 'FontName', 'Times New Roman');

title('Aggregate', 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel(['WARP violations (n1 = ', num2str(n1), ', n2 = ', num2str(n2), ')'], 'FontName', 'Times New Roman');
ylabel('Raven', 'FontName', 'Times New Roman');
hold off;


id_heu = ti_exp' | deliberate_time ==1 | tp_exp';
id_rum = ~ti_exp' & ~tp_exp' & deliberate_time ==0;

n1 = sum(id_heu);
n2 = sum(id_rum);


% Smaller panels for "Time" and "Risk"
% "Time" panel (WARP_D)
nexttile
scatter(WARP_D(id_heu), Raven(id_heu));
hold on;
scatter(WARP_D(id_rum), Raven(id_rum));

% Regression lines for "Time" panel
% Heuristic group regression line
p_heu_time = polyfit(WARP_D(id_heu), Raven(id_heu), 1);
yfit_heu_time = polyval(p_heu_time, WARP_D(id_heu));
plot(WARP_D(id_heu), yfit_heu_time, 'b-', 'LineWidth', 1.5);

% Rumination group regression line
p_rum_time = polyfit(WARP_D(id_rum), Raven(id_rum), 1);
yfit_rum_time = polyval(p_rum_time, WARP_D(id_rum));
plot(WARP_D(id_rum), yfit_rum_time, 'r-', 'LineWidth', 1.5);

% Correlation coefficients for "Time" panel
[r1_time, p1_time] = corr(WARP_D(id_heu), Raven(id_heu));
[r2_time, p2_time] = corr(WARP_D(id_rum), Raven(id_rum));

% Display correlation coefficients on "Time" panel
text(1, max(Raven)*0.95, ['Ext & Rand: r = ', num2str(r1_time, '%.2f'), ', p = ', num2str(p1_time, '%.3f')], 'Color', 'blue', 'FontName', 'Times New Roman');
text(1, max(Raven)*0.85, ['Others: r = ', num2str(r2_time, '%.2f'), ', p = ', num2str(p2_time, '%.3f')], 'Color', 'red', 'FontName', 'Times New Roman');

title('Time', 'FontName', 'Times New Roman');
xlabel(['WARP violations (n1 = ', num2str(n1), ', n2 = ', num2str(n2), ')'], 'FontName', 'Times New Roman');
ylabel('Raven', 'FontName', 'Times New Roman');
hold off;


id_heu =  ra_crra' |  deliberate ==1 | rl_crra';
id_rum =  ~ra_crra' & ~rl_crra' & deliberate ==0;

n1 = sum(id_heu);
n2 = sum(id_rum);




% "Risk" panel (WARP_L)
nexttile
scatter(WARP_L(id_heu), Raven(id_heu));
hold on;
scatter(WARP_L(id_rum), Raven(id_rum));

% Regression lines for "Risk" panel
% Heuristic group regression line
p_heu_risk = polyfit(WARP_L(id_heu), Raven(id_heu), 1);
yfit_heu_risk = polyval(p_heu_risk, WARP_L(id_heu));
plot(WARP_L(id_heu), yfit_heu_risk, 'b-', 'LineWidth', 1.5);

% Rumination group regression line
p_rum_risk = polyfit(WARP_L(id_rum), Raven(id_rum), 1);
yfit_rum_risk = polyval(p_rum_risk, WARP_L(id_rum));
plot(WARP_L(id_rum), yfit_rum_risk, 'r-', 'LineWidth', 1.5);

% Correlation coefficients for "Risk" panel
[r1_risk, p1_risk] = corr(WARP_L(id_heu), Raven(id_heu));
[r2_risk, p2_risk] = corr(WARP_L(id_rum), Raven(id_rum));

% Display correlation coefficients on "Risk" panel
text(1, max(Raven)*0.95, ['Ext & Rand: r = ', num2str(r1_risk, '%.2f'), ', p = ', num2str(p1_risk, '%.3f')], 'Color', 'blue', 'FontName', 'Times New Roman');
text(1, max(Raven)*0.85, ['Others: r = ', num2str(r2_risk, '%.2f'), ', p = ', num2str(p2_risk, '%.3f')], 'Color', 'red', 'FontName', 'Times New Roman');

title('Risk', 'FontName', 'Times New Roman');
xlabel(['WARP violations (n1 = ', num2str(n1), ', n2 = ', num2str(n2), ')'], 'FontName', 'Times New Roman');
ylabel('Raven', 'FontName', 'Times New Roman');
hold off;



