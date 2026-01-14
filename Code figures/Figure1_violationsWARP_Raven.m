clear 
clc

%% Figure 1 %%

load warpD.mat
load warpL.mat
Raven = xlsread('DATI.xlsx','Variables','S1:S146');


% Set default font to Times New Roman for the entire figure
set(gca,'FontName','Times New Roman')
subplot(1,2,1)
scatter(WARPxyzw, Raven)
xlim([0 12])
[r1, p1] = corr(WARPxyzw, Raven);  % Calculate correlation and p-value
text(1, max(Raven)*0.9, ['r = ', num2str(r1,'%.2f'), ', p = ', num2str(p1,'%.3f')], 'FontName', 'Times New Roman')  % Display correlation and p-value
xlabel('WARP violations', 'FontName', 'Times New Roman')
ylabel('Raven scores', 'FontName', 'Times New Roman')

subplot(1,2,2)
scatter(WARPabcd, Raven)
xlim([0 12])
[r2, p2] = corr(WARPabcd, Raven);  % Calculate correlation and p-value
text(1, max(Raven)*0.9, ['r = ', num2str(r2,'%.2f'), ', p = ', num2str(p2,'%.3f')], 'FontName', 'Times New Roman')  % Display correlation and p-value
xlabel('WARP violations', 'FontName', 'Times New Roman')
ylabel('Raven scores', 'FontName', 'Times New Roman')
