%% replication with all the methodologies

%% loading all the methodologies


XYZW_OW = xlsread('DATI.xlsx','Variables','CH1:CH146');
WZYX_OW = xlsread('DATI.xlsx','Variables','CG1:CG146');

ABCD_OW = xlsread('DATI.xlsx','Variables','CE1:CE146');
DCBA_OW = xlsread('DATI.xlsx','Variables','CD1:CD146');


XYZW_SEQ = xlsread('DATI.xlsx','Variables','EB1:EB146');
WZYX_SEQ = xlsread('DATI.xlsx','Variables','EA1:EA146');

ABCD_SEQ = xlsread('DATI.xlsx','Variables','ED1:ED146');
DCBA_SEQ = xlsread('DATI.xlsx','Variables','EC1:EC146');


XYZW_rep = xlsread('DATI.xlsx','Variables','AA1:AA146');
WZYX_rep = xlsread('DATI.xlsx','Variables','Z1:Z146');

ABCD_rep = xlsread('DATI.xlsx','Variables','X1:X146');
DCBA_rep = xlsread('DATI.xlsx','Variables','W1:W146');


XYZW_MS = xlsread('DATI.xlsx','Variables','EI1:EI146');
WZYX_MS = xlsread('DATI.xlsx','Variables','EH1:EH146');

ABCD_MS = xlsread('DATI.xlsx','Variables','EG1:EG146');
DCBA_MS = xlsread('DATI.xlsx','Variables','EF1:EF146');


XYZW_que = xlsread('DATI.xlsx','Variables','DJ1:DJ146');
WZYX_que = xlsread('DATI.xlsx','Variables','DI1:DI146');

ABCD_que = xlsread('DATI.xlsx','Variables','FM1:FM146');
DCBA_que = xlsread('DATI.xlsx','Variables','FL1:FL146');

heu_time = [ti_exp' | tp_exp', XYZW_OW | WZYX_OW, XYZW_MS | WZYX_MS, XYZW_SEQ | WZYX_SEQ, XYZW_rep | WZYX_rep, XYZW_que | WZYX_que];
heu_risk = [ra_crra' | rl_crra', ABCD_OW | DCBA_OW, ABCD_MS | DCBA_MS, ABCD_SEQ | DCBA_SEQ, ABCD_rep | DCBA_rep, ABCD_que | DCBA_que];

cat = categorical(["Structural", "Optimal Weighting", "Minimum Swaps", "Sequential", "Reported", "Questionnaire"]);
cat = reordercats(cat,["Structural", "Optimal Weighting", "Minimum Swaps", "Sequential", "Reported", "Questionnaire"]);



%% WARP


WARP_D = xlsread('DATI.xlsx','Variables','C1:C146');
WARP_L = xlsread('DATI.xlsx','Variables','E1:E146');



for i=1:size(heu_time,2)
mean_WARP(1,i) = mean(WARP_D(heu_time(:,i))); 
mean_WARP(2,i)= mean(WARP_D(~heu_time(:,i)));
end


% Assuming mean_WARP is already computed as provided in the previous code

% Perform t-tests for each comparison and store p-values
p_values = zeros(1, 6); % Preallocate p-values array for required comparisons

% Conduct t-tests for specific pairs between decisions with and without heuristics
% Replace these comparisons with the relevant WARP_D data as needed
[~, p_values(1)] = ttest2(WARP_D(heu_time(:, 1)), WARP_D(~heu_time(:, 1)),'Tail','left');
[~, p_values(2)] = ttest2(WARP_D(heu_time(:, 2)), WARP_D(~heu_time(:, 2)),'Tail','left'); 
[~, p_values(3)] = ttest2(WARP_D(heu_time(:, 3)), WARP_D(~heu_time(:, 3)),'Tail','left');
[~, p_values(4)] = ttest2(WARP_D(heu_time(:, 4)), WARP_D(~heu_time(:, 4)),'Tail','left');
[~, p_values(5)] = ttest2(WARP_D(heu_time(:, 5)), WARP_D(~heu_time(:, 5)),'Tail','left');
[~, p_values(6)] = ttest2(WARP_D(heu_time(:, 6)), WARP_D(~heu_time(:, 6)),'Tail','left');

% Plot the bar graph
figure;
bar(cat, mean_WARP);

% Set the x-tick labels' font to Times New Roman
set(gca, 'FontName', 'Times New Roman');

% Set the y-axis label with the desired font settings
ylabel('WARP violations', 'FontName', 'Times New Roman');

% Set the limits for the y-axis; adjust as per the expected range of your data
ylim([0 8]);

% Define x positions for the bars
x_positions = 1:length(mean_WARP);

% Define y positions for each comparison (adjust these as needed)
y_bracket_12 = 6.5;   % Position for comparison Bin 1 vs Bin 2
y_bracket_13 = 6.5;   % Position for comparison Bin 1 vs Bin 3
y_bracket_14 = 6.5;   % Position for comparison Bin 1 vs Bin 4
y_bracket_23 = 6.5;   % Position for comparison Bin 2 vs Bin 3
y_bracket_24 = 6.5;   % Position for comparison Bin 2 vs Bin 4
y_bracket_34 = 6.5;   % Position for comparison Bin 3 vs Bin 4


hold on;

% Add comparisons and p-values for each

% First comparison (Bin 1 vs Bin 2)
plot([x_positions(1)-0.15, x_positions(1)-.15, x_positions(1)+.15, x_positions(1)+.15], ...
     [y_bracket_12, y_bracket_12 + 0.25, y_bracket_12 + 0.25, y_bracket_12], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(1)]), y_bracket_12 + 0.5, sprintf('p = %.3f', p_values(1)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Second comparison (Bin 1 vs Bin 3)
plot([x_positions(2)-0.15, x_positions(2)-.15, x_positions(2)+.15, x_positions(2)+.15], ...
     [y_bracket_13, y_bracket_13 + 0.25, y_bracket_13 + 0.25, y_bracket_13], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(2)]), y_bracket_13 + 0.5, sprintf('p = %.3f', p_values(2)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Third comparison (Bin 1 vs Bin 4)
plot([x_positions(3)-0.15, x_positions(3)-.15, x_positions(3)+.15, x_positions(3)+.15], ...
     [y_bracket_14, y_bracket_14 + 0.25, y_bracket_14 + 0.25, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(3), x_positions(3)]), y_bracket_14 + 0.5, sprintf('p = %.3f', p_values(3)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Fourth comparison (Bin 2 vs Bin 3)
plot([x_positions(4)-0.15, x_positions(4)-.15, x_positions(4)+.15, x_positions(4)+.15], ...
     [y_bracket_23, y_bracket_23 + 0.25, y_bracket_23 + 0.25, y_bracket_23], 'k', 'LineWidth', 1.5);
text(mean([x_positions(4), x_positions(4)]), y_bracket_23 + 0.5, sprintf('p = %.3f', p_values(4)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Fifth comparison (Bin 2 vs Bin 4)
plot([x_positions(5)-0.15, x_positions(5)-.15, x_positions(5)+.15, x_positions(5)+.15], ...
     [y_bracket_24, y_bracket_24 + 0.25, y_bracket_24 + 0.25, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(5), x_positions(5)]), y_bracket_24 + 0.5, sprintf('p = %.3f', p_values(5)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Sixth comparison (Bin 3 vs Bin 4)
plot([x_positions(6)-0.15, x_positions(6)-.15, x_positions(6)+.15, x_positions(6)+.15], ...
     [y_bracket_34, y_bracket_34 + 0.25, y_bracket_34 + 0.25, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean([x_positions(6), x_positions(6)]), y_bracket_34 + 0.5, sprintf('p = %.3f', p_values(6)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');




for i=1:size(heu_risk,2)
mean_WARP(1,i) = mean(WARP_L(heu_risk(:,i))); 
mean_WARP(2,i)= mean(WARP_L(~heu_risk(:,i)));
end


% Assuming mean_WARP is already computed as provided in the previous code

% Perform t-tests for each comparison and store p-values
p_values = zeros(1, 6); % Preallocate p-values array for required comparisons

% Conduct t-tests for specific pairs between decisions with and without heuristics
% Replace these comparisons with the relevant WARP_L data as needed
[~, p_values(1)] = ttest2(WARP_L(heu_risk(:, 1)), WARP_L(~heu_risk(:, 1)),'Tail','left');
[~, p_values(2)] = ttest2(WARP_L(heu_risk(:, 2)), WARP_L(~heu_risk(:, 2)),'Tail','left');
[~, p_values(3)] = ttest2(WARP_L(heu_risk(:, 3)), WARP_L(~heu_risk(:, 3)),'Tail','left');
[~, p_values(4)] = ttest2(WARP_L(heu_risk(:, 4)), WARP_L(~heu_risk(:, 4)),'Tail','left');
[~, p_values(5)] = ttest2(WARP_L(heu_risk(:, 5)), WARP_L(~heu_risk(:, 5)),'Tail','left');
[~, p_values(6)] = ttest2(WARP_L(heu_risk(:, 6)), WARP_L(~heu_risk(:, 6)),'Tail','left');

% Plot the bar graph
figure;
bar(cat, mean_WARP);

% Set the x-tick labels' font to Times New Roman
set(gca, 'FontName', 'Times New Roman');

% Set the y-axis label with the desired font settings
ylabel('WARP violations', 'FontName', 'Times New Roman');

% Set the limits for the y-axis; adjust as per the expected range of your data
ylim([0 8]);

% Define x positions for the bars
x_positions = 1:length(mean_WARP);

% Define y positions for each comparison (adjust these as needed)
y_bracket_12 = 6.5;   % Position for comparison Bin 1 vs Bin 2
y_bracket_13 = 6.5;   % Position for comparison Bin 1 vs Bin 3
y_bracket_14 = 6.5;   % Position for comparison Bin 1 vs Bin 4
y_bracket_23 = 6.5;   % Position for comparison Bin 2 vs Bin 3
y_bracket_24 = 6.5;   % Position for comparison Bin 2 vs Bin 4
y_bracket_34 = 6.5;   % Position for comparison Bin 3 vs Bin 4

hold on;

% Add comparisons and p-values for each

% First comparison (Bin 1 vs Bin 2)
plot([x_positions(1)-0.15, x_positions(1)-.15, x_positions(1)+.15, x_positions(1)+.15], ...
     [y_bracket_12, y_bracket_12 + 0.25, y_bracket_12 + 0.25, y_bracket_12], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(1)]), y_bracket_12 + 0.5, sprintf('p = %.3f', p_values(1)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Second comparison (Bin 1 vs Bin 3)
plot([x_positions(2)-0.15, x_positions(2)-.15, x_positions(2)+.15, x_positions(2)+.15], ...
     [y_bracket_13, y_bracket_13 + 0.25, y_bracket_13 + 0.25, y_bracket_13], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(2)]), y_bracket_13 + 0.5, sprintf('p = %.3f', p_values(2)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Third comparison (Bin 1 vs Bin 4)
plot([x_positions(3)-0.15, x_positions(3)-.15, x_positions(3)+.15, x_positions(3)+.15], ...
     [y_bracket_14, y_bracket_14 + 0.25, y_bracket_14 + 0.25, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(3), x_positions(3)]), y_bracket_14 + 0.5, sprintf('p = %.3f', p_values(3)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Fourth comparison (Bin 2 vs Bin 3)
plot([x_positions(4)-0.15, x_positions(4)-.15, x_positions(4)+.15, x_positions(4)+.15], ...
     [y_bracket_23, y_bracket_23 + 0.25, y_bracket_23 + 0.25, y_bracket_23], 'k', 'LineWidth', 1.5);
text(mean([x_positions(4), x_positions(4)]), y_bracket_23 + 0.5, sprintf('p = %.3f', p_values(4)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Fifth comparison (Bin 2 vs Bin 4)
plot([x_positions(5)-0.15, x_positions(5)-.15, x_positions(5)+.15, x_positions(5)+.15], ...
     [y_bracket_24, y_bracket_24 + 0.25, y_bracket_24 + 0.25, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(5), x_positions(5)]), y_bracket_24 + 0.5, sprintf('p = %.3f', p_values(5)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Sixth comparison (Bin 3 vs Bin 4)
plot([x_positions(6)-0.15, x_positions(6)-.15, x_positions(6)+.15, x_positions(6)+.15], ...
     [y_bracket_34, y_bracket_34 + 0.25, y_bracket_34 + 0.25, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean([x_positions(6), x_positions(6)]), y_bracket_34 + 0.5, sprintf('p = %.3f', p_values(6)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

hold off;




%% Raven

Raven = xlsread('DATI.xlsx','Variables','S1:S146');



heu_time1 = [XYZW_str , XYZW_OW , XYZW_MS , XYZW_SEQ , XYZW_rep, XYZW_que];
heu_time2 = [WZYX_str, WZYX_OW, WZYX_MS, WZYX_SEQ,  WZYX_rep,  WZYX_que];

heu_time1 = logical(heu_time1);
heu_time2 = logical(heu_time2);

heu_risk1 = [ABCD_str , ABCD_OW , ABCD_MS, ABCD_SEQ, ABCD_rep , ABCD_que];
heu_risk2 = [DCBA_str,  DCBA_OW,  DCBA_MS, DCBA_SEQ,  DCBA_rep, DCBA_que];

heu_risk1 = logical(heu_risk1);
heu_risk2 = logical(heu_risk2);


for i=1:size(heu_time,2)
mean_raven(1,i) = mean(Raven(heu_time1(:,i))); 
mean_raven(2,i)= mean(Raven(heu_time2(:,i)));
end

% Assuming mean_WARP is already computed as provided in the previous code

% Perform t-tests for each comparison and store p-values
p_values = zeros(1, 6); % Preallocate p-values array for required comparisons

% Conduct t-tests for specific pairs between decisions with and without heuristics
% Replace these comparisons with the relevant Raven data as needed
[~, p_values(1)] = ttest2(Raven(heu_time1(:, 1)), Raven(heu_time2(:, 1)),'Tail','left'); 
[~, p_values(2)] = ttest2(Raven(heu_time1(:, 2)), Raven(heu_time2(:, 2)),'Tail','left');
[~, p_values(3)] = ttest2(Raven(heu_time1(:, 3)), Raven(heu_time2(:, 3)),'Tail','left');
[~, p_values(4)] = ttest2(Raven(heu_time1(:, 4)), Raven(heu_time2(:, 4)),'Tail','left');
[~, p_values(5)] = ttest2(Raven(heu_time1(:, 5)), Raven(heu_time2(:, 5)),'Tail','left');
[~, p_values(6)] = ttest2(Raven(heu_time1(:, 6)), Raven(heu_time2(:, 6)),'Tail','left');

% Plot the bar graph
figure;
bar(cat, mean_raven);

% Set the x-tick labels' font to Times New Roman
set(gca, 'FontName', 'Times New Roman');

% Set the y-axis label with the desired font settings
ylabel('Raven scores', 'FontName', 'Times New Roman');

% Set the limits for the y-axis; adjust as per the expected range of your data
ylim([5 15]);

% Define x positions for the bars
x_positions = 1:length(mean_WARP);

% Define y positions for each comparison (adjust these as needed)
y_bracket_12 = 11.5;   % Position for comparison Bin 1 vs Bin 2
y_bracket_13 = 11.5;   % Position for comparison Bin 1 vs Bin 3
y_bracket_14 = 11.5;   % Position for comparison Bin 1 vs Bin 4
y_bracket_23 = 11.5;   % Position for comparison Bin 2 vs Bin 3
y_bracket_24 = 11.5;   % Position for comparison Bin 2 vs Bin 4
y_bracket_34 = 11.5;   % Position for comparison Bin 3 vs Bin 4


hold on;

% Add comparisons and p-values for each

% First comparison (Bin 1 vs Bin 2)
plot([x_positions(1)-0.15, x_positions(1)-.15, x_positions(1)+.15, x_positions(1)+.15], ...
     [y_bracket_12, y_bracket_12 + 0.25, y_bracket_12 + 0.25, y_bracket_12], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(1)]), y_bracket_12 + 0.5, sprintf('p = %.3f', p_values(1)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Second comparison (Bin 1 vs Bin 3)
plot([x_positions(2)-0.15, x_positions(2)-.15, x_positions(2)+.15, x_positions(2)+.15], ...
     [y_bracket_13, y_bracket_13 + 0.25, y_bracket_13 + 0.25, y_bracket_13], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(2)]), y_bracket_13 + 0.5, sprintf('p = %.3f', p_values(2)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Third comparison (Bin 1 vs Bin 4)
plot([x_positions(3)-0.15, x_positions(3)-.15, x_positions(3)+.15, x_positions(3)+.15], ...
     [y_bracket_14, y_bracket_14 + 0.25, y_bracket_14 + 0.25, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(3), x_positions(3)]), y_bracket_14 + 0.5, sprintf('p = %.3f', p_values(3)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Fourth comparison (Bin 2 vs Bin 3)
plot([x_positions(4)-0.15, x_positions(4)-.15, x_positions(4)+.15, x_positions(4)+.15], ...
     [y_bracket_23, y_bracket_23 + 0.25, y_bracket_23 + 0.25, y_bracket_23], 'k', 'LineWidth', 1.5);
text(mean([x_positions(4), x_positions(4)]), y_bracket_23 + 0.5, sprintf('p = %.3f', p_values(4)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Fifth comparison (Bin 2 vs Bin 4)
plot([x_positions(5)-0.15, x_positions(5)-.15, x_positions(5)+.15, x_positions(5)+.15], ...
     [y_bracket_24, y_bracket_24 + 0.25, y_bracket_24 + 0.25, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(5), x_positions(5)]), y_bracket_24 + 0.5, sprintf('p = %.3f', p_values(5)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Sixth comparison (Bin 3 vs Bin 4)
plot([x_positions(6)-0.15, x_positions(6)-.15, x_positions(6)+.15, x_positions(6)+.15], ...
     [y_bracket_34, y_bracket_34 + 0.25, y_bracket_34 + 0.25, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean([x_positions(6), x_positions(6)]), y_bracket_34 + 0.5, sprintf('p = %.3f', p_values(6)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
legend('Impatient', 'Patient')





for i=1:size(heu_time,2)
mean_raven(1,i) = mean(Raven(heu_risk1(:,i))); 
mean_raven(2,i)= mean(Raven(heu_risk2(:,i)));
end

% Assuming mean_WARP is already computed as provided in the previous code

% Perform t-tests for each comparison and store p-values
p_values = zeros(1, 6); % Preallocate p-values array for required comparisons

% Conduct t-tests for specific pairs between decisions with and without heuristics
% Replace these comparisons with the relevant Raven data as needed
[~, p_values(1)] = ttest2(Raven(heu_risk1(:, 1)), Raven(heu_risk2(:, 1)),'Tail','left');
[~, p_values(2)] = ttest2(Raven(heu_risk1(:, 2)), Raven(heu_risk2(:, 2)),'Tail','left');
[~, p_values(3)] = ttest2(Raven(heu_risk1(:, 3)), Raven(heu_risk2(:, 3)),'Tail','left');
[~, p_values(4)] = ttest2(Raven(heu_risk1(:, 4)), Raven(heu_risk2(:, 4)),'Tail','left');
[~, p_values(5)] = ttest2(Raven(heu_risk1(:, 5)), Raven(heu_risk2(:, 5)),'Tail','left');
[~, p_values(6)] = ttest2(Raven(heu_risk1(:, 6)), Raven(heu_risk2(:, 6)),'Tail','left');

% Plot the bar graph
figure;
bar(cat, mean_raven);

% Set the x-tick labels' font to Times New Roman
set(gca, 'FontName', 'Times New Roman');

% Set the y-axis label with the desired font settings
ylabel('Raven scores', 'FontName', 'Times New Roman');

% Set the limits for the y-axis; adjust as per the expected range of your data
ylim([5 15]);

% Define x positions for the bars
x_positions = 1:length(mean_WARP);

% Define y positions for each comparison (adjust these as needed)
y_bracket_12 = 12.5;   % Position for comparison Bin 1 vs Bin 2
y_bracket_13 = 12.5;   % Position for comparison Bin 1 vs Bin 3
y_bracket_14 = 12.5;   % Position for comparison Bin 1 vs Bin 4
y_bracket_23 = 12.5;   % Position for comparison Bin 2 vs Bin 3
y_bracket_24 = 12.5;   % Position for comparison Bin 2 vs Bin 4
y_bracket_34 = 12.5;   % Position for comparison Bin 3 vs Bin 4


hold on;

% Add comparisons and p-values for each

% First comparison (Bin 1 vs Bin 2)
plot([x_positions(1)-0.15, x_positions(1)-.15, x_positions(1)+.15, x_positions(1)+.15], ...
     [y_bracket_12, y_bracket_12 + 0.25, y_bracket_12 + 0.25, y_bracket_12], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(1)]), y_bracket_12 + 0.5, sprintf('p = %.3f', p_values(1)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Second comparison (Bin 1 vs Bin 3)
plot([x_positions(2)-0.15, x_positions(2)-.15, x_positions(2)+.15, x_positions(2)+.15], ...
     [y_bracket_13, y_bracket_13 + 0.25, y_bracket_13 + 0.25, y_bracket_13], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(2)]), y_bracket_13 + 0.5, sprintf('p = %.3f', p_values(2)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Third comparison (Bin 1 vs Bin 4)
plot([x_positions(3)-0.15, x_positions(3)-.15, x_positions(3)+.15, x_positions(3)+.15], ...
     [y_bracket_14, y_bracket_14 + 0.25, y_bracket_14 + 0.25, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(3), x_positions(3)]), y_bracket_14 + 0.5, sprintf('p = %.3f', p_values(3)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Fourth comparison (Bin 2 vs Bin 3)
plot([x_positions(4)-0.15, x_positions(4)-.15, x_positions(4)+.15, x_positions(4)+.15], ...
     [y_bracket_23, y_bracket_23 + 0.25, y_bracket_23 + 0.25, y_bracket_23], 'k', 'LineWidth', 1.5);
text(mean([x_positions(4), x_positions(4)]), y_bracket_23 + 0.5, sprintf('p = %.3f', p_values(4)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Fifth comparison (Bin 2 vs Bin 4)
plot([x_positions(5)-0.15, x_positions(5)-.15, x_positions(5)+.15, x_positions(5)+.15], ...
     [y_bracket_24, y_bracket_24 + 0.25, y_bracket_24 + 0.25, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(5), x_positions(5)]), y_bracket_24 + 0.5, sprintf('p = %.3f', p_values(5)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Sixth comparison (Bin 3 vs Bin 4)
plot([x_positions(6)-0.15, x_positions(6)-.15, x_positions(6)+.15, x_positions(6)+.15], ...
     [y_bracket_34, y_bracket_34 + 0.25, y_bracket_34 + 0.25, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean([x_positions(6), x_positions(6)]), y_bracket_34 + 0.5, sprintf('p = %.3f', p_values(6)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
legend('Risk-averse', 'Risk-neutral')





%% Response Times

Times = xlsread('DATI.xlsx','Variables','BN1:BO146');
Response_D = Times(:,1);
Response_L = Times(:,2);


for i=1:size(heu_time,2)
mean_resp(1,i) = mean(Response_D(heu_time1(:,i))); 
mean_resp(2,i)= mean(Response_D(heu_time2(:,i)));
end



% Perform t-tests for each comparison and store p-values
p_values = zeros(1, 6); % Preallocate p-values array for required comparisons

% Conduct t-tests for specific pairs between decisions with and without heuristics
% Replace these comparisons with the relevant Raven data as needed
[~, p_values(1)] = ttest2(Response_D(heu_time1(:, 1)), Response_D(heu_time2(:, 1)),'Tail','left');
[~, p_values(2)] = ttest2(Response_D(heu_time1(:, 2)), Response_D(heu_time2(:, 2)),'Tail','left');
[~, p_values(3)] = ttest2(Response_D(heu_time1(:, 3)), Response_D(heu_time2(:, 3)),'Tail','left'); 
[~, p_values(4)] = ttest2(Response_D(heu_time1(:, 4)), Response_D(heu_time2(:, 4)),'Tail','left');
[~, p_values(5)] = ttest2(Response_D(heu_time1(:, 5)), Response_D(heu_time2(:, 5)),'Tail','left');
[~, p_values(6)] = ttest2(Response_D(heu_time1(:, 6)), Response_D(heu_time2(:, 6)),'Tail','left');

% Plot the bar graph
figure;
bar(cat, mean_resp);

% Set the x-tick labels' font to Times New Roman
set(gca, 'FontName', 'Times New Roman');

% Set the y-axis label with the desired font settings
ylabel('Response times', 'FontName', 'Times New Roman');

% Set the limits for the y-axis; adjust as per the expected range of your data
ylim([0 650]);

% Define x positions for the bars
x_positions = 1:length(mean_WARP);

% Define y positions for each comparison (adjust these as needed)
y_bracket_12 = 500;   % Position for comparison Bin 1 vs Bin 2
y_bracket_13 = 500;   % Position for comparison Bin 1 vs Bin 3
y_bracket_14 = 500;   % Position for comparison Bin 1 vs Bin 4
y_bracket_23 = 500;   % Position for comparison Bin 2 vs Bin 3
y_bracket_24 = 500;   % Position for comparison Bin 2 vs Bin 4
y_bracket_34 = 500;   % Position for comparison Bin 3 vs Bin 4


hold on;

% Add comparisons and p-values for each

% First comparison (Bin 1 vs Bin 2)
plot([x_positions(1)-0.15, x_positions(1)-.15, x_positions(1)+.15, x_positions(1)+.15], ...
     [y_bracket_12, y_bracket_12 + 10, y_bracket_12 + 10, y_bracket_12], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(1)]), y_bracket_12 + 35, sprintf('p = %.3f', p_values(1)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Second comparison (Bin 1 vs Bin 3)
plot([x_positions(2)-0.15, x_positions(2)-.15, x_positions(2)+.15, x_positions(2)+.15], ...
     [y_bracket_13, y_bracket_13 + 10, y_bracket_13 + 10, y_bracket_13], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(2)]), y_bracket_13 + 35, sprintf('p = %.3f', p_values(2)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Third comparison (Bin 1 vs Bin 4)
plot([x_positions(3)-0.15, x_positions(3)-.15, x_positions(3)+.15, x_positions(3)+.15], ...
     [y_bracket_14, y_bracket_14 + 10, y_bracket_14 + 10, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(3), x_positions(3)]), y_bracket_14 + 35, sprintf('p = %.3f', p_values(3)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Fourth comparison (Bin 2 vs Bin 3)
plot([x_positions(4)-0.15, x_positions(4)-.15, x_positions(4)+.15, x_positions(4)+.15], ...
     [y_bracket_23, y_bracket_23 + 10, y_bracket_23 + 10, y_bracket_23], 'k', 'LineWidth', 1.5);
text(mean([x_positions(4), x_positions(4)]), y_bracket_23 + 35, sprintf('p = %.3f', p_values(4)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Fifth comparison (Bin 2 vs Bin 4)
plot([x_positions(5)-0.15, x_positions(5)-.15, x_positions(5)+.15, x_positions(5)+.15], ...
     [y_bracket_24, y_bracket_24 + 10, y_bracket_24 + 10, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(5), x_positions(5)]), y_bracket_24 + 35, sprintf('p = %.3f', p_values(5)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Sixth comparison (Bin 3 vs Bin 4)
plot([x_positions(6)-0.15, x_positions(6)-.15, x_positions(6)+.15, x_positions(6)+.15], ...
     [y_bracket_34, y_bracket_34 + 10, y_bracket_34 + 10, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean([x_positions(6), x_positions(6)]), y_bracket_34 + 35, sprintf('p = %.3f', p_values(6)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
legend('Impatient', 'Patient')





for i=1:size(heu_time,2)
mean_resp(1,i) = mean(Response_L(heu_risk1(:,i))); 
mean_resp(2,i)= mean(Response_L(heu_risk2(:,i)));
end


% Assuming mean_WARP is already computed as provided in the previous code

% Perform t-tests for each comparison and store p-values
p_values = zeros(1, 6); % Preallocate p-values array for required comparisons

% Conduct t-tests for specific pairs between decisions with and without heuristics
% Replace these comparisons with the relevant Raven data as needed
[~, p_values(1)] = ttest2(Response_L(heu_risk1(:, 1)), Response_L(heu_risk2(:, 1)),'Tail','left');
[~, p_values(2)] = ttest2(Response_L(heu_risk1(:, 2)), Response_L(heu_risk2(:, 2)),'Tail','left');
[~, p_values(3)] = ttest2(Response_L(heu_risk1(:, 3)), Response_L(heu_risk2(:, 3)),'Tail','left');
[~, p_values(4)] = ttest2(Response_L(heu_risk1(:, 4)), Response_L(heu_risk2(:, 4)),'Tail','left');
[~, p_values(5)] = ttest2(Response_L(heu_risk1(:, 5)), Response_L(heu_risk2(:, 5)),'Tail','left');
[~, p_values(6)] = ttest2(Response_L(heu_risk1(:, 6)), Response_L(heu_risk2(:, 6)),'Tail','left');

% Plot the bar graph
figure;
bar(cat, mean_resp);

% Set the x-tick labels' font to Times New Roman
set(gca, 'FontName', 'Times New Roman');

% Set the y-axis label with the desired font settings
ylabel('Response times', 'FontName', 'Times New Roman');

% Set the limits for the y-axis; adjust as per the expected range of your data
ylim([0 1400]);

% Define x positions for the bars
x_positions = 1:length(mean_WARP);

% Define y positions for each comparison (adjust these as needed)
y_bracket_12 = 1100;   % Position for comparison Bin 1 vs Bin 2
y_bracket_13 = 1100;   % Position for comparison Bin 1 vs Bin 3
y_bracket_14 = 1100;   % Position for comparison Bin 1 vs Bin 4
y_bracket_23 = 1100;   % Position for comparison Bin 2 vs Bin 3
y_bracket_24 = 1100;   % Position for comparison Bin 2 vs Bin 4
y_bracket_34 = 1100;   % Position for comparison Bin 3 vs Bin 4

hold on;

% Add comparisons and p-values for each

% First comparison (Bin 1 vs Bin 2)
plot([x_positions(1)-0.15, x_positions(1)-.15, x_positions(1)+.15, x_positions(1)+.15], ...
     [y_bracket_12, y_bracket_12 + 20, y_bracket_12 + 20, y_bracket_12], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(1)]), y_bracket_12 + 70, sprintf('p = %.3f', p_values(1)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Second comparison (Bin 1 vs Bin 3)
plot([x_positions(2)-0.15, x_positions(2)-.15, x_positions(2)+.15, x_positions(2)+.15], ...
     [y_bracket_13, y_bracket_13 + 20, y_bracket_13 + 20, y_bracket_13], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(2)]), y_bracket_13 + 70, sprintf('p = %.3f', p_values(2)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Third comparison (Bin 1 vs Bin 4)
plot([x_positions(3)-0.15, x_positions(3)-.15, x_positions(3)+.15, x_positions(3)+.15], ...
     [y_bracket_14, y_bracket_14 + 20, y_bracket_14 + 20, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(3), x_positions(3)]), y_bracket_14 + 70, sprintf('p = %.3f', p_values(3)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Fourth comparison (Bin 2 vs Bin 3)
plot([x_positions(4)-0.15, x_positions(4)-.15, x_positions(4)+.15, x_positions(4)+.15], ...
     [y_bracket_23, y_bracket_23 + 20, y_bracket_23 + 20, y_bracket_23], 'k', 'LineWidth', 1.5);
text(mean([x_positions(4), x_positions(4)]), y_bracket_23 + 70, sprintf('p = %.3f', p_values(4)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Fifth comparison (Bin 2 vs Bin 4)
plot([x_positions(5)-0.15, x_positions(5)-.15, x_positions(5)+.15, x_positions(5)+.15], ...
     [y_bracket_24, y_bracket_24 + 20, y_bracket_24 + 20, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(5), x_positions(5)]), y_bracket_24 + 70, sprintf('p = %.3f', p_values(5)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

% Sixth comparison (Bin 3 vs Bin 4)
plot([x_positions(6)-0.15, x_positions(6)-.15, x_positions(6)+.15, x_positions(6)+.15], ...
     [y_bracket_34, y_bracket_34 + 20, y_bracket_34 + 20, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean([x_positions(6), x_positions(6)]), y_bracket_34 + 70, sprintf('p = %.3f', p_values(6)), ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

legend('Risk-averse', 'Risk-neutral')





