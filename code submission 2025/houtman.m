

% Set default font to Times New Roman for the entire figure
set(gca,'FontName','Times New Roman')
subplot(1,2,1)
scatter(HMTime, Raven)
xlim([0 6])
[r1, p1] = corr(HMTime, Raven);  % Calculate correlation and p-value
text(1, max(Raven)*0.9, ['r = ', num2str(r1,'%.2f'), ', p = ', num2str(p1,'%.2f')], 'FontName', 'Times New Roman')  % Display correlation and p-value
xlabel('WARP violations', 'FontName', 'Times New Roman')
ylabel('Raven', 'FontName', 'Times New Roman')

subplot(1,2,2)
scatter(HMRisk, Raven)
xlim([0 6])
[r2, p2] = corr(HMRisk, Raven);  % Calculate correlation and p-value
text(1, max(Raven)*0.9, ['r = ', num2str(r2,'%.2f'), ', p = ', num2str(p2,'%.2f')], 'FontName', 'Times New Roman')  % Display correlation and p-value
xlabel('WARP violations', 'FontName', 'Times New Roman')
ylabel('Raven', 'FontName', 'Times New Roman')





WARP_T = HMTime + HMRisk;

id_heu = time_ind < 0.91 | risk_ind >2.2 | deliberate_time' ==1 |  deliberate' ==1 | risk_ind< 0.2 | time_ind>0.99;
id_rum = time_ind >= 0.91 & time_ind <= 0.99 & risk_ind <= 2.2 & risk_ind >=0.2 & deliberate' ==0 & deliberate_time' ==0;

n1 = sum(id_heu);
n2 = sum(id_rum);

prec_ind_risk = log(prec_ind_risk);


% Set default font to Times New Roman for the entire figure
set(gca,'FontName','Times New Roman')
subplot(1,2,1)
scatter(WARP_T(id_heu), Raven(id_heu))
xlim([0 24])
[r1, p1] = corr(WARP_T(id_heu), Raven(id_heu));  % Calculate correlation and p-value
text(1, max(Raven)*0.9, ['r = ', num2str(r1,'%.2f'), ', p = ', num2str(p1,'%.2f')], 'FontName', 'Times New Roman')  % Display correlation and p-value
xlabel(['WARP violations (n = ', num2str(n1), ')'], 'FontName', 'Times New Roman');  % For the first subplot
ylabel('Raven', 'FontName', 'Times New Roman')

subplot(1,2,2)
scatter(WARP_T(id_rum), Raven(id_rum))
xlim([0 24])
[r2, p2] = corr(WARP_T(id_rum), Raven(id_rum));  % Calculate correlation and p-value
text(1, max(Raven)*0.9, ['r = ', num2str(r2,'%.2f'), ', p = ', num2str(p2,'%.2f')], 'FontName', 'Times New Roman')  % Display correlation and p-value
ylabel('Raven', 'FontName', 'Times New Roman')
xlabel(['WARP violations (n = ', num2str(n2), ')'], 'FontName', 'Times New Roman');  % For the second subplot

compare_correlation_coefficients(  r1, r2, n1, n2  )




id_heu = time_ind < 0.91 | deliberate_time' ==1 |  time_ind>0.99;
id_rum = time_ind >= 0.91 & time_ind <= 0.99 & deliberate_time' ==0;

n1 = sum(id_heu);
n2 = sum(id_rum);


% Set default font to Times New Roman for the entire figure
set(gca,'FontName','Times New Roman')
subplot(1,2,1)
scatter(HMTime(id_heu), Raven(id_heu))
xlim([0 12])
[r1, p1] = corr(HMTime(id_heu), Raven(id_heu));  % Calculate correlation and p-value
text(1, max(Raven)*0.9, ['r = ', num2str(r1,'%.2f'), ', p = ', num2str(p1,'%.2f')], 'FontName', 'Times New Roman')  % Display correlation and p-value
xlabel(['WARP violations (n = ', num2str(n1), ')'], 'FontName', 'Times New Roman');  % For the first subplot
ylabel('Raven', 'FontName', 'Times New Roman')

subplot(1,2,2)
scatter(HMTime(id_rum), Raven(id_rum))
xlim([0 12])
[r2, p2] = corr(HMTime(id_rum), Raven(id_rum));  % Calculate correlation and p-value
text(1, max(Raven)*0.9, ['r = ', num2str(r2,'%.2f'), ', p = ', num2str(p2,'%.2f')], 'FontName', 'Times New Roman')  % Display correlation and p-value
ylabel('Raven', 'FontName', 'Times New Roman')
xlabel(['WARP violations (n = ', num2str(n2), ')'], 'FontName', 'Times New Roman');  % For the second subplot

compare_correlation_coefficients(  r1, r2, n1, n2  )




id_heu =  risk_ind >2.2 |  deliberate' ==1 ;
id_rum =  risk_ind <= 2.2 & risk_ind >=0.2 & deliberate' ==0;

n1 = sum(id_heu);
n2 = sum(id_rum);



% Set default font to Times New Roman for the entire figure
set(gca,'FontName','Times New Roman')
subplot(1,2,1)
scatter(HMRisk(id_heu), Raven(id_heu))
xlim([0 12])
[r1, p1] = corr(HMRisk(id_heu), Raven(id_heu));  % Calculate correlation and p-value
text(1, max(Raven)*0.9, ['r = ', num2str(r1,'%.2f'), ', p = ', num2str(p1,'%.2f')], 'FontName', 'Times New Roman')  % Display correlation and p-value
xlabel(['WARP violations (n = ', num2str(n1), ')'], 'FontName', 'Times New Roman');  % For the first subplot
ylabel('Raven', 'FontName', 'Times New Roman')

subplot(1,2,2)
scatter(HMRisk(id_rum), Raven(id_rum))
xlim([0 12])
[r2, p2] = corr(HMRisk(id_rum), Raven(id_rum));  % Calculate correlation and p-value
text(1, max(Raven)*0.9, ['r = ', num2str(r2,'%.2f'), ', p = ', num2str(p2,'%.2f')], 'FontName', 'Times New Roman')  % Display correlation and p-value
ylabel('Raven', 'FontName', 'Times New Roman')
xlabel(['WARP violations (n = ', num2str(n2), ')'], 'FontName', 'Times New Roman');  % For the second subplot

compare_correlation_coefficients(  r1, r2, n1, n2  )



%% ALL TOGETHER


id_heu = time_ind < 0.91 | risk_ind >2.2 | deliberate_time' ==1 |  deliberate' ==1 | risk_ind< 0.2 | time_ind>0.99;
id_rum = time_ind >= 0.91 & time_ind <= 0.99 & risk_ind <= 2.2 & risk_ind >=0.2 & deliberate' ==0 & deliberate_time' ==0;

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
text(1, max(Raven)*0.9, ['Heu & Rand: r = ', num2str(r1, '%.2f'), ', p = ', num2str(p1, '%.2f')], 'Color', 'blue', 'FontName', 'Times New Roman');
text(1, max(Raven)*0.8, ['Rum: r = ', num2str(r2, '%.2f'), ', p = ', num2str(p2, '%.2f')], 'Color', 'red', 'FontName', 'Times New Roman');

title('Aggregate', 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel(['WARP violations (n1 = ', num2str(n1), ', n2 = ', num2str(n2), ')'], 'FontName', 'Times New Roman');
ylabel('Raven', 'FontName', 'Times New Roman');
legend('Heu & Rand', 'Rum', 'Location', 'best');
hold off;


id_heu = time_ind < 0.91 | deliberate_time' ==1 |  time_ind>0.99;
id_rum = time_ind >= 0.91 & time_ind <= 0.99 & deliberate_time' ==0;

n1 = sum(id_heu);
n2 = sum(id_rum);


% Smaller panels for "Time" and "Risk"
% "Time" panel (HMTime)
nexttile
scatter(HMTime(id_heu), Raven(id_heu));
hold on;
scatter(HMTime(id_rum), Raven(id_rum));

% Regression lines for "Time" panel
% Heuristic group regression line
p_heu_time = polyfit(HMTime(id_heu), Raven(id_heu), 1);
yfit_heu_time = polyval(p_heu_time, HMTime(id_heu));
plot(HMTime(id_heu), yfit_heu_time, 'b-', 'LineWidth', 1.5);

% Rumination group regression line
p_rum_time = polyfit(HMTime(id_rum), Raven(id_rum), 1);
yfit_rum_time = polyval(p_rum_time, HMTime(id_rum));
plot(HMTime(id_rum), yfit_rum_time, 'r-', 'LineWidth', 1.5);

% Correlation coefficients for "Time" panel
[r1_time, p1_time] = corr(HMTime(id_heu), Raven(id_heu));
[r2_time, p2_time] = corr(HMTime(id_rum), Raven(id_rum));

% Display correlation coefficients on "Time" panel
text(1, max(Raven)*0.95, ['Heu & Rand: r = ', num2str(r1_time, '%.2f'), ', p = ', num2str(p1_time, '%.2f')], 'Color', 'blue', 'FontName', 'Times New Roman');
text(1, max(Raven)*0.85, ['Rum: r = ', num2str(r2_time, '%.2f'), ', p = ', num2str(p2_time, '%.2f')], 'Color', 'red', 'FontName', 'Times New Roman');

title('Time', 'FontName', 'Times New Roman');
xlabel(['WARP violations (n1 = ', num2str(n1), ', n2 = ', num2str(n2), ')'], 'FontName', 'Times New Roman');
ylabel('Raven', 'FontName', 'Times New Roman');
hold off;


id_heu =  risk_ind >2.2 |  deliberate' ==1 | risk_ind< 0.2;
id_rum =  risk_ind <= 2.2 & risk_ind >=0.2 & deliberate' ==0;

n1 = sum(id_heu);
n2 = sum(id_rum);




% "Risk" panel (HMRisk)
nexttile
scatter(HMRisk(id_heu), Raven(id_heu));
hold on;
scatter(HMRisk(id_rum), Raven(id_rum));

% Regression lines for "Risk" panel
% Heuristic group regression line
p_heu_risk = polyfit(HMRisk(id_heu), Raven(id_heu), 1);
yfit_heu_risk = polyval(p_heu_risk, HMRisk(id_heu));
plot(HMRisk(id_heu), yfit_heu_risk, 'b-', 'LineWidth', 1.5);

% Rumination group regression line
p_rum_risk = polyfit(HMRisk(id_rum), Raven(id_rum), 1);
yfit_rum_risk = polyval(p_rum_risk, HMRisk(id_rum));
plot(HMRisk(id_rum), yfit_rum_risk, 'r-', 'LineWidth', 1.5);

% Correlation coefficients for "Risk" panel
[r1_risk, p1_risk] = corr(HMRisk(id_heu), Raven(id_heu));
[r2_risk, p2_risk] = corr(HMRisk(id_rum), Raven(id_rum));

% Display correlation coefficients on "Risk" panel
text(1, max(Raven)*0.95, ['Heu & Rand: r = ', num2str(r1_risk, '%.2f'), ', p = ', num2str(p1_risk, '%.2f')], 'Color', 'blue', 'FontName', 'Times New Roman');
text(1, max(Raven)*0.85, ['Rum: r = ', num2str(r2_risk, '%.2f'), ', p = ', num2str(p2_risk, '%.2f')], 'Color', 'red', 'FontName', 'Times New Roman');

title('Risk', 'FontName', 'Times New Roman');
xlabel(['WARP violations (n1 = ', num2str(n1), ', n2 = ', num2str(n2), ')'], 'FontName', 'Times New Roman');
ylabel('Raven', 'FontName', 'Times New Roman');
hold off;












% Define bin edges and labels
bin_edges = [0.9, 0.91, 0.99, 1.0];
bin_labels = {'Impatient', 'Patient', 'Randomize', 'Others'};

% Initialize array to store mean Response_D scores for each bin
mean_warp = zeros(1, length(bin_labels));

% Calculate mean "Response_D" score for each bin

bin(:,3) = deliberate_time;
bin(:,1) = ( time_ind >= bin_edges(1) & time_ind < bin_edges(2) & deliberate_time'==0 ) ;
bin(:,2) = (  time_ind >= bin_edges(3) & time_ind < bin_edges(4) & deliberate_time'==0  ) ;
bin(:,4) = bin(:,1)+bin(:,2)+bin(:,3)==0;


bin=logical(bin);

for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    % Calculate mean "Response_D" for entries within this bin
    mean_warp(i) = mean(HMTime(bin(:,i)));
end

% Plot the bar graph
figure;

subplot(1,2,1)

bar(mean_warp);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('Houtman-Maks index', 'FontName', 'Times New Roman');
ylim([0 7])

% Perform t-tests and store results for p-value annotations
[h12, p12, ~] = ttest2(HMTime(bin(:,1)), HMTime(bin(:,2)));
[h23, p23, ~] = ttest2(HMTime(bin(:,2)), HMTime(bin(:,3)));
[h34, p34, ~] = ttest2(HMTime(bin(:,3)), HMTime(bin(:,4)),'Tail','right');
[h13, p13, ~] = ttest2(HMTime(bin(:,1)), HMTime(bin(:,3)));
[h14, p14, ~] = ttest2(HMTime(bin(:,1)), HMTime(bin(:,4)),'Tail','left');
[h24, p24, ~] = ttest2(HMTime(bin(:,2)), HMTime(bin(:,4)),'Tail','left');

% Define x positions for the bars
x_positions = 1:length(mean_warp);

% Y positions for each comparison, slightly above the bar heights
% y_bracket_12 = 7.5;  % Bin 1 vs Bin 2
% y_bracket_23 = 8.5;  % Bin 2 vs Bin 3
y_bracket_34 = 5.5;  % Bin 3 vs Bin 4
% y_bracket_13 = 10.5;  % Bin 1 vs Bin 3
y_bracket_14 = 3.5;  % Bin 1 vs Bin 4
y_bracket_24 = 4.5; % Bin 2 vs Bin 4


hold on

% % First comparison (Bin 1 vs Bin 2)
% plot([x_positions(1), x_positions(1), x_positions(2), x_positions(2)], [y_bracket_12, y_bracket_12+0.25, y_bracket_12+0.25, y_bracket_12], 'k', 'LineWidth', 1.5);
% text(mean(x_positions(1:2)), y_bracket_12 + 0.75, sprintf('p = %.3f', p12), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
% hold on
% % Second comparison (Bin 2 vs Bin 3)
% plot([x_positions(2), x_positions(2), x_positions(3), x_positions(3)], [y_bracket_23, y_bracket_23+0.25, y_bracket_23+0.25, y_bracket_23], 'k', 'LineWidth', 1.5);
% text(mean(x_positions(2:3)), y_bracket_13 + 0.75, sprintf('p = %.3f', p23), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
% hold on
% Fifth comparison (Bin 1 vs Bin 4)
plot([x_positions(1), x_positions(1), x_positions(4), x_positions(4)], [y_bracket_14, y_bracket_14+0.25, y_bracket_14+0.25, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(4)]), y_bracket_14 + 0.5, sprintf('p = %.3f', p14), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Sixth comparison (Bin 2 vs Bin 4)
plot([x_positions(2), x_positions(2), x_positions(4), x_positions(4)], [y_bracket_24, y_bracket_24+0.25, y_bracket_24+0.25, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(4)]), y_bracket_24 + 0.5, sprintf('p = %.3f', p24), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
% Third comparison (Bin 3 vs Bin 4)
plot([x_positions(3), x_positions(3), x_positions(4), x_positions(4)], [y_bracket_34, y_bracket_34+0.25, y_bracket_34+0.25, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean(x_positions(3:4)), y_bracket_34 + 0.5, sprintf('p = %.3f', p34), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% % Fourth comparison (Bin 1 vs Bin 3)
% plot([x_positions(1), x_positions(1), x_positions(3), x_positions(3)], [y_bracket_13, y_bracket_13+0.25, y_bracket_13+0.25, y_bracket_13], 'k', 'LineWidth', 1.5);
% text(mean([x_positions(1), x_positions(3)]), y_bracket_13 + 0.75, sprintf('p = %.3f', p13), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
% hold on






% Define bin edges and labels
bin_edges = [0, 0.2, 2.2, 2.4];
bin_labels = {'Risk averse', 'Risk neutral', 'Randomize', 'Others'};

% Initialize array to store mean Response_D scores for each bin
mean_warp = zeros(1, length(bin_labels));

% Calculate mean "Response_D" score for each bin

bin(:,3) = deliberate;
bin(:,2) = ( risk_ind >= bin_edges(1) & risk_ind < bin_edges(2) & deliberate' ==0 ) ;
bin(:,1) = (  risk_ind >= bin_edges(3) & risk_ind < bin_edges(4)  & deliberate' ==0 ) ;
bin(:,4) = bin(:,1)+bin(:,2)+bin(:,3)==0;


bin=logical(bin);

for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    % Calculate mean "Response_D" for entries within this bin
    mean_warp(i) = mean(HMRisk(bin(:,i)));
end

% Plot the bar graph
subplot(1,2,2)
bar(mean_warp);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('Houtman-Maks index', 'FontName', 'Times New Roman');
ylim([0 7])

% Perform t-tests and store results for p-value annotations
[h12, p12, ~] = ttest2(HMRisk(bin(:,1)), HMRisk(bin(:,2)));
[h23, p23, ~] = ttest2(HMRisk(bin(:,2)), HMRisk(bin(:,3)));
[h34, p34, ~] = ttest2(HMRisk(bin(:,3)), HMRisk(bin(:,4)),'Tail','right');
[h13, p13, ~] = ttest2(HMRisk(bin(:,1)), HMRisk(bin(:,3)));
[h14, p14, ~] = ttest2(HMRisk(bin(:,1)), HMRisk(bin(:,4)),'Tail','left');
[h24, p24, ~] = ttest2(HMRisk(bin(:,2)), HMRisk(bin(:,4)),'Tail','left');

% Define x positions for the bars
x_positions = 1:length(mean_warp);

% Y positions for each comparison, slightly above the bar heights
% y_bracket_12 = 7.5;  % Bin 1 vs Bin 2
% y_bracket_23 = 8.5;  % Bin 2 vs Bin 3
y_bracket_34 = 5.5;  % Bin 3 vs Bin 4
% y_bracket_13 = 10.5;  % Bin 1 vs Bin 3
y_bracket_14 = 3.5;  % Bin 1 vs Bin 4
y_bracket_24 = 4.5; % Bin 2 vs Bin 4


hold on

% % First comparison (Bin 1 vs Bin 2)
% plot([x_positions(1), x_positions(1), x_positions(2), x_positions(2)], [y_bracket_12, y_bracket_12+0.25, y_bracket_12+0.25, y_bracket_12], 'k', 'LineWidth', 1.5);
% text(mean(x_positions(1:2)), y_bracket_12 + 0.75, sprintf('p = %.3f', p12), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
% hold on
% % Second comparison (Bin 2 vs Bin 3)
% plot([x_positions(2), x_positions(2), x_positions(3), x_positions(3)], [y_bracket_23, y_bracket_23+0.25, y_bracket_23+0.25, y_bracket_23], 'k', 'LineWidth', 1.5);
% text(mean(x_positions(2:3)), y_bracket_13 + 0.75, sprintf('p = %.3f', p23), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
% hold on
% Fifth comparison (Bin 1 vs Bin 4)
plot([x_positions(1), x_positions(1), x_positions(4), x_positions(4)], [y_bracket_14, y_bracket_14+0.25, y_bracket_14+0.25, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(4)]), y_bracket_14 + 0.5, sprintf('p = %.3f', p14), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Sixth comparison (Bin 2 vs Bin 4)
plot([x_positions(2), x_positions(2), x_positions(4), x_positions(4)], [y_bracket_24, y_bracket_24+0.25, y_bracket_24+0.25, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(4)]), y_bracket_24 + 0.5, sprintf('p = %.3f', p24), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
% Third comparison (Bin 3 vs Bin 4)
plot([x_positions(3), x_positions(3), x_positions(4), x_positions(4)], [y_bracket_34, y_bracket_34+0.25, y_bracket_34+0.25, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean(x_positions(3:4)), y_bracket_34 + 0.5, sprintf('p = %.3f', p34), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% % Fourth comparison (Bin 1 vs Bin 3)
% plot([x_positions(1), x_positions(1), x_positions(3), x_positions(3)], [y_bracket_13, y_bracket_13+0.25, y_bracket_13+0.25, y_bracket_13], 'k', 'LineWidth', 1.5);
% text(mean([x_positions(1), x_positions(3)]), y_bracket_13 + 0.75, sprintf('p = %.3f', p13), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
% hold on
hold off





