clear
clc

structural_estimation_3rdSubmission_CRRACARAEXPBETA  % estimates under the different specifications

%% Figures 11 and 12 - WARP violations and Raven scores %%


deliberate_time =  xlsread('DATI.xlsx','Variables','FN1:FN146');
deliberate = xlsread('DATI.xlsx','Variables','EE1:EE146');

% McNemar test , differences in proportions

x = [sum(deliberate_time==1 & deliberate==1), sum(deliberate_time==1 & deliberate==0);
    sum(deliberate_time==0 & deliberate==1), sum(deliberate_time==0 & deliberate==0)];

mcnemar(x) % Test for the difference in proportion of deliberate randomization behavior between Time and Risk

WARP_D = xlsread('DATI.xlsx','Variables','C1:C146');
WARP_L = xlsread('DATI.xlsx','Variables','E1:E146');
Raven = xlsread('DATI.xlsx','Variables','S1:S146');


%% FIGURE 11 

ti = ti_exp';
tp = tp_exp';



% Define bin edges and labels
bin_labels = {'Impatient', 'Patient', 'Randomize', 'Others'};

% Initialize array to store mean Response_D scores for each bin
mean_warp = zeros(1, length(bin_labels));

% Calculate mean "Response_D" score for each bin

bin(:,3) = deliberate_time;
bin(:,1) = ( ti & deliberate_time==0 ) ;
bin(:,2) = (  tp & deliberate_time==0  ) ;
bin(:,4) = bin(:,1)+bin(:,2)+bin(:,3)==0;


bin=logical(bin);

for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    % Calculate mean "Response_D" for entries within this bin
    mean_warp(i) = mean(WARP_D(bin(:,i)));
    se_warp(i) = std(WARP_D(bin(:,i))) / sqrt(length(WARP_D(bin(:,i))));
end

% Plot the bar graph

subplot(1,2,1)

bar(mean_warp);
hold on
errorbar([1 2 3 4], mean_warp, se_warp, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('WARP violations', 'FontName', 'Times New Roman');
ylim([0 12])

% Perform t-tests and store results for p-value annotations

[~, p34, ~] = ttest2(WARP_D(bin(:,3)), WARP_D(bin(:,4)),'Tail','right');
[~, p14, ~] = ttest2(WARP_D(bin(:,1)), WARP_D(bin(:,4)),'Tail','left');
[~, p24, ~] = ttest2(WARP_D(bin(:,2)), WARP_D(bin(:,4)),'Tail','left');

% Define x positions for the bars
x_positions = 1:length(mean_warp);

y_bracket_34 = 10;  % Bin 3 vs Bin 4
y_bracket_14 = 8;  % Bin 1 vs Bin 4
y_bracket_24 = 9; % Bin 2 vs Bin 4

hold on

plot([x_positions(1), x_positions(1), x_positions(4), x_positions(4)], [y_bracket_14, y_bracket_14+0.15, y_bracket_14+0.15, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(4)]), y_bracket_14 + 0.5, sprintf('p = %.3f', p14), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Sixth comparison (Bin 2 vs Bin 4)
plot([x_positions(2), x_positions(2), x_positions(4), x_positions(4)], [y_bracket_24, y_bracket_24+0.15, y_bracket_24+0.15, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(4)]), y_bracket_24 + 0.5, sprintf('p = %.3f', p24), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Third comparison (Bin 3 vs Bin 4)
plot([x_positions(3), x_positions(3), x_positions(4), x_positions(4)], [y_bracket_34, y_bracket_34+0.15, y_bracket_34+0.15, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean(x_positions(3:4)), y_bracket_34 + 0.5, sprintf('p = %.3f', p34), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on



subplot(1,2,2)

bin_labels = {'Risk averse', 'Risk neutral', 'Randomize', 'Others'};

% Initialize array to store mean Response_D scores for each bin
mean_warp = zeros(1, length(bin_labels));

% Calculate mean "Response_D" score for each bin

ra = ra_crra';
rl = rl_crra';

bin(:,3) = deliberate;
bin(:,2) = ( rl & deliberate ==0 ) ;
bin(:,1) = (  ra  & deliberate ==0 ) ;
bin(:,4) = bin(:,1)+bin(:,2)+bin(:,3)==0;

bin=logical(bin);

for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    % Calculate mean "Response_D" for entries within this bin
    mean_warp(i) = mean(WARP_L(bin(:,i)));
    se_warp(i) = std(WARP_L(bin(:,i))) / sqrt(length(WARP_L(bin(:,i))));
end


bar(mean_warp);
hold on
errorbar([1 2 3 4], mean_warp, se_warp, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('WARP violations', 'FontName', 'Times New Roman');
ylim([0 12])

% Perform t-tests and store results for p-value annotations

[~, p34, ~] = ttest2(WARP_L(bin(:,3)), WARP_L(bin(:,4)),'Tail','right');
[~, p14, ~] = ttest2(WARP_L(bin(:,1)), WARP_L(bin(:,4)),'Tail','left');
[~, p24, ~] = ttest2(WARP_L(bin(:,2)), WARP_L(bin(:,4)),'Tail','left');

% Define x positions for the bars
x_positions = 1:length(mean_warp);

y_bracket_34 = 10;  % Bin 3 vs Bin 4
y_bracket_14 = 8;  % Bin 1 vs Bin 4
y_bracket_24 = 9; % Bin 2 vs Bin 4

hold on

plot([x_positions(1), x_positions(1), x_positions(4), x_positions(4)], [y_bracket_14, y_bracket_14+0.15, y_bracket_14+0.15, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(4)]), y_bracket_14 + 0.5, sprintf('p = %.3f', p14), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Sixth comparison (Bin 2 vs Bin 4)
plot([x_positions(2), x_positions(2), x_positions(4), x_positions(4)], [y_bracket_24, y_bracket_24+0.15, y_bracket_24+0.15, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(4)]), y_bracket_24 + 0.5, sprintf('p = %.3f', p24), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Third comparison (Bin 3 vs Bin 4)
plot([x_positions(3), x_positions(3), x_positions(4), x_positions(4)], [y_bracket_34, y_bracket_34+0.15, y_bracket_34+0.15, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean(x_positions(3:4)), y_bracket_34 + 0.5, sprintf('p = %.3f', p34), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on


%% FIGURE 12

clearvars bin_labels
% Define bin edges and labels
bin_labels = {'Impatient', 'Patient', 'Randomize', 'Others'};

bin(:,3) = deliberate_time;
bin(:,1) = ( ti & deliberate_time==0 ) ;
bin(:,2) = (  tp & deliberate_time==0  ) ;
bin(:,4) = bin(:,1)+bin(:,2)+bin(:,3)==0;


bin=logical(bin);

% Initialize array to store mean Response_D scores for each bin
mean_raven = zeros(1, length(bin_labels));


for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    % Calculate mean "Response_D" for entries within this bin
    mean_raven(i) = mean(Raven(bin(:,i)));
    se_raven(i) = std(Raven(bin(:,i))) / sqrt(length(Raven(bin(:,i))));
end




subplot(1,2,1)

bar(mean_raven);
hold on
errorbar([1 2 3 4], mean_raven, se_raven, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('Raven Scores', 'FontName', 'Times New Roman');
ylim([0 18])

% Perform t-tests and store results for p-value annotations
% [~, p12, ~] = ttest2(Raven(bin(:,1)), Raven(bin(:,2)),'Tail','left');
%[h23, p23, ~] = ttest2(Raven(bin(:,2)), Raven(bin(:,3)));
[~, p34, ~] = ttest2(Raven(bin(:,3)), Raven(bin(:,4)),'Tail','right');
%[h13, p13, ~] = ttest2(Raven(bin(:,1)), Raven(bin(:,3)));
[~, p14, ~] = ttest2(Raven(bin(:,1)), Raven(bin(:,4)),'Tail','left');
[~, p24, ~] = ttest2(Raven(bin(:,2)), Raven(bin(:,4)),'Tail','left');

% Define x positions for the bars
x_positions = 1:length(mean_raven);

y_bracket_34 = 16;  % Bin 3 vs Bin 4
y_bracket_14 = 12;  % Bin 1 vs Bin 4
y_bracket_24 = 14; % Bin 2 vs Bin 4

hold on
plot([x_positions(3), x_positions(3), x_positions(4), x_positions(4)], [y_bracket_34, y_bracket_34+0.25, y_bracket_34+0.25, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean(x_positions(3:4)), y_bracket_34 + 0.75, sprintf('p = %.3f', p34), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
plot([x_positions(1), x_positions(1), x_positions(4), x_positions(4)], [y_bracket_14, y_bracket_14+0.25, y_bracket_14+0.25, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(4)]), y_bracket_14 + 0.75, sprintf('p = %.3f', p14), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Sixth comparison (Bin 2 vs Bin 4)
plot([x_positions(2), x_positions(2), x_positions(4), x_positions(4)], [y_bracket_24, y_bracket_24+0.25, y_bracket_24+0.25, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(4)]), y_bracket_24 + 0.75, sprintf('p = %.3f', p24), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on



bin_labels = {'Risk averse', 'Risk neutral', 'Randomize', 'Others'};

% Initialize array to store mean Response_D scores for each bin
mean_raven = zeros(1, length(bin_labels));


bin(:,3) = deliberate;
bin(:,2) = ( rl & deliberate ==0 ) ;
bin(:,1) = (  ra  & deliberate ==0 ) ;
bin(:,4) = bin(:,1)+bin(:,2)+bin(:,3)==0;

bin=logical(bin);

for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    % Calculate mean "Response_D" for entries within this bin
    mean_raven(i) = mean(Raven(bin(:,i)));
    se_raven(i) = std(Raven(bin(:,i))) / sqrt(length(Raven(bin(:,i))));
end



subplot(1,2,2)

bar(mean_raven);
hold on
errorbar([1 2 3 4], mean_raven, se_raven, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('Raven Scores', 'FontName', 'Times New Roman');
ylim([0 18])

% Perform t-tests and store results for p-value annotations
[h12, p12, ~] = ttest2(Raven(bin(:,1)), Raven(bin(:,2)),'Tail','left');
%[h23, p23, ~] = ttest2(Raven(bin(:,2)), Raven(bin(:,3)));
[h34, p34, ~] = ttest2(Raven(bin(:,3)), Raven(bin(:,4)),'Tail','right');
%[h13, p13, ~] = ttest2(Raven(bin(:,1)), Raven(bin(:,3)));
[h14, p14, ~] = ttest2(Raven(bin(:,1)), Raven(bin(:,4)),'Tail','left');
[h24, p24, ~] = ttest2(Raven(bin(:,2)), Raven(bin(:,4)),'Tail','left');

% Define x positions for the bars
x_positions = 1:length(mean_raven);


y_bracket_34 = 16.2;  % Bin 3 vs Bin 4
y_bracket_14 = 12.2;  % Bin 1 vs Bin 4
y_bracket_24 = 14.2; % Bin 2 vs Bin 4

hold on
plot([x_positions(3), x_positions(3), x_positions(4), x_positions(4)], [y_bracket_34, y_bracket_34+0.25, y_bracket_34+0.25, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean(x_positions(3:4)), y_bracket_34 + 0.75, sprintf('p = %.3f', p34), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
plot([x_positions(1), x_positions(1), x_positions(4), x_positions(4)], [y_bracket_14, y_bracket_14+0.25, y_bracket_14+0.25, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(4)]), y_bracket_14 + 0.75, sprintf('p = %.3f', p14), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Sixth comparison (Bin 2 vs Bin 4)
plot([x_positions(2), x_positions(2), x_positions(4), x_positions(4)], [y_bracket_24, y_bracket_24+0.25, y_bracket_24+0.25, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(4)]), y_bracket_24 + 0.75, sprintf('p = %.3f', p24), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on

