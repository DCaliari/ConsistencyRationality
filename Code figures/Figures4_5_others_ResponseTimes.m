clear
clc

structural_estimation_3rdSubmission_CRRACARAEXPBETA  % estimates under the different specifications

%% Figures 4 and 5 - Response Times %%


Times = xlsread('DATI.xlsx','Variables','BN1:BO146');
Response_D = Times(:,1);
Response_L = Times(:,2);



%% FIGURE 4

ti = ti_exp';
tp = tp_exp';


% Define bin edges and labels
bin_labels = {'Impatient', 'Others', 'Patient'};

% Initialize array to store mean Raven scores for each bin
mean_times = zeros(1, length(bin_labels));
se_times = zeros(1, length(bin_labels));

bin(:,2) = ( ~ti&~tp );
bin(:,1) = ( ti ) ;
bin(:,3) = (  tp ) ;


% Calculate mean "Raven" score for each bin
for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    % Calculate mean "Response_D" for entries within this bin
    mean_times(i) = mean(Response_D(bin(:,i)));
    se_times(i) = std(Response_D(bin(:,i))) / sqrt(length(Response_D(bin(:,i))));
end

% Plot the bar graph
bar(mean_times)
hold on
errorbar([1 2 3], mean_times, se_times, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('Response Times', 'FontName', 'Times New Roman');
ylim([300 530])

% Perform t-tests and store results for p-value annotations
[~, p12, ~] = ttest2(Response_D(bin(:,1)), Response_D(bin(:,2)));
[~, p23, ~] = ttest2(Response_D(bin(:,2)), Response_D(bin(:,3)));
[~, p13, ~] = ttest2(Response_D(bin(:,1)), Response_D(bin(:,3)));
% Define x positions for the bars
x_positions = 1:length(mean_times); 

% Y positions for each comparison, slightly above the bar heights
y_bracket_12 = 470; % Bin 1 vs Bin 2
y_bracket_23 = 490; % Bin 2 vs Bin 3
y_bracket_13 = 510; % Bin 1 vs Bin 3
hold on
% First comparison (Bin 1 vs Bin 2)
plot([x_positions(1), x_positions(1), x_positions(2), x_positions(2)], [y_bracket_12, y_bracket_12+5, y_bracket_12+5, y_bracket_12], 'k', 'LineWidth', 1.5);
text(mean(x_positions(1:2)), y_bracket_12 + 10, sprintf('p = %.3f', p12), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Second comparison (Bin 2 vs Bin 3)
plot([x_positions(2), x_positions(2), x_positions(3), x_positions(3)], [y_bracket_23, y_bracket_23+5, y_bracket_23+5, y_bracket_23], 'k', 'LineWidth', 1.5);
text(mean(x_positions(2:3)), y_bracket_23 + 10, sprintf('p = %.3f', p23), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Third comparison (Bin 1 vs Bin 3)
plot([x_positions(1), x_positions(1), x_positions(3), x_positions(3)], [y_bracket_13, y_bracket_13+5, y_bracket_13+5, y_bracket_13], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(3)]), y_bracket_13 + 10, sprintf('p = %.3f', p13), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold off


%% FIGURE 5

ra = ra_crra';
rl = rl_crra';


% Define bin edges and labels
bin_labels = {'Risk-neutral', 'Others', 'Risk-averse'};

% Initialize array to store mean Raven scores for each bin
mean_times = zeros(1, length(bin_labels));
se_times = zeros(1, length(bin_labels));

bin(:,2) = ( ~ra&~rl );
bin(:,3) = ( ra ) ;
bin(:,1) = (  rl ) ;


% Calculate mean "Raven" score for each bin
for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    mean_times(i) = mean(Response_L(bin(:,i)));
    se_times(i) = std(Response_L(bin(:,i))) / sqrt(length(Response_L(bin(:,i))));
end

% Plot the bar graph
bar(mean_times)
hold on
errorbar([1 2 3], mean_times, se_times, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('Response times', 'FontName', 'Times New Roman');
ylim([400 1000])


% Perform t-tests and store results for p-value annotations
[h12, p12, ~] = ttest2(Response_L(bin(:,1)), Response_L(bin(:,2)));
[h23, p23, ~] = ttest2(Response_L(bin(:,2)), Response_L(bin(:,3)));
[h13, p13, ~] = ttest2(Response_L(bin(:,1)), Response_L(bin(:,3)));

% Define x positions for the bars
x_positions = 1:length(mean_times); 

% Y positions for each comparison, slightly above the bar heights
y_bracket_12 = 860; % Bin 1 vs Bin 2
y_bracket_23 = 895; % Bin 2 vs Bin 3
y_bracket_13 = 930; % Bin 1 vs Bin 3
hold on
% First comparison (Bin 1 vs Bin 2)
plot([x_positions(1), x_positions(1), x_positions(2), x_positions(2)], [y_bracket_12, y_bracket_12+7.5, y_bracket_12+7.5, y_bracket_12], 'k', 'LineWidth', 1.5);
text(mean(x_positions(1:2)), y_bracket_12 + 25, sprintf('p = %.3f', p12), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Second comparison (Bin 2 vs Bin 3)
plot([x_positions(2), x_positions(2), x_positions(3), x_positions(3)], [y_bracket_23, y_bracket_23+7.5, y_bracket_23+7.5, y_bracket_23], 'k', 'LineWidth', 1.5);
text(mean(x_positions(2:3)), y_bracket_23 + 25, sprintf('p = %.3f', p23), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Third comparison (Bin 1 vs Bin 3)
plot([x_positions(1), x_positions(1), x_positions(3), x_positions(3)], [y_bracket_13, y_bracket_13+7.5, y_bracket_13+7.5, y_bracket_13], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(3)]), y_bracket_13 + 25, sprintf('p = %.3f', p13), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold off

