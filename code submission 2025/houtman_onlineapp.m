
%% FIGURE 11 



load houtman_swaps_stat_ind.mat



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
    mean_warp(i) = mean(HMTime(bin(:,i)));
    se_warp(i) = std(HMTime(bin(:,i))) / sqrt(length(HMTime(bin(:,i))));
end

% Plot the bar graph

subplot(1,2,1)

bar(mean_warp);
hold on
errorbar([1 2 3 4], mean_warp, se_warp, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('WARP violations', 'FontName', 'Times New Roman');
ylim([0 7])

% Perform t-tests and store results for p-value annotations

[~, p34, ~] = ttest2(HMTime(bin(:,3)), HMTime(bin(:,4)),'Tail','right');
[~, p14, ~] = ttest2(HMTime(bin(:,1)), HMTime(bin(:,4)),'Tail','left');
[~, p24, ~] = ttest2(HMTime(bin(:,2)), HMTime(bin(:,4)),'Tail','left');

% Define x positions for the bars
x_positions = 1:length(mean_warp);

y_bracket_34 = 6;  % Bin 3 vs Bin 4
y_bracket_14 = 4;  % Bin 1 vs Bin 4
y_bracket_24 = 5; % Bin 2 vs Bin 4

hold on


plot([x_positions(1), x_positions(1), x_positions(4), x_positions(4)], [y_bracket_14, y_bracket_14+0.10, y_bracket_14+0.10, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(4)]), y_bracket_14 + 0.35, sprintf('p = %.3f', p14), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Sixth comparison (Bin 2 vs Bin 4)
plot([x_positions(2), x_positions(2), x_positions(4), x_positions(4)], [y_bracket_24, y_bracket_24+0.10, y_bracket_24+0.10, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(4)]), y_bracket_24 + 0.35, sprintf('p = %.3f', p24), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Third comparison (Bin 3 vs Bin 4)
plot([x_positions(3), x_positions(3), x_positions(4), x_positions(4)], [y_bracket_34, y_bracket_34+0.10, y_bracket_34+0.10, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean(x_positions(3:4)), y_bracket_34 + 0.35, sprintf('p = %.3f', p34), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
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
    mean_warp(i) = mean(HMRisk(bin(:,i)));
    se_warp(i) = std(HMRisk(bin(:,i))) / sqrt(length(HMRisk(bin(:,i))));
end


bar(mean_warp);
hold on
errorbar([1 2 3 4], mean_warp, se_warp, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('WARP violations', 'FontName', 'Times New Roman');
ylim([0 7])

% Perform t-tests and store results for p-value annotations

[~, p34, ~] = ttest2(HMRisk(bin(:,3)), HMRisk(bin(:,4)),'Tail','right');
[~, p14, ~] = ttest2(HMRisk(bin(:,1)), HMRisk(bin(:,4)),'Tail','left');
[~, p24, ~] = ttest2(HMRisk(bin(:,2)), HMRisk(bin(:,4)),'Tail','left');

% Define x positions for the bars
x_positions = 1:length(mean_warp);

y_bracket_34 = 6;  % Bin 3 vs Bin 4
y_bracket_14 = 4;  % Bin 1 vs Bin 4
y_bracket_24 = 5; % Bin 2 vs Bin 4

hold on

plot([x_positions(1), x_positions(1), x_positions(4), x_positions(4)], [y_bracket_14, y_bracket_14+0.10, y_bracket_14+0.10, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(4)]), y_bracket_14 + 0.35, sprintf('p = %.3f', p14), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Sixth comparison (Bin 2 vs Bin 4)
plot([x_positions(2), x_positions(2), x_positions(4), x_positions(4)], [y_bracket_24, y_bracket_24+0.10, y_bracket_24+0.10, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(4)]), y_bracket_24 + 0.35, sprintf('p = %.3f', p24), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Third comparison (Bin 3 vs Bin 4)
plot([x_positions(3), x_positions(3), x_positions(4), x_positions(4)], [y_bracket_34, y_bracket_34+0.10, y_bracket_34+0.10, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean(x_positions(3:4)), y_bracket_34 + 0.35, sprintf('p = %.3f', p34), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on

