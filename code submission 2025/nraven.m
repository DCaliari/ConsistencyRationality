

deliberate_time =  xlsread('DATI.xlsx','Variables','FN1:FN146');
deliberate = xlsread('DATI.xlsx','Variables','EE1:EE146');

% McNemar test , differences in proportions

x = [sum(deliberate_time==1 & deliberate==1), sum(deliberate_time==1 & deliberate==0);
    sum(deliberate_time==0 & deliberate==1), sum(deliberate_time==0 & deliberate==0)];

mcnemar(x) % Test for the difference in proportion of deliberate randomization behavior between Time and Risk

WARP_D = xlsread('DATI.xlsx','Variables','C1:C146');
WARP_L = xlsread('DATI.xlsx','Variables','E1:E146');
Raven = xlsread('DATI.xlsx','Variables','S1:S146');



Nraven = xlsread('DATI.xlsx','Variables','P1:P146');


ti = ti_exp';
tp = tp_exp';

ra= ra_crra';
rl = rl_crra';




% Define bin edges and labels

bin_labels = {'Impatient', 'Patient', 'Randomize', 'Others'};

% Initialize array to store mean Response_D scores for each bin
mean_raven = zeros(1, length(bin_labels));

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
    mean_raven(i) = mean(Nraven(bin(:,i)));
    se_raven(i) = std(Nraven(bin(:,i))) / sqrt(length(Nraven(bin(:,i))));
end

% Plot the bar graph
figure;

subplot(1,2,1)

bar(mean_raven);
hold on
errorbar([1 2 3 4], mean_raven, se_raven, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('Raven Scores', 'FontName', 'Times New Roman');
ylim([4 11])

% Perform t-tests and store results for p-value annotations
[h12, p12, ~] = ttest2(Nraven(bin(:,1)), Nraven(bin(:,2)));
[h23, p23, ~] = ttest2(Nraven(bin(:,2)), Nraven(bin(:,3)));
[h34, p34, ~] = ttest2(Nraven(bin(:,3)), Nraven(bin(:,4)),'Tail','right');
[h13, p13, ~] = ttest2(Nraven(bin(:,1)), Nraven(bin(:,3)));
[h14, p14, ~] = ttest2(Nraven(bin(:,1)), Nraven(bin(:,4)),'Tail','left');
[h24, p24, ~] = ttest2(Nraven(bin(:,2)), Nraven(bin(:,4)),'Tail','left');

% Define x positions for the bars
x_positions = 1:length(mean_raven);

% Y positions for each comparison, slightly above the bar heights
% y_bracket_12 = 12;  % Bin 1 vs Bin 2
% y_bracket_23 = 13.5;  % Bin 2 vs Bin 3
y_bracket_34 = 9.5;  % Bin 3 vs Bin 4
% y_bracket_13 = 16.5;  % Bin 1 vs Bin 3
y_bracket_14 = 7.5;  % Bin 1 vs Bin 4
y_bracket_24 = 8.5; % Bin 2 vs Bin 4
hold on

% First comparison (Bin 1 vs Bin 2)
%plot([x_positions(1), x_positions(1), x_positions(2), x_positions(2)], [y_bracket_12, y_bracket_12+0.25, y_bracket_12+0.25, y_bracket_12], 'k', 'LineWidth', 1.5);
%text(mean(x_positions(1:2)), y_bracket_12 + 0.75, sprintf('p = %.3f', p12), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
%hold on
% Second comparison (Bin 2 vs Bin 3)
%plot([x_positions(2), x_positions(2), x_positions(3), x_positions(3)], [y_bracket_23, y_bracket_23+0.25, y_bracket_23+0.25, y_bracket_23], 'k', 'LineWidth', 1.5);
%text(mean(x_positions(2:3)), y_bracket_23 + 0.75, sprintf('p = %.3f', p23), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
%hold on
% Third comparison (Bin 3 vs Bin 4)
plot([x_positions(3), x_positions(3), x_positions(4), x_positions(4)], [y_bracket_34, y_bracket_34+0.25, y_bracket_34+0.25, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean(x_positions(3:4)), y_bracket_34 + 0.75, sprintf('p = %.3f', p34), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Fourth comparison (Bin 1 vs Bin 3)
%plot([x_positions(1), x_positions(1), x_positions(3), x_positions(3)], [y_bracket_13, y_bracket_13+0.25, y_bracket_13+0.25, y_bracket_13], 'k', 'LineWidth', 1.5);
%text(mean([x_positions(1), x_positions(3)]), y_bracket_13 + 0.75, sprintf('p = %.3f', p13), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
%hold on
% Fifth comparison (Bin 1 vs Bin 4)
plot([x_positions(1), x_positions(1), x_positions(4), x_positions(4)], [y_bracket_14, y_bracket_14+0.25, y_bracket_14+0.25, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(4)]), y_bracket_14 + 0.75, sprintf('p = %.3f', p14), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Sixth comparison (Bin 2 vs Bin 4)
plot([x_positions(2), x_positions(2), x_positions(4), x_positions(4)], [y_bracket_24, y_bracket_24+0.25, y_bracket_24+0.25, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(4)]), y_bracket_24 + 0.75, sprintf('p = %.3f', p24), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');





% Define bin edges and labels
bin_edges = [0, 0.2, 2.2, 2.4];
bin_labels = {'Risk averse', 'Risk neutral', 'Randomize', 'Others'};

% Initialize array to store mean Response_D scores for each bin
mean_raven = zeros(1, length(bin_labels));

% Calculate mean "Response_D" score for each bin

bin(:,3) = deliberate;
bin(:,2) = ( rl & deliberate ==0 ) ;
bin(:,1) = (  ra  & deliberate ==0 ) ;
bin(:,4) = bin(:,1)+bin(:,2)+bin(:,3)==0;


bin=logical(bin);

for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    % Calculate mean "Response_D" for entries within this bin
    mean_raven(i) = mean(Nraven(bin(:,i)));
    se_raven(i) = std(Nraven(bin(:,i))) / sqrt(length(Nraven(bin(:,i))));
end

% Plot the bar graph
subplot(1,2,2)
bar(mean_raven);
hold on
errorbar([1 2 3 4], mean_raven, se_raven, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('Raven Scores', 'FontName', 'Times New Roman');
ylim([4 11])

% Perform t-tests and store results for p-value annotations
[h12, p12, ~] = ttest2(Nraven(bin(:,1)), Nraven(bin(:,2)));
[h23, p23, ~] = ttest2(Nraven(bin(:,2)), Nraven(bin(:,3)));
[h34, p34, ~] = ttest2(Nraven(bin(:,3)), Nraven(bin(:,4)),'Tail','right');
[h13, p13, ~] = ttest2(Nraven(bin(:,1)), Nraven(bin(:,3)));
[h14, p14, ~] = ttest2(Nraven(bin(:,1)), Nraven(bin(:,4)),'Tail','left');
[h24, p24, ~] = ttest2(Nraven(bin(:,2)), Nraven(bin(:,4)),'Tail','left');

% Define x positions for the bars
x_positions = 1:length(mean_raven);

% Y positions for each comparison, slightly above the bar heights
% y_bracket_12 = 12;  % Bin 1 vs Bin 2
% y_bracket_23 = 13.5;  % Bin 2 vs Bin 3
y_bracket_34 = 9.5;  % Bin 3 vs Bin 4
% y_bracket_13 = 16.5;  % Bin 1 vs Bin 3
y_bracket_14 = 7.5;  % Bin 1 vs Bin 4
y_bracket_24 = 8.5; % Bin 2 vs Bin 4


hold on

% First comparison (Bin 1 vs Bin 2)
%plot([x_positions(1), x_positions(1), x_positions(2), x_positions(2)], [y_bracket_12, y_bracket_12+0.25, y_bracket_12+0.25, y_bracket_12], 'k', 'LineWidth', 1.5);
%text(mean(x_positions(1:2)), y_bracket_12 + 0.75, sprintf('p = %.3f', p12), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
%hold on
% Second comparison (Bin 2 vs Bin 3)
%plot([x_positions(2), x_positions(2), x_positions(3), x_positions(3)], [y_bracket_23, y_bracket_23+0.25, y_bracket_23+0.25, y_bracket_23], 'k', 'LineWidth', 1.5);
%text(mean(x_positions(2:3)), y_bracket_23 + 0.75, sprintf('p = %.3f', p23), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
%hold on
% Third comparison (Bin 3 vs Bin 4)
plot([x_positions(3), x_positions(3), x_positions(4), x_positions(4)], [y_bracket_34, y_bracket_34+0.25, y_bracket_34+0.25, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean(x_positions(3:4)), y_bracket_34 + 0.75, sprintf('p = %.3f', p34), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Fourth comparison (Bin 1 vs Bin 3)
%plot([x_positions(1), x_positions(1), x_positions(3), x_positions(3)], [y_bracket_13, y_bracket_13+0.25, y_bracket_13+0.25, y_bracket_13], 'k', 'LineWidth', 1.5);
%text(mean([x_positions(1), x_positions(3)]), y_bracket_13 + 0.75, sprintf('p = %.3f', p13), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
%hold on
% Fifth comparison (Bin 1 vs Bin 4)
plot([x_positions(1), x_positions(1), x_positions(4), x_positions(4)], [y_bracket_14, y_bracket_14+0.25, y_bracket_14+0.25, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(4)]), y_bracket_14 + 0.75, sprintf('p = %.3f', p14), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Sixth comparison (Bin 2 vs Bin 4)
plot([x_positions(2), x_positions(2), x_positions(4), x_positions(4)], [y_bracket_24, y_bracket_24+0.25, y_bracket_24+0.25, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(4)]), y_bracket_24 + 0.75, sprintf('p = %.3f', p24), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');

hold off
















