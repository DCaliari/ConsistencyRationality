clear
clc

load variables.mat

Times = xlsread('DATI.xlsx','Variables','BN1:BO146');
Response_D = Times(:,1);
Response_L = Times(:,2);




%% Figure 4 for all cardinal specifications in time preferences (online appendix)

ti_all = [ti_hyp', ti_betadelta'];
tp_all = [tp_hyp', tp_betadelta'];

for z=1:size(ti_all,2)

ti = ti_all(:,z);
tp = tp_all(:,z);


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


subplot(1,size(ti_all,2),z)

% Plot the bar graph
bar(mean_times)
hold on
errorbar([1 2 3], mean_times, se_times, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('Response Times', 'FontName', 'Times New Roman');
ylim([300 550])

% Perform t-tests and store results for p-value annotations
[h12, p12, ~] = ttest2(Response_D(bin(:,1)), Response_D(bin(:,2)));
[h23, p23, ~] = ttest2(Response_D(bin(:,2)), Response_D(bin(:,3)));
[h13, p13, ~] = ttest2(Response_D(bin(:,1)), Response_D(bin(:,3)));
% Define x positions for the bars
x_positions = 1:length(mean_times); 

% Y positions for each comparison, slightly above the bar heights
y_bracket_12 = 470; % Bin 1 vs Bin 2
y_bracket_23 = 495; % Bin 2 vs Bin 3
y_bracket_13 = 520; % Bin 1 vs Bin 3
hold on
% First comparison (Bin 1 vs Bin 2)
plot([x_positions(1), x_positions(1), x_positions(2), x_positions(2)], [y_bracket_12, y_bracket_12+5, y_bracket_12+5, y_bracket_12], 'k', 'LineWidth', 1.5);
text(mean(x_positions(1:2)), y_bracket_12 + 14, sprintf('p = %.3f', p12), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Second comparison (Bin 2 vs Bin 3)
plot([x_positions(2), x_positions(2), x_positions(3), x_positions(3)], [y_bracket_23, y_bracket_23+5, y_bracket_23+5, y_bracket_23], 'k', 'LineWidth', 1.5);
text(mean(x_positions(2:3)), y_bracket_23 + 14, sprintf('p = %.3f', p23), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Third comparison (Bin 1 vs Bin 3)
plot([x_positions(1), x_positions(1), x_positions(3), x_positions(3)], [y_bracket_13, y_bracket_13+5, y_bracket_13+5, y_bracket_13], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(3)]), y_bracket_13 + 14, sprintf('p = %.3f', p13), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold off


end




%% Figure 5 for all cardinal specifications in risk preferences (online appendix)


ra_all = [ra_cara', ra_crrace', ra_carace'];
rl_all = [rl_cara', rl_crrace', rl_carace'];


for z=1:size(ra_all,2)

ra = ra_all(:,z);
rl = rl_all(:,z);


% Define bin edges and labels
bin_labels = {'Risk neutral', 'Others', 'Risk averse'};

% Initialize array to store mean Raven scores for each bin
mean_times = zeros(1, length(bin_labels));
se_times = zeros(1, length(bin_labels));

bin(:,2) = ( ~ra&~rl );
bin(:,1) = ( rl ) ;
bin(:,3) = (  ra ) ;
% Calculate mean "Raven" score for each bin
for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    % Calculate mean "Response_D" for entries within this bin
    mean_times(i) = mean(Response_L(bin(:,i)));
    se_times(i) = std(Response_L(bin(:,i))) / sqrt(length(Response_L(bin(:,i))));
end

subplot(1,size(ra_all,2),z)

% Plot the bar graph
bar(mean_times)
hold on
errorbar([1 2 3], mean_times, se_times, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('Response times', 'FontName', 'Times New Roman');
ylim([400 1100])


% Perform t-tests and store results for p-value annotations
[h12, p12, ~] = ttest2(Response_L(bin(:,1)), Response_L(bin(:,2)));
[h23, p23, ~] = ttest2(Response_L(bin(:,2)), Response_L(bin(:,3)));
[h13, p13, ~] = ttest2(Response_L(bin(:,1)), Response_L(bin(:,3)));

% Define x positions for the bars
x_positions = 1:length(mean_times); 

% Y positions for each comparison, slightly above the bar heights
y_bracket_12 = 840; % Bin 1 vs Bin 2
y_bracket_23 = 920; % Bin 2 vs Bin 3
y_bracket_13 = 1000; % Bin 1 vs Bin 3
hold on
% First comparison (Bin 1 vs Bin 2)
plot([x_positions(1), x_positions(1), x_positions(2), x_positions(2)], [y_bracket_12, y_bracket_12+7.5, y_bracket_12+7.5, y_bracket_12], 'k', 'LineWidth', 1.5);
text(mean(x_positions(1:2)), y_bracket_12 + 40, sprintf('p = %.3f', p12), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Second comparison (Bin 2 vs Bin 3)
plot([x_positions(2), x_positions(2), x_positions(3), x_positions(3)], [y_bracket_23, y_bracket_23+7.5, y_bracket_23+7.5, y_bracket_23], 'k', 'LineWidth', 1.5);
text(mean(x_positions(2:3)), y_bracket_23 + 40, sprintf('p = %.3f', p23), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Third comparison (Bin 1 vs Bin 3)
plot([x_positions(1), x_positions(1), x_positions(3), x_positions(3)], [y_bracket_13, y_bracket_13+7.5, y_bracket_13+7.5, y_bracket_13], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(3)]), y_bracket_13 + 40, sprintf('p = %.3f', p13), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold off

end



%% Other deciles for extreme preferences (online appendix)


ti_exp2 = time_ind_exp>= 0.08/0.92;
ti_exp3 = time_ind_exp>= 0.07/0.93;

tp_exp2 = time_ind_exp<= 0.02/0.98;
tp_exp3 = time_ind_exp<= 0.03/0.97;

ra_crra2 = risk_ind_crra>=2.0;
ra_crra3 = risk_ind_crra>=1.8;

rl_crra2 = risk_ind_crra<=0.4;
rl_crra3 = risk_ind_crra<=0.6;


ti_n = [ti_exp2', ti_exp3'];
tp_n = [tp_exp2', tp_exp3'];
ra_n = [ra_crra2', ra_crra3'];
rl_n = [rl_crra2', rl_crra3'];




for z=1:size(ti_n,2)

ti = ti_n(:,z);
tp = tp_n(:,z);


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


subplot(1,size(ti_n,2),z)

% Plot the bar graph
bar(mean_times)
hold on
errorbar([1 2 3], mean_times, se_times, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('Response Times', 'FontName', 'Times New Roman');
ylim([300 530])

% Perform t-tests and store results for p-value annotations
[h12, p12, ~] = ttest2(Response_D(bin(:,1)), Response_D(bin(:,2)));
[h23, p23, ~] = ttest2(Response_D(bin(:,2)), Response_D(bin(:,3)));
[h13, p13, ~] = ttest2(Response_D(bin(:,1)), Response_D(bin(:,3)));
% Define x positions for the bars
x_positions = 1:length(mean_times); 

% Y positions for each comparison, slightly above the bar heights
y_bracket_12 = 470; % Bin 1 vs Bin 2
y_bracket_23 = 490; % Bin 2 vs Bin 3
y_bracket_13 = 510; % Bin 1 vs Bin 3
hold on
% First comparison (Bin 1 vs Bin 2)
plot([x_positions(1), x_positions(1), x_positions(2), x_positions(2)], [y_bracket_12, y_bracket_12+5, y_bracket_12+5, y_bracket_12], 'k', 'LineWidth', 1.5);
text(mean(x_positions(1:2)), y_bracket_12 + 12, sprintf('p = %.3f', p12), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Second comparison (Bin 2 vs Bin 3)
plot([x_positions(2), x_positions(2), x_positions(3), x_positions(3)], [y_bracket_23, y_bracket_23+5, y_bracket_23+5, y_bracket_23], 'k', 'LineWidth', 1.5);
text(mean(x_positions(2:3)), y_bracket_23 + 12, sprintf('p = %.3f', p23), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Third comparison (Bin 1 vs Bin 3)
plot([x_positions(1), x_positions(1), x_positions(3), x_positions(3)], [y_bracket_13, y_bracket_13+5, y_bracket_13+5, y_bracket_13], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(3)]), y_bracket_13 + 12, sprintf('p = %.3f', p13), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold off


end






for z=1:size(ra_n,2)

ra = ra_n(:,z);
rl = rl_n(:,z);


% Define bin edges and labels
bin_labels = {'Risk neutral', 'Others', 'Risk averse'};

% Initialize array to store mean Raven scores for each bin
mean_times = zeros(1, length(bin_labels));
se_times = zeros(1, length(bin_labels));

bin(:,2) = ( ~ra&~rl );
bin(:,1) = ( rl ) ;
bin(:,3) = (  ra ) ;
% Calculate mean "Raven" score for each bin
for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    % Calculate mean "Response_D" for entries within this bin
    mean_times(i) = mean(Response_L(bin(:,i)));
    se_times(i) = std(Response_L(bin(:,i))) / sqrt(length(Response_L(bin(:,i))));
end

subplot(1,size(ra_n,2),z)

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
y_bracket_12 = 820; % Bin 1 vs Bin 2
y_bracket_23 = 870; % Bin 2 vs Bin 3
y_bracket_13 = 920; % Bin 1 vs Bin 3
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

end

