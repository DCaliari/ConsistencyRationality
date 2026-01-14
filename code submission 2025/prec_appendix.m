

prec_ind_time = (prec_ind_exp - min(prec_ind_exp))./(max(prec_ind_exp)- min(prec_ind_exp));


ti = ti_exp';
tp = tp_exp';



% Define bin edges and labels
bin_edges = [0.9, 0.91, 0.99, 1.0];
bin_labels = {'Impatient', 'Patient', 'Randomize', 'Others'};

% Initialize array to store mean Response_D scores for each bin
mean_prec = zeros(1, length(bin_labels));

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
    mean_prec(i) = mean(prec_ind_time(bin(:,i)));
    se_prec(i) = std(prec_ind_time(bin(:,i))) / sqrt(length(prec_ind_time(bin(:,i))));
end

% Plot the bar graph

subplot(1,2,1)

bar(mean_prec);
hold on
errorbar([1 2 3 4], mean_prec, se_prec, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('Precision', 'FontName', 'Times New Roman');
ylim([0 1.2])
yticks([0:0.2:1])

% Perform t-tests and store results for p-value annotations

[~, p34, ~] = ttest2(prec_ind_time(bin(:,3)), prec_ind_time(bin(:,4)),'Tail','left');
[~, p14, ~] = ttest2(prec_ind_time(bin(:,1)), prec_ind_time(bin(:,4)),'Tail','right');
[~, p24, ~] = ttest2(prec_ind_time(bin(:,2)), prec_ind_time(bin(:,4)),'Tail','right');


% Define x positions for the bars
x_positions = 1:length(mean_prec);


% Y positions for each comparison, slightly above the bar heights
y_bracket_34 = 1.1;  % Bin 3 vs Bin 4
y_bracket_14 = 0.9;  % Bin 1 vs Bin 4
y_bracket_24 = 1; % Bin 2 vs Bin 4

hold on

plot([x_positions(1), x_positions(1), x_positions(4), x_positions(4)], [y_bracket_14, y_bracket_14+0.015, y_bracket_14+0.015, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(4)]), y_bracket_14 + 0.05, sprintf('p = %.3f', p14), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Sixth comparison (Bin 2 vs Bin 4)
plot([x_positions(2), x_positions(2), x_positions(4), x_positions(4)], [y_bracket_24, y_bracket_24+0.015, y_bracket_24+0.015, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(4)]), y_bracket_24 + 0.05, sprintf('p = %.3f', p24), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Third comparison (Bin 3 vs Bin 4)
plot([x_positions(3), x_positions(3), x_positions(4), x_positions(4)], [y_bracket_34, y_bracket_34+0.015, y_bracket_34+0.015, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean(x_positions(3:4)), y_bracket_34 + 0.05, sprintf('p = %.3f', p34), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on





ra = ra_crra';
rl = rl_crra';



prec_ind_risk = (prec_ind_crra - min(prec_ind_crra))./(max(prec_ind_crra)- min(prec_ind_crra));

% Define bin edges and labels
bin_edges = [0, 0.2, 2.2, 2.4];
bin_labels = {'Risk averse', 'Risk neutral', 'Randomize', 'Others'};

% Initialize array to store mean Response_D scores for each bin
mean_warp = zeros(1, length(bin_labels));

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
    mean_prec(i) = mean(prec_ind_risk(bin(:,i)));
    se_prec(i) = std(prec_ind_risk(bin(:,i))) / sqrt(length(prec_ind_risk(bin(:,i))));
end

% Plot the bar graph

subplot(1,2,2)

bar(mean_prec);
hold on
errorbar([1 2 3 4], mean_prec, se_prec, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('Precision', 'FontName', 'Times New Roman');
ylim([0 1.2])
yticks([0:0.2:1])

[~, p34, ~] = ttest2(prec_ind_risk(bin(:,3)), prec_ind_risk(bin(:,4)),'Tail','left');
[~, p14, ~] = ttest2(prec_ind_risk(bin(:,1)), prec_ind_risk(bin(:,4)),'Tail','right');
[~, p24, ~] = ttest2(prec_ind_risk(bin(:,2)), prec_ind_risk(bin(:,4)),'Tail','right');


% Define x positions for the bars
x_positions = 1:length(mean_prec);


% Y positions for each comparison, slightly above the bar heights
y_bracket_34 = 1.1;  % Bin 3 vs Bin 4
y_bracket_14 = 0.9;  % Bin 1 vs Bin 4
y_bracket_24 = 1; % Bin 2 vs Bin 4

hold on

plot([x_positions(1), x_positions(1), x_positions(4), x_positions(4)], [y_bracket_14, y_bracket_14+0.015, y_bracket_14+0.015, y_bracket_14], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(4)]), y_bracket_14 + 0.05, sprintf('p = %.3f', p14), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Sixth comparison (Bin 2 vs Bin 4)
plot([x_positions(2), x_positions(2), x_positions(4), x_positions(4)], [y_bracket_24, y_bracket_24+0.015, y_bracket_24+0.015, y_bracket_24], 'k', 'LineWidth', 1.5);
text(mean([x_positions(2), x_positions(4)]), y_bracket_24 + 0.05, sprintf('p = %.3f', p24), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Third comparison (Bin 3 vs Bin 4)
plot([x_positions(3), x_positions(3), x_positions(4), x_positions(4)], [y_bracket_34, y_bracket_34+0.015, y_bracket_34+0.015, y_bracket_34], 'k', 'LineWidth', 1.5);
text(mean(x_positions(3:4)), y_bracket_34 + 0.05, sprintf('p = %.3f', p34), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');


hold off







