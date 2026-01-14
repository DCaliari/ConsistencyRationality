clear
clc


load variables.mat

load warpD.mat

WARP_D = WARPxyzw;

Impatience = xlsread('DATI.xlsx','Variables','AN1:AN146');
Monotonicity = xlsread('DATI.xlsx','Variables','AO1:AO146');


FOSD = xlsread('DATI.xlsx','Variables','AP1:AQ146');
FOSD = sum(FOSD,2)>0;
SOSD = xlsread('DATI.xlsx','Variables','AS1:AS146');


load warpL.mat

WARP_L = WARPabcd;



%% For all other cardinal specifications (online appendix)


ti_all = [ti_hyp', ti_betadelta'];
tp_all = [tp_hyp', tp_betadelta'];

for z=1:size(ti_all,2)

ti = ti_all(:,z);
tp = tp_all(:,z);


% Define bin edges and labels
bin_labels = {'Impatient', 'Others', 'Others & consistent', 'Patient'};


% Calculate mean "Response_D" score for each bin

clearvars bin

bin(:,2) = ( ~ti&~tp );
bin(:,1) = ( ti ) ;
bin(:,4) = (  tp ) ;

bin = logical(bin);

bin(:,3) = ( ~ti&~tp & WARP_D <= 1 );


for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    % Calculate mean "Response_D" for entries within this bin
    mean_impatience(i) = mean(Impatience(bin(:,i)));
    se_impatience(i) = std(Impatience(bin(:,i))) / sqrt(length(Impatience(bin(:,i))));
end



% Calculate the mean WARP_D for each bin
mean_warp_d = zeros(1, length(bin_labels));
for i = 1:length(bin_labels)
    mean_warp_d(i) = mean(WARP_D(bin(:, i)));
end


subplot(1,size(ti_all,2),z)


yyaxis left; % Use left y-axis for mean impatience
bar(mean_impatience)
hold on
errorbar([1 2 3 4], mean_impatience, se_impatience, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman');
ylabel('Share of DMs violating Impatience', 'FontName', 'Times New Roman');
ylim([0 1]); yticks(0:0.2:1);

% Plot mean WARP_D on the right axis
yyaxis right;
plot(1:length(mean_warp_d), mean_warp_d, '-o', 'LineWidth', 1, ...
    'MarkerSize', 6, 'Color', [0.8 0.2 0.2]);
ylabel('Mean WARP violations', 'FontName', 'Times New Roman');
ylim([min(mean_warp_d) - 0.1, max(mean_warp_d) + 1]); % Adjust y-limits dynamically

% Enhance the plot appearance
set(gca, 'XTick', 1:length(bin_labels)); % Align ticks with bin positions
legend({'Mean Impatience', 'Mean WARP violations'}, 'FontName', 'Times New Roman', ...
       'Location', 'Best', 'Box', 'off');



x = [sum(Impatience(bin(:,1))), sum(~Impatience(bin(:,1)));
    sum(Impatience(bin(:,2) )), sum(~Impatience(bin(:,2) ))]

[h12,p12,stats] = fishertest(x)

x = [sum(Impatience(bin(:,1))), sum(~Impatience(bin(:,1)));
    sum(Impatience(bin(:,3) )), sum(~Impatience(bin(:,3) ))]
[h13,p13,stats] = fishertest(x)

x = [sum(Impatience(bin(:,2))), sum(~Impatience(bin(:,2)));
    sum(Impatience(bin(:,4) )), sum(~Impatience(bin(:,4) ))]
[h24,p24,stats] = fishertest(x)

x = [sum(Impatience(bin(:,3))), sum(~Impatience(bin(:,3)));
    sum(Impatience(bin(:,4) )), sum(~Impatience(bin(:,4) ))]
[h34,p34,stats] = fishertest(x)



% Define x positions for the bars
x_positions = 1:length(mean_impatience);

% Get current limits of the left y-axis
left_ylim = ylim; % Save the left y-axis limits

% Define Y positions for each comparison as fractions of the left y-axis range
y_bracket_fractions = [0.6, 0.7, 0.8, 0.9]; % Relative heights (fractions of the range)
y_brackets = left_ylim(1) + (left_ylim(2) - left_ylim(1)) * y_bracket_fractions; % Convert to actual y-values

% Define x positions for the bars
x_positions = 1:length(mean_impatience);

% P-values and comparisons
p_values = [p12, p13, p24, p34]; % Corresponding p-values
comparisons = [1 2; 1 3; 2 4; 3 4]; % Pairs of bins being compared

hold on;

% Loop through each comparison
for i = 1:size(comparisons, 1)
    % Get x positions and y bracket height for the comparison
    x_pair = x_positions(comparisons(i, :));
    y_bracket = y_brackets(i); % Base y-position for the bracket

    % Define the height of the bracket consistently
    bracket_height = 0.02 * diff(left_ylim); % Fixed offset for height

    % Plot the comparison bracket using continuous lines
    plot([x_pair(1), x_pair(1)], [y_bracket, y_bracket + bracket_height], 'k-', 'LineWidth', 1); % Left vertical line
    plot([x_pair(1), x_pair(2)], [y_bracket + bracket_height, y_bracket + bracket_height], 'k-', 'LineWidth', 1); % Horizontal line
    plot([x_pair(2), x_pair(2)], [y_bracket + bracket_height, y_bracket], 'k-', 'LineWidth', 1); % Right vertical line

    % Add the p-value annotation
    text(mean(x_pair), y_bracket + bracket_height + 0.03 * diff(left_ylim), ...
         sprintf('p = %.3f', p_values(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontName', 'Times New Roman');
end


legend('off')

hold off;

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


ti_all = [ti_exp2', ti_exp3'];
tp_all = [tp_exp2', tp_exp3'];

for z=1:size(ti_all,2)

ti = ti_all(:,z);
tp = tp_all(:,z);


% Define bin edges and labels
bin_labels = {'Impatient', 'Others', 'Others & consistent', 'Patient'};


% Calculate mean "Response_D" score for each bin

clearvars bin

bin(:,2) = ( ~ti&~tp );
bin(:,1) = ( ti ) ;
bin(:,4) = (  tp ) ;

bin = logical(bin);

bin(:,3) = ( ~ti&~tp & WARP_D <= 1 );


for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    % Calculate mean "Response_D" for entries within this bin
    mean_impatience(i) = mean(Impatience(bin(:,i)));
    se_impatience(i) = std(Impatience(bin(:,i))) / sqrt(length(Impatience(bin(:,i))));
end



% Calculate the mean WARP_D for each bin
mean_warp_d = zeros(1, length(bin_labels));
for i = 1:length(bin_labels)
    mean_warp_d(i) = mean(WARP_D(bin(:, i)));
end


subplot(1,size(ti_all,2),z)


yyaxis left; % Use left y-axis for mean impatience
bar(mean_impatience)
hold on
errorbar([1 2 3 4], mean_impatience, se_impatience, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman');
ylabel('Share of DMs violating Impatience', 'FontName', 'Times New Roman');
ylim([0 1]); yticks(0:0.2:1);

% Plot mean WARP_D on the right axis
yyaxis right;
plot(1:length(mean_warp_d), mean_warp_d, '-o', 'LineWidth', 1, ...
    'MarkerSize', 6, 'Color', [0.8 0.2 0.2]);
ylabel('Mean WARP violations', 'FontName', 'Times New Roman');
ylim([min(mean_warp_d) - 0.1, max(mean_warp_d) + 1]); % Adjust y-limits dynamically

% Enhance the plot appearance
set(gca, 'XTick', 1:length(bin_labels)); % Align ticks with bin positions
legend({'Mean Impatience', 'Mean WARP violations'}, 'FontName', 'Times New Roman', ...
       'Location', 'Best', 'Box', 'off');



x = [sum(Impatience(bin(:,1))), sum(~Impatience(bin(:,1)));
    sum(Impatience(bin(:,2) )), sum(~Impatience(bin(:,2) ))]

[h12,p12,stats] = fishertest(x)

x = [sum(Impatience(bin(:,1))), sum(~Impatience(bin(:,1)));
    sum(Impatience(bin(:,3) )), sum(~Impatience(bin(:,3) ))]
[h13,p13,stats] = fishertest(x)

x = [sum(Impatience(bin(:,2))), sum(~Impatience(bin(:,2)));
    sum(Impatience(bin(:,4) )), sum(~Impatience(bin(:,4) ))]
[h24,p24,stats] = fishertest(x)

x = [sum(Impatience(bin(:,3))), sum(~Impatience(bin(:,3)));
    sum(Impatience(bin(:,4) )), sum(~Impatience(bin(:,4) ))]
[h34,p34,stats] = fishertest(x)



% Define x positions for the bars
x_positions = 1:length(mean_impatience);

% Get current limits of the left y-axis
left_ylim = ylim; % Save the left y-axis limits

% Define Y positions for each comparison as fractions of the left y-axis range
y_bracket_fractions = [0.6, 0.7, 0.8, 0.9]; % Relative heights (fractions of the range)
y_brackets = left_ylim(1) + (left_ylim(2) - left_ylim(1)) * y_bracket_fractions; % Convert to actual y-values

% Define x positions for the bars
x_positions = 1:length(mean_impatience);

% P-values and comparisons
p_values = [p12, p13, p24, p34]; % Corresponding p-values
comparisons = [1 2; 1 3; 2 4; 3 4]; % Pairs of bins being compared

hold on;

% Loop through each comparison
for i = 1:size(comparisons, 1)
    % Get x positions and y bracket height for the comparison
    x_pair = x_positions(comparisons(i, :));
    y_bracket = y_brackets(i); % Base y-position for the bracket

    % Define the height of the bracket consistently
    bracket_height = 0.02 * diff(left_ylim); % Fixed offset for height

    % Plot the comparison bracket using continuous lines
    plot([x_pair(1), x_pair(1)], [y_bracket, y_bracket + bracket_height], 'k-', 'LineWidth', 1); % Left vertical line
    plot([x_pair(1), x_pair(2)], [y_bracket + bracket_height, y_bracket + bracket_height], 'k-', 'LineWidth', 1); % Horizontal line
    plot([x_pair(2), x_pair(2)], [y_bracket + bracket_height, y_bracket], 'k-', 'LineWidth', 1); % Right vertical line

    % Add the p-value annotation
    text(mean(x_pair), y_bracket + bracket_height + 0.03 * diff(left_ylim), ...
         sprintf('p = %.3f', p_values(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontName', 'Times New Roman');
end


legend('off')

hold off;

end








%% For all other cardinal specifications (online appendix)



ra_all = [ra_cara', ra_crrace', ra_carace'];
rl_all = [rl_cara', rl_crrace', rl_carace'];

for z=1:size(ra_all, 2)

clearvars bin_labels
ra = ra_all(:,z);
rl = rl_all(:,z);


bin_labels = {'Risk neutral', 'Others', 'Others & consistent' , 'Risk averse'};


% Calculate mean "Response_D" score for each bin

bin(:,2) = ( ~rl&~ra );
bin(:,1) = ( rl ) ;
bin(:,4) = ( ra ) ;

bin = logical(bin);

bin(:,3) = (~rl&~ra & WARP_L <= 7 );


mean_SOSD = zeros(1, length(bin_labels));


for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    % Calculate mean "Response_D" for entries within this bin
    mean_SOSD(i) = mean(SOSD(bin(:,i)));
    se_SOSD(i) = std(SOSD(bin(:,i))) / sqrt(length(SOSD(bin(:,i))));
end


% Calculate the mean WARP_D for each bin
mean_warp_d = zeros(1, length(bin_labels));
for i = 1:length(bin_labels)
    mean_warp_d(i) = mean(WARP_L(bin(:, i)));
end

% Plot the bar graph for mean SOSD
subplot(1,size(ra_all, 2),z)
yyaxis left; % Use left y-axis for mean SOSD
bar(mean_SOSD);
hold on
errorbar([1 2 3 4], mean_SOSD, se_SOSD, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman');
ylabel('Share of DMs violating SOSD', 'FontName', 'Times New Roman');
ylim([0 1.8]); yticks(0:0.2:1);

% Plot mean WARP_D on the right axis
yyaxis right;
plot(1:length(mean_warp_d), mean_warp_d, '-o', 'LineWidth', 1, ...
    'MarkerSize', 6, 'Color', [0.8 0.2 0.2]);
ylabel('Mean WARP violations', 'FontName', 'Times New Roman');
ylim([min(mean_warp_d) - 1, max(mean_warp_d) + 3]); % Adjust y-limits dynamically

% Enhance the plot appearance
set(gca, 'XTick', 1:length(bin_labels)); % Align ticks with bin positions
legend({'Mean SOSD', 'Mean WARP violations'}, 'FontName', 'Times New Roman', ...
       'Location', 'Best', 'Box', 'off');



x = [sum(SOSD(bin(:,1))), sum(~SOSD(bin(:,1)));
    sum(SOSD(bin(:,2) )), sum(~SOSD(bin(:,2) ))]

[h12,p12,stats] = fishertest(x)

x = [sum(SOSD(bin(:,1))), sum(~SOSD(bin(:,1)));
    sum(SOSD(bin(:,3) )), sum(~SOSD(bin(:,3) ))]
[h13,p13,stats] = fishertest(x)

x = [sum(SOSD(bin(:,2))), sum(~SOSD(bin(:,2)));
    sum(SOSD(bin(:,4) )), sum(~SOSD(bin(:,4) ))]
[h24,p24,stats] = fishertest(x)

x = [sum(SOSD(bin(:,3))), sum(~SOSD(bin(:,3)));
    sum(SOSD(bin(:,4) )), sum(~SOSD(bin(:,4) ))]
[h34,p34,stats] = fishertest(x)



% Define x positions for the bars
x_positions = 1:length(mean_SOSD);

% Get current limits of the left y-axis
left_ylim = ylim; % Save the left y-axis limits

% Define Y positions for each comparison as fractions of the left y-axis range
y_bracket_fractions = [0.6, 0.7, 0.8, 0.9]; % Relative heights (fractions of the range)
y_brackets = left_ylim(1) + (left_ylim(2) - left_ylim(1)) * y_bracket_fractions; % Convert to actual y-values

% Define x positions for the bars
x_positions = 1:length(mean_SOSD);

% P-values and comparisons
p_values = [p12, p13, p24, p34]; % Corresponding p-values
comparisons = [1 2; 1 3; 2 4; 3 4]; % Pairs of bins being compared

hold on;

% Loop through each comparison
for i = 1:size(comparisons, 1)
    % Get x positions and y bracket height for the comparison
    x_pair = x_positions(comparisons(i, :));
    y_bracket = y_brackets(i); % Base y-position for the bracket

    % Define the height of the bracket consistently
    bracket_height = 0.02 * diff(left_ylim); % Fixed offset for height

    % Plot the comparison bracket using continuous lines
    plot([x_pair(1), x_pair(1)], [y_bracket, y_bracket + bracket_height], 'k-', 'LineWidth', 1); % Left vertical line
    plot([x_pair(1), x_pair(2)], [y_bracket + bracket_height, y_bracket + bracket_height], 'k-', 'LineWidth', 1); % Horizontal line
    plot([x_pair(2), x_pair(2)], [y_bracket + bracket_height, y_bracket], 'k-', 'LineWidth', 1); % Right vertical line

    % Add the p-value annotation
    text(mean(x_pair), y_bracket + bracket_height + 0.03 * diff(left_ylim), ...
         sprintf('p = %.3f', p_values(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontName', 'Times New Roman');
end


legend('off')

hold off;


end


%% other extreme preferences definitions (online appendix)


ra_all = [ra_crra2', ra_crra3'];
rl_all = [rl_crra2', rl_crra3'];

for z=1:size(ra_all, 2)

clearvars bin_labels
ra = ra_all(:,z);
rl = rl_all(:,z);


bin_labels = {'Risk neutral', 'Others', 'Others & consistent' , 'Risk averse'};


% Calculate mean "Response_D" score for each bin

bin(:,2) = ( ~rl&~ra );
bin(:,1) = ( rl ) ;
bin(:,4) = ( ra ) ;

bin = logical(bin);

bin(:,3) = (~rl&~ra & WARP_L <= 7 );


mean_SOSD = zeros(1, length(bin_labels));


for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    % Calculate mean "Response_D" for entries within this bin
    mean_SOSD(i) = mean(SOSD(bin(:,i)));
    se_SOSD(i) = std(SOSD(bin(:,i))) / sqrt(length(SOSD(bin(:,i))));
end


% Calculate the mean WARP_D for each bin
mean_warp_d = zeros(1, length(bin_labels));
for i = 1:length(bin_labels)
    mean_warp_d(i) = mean(WARP_L(bin(:, i)));
end

% Plot the bar graph for mean SOSD
subplot(1,size(ra_all, 2),z)
yyaxis left; % Use left y-axis for mean SOSD
bar(mean_SOSD);
hold on
errorbar([1 2 3 4], mean_SOSD, se_SOSD, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman');
ylabel('Share of DMs violating SOSD', 'FontName', 'Times New Roman');
ylim([0 1.5]); yticks(0:0.2:1);

% Plot mean WARP_D on the right axis
yyaxis right;
plot(1:length(mean_warp_d), mean_warp_d, '-o', 'LineWidth', 1, ...
    'MarkerSize', 6, 'Color', [0.8 0.2 0.2]);
ylabel('Mean WARP violations', 'FontName', 'Times New Roman');
ylim([min(mean_warp_d) - 1, max(mean_warp_d) + 1]); % Adjust y-limits dynamically

% Enhance the plot appearance
set(gca, 'XTick', 1:length(bin_labels)); % Align ticks with bin positions
legend({'Mean SOSD', 'Mean WARP violations'}, 'FontName', 'Times New Roman', ...
       'Location', 'Best', 'Box', 'off');



x = [sum(SOSD(bin(:,1))), sum(~SOSD(bin(:,1)));
    sum(SOSD(bin(:,2) )), sum(~SOSD(bin(:,2) ))]

[h12,p12,stats] = fishertest(x)

x = [sum(SOSD(bin(:,1))), sum(~SOSD(bin(:,1)));
    sum(SOSD(bin(:,3) )), sum(~SOSD(bin(:,3) ))]
[h13,p13,stats] = fishertest(x)

x = [sum(SOSD(bin(:,2))), sum(~SOSD(bin(:,2)));
    sum(SOSD(bin(:,4) )), sum(~SOSD(bin(:,4) ))]
[h24,p24,stats] = fishertest(x)

x = [sum(SOSD(bin(:,3))), sum(~SOSD(bin(:,3)));
    sum(SOSD(bin(:,4) )), sum(~SOSD(bin(:,4) ))]
[h34,p34,stats] = fishertest(x)



% Define x positions for the bars
x_positions = 1:length(mean_SOSD);

% Get current limits of the left y-axis
left_ylim = ylim; % Save the left y-axis limits

% Define Y positions for each comparison as fractions of the left y-axis range
y_bracket_fractions = [0.6, 0.7, 0.8, 0.9]; % Relative heights (fractions of the range)
y_brackets = left_ylim(1) + (left_ylim(2) - left_ylim(1)) * y_bracket_fractions; % Convert to actual y-values

% Define x positions for the bars
x_positions = 1:length(mean_SOSD);

% P-values and comparisons
p_values = [p12, p13, p24, p34]; % Corresponding p-values
comparisons = [1 2; 1 3; 2 4; 3 4]; % Pairs of bins being compared

hold on;

% Loop through each comparison
for i = 1:size(comparisons, 1)
    % Get x positions and y bracket height for the comparison
    x_pair = x_positions(comparisons(i, :));
    y_bracket = y_brackets(i); % Base y-position for the bracket

    % Define the height of the bracket consistently
    bracket_height = 0.02 * diff(left_ylim); % Fixed offset for height

    % Plot the comparison bracket using continuous lines
    plot([x_pair(1), x_pair(1)], [y_bracket, y_bracket + bracket_height], 'k-', 'LineWidth', 1); % Left vertical line
    plot([x_pair(1), x_pair(2)], [y_bracket + bracket_height, y_bracket + bracket_height], 'k-', 'LineWidth', 1); % Horizontal line
    plot([x_pair(2), x_pair(2)], [y_bracket + bracket_height, y_bracket], 'k-', 'LineWidth', 1); % Right vertical line

    % Add the p-value annotation
    text(mean(x_pair), y_bracket + bracket_height + 0.03 * diff(left_ylim), ...
         sprintf('p = %.3f', p_values(i)), ...
         'HorizontalAlignment', 'center', 'FontSize', 11, 'FontName', 'Times New Roman');
end


legend('off')

hold off;


end


















