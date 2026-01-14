%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TIME Preferences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc



loadtimedata
Sets=Sets==0;

%% EXPONENTIAL PREFERENCES

%% now we set the estimated values of utilities v at the relevant types (TIME)

X1 = [160 0 0 0 0 ];
X2 = [110 50 25 0 0 ];
X3 = [50 50 50 50 0 ];
X4 = [0 15 40 170 0];

X = [X1; X2; X3; X4];

T = [0 3 6 9 12];

T = [T;T;T;T];

clearvars X1 X2 X3 X4

beta = linspace(0.2, 1, 10000);


[~,PERexp] = exptypes(X,T,beta,Na);




%% Logit estimation of utilities v in Time

options=optimoptions('fmincon','display','off','MaxFunctionEvaluations',1e7,...
    'StepTolerance',1e-9,'FunctionTolerance',1e-9,'OptimalityTolerance',1e-9,'FiniteDifferenceType','central' );


for ind = 1:N_subjects
x0 = ones(1,size(PERexp,1));
x0= x0./(sum(x0));

Aeq = ones(1,size(PERexp,1));
beq = 1;
lb = zeros(1,size(PERexp,1));
ub = ones(1,size(PERexp,1));

% [utilities_risk(:,ind), Lrisk(ind,1)] = fmincon(@(x) loglikelihood_logit(x, Data, Na, Main_sets, Sets, ind), x0, [], [], Aeq, beq, lb, ub, [], options);
[probabilities_ordinal_time(:,ind), Ltime_ordinal(ind,1)] = fmincon(@(x) loglikelihood_ordinal(x, PERexp, Data, Na, Main_sets, Sets, ind), x0, [], [], Aeq, beq, lb, ub, [], options);
end

ti_ord = (probabilities_ordinal_time(1,:)>=0.5)';
tp_ord = (probabilities_ordinal_time(7,:)>=0.5)';


DataTime=Data;
clearvars -except probabilities_ordinal_time Na Main_sets Sets N_subjects Ltime_ordinal DataTime ti_ord tp_ord


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RISK Preferences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



loadriskdata
Sets=Sets==0;

%% CRRA PREFERENCES

PER = flipud(perms(1:4)); % ordinal preferences

% select the ordinal preferences that are in line with CRRA utility
% function

crra = linspace(-3, 4, 10000);

X1 = [50 0.00001];
X2 = [65 25];
X3 = [90 25];
X4 = [300 5];

X = [X1;X2;X3;X4];

P1 = [1 0];
P2 = [0.8 0.2];
P3 = [0.5 0.5];
P4 = [0.2 0.8];

P = [P1;P2;P3;P4];

clearvars X1 X2 X3 X4 P1 P2 P3 P4

[~,PERcrra] = crratypes(X,P,crra,Na);



%% CARA PREFERENCES

% select the ordinal preferences that are in line with CRRA utility
% function

cara = linspace(-3, 4, 10000);

X1 = [50 0.00001];
X2 = [65 25];
X3 = [90 25];
X4 = [300 5];

X = [X1;X2;X3;X4];

P1 = [1 0];
P2 = [0.8 0.2];
P3 = [0.5 0.5];
P4 = [0.2 0.8];

P = [P1;P2;P3;P4];

clearvars X1 X2 X3 X4 P1 P2 P3 P4

[~,PERcara] = caratypes(X,P,cara,Na);



%% Estimation


options=optimoptions('fmincon','display','off','MaxFunctionEvaluations',1e7,...
    'StepTolerance',1e-9,'FunctionTolerance',1e-9,'OptimalityTolerance',1e-9,'FiniteDifferenceType','central' );

for ind = 1:N_subjects
x0 = ones(1,size(PERcrra,1));
x0= x0./(sum(x0));

Aeq = ones(1,size(PERcrra,1));
beq = 1;
lb = zeros(1,size(PERcrra,1));
ub = ones(1,size(PERcrra,1));

% [utilities_risk(:,ind), Lrisk(ind,1)] = fmincon(@(x) loglikelihood_logit(x, Data, Na, Main_sets, Sets, ind), x0, [], [], Aeq, beq, lb, ub, [], options);
[probabilities_ordinal(:,ind), Lrisk_ordinal(ind,1)] = fmincon(@(x) loglikelihood_ordinal(x, PERcrra, Data, Na, Main_sets, Sets, ind), x0, [], [], Aeq, beq, lb, ub, [], options);
end

ra_ord = (probabilities_ordinal(1,:)>=0.5)';
rl_ord = (probabilities_ordinal(7,:)>=0.5)';

DataRisk=Data;


% ESTIMATION CARA

options=optimoptions('fmincon','display','off','MaxFunctionEvaluations',1e7,...
    'StepTolerance',1e-9,'FunctionTolerance',1e-9,'OptimalityTolerance',1e-9,'FiniteDifferenceType','central' );

for ind = 1:N_subjects
x0 = ones(1,size(PERcara,1));
x0= x0./(sum(x0));

Aeq = ones(1,size(PERcara,1));
beq = 1;
lb = zeros(1,size(PERcara,1));
ub = ones(1,size(PERcara,1));

% [utilities_risk(:,ind), Lrisk(ind,1)] = fmincon(@(x) loglikelihood_logit(x, Data, Na, Main_sets, Sets, ind), x0, [], [], Aeq, beq, lb, ub, [], options);
[probabilities_ordinal(:,ind), Lrisk_ordinal(ind,1)] = fmincon(@(x) loglikelihood_ordinal(x, PERcara, Data, Na, Main_sets, Sets, ind), x0, [], [], Aeq, beq, lb, ub, [], options);
end

ra_ord_cara = (probabilities_ordinal(1,:)>=0.5)';
rl_ord_cara = (probabilities_ordinal(7,:)>=0.5)';

DataRisk=Data;

clearvars -except probabilities_ordinal_time probabilities_ordinal Na Main_sets Sets N_subjects Ltime_ordinal Lrisk_ordinal DataTime ti_ord tp_ord ra_ord rl_ord ra_ord_cara rl_ord_cara


%% CORRELATION  

load variables.mat
ti_exp = time_ind_exp>= 0.09/0.91;
tp_exp = time_ind_exp<= 0.01/0.99;

ra_crra = risk_ind_crra>=2.2;
rl_crra = risk_ind_crra<=0.2;


corr_matrix = corr([ti_exp', tp_exp', ti_ord, tp_ord]);
    
% Define the column and row names
column_names = {'EXP_ti', 'EXP_tp', 'EXP_ti ordinal', 'EXP_tp ordinal'};
row_names = {'EXP_ti', 'EXP_tp', 'EXP_ti ordinal', 'EXP_tp ordinal'};

% Create a table for the correlation coefficients
corr_table = array2table(corr_matrix, 'VariableNames', column_names, 'RowNames', row_names)


corr_matrix = corr([ra_crra', rl_crra', ra_ord, rl_ord]);

% Define the column and row names
column_names = {'CRRA_ra', 'CRRA_rl', 'CRRA_ra ordinal', 'CRRA_rl ordinal'};
row_names = {'CRRA_ra', 'CRRA_rl', 'CRRA_ra ordinal', 'CRRA_rl ordinal'};

% Create a table for the correlation coefficients
corr_table = array2table(corr_matrix, 'VariableNames', column_names, 'RowNames', row_names)





%% FIGURES 4-5-6-7

ti = ti_ord;
tp = tp_ord;

ra = ra_ord;
rl = rl_ord;


Times = xlsread('DATI.xlsx','Variables','BN1:BO146');
Response_D = Times(:,1);



% Define bin edges and labels
bin_labels = {'Impatient', 'Others', 'Patient'};


bin(:,2) = ( ~ti&~tp );
bin(:,1) = ( ti ) ;
bin(:,3) = (  tp ) ;

% Calculate mean "Raven" score for each bin
for i = 1:3
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);    
    mean_times(i) = mean(Response_D(bin(:,i)));
    se_times(i) = std(Response_D(bin(:,i))) / sqrt(length(Response_D(bin(:,i))));
end

% Plot the bar graph

bar(mean_times);
hold on
errorbar([1 2 3], mean_times, se_times, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('Response Times', 'FontName', 'Times New Roman');
ylim([300 540])

% Perform t-tests and store results for p-value annotations
[~, p12, ~] = ttest2(Response_D(bin(:,1)), Response_D(bin(:,2)));
[~, p23, ~] = ttest2(Response_D(bin(:,2)), Response_D(bin(:,3)));
[~, p13, ~] = ttest2(Response_D(bin(:,1)), Response_D(bin(:,3)));
% Define x positions for the bars
x_positions = 1:length(mean_times); 

% Y positions for each comparison, slightly above the bar heights
y_bracket_12 = 475; % Bin 1 vs Bin 2
y_bracket_23 = 495; % Bin 2 vs Bin 3
y_bracket_13 = 515; % Bin 1 vs Bin 3
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





Response_L = Times(:,2);



% Define bin edges and labels
bin_labels = {'Risk neutral', 'Others', 'Risk averse'};

% Initialize array to store mean Raven scores for each bin
mean_raven = zeros(1, length(bin_labels));

bin(:,2) = ( ~ra&~rl );
bin(:,1) = ( rl ) ;
bin(:,3) = (  ra ) ;

% Calculate mean "Raven" score for each bin
for i = 1:3
        bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);    
    mean_times(i) = mean(Response_L(bin(:,i)));
    se_times(i) = std(Response_L(bin(:,i))) / sqrt(length(Response_L(bin(:,i))));
end

% Plot the bar graph

bar(mean_times);
hold on
errorbar([1 2 3], mean_times, se_times, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('Response times', 'FontName', 'Times New Roman');
ylim([400 1200])


% Perform t-tests and store results for p-value annotations
[h12, p12, ~] = ttest2(Response_L(bin(:,1)), Response_L(bin(:,2)));
[h23, p23, ~] = ttest2(Response_L(bin(:,2)), Response_L(bin(:,3)));
[h13, p13, ~] = ttest2(Response_L(bin(:,1)), Response_L(bin(:,3)));

% Define x positions for the bars
x_positions = 1:length(mean_times); 

% Y positions for each comparison, slightly above the bar heights
y_bracket_12 = 1000; % Bin 1 vs Bin 2
y_bracket_23 = 1050; % Bin 2 vs Bin 3
y_bracket_13 = 1100; % Bin 1 vs Bin 3
hold on
% First comparison (Bin 1 vs Bin 2)
plot([x_positions(1), x_positions(1), x_positions(2), x_positions(2)], [y_bracket_12, y_bracket_12+10, y_bracket_12+10, y_bracket_12], 'k', 'LineWidth', 1.5);
text(mean(x_positions(1:2)), y_bracket_12 + 30, sprintf('p = %.3f', p12), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Second comparison (Bin 2 vs Bin 3)
plot([x_positions(2), x_positions(2), x_positions(3), x_positions(3)], [y_bracket_23, y_bracket_23+10, y_bracket_23+10, y_bracket_23], 'k', 'LineWidth', 1.5);
text(mean(x_positions(2:3)), y_bracket_23 + 30, sprintf('p = %.3f', p23), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Third comparison (Bin 1 vs Bin 3)
plot([x_positions(1), x_positions(1), x_positions(3), x_positions(3)], [y_bracket_13, y_bracket_13+10, y_bracket_13+10, y_bracket_13], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(3)]), y_bracket_13 + 30, sprintf('p = %.3f', p13), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold off





load warpD.mat

WARP_D = WARPxyzw;

Impatience = xlsread('DATI.xlsx','Variables','AN1:AN146');
Monotonicity = xlsread('DATI.xlsx','Variables','AO1:AO146');



% Define bin edges and labels
bin_labels = {'Impatient', 'Others', 'Others & consistent', 'Patient'};


% Calculate mean "Response_D" score for each bin

clearvars bin

bin(:,2) = ( ~ti&~tp );
bin(:,1) = ( ti ) ;
bin(:,4) = (  tp ) ;

bin = logical(bin);

bin(:,3) = ( ~ti&~tp & WARP_D <= 3 );


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



yyaxis left; % Use left y-axis for mean impatience
bar(mean_impatience)
hold on
errorbar([1 2 3 4], mean_impatience, se_impatience, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman');
ylabel('Share of DMs violating Impatience', 'FontName', 'Times New Roman');
ylim([0 0.8]); yticks(0:0.2:1);

% Plot mean WARP_D on the right axis
yyaxis right;
plot(1:length(mean_warp_d), mean_warp_d, '-o', 'LineWidth', 1, ...
    'MarkerSize', 6, 'Color', [0.8 0.2 0.2]);
ylabel('Mean WARP violations', 'FontName', 'Times New Roman');
ylim([min(mean_warp_d) - 0.1, max(mean_warp_d) + 4]); % Adjust y-limits dynamically

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




SOSD = xlsread('DATI.xlsx','Variables','AS1:AS146');


load warpL.mat

WARP_L = WARPabcd;




bin_labels = {'Risk neutral', 'Others', 'Others & consistent' , 'Risk averse'};


% Calculate mean "Response_D" score for each bin

bin(:,2) = ( ~rl&~ra );
bin(:,1) = ( rl ) ;
bin(:,4) = ( ra ) ;

bin = logical(bin);

bin(:,3) = (~rl&~ra & WARP_L <= 5 );


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





%% main figure

deliberate_time =  xlsread('DATI.xlsx','Variables','FN1:FN146');
deliberate = xlsread('DATI.xlsx','Variables','EE1:EE146');

WARP_D = xlsread('DATI.xlsx','Variables','C1:C146');
WARP_L = xlsread('DATI.xlsx','Variables','E1:E146');


ti = ti_ord;
tp = tp_ord;

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

ra = ra_ord;
rl = rl_ord;

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







Raven = xlsread('DATI.xlsx','Variables','S1:S146');



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




