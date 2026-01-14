clear
clc

%% Figure 11 - Violations of PrAx and Raven scores %%

Impatience = xlsread('DATI.xlsx','Variables','AN1:AN146');
Monotonicity = xlsread('DATI.xlsx','Variables','AO1:AO146');

FOSD = xlsread('DATI.xlsx','Variables','AP1:AQ146');
FOSD = sum(FOSD,2)>0;
SOSD = xlsread('DATI.xlsx','Variables','AS1:AS146');

Raven = xlsread('DATI.xlsx','Variables','S1:S146');


% Separate Raven values based on Impatience categories
Raven_0 = Raven(Monotonicity == 0);
Raven_1 = Raven(Monotonicity == 1);

% Calculate means for each group
mean_Raven_0 = mean(Raven_0);
mean_Raven_1 = mean(Raven_1);

% Calculate standard errors for each group
se_Raven_0 = std(Raven_0) / sqrt(length(Raven_0));
se_Raven_1 = std(Raven_1) / sqrt(length(Raven_1));

% Perform the t-test
[~, p_value] = ttest2(Raven_0, Raven_1);


% Create the bar plot
subplot(1,4,1)

bar([mean_Raven_0, mean_Raven_1]); 
hold on;

% Add error bars
errorbar([1, 2], [mean_Raven_0, mean_Raven_1], [se_Raven_0, se_Raven_1], '.k', 'LineWidth', 1);

% Set the x-axis labels and title
set(gca, 'XTickLabel', {'Satisfied', 'Violated'});
ylabel('Raven Score');
title('Monotonicity')
ylim([0 14])


% Add p-value to the plot
if p_value < 0.05
    sig_text = sprintf('p = %.3f*', p_value); % Add * for significance
else
    sig_text = sprintf('p = %.3f', p_value);
end
text(1.5, 13, sig_text, ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');


% Separate Raven values based on Impatience categories
Raven_0 = Raven(Impatience == 0);
Raven_1 = Raven(Impatience == 1);

% Calculate means for each group
mean_Raven_0 = mean(Raven_0);
mean_Raven_1 = mean(Raven_1);

% Calculate standard errors for each group
se_Raven_0 = std(Raven_0) / sqrt(length(Raven_0));
se_Raven_1 = std(Raven_1) / sqrt(length(Raven_1));

% Perform the t-test
[~, p_value] = ttest2(Raven_0, Raven_1);


% Create the bar plot
subplot(1,4,2)

bar([mean_Raven_0, mean_Raven_1]); 
hold on;

% Add error bars
errorbar([1, 2], [mean_Raven_0, mean_Raven_1], [se_Raven_0, se_Raven_1], '.k', 'LineWidth', 1);

% Set the x-axis labels and title
set(gca, 'XTickLabel', {'Satisfied', 'Violated'});
ylabel('Raven Score');
title('Impatience')
ylim([0 14])

% Add p-value to the plot
if p_value < 0.05
    sig_text = sprintf('p = %.3f*', p_value); % Add * for significance
else
    sig_text = sprintf('p = %.3f', p_value);
end
text(1.5, 13, sig_text, ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');


% Separate Raven values based on Impatience categories
Raven_0 = Raven(FOSD == 0);
Raven_1 = Raven(FOSD == 1);

% Calculate means for each group
mean_Raven_0 = mean(Raven_0);
mean_Raven_1 = mean(Raven_1);

% Calculate standard errors for each group
se_Raven_0 = std(Raven_0) / sqrt(length(Raven_0));
se_Raven_1 = std(Raven_1) / sqrt(length(Raven_1));

% Perform the t-test
[~, p_value] = ttest2(Raven_0, Raven_1);


% Create the bar plot
subplot(1,4,3)

bar([mean_Raven_0, mean_Raven_1]); 
hold on;

% Add error bars
errorbar([1, 2], [mean_Raven_0, mean_Raven_1], [se_Raven_0, se_Raven_1], '.k', 'LineWidth', 1);

% Set the x-axis labels and title
set(gca, 'XTickLabel', {'Satisfied', 'Violated'});
ylabel('Raven Score');
title('FOSD')
ylim([0 14])

% Add p-value to the plot
if p_value < 0.05
    sig_text = sprintf('p = %.3f*', p_value); % Add * for significance
else
    sig_text = sprintf('p = %.3f', p_value);
end
text(1.5, 13, sig_text, ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');


% Separate Raven values based on Impatience categories
Raven_0 = Raven(SOSD == 0);
Raven_1 = Raven(SOSD == 1);

% Calculate means for each group
mean_Raven_0 = mean(Raven_0);
mean_Raven_1 = mean(Raven_1);

% Calculate standard errors for each group
se_Raven_0 = std(Raven_0) / sqrt(length(Raven_0));
se_Raven_1 = std(Raven_1) / sqrt(length(Raven_1));

% Perform the t-test
[~, p_value] = ttest2(Raven_0, Raven_1);


% Create the bar plot
subplot(1,4,4)

bar([mean_Raven_0, mean_Raven_1]); 
hold on;

% Add error bars
errorbar([1, 2], [mean_Raven_0, mean_Raven_1], [se_Raven_0, se_Raven_1], '.k', 'LineWidth', 1);

% Set the x-axis labels and title
set(gca, 'XTickLabel', {'Satisfied', 'Violated'});
ylabel('Raven Score');
title('SOSD')
ylim([0 14])

% Add p-value to the plot
if p_value < 0.05
    sig_text = sprintf('p = %.3f*', p_value); % Add * for significance
else
    sig_text = sprintf('p = %.3f', p_value);
end
text(1.5, 13, sig_text, ...
     'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');

set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');





