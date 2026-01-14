

load('C:\Users\caliari\Dropbox\JOB MARKET PAPER\Codes\Creating all datasets\Time Datasets\Datasets_Time.mat')

stationarityData = choiceDelayed(:,[11 16]);

choicesStationarity = zeros(size(stationarityData,1),size(stationarityData,2));

choicesStationarity(stationarityData=="x")=1;
choicesStationarity(stationarityData=="y")=2;
choicesStationarity(stationarityData=="z")=3;
choicesStationarity(stationarityData=="w")=4;

choicesStationarity(stationarityData=="ix")=5;
choicesStationarity(stationarityData=="iy")=6;
choicesStationarity(stationarityData=="iz")=7;
choicesStationarity(stationarityData=="iw")=8;

clearvars -except choicesStationarity

N = size(choicesStationarity,1);

for i=1:N
    testStationarity(i,1) = choicesStationarity(i,2) - choicesStationarity(i,1) == 4; 
end

%% TEST FOR INDEPENDENCE


load('C:\Users\caliari\Dropbox\JOB MARKET PAPER\Codes\Creating all datasets\Risk Datasets\Datasets_Risk.mat')

IndependenceData = choiceLotteries(:,[10 18]);

choicesIndependence = zeros(size(IndependenceData,1),size(IndependenceData,2));

choicesIndependence(IndependenceData=="b")=1;
choicesIndependence(IndependenceData=="c")=2;
choicesIndependence(IndependenceData=="d")=3;

choicesIndependence(IndependenceData=="f2b")=5;
choicesIndependence(IndependenceData=="f2c")=6;
choicesIndependence(IndependenceData=="f2d")=7;


clearvars -except choicesIndependence testStationarity

N = size(choicesIndependence,1);

for i=1:N
    testIndependence(i,1) = choicesIndependence(i,2) - choicesIndependence(i,1) == 4; 
end


%%

Raven = xlsread('DATI.xlsx','Variables','S1:S146');
Impatience = xlsread('DATI.xlsx','Variables','AN1:AN146');
FOSD = xlsread('DATI.xlsx','Variables','AP1:AQ146');
FOSD = sum(FOSD,2)>0;
WARP_D = xlsread('DATI.xlsx','Variables','C1:C146');
WARP_L = xlsread('DATI.xlsx','Variables','E1:E146');
deliberate_time =  xlsread('DATI.xlsx','Variables','FN1:FN146');
deliberate = xlsread('DATI.xlsx','Variables','EE1:EE146');






% Separate Raven values based on Impatience categories
Raven_0 = Raven(testStationarity(Impatience==0));
Raven_1 = Raven(~testStationarity(Impatience==0));

% Calculate means for each group
mean_Raven_0 = mean(Raven_0);
mean_Raven_1 = mean(Raven_1);

% Calculate standard errors for each group
se_Raven_0 = std(Raven_0) / sqrt(length(Raven_0));
se_Raven_1 = std(Raven_1) / sqrt(length(Raven_1));

% Perform the t-test
[~, p_value] = ttest2(Raven_0, Raven_1);


% Create the bar plot
subplot(1,2,1)

bar([mean_Raven_0, mean_Raven_1]); 
hold on;

% Add error bars
errorbar([1, 2], [mean_Raven_0, mean_Raven_1], [se_Raven_0, se_Raven_1], '.k', 'LineWidth', 1);

% Set the x-axis labels and title
set(gca, 'XTickLabel', {'Satisfied', 'Violated'});
ylabel('Raven Score');
title('Stationarity')
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
Raven_0 = Raven(testIndependence(FOSD==0));
Raven_1 = Raven(~testIndependence(FOSD==0));

% Calculate means for each group
mean_Raven_0 = mean(Raven_0);
mean_Raven_1 = mean(Raven_1);

% Calculate standard errors for each group
se_Raven_0 = std(Raven_0) / sqrt(length(Raven_0));
se_Raven_1 = std(Raven_1) / sqrt(length(Raven_1));

% Perform the t-test
[~, p_value] = ttest2(Raven_0, Raven_1);


% Create the bar plot
subplot(1,2,2)

bar([mean_Raven_0, mean_Raven_1]); 
hold on;

% Add error bars
errorbar([1, 2], [mean_Raven_0, mean_Raven_1], [se_Raven_0, se_Raven_1], '.k', 'LineWidth', 1);

% Set the x-axis labels and title
set(gca, 'XTickLabel', {'Satisfied', 'Violated'});
ylabel('Raven Score');
title('Independence')
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




x = [sum(testIndependence(deliberate==1)),sum(~testIndependence(deliberate==1));
    sum(testIndependence(deliberate==0)),sum(~testIndependence(deliberate==0))];

[h,p] = fishertest(x)

