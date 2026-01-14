clear
clc

structural_estimation_3rdSubmission_CRRACARAEXPBETA  % estimates under the different specifications

%% Figures 2 and 3 - Joint distribution preference and precision parameters %%


time_ind =  1./(1+time_ind_exp);
prec_ind_time = prec_ind_exp;

% Joint distribution scatter plot
bubble_time = [time_ind', prec_ind_time'];
bubble_time = sortrows(bubble_time, 1);
bubble_time(:, 1) = round(bubble_time(:, 1), 6);

% Count occurrences for point sizes
for i = 1:size(time_ind, 2)
    bubble_time(i, 3) = sum(bubble_time(:, 1) == bubble_time(i, 1));
end
bubble_time(:, 3) = bubble_time(:, 3) * 20; % Scale factor for scatter sizes

% Create a figure with two subplots: one for the scatter plot and one for the histogram
figure;

% Top subplot: Scatter plot
subplot(2, 1, 1); % Creates a 2-row, 1-column layout and selects the first subplot
scatter(bubble_time(:, 1)', bubble_time(:, 2)', bubble_time(:, 3)');
xlabel('Discount factor');
ylabel('Precision');

% Bottom subplot: Histogram for the marginal distribution of time_ind
subplot(2, 1, 2); % Selects the second subplot in the layout
histogram(time_ind, 'BinWidth', 0.01, 'Normalization', 'probability', 'FaceAlpha', 0.6);
xlabel('Discount factor');
ylabel('Probability');
xlim([min(time_ind)-0.001, max(time_ind)+0.001]); % Align x-axis with scatter plot for consistency



risk_ind = risk_ind_crra;
prec_ind = prec_ind_crra;

% Joint distribution scatter plot
bubble_risk = [risk_ind', log(prec_ind)'];
bubble_risk = sortrows(bubble_risk, 1);
bubble_risk(:, 1) = round(bubble_risk(:, 1), 6);

% Count occurrences for point sizes
for i = 1:size(risk_ind, 2)
    bubble_risk(i, 3) = sum(bubble_risk(:, 1) == bubble_risk(i, 1));
end
bubble_risk(:, 3) = bubble_risk(:, 3) * 20; % Scale factor for scatter sizes

% Create a figure with two subplots: one for the scatter plot and one for the histogram
figure;

% Top subplot: Scatter plot
subplot(2, 1, 1); % Creates a 2-row, 1-column layout and selects the first subplot
scatter(bubble_risk(:, 1)', bubble_risk(:, 2)', bubble_risk(:, 3)');
xlabel('Risk parameter');
ylabel('Precision');
xlim([min(risk_ind)-0.001, max(risk_ind)+0.001]); 

% Bottom subplot: Histogram for the marginal distribution of time_ind
subplot(2, 1, 2); % Selects the second subplot in the layout
histogram(risk_ind, 'BinWidth', 0.2, 'Normalization', 'probability', 'FaceAlpha', 0.6);
xlabel('Risk parameter');
ylabel('Probability');
xlim([min(risk_ind)-0.001, max(risk_ind)+0.001]); % Align x-axis with scatter plot for consistency


