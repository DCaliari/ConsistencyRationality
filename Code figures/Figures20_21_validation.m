clear
clc

load 'Time Ranking.mat' % load the full reported ranking
load 'Results Time MAIN.mat' % load the full sequential, swaps relations

load variables.mat

Na = 4; % number of alternatives
[any,irreflexive,asymmetric,acyclic,transitive,linear,intervalorder,semiorder] = binaryrelations(Na);

T = linear;

% these are indicators of the which element is the best in the linear
% orders indexed in the matrix T
indA = [1 2 3 4 5 6];
indB = [13 14 15 17 19 20];
indC = [9 10 11 18 22 23];
indD = [7 8 12 16 21 24];

type_D = zeros(1,145);



for i=1:145
    for j=1:24
    if isequal(seqsolDelayed(:,:,i),T(:,:,j))
        type_D(1,i) = j;
    else
    end
    end
end



beta = 1./(1+time_ind_exp');

beta_type = [type_D',beta];



% Calculate the means for each group
meanA = mean(beta_type(ismember(beta_type(:,1), indA), 2));
meanB = mean(beta_type(ismember(beta_type(:,1), indB), 2));
meanC = mean(beta_type(ismember(beta_type(:,1), indC), 2));
meanD = mean(beta_type(ismember(beta_type(:,1), indD), 2));

% Store the means in a vector
means = [meanA, meanB, meanC, meanD];


std_devs = [
    std(beta_type(ismember(beta_type(:,1), indA), 2)), ...
    std(beta_type(ismember(beta_type(:,1), indB), 2)), ...
    std(beta_type(ismember(beta_type(:,1), indC), 2)), ...
    std(beta_type(ismember(beta_type(:,1), indD), 2))
];

n = [
    sum(ismember(beta_type(:, 1), indA)), ...
    sum(ismember(beta_type(:, 1), indB)), ...
    sum(ismember(beta_type(:, 1), indC)), ...
    sum(ismember(beta_type(:, 1), indD))
];

% Calculate standard errors
std_errors = std_devs ./ sqrt(n);


means_sequential = means;
std_errors_sequential = std_errors;


for i=1:145
    for j=1:24
    if isequal(Time_Ranking(:,:,i),T(:,:,j))
        type_D(1,i) = j;
    else
    end
    end
end



beta_type = [type_D',beta];



% Calculate the means for each group
meanA = mean(beta_type(ismember(beta_type(:,1), indA), 2));
meanB = mean(beta_type(ismember(beta_type(:,1), indB), 2));
meanC = mean(beta_type(ismember(beta_type(:,1), indC), 2));
meanD = mean(beta_type(ismember(beta_type(:,1), indD), 2));

% Store the means in a vector
means = [meanA, meanB, meanC, meanD];


std_devs = [
    std(beta_type(ismember(beta_type(:,1), indA), 2)), ...
    std(beta_type(ismember(beta_type(:,1), indB), 2)), ...
    std(beta_type(ismember(beta_type(:,1), indC), 2)), ...
    std(beta_type(ismember(beta_type(:,1), indD), 2))
];

n = [
    sum(ismember(beta_type(:, 1), indA)), ...
    sum(ismember(beta_type(:, 1), indB)), ...
    sum(ismember(beta_type(:, 1), indC)), ...
    sum(ismember(beta_type(:, 1), indD))
];

% Calculate standard errors
std_errors = std_devs ./ sqrt(n);


means_reported = means;
std_errors_reported = std_errors;



type_D = zeros(1,145);
for i=1:145
    temp = zeros(1,4);
    for j=1:4
    temp(j) = sum(minswapssolDelayed(j,:,i));
    end
    for j=1:4
    if temp(j)==3 && sum(temp(j)>temp)==3
        type_D(1,i) = j;
    else
    end
    end
end



beta_type = [type_D',beta];



% Calculate the means for each group
meanA = mean(beta_type(ismember(beta_type(:,1), [1]), 2));
meanB = mean(beta_type(ismember(beta_type(:,1), [2]), 2));
meanC = mean(beta_type(ismember(beta_type(:,1), [3]), 2));
meanD = mean(beta_type(ismember(beta_type(:,1), [4]), 2));

% Store the means in a vector
means = [meanA, meanB, meanC, meanD];


std_devs = [
    std(beta_type(ismember(beta_type(:,1), [1]), 2)), ...
    std(beta_type(ismember(beta_type(:,1), [2]), 2)), ...
    std(beta_type(ismember(beta_type(:,1), [3]), 2)), ...
    std(beta_type(ismember(beta_type(:,1), [4]), 2))
];

n = [
    sum(ismember(beta_type(:, 1), [1])), ...
    sum(ismember(beta_type(:, 1), [2])), ...
    sum(ismember(beta_type(:, 1), [3])), ...
    sum(ismember(beta_type(:, 1), [4]))
];

z = sum(beta_type(:,1)==0);

% Calculate standard errors
std_errors = std_devs ./ sqrt(n);


means_swaps = means;
std_errors_swaps = std_errors;



% Combine means and errors
means = [means_swaps; means_sequential; means_reported]; % 2 rows for 2 groups
std_errors = [std_errors_swaps; std_errors_sequential; std_errors_reported]; % 2 rows for 2 groups




% Bar width adjustments
bar_width = 1; % Bar width
x_positions = 1:size(means, 2); % X positions for the bars

% Create grouped bar plot with two sets
hBar = bar(x_positions, means', 'grouped'); % Use 'grouped' for clarity
hold on; % Keep the current figure

% Now, add error bars
% Calculate the positions directly from the bar objects
for i = 1:size(means, 1) % Loop over the two groups
    % Get the X data for the current group from the bar object
    % Since the bars are grouped, we can access their data easily
    xCenter = hBar(i).XData + (i - 2) * bar_width / 4.5; % Correctly calculate x-center positions
    errorbar(xCenter, means(i, :), std_errors(i, :), 'k', 'linestyle', 'none', 'LineWidth', 1);
end

% Set x-tick labels
set(gca, 'xtick', x_positions, 'xticklabel', {'OS', 'D', 'K', 'I'}, 'FontName', 'Times New Roman', 'FontSize', 12);

% Set y-axis limits
ylim([0.9, 1]); % Adjust as necessary

% Add title and labels
xlabel('Elicited best alternative', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Estimated Discount Factor', 'FontName', 'Times New Roman', 'FontSize', 12);

% Add legend
legend({'Minimum Swaps','Sequential', 'Reported'}, 'Location', 'northwest', 'FontName', 'Times New Roman', 'FontSize', 12);

% Set font for x and y axes and other text elements
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
hold off; % Release the current figure











load 'Risk Ranking.mat' % load the full reported ranking
load 'Results Risk MAIN.mat' % load the full sequential, swaps relations


% these are indicators of the which element is the best in the linear
% orders indexed in the matrix T
indA = [1 2 3 4 5 6];
indB = [13 14 15 17 19 20];
indC = [9 10 11 18 22 23];
indD = [7 8 12 16 21 24];


type_L = zeros(1,145);


for i=1:145
    for j=1:24
    if isequal(seqsolLotteries(:,:,i),T(:,:,j))
        type_L(1,i) = j;
    else
    end
    end
end

rho = risk_ind_crra';

rho_type = [type_L',rho];



% Calculate the means for each group
meanA = mean(rho_type(ismember(rho_type(:,1), indA), 2));
meanB = mean(rho_type(ismember(rho_type(:,1), indB), 2));
meanC = mean(rho_type(ismember(rho_type(:,1), indC), 2));
meanD = mean(rho_type(ismember(rho_type(:,1), indD), 2));

% Store the means in a vector
means = [meanA, meanB, meanC, meanD];


std_devs = [
    std(rho_type(ismember(rho_type(:,1), indA), 2)), ...
    std(rho_type(ismember(rho_type(:,1), indB), 2)), ...
    std(rho_type(ismember(rho_type(:,1), indC), 2)), ...
    std(rho_type(ismember(rho_type(:,1), indD), 2))
];

n = [
    sum(ismember(rho_type(:, 1), indA)), ...
    sum(ismember(rho_type(:, 1), indB)), ...
    sum(ismember(rho_type(:, 1), indC)), ...
    sum(ismember(rho_type(:, 1), indD))
];

% Calculate standard errors
std_errors = std_devs ./ sqrt(n);


means_sequential = means;
std_errors_sequential = std_errors;


for i=1:145
    for j=1:24
    if isequal(Risk_Ranking(:,:,i),T(:,:,j))
        type_L(1,i) = j;
    else
    end
    end
end


rho_type = [type_L',rho];



% Calculate the means for each group
meanA = mean(rho_type(ismember(rho_type(:,1), indA), 2));
meanB = mean(rho_type(ismember(rho_type(:,1), indB), 2));
meanC = mean(rho_type(ismember(rho_type(:,1), indC), 2));
meanD = mean(rho_type(ismember(rho_type(:,1), indD), 2));

% Store the means in a vector
means = [meanA, meanB, meanC, meanD];


std_devs = [
    std(rho_type(ismember(rho_type(:,1), indA), 2)), ...
    std(rho_type(ismember(rho_type(:,1), indB), 2)), ...
    std(rho_type(ismember(rho_type(:,1), indC), 2)), ...
    std(rho_type(ismember(rho_type(:,1), indD), 2))
];

n = [
    sum(ismember(rho_type(:, 1), indA)), ...
    sum(ismember(rho_type(:, 1), indB)), ...
    sum(ismember(rho_type(:, 1), indC)), ...
    sum(ismember(rho_type(:, 1), indD))
];

% Calculate standard errors
std_errors = std_devs ./ sqrt(n);


means_reported = means;
std_errors_reported = std_errors;



type_L = zeros(1,145);
for i=1:145
    temp = zeros(1,4);
    for j=1:4
    temp(j) = sum(minswapssolLotteries(j,:,i));
    end
    for j=1:4
    if temp(j)==3 && sum(temp(j)>temp)==3
        type_L(1,i) = j;
    else
    end
    end
end


rho_type = [type_L',rho];



% Calculate the means for each group
meanA = mean(rho_type(ismember(rho_type(:,1), [1]), 2));
meanB = mean(rho_type(ismember(rho_type(:,1), [2]), 2));
meanC = mean(rho_type(ismember(rho_type(:,1), [3]), 2));
meanD = mean(rho_type(ismember(rho_type(:,1), [4]), 2));

% Store the means in a vector
means = [meanA, meanB, meanC, meanD];


std_devs = [
    std(rho_type(ismember(rho_type(:,1), [1]), 2)), ...
    std(rho_type(ismember(rho_type(:,1), [2]), 2)), ...
    std(rho_type(ismember(rho_type(:,1), [3]), 2)), ...
    std(rho_type(ismember(rho_type(:,1), [4]), 2))
];

n = [
    sum(ismember(rho_type(:, 1), [1])), ...
    sum(ismember(rho_type(:, 1), [2])), ...
    sum(ismember(rho_type(:, 1), [3])), ...
    sum(ismember(rho_type(:, 1), [4]))
];

z = sum(rho_type(:,1)==0);

% Calculate standard errors
std_errors = std_devs ./ sqrt(n);


means_swaps = means;
std_errors_swaps = std_errors;



% Combine means and errors
means = [means_swaps; means_sequential; means_reported]; % 2 rows for 2 groups
std_errors = [std_errors_swaps; std_errors_sequential; std_errors_reported]; % 2 rows for 2 groups



% Assuming means and std_errors have already been calculated as before

% Create bar plot with error bars
figure; % Create a new figure

% Bar width adjustments
bar_width = 1; % Bar width
x_positions = 1:size(means, 2); % X positions for the bars

% Create grouped bar plot with two sets
hBar = bar(x_positions, means', 'grouped'); % Use 'grouped' for clarity
hold on; % Keep the current figure

% Now, add error bars
% Calculate the positions directly from the bar objects
for i = 1:size(means, 1) % Loop over the two groups
    % Get the X data for the current group from the bar object
    % Since the bars are grouped, we can access their data easily
    xCenter = hBar(i).XData + (i - 2) * bar_width / 4.5; % Correctly calculate x-center positions
    errorbar(xCenter, means(i, :), std_errors(i, :), 'k', 'linestyle', 'none', 'LineWidth', 1);
end

% Set x-tick labels
set(gca, 'xtick', x_positions, 'xticklabel', {'D', 'S', '50', 'R'}, 'FontName', 'Times New Roman', 'FontSize', 12);

% Set y-axis limits
ylim([0, 2.4]); % Adjust as necessary

% Add title and labels
xlabel('Elicited best alternative', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Estimated Risk Parameter', 'FontName', 'Times New Roman', 'FontSize', 12);

% Add legend
legend({'Minimum Swaps','Sequential', 'Reported'}, 'Location', 'northeast', 'FontName', 'Times New Roman', 'FontSize', 12);

% Set font for x and y axes and other text elements
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
hold off; % Release the current figure







