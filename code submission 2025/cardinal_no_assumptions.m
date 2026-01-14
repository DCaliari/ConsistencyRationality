%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TIME Preferences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

loadtimedata
Sets=Sets==0;

%% Logit estimation of utilities v in Time

options=optimoptions('fmincon','display','off','MaxFunctionEvaluations',1e7,...
    'StepTolerance',1e-9,'FunctionTolerance',1e-9,'OptimalityTolerance',1e-9,'FiniteDifferenceType','central' );

utilities_time = zeros(Na,N_subjects);

for ind = 1:N_subjects
x0 = [0.25 0.25 0.25 0.25]';

Aeq = [1 1 1 1];
beq = 1;
lb = [0; 0; 0; 0];
ub = [1; 1; 1; 1];

% [utilities_time(:,ind),Ltime(ind,1)] = fmincon(@(x) loglikelihood_logit(x, Data, Na, Main_sets, Sets, ind), x0, [], [], Aeq, beq, lb, ub, [], options);
[utilities_time(:,ind),Ltime(ind,1)] = fmincon(@(x) loglikelihood_logit(x, Data, Na, Main_sets, Sets, ind), x0, [], [], Aeq, beq, [], [], [], options);

end

DataTime=Data;
clearvars -except utilities_time Na Main_sets Sets N_subjects Ltime DataTime


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RISK Preferences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



loadriskdata
Sets=Sets==0;

%% Description of the lotteries


options=optimoptions('fmincon','display','off','MaxFunctionEvaluations',1e7,...
    'StepTolerance',1e-9,'FunctionTolerance',1e-9,'OptimalityTolerance',1e-9,'FiniteDifferenceType','central' );

utilities_risk = zeros(Na, N_subjects);

for ind = 1:N_subjects
x0 = [0.25 0.25 0.25 0.25]';

Aeq = [1 1 1 1];
beq = 1;
lb = [0; 0; 0; 0];
ub = [1; 1; 1; 1];

% [utilities_risk(:,ind), Lrisk(ind,1)] = fmincon(@(x) loglikelihood_logit(x, Data, Na, Main_sets, Sets, ind), x0, [], [], Aeq, beq, lb, ub, [], options);
[utilities_risk(:,ind), Lrisk(ind,1)] = fmincon(@(x) loglikelihood_logit(x, Data, Na, Main_sets, Sets, ind), x0, [], [], Aeq, beq, [], [], [], options);
end


DataRisk=Data;

clearvars -except utilities_time utilities_risk Na Main_sets Sets N_subjects Ltime Lrisk DataTime DataRisk options


RatIndexTime = zeros(N_subjects,3);
RatIndexRisk = zeros(N_subjects,3);

for i=1:N_subjects
    temp = sort(utilities_risk(:,i),'descend');
    RatIndexRisk(i,1) = exp(temp(1))/sum(exp(temp));
    RatIndexRisk(i,2) = exp(temp(2))/sum(exp(temp(2))+exp(temp(3))+exp(temp(4)));
    RatIndexRisk(i,3) = exp(temp(3))/sum(exp(temp(3))+exp(temp(4)));
    temp = sort(utilities_time(:,i),'descend');
    RatIndexTime(i,1) = exp(temp(1))/sum(exp(temp));
    RatIndexTime(i,2) = exp(temp(2))/sum(exp(temp(2))+exp(temp(3))+exp(temp(4)));
    RatIndexTime(i,3) = exp(temp(3))/sum(exp(temp(3))+exp(temp(4)));
end

% utility differences
for i=1:N_subjects
    temp = sort(utilities_risk(:,i),'descend');
    UdiffRisk(i,1) = 3*(temp(1)-temp(2))+4*(temp(2)-temp(3))+3*(temp(3)-temp(4));
    MdiffRisk(i,1) = UdiffRisk(i,1)/6;
    temp = sort(utilities_time(:,i),'descend');
    UdiffTime(i,1) = 3*(temp(1)-temp(2))+4*(temp(2)-temp(3))+3*(temp(3)-temp(4));
    MdiffTime(i,1) = UdiffTime(i,1)/6;
end


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


%% HOW MANY PARTICIPANTS AGREE WITH CRRA AND CARA ORDERINGS

    
clearvars type type2 

indifutilities_risk = round(utilities_risk,2);

for i=1:N_subjects
[~,rank_risk(i,:)]=sort(indifutilities_risk(:,i),'descend');
[~, rank_risk2(i,:)] = sortrows([(1:numel(indifutilities_risk(:,i)))' indifutilities_risk(:,i)], [-2, -1]);
end


ra = zeros(N_subjects,1);
ra2 = zeros(N_subjects,1);
rl = zeros(N_subjects,1);
rl2 = zeros(N_subjects,1);

for i=1:N_subjects
     for j=1:length(PERcrra)
     if isequal(rank_risk(i,:),PERcrra(j,:))
        type(i,j)=1; 
     end
     if isequal(rank_risk(i,:),PERcrra(1,:))
        ra(i,1)=1; 
     end   
     if isequal(rank_risk(i,:),PERcrra(end,:))
        rl(i,1)=1; 
     end        
     if isequal(rank_risk2(i,:),PERcrra(j,:))
        type2(i,j)=1;
     end     
     if isequal(rank_risk2(i,:),PERcrra(1,:))
        ra2(i,1)=1; 
     end
     if isequal(rank_risk2(i,:),PERcrra(end,:))
        rl2(i,1)=1; 
     end
     end
end


temp1 = [sum(type,2), sum(type2,2)];
temp1 = sum(temp1,2)>0;
proportion_crra = sum(temp1); % proportion the is aligned with crra
ra_crra = sum([ra,ra2],2)>0;
rl_crra = sum([rl,rl2],2)>0;
sum(ra_crra)
sum(rl_crra)


clearvars type type2

ra = zeros(N_subjects,1);
ra2 = zeros(N_subjects,1);
rl = zeros(N_subjects,1);
rl2 = zeros(N_subjects,1);

for i=1:N_subjects
     for j=1:length(PERcara)
     if isequal(rank_risk(i,:),PERcara(j,:))
        type(i,j)=1; 
     end
     if isequal(rank_risk(i,:),PERcara(1,:))
        ra(i,1)=1; 
     end   
     if isequal(rank_risk(i,:),PERcara(end,:))
        rl(i,1)=1; 
     end        
     if isequal(rank_risk2(i,:),PERcara(j,:))
        type2(i,j)=1;
     end     
     if isequal(rank_risk2(i,:),PERcara(1,:))
        ra2(i,1)=1; 
     end
     if isequal(rank_risk2(i,:),PERcara(end,:))
        rl2(i,1)=1; 
     end
     end
end


temp1 = [sum(type,2), sum(type2,2)];
temp1 = sum(temp1,2)>0;
proportion_cara = sum(temp1); % proportion the is aligned with crra
ra_cara = sum([ra,ra2],2)>0;
rl_cara = sum([rl,rl2],2)>0;
sum(ra_cara)
sum(rl_cara)





indifutilities_time = round(utilities_time,2);

for i=1:N_subjects
[~,rank_time(i,:)]=sort(indifutilities_time(:,i),'descend');
[~, rank_time2(i,:)] = sortrows([(1:numel(indifutilities_time(:,i)))' indifutilities_time(:,i)], [-2, -1]);
end

clearvars type type2

ti = zeros(N_subjects,1);
ti2 = zeros(N_subjects,1);
tp = zeros(N_subjects,1);
tp2 = zeros(N_subjects,1);


for i=1:N_subjects
     for j=1:length(PERcrra)
     if isequal(rank_time(i,:),PERexp(j,:))
        type(i,j)=1; 
     end
     if isequal(rank_risk(i,:),PERexp(1,:))
        ti(i,1)=1; 
     end   
     if isequal(rank_risk(i,:),PERexp(end,:))
        tp(i,1)=1; 
     end         
     if isequal(rank_time2(i,:),PERexp(j,:))
        type2(i,j)=1;
     end
     if isequal(rank_time2(i,:),PERexp(1,:))
        ti2(i,1)=1; 
     end
     if isequal(rank_time2(i,:),PERexp(end,:))
        tp2(i,1)=1; 
     end     
     end
end


temp3 = [sum(type,2), sum(type2,2)];
temp3 = sum(temp3,2)>0;
sum(temp3) % proportion the is aligned with crra
ti_exp = sum([ti,ti2],2)>0;
tp_exp = sum([tp,tp2],2)>0;
sum(ti_exp)
sum(tp_exp)


%% correlation with WARP

load warpD.mat
load warpL.mat

corr(MdiffTime, WARPxyzw)
corr(MdiffRisk, WARPabcd)

%% FIGURE DIFFERENCE UTILITIES

deliberate_time =  xlsread('DATI.xlsx','Variables','FN1:FN146');
deliberate = xlsread('DATI.xlsx','Variables','EE1:EE146');



ti = ti_exp;
tp = tp_exp;


subplot(1,2,1)

% Define bin edges and labels
bin_labels = {'Impatient', 'Patient', 'Randomize', 'Others'};


% Calculate mean "Response_D" score for each bin

bin(:,3) = deliberate_time;
bin(:,1) = ( ti & deliberate_time==0 ) ;
bin(:,2) = (  tp & deliberate_time==0  ) ;
bin(:,4) = bin(:,1)+bin(:,2)+bin(:,3)==0;



% Calculate mean "Raven" score for each bin
for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    mean_diff(i) = mean(WARPxyzw(bin(:,i)));
    se_diff(i) = std(WARPxyzw(bin(:,i))) / sqrt(length(WARPxyzw(bin(:,i))));
end


% Plot the bar graph
bar(mean_diff);
hold on
errorbar([1 2 3 4], mean_diff, se_diff, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('WARP violations', 'FontName', 'Times New Roman');
ylim([0 12])

% Perform t-tests and store results for p-value annotations
[~, p14, ~] = ttest2(WARPxyzw(bin(:,1)), WARPxyzw(bin(:,4)));
[~, p24, ~] = ttest2(WARPxyzw(bin(:,2)), WARPxyzw(bin(:,4)));
[~, p34, ~] = ttest2(WARPxyzw(bin(:,3)), WARPxyzw(bin(:,4)));

% Define x positions for the bars
x_positions = 1:length(mean_diff); 

% Y positions for each comparison, slightly above the bar heights
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


ra = ra_crra;
rl = rl_crra;

bin(:,3) = deliberate;
bin(:,2) = ( rl & deliberate ==0 ) ;
bin(:,1) = (  ra  & deliberate ==0 ) ;
bin(:,4) = bin(:,1)+bin(:,2)+bin(:,3)==0;


% Initialize array to store mean Raven scores for each bin
mean_diff = zeros(1, length(bin_labels));

for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);    
    mean_diff(i) = mean(WARPabcd(bin(:,i)));
    se_diff(i) = std(WARPabcd(bin(:,i))) / sqrt(length(WARPabcd(bin(:,i))));
end


% Plot the bar graph
bar(mean_diff);
hold on
errorbar([1 2 3 4], mean_diff, se_diff, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('WARP violations', 'FontName', 'Times New Roman');
ylim([0 12])

% Perform t-tests and store results for p-value annotations
[~, p14, ~] = ttest2(WARPabcd(bin(:,1)), WARPabcd(bin(:,4)));
[~, p24, ~] = ttest2(WARPabcd(bin(:,2)), WARPabcd(bin(:,4)));
[~, p34, ~] = ttest2(WARPabcd(bin(:,3)), WARPabcd(bin(:,4)));

% Define x positions for the bars
x_positions = 1:length(mean_diff); 

% Y positions for each comparison, slightly above the bar heights
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
hold off



%%




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


y_bracket_34 = 16.5;  % Bin 3 vs Bin 4
y_bracket_14 = 12.5;  % Bin 1 vs Bin 4
y_bracket_24 = 14.5; % Bin 2 vs Bin 4

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





%%

Times = xlsread('DATI.xlsx','Variables','BN1:BO146');
Response_D = Times(:,1);
Response_L = Times(:,2);




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

% Plot the bar graph
bar(mean_times)
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


%% FIGURE 7



% Define bin edges and labels
bin_labels = {'Risk neutral', 'Others', 'Risk averse'};

% Initialize array to store mean Raven scores for each bin
mean_times = zeros(1, length(bin_labels));
se_times = zeros(1, length(bin_labels));

bin(:,2) = ( ~ra&~rl );
bin(:,3) = ( ra ) ;
bin(:,1) = (  rl ) ;


% Calculate mean "Raven" score for each bin
for i = 1:length(bin_labels)
    bin_count = sum(bin(:,i));
    bin_labels{i} = sprintf('%s (n=%d)', bin_labels{i}, bin_count);
    mean_times(i) = mean(Response_L(bin(:,i)));
    se_times(i) = std(Response_L(bin(:,i))) / sqrt(length(Response_L(bin(:,i))));
end

% Plot the bar graph
bar(mean_times)
hold on
errorbar([1 2 3], mean_times, se_times, '.k', 'LineWidth', 1);
set(gca, 'XTickLabel', bin_labels, 'FontName', 'Times New Roman'); % Set x-axis labels to bin ranges
ylabel('Response times', 'FontName', 'Times New Roman');
ylim([400 1500])


% Perform t-tests and store results for p-value annotations
[h12, p12, ~] = ttest2(Response_L(bin(:,1)), Response_L(bin(:,2)));
[h23, p23, ~] = ttest2(Response_L(bin(:,2)), Response_L(bin(:,3)));
[h13, p13, ~] = ttest2(Response_L(bin(:,1)), Response_L(bin(:,3)));

% Define x positions for the bars
x_positions = 1:length(mean_times); 

% Y positions for each comparison, slightly above the bar heights
y_bracket_12 = 1240; % Bin 1 vs Bin 2
y_bracket_23 = 1300; % Bin 2 vs Bin 3
y_bracket_13 = 1360; % Bin 1 vs Bin 3
hold on
% First comparison (Bin 1 vs Bin 2)
plot([x_positions(1), x_positions(1), x_positions(2), x_positions(2)], [y_bracket_12, y_bracket_12+10, y_bracket_12+10, y_bracket_12], 'k', 'LineWidth', 1.5);
text(mean(x_positions(1:2)), y_bracket_12 + 40, sprintf('p = %.3f', p12), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Second comparison (Bin 2 vs Bin 3)
plot([x_positions(2), x_positions(2), x_positions(3), x_positions(3)], [y_bracket_23, y_bracket_23+10, y_bracket_23+10, y_bracket_23], 'k', 'LineWidth', 1.5);
text(mean(x_positions(2:3)), y_bracket_23 + 40, sprintf('p = %.3f', p23), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
hold on
% Third comparison (Bin 1 vs Bin 3)
plot([x_positions(1), x_positions(1), x_positions(3), x_positions(3)], [y_bracket_13, y_bracket_13+10, y_bracket_13+10, y_bracket_13], 'k', 'LineWidth', 1.5);
text(mean([x_positions(1), x_positions(3)]), y_bracket_13 + 40, sprintf('p = %.3f', p13), 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Times New Roman');
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
ylim([0 0.6]); yticks(0:0.2:1);

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










