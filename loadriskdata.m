
[numbers, choiceLotteries, everything] = xlsread('Dataset Creation - Sep','Lotteries Choices','A2:Y146');

choice_L_MAIN = choiceLotteries(:,1:11); % Main 11 sets - Binary, Ternary, Quaternary

N_subjects = 145;
Main_sets = 11;
Na = 4;
Data_DC_exp = zeros(N_subjects*Main_sets,Na+1);

for i=1:N_subjects
Data_DC_exp(i+10*(i-1):10+i+10*(i-1),1) = ones()*i; % creates the line of ID
end

Sets = [1,1,0,0;
    1,0,1,0;
    1,0,0,1;
    0,1,1,0;
    0,1,0,1;
    0,0,1,1;
    1,1,1,0;
    1,1,0,1;
    1,0,1,1;
    0,1,1,1;
    1,1,1,1];

Sets(Sets==0)=-99;
Sets(Sets==1)=0;
Data_DC_exp(:,2:5)= Data_DC_exp(:,2:5) + (repmat(Sets,N_subjects,1));

for i=1:N_subjects
    for j=1:Main_sets
        if choice_L_MAIN{i,j} == 'a'
        Data_DC_exp(i+10*(i-1)-1+j,2) = 1;
        elseif choice_L_MAIN{i,j} == 'b'
        Data_DC_exp(i+10*(i-1)-1+j,3) = 1;
        elseif choice_L_MAIN{i,j} == 'c'
        Data_DC_exp(i+10*(i-1)-1+j,4) = 1;
        elseif choice_L_MAIN{i,j} == 'd'
        Data_DC_exp(i+10*(i-1)-1+j,5) = 1;
        end
    end
end

Data = zeros(N_subjects,Main_sets);
for i=1:N_subjects
    for j=1:Main_sets
        if choice_L_MAIN{i,j} == 'a'
        Data(i,j)=1;
        elseif choice_L_MAIN{i,j} == 'b'
        Data(i,j)=2;
        elseif choice_L_MAIN{i,j} == 'c'
        Data(i,j)=3;
        elseif choice_L_MAIN{i,j} == 'd'
        Data(i,j)=4;
        end
    end
end

Data = Data';

%% Description of the lotteries

% X1, X2, X3, X4: are related to the numbers 1,2,3,4 in the matrix data.

scale =100;

X1 = repmat([50 0.001], Main_sets,1);
X2 = repmat([65 25], Main_sets,1);
X3 = repmat([95 20], Main_sets,1);
X4 = repmat([300 5], Main_sets,1);

X1=X1./scale;
X2=X2./scale;
X3=X3./scale;
X4=X4./scale;

P1 = repmat([1 0], Main_sets,1);
P2 = repmat([0.8 0.2], Main_sets,1);
P3 = repmat([0.5 0.5], Main_sets,1);
P4 = repmat([0.2 0.8], Main_sets,1);

Ns =N_subjects;

%% Create four dummy vectors if alternatives are in or out

D1 = double(Sets(:,1)==0);
D2 = double(Sets(:,2)==0);
D3 = double(Sets(:,3)==0);
D4 = double(Sets(:,4)==0);