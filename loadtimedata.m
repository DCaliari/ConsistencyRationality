
[~, choiceDelayed, ~] = xlsread('Dataset Creation - Sep','Delayed Plans Choices','A2:Y146');

choice_D_MAIN = choiceDelayed(:,1:11); % Main 11 sets - Binary, Ternary, Quaternary

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
        if choice_D_MAIN{i,j} == 'x'
        Data_DC_exp(i+10*(i-1)-1+j,2) = 1;
        elseif choice_D_MAIN{i,j} == 'y'
        Data_DC_exp(i+10*(i-1)-1+j,3) = 1;
        elseif choice_D_MAIN{i,j} == 'z'
        Data_DC_exp(i+10*(i-1)-1+j,4) = 1;
        elseif choice_D_MAIN{i,j} == 'w'
        Data_DC_exp(i+10*(i-1)-1+j,5) = 1;
        end
    end
end

Data = zeros(N_subjects,Main_sets);
for i=1:N_subjects
    for j=1:Main_sets
        if choice_D_MAIN{i,j} == 'x'
        Data(i,j)=1;
        elseif choice_D_MAIN{i,j} == 'y'
        Data(i,j)=2;
        elseif choice_D_MAIN{i,j} == 'z'
        Data(i,j)=3;
        elseif choice_D_MAIN{i,j} == 'w'
        Data(i,j)=4;
        end
    end
end

Data = Data';

%% Description of the delayed payment plans

% X1, X2, X3, X4: are related to the numbers 1,2,3,4 in the matrix data.

X1 = repmat([160 0 0 0 0 ], Main_sets,1);
X2 = repmat([110 50 25 0 0 ], Main_sets,1);
X3 = repmat([50 50 50 50 0 ], Main_sets,1);
X4 = repmat([0 15 40 170 0], Main_sets,1);

T = repmat([0 3 6 9 12], Main_sets,1);

Ns =N_subjects;

%% Create four dummy vectors if alternatives are in or out

D1 = double(Sets(:,1)==0);
D2 = double(Sets(:,2)==0);
D3 = double(Sets(:,3)==0);
D4 = double(Sets(:,4)==0);
