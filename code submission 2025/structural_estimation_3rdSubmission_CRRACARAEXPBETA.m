%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TIME Preferences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

loadtimedata

%% Exponential discounting

sr=10;
sp=5;

time=linspace(0,1/9,sr);
prec=linspace(1e-6,20,sp);

V = {time, prec};
C = cell(1,numel(V));
[C{:}] = ndgrid(V{:});
C = cellfun(@(X) reshape(X,[],1),C,'UniformOutput',false);
C = horzcat(C{:});
p_par = C';

lc = size(p_par,2);

pi_par_old=ones(lc,1)/lc;
F=zeros(Ns,lc);

parfor t=1:lc

% Recover parameters
omega        = p_par(1,t); % Risk Aversion Coefficient
lambda       = p_par(2,t); % Precision Parameter
% kappa      = p_par(3,t);
%
% Transform model parameters to restrict their domains
% lambda = exp(LN_lambda);
%

DU1 = sum(((1./(1+omega)).^T).*X1,2);
DU2 = sum(((1./(1+omega)).^T).*X2,2);
DU3 = sum(((1./(1+omega)).^T).*X3,2);
DU4 = sum(((1./(1+omega)).^T).*X4,2);

temp = [DU1, DU2, DU3, DU4];

D = [D1, D2, D3, D4];

LL = loglikelihood_contextual(temp, D, Ns, Data, lambda);

F(:,t)=exp(sum(LL));

end

tol_par = 01e-08;
par_dist=1; it=1; 
%
%
while   par_dist>tol_par
%    
[a, pi_par_new] = E_fun_nonpar_simplified(pi_par_old,Ns,F); % This is the E-STEP in the EM algorithm
%,
% norm(pi_par_old-pi_par_new,Inf)
par_dist=norm(pi_par_old-pi_par_new)
pi_par_old=pi_par_new;
%
it = it+1;
% save pi_par_new pi_par_new
end


% Using Bayes Theorem we recover the individual posterior distribution

for i=1:Ns

ind_par = a(i,:)';

for j=1:sr
    t = find(p_par(1,:)==time(j));
    t = ind_par(t);
    time_marg_ind(j) = sum(t);
end
time_marg_ind = time_marg_ind';

for j=1:sp
    t = find(p_par(2,:)==prec(j));
    t = ind_par(t);
    prec_marg_ind(j) = sum(t);
end
prec_marg_ind = prec_marg_ind';

time_ind(i)=time*time_marg_ind;
prec_ind(i)=prec*prec_marg_ind;

clearvars time_marg_ind prec_marg_ind
end

time_ind_exp = time_ind;
prec_ind_exp = prec_ind;

clearvars -except time_ind_exp prec_ind_exp

%% HYPERBOLIC DISCOUNTING 

loadtimedata

sr=10;
sp=5;

time=linspace(0,1/9,sr);
prec=linspace(1e-6,20,sp);

V = {time, prec};
C = cell(1,numel(V));
[C{:}] = ndgrid(V{:});
C = cellfun(@(X) reshape(X,[],1),C,'UniformOutput',false);
C = horzcat(C{:});
p_par = C';

lc = size(p_par,2);

pi_par_old=ones(lc,1)/lc;
F=zeros(Ns,lc);

parfor t=1:lc

% Recover parameters
omega        = p_par(1,t); % Risk Aversion Coefficient
lambda    = p_par(2,t); % Precision Parameter
% kappa      = p_par(3,t);
%
% Transform model parameters to restrict their domains
% lambda = exp(LN_lambda);
%

DU1 = sum(((1./(1+omega.*T))).*X1,2);
DU2 = sum(((1./(1+omega.*T))).*X2,2);
DU3 = sum(((1./(1+omega.*T))).*X3,2);
DU4 = sum(((1./(1+omega.*T))).*X4,2);

temp = [DU1, DU2, DU3, DU4];

D = [D1, D2, D3, D4];

LL = loglikelihood_contextual(temp, D, Ns, Data, lambda);

F(:,t)=exp(sum(LL));

end

tol_par = 01e-08;
par_dist=1; it=1; 
%
%
while   par_dist>tol_par
%    
[a, pi_par_new] = E_fun_nonpar_simplified(pi_par_old,Ns,F); % This is the E-STEP in the EM algorithm
%,
% norm(pi_par_old-pi_par_new,Inf)
par_dist=norm(pi_par_old-pi_par_new)
pi_par_old=pi_par_new;
%
it = it+1;
% save pi_par_new pi_par_new
end


for i=1:Ns

ind_par = a(i,:)';

for j=1:sr
    t = find(p_par(1,:)==time(j));
    t = ind_par(t);
    time_marg_ind(j) = sum(t);
end
time_marg_ind = time_marg_ind';

for j=1:sp
    t = find(p_par(2,:)==prec(j));
    t = ind_par(t);
    prec_marg_ind(j) = sum(t);
end
prec_marg_ind = prec_marg_ind';

time_ind(i)=time*time_marg_ind;
prec_ind(i)=prec*prec_marg_ind;

clearvars time_marg_ind prec_marg_ind
end


time_ind_hyp = time_ind;
prec_ind_hyp = prec_ind;

clearvars -except time_ind_exp prec_ind_exp time_ind_hyp prec_ind_hyp



%% BETA-DELTA DISCOUNTING 

loadtimedata


sr=10;
sp=5;
sb=3;

beta=linspace(0.8,1,sb);
time=linspace(0,1/9,sr);
prec=linspace(1e-6,20,sp);

V = {time, prec, beta};
C = cell(1,numel(V));
[C{:}] = ndgrid(V{:});
C = cellfun(@(X) reshape(X,[],1),C,'UniformOutput',false);
C = horzcat(C{:});
p_par = C';

lc = size(p_par,2);

pi_par_old=ones(lc,1)/lc;
F=zeros(Ns,lc);

parfor t=1:lc

% Recover parameters
omega        = p_par(1,t); % Time Coefficient
lambda    = p_par(2,t); % Precision Parameter
beta = p_par(3,t); % Present Bias
% kappa      = p_par(3,t);
%
% Transform model parameters to restrict their domains
% lambda = exp(LN_lambda);
%

DU1 =  ((1./(1+omega)).^T(:,1)).*X1(:,1) +  sum( beta.*((1./(1+omega)).^T(:,[2:5])).*X1(:,[2:5]),2);
DU2 =  ((1./(1+omega)).^T(:,1)).*X2(:,1) +  sum( beta.*((1./(1+omega)).^T(:,[2:5])).*X2(:,[2:5]),2);
DU3 =  ((1./(1+omega)).^T(:,1)).*X3(:,1) +  sum( beta.*((1./(1+omega)).^T(:,[2:5])).*X3(:,[2:5]),2);
DU4 =  ((1./(1+omega)).^T(:,1)).*X4(:,1) +  sum( beta.*((1./(1+omega)).^T(:,[2:5])).*X4(:,[2:5]),2);

temp = [DU1, DU2, DU3, DU4];

D = [D1, D2, D3, D4];

LL = loglikelihood_contextual(temp, D, Ns, Data, lambda);

F(:,t)=exp(sum(LL));

end

tol_par = 01e-08;
par_dist=1; it=1; 
%
%
while   par_dist>tol_par
%    
[a, pi_par_new] = E_fun_nonpar_simplified(pi_par_old,Ns,F); % This is the E-STEP in the EM algorithm
%,
% norm(pi_par_old-pi_par_new,Inf)
par_dist=norm(pi_par_old-pi_par_new)
pi_par_old=pi_par_new;
%
it = it+1;
% save pi_par_new pi_par_new
end



for i=1:Ns

ind_par = a(i,:)';

for j=1:sr
    t = find(p_par(1,:)==time(j));
    t = ind_par(t);
    time_marg_ind(1,j) = sum(t);
end
time_marg_ind = time_marg_ind';

for j=1:sp
    t = find(p_par(2,:)==prec(j));
    t = ind_par(t);
    prec_marg_ind(1,j) = sum(t);
end
prec_marg_ind = prec_marg_ind';

for j=1:sb
    t = find(p_par(3,:)==beta(j));
    t = ind_par(t);
    beta_marg_ind(1,j) = sum(t);
end
beta_marg_ind = beta_marg_ind';

time_ind(i)=time*time_marg_ind;
prec_ind(i)=prec*prec_marg_ind;
beta_ind(i)=beta*beta_marg_ind;

clearvars time_marg_ind prec_marg_ind beta_marg_ind
end


time_ind_betadelta = time_ind;
prec_ind_betadelta = prec_ind;

clearvars -except time_ind_betadelta prec_ind_betadelta time_ind_exp prec_ind_exp time_ind_hyp prec_ind_hyp beta_ind




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RISK Preferences %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CRRA 

loadriskdata

sr = 8;
sp = 5;

% risk=linspace(-0.1,2.4,sr);
risk=linspace(0,2.4,sr);
prec=exp(linspace(1e-6,4,sp));

V = {risk, prec};
C = cell(1,numel(V));
[C{:}] = ndgrid(V{:});
C = cellfun(@(X) reshape(X,[],1),C,'UniformOutput',false);
C = horzcat(C{:});
p_par = C';

lc = size(p_par,2);

pi_par_old=ones(lc,1)/lc;
F=zeros(Ns,lc);
%
parfor t=1:lc

% Recover parameters
r         = p_par(1,t); % Risk Aversion Coefficient
lambda    = p_par(2,t); % Precision Parameter
%
% Transform model parameters to restrict their domains
% lambda = exp(LN_lambda);
%

% Compute expected utility and then standardization of the utility scale

EUX1 = sum(( (X1).^(1-r)./(1-r) ).*(P1),2);
EUX2 = sum(( (X2).^(1-r)./(1-r) ).*(P2),2);
EUX3 = sum(( (X3).^(1-r)./(1-r) ).*(P3),2);
EUX4 = sum(( (X4).^(1-r)./(1-r) ).*(P4),2);

temp = [EUX1, EUX2, EUX3, EUX4];

D = [D1, D2, D3, D4];

LL = loglikelihood_contextual(temp, D, Ns, Data, lambda);

F(:,t)=exp(sum(LL));


end


tol_par =  01e-08;
par_dist=1; it=1; 
%
%
while   par_dist>tol_par
%    
[a, pi_par_new] = E_fun_nonpar_simplified(pi_par_old,Ns,F); % This is the E-STEP in the EM algorithm
%,
% norm(pi_par_old-pi_par_new,Inf)
par_dist=norm(pi_par_old-pi_par_new)
pi_par_old=pi_par_new;
%
it = it+1;
% save pi_par_new pi_par_new
end


for i=1:Ns

ind_par = a(i,:)';

for j=1:sr
    t = find(p_par(1,:)==risk(j));
    t = ind_par(t);
    risk_marg_ind(j) = sum(t);
end
risk_marg_ind = risk_marg_ind';

for j=1:sp
    t = find(p_par(2,:)==prec(j));
    t = ind_par(t);
    prec_marg_ind(j) = sum(t);
end
prec_marg_ind = prec_marg_ind';

risk_ind(i)=risk*risk_marg_ind;
prec_ind(i)=prec*prec_marg_ind;

clearvars risk_marg_ind prec_marg_ind
end

risk_ind_crra = risk_ind;
prec_ind_crra = prec_ind;

clearvars -except time_ind_betadelta prec_ind_betadelta time_ind_exp prec_ind_exp time_ind_hyp prec_ind_hyp beta_ind risk_ind_crra prec_ind_crra

%% CARA

loadriskdata

sr = 8;
sp = 5;

% risk=linspace(-0.1,2.4,sr);
risk=linspace(-0.1,5,sr);
prec=exp(linspace(1e-6,4,sp));

V = {risk, prec};
C = cell(1,numel(V));
[C{:}] = ndgrid(V{:});
C = cellfun(@(X) reshape(X,[],1),C,'UniformOutput',false);
C = horzcat(C{:});
p_par = C';

lc = size(p_par,2);

pi_par_old=ones(lc,1)/lc;
F=zeros(Ns,lc);
%
parfor t=1:lc

% Recover parameters
r         = p_par(1,t); % Risk Aversion Coefficient
lambda    = p_par(2,t); % Precision Parameter
%
% Transform model parameters to restrict their domains
% lambda = exp(LN_lambda);
%

% Compute expected utility and then standardization of the utility scale

EUX1 = sum(( (1 - exp(- r.*X1) )./(r) ).*(P1),2);
EUX2 = sum(( (1 - exp(- r.*X2) )./(r) ).*(P2),2);
EUX3 = sum(( (1 - exp(- r.*X3) )./(r) ).*(P3),2);
EUX4 = sum(( (1 - exp(- r.*X4) )./(r) ).*(P4),2);

temp = [EUX1, EUX2, EUX3, EUX4];

D = [D1, D2, D3, D4];

LL = loglikelihood_contextual(temp, D, Ns, Data, lambda);

F(:,t)=exp(sum(LL));


end


tol_par =  01e-08;
par_dist=1; it=1; 
%
%
while   par_dist>tol_par
%    
[a, pi_par_new] = E_fun_nonpar_simplified(pi_par_old,Ns,F); % This is the E-STEP in the EM algorithm
%,
% norm(pi_par_old-pi_par_new,Inf)
par_dist=norm(pi_par_old-pi_par_new)
pi_par_old=pi_par_new;
%
it = it+1;
% save pi_par_new pi_par_new
end


for i=1:Ns

ind_par = a(i,:)';

for j=1:sr
    t = find(p_par(1,:)==risk(j));
    t = ind_par(t);
    risk_marg_ind(j) = sum(t);
end
risk_marg_ind = risk_marg_ind';

for j=1:sp
    t = find(p_par(2,:)==prec(j));
    t = ind_par(t);
    prec_marg_ind(j) = sum(t);
end
prec_marg_ind = prec_marg_ind';

risk_ind(i)=risk*risk_marg_ind;
prec_ind(i)=prec*prec_marg_ind;

clearvars risk_marg_ind prec_marg_ind
end

risk_ind_cara = risk_ind;
prec_ind_cara = prec_ind;

clearvars -except time_ind_betadelta prec_ind_betadelta time_ind_exp prec_ind_exp time_ind_hyp prec_ind_hyp beta_ind ...
    risk_ind_crra prec_ind_crra risk_ind_cara prec_ind_cara



%% CERTAINTY EQUIVALENT - CRRA

loadriskdata

sr = 8;
sp = 5;

% risk=linspace(-0.1,2.4,sr);
risk=linspace(0,2.4,sr);
prec=exp(linspace(1e-6,4,sp));

V = {risk, prec};
C = cell(1,numel(V));
[C{:}] = ndgrid(V{:});
C = cellfun(@(X) reshape(X,[],1),C,'UniformOutput',false);
C = horzcat(C{:});
p_par = C';

lc = size(p_par,2);

pi_par_old=ones(lc,1)/lc;
F=zeros(Ns,lc);
%
parfor t=1:lc

% Recover parameters
r         = p_par(1,t); % Risk Aversion Coefficient
lambda    = p_par(2,t); % Precision Parameter
%
% Transform model parameters to restrict their domains
% lambda = exp(LN_lambda);
%

% Compute expected utility and then standardization of the utility scale

EUX1 = sum(( (X1).^(1-r)./(1-r) ).*(P1),2);
EUX2 = sum(( (X2).^(1-r)./(1-r) ).*(P2),2);
EUX3 = sum(( (X3).^(1-r)./(1-r) ).*(P3),2);
EUX4 = sum(( (X4).^(1-r)./(1-r) ).*(P4),2);

CEU1 = (EUX1.*(1-r)).^(1/(1-r));
CEU2 = (EUX2.*(1-r)).^(1/(1-r));
CEU3 = (EUX3.*(1-r)).^(1/(1-r));
CEU4 = (EUX4.*(1-r)).^(1/(1-r));

temp = [CEU1, CEU2, CEU3, CEU4];

D = [D1, D2, D3, D4];

LL = loglikelihood_contextual(temp, D, Ns, Data, lambda);

F(:,t)=exp(sum(LL));


end


tol_par =  01e-08;
par_dist=1; it=1; 
%
%
while   par_dist>tol_par
%    
[a, pi_par_new] = E_fun_nonpar_simplified(pi_par_old,Ns,F); % This is the E-STEP in the EM algorithm
%,
% norm(pi_par_old-pi_par_new,Inf)
par_dist=norm(pi_par_old-pi_par_new)
pi_par_old=pi_par_new;
%
it = it+1;
% save pi_par_new pi_par_new
end



for i=1:Ns

ind_par = a(i,:)';

for j=1:sr
    t = find(p_par(1,:)==risk(j));
    t = ind_par(t);
    risk_marg_ind(j) = sum(t);
end
risk_marg_ind = risk_marg_ind';

for j=1:sp
    t = find(p_par(2,:)==prec(j));
    t = ind_par(t);
    prec_marg_ind(j) = sum(t);
end
prec_marg_ind = prec_marg_ind';

risk_ind(i)=risk*risk_marg_ind;
prec_ind(i)=prec*prec_marg_ind;

clearvars risk_marg_ind prec_marg_ind
end

risk_ind_crrace = risk_ind;
prec_ind_crrace = prec_ind;

clearvars -except time_ind_betadelta prec_ind_betadelta time_ind_exp prec_ind_exp time_ind_hyp prec_ind_hyp beta_ind ...
    risk_ind_crra prec_ind_crra risk_ind_cara prec_ind_cara risk_ind_crrace prec_ind_crrace


%% CERTAINTY EQUIVALENT - CARA

loadriskdata

sr = 8;
sp = 5;

% risk=linspace(-0.1,2.4,sr);
risk=linspace(-0.1,5,sr);
prec=exp(linspace(1e-6,4,sp));

V = {risk, prec};
C = cell(1,numel(V));
[C{:}] = ndgrid(V{:});
C = cellfun(@(X) reshape(X,[],1),C,'UniformOutput',false);
C = horzcat(C{:});
p_par = C';

lc = size(p_par,2);

pi_par_old=ones(lc,1)/lc;
F=zeros(Ns,lc);
%
parfor t=1:lc

% Recover parameters
r         = p_par(1,t); % Risk Aversion Coefficient
lambda    = p_par(2,t); % Precision Parameter
%
% Transform model parameters to restrict their domains
% lambda = exp(LN_lambda);
%

% Compute expected utility and then standardization of the utility scale

EUX1 = sum(( (1 - exp(- r.*X1) )./(r) ).*(P1),2);
EUX2 = sum(( (1 - exp(- r.*X2) )./(r) ).*(P2),2);
EUX3 = sum(( (1 - exp(- r.*X3) )./(r) ).*(P3),2);
EUX4 = sum(( (1 - exp(- r.*X4) )./(r) ).*(P4),2);

CEU1 = -log(1-r.*EUX1)./(r);
CEU2 = -log(1-r.*EUX2)./(r);
CEU3 = -log(1-r.*EUX3)./(r);
CEU4 = -log(1-r.*EUX4)./(r);

temp = [CEU1, CEU2, CEU3, CEU4];

D = [D1, D2, D3, D4];

LL = loglikelihood_contextual(temp, D, Ns, Data, lambda);
F(:,t)=exp(sum(LL));


end


tol_par =  01e-08;
par_dist=1; it=1; 
%
%
while   par_dist>tol_par
%    
[a, pi_par_new] = E_fun_nonpar_simplified(pi_par_old,Ns,F); % This is the E-STEP in the EM algorithm
%,
% norm(pi_par_old-pi_par_new,Inf)
par_dist=norm(pi_par_old-pi_par_new)
pi_par_old=pi_par_new;
%
it = it+1;
% save pi_par_new pi_par_new
end


for i=1:Ns

ind_par = a(i,:)';

for j=1:sr
    t = find(p_par(1,:)==risk(j));
    t = ind_par(t);
    risk_marg_ind(j) = sum(t);
end
risk_marg_ind = risk_marg_ind';

for j=1:sp
    t = find(p_par(2,:)==prec(j));
    t = ind_par(t);
    prec_marg_ind(j) = sum(t);
end
prec_marg_ind = prec_marg_ind';

risk_ind(i)=risk*risk_marg_ind;
prec_ind(i)=prec*prec_marg_ind;

clearvars risk_marg_ind prec_marg_ind
end

risk_ind_carace = risk_ind;
prec_ind_carace = prec_ind;

clearvars -except time_ind_betadelta prec_ind_betadelta time_ind_exp prec_ind_exp time_ind_hyp prec_ind_hyp beta_ind ...
    risk_ind_crra prec_ind_crra risk_ind_cara prec_ind_cara risk_ind_crrace prec_ind_crrace risk_ind_carace prec_ind_carace

%% EXPONENTIAL - RANDOM PARAMETERS MODEL

loadtimedata

Sets=Sets==0;

%% EXPONENTIAL PREFERENCES

%% now we set the estimated values of utilities v at the relevant types (TIME)

PER = flipud(perms(1:4)); % ordinal preferences

X1 = [160 0 0 0 0 ];
X2 = [110 50 25 0 0 ];
X3 = [50 50 50 50 0 ];
X4 = [0 15 40 170 0];

X = [X1; X2; X3; X4];

T = [0 3 6 9 12];

T = [T;T;T;T];

clearvars X1 X2 X3 X4

beta = linspace(0.9,1, 10000);
betax = (1-beta)./(beta);

[t,PERexp] = exptypes(X,T,beta,Na);

for i=1:size(PERexp,1)
    for j=1:length(beta)
        tempr(j,1)= isequal(t(j,:),PERexp(i,:));
    end
    omega(i,1) = max(betax(tempr));
end




%% Exponential discounting

sr=10;
sp=8;
st=5;


time=linspace(0,1/9,sr);
prec=linspace(1e-6,5,sp);
trem=linspace(1e-6,.25,st);
% 
V = {time, prec, trem};
C = cell(1,numel(V));
[C{:}] = ndgrid(V{:});
C = cellfun(@(X) reshape(X,[],1),C,'UniformOutput',false);
C = horzcat(C{:});
p_par = C';
lc = size(p_par,2);


pi_par_old=ones(lc,1)/lc;
F=zeros(Ns,lc);

parfor t=1:lc

%
% Recover parameters
r         = p_par(1,t); % Risk Aversion Coefficient
LN_lambda = p_par(2,t); % Precision Parameter
kappa     = p_par(3,t); % Tremble Parameter
%
% Transform model parameters to restrict their domains
lambda = exp(LN_lambda);

% Compute expected utility and then standardization of the utility scale

LL = loglikelihood_rpm(omega, r, lambda, kappa, Data, PERexp);

F(:,t)=exp(sum(LL));


end

tol_par = 01e-08;
par_dist=1; it=1; 
%
%
while   par_dist>tol_par
%    
[a, pi_par_new] = E_fun_nonpar_simplified(pi_par_old,Ns,F); % This is the E-STEP in the EM algorithm
%,
% norm(pi_par_old-pi_par_new,Inf)
par_dist=norm(pi_par_old-pi_par_new)
pi_par_old=pi_par_new;
%
it = it+1;
% save pi_par_new pi_par_new
end



% Using Bayes Theorem we recover the individual posterior distribution

for i=1:Ns

ind_par = a(i,:)';

for j=1:sr
    t = find(p_par(1,:)==time(j));
    t = ind_par(t);
    time_marg_ind(j) = sum(t);
end
time_marg_ind = time_marg_ind';

for j=1:sp
    t = find(p_par(2,:)==prec(j));
    t = ind_par(t);
    prec_marg_ind(j) = sum(t);
end
prec_marg_ind = prec_marg_ind';

time_ind(i)=time*time_marg_ind;
prec_ind(i)=prec*prec_marg_ind;

clearvars time_marg_ind prec_marg_ind
end

time_ind_rpm = time_ind;
prec_ind_rpm_time = prec_ind;

ti_rpm = (1)./(1+time_ind_rpm)<=0.91;
tp_rpm = (1)./(1+time_ind_rpm)>=0.99;


clearvars -except time_ind_betadelta prec_ind_betadelta time_ind_exp prec_ind_exp time_ind_hyp prec_ind_hyp beta_ind ...
    risk_ind_crra prec_ind_crra risk_ind_cara prec_ind_cara risk_ind_crrace prec_ind_crrace risk_ind_carace prec_ind_carace ...
    time_ind_rpm ti_rpm tp_rpm prec_ind_rpm_time

%% CRRA RANDOM PARAMETERS MODEL

loadriskdata

Sets=Sets==0;

%% Calculate the Omega

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

[r,PERcrra] = crratypes(X,P,crra,Na);

for i=1:size(PERcrra,1)
    for j=1:length(crra)
        tempr(j,1)= isequal(r(j,:),PERcrra(i,:));
    end
    omega(i,1) = max(crra(tempr));
end



sr=10;
sp=8;
st=5;

risk=linspace(-0.1,2.4,sr);
prec=linspace(1e-6,2.1,sp);
trem=linspace(1e-6,.25,st);
% 
V = {risk, prec, trem};
C = cell(1,numel(V));
[C{:}] = ndgrid(V{:});
C = cellfun(@(X) reshape(X,[],1),C,'UniformOutput',false);
C = horzcat(C{:});
p_par = C';
lc = size(p_par,2);

pi_par_old=ones(lc,1)/lc;

F=zeros(N_subjects,lc);


parfor t=1:lc

%
% Recover parameters
r         = p_par(1,t); % Risk Aversion Coefficient
LN_lambda = p_par(2,t); % Precision Parameter
kappa     = p_par(3,t); % Tremble Parameter
%
% Transform model parameters to restrict their domains
lambda = exp(LN_lambda);

% Compute expected utility and then standardization of the utility scale

LL = loglikelihood_rpm(omega, r, lambda, kappa, Data, PERcrra);

F(:,t)=exp(sum(LL));


end


tol_par =  01e-08;
par_dist=1; it=1; 
%
%
while   par_dist>tol_par
%    
[a, pi_par_new] = E_fun_nonpar_simplified(pi_par_old,Ns,F); % This is the E-STEP in the EM algorithm
%,
% norm(pi_par_old-pi_par_new,Inf)
par_dist=norm(pi_par_old-pi_par_new)
pi_par_old=pi_par_new;
%
it = it+1;
% save pi_par_new pi_par_new
end


for i=1:Ns

ind_par = a(i,:)';

for j=1:sr
    t = find(p_par(1,:)==risk(j));
    t = ind_par(t);
    risk_marg_ind(j) = sum(t);
end
risk_marg_ind = risk_marg_ind';

for j=1:sp
    t = find(p_par(2,:)==prec(j));
    t = ind_par(t);
    prec_marg_ind(j) = sum(t);
end
prec_marg_ind = prec_marg_ind';

risk_ind(i)=risk*risk_marg_ind;
prec_ind(i)=prec*prec_marg_ind;

clearvars risk_marg_ind prec_marg_ind
end

risk_ind_rpm = risk_ind;
prec_ind_rpm_risk = prec_ind;

ra_rpm = risk_ind_rpm>=2.2;
rl_rpm = risk_ind_rpm<=0.2;


%% grouping participants

ti_exp = time_ind_exp>= 0.09/0.91;
ti_hyp = time_ind_hyp>= 0.09/0.91;
ti_betadelta = time_ind_betadelta>= 0.09/0.91;

tp_exp = time_ind_exp<= 0.01/0.99;
tp_hyp = time_ind_hyp<= 0.01/0.99;
tp_betadelta = time_ind_betadelta<= 0.01/0.99;

ra_carace = risk_ind_carace>=4.5;
ra_crrace = risk_ind_crrace>=2.2;
ra_crra = risk_ind_crra>=2.2;
ra_cara = risk_ind_cara>=4.5;

rl_carace = risk_ind_carace<=0.5;
rl_crrace = risk_ind_crrace<=0.2;
rl_crra = risk_ind_crra<=0.2;
rl_cara = risk_ind_cara<=0.5;

clearvars -except time_ind_betadelta prec_ind_betadelta time_ind_exp prec_ind_exp time_ind_hyp prec_ind_hyp beta_ind ...
    risk_ind_crra prec_ind_crra risk_ind_cara prec_ind_cara risk_ind_crrace prec_ind_crrace risk_ind_carace prec_ind_carace ...
    time_ind_rpm ti_rpm tp_rpm risk_ind_rpm ra_rpm rl_rpm prec_ind_rpm_time prec_ind_rpm_risk ...
    ti_exp ti_hyp ti_betadelta tp_exp tp_hyp tp_betadelta ...
    ra_carace ra_crrace ra_crra ra_cara rl_carace rl_crrace rl_crra rl_cara
    
save 'variables.mat' time_ind_betadelta prec_ind_betadelta time_ind_exp prec_ind_exp time_ind_hyp prec_ind_hyp beta_ind ...
    risk_ind_crra prec_ind_crra risk_ind_cara prec_ind_cara risk_ind_crrace prec_ind_crrace risk_ind_carace prec_ind_carace ...
    time_ind_rpm ti_rpm tp_rpm risk_ind_rpm ra_rpm rl_rpm prec_ind_rpm_time prec_ind_rpm_risk ...
    ti_exp ti_hyp ti_betadelta tp_exp tp_hyp tp_betadelta ...
    ra_carace ra_crrace ra_crra ra_cara rl_carace rl_crrace rl_crra rl_cara

%% comparisons between EXP and HYP and BETADELTA

corrcoef([ti_exp',ti_hyp', ti_betadelta', ti_rpm'])
corrcoef([tp_exp',tp_hyp', tp_betadelta', tp_rpm'])

time_ind = [time_ind_exp',time_ind_hyp',time_ind_betadelta', time_ind_rpm'];

corrcoef(time_ind)

time_ind = [ 1./(1+time_ind(:,1)), 1./(1+time_ind(:,2)), beta_ind'.*(1./(1+time_ind(:,1))), 1./(1+time_ind(:,4)) ];

% Assuming time_ind is already defined and has the necessary calculations
% Example correlation calculation
corr_matrix = corr(time_ind);

% Define the column and row names
column_names = {'Exponential', 'Hyperbolic', 'Beta-Delta', 'Exp/rpm'};
row_names = {'Exponential', 'Hyperbolic', 'Beta-Delta', 'Exp/rpm'};

% Create a table for the correlation coefficients
corr_table = array2table(corr_matrix, 'VariableNames', column_names, 'RowNames', row_names);

% Display the table
disp(corr_table);

% % Write the table to an Excel file
% filename = 'Replication Package\figures and tables\correlation_time.xlsx';
% writetable(corr_table, filename, 'WriteRowNames', true);
% 



prec_ind = [ prec_ind_exp', prec_ind_hyp', prec_ind_betadelta', prec_ind_rpm_time' ];

% Assuming time_ind is already defined and has the necessary calculations
% Example correlation calculation
corr_matrix = corr(prec_ind);

% Define the column and row names
column_names = {'Exponential', 'Hyperbolic', 'Beta-Delta', 'Exp/rpm'};
row_names = {'Exponential', 'Hyperbolic', 'Beta-Delta', 'Exp/rpm'};

% Create a table for the correlation coefficients
corr_table = array2table(corr_matrix, 'VariableNames', column_names, 'RowNames', row_names);

% Display the table
disp(corr_table);

% % Write the table to an Excel file
% filename = 'Replication Package\figures and tables\correlation_time_prec.xlsx';
% writetable(corr_table, filename, 'WriteRowNames', true);


%% comparisons between CRRA and CARA and CE



corrcoef([ra_crra',ra_cara', ra_crrace', ra_carace'])
corrcoef([rl_crra',rl_cara', rl_crrace', rl_carace'])

risk_ind = [risk_ind_crra',risk_ind_cara',risk_ind_crrace',risk_ind_carace', risk_ind_rpm'];

% Assuming time_ind is already defined and has the necessary calculations
% Example correlation calculation
corr_matrix = corr(risk_ind);

% Define the column and row names
column_names = {'crra', 'cara', 'crra/ce', 'cara/ce', 'crra/rpm'};
row_names = {'crra', 'cara', 'crra/ce', 'cara/ce', 'crra/rpm'};

% Create a table for the correlation coefficients
corr_table = array2table(corr_matrix, 'VariableNames', column_names, 'RowNames', row_names);

% Display the table
disp(corr_table);
% 
% % Write the table to an Excel file
% filename = 'Replication Package\figures and tables\correlation_risk.xlsx';
% writetable(corr_table, filename, 'WriteRowNames', true);
% 
% 

prec_ind = [ prec_ind_crra', prec_ind_cara', prec_ind_crrace', prec_ind_carace', prec_ind_rpm_risk' ];

% Assuming time_ind is already defined and has the necessary calculations
% Example correlation calculation
corr_matrix = corr(prec_ind);

% Define the column and row names
column_names = {'crra', 'cara', 'crra/ce', 'cara/ce', 'crra/rpm'};
row_names = {'crra', 'cara', 'crra/ce', 'cara/ce', 'crra/rpm'};
% Create a table for the correlation coefficients
corr_table = array2table(corr_matrix, 'VariableNames', column_names, 'RowNames', row_names);

% Display the table
disp(corr_table);

% % Write the table to an Excel file
% filename = 'Replication Package\figures and tables\correlation_risk_prec.xlsx';
% writetable(corr_table, filename, 'WriteRowNames', true);
