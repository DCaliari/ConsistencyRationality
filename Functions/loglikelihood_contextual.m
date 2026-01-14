function LL = loglikelihood_contextual(temp, D, Ns, Data, lambda)

a=temp(:,1);
b=temp(:,2);
c=temp(:,3);
d=temp(:,4);

D1=D(:,1);
D2=D(:,2);
D3=D(:,3);
D4=D(:,4);

temp_min = min(temp,[],'all');
temp_max = max(temp,[],'all');

a = (a - temp_min)./(temp_max - temp_min);
b = (b - temp_min)./(temp_max - temp_min);
c = (c - temp_min)./(temp_max - temp_min);
d = (d - temp_min)./(temp_max - temp_min);

% Compute probability of choosing lottery X (0) over lottery  Y (1)

q1 = (exp((a).* lambda).*D1)./...
    (exp((a).* lambda).*D1+ ...
    exp(( b).* lambda).*D2+ ...
    exp(( c).* lambda).*D3+ ...
    exp(( d).* lambda).*D4);

q2 = (exp((b).* lambda).*D2)./...
    (exp((a).* lambda).*D1+ ...
    exp(( b).* lambda).*D2+ ...
    exp(( c).* lambda).*D3+ ...
    exp(( d).* lambda).*D4);

q3 = (exp((c).* lambda).*D3)./...
    (exp((a).* lambda).*D1+ ...
    exp(( b).* lambda).*D2+ ...
    exp(( c).* lambda).*D3+ ...
    exp(( d).* lambda).*D4);

q4 = (exp((d).* lambda).*D4)./...
    (exp((a).* lambda).*D1+ ...
    exp(( b).* lambda).*D2+ ...
    exp(( c).* lambda).*D3+ ...
    exp(( d).* lambda).*D4);

q1 = repmat(q1,1,Ns);
q2 = repmat(q2,1,Ns);
q3 = repmat(q3,1,Ns);
q4 = repmat(q4,1,Ns);

% Create indicators for each case
i1 = Data == 1 ;
i2 = Data == 2 ; 
i3 = Data == 3 ;
i4 = Data == 4 ;

LL1= (log(q1).*i1);
LL1(isnan(LL1))=0;

LL2= (log(q2).*i2);
LL2(isnan(LL2))=0;

LL3= (log(q3).*i3);
LL3(isnan(LL3))=0;

LL4= (log(q4).*i4);
LL4(isnan(LL4))=0;

LL = LL1 + LL2 + LL3 + LL4;

end