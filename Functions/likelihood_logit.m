function L = likelihood_logit(x, Data, Na, Main_sets, Sets, ind)

x=x';

for j=1:Main_sets
for i=1:Na
P(j,i) = x(i)/sum(x*Sets(j,:)','all');
end
end

P = P.*Sets;

ChoiceInd = zeros(Main_sets,Na);
for i=1:Main_sets
ChoiceInd(i,Data(i,ind)) = 1;
end

L = P.*ChoiceInd;
L = L(L~=0);
L = -prod(L);

end