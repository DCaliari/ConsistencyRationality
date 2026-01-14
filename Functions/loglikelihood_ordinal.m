function LL = loglikelihood_ordinal(x,PER,Data,Na,Main_sets,Sets,ind)

x=x';

MaxChoices = zeros(Main_sets,Na,size(PER,1));

for i=1: size(PER,1)
    temp = PER(i,:);
    MaxChoices(:,:,i) = Sets;
    for j=1:Main_sets
        if MaxChoices(j,temp(1),i)==1
           MaxChoices(j,[temp(2) temp(3) temp(4)], i) = 0;
        
        elseif MaxChoices(j,temp(1),i)==0 && MaxChoices(j,temp(2),i)==1
           MaxChoices(j,[temp(3) temp(4)], i) = 0;
        
        elseif MaxChoices(j,temp(1),i)==0 && MaxChoices(j,temp(2),i)==0 
           MaxChoices(j,[temp(4)], i) = 0;
        end
    end
end

P = zeros(Main_sets,Na,size(PER,1));

for i=1:size(PER,1)
P(:,:,i) = x(i).*MaxChoices(:,:,i);
end

P = sum(P,3);

ChoiceInd = zeros(Main_sets,Na);
for i=1:Main_sets
ChoiceInd(i,Data(i,ind)) = 1;
end

L = P.*ChoiceInd;
L = L(L~=0);
LL = - sum(log(L));

end