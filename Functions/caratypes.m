function [r,PERcara] = caratypes(X,P,cara,Na)

for i=1:length(cara)
    for j=1:Na
        ucara(i,j) = sum( P(j,:).* (1 - exp(- X(j,:).*cara(i)))./cara(i) );
    end
end


for i=1:length(cara)
    [~,idx]=sort(ucara(i,:),'descend');
    r(i,:) = idx;
end

[caratypes, ~, ~] = unique(r, 'rows');

PERcara = caratypes;

end