function [r,PERcrra] = crratypes(X,P,crra,Na)

for i=1:length(crra)
    for j=1:Na
        ucrra(i,j) = sum( P(j,:).*(X(j,:).^(1-crra(i)))./(1-crra(i)) );
    end
end


for i=1:length(crra)
    [~,idx]=sort(ucrra(i,:),'descend');
    r(i,:) = idx;
end

[crratypes, ~, ~] = unique(r, 'rows');

PERcrra = crratypes;

end