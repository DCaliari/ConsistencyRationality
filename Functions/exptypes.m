function [t,PERexp] = exptypes(X,T,beta,Na)


for i=1:length(beta)
    for j=1:Na
        uexp(i,j) = sum( X(j,:).*(beta(i)).^(T(j,:)) );
    end
end


for i=1:length(beta)
    [~,idx]=sort(uexp(i,:),'descend');
    t(i,:) = idx;
end

[exptypes, ~, ~] = unique(t, 'rows');

PERexp = exptypes;

end