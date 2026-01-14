function [any,irreflexive,asymmetric,acyclic,transitive,linear,intervalorder,semiorder] = binaryrelations(Na)

%% ANY

matrices = HHH(Na);
A = cat(3, matrices);
[n,m,~]=size(A);
a=reshape(A,n,[],1);
b=reshape(a(:),n*m,[])';
c=unique(b,'rows','stable')';
any = reshape(c,n,m,[]);


%% IRREFLEXIVE

for z=1:size(matrices,3)  %transform into irreflexive relation
for i = 1:size(matrices,1)
    for j= 1:size(matrices,2)
        if i==j
            matrices(i,j,z) = 0;
        else
        end
    end
end
end

A = cat(3, matrices);
[n,m,~]=size(A);
a=reshape(A,n,[],1);
b=reshape(a(:),n*m,[])';
c=unique(b,'rows','stable')';
irreflexive = reshape(c,n,m,[]);


%% ASYMMETRIC

for z=1:size(irreflexive,3)  %transform into asymmetric relation
for i = 1:size(irreflexive,1)
    for j= 1:size(irreflexive,2)
        if irreflexive(i,j,z) == irreflexive(j,i,z) && irreflexive(i,j,z) == 1
         irreflexive(i,j,z) = 0;
        else
        end
    end
end
end

A = cat(3, irreflexive);
[n,m,~]=size(A);
a=reshape(A,n,[],1);
b=reshape(a(:),n*m,[])';
c=unique(b,'rows','stable')';
asymmetric = reshape(c,n,m,[]);


%% ACYCLIC

A = cell(1,size(asymmetric,3));
 
for i=1:size(asymmetric,3)
 A{1,i} = digraph(asymmetric(:,:,i)); 
end %set of digraphs 

ttt = zeros(size(asymmetric,3),1);
for i=1:size(asymmetric,3)
    ttt(i) = isdag(A{1,i});
end
ttt=ttt==1;
acyclic = asymmetric(:,:,ttt);

%% TRANSITIVE

A = cell(1,size(acyclic,3));
B = cell(1,size(acyclic,3));
C = zeros(Na,Na,size(acyclic,3));
for i=1:size(acyclic,3)
 A{1,i} = digraph(acyclic(:,:,i));
 A{1,i} = transclosure(A{1,i});
 B{1,i} = adjacency(A{1,i});
 C(:,:,i) = B{1,i};
end %set of digraphs 

A = cat(3, C);
[n,m,~]=size(A);
a=reshape(A,n,[],1);
b=reshape(a(:),n*m,[])';
c=unique(b,'rows','stable')';
transitive = reshape(c,n,m,[]);

%% LINEAR

temp_asy = asymmetric;

for z=1:size(asymmetric,3) %transform into complete relation
for i = 1:size(asymmetric,1)
    for j= 1:size(asymmetric,2)
      r = sum(temp_asy,1);
      q =  sum(r,2);
      if q(1,1,z) < nchoosek(Na,2)
          temp_asy(i,j,z) = 0;
      else
      end
    end
end
end

A = cat(3, temp_asy);
[n,m,~]=size(A);
a=reshape(A,n,[],1);
b=reshape(a(:),n*m,[])';
c=unique(b,'rows','stable')';
transitivecomplete = reshape(c,n,m,[]);

transitivecompleteadj= cell(1,size(transitivecomplete,3));

transitivecompletegraph = cell(1,size(transitivecomplete,3));

for i = 1:size(transitivecomplete,3)
transitivecompleteadj{1,i} = transitivecomplete(:,:,i);
end

for z=1:size(transitivecomplete,3)
transitivecompletegraph{1,z} = digraph(transitivecompleteadj{1,z});
end
   

for z = 1:size(transitivecomplete,3)  %transform into acyclic relation
for i = 1:size(transitivecomplete,1)
   for j= 1:size(transitivecomplete,2)
   if isdag(transitivecompletegraph{1,z}) == 1
      transitivecomplete(i,j,z) = transitivecomplete(i,j,z);
   elseif isdag(transitivecompletegraph{1,z}) == 0
     transitivecomplete(i,j,z) = 0;
   end
   end
end
end

A = cat(3, transitivecomplete);
[n,m,~]=size(A);
a=reshape(A,n,[],1);
b=reshape(a(:),n*m,[])';
c=unique(b,'rows','stable')';
transitivecomplete = reshape(c,n,m,[]);
linear = transitivecomplete(:,:,2:size(transitivecomplete,3));

% STRONG INTERVALITY CONDITION 
% xPy and zPt implies xPt or zPy

STRONGINTERVALITY = ones(size(transitive,3),1);
for i=1:size(transitive,3)
for x=1:4
    for y=1:4
        if transitive(x,y,i)==1
            for z=1:4
                for t=1:4
                    if transitive(z,t,i)==1 && transitive(x,t,i) ==0 && transitive(z,y,i)==0
                        STRONGINTERVALITY(i) = 0;
                    end
                end
            end
        end
    end
end
end

% SEMITRANSITIVITY CONDITION
% xPy and yPz implies xPt or tPz

SEMITRANSITIVITY = ones(size(transitive,3),1);
for i=1:size(transitive,3)
for x=1:4
    for y=1:4
        if transitive(x,y,i)==1
            for z=1:4
                for t=1:4
                    if transitive(y,z,i)==1 && transitive(x,t,i) ==0 && transitive(t,z,i)==0
                        SEMITRANSITIVITY(i) = 0;
                    end
                end
            end
        end
    end
end
end

intervalorder = transitive(:,:,logical(STRONGINTERVALITY));
semiorder = transitive(:,:,logical(SEMITRANSITIVITY.*STRONGINTERVALITY));


end