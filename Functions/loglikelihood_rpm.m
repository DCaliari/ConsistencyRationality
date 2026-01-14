function LL = loglikelihood_rpm(omega, r, lambda, kappa, Data,PER)

for j=1:size(PER,1)
    tempab(j,1) = find(PER(j,:)==1);
    tempab(j,2) = find(PER(j,:)==2);
    ab = sum(tempab(:,1)<tempab(:,2))+1;
    tempac(j,1) = find(PER(j,:)==1);
    tempac(j,2) = find(PER(j,:)==3);
    ac = sum(tempac(:,1)<tempac(:,2))+1;
    tempad(j,1) = find(PER(j,:)==1);
    tempad(j,2) = find(PER(j,:)==4);
    ad = sum(tempad(:,1)<tempad(:,2))+1;
    tempbc(j,1) = find(PER(j,:)==2);
    tempbc(j,2) = find(PER(j,:)==3);
    bc = sum(tempbc(:,1)<tempbc(:,2))+1;
    tempbd(j,1) = find(PER(j,:)==2);
    tempbd(j,2) = find(PER(j,:)==4);
    bd = sum(tempbd(:,1)<tempbd(:,2))+1;
    tempcd(j,1) = find(PER(j,:)==3);
    tempcd(j,2) = find(PER(j,:)==4);
    cd = sum(tempcd(:,1)<tempcd(:,2))+1;
    tempabc(j,1) = find(PER(j,:)==1);
    tempabc(j,2) = find(PER(j,:)==2);
    tempabc(j,3) = find(PER(j,:)==3);
    abc(j,1) = find(min(tempabc(j,:)) == tempabc(j,:));
    tempabd(j,1) = find(PER(j,:)==1);
    tempabd(j,2) = find(PER(j,:)==2);
    tempabd(j,3) = find(PER(j,:)==4);
    abd(j,1) = find(min(tempabd(j,:)) == tempabd(j,:));
    tempacd(j,1) = find(PER(j,:)==1);
    tempacd(j,2) = find(PER(j,:)==3);
    tempacd(j,3) = find(PER(j,:)==4);
    acd(j,1) = find(min(tempacd(j,:)) == tempacd(j,:));
    tempbcd(j,1) = find(PER(j,:)==2);
    tempbcd(j,2) = find(PER(j,:)==3);
    tempbcd(j,3) = find(PER(j,:)==4);
    bcd(j,1) = find(min(tempbcd(j,:)) == tempbcd(j,:));
    tempabcd(j,1) = find(PER(j,:)==1);
    tempabcd(j,2) = find(PER(j,:)==2);
    tempabcd(j,3) = find(PER(j,:)==3);
    tempabcd(j,4) = find(PER(j,:)==4);
    abcd(j,1) = find(min(tempabcd(j,:)) == tempabcd(j,:));
end


% ab = omega(2);
% ac = omega(3);
% ad = omega(4);
% bc = omega(7);
% bd = omega(6);
% cd = omega(5);
% 
% aabc = omega(2);
% babc = omega(7);
% 
% aabd = omega(2);
% babd = omega(6);
% 
% aacd = omega(3);
% cacd = omega(5);
% 
% bbcd = omega(6);
% cbcd = 0;
% 
% aabcd = omega(2);
% babcd = omega(6);
% cabcd = 0;

ab = omega(ab);
ac = omega(ac);
ad = omega(ad);
bc = omega(bc);
bd = omega(bd);
cd = omega(cd);

% Triplet & quadruplet assignments with empty checks
idx = find(abc == 1, 1, 'last');
if ~isempty(idx)
    aabc = omega(idx + 1);
else
    aabc = 0;
end


idx0 = find(abc == 2, 1, 'last');
if ~isempty(idx0)
    babc = omega(idx0 + 1);
else
    babc = 0;
end
idx1 = find(abd == 1, 1, 'last');
if ~isempty(idx1)
    aabd = omega(idx1);
else
    aabd = 0;
end

idx2 = find(abd == 2, 1, 'last');
if ~isempty(idx2)
    babd = omega(idx2);
else
    babd = 0;
end

idx3 = find(acd == 1, 1, 'last');
if ~isempty(idx3)
    aacd = omega(idx3);
else
    aacd = 0;
end

idx4 = find(acd == 3, 1, 'last');
if ~isempty(idx4)
    cacd = omega(idx4);
else
    cacd = 0;
end

idx5 = find(bcd == 2, 1, 'last');
if ~isempty(idx5)
    bbcd = omega(idx5);
else
    bbcd = 0;
end

idx6 = find(bcd == 3, 1, 'last');
if ~isempty(idx6)
    cbcd = omega(idx6);
else
    cbcd = 0;
end

idx7 = find(abcd == 1, 1, 'last');
if ~isempty(idx7)
    aabcd = omega(idx7);
else
    aabcd = 0;
end

idx8 = find(abcd == 2, 1, 'last');
if ~isempty(idx8)
    babcd = omega(idx8);
else
    babcd = 0;
end

idx9 = find(abcd == 3, 1, 'last');
if ~isempty(idx9)
    cabcd = omega(idx9);
else
    cabcd = 0;
end



%
% Compute probability of choosing lottery X (0) over lottery  Y (1)

pab = (1-2.*kappa).*( exp ( (r).*lambda ) )./(exp ( (r).*lambda ) + exp ( (ab).*lambda )  ) + kappa;
pba = 1- pab;
pac = (1-2.*kappa).*( exp ( (r).*lambda ) )./(exp ( (r).*lambda ) + exp ( (ac).*lambda )  ) + kappa;
pca = 1 - pac;
pad = (1-2.*kappa).*( exp ( (r).*lambda ) )./(exp ( (r).*lambda ) + exp ( (ad).*lambda )  ) + kappa;
pda = 1 - pad;
pbc = (1-2.*kappa).*( exp ( (r).*lambda ) )./(exp ( (r).*lambda ) + exp ( (bc).*lambda )  ) + kappa;
pcb = 1 - pbc;
pbd = (1-2.*kappa).*( exp ( (r).*lambda ) )./(exp ( (r).*lambda ) + exp ( (bd).*lambda )  ) + kappa;
pdb = 1 - pbd;
pcd = (1-2.*kappa).*( exp ( (r).*lambda ) )./(exp ( (r).*lambda ) + exp ( (cd).*lambda )  ) + kappa;
pdc = 1 - pcd;
paabc = (1-3.*kappa).*( exp ( (r).*lambda ) )./(exp ( (r).*lambda ) + exp ( (aabc).*lambda ) + exp ( (babc).*lambda ) ) + kappa;
pbabc = (1-3.*kappa).*( exp ( (babc).*lambda ) )./(exp ( (r).*lambda ) + exp ( (aabc).*lambda ) + exp ( (babc).*lambda ) ) + kappa;
pcabc = 1 - paabc - pbabc;
paabd = (1-3.*kappa).*( exp ( (r).*lambda ) )./(exp ( (r).*lambda ) + exp ( (aabd).*lambda ) + exp ( (babd).*lambda ) ) + kappa;
pbabd = (1-3.*kappa).*( exp ( (babd).*lambda ) )./(exp ( (r).*lambda ) + exp ( (aabd).*lambda ) + exp ( (babd).*lambda ) ) + kappa;
pdabd = 1 - paabd - pbabd;
paacd = (1-3.*kappa).*( exp ( (r).*lambda ) )./(exp ( (r).*lambda ) + exp ( (aacd).*lambda ) + exp ( (cacd).*lambda ) ) + kappa;
pcacd = (1-3.*kappa).*( exp ( (cacd).*lambda ) )./(exp ( (r).*lambda ) + exp ( (aacd).*lambda ) + exp ( (cacd).*lambda ) ) + kappa;
pdacd = 1 - paacd - pcacd;
pbbcd = (1-3.*kappa).*( exp ( (r).*lambda ) )./(exp ( (r).*lambda ) + exp ( (bbcd).*lambda ) + exp ( (cbcd).*lambda ) ) + kappa;
pcbcd = (1-3.*kappa).*( exp ( (cbcd).*lambda ) )./(exp ( (r).*lambda ) + exp ( (bbcd).*lambda ) + exp ( (cbcd).*lambda ) ) + kappa;
pdbcd = 1 - pbbcd - pcbcd;
paabcd =  (1-4.*kappa).*( exp ( (r).*lambda ) )./(exp ( (r).*lambda ) + exp ( (aabcd).*lambda ) + exp ( (babcd).*lambda ) + exp ( (cabcd).*lambda ) ) + kappa;
pbabcd =  (1-4.*kappa).*( exp ( (babcd).*lambda ) )./(exp ( (r).*lambda ) + exp ( (aabcd).*lambda ) + exp ( (babcd).*lambda ) + exp ( (cabcd).*lambda ) ) + kappa;
pcabcd =  (1-4.*kappa).*( exp ( (cabcd).*lambda ) )./(exp ( (r).*lambda ) + exp ( (aabcd).*lambda ) + exp ( (babcd).*lambda ) + exp ( (cabcd).*lambda ) ) + kappa;
pdabcd =  1- paabcd - pbabcd - pcabcd;

prob = [pab  pba  0  0; 
    pac 0 pca 0;
    pad 0 0 pda;
    0 pbc pcb 0;
    0 pbd 0 pdb;
    0 0 pcd pdc;
    paabc pbabc pcabc 0;
    paabd pbabd 0 pdabd;
    paacd 0 pcacd pdacd;
    0 pbbcd pcbcd pdbcd;
    paabcd pbabcd pcabcd pdabcd];

i1 = Data == 1 ;
i2 = Data == 2 ; 
i3 = Data == 3 ;
i4 = Data == 4 ;
%

LL1= (log(prob(:,1)).*i1);
LL1(isnan(LL1))=0;

LL2= (log(prob(:,2)).*i2);
LL2(isnan(LL2))=0;

LL3= (log(prob(:,3)).*i3);
LL3(isnan(LL3))=0;

LL4= (log(prob(:,4)).*i4);
LL4(isnan(LL4))=0;

LL = LL1 + LL2 + LL3 + LL4;


% Create indicators for each case

end