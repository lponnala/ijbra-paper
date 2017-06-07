clear; % clc;

% -------------------------------
% Do rare codons really cluster?
% -------------------------------
load EcoliGBK.mat;
load ecoli_ta_am.mat;

% ---- make selection of gene-set here ----
load EcoliPhysRoleGenes.mat; selectgenes=Name;

numgenes=length(selectgenes);
rc_wt_limit=0.1; % zones will be centred on codons whose wt is more than this
stepsize=1; 
ntctm=100; % number of top clusters to merge
numinitclust=[]; % number of initial clusters found (before merge-processing)
Laa=[]; X={}; H=[]; Nc=[]; Name={}; SignifLambda={}; SignifCz={}; SignifNz={}; SignifDz={}; SS={};
count=0;
for i=1:length(S)
    ind=strmatch(S(i).name,selectgenes,'exact');
    if isempty(ind), continue, end
    if length(ind)>1, error('multiple matches for %s\n',S(i).name); end        
    seq=S(i).sequence(1:end-3); % ignore the last (stop) codon
    if rem(length(seq),3)~=0, continue, end
    numcodons=length(seq)/3;
    count=count+1;
    seq(find(seq=='T'))='U'; % convert DNA to RNA
    x=nan(1,numcodons);
    for j=1:numcodons 
        codon=seq(3*(j-1)+1:3*j); x(j)=1/ta.(codon);
    end
    if ~isempty(find(isnan(x))), error('NaNs remain in the x vector!'); end
    if (length(find(x>rc_wt_limit))<2), continue, end % too few rare codons!
    fprintf('\nProcessing %d of %d: i=%d, name=%s, numcodons=%d, numrare = %d...\n',count,numgenes,i,S(i).name,numcodons,length(find(x>rc_wt_limit)));    
    tic
    [Lambda,Start,Stop,Cz,Nz]=kscanstat(x,rc_wt_limit,stepsize);
    Laa=[Laa; length(x)]; X=[X; x]; numinitclust=[numinitclust; length(Lambda)];
    [sLambda,sStart,sStop,sCz,sNz]=mergeclust(Lambda(1:min(ntctm,length(Lambda))),Start(1:min(ntctm,length(Lambda))),Stop(1:min(ntctm,length(Lambda))),Cz(1:min(ntctm,length(Lambda))),Nz(1:min(ntctm,length(Lambda))));
    if isempty(sLambda), error('no clusters left after merge-processing!'), else fprintf('number of clusters left after merge-processing = %d\n',length(sLambda)); end
    ss=[sStart'; sStop'];
    H=[H; 1]; Nc=[Nc; length(sLambda)]; Name=[Name; S(i).name]; SS=[SS; ss];
    SignifLambda=[SignifLambda; sLambda']; SignifCz=[SignifCz; sCz']; SignifNz=[SignifNz; sNz']; sDz=sCz./sNz; SignifDz=[SignifDz; sDz'];
    toc
    save Wc_am_Laa.mat Laa; % Laa excludes the last (stop) codon
    save Wc_am_X.mat X;
    save Wc_am_Name.mat Name;
    save Wc_am_H.mat H;
    save Wc_am_Nc.mat Nc;
    save Wc_am_SignifLambda.mat SignifLambda;
    save Wc_am_SignifCz.mat SignifCz;
    save Wc_am_SignifNz.mat SignifNz;
    save Wc_am_SignifDz.mat SignifDz;
    save Wc_am_SS.mat SS;
    save Wc_am_numinitclust.mat numinitclust;
end
% -------------------------------
hist(numinitclust), xlabel('number of initial clusters'), ylabel('number of genes')
