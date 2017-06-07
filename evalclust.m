
% //// This script answers the questions ////
% NOTE: the word "cluster" means a statistically significant cluster of
% rare codons, found by the method of [SpScanStatSurvivalData]

% -----------------------------------------------
% mRNA-protein Vs cluster-features
% -----------------------------------------------
clear;
load 'ChurchData.mat';
% columns contain, in order:
% name,start,stop,strand,func,seq,signal,nabd,spint,rnacopies/cell
name=D(:,1); start=cell2mat(D(:,2)); stop=cell2mat(D(:,3)); strand=D(:,4); func=D(:,5); seq=D(:,6); signal=D(:,7);
prot=cell2mat(D(:,8)); % protein (N-abd in molecules/cell)
spint=cell2mat(D(:,9)); % rna "SP Intensity"
rna=cell2mat(D(:,10)); % rna (estimated copies/cell)
ratio=prot./rna; % hist(ratio), title('prot/rna ratio')

load Wc_am_Name.mat; load Wc_am_H.mat; load Wc_am_Nc.mat;
load Wc_am_SignifNz.mat; load Wc_am_SignifCz.mat; load Wc_am_SignifDz.mat;
myname={}; myratio=[]; myprot=[]; myrna=[]; myNc=[]; myNz=[]; myCz=[]; myDz=[];
for i=1:length(name)
    ind=strmatch(name{i},Name,'exact');
    if isempty(ind), continue, end
    if length(ind)>1, error('multiple matches for %s\n',name{i}); end    
    if Nc(ind)==0
        % myNc=[myNc; 0];
        % myNz=[myNz; NaN]; myDz=[myDz; NaN];
    elseif Nc(ind)>1
        myname=[myname; name{i}]; myratio=[myratio; ratio(i)]; myprot=[myprot; prot(i)]; myrna=[myrna; rna(i)]; myNc=[myNc; Nc(ind)]; 
        myDz=[myDz; prod(SignifDz{ind})^(1/Nc(ind))]; 
    else
        myname=[myname; name{i}]; myratio=[myratio; ratio(i)]; myprot=[myprot; prot(i)]; myrna=[myrna; rna(i)]; myNc=[myNc; Nc(ind)]; 
        % myNz=[myNz; SignifNz{ind}]; myCz=[myCz; SignifCz{ind}]; 
        myDz=[myDz; SignifDz{ind}];
    end
end

% compute Nc-level-wise total prot, rna
fprintf('\n---- Computing Nc-level-wise avg(prot),avg(rna) levels ----\n');
Nc_level=[]; total_prot=[]; total_rna=[]; avg_prot=[]; avg_rna=[];
for i=1:max(myNc)
    I=find(myNc==i);
    if ~isempty(I)
        Nc_level=[Nc_level; i]; total_prot=[total_prot; sum(myprot(I))]; total_rna=[total_rna; sum(myrna(I))];  
        avg_prot=[avg_prot; mean(myprot(I))]; avg_rna=[avg_rna; mean(myrna(I))];        
    end
end
% figure, plot(Nc_level,log10(avg_prot),'k*'), xlabel('number of slow-translating clusters'), ylabel('average protein abundance')

% [pp_rho,pp_pval]=partialcorr(total_prot,Nc_level,total_rna); fprintf('partial corr between total_prot and number-of-clusters (while controlling for total_rna) = %f, p-value = %f\n',pp_rho,pp_pval);
[pp_rho,pp_pval]=partialcorr(avg_prot,Nc_level,avg_rna); fprintf('partial corr between avg_prot and number-of-clusters (while controlling for avg_rna) = %f, p-value = %f\n',pp_rho,pp_pval);

[R,P]=corrcoef(avg_prot,Nc_level); fprintf('simple correlation between avg_prot and number-of-clusters = %f, p-value = %f\n',R(1,2),P(1,2));
% pause

% -----------------------------------------------------
% how many genes in each category contain rare codons?
% -----------------------------------------------------
clear;
load Wc_am_Name.mat; load Wc_am_H.mat;
% load genes in each category
enzyme=load('Enzyme.mat');
metab=load('Metabolism.mat');
cellstr=load('CellStructure.mat');
logp=load('LocationOfGeneProducts.mat');
inftr=load('InformationTransfer.mat');
regul=load('Regulation.mat');
transp=load('Transport.mat');
% create blank index-vectors for each category of genes
I_enzyme=[]; I_metab=[]; I_cellstr=[]; I_logp=[]; I_inftr=[]; I_regul=[]; I_transp=[];
for i=1:length(Name)
    % enzyme
    ind=strmatch(Name{i},enzyme.G,'exact'); if length(ind)==1, I_enzyme=[I_enzyme; i]; elseif length(ind)>1, error('multiple matches found for gene %s among enzyme genes\n',Name{i}), end
    % metab
    ind=strmatch(Name{i},metab.G,'exact'); if length(ind)==1, I_metab=[I_metab; i]; elseif length(ind)>1, error('multiple matches found for gene %s among metab genes\n',Name{i}), end
    % cellstr
    ind=strmatch(Name{i},cellstr.G,'exact'); if length(ind)==1, I_cellstr=[I_cellstr; i]; elseif length(ind)>1, error('multiple matches found for gene %s among cellstr genes\n',Name{i}), end    
    % logp
    ind=strmatch(Name{i},logp.G,'exact'); if length(ind)==1, I_logp=[I_logp; i]; elseif length(ind)>1, error('multiple matches found for gene %s among logp genes\n',Name{i}), end
    % inftr
    ind=strmatch(Name{i},inftr.G,'exact'); if length(ind)==1, I_inftr=[I_inftr; i]; elseif length(ind)>1, error('multiple matches found for gene %s among inftr genes\n',Name{i}), end
    % regul
    ind=strmatch(Name{i},regul.G,'exact'); if length(ind)==1, I_regul=[I_regul; i]; elseif length(ind)>1, error('multiple matches found for gene %s among regul genes\n',Name{i}), end
    % transp
    ind=strmatch(Name{i},transp.G,'exact'); if length(ind)==1, I_transp=[I_transp; i]; elseif length(ind)>1, error('multiple matches found for gene %s among transp genes\n',Name{i}), end    
end
fprintf('\n\n');
fprintf('Number of *enzyme* genes examined = %d (%f percent)\n',length(I_enzyme),100*length(I_enzyme)/length(enzyme.G));
fprintf('Number of them having cluster(s) = %d (%f percent)\n\n',length(find(H(I_enzyme))),100*length(find(H(I_enzyme)))/length(I_enzyme));
fprintf('Number of *metab* genes examined = %d (%f percent)\n',length(I_metab),100*length(I_metab)/length(metab.G));
fprintf('Number of them having cluster(s) = %d (%f percent)\n\n',length(find(H(I_metab))),100*length(find(H(I_metab)))/length(I_metab));
fprintf('Number of *cellstr* genes examined = %d (%f percent)\n',length(I_cellstr),100*length(I_cellstr)/length(cellstr.G));
fprintf('Number of them having cluster(s) = %d (%f percent)\n\n',length(find(H(I_cellstr))),100*length(find(H(I_cellstr)))/length(I_cellstr));
fprintf('Number of *logp* genes examined = %d (%f percent)\n',length(I_logp),100*length(I_logp)/length(logp.G));
fprintf('Number of them having cluster(s) = %d (%f percent)\n\n',length(find(H(I_logp))),100*length(find(H(I_logp)))/length(I_logp));
fprintf('Number of *inftr* genes examined = %d (%f percent)\n',length(I_inftr),100*length(I_inftr)/length(inftr.G));
fprintf('Number of them having cluster(s) = %d (%f percent)\n\n',length(find(H(I_inftr))),100*length(find(H(I_inftr)))/length(I_inftr));
fprintf('Number of *regul* genes examined = %d (%f percent)\n',length(I_regul),100*length(I_regul)/length(regul.G));
fprintf('Number of them having cluster(s) = %d (%f percent)\n\n',length(find(H(I_regul))),100*length(find(H(I_regul)))/length(I_regul));
fprintf('Number of *transp* genes examined = %d (%f percent)\n',length(I_transp),100*length(I_transp)/length(transp.G));
fprintf('Number of them having cluster(s) = %d (%f percent)\n\n',length(find(H(I_transp))),100*length(find(H(I_transp)))/length(I_transp));

% --------------------------------------------------
% percentage of genes that have atleast one cluster
% --------------------------------------------------
clear;
load Wc_am_H.mat;
fprintf('percentage of genes that have atleast one cluster = %f\n\n\n',100*length(find(H==1))/length(H));

% --------------------------------------------------
% some cluster-size evaluation metrics
% --------------------------------------------------
clear;
load Wc_am_H.mat; load Wc_am_Nc.mat; load Wc_am_Laa.mat; load Wc_am_SignifCz.mat; load Wc_am_SignifNz.mat; load Wc_am_SignifDz.mat; load Wc_am_Name.mat;
nz=[]; cz=[]; dz=[]; gene={};
for i=1:length(H)
    if (H(i)==1)
        signifNz=SignifNz{i}; signifCz=SignifCz{i}; signifDz=SignifDz{i};
        for j=1:Nc(i)
            nz=[nz; signifNz(j)];
            cz=[cz; signifCz(j)];
            dz=[dz; signifDz(j)];
            gene=[gene; Name{i}];
        end
    end
end
% max cluster size
fprintf('max cluster size = %d\n',max(nz));
N=find(nz==max(nz)); fprintf('number of clusters having max size = %d\n',length(N));
fprintf('gene(s) having maximum-size clusters:\n'); disp(gene(N))
% max cluster content (i.e. sum(wt))
fprintf('max cluster content = %f\n',max(cz));
C=find(cz==max(cz)); fprintf('number of clusters having max content = %d\n',length(C));
fprintf('gene(s) having maximum-content clusters:\n'); disp(gene(C))
% max cluster density
fprintf('max cluster density = %f\n',max(dz));
D=find(dz==max(dz)); fprintf('number of clusters having max density = %d\n',length(D));
fprintf('gene(s) having maximum-density clusters:\n'); disp(gene(D))
fprintf('\n\n');

% --------------------------------------------------
% gap between clusters
% --------------------------------------------------
clear;
load Wc_am_Nc.mat; load Wc_am_SS.mat; load Wc_am_Name.mat;
gap=[]; gene={};
for i=1:length(Nc)
    if (Nc(i)>1)
        ss=SS{i}; start=ss(1,:); stop=ss(2,:);
        [start,I]=sort(start); stop=stop(I);
        for j=2:length(start)
            if (start(j)>stop(j-1))
                gap=[gap; (start(j)-stop(j-1)-1)];
                gene=[gene; Name{i}];
            end
        end
    end
end
figure, hist(gap), title('gap between non-overlapping clusters')
fprintf('Average gap (in codons) between non-overlapping clusters = %f\n\n\n',mean(gap));

% --------------------------------------------------
% percentage of all clusters in short proteins
% --------------------------------------------------
clear;
load Wc_am_Laa.mat; load Wc_am_H.mat; load Wc_am_Nc.mat;
% lothresh=300; I=find(Laa<=lothresh); fprintf('percentage of short (<=%d aa) proteins that have atleast one cluster = %f\n\n\n',lothresh,100*length(find(H(I)))/length(I));
% lothresh=[250 350]; I2=find((Laa>=lothresh(1))&(Laa<=lothresh(2))); fprintf('percentage of short (between %d and %d aa) proteins that have atleast one cluster = %f\n\n\n',lothresh(1),lothresh(2),100*length(find(H(I2)))/length(I2));
% hithresh=500; J=find(Laa>hithresh); fprintf('percentage of long (>%d aa) proteins that have atleast one cluster = %f\n\n\n',hithresh,100*length(find(H(J)))/length(J));
lothresh=300; I=find(Laa<=lothresh); fprintf('percentage of all clusters that occur in short (<=%d aa) proteins = %f\n',lothresh,100*sum(Nc(I))/sum(Nc));
lothresh=[250 350]; I2=find((Laa>=lothresh(1))&(Laa<=lothresh(2))); fprintf('percentage of all clusters that occur in short (between %d and %d aa) proteins = %f\n',lothresh(1),lothresh(2),100*sum(Nc(I2))/sum(Nc));
hithresh=500; J=find(Laa>hithresh); fprintf('percentage of all clusters that occur in long (>%d aa) proteins = %f\n\n\n',hithresh,100*sum(Nc(J))/sum(Nc));

% --------------------------------------------------
% number of proteins having atleast one slow-translating cluster Vs. length
% of protein
% --------------------------------------------------
clear;
load Wc_am_Laa.mat; load Wc_am_H.mat;
I=find(H==1); L=Laa(I); figure, hist(L), title('length of proteins that contain atleast one cluster')
% disp(mean(L))
thresh=300;
fprintf('Out of %d proteins having atleast one cluster, %d have length>=%d\n\n\n',length(L),length(find(L>=thresh)),thresh);

% //// cluster analysis ////
clear;
load Wc_am_H.mat; load Wc_am_SS.mat; load Wc_am_Name.mat;
S=[]; gene={};
for i=1:length(H)
    if (H(i)==1)
        ss=SS{i};
        for j=1:size(ss,2)
            S=[S; [ss(1,j),ss(2,j)]];
            gene=[gene; Name{i}];
        end
    end
end
% what percentage of the clusters are completely inside the first 50 amino acids (150 bases)?
fprintf('percentage of clusters that are completely inside the first 50 amino acids = %f\n',100*length(find(S(:,2)<=50))/size(S,1));
% what percentage of the clusters have atleast one-aa inside the first 50 amino acids (150 bases)?
fprintf('percentage of clusters that have atleast one-aa inside the first 50 amino acids = %f\n',100*length(find(S(:,1)<=50))/size(S,1));
% what percentage of the clusters have atleast 5-aa inside the first 50 amino acids (150 bases)?
fprintf('percentage of clusters that have atleast 5-aa inside the first 50 amino acids and the remaining outside = %f\n',100*length(find((S(:,1)<=46)&(S(:,2)>50)))/size(S,1));
% what percentage of the clusters have atleast 10-aa inside the first 50 amino acids (150 bases)?
fprintf('percentage of clusters that have atleast 10-aa inside the first 50 amino acids and the remaining outside = %f\n',100*length(find((S(:,1)<=41)&(S(:,2)>50)))/size(S,1));
% how many of the clusters are atleast half-inside the first 50 amino acids (150 bases)?
count=0; fcount=0;
for i=1:size(S,1)
    if S(i,2)<=50, fcount=fcount+1; end
    midpointindex=S(i,1)+floor((S(i,2)-S(i,1))/2);
    if (midpointindex<=50), count=count+1; end
end
fprintf('percentage of clusters that are atleast half-inside the first 50 amino acids = %f\n',100*count/size(S,1));
fprintf('percentage of clusters that are completely inside the first 50 amino acids = %f\n',100*fcount/size(S,1));

% --------------------------------------------------
% Plot the histogram of estimated waiting time
% --------------------------------------------------
clear;
load ecoli_ta_am.mat;
F=fields(ta); F(strmatch('UGA',F))=[]; for i=1:length(F), tav(i,1)=ta.(F{i}); end
wt=1./tav;
fprintf('\n\nfit distr to waiting-time: negative log-likelihood (min is best)\n');
fprintf('exp: %f\n',explike(expfit(wt),wt))
fprintf('gamma: %f\n',gamlike(gamfit(wt),wt))
[mu,sig]=normfit(wt); fprintf('normal: %f\n',normlike([mu,sig],wt))
waiting_time=wt; % dfittool
figure, wthist(waiting_time)

% --------------------------------------------------
% Print table of estimated waiting time for each codon
% --------------------------------------------------
clear;
load ecoli_ta_am.mat;
F=fields(ta); F(strmatch('UGA',F))=[];
fprintf('\n\nCodon\tEstWt\n');
for i=1:length(F)
    fprintf('%s\t%f\n',F{i},1/ta.(F{i}));
end
