% Reconstruction of zonation expression patterns of enterocytes along the villus axis
% This matlab code loads the Laser-Capture-Microdissection dataset, extracts landmark genes, determines the villus zones
% of all cells in the scRNAseq dataset and produces the zonation table. The code also clusters the data and analyzed enriched GO annotations in each zone.

%% Load the RNAseq data
seq_data=importdata('./raw_data/CSC_dataset_seurat_raw.data_export_171212.txt');
[dd,tSNE_1,tSNE_2,cluster]=textread('./raw_data/CSC_dataset_seurat_tsne_export_171212.txt','%s%f%f%s','headerlines',1);
gene_names_scRNAseq=textread('./raw_data/CSC_dataset_seurat_genes_export_171212.txt','%s');


cid=zeros(length(cluster),1);
cid(find(strcmp(cluster,'enterocyte_1')))=1;
cid(find(strcmp(cluster,'enterocyte_2')))=2;
cid(find(strcmp(cluster,'enterocyte_3')))=3;
cid(find(strcmp(cluster,'enterocyte_4')))=4;
cid(find(strcmp(cluster,'enterocyte_5')))=5;
cid(find(strcmp(cluster,'enterocyte_6')))=6;
cid(find(strcmp(cluster,'transient_amp_1')))=8;
cid(find(strcmp(cluster,'transient_amp_2')))=9;

scRNAseq_mat=seq_data.data;

% Normalize UMI table by the sum of all UMIs per cell
scRNAseq_mat_norm=scRNAseq_mat./repmat(sum(scRNAseq_mat),size(scRNAseq_mat,1),1);

addpath('/Volumes/shalevi/amoor5/Matlab code/export_fig-master'); %export_fig function from https://github.com/altmany/export_fig


%% Handle the LCM data
% Mat_full is the raw TPM for three mice and 5 segments
load('./raw_data/data_LCM.mat')

% Normalize each segment to a sum of 1
mat_norm=mat_full./repmat(sum(mat_full),size(mat_full,1),1);
Ma=mat_norm(:,1:3:end);
Mb=mat_norm(:,2:3:end);
Mc=mat_norm(:,3:3:end);
Mall=zeros(size(Ma,1),size(Ma,2),3);
Mall(:,:,1)=Ma;Mall(:,:,2)=Mb;Mall(:,:,3)=Mc;

% Generate the mean and SEM of expression for each segment across the three mice
mat=mean(Mall,3);
SE=std(Mall,[],3)/sqrt(3);

% Analyze only genes with a maximal mean expression above 10^-4
TPM_THRESHOLD=5*10^-4;
m=max(mat,[],2);
indin=find(m>TPM_THRESHOLD);
mat=mat(indin,:);
SE=SE(indin,:);
mat_full=mat_full(indin,:);
gene_name=gene_name(indin);
m=m(indin);

% calculate profile center of mass and the maximal expression across zones
com=zeros(size(mat,1),1);
for i=1:size(mat,1),
    com(i)=sum((1:5).*mat(i,:)/sum(mat(i,:)));
    [~,mx(i)]=max(mat(i,:));
end

% Display an LCM table with each LCM profile normalized to a maximum of 1
mx=mx(:);
M=mat./repmat(max(mat,[],2),1,size(mat,2));
[y,ord]=sort(com);
figure;
imagesc(M(ord,:));
for i=0:10, line([i+0.5 i+0.5],ylim,'color','k','linewidth',1);end
set(gca,'ytick',0:100:500);
set(gca,'fontsize',15);
xlabel('zones','fontsize',21);
ylabel('genes','fontsize',21);
set(gca,'xtick',1:10);
colorbar;box on;
ylabel('Genes');
xlabel('Zones');
set(gcf,'position',[440 69 330 729]);
export_fig('./figures/supplement/FigS1A_carpet_lcm', '-eps','-painters');


%% selection of signature genes
% choose a subset of  highly expressed and zonated genes to use for scRNAseq reconstruction
THRESH=10^-3;
ind_low_LCM=find(m>THRESH & com<2.5 & mx==1);
ind_high_LCM=find(m>THRESH &  com>3.5 & mx==5);

% Filter out landmark genes that do not appear in the scRNAseq data
genes_low=gene_name(ind_low_LCM);
genes_high=gene_name(ind_high_LCM);
genes_low=intersect(gene_names_scRNAseq,genes_low);
genes_high=intersect(gene_names_scRNAseq,genes_high);

% Find the indices of the bottom and top landmark genes in the LCM dataset
% (after filtering for the genes that did not appear in the scRNAseq
% dataset)
ind_low_LCM=[];
for i=1:length(genes_low), 
    indd=find(strcmpi(gene_name,genes_low{i}));
    ind_low_LCM=[ind_low_LCM indd(1)];
end
ind_high_LCM=[];
for i=1:length(genes_high),
    indd=find(strcmpi(gene_name,genes_high{i})); 
    ind_high_LCM=[ind_high_LCM indd(1)];
end

% Find the indices of the bottom and top landmark genes in the scRNAseq
% gene names
ind_low=[];
for i=1:length(genes_low),
    indd=find(strcmpi(gene_names_scRNAseq,genes_low{i}));
    ind_low=[ind_low indd(1)];
end
ind_high=[];
for i=1:length(genes_high),
    indd=find(strcmpi(gene_names_scRNAseq,genes_high{i}));
    ind_high=[ind_high indd(1)]; 
end

%% Display LCM zonation of low signature genes
figure;
imagesc(M(ind_low_LCM,:));
for i=0:10, line([i+0.5 i+0.5],ylim,'color','k','linewidth',1);end
for i=0:length(ind_low_LCM), line(xlim,[i+0.5 i+0.5],'color','k','linewidth',1);end

set(gca,'ytick',1:length(ind_low_LCM));
set(gca,'yticklabel',gene_name(ind_low_LCM));
set(gca,'fontsize',12);
xlabel('zones','fontsize',21);
ylabel('genes','fontsize',21);
set(gca,'xtick',1:10);
colorbar;box on;
ylabel('Genes');
xlabel('Zones');
set(gcf,'position',[440    46   411   752]);
title('Low signature gene set','fontsize',21);
set(gcf, 'Color', 'w');
export_fig('./figures/main/Fig2B_carpet_lcm_signature_low', '-eps','-painters');



%% Display LCM zonation of high signature genes
figure;
imagesc(M(ind_high_LCM,:));
for i=0:10, line([i+0.5 i+0.5],ylim,'color','k','linewidth',1);end
for i=0:length(ind_high_LCM), line(xlim,[i+0.5 i+0.5],'color','k','linewidth',1);end

set(gca,'ytick',1:length(ind_high_LCM));
set(gca,'yticklabel',gene_name(ind_high_LCM));
set(gca,'fontsize',12);
xlabel('zones','fontsize',21);
ylabel('genes','fontsize',21);
set(gca,'xtick',1:10);
colorbar;box on;
ylabel('Genes');
xlabel('Zones');
set(gcf,'position',[440    46   411   752]);
title('High signature gene set','fontsize',24);
set(gcf, 'Color', 'w');
export_fig('./figures/main/Fig2C_carpet_lcm_signature_high', '-eps','-painters');


%% For each cell calculate the summed normalized expression of the bottom and top landmark genes
vec_low=sum((scRNAseq_mat_norm(ind_low,:)));
vec_high=sum((scRNAseq_mat_norm(ind_high,:)));

% Calculate the zonation variable ("x" in the paper, denoted here as indicator)
% This is the ratio between the summed expression of the top and bottom+top landmark genes 
ind_enterocyte_cells=find(cid<=6); % These are the enterocyte cells
indicator=vec_high(ind_enterocyte_cells)./(vec_low(ind_enterocyte_cells)+vec_high(ind_enterocyte_cells));

%% Color  enterocytes according to the spatial coordinate indicator 
figure;
subplot(1,3,1);
scatter(tSNE_1(ind_enterocyte_cells),tSNE_2(ind_enterocyte_cells),40,vec_low(ind_enterocyte_cells),'filled');
box on;axis square;colorbar;axis tight;
title('Low villus gene sum');
xlabel('tSNE1');ylabel('tSNE2');
set(gca,'fontsize',16);grid on

subplot(1,3,2);
scatter(tSNE_1(ind_enterocyte_cells),tSNE_2(ind_enterocyte_cells),40,vec_high(ind_enterocyte_cells),'filled');
box on;axis square;colorbar;axis tight;
title('High villus gene sum');
xlabel('tSNE1');ylabel('tSNE2');
set(gca,'fontsize',16);grid on

subplot(1,3,3);
scatter(tSNE_1(ind_enterocyte_cells),tSNE_2(ind_enterocyte_cells),40,indicator,'filled');
box on;
axis square;
colorbar;
axis tight;
title('High/(High+Low)');
xlabel('tSNE1');ylabel('tSNE2');
set(gcf,'position',[18 102 1860 851]);
set(gca,'fontsize',16);
grid on
export_fig('./figures/supplement/FigS1HtoJ', '-eps','-painters');


%% Calculate the zonation variable for the bulk LCM populations - this will be used to convert indicator to villus zones 
vec_low_LCM=sum(mat(ind_low_LCM,:));
vec_high_LCM=sum(mat(ind_high_LCM,:));
indicator_LCM=vec_high_LCM./(vec_low_LCM+vec_high_LCM);
figure;
subplot(1,2,1);
bar(1:5,indicator_LCM,'linewidth',2);
xlabel('LCM population');
ylabel('Spatial coordinate (x)');
set(gca,'fontsize',24);

% Examine the range of villus coordinates spanned by the sequenced single enterocytes
subplot(1,2,2);
hist(indicator,30);
hold on;
% overlay the zonation coordinates of the LCM populations
for i=1:length(indicator_LCM),
    line([indicator_LCM(i) indicator_LCM(i)],ylim,'linewidth',2,'color','k','linestyle','--');
end
xlabel('Spatial coordinate');ylabel('Number of cells');
set(gca,'fontsize',24);
axis tight;
xlim([0 1]);

set(gcf,'position',[ 397 296 1153 478]);
export_fig('./figures/supplement/FigSDE', '-eps','-painters');


%% Convert zonation coordinate to one of 6 villus zones using to the LCM zonation coordinates
zone=zeros(length(indicator),1);
boundaries=[0 indicator_LCM 1];
for i=1:length(boundaries)-1
    zone(find(indicator>=boundaries(i) & indicator<boundaries(i+1)))=i;
end

% Present a tSNE plot colored by the inferred villus zone
figure;
scatter(tSNE_1(ind_enterocyte_cells),tSNE_2(ind_enterocyte_cells),100,zone,'filled');
box on;
axis square;
colorbar;
axis tight;
num_zones=max(zone);
yy=linspace(0,1,num_zones);
xlabel('tSNE1');ylabel('tSNE2');
set(gcf,'position',[40 6 1289 792]);
set(gca,'fontsize',16);grid on

% create a zonation table containing mean expression of all genes in each zone (this table is called mn)
% and a table containing for each gene the standard error of the mean for each zone (this table is called se)
mn=zeros(length(gene_names_scRNAseq),max(zone)+1);
se=zeros(length(gene_names_scRNAseq),max(zone)+1);
% Add a first column for the crypt cells, based on the cells annotated as transit amplifying cells
crypt=zeros(length(gene_names_scRNAseq),1);
crypt_cells=find(cid>=8 & cid<=9);
mn(:,1)=mean(scRNAseq_mat_norm(:,crypt_cells),2);
se(:,1)=std(scRNAseq_mat_norm(:,crypt_cells),[],2)/sqrt(length(crypt_cells));

% fill up the zonation table for the 6 villus zones
for j=1:num_zones,
    indd=find(zone==j);    
    mn(:,j+1)=mean(scRNAseq_mat_norm(:,ind_enterocyte_cells(indd)),2);
    se(:,j+1)=std(scRNAseq_mat_norm(:,ind_enterocyte_cells(indd)),[],2)/sqrt(length(ind_enterocyte_cells(indd)));
end
yy=[-0.2 yy];

%% Display the scRNAseq-based zonation table

comSC=zeros(size(mn,1),1);
for i=1:length(comSC)
    comSC(i)=sum((1:size(mn,2)).*mn(i,:)./sum(mn(i,:)));  %calculate the center of mass for each gene (results being between 1-7)
end
indin=find(max(mn,[],2)>1*10^-5); 
[y,ord]=sort(comSC(indin),'descend');
figure;
imagesc(mn(indin(ord),:)./repmat(max(mn(indin(ord),:),[],2),1,size(mn(indin(ord),:),2)));
for i=0:size(mn,2), line([i+0.5 i+0.5],ylim,'color','k','linewidth',1);end
set(gca,'ytick',0:1000:11000);
set(gca,'fontsize',15);
xlabel('zones','fontsize',21);
ylabel('genes','fontsize',21);
set(gca,'xtick',1:size(mn,2));
set(gcf,'position',[440    46   411   752]);
colorbar;box on;
set(gcf, 'Color', 'w');
set(gca,'xticklabel',{'Crypt','V1','V2','V3','V4','V5','V6'});
export_fig('./figures/main/Fig2E_carpet', '-eps','-painters');


%% Assign pvalues for zonation profiles - compare for each profile the dynamic range (max-min)/mean 
% with that obtained in datasets where cell zone is randomly permuted. Zonation significance is generated
% for the 6 villus zones (without the crypt).

val_norm=mean(mn(:,2:end),2);
stat_real=max(mn(:,2:end),[],2)./val_norm-min(mn(:,2:end),[],2)./val_norm;

% Create permutation pvalues
rng(3)  %set the seed for randn
ITER=1000;
stat=zeros(size(mn,1),ITER);
for iter=1:ITER,
    zonesr=zone(randperm(length(zone)));
    mnr=zeros(length(gene_names_scRNAseq),max(zone));
    for j=1:num_zones,
        indd=find(zonesr==j);
        mnr(:,j)=mean(scRNAseq_mat_norm(:,ind_enterocyte_cells(indd)),2);
    end
    val_norm=mean(mnr,2);
    stat(:,iter)=max(mnr,[],2)./val_norm-min(mnr,[],2)./val_norm;
end

Z=(stat_real-mean(stat,2))./std(stat,[],2);
pval=1-normcdf(Z);
pval(pval==0)=1.102*10^(-16);
% Perform multiple hypotheses correction
[~,qval]=mafdr(pval); 

%% Export zonation table
fid=fopen('./raw_data/zonation_table.txt','w');
fprintf(fid,'Gene name\tCrypt_mean\tV1_mean\tV2_mean\tV3_mean\tV4_mean\tV5_mean\tV6_mean\tCrypt_SE\tV1_SE\tV2_SE\tV3_SE\tV4_SE\tV5_SE\tV6_SE\tpval\tqval\n');
for i=1:length(gene_names_scRNAseq)
    fprintf(fid,'%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%e\t%e\n', ...
        gene_names_scRNAseq{i},mn(i,1),mn(i,2),mn(i,3),mn(i,4),mn(i,5),mn(i,6),mn(i,7), ...
        se(i,1),se(i,2),se(i,3),se(i,4),se(i,5),se(i,6),se(i,7), pval(i),qval(i));
end
fclose all;

cell_zone_table(:,1)=dd(ind_enterocyte_cells);
cell_zone_table(:,2)=num2cell(zone);
T = cell2table(cell_zone_table,'VariableNames',{'cell_id','zone'});
writetable(T,'./raw_data/cell_zone_table.txt','Delimiter','\t');
save('./raw_data/gene_zone_matrix.mat','mn','gene_names_scRNAseq','se');


%% Calculate percentage of expressed and signficiantly zonated genes
mx=max(mn,[],2);
indin=find(mx>5*10^-6);
indin2=find(mx>5*10^-6 & qval<0.05);
fprintf('number of genes above threshold is %.f\n',length(indin));
fprintf('number of significantly zonated genes is %.f\n',length(indin2));
fprintf('proportion of significant genes %.2f\n',length(indin2)/length(indin));

%% kmeans clustering
% Use genes with max zonation value > 10^-4
indin=find(max(mn,[],2)>0.0001);
% normalize zonation table to a maximum of 1 for each zonation profile
Mnorm=mn(indin,:)./repmat(max(mn(indin,:),[],2),1,size(mn(indin,:),2));
num_clusters=5;

rng(1);
[idx,C]=kmeans(Mnorm,num_clusters,'Distance','correlation');
%reorder clusters by increasing center of mass of their average zonation
clear p;
clear comp;
for i=1:num_clusters,
    index=find(idx==i); 
    p=mean(Mnorm(index,:));
    p=p/sum(p); 
    comp(i)=sum((1:length(p)).*p);
end
[~,cluster_ord]=sort(comp);

figure;
for i=1:num_clusters
    ind=find(idx==cluster_ord(i));
    subplot(num_clusters,1,i);
    plot(Mnorm(ind,:)','color',[0.7 0.7 0.7]);
    hold on;
    plot(1:size(Mnorm,2),mean(Mnorm(ind,:)),'b','linewidth',5);
    title(['Cluster ' num2str(i)],'fontsize',18);
    set(gca,'xticklabel',{'Crypt','V1','V2','V3','V4','V5','V6'});
    set(gca,'ytick',[]);
    set(gca,'fontsize',12);
end
set(gcf,'position',[1244 6 592 979],'Color', 'w');
fclose all;
export_fig('./figures/main/Fig3_clusters', '-eps','-painters');

