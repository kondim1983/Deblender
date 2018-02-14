function [A_under,high_variable_genes,high_variable_data, cluster_set]= calc_A_unsupervised_underdetermined(mixed_genes, mixed_data, low_cutoff, upper_cutoff, coef_var, K_source, limit_cent_neigh, clust_algo, log_option, NMF_option)
% % % Input:
% % % mixed_genes: <column vector> unique gene names, each name corresponding to a row of mixed_data
% % % mixed_data: expression data matrix (class double) of mixed samples (non-logged)[genes x samples]
% % % low_cutoff: filter proportion of genes with the lowest vector norm 
% % % upper_cutoff: filter proportion of genes genes with the highest vector norm 
% % % coef_var: filter genes based on coefficient of variation
% % % K_source: integer number of cell types present in the mixture (>=2)
% % % limit_cent_neigh: percentage of closest profiles to centroids
% % % clust_algo:choose between kmeans (1) and kmedoids (2)
% % % log_option: choose to cluster in log (1) or linear space(2)
% % % NMF_option: three adapted NMF schemes,(1)run Matlab 'nnmf' function that does normalization 
% % % on the rows of H first and then on cols of W,(2)with normalization
% % % on cols of W first and then on rows of H (with option to disable normalization) , 
% % % or(3) apply UPSO-NMF adapted scheme with space bounds and with normalization
% % % on cols of W first and then on rows of H (with option to disable normalization)
 
% % % Output:
% % % A_under:proportion matrix [samples x cell types]
% % % high_variable_genes,high_variable_data: the filtered dataset
% % % cluster_set: the gene set assigned to each cell type after choosing
% % % the closest to centroid genes (as defined by the limit_cent_neigh cutoff )

% % ######################################################################
addpath(genpath('lib_Deblender'));

if ~isnumeric(mixed_data);
       error('The mixed_data are not a valid matrix');
end

if sum(sum(isnan(mixed_data)))>0
       error('The mixed_data have NaNs');
end

if min(mixed_data(:))<=0
       error('The mixed_data contain zero or negative values');
end

% % #################PREPROCESSING########################################
% % calculate the frobenius norm of each gene expression vector
rows=size(mixed_data,1);
gene_norm=zeros(rows,1);
for i=1:rows
    gene_norm(i)=norm(mixed_data(i,:),'fro');
end

% % define the limits to filter the mixed_data 
gene_norm_sorted=sort(gene_norm,'ascend');
low_filter=gene_norm_sorted((round(numel(gene_norm_sorted)*low_cutoff)+1));
upper_filter=gene_norm_sorted(rows-(round(numel(gene_norm_sorted)*upper_cutoff)));

% % keep data that passed the high/low expression filter
[p, ~]=find(gene_norm>=low_filter & gene_norm<=upper_filter);
filtered_data=mixed_data(p,:);
filtered_genes=mixed_genes(p);

% % calculate coefficient of variation for the filtered_data
if isempty(filtered_data)==0;
    rows1=size(filtered_data,1);
    gene_cv=zeros(rows1,1);
    for j=1:rows1
        gene_cv(j)=std(filtered_data(j,:))/mean(filtered_data(j,:));
    end
else
   error('Decrease low_cutoff and/or upper_cutoff!'); 
end

% % filter further the filtered_data based on coefficient of variation
[h, ~]=find(gene_cv>=coef_var);
high_variable_data=filtered_data(h,:);
high_variable_genes=filtered_genes(h);

% % decide to cluster in log or linear space
if log_option==1
    clust_high_variable_data=log2(high_variable_data);
else
    clust_high_variable_data=high_variable_data;
end

% % ######################################################################
% % #################CLUSTERING###########################################

if isempty(high_variable_data)==0 || size(high_variable_data,1)>=K_source;
    if clust_algo==1
        try
            [cidx,~,~,D]=kmeans(clust_high_variable_data,K_source,'Distance','correlation','MaxIter',300000,'Replicates',500,'OnlinePhase','on');
        catch
            warning('Problem using kmeans, it will be rerun'); 
            [cidx,~,~,D]=kmeans(clust_high_variable_data,K_source,'Distance','correlation','MaxIter',300000,'Replicates',500,'OnlinePhase','on');
        end
    else
        opts = statset('MaxIter',300000);
        try
            [cidx,~,~,D] =kmedoids(clust_high_variable_data,K_source,'Distance','correlation','Algorithm','large','Replicates',500,'OnlinePhase','on','Option',opts);
        catch
            warning('Problem using kmedoids, it will be rerun'); 
            [cidx,~,~,D] =kmedoids(clust_high_variable_data,K_source,'Distance','correlation','Algorithm','large','Replicates',500,'OnlinePhase','on','Option',opts);
        end
        
    end

else
    error('Decrease coef_var!'); 
end
% % ######################################################################
% % #################CALCULATING_A WITH NMF################################

if max(cidx)==K_source
    cluster_set=cell(K_source,1);
    A=[];
    for m=1:K_source
        % %  find genes per cluster
        tf=ismember(cidx,m);
        temp_clust_data=high_variable_data(tf==1,:);
        temp_clust_genes=high_variable_genes(tf==1);
        % %  find distance of each gene in the cluster to centroid
        temp_D=D(tf==1,m);
        % %  sort distances and find the top 'limit_cent_neigh' closest to centroid genes
        sort_temp_D=sort(temp_D,'ascend');
        max_cent_dist=sort_temp_D(ceil(limit_cent_neigh*length(temp_clust_genes)));
        [p, ~]=find(temp_D<=max_cent_dist);
        marker_mixed_data=temp_clust_data(p,:);
        cluster_set{m,1}=temp_clust_genes(p);
        if NMF_option==1
        % % %     apply nnmf matlab function with factorization rank equal to 1
            opt = statset('MaxIter',100000,'TolFun',1e-4,'TolX',1e-4 );   
            [~, A_m] = nnmf(marker_mixed_data,1,'replicates',100,'algorithm','mult','options',opt);
            A=[A;A_m];
        elseif NMF_option==2
            repl=100;
            tol=1e-4;%NMF tolerance
            num_iter=100000;
            norm_nmf=1;
            A_m = adapted_nmf2(marker_mixed_data,num_iter,tol,repl,norm_nmf);  
            A=[A;A_m];
        else
            iter=10000;
            parts=300;
            repl=100;
            tol=1e-4;
            norm_nmf=1;
            A_m = UPSO_NMF(mixed_data,marker_mixed_data,iter,parts,repl,tol,norm_nmf);
            A=[A; A_m];
        end
    end
else
    error('The data were not properly clustered to proceed, try different preprocessing!');
end
% % #########################################################################

% scale the proportions so that for each sample the sum is equal to 1
A_under=A';
for qf=1:size(A_under,1);
    A_under(qf,:)=A_under(qf,:)/sum(A_under(qf,:));
end