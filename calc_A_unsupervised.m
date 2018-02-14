function [A_proportions, high_variable_genes, high_variable_data, cluster_set] = calc_A_unsupervised(mixed_genes, mixed_data, full_mixed_data, K_source, low_cutoff, upper_cutoff, coef_var, clust_algo, log_option, algo_solver, call_NMF, limit_cent_neigh)

% % % Input:
% % % mixed_genes: <column vector> unique gene names, each name corresponding to a row of mixed_data
% % % mixed_data: mixed expression matrix class (gene or probe) double[genes x samples] containing positive non-log expression values 
% % % full_mixed_data: all mixed expression matrix containing positive non-log expression values, including all recorded probes in the experiment
% % % K_source: integer number of cell types present in the mixture (>=2)
% % % low_cutoff: filter proportion of genes with the lowest vector norm 
% % % upper_cutoff: filter proportion of genes with the highest vector norm 
% % % coef_var: filter genes based on coefficient of variation 
% % % clust_algo:choose between kmeans-using the centroids (1), kmedoids-using the medoids(2)and kmedoids-using the centroids (3)
% % % log_option: choose to cluster in log (1) or linear space(2)
% % % algo_solver: choose between the lsqnnoneg solver (1)  or UPSO (2)
% % % call_NMF: apply stage II with adapted NMF scheme, 
% % % (0)not apply 
% % % (1)apply on all mixed data with H prior only, 
% % % (2)apply on all mixed data with H and W prior (calculated with lsqlin),
% % % (3)apply on all mixed data with H and W prior (calculated with quadprog), 
% % % (4)apply on all mixed data with H and W prior (calculated with UPSO)
% % % (5)apply on the data of the neighborhood of the cluster centroid/medoid with H prior only 
% % % limit_cent_neigh: percentage of closest to centroid/medoid genes [0,1],input only if call_NMF >=1

% % % Output:
% % % A_proportions: the proportion matrix [samples x cell types]
% % % high_variable_genes, high_variable_data: the filtered gene dataset
% % % cluster_set: the set of genes per cluster

addpath(genpath('lib_Deblender'));

if nargin<12 && call_NMF~=0
    error('You want to apply both Deblender S1 and S2, insert argument limit_cent_neigh');
end

if ~isnumeric(mixed_data) || ~isequal(length(mixed_data),length(mixed_genes))
       error('The mixed_data is not a valid matrix');
end

if sum(isnan(mixed_data))
       error('The mixed_data have NaNs');
end

if min(mixed_data(:))<=0
       error('The mixed_data contain zero or negative values');
end
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
[p,~]=find(gene_norm>=low_filter & gene_norm<=upper_filter);
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
[h,~]=find(gene_cv>=coef_var);
high_variable_data=filtered_data(h,:);
high_variable_genes=filtered_genes(h);

% % decide to cluster in log or linear space
if log_option==1
    clust_high_variable_data=log2(high_variable_data);
else
    clust_high_variable_data=high_variable_data;
end

% % apply kmeans or kmedoids clustering on the high_variable_data
cluster_set=cell(K_source,1);
if isempty(high_variable_data)==0 || size(high_variable_data,1)>=K_source;
    if clust_algo==1
        try
            [cidx, ~ , ~ , D]=kmeans(clust_high_variable_data,K_source,'Distance','correlation','MaxIter',300000,'Replicates',500,'OnlinePhase','on');
        catch
            warning('Problem using kmeans, it will be rerun'); 
            [cidx, ~ , ~ , D]=kmeans(clust_high_variable_data,K_source,'Distance','correlation','MaxIter',300000,'Replicates',500,'OnlinePhase','on');
        end

        cluster_data=[];
        for m=1:K_source
            tf=ismember(cidx,m);
            cluster_m=mean(high_variable_data(tf==1,:),1);
            cluster_data=[cluster_data;cluster_m];
            cluster_set{m,1}=high_variable_genes(tf==1);
        end
    else
        opts = statset('MaxIter',300000);
        try
            [cidx, cmedoids, ~, D] =kmedoids(clust_high_variable_data,K_source,'Distance','correlation','Algorithm','large','Replicates',500,'OnlinePhase','on','Option',opts);
        catch
            warning('Problem using kmedoids, it will be rerun'); 
            [cidx, cmedoids, ~, D] =kmedoids(clust_high_variable_data,K_source,'Distance','correlation','Algorithm','large','Replicates',500,'OnlinePhase','on','Option',opts);
        end
        
        if log_option==1 && clust_algo==2
            cluster_data=2.^cmedoids;
            for m=1:K_source
                tf=ismember(cidx,m);
                cluster_set{m,1}=high_variable_genes(tf==1);
            end
        elseif log_option ~=1 && clust_algo==2
            cluster_data= cmedoids;
            for m=1:K_source
                tf=ismember(cidx,m);
                cluster_set{m,1}=high_variable_genes(tf==1);
            end
        else
            cluster_data=[];
            for m=1:K_source
                tf=ismember(cidx,m);
                cluster_m=mean(high_variable_data(tf==1,:),1);
                cluster_data=[cluster_data;cluster_m];
                cluster_set{m,1}=high_variable_genes(tf==1);
            end
        end       
        
    end

else
    error('Decrease coef_var!'); 
end

% % calculate proportions based on the desired solver
C=cluster_data';
if algo_solver==1
    d=ones(size(high_variable_data,2),1);
    f=lsqnonneg(C,d);% call Matlab non-negative least squares solver
    A_proportions=C*diag(f);
else
    particles=100;
    swarm=rand(size(C,2),particles)*(0.01-0.0000000001)+0.0000000001;% call UPSO solver and define number of particles(=solutions)
    Xmin=1e-10;
    Xmax=0.01;
    vmax=0.0001;
    f = opt_UPSO(C,ones(size(C,1),1), swarm,Xmin,Xmax,vmax);% call UPSO solver
    A_proportions=C*diag(f);
end

% % rescale A_proporions so that for each sample the sum of proportions
% equals to 1.
for z=1:size(A_proportions,1);
    A_proportions(z,:)=A_proportions(z,:)/sum(A_proportions(z,:));
end

% % call NMF as stage II to estimate the proportion matrix
if call_NMF>0
    
   subset_cluster_set=cell(K_source,1);
   for m=1:K_source
                % %  find genes per cluster
       tf=ismember(cidx,m);
       temp_clust_genes=high_variable_genes(tf==1);
                % %  find distance of each gene in the cluster to centroid
       temp_D=D(tf==1,m);
                % %  sort distances and find the top 'limit_cent_neigh'
                % % percent of genes closest to centroids
       sort_temp_D=sort(temp_D,'ascend');
       max_cent_dist=sort_temp_D(ceil(limit_cent_neigh*length(temp_clust_genes)));
       [p,~]=find(temp_D<=max_cent_dist);
       subset_cluster_set{m,1}=temp_clust_genes(p);
    end

    if call_NMF==1 

        marker_cell_type_index=zeros(size(mixed_genes,1),1);
        for qf=1:K_source
            tf=ismember(mixed_genes,subset_cluster_set{qf,:});
            marker_cell_type_index(tf==1,1)=qf;
        end


        tol=1e-4;%NMF tolerance
        num_iter=100000; 
        repl=100; % run NMF multiple replicates
        norm_nmf=1; % % if set to (0) the nmf solution is not normalized, if set to (1) is normalized by normalizing 
        % % the columns of W to h have L2-norm equal to 1 and then H is
        % % updated accordingly and if set to (2) the rows of H are normalized to have
        % % L2-norm equal to 1 and then W is updated accordingly (as in 'nnmf' Matlab function)
        H = adapted_nmf1(mixed_data,K_source,marker_cell_type_index,tol,num_iter,repl,norm_nmf,A_proportions',[]);
        A_proportions=H';
        
    elseif call_NMF>=2 && call_NMF<=4
        marker_cell_type_index=zeros(size(mixed_genes,1),1);
        for qf=1:K_source
            tf=ismember(mixed_genes,subset_cluster_set{qf,:});
            marker_cell_type_index(tf==1,1)=qf;
        end

        S_data = calculate_S (A_proportions, full_mixed_data, mixed_data,call_NMF-1);

        tol=1e-4;%NMF tolerance
        num_iter=100000;
        repl=100; % run NMF multiple replicates
        norm_nmf=1; % if set to 1 the nmf solution is normalized
        H =adapted_nmf1(mixed_data,K_source,marker_cell_type_index,tol,num_iter,repl,norm_nmf,A_proportions',S_data);
        A_proportions=H';

    elseif call_NMF>4

        subset_cluster_set_data=[];    
        marker_cell_type_index=[];
        for qf=1:K_source
            tf=ismember(mixed_genes,subset_cluster_set{qf,:});
            subset_cluster_set_data_temp=mixed_data(tf==1,:);
            subset_cluster_set_data=[subset_cluster_set_data;subset_cluster_set_data_temp];
            marker_cell_type_index_temp=qf*ones(size(subset_cluster_set_data_temp,1),1);
            marker_cell_type_index=[marker_cell_type_index;marker_cell_type_index_temp];
        end

        tol=1e-4;%NMF tolerance
        num_iter=100000; 
        repl=100; % run NMF multiple replicates
        norm_nmf=1; % if set to (0) the nmf solution is not normalized, if set to (1) is normalized by normalizing 
        % % the columns of W to have L2-norm equal to 1 and then H is updated accordingly and if set to (2) the rows of H are normalized to have
        % % L2-norm equal to 1 and then W is updated accordingly (as in 'nnmf' Matlab function)
        H =adapted_nmf1(subset_cluster_set_data,K_source,marker_cell_type_index,tol,num_iter,repl,norm_nmf,A_proportions',[]);
        A_proportions=H';

    end
end
