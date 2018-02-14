function A_proportions = calc_A_known_marker_genes (mixed_genes, mixed_data,full_mixed_data, marker_genes, marker_cell_type_index, algo_solver1, call_NMF)
% % % Input:
% % % mixed_genes:<column vector> all gene names (class: cell or double) of mixed dataset with 
% % % each each name corresponding to the row of "mixed_data"
% % % mixed_data: expression data matrix (class double) of mixed samples(non-logged)
% % % marker_genes: <column vector> genes (class: cell or double) with same nomenclature like "mixed_genes" 
% % % marker_cell_type_index: column vector (same length as marker_genes) with integer values 1...N where
% % % N is the number of cell types, which designates the cell type where each marker gene belongs to
% % % algo_solver1: choose between lsqnonneg (1) or UPSO (2)
% % % call_NMF: apply stage II with adapted NMF scheme, (0) not apply 
% % % (1)apply with H prior only, (2)apply with H and W prior (calculated with
% % % lsqlin),(3)apply with H and W prior (calculated with quadprog), (4)apply with H and W prior (calculated with UPSO)

% % % Output:
% % % A_proportions:proportion matrix [samples x sources]


% % ########################################################################
addpath(genpath('lib_Deblender'));

if ~isequal(class(mixed_genes), class(marker_genes))
       error('The nomenclature of genes and marker genes is not the same');
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

% % construct marker meta-profile data based on marker gene lists specific
% % for each cell type
marker_data=create_marker_data(mixed_genes, mixed_data, marker_genes, marker_cell_type_index);

% % calculate proportions with the use of lsqnonneg or UPSO solver
if isempty(marker_data)==0 && size(marker_data,1)==max(marker_cell_type_index) 
    if algo_solver1==1
        d=ones(size(marker_data,2),1);
        fS=lsqnonneg(marker_data',d);
        A_proportions=marker_data'*diag(fS);
    else
        particles=100;
        swarm=rand(size(marker_data',2),particles)*(0.01-0.0000000001)+0.0000000001;% call UPSO solver and define number of particles(=solutions)
        Xmin=1e-10;
        Xmax=0.01;
        vmax=0.0001;
        fS = opt_UPSO(marker_data',ones(size(marker_data',1),1), swarm,Xmin,Xmax,vmax);% call UPSO solver
        A_proportions=marker_data'*diag(fS);
    end
else
    error('Error. The marker data are not valid, check input')
end

% % scale A_proportions so that sum equals to 1
for j=1:size(A_proportions,1)
    A_proportions(j,:)=A_proportions(j,:)/sum(A_proportions(j,:));
end

if call_NMF==1

    % % % run NMF
    nmf_marker_cell_type_index=zeros(size(mixed_data,1),1);
    for qf=1:max(marker_cell_type_index)
        tf=ismember(marker_cell_type_index,qf);
        temp=marker_genes(tf==1);
        tf=ismember(mixed_genes,temp);
        nmf_marker_cell_type_index(tf==1,1)=qf;
    end
    
    tol=1e-4;%NMF tolerance
    maxiter=100000; 
    repl=500; % run NMF multiple replicates
    norm_nmf=1; % % if set to (0) the nmf solution is not normalized, if set to (1) is normalized by normalizing 
    % % the columns of W to have L2-norm equal to 1 and then H is updated accordingly so that objective function 
    % % does not change  and if set to (2) the rows of H are normalized to have
    % % L2-norm equal to 1 and then W is updated accordingly (as in 'nnmf' Matlab function)
    h =adapted_nmf1(mixed_data,max(marker_cell_type_index),nmf_marker_cell_type_index,tol,maxiter,repl,norm_nmf,A_proportions',[]);
    A_proportions=h';

elseif call_NMF>=2 && call_NMF<=4
    S_data = calculate_S (A_proportions, full_mixed_data, mixed_data,call_NMF-1);
    % % % run NMF
    nmf_marker_cell_type_index=zeros(size(mixed_data,1),1);
    for qf=1:max(marker_cell_type_index)
        tf=ismember(marker_cell_type_index,qf);
        temp=marker_genes(tf==1);
        tf=ismember(mixed_genes,temp);
        nmf_marker_cell_type_index(tf==1,1)=qf;
    end

    tol=1e-4;%NMF tolerance
    maxiter=100000; 
    repl=500; % run NMF multiple replicates
    norm_nmf=1; % if set to 1 the nmf solution is normalized
    h = adapted_nmf1(mixed_data,max(marker_cell_type_index),nmf_marker_cell_type_index,tol,maxiter,repl,norm_nmf,A_proportions',S_data);
    A_proportions=h';
elseif call_NMF>4

        subset_cluster_set_data=[];    
        nmf_marker_cell_type_index=[];
        for qf=1:max(marker_cell_type_index)
            tf=ismember(marker_cell_type_index,qf);
            temp=marker_genes(tf==1);
            tf=ismember(mixed_genes,temp);
            subset_cluster_set_data_temp=mixed_data(tf==1,:);
            subset_cluster_set_data=[subset_cluster_set_data;subset_cluster_set_data_temp];
            nmf_marker_cell_type_index_temp=qf*ones(size(subset_cluster_set_data_temp,1),1);
            nmf_marker_cell_type_index=[nmf_marker_cell_type_index;nmf_marker_cell_type_index_temp];
        end

        tol=1e-4;%NMF tolerance
        num_iter=100000; 
        repl=500; % run NMF multiple replicates
        norm_nmf=1;
        h =adapted_nmf1(subset_cluster_set_data,max(marker_cell_type_index),nmf_marker_cell_type_index,tol,num_iter,repl,norm_nmf,A_proportions',[]);
        A_proportions=h';
end



