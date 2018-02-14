function A_under = calc_A_known_marker_underdetermined(mixed_genes, mixed_data,marker_genes,marker_cell_type_index, NMF_option)

% % % Input:
% % % mixed_genes:<column vector> all gene names (class: cell or double) of mixed dataset with 
% % % each each name corresponding to the row of "mixed_data"
% % % mixed_data: expression data matrix (class double) of mixed samples (non-logged)
% % % marker_genes: <column vector> genes (class: cell or double) with same nomenclature like "genes", all marker genes are included in the mixed_genes 
% % % marker_cell_type_index: column vector (same length as marker_genes) with integer values 1...N where
% % % N is the number of cell types, which designates the source where each marker gene belongs to
% % % NMF_option: three adapted NMF schemes,(1)run multiplicative formula either with Matlab 'nnmf'
% % % normalization,(2)with other normalization or (3) apply UPSO-NMF

% % % Output:
% % % A_under:proportion matrix [samples x sources]
% ##############################################################
% % % Compute proportions with non-negative matrix factorization by solving
% % %  each cell-type marker set independently

addpath(genpath('lib_Deblender'));

if ~isequal(class(mixed_genes), class(marker_genes))
       error('The nomenclature of genes and marker genes is not the same');
end

if ~isnumeric(mixed_data);
       error('The mixed_data are not a valid matrix');
end

if sum(sum(isnan(mixed_data)))>0
       error('The mixed_data have NaNs');
end

if min(mixed_data(:))<=0
       error('The mixed_data contain zero or negative values');
end

% %  Calculate cell type-specific proportion vectors and then combine
A=[];
for i=1:max(marker_cell_type_index)
    tf=ismember(marker_cell_type_index,i);
    temp_markers=marker_genes(tf==1);
    tf=ismember(mixed_genes,temp_markers);
    temp_marker_mixed_data=mixed_data(tf==1,:);
    if size(temp_marker_mixed_data,1)*size(temp_marker_mixed_data,2)>= 2*size(temp_marker_mixed_data,1)
        if NMF_option==1
            opt = statset('MaxIter',100000,'TolFun',1e-4,'TolX',1e-4);
            [~,A_m] = nnmf(temp_marker_mixed_data,1,'replicates',500,'algorithm','mult','options',opt);
            A=[A;A_m];
        elseif NMF_option==2         
            repl=500;
            tol=1e-4;%NMF tolerance
            num_iter=10000;
            norm_nmf=1;
            A_m=adapted_nmf2(temp_marker_mixed_data,num_iter,tol,repl,norm_nmf);  
            A=[A;A_m];
        else
            iter=10000;
            parts=300;
            repl=500;
            tol=1e-4;
            norm_nmf=1;
            A_m = UPSO_NMF(mixed_data,temp_marker_mixed_data,iter,parts,repl,tol,norm_nmf);
            A=[A;A_m];
        end
    else
        error('You need more samples!');
    end
 
end
% scale the proportions so that for each sample the sum equal to 1
A_under=A';
for k=1:size(A_under,1);
    A_under(k,:)=A_under(k,:)/sum(A_under(k,:));
end


