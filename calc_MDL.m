function MDL  = calc_MDL(mixed_genes,mixed_data,full_mixed_data,Kmin, Kmax, low_cutoff, upper_cutoff, coef_var, clust_algo, log_option, algo_solver1,algo_solver2)

% % % Input:
% % % mixed_genes: <column vector> the names of rows in mixed_data
% % % mixed_data: mixed expression matrix class double[genes x samples] containing positive non-log expression values
% % % full_mixed_data: all mixed expression matrix containing positive
% % % non-log expression values, including all recorded probes in the
% % % experiment (will be used to define lower/upper space boundaries)
% % % Kmin, Kmax: integer number of cell types present in the mixture (Kmin>=2)
% % % low_cutoff: filter percentage for genes with the lowest vector norm 
% % % upper_cutoff: filter percentage for genes with the highest vector norm 
% % % coef_var: futher filtering of genes based on coefficient of variation
% % % clust_algo:choose between kmeans-using the centroids (1), kmedoids-using the medoids(2)and kmedoids-using the centroids (3)
% % % log_option: choose to cluster in log (1) or linear space(2)
% % % algo_solver1: choose between the lsqnnoneg solver (1)  or UPSO (2)
% % % algo_solver2: choose between the lsqlin solver (1), quadprog (2) or UPSO (3)

% % % Output:
% % % MDL

addpath(genpath('lib_Deblender'));

if (~isnumeric(mixed_data)) || sum(sum(isnan(mixed_data)))>0 || (min(mixed_data(:))<=0)
       error('The genes_data is not a valid matrix');
end

% % parameters used for calculating proportions
call_NMF=5;% by default we use the case that calls S1&S2 unsupervised mode and applies NMF scheme only in the neighborhood around cluster exemplars
limit_cent_neigh=0.3;% percent of neighbors around cluster exemplars


MDL=zeros(length(Kmin:1:Kmax),1);
index=1;
for i=Kmin:1:Kmax
    % % calculate proportions based on number of sources specified    
    [A_proportions, ~, high_variable_data, ~] = calc_A_unsupervised(mixed_genes, mixed_data,full_mixed_data, i, low_cutoff, upper_cutoff, coef_var, clust_algo, log_option, algo_solver1, call_NMF, limit_cent_neigh);
    
    
    % % calculate the real pure profiles of the highly variable genes based on the estimated proportions
    S_data = calculate_S (A_proportions, full_mixed_data, high_variable_data,algo_solver2);

    % % estimate Minimum Description length (MDL) for the specific number of sources
    % % code for MDL is based on the MDL.R function of the CAM tool of Wang et al. (2016)"Mathematical modelling of transcriptional 
    % % heterogeneity identifies novel markers and subpopulations in complex tissues".DOI: 10.1038/srep18909    
    N1=size(high_variable_data,1);
    N2=size(high_variable_data,2);
    % calculate the error variance
    var_error=var(reshape(high_variable_data'-A_proportions*S_data',N1*N2,1));
    
    % set the number of parameters, calculate likelihood and penalty and then add
    sources=i;
    likelihood=(N1*N2)/2*log(var_error);
    penalty =((sources-1)*N2)/2*log(N1)+((sources*N1)/2*log(N2));
    MDL(index) = likelihood + penalty;
    index=index+1;
end