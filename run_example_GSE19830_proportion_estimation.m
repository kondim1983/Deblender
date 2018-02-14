% % % test unsupervised mode of Deblender tool with GSE19830 dataset. It has 33 mixed
% % % samples including 3 tissues. In this example, all probes of the dataset are included except the control probes starting with AFFX-. 
% % % the input is the mixed expression data 'mixed_data', the probe names 'mixed_genes'. 
% % % the output is the proportion matrix 'A_estimated'
% % % the ground truth proportion matrix is 'A_real'


load('GSE19830.mat')

% % run mixture proportion estimation with Stage I
% % How to set these cutoffs check the manuscript.

K_source=3; low_cutoff=0.1; upper_cutoff=0.1; coef_var=0.1;clust_algo=1; log_option=1;algo_solver=1;call_NMF=0;limit_cent_neigh=0;

[A_estimated , ~, ~, ~]  = calc_A_unsupervised(mixed_genes, mixed_data, mixed_data, K_source, low_cutoff, upper_cutoff, coef_var, clust_algo, log_option, algo_solver, call_NMF, limit_cent_neigh);


% % try all possible order configurations relative to ground truth in order
% % to find the one with maximum correlation
v=perms((1:1:K_source));
index_corr=zeros(size(v,1),1);
A_real_vectorized=reshape(A_real,size(A_real,1)*size(A_real,2),1);     
for i=1:size(v,1)
        temp_A_est=[A_estimated(:,v(i,1)) A_estimated(:,v(i,2)) A_estimated(:,v(i,3))];
        temp_A_est_vectorized=reshape(temp_A_est,size(temp_A_est,1)*size(temp_A_est,2),1);
        index_corr(i,1)=corr(temp_A_est_vectorized,A_real_vectorized);
end
[p, ~]=find(index_corr==max(index_corr));


% % define the best configuration for proportions and plot against ground
% truth
A_estimated=[A_estimated(:,v(p,1)) A_estimated(:,v(p,2)) A_estimated(:,v(p,3))];

figure('Name', 'Mixture proportions')
subplot(2,2,1); scatter(A_estimated(:,1),A_real(:,1));hold on;plot((0:0.000001:1),(0:0.000001:1)); title('Liver')
subplot(2,2,2); scatter(A_estimated(:,2),A_real(:,2));hold on;plot((0:0.000001:1),(0:0.000001:1)); title('Brain')
subplot(2,2,3); scatter(A_estimated(:,3),A_real(:,3));hold on;plot((0:0.000001:1),(0:0.000001:1)); title('Lung')

fprintf('The correlation for Liver is %i.\n',corr(A_estimated(:,1),A_real(:,1)))
fprintf('The correlation for Brain is %i.\n',corr(A_estimated(:,2),A_real(:,2)))
fprintf('The correlation for Lung is %i.\n',corr(A_estimated(:,3),A_real(:,3)))