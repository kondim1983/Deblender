
% % % test unsupervised mode of Deblender tool, i.e. calculating number of tissue types and mixture proportions,
% % % with RNA-Seq dataset downloaded by DECONRNASeq tool (Gong T, Szustakowski JD. DeconRNASeq: a statistical framework for deconvolution 
% % % of heterogeneous tissue samples based on mRNA-Seq data. Bioinformatics. 2013;29: 1083-5.).
% % % The dataset has 10 mixed samples containing 5 different tissues. In the 'summarized_mixed_data' one RefSeq is chosen per gene name in case of multiple RefSeq Ids corresponding to the same gene name.
% % % In few cases RefSeq Ids corresponding to the same gene names displayed the same variance, so the mean profile and one Refseq ID were retained. 
% % % An offset of 0.0001 was added to all values
% % % the output is MDL and the proportion matrix 'A_estimated'
% % % the ground truth proportion matrix is 'A_real'


load('RNA_Seq.mat')

% % test the automatic calculation of the number of cell types without any prior
% % info. For this, calculate MDL criterion in a range of cell types, here we chose {2,8} (the real number is 5).
% % How to set preprocessing cutoffs check the manuscript and the function 'calc_MDL.m'.

Kmin=2;Kmax=8;
low_cutoff=0.05; upper_cutoff=0.05;coef_var=0.4;
algo_solver1=1; algo_solver2=2; clust_algo=1;log_option=1;

MDL  = calc_MDL(summarized_genes,summarized_mixed_data,full_mixed_data,Kmin, Kmax, low_cutoff, upper_cutoff, coef_var, clust_algo, log_option, algo_solver1,algo_solver2);

figure(1)
h=plot((Kmin:1:Kmax),MDL);hold on;set(h, 'Marker', 'o');

% % calculate the proportions for number of cell type equal to 5 based on
% % stageI&stageII
% % How to set these cutoffs check the manuscript.

K_source=5; low_cutoff=0.1; upper_cutoff=0.1; coef_var=0.3; algo_solver=1; clust_algo=1;log_option=1;call_NMF=5;limit_cent_neigh=0.3;

[A_estimated , ~, ~, ~]  = calc_A_unsupervised(summarized_genes, summarized_mixed_data, full_mixed_data, K_source, low_cutoff, upper_cutoff, coef_var, clust_algo, log_option, algo_solver, call_NMF, limit_cent_neigh);

% % try all possible order configurations relative to ground truth in order
% % to find the one with maximum correlation

v=perms((1:1:K_source));
index_corr=zeros(size(v,1),1);
A_real_vectorized=reshape(A_real,size(A_real,1)*size(A_real,2),1);     
for i=1:size(v,1)
        temp_A_est=[A_estimated(:,v(i,1)) A_estimated(:,v(i,2)) A_estimated(:,v(i,3)) A_estimated(:,v(i,4))  A_estimated(:,v(i,5))];
        temp_A_est_vectorized=reshape(temp_A_est,size(temp_A_est,1)*size(temp_A_est,2),1);
        index_corr(i,1)=corr(temp_A_est_vectorized,A_real_vectorized);
end
[p,~]=find(index_corr==max(index_corr));


% % define the best configuration for proportions and plot against ground
% % truth
A_estimated=[A_estimated(:,v(p,1)) A_estimated(:,v(p,2)) A_estimated(:,v(p,3)) A_estimated(:,v(p,4)) A_estimated(:,v(p,5))];

figure(2)
subplot(3,2,1); scatter(A_estimated(:,1),A_real(:,1));hold on;plot((0:0.000001:1),(0:0.000001:1)); title('Brain')
subplot(3,2,2); scatter(A_estimated(:,2),A_real(:,2));hold on;plot((0:0.000001:1),(0:0.000001:1)); title('Muscle')
subplot(3,2,3); scatter(A_estimated(:,3),A_real(:,3));hold on;plot((0:0.000001:1),(0:0.000001:1)); title('Lung')
subplot(3,2,4); scatter(A_estimated(:,4),A_real(:,4));hold on;plot((0:0.000001:1),(0:0.000001:1)); title('Liver')
subplot(3,2,5); scatter(A_estimated(:,5),A_real(:,5));hold on;plot((0:0.000001:1),(0:0.000001:1)); title('Heart')

fprintf('The correlation for Brain is %i.\n',corr(A_estimated(:,1),A_real(:,1)))
fprintf('The correlation for Muscle is %i.\n',corr(A_estimated(:,2),A_real(:,2)))
fprintf('The correlation for Lung is %i.\n',corr(A_estimated(:,3),A_real(:,3)))
fprintf('The correlation for Liver is %i.\n',corr(A_estimated(:,4),A_real(:,4)))
fprintf('The correlation for Heart is %i.\n',corr(A_estimated(:,5),A_real(:,5)))

