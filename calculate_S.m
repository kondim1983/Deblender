function S_data = calculate_S(A_proportions, full_mixed_data, mixed_data,algo_solver)
% % A_proportions: [samples x cell types] mixture proportion matrix
% % full_mixed_data: [genes x samples] gene expression matrix
% % mixed_data: [genes x samples] preprocessed gene expression matrix
% % algo_solver: choose between the lsqlin (1)  or quadprog (2)  or UPSO (3)
% % calculate cell-type-specific gene expression profiles with the use of lsqlin, quadprog or UPSO solver

addpath(genpath('lib_Deblender'));

C=A_proportions;
minimal_val= min(full_mixed_data(:));
maximal_val=max(full_mixed_data(:));

d=mixed_data;
S_data=zeros(size(d,1),size(A_proportions,2));
if algo_solver==1
   % lsqlin Matlab solver
    options = optimset('Algorithm','trust-region-reflective' );
    parfor i=1:size(d,1)
        S_data(i,:)=lsqlin(C,d(i,:),[],[],[],[],repmat(minimal_val,size(C,2),1),repmat(maximal_val,size(C,2),1),[],options);
    end
elseif algo_solver==2
    % quadratic programming Matlab solver
    options = optimset('Algorithm','interior-point-convex');
    parfor i=1:size(d,1)
        S_data(i,:)=quadprog(C'*C,-C'*d(i,:)',[],[],[],[],repmat(minimal_val,size(C,2),1),repmat(maximal_val,size(C,2),1),[],options);  
    end
else
    % UPSO solver
    particles=100;
    Xmax =max(full_mixed_data(:)); % Upper bound of particles per coordinate direction
    Xmid2=prctile(full_mixed_data(:),80);% Middle bound2 of particles per coordinate direction
    Xmid1=prctile(full_mixed_data(:),40);% Middle bound1 of particles per coordinate direction
    Xmin = min(full_mixed_data(:)); % Lower bound of particles per coordinate direction
    vmax=10;
    
    parfor i=1:size(d,1)
        swarm1 = rand(size(C,2), particles)*(Xmid1-Xmin) + Xmin;
        swarm2 = rand(size(C,2), particles)*(Xmid2-Xmid1) + Xmid1;
        swarm3 = rand(size(C,2), particles)*(Xmax-Xmid2) + Xmid2;
        swarm=[swarm1(:,1:round(0.4*particles)) swarm2(:,1:round(0.3*particles)) swarm3(:,1:(particles-(round(0.4*particles)+round(0.3*particles))))];
        S_data(i,:) = opt_UPSO(C,d(i,:)', swarm,Xmin,Xmax,vmax);% call UPSO solver
    end         
end
