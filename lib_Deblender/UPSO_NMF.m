function A = UPSO_NMF(mixed_data,marker_mixed_data,iter,parts,repl,tol,norm_nmf)

X=marker_mixed_data;
sqrteps = sqrt(eps);
dnorm_repl=zeros(repl,1);% gather results from replicates to check in the end which h to choose as final
dnorm_repl_h=cell(repl,2);

for i=1:repl
    S0 = rand(size(marker_mixed_data,1),1)*(prctile(mixed_data(:),99)-min(mixed_data(:))) + min(mixed_data(:));
    A0 = rand(1, size(marker_mixed_data,2));
    for j=1:iter
        upso_iter=4000;
        % % set upso parameters and initializations
        swarm = rand(size(marker_mixed_data,2), parts);
        Xmin=0;
        Xmax=1;
        vmax = 0.01; 
        nmf_objective_function=1;
        % % S fixed, calculate A
        A=nmf_upso(S0, marker_mixed_data,Xmin,Xmax,swarm,upso_iter,vmax,nmf_objective_function);
        A=A';
        
        % A fixed, calculate S
        % % set upso parameters and initializations
        Xmax=prctile(mixed_data(:),99);%  almost maximum value
        Xmid2=prctile(mixed_data(:),80);% Middle bound of particles per coordinate direction
        Xmid1=prctile(mixed_data(:),40);% Middle bound of particles per coordinate direction
        Xmin=min(mixed_data(:));
        
        swarm1 = rand(size(marker_mixed_data,1), parts)*(Xmid1-Xmin) + Xmin;
        swarm2 = rand(size(marker_mixed_data,1), parts)*(Xmid2-Xmid1) + Xmid1;
        swarm3 = rand(size(marker_mixed_data,1), parts)*(Xmax-Xmid2) + Xmid2;
        swarm=[swarm1(:,1:round(0.3*parts)) swarm2(:,1:round(0.4*parts)) swarm3(:,1:(parts-(round(0.3*parts)+round(0.4*parts))))];      
        vmax=Xmax-Xmin/100;
        nmf_objective_function=2;        
        % %                
        S=nmf_upso(A', marker_mixed_data,Xmin, Xmax,swarm,upso_iter,vmax,nmf_objective_function);
            
        % Get norm of difference and max change in factors
        dnorm = sqrt(sum(sum((X-S*A).^2))/(size(X,1)*size(X,2)));
        dS = max(max(abs(S-S0) / (sqrteps+max(max(abs(S0))))));
        dA = max(max(abs(A-A0) / (sqrteps+max(max(abs(A0))))));
        delta = max(dS,dA);

    % Check for convergence
    if j>1
        if delta <= tol
            break;
        elseif dnorm0-dnorm <= tol*max(1,dnorm0)
            break;
        elseif j==iter
            break
        end
    end
    %     initialize for next iteration      
        dnorm0 = dnorm;
        S0 = S;
        A0 = A;
    end
    dnorm_repl(i,1)=dnorm;
    dnorm_repl_h{i,1}= A;
    dnorm_repl_h{i,2}= S;
end

[k1, ~]=find(dnorm_repl==min(dnorm_repl));

A=dnorm_repl_h{k1(1),1};
S=dnorm_repl_h{k1(1),2};

%     normalize nmf A solution based on the l2-norm of S
if norm_nmf==1
    Slen=norm(S,2);
    % Slen=norm(S,1); % use this if you want l1-norm
%     S=S./Slen;
    A=A.*Slen;
end

