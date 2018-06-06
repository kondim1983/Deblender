function H = adapted_nmf1(mixed_data,k,marker_cell_type_index,tol,maxiter,repl,norm_nmf,h0_prior,w0_prior)
% % the function is based on the nnmf function from mathworks (https://se.mathworks.com/help/stats/nnmf.html)

a=mixed_data;
[n,m]=size(a);
sqrteps = sqrt(eps);

dnorm_repl=zeros(repl,1);% gather results from replicates to find the best H
dnorm_repl_h=cell(repl,2);
% % % MULTIPLICATIVE UPDATE% % % % % 
for i=1:repl
    if i==1 && isempty(w0_prior)==1
        h0=h0_prior;
        w0=rand(n,k);

        for qz=1:k
            tf=ismember(marker_cell_type_index,qz);
            w0(tf==1,1:end ~= qz)=0;
        end
    elseif i==1 && isempty(w0_prior)==0
        h0=h0_prior;
        w0=w0_prior;
        
        for qz=1:k
            tf=ismember(marker_cell_type_index,qz);
            w0(tf==1,1:end ~= qz)=0;
        end
    else
        h0=rand(k,m);
        w0=rand(n,k);

        for qz=1:k
            tf=ismember(marker_cell_type_index,qz);
            w0(tf==1,1:end ~= qz)=0;
        end
    end
    
    for j=1:maxiter
        % Multiplicative update 
    %     calculate h
        numer = w0'*a;
        h = max(0,h0 .* (numer ./ ((w0'*w0)*h0 + eps(numer))));
    %     calculate w
        numer = a*h';
        w = max(0,w0 .* (numer ./ (w0*(h*h') + eps(numer))));

    % Get norm of difference and max change in factors
        dnorm = sqrt(sum(sum((a-w*h).^2))/(n*m));
        dw = max(max(abs(w-w0) / (sqrteps+max(max(abs(w0))))));
        dh = max(max(abs(h-h0) / (sqrteps+max(max(abs(h0))))));
        delta = max(dw,dh);

    % Check for convergence
        if j>1
            if delta <= tol
                break;
            elseif dnorm0-dnorm <= tol*max(1,dnorm0)
                break;
            elseif j==maxiter
                break
            end
        end

    %     initialize for next iteration      
        dnorm0 = dnorm;
    %     ensure that marker genes are equal to zero in the
    %     non-specific cell types    
        for qw=1:k
            tf=ismember(marker_cell_type_index,qw);
            w(tf==1,1:end ~= qw)=0;
        end
        w0 = w;  
        h0 = h;   
    end

    %     return best result 
    dnorm_repl(i,1)=dnorm;
    dnorm_repl_h{i,1}=h;
    dnorm_repl_h{i,2}=w;
end

% % find the best result across replicates
[k1, ~]=find(dnorm_repl==min(dnorm_repl));
H=dnorm_repl_h{(k1(1)),1};
W=dnorm_repl_h{(k1(1)),2};

% % % % % % % % 
% % if norm_nmf is set to (0) the nmf solution is not normalized, if set to (1) is normalized by normalizing 
% % the columns of W to have L2-norm equal to 1 and then H is updated accordingly so that objective function 
% % does not change and if set to (2) the rows of H are normalized to have L2-norm equal to 1 and 
% % then W is updated accordingly (as in 'nnmf' Matlab function)
if norm_nmf==1
    
    wlen = sqrt(sum(W.^2,1));
%     wlen = sum(W,1); % use this if you want l1-norm
%     wbest = bsxfun(@times,W,1./wlen);% unused
    H = bsxfun(@times,H,wlen');
    for g=1:size(H,2)
        H(:,g)=H(:,g)./sum(H(:,g));
    end
elseif norm_nmf==2
    hlen=sqrt(sum(H.^2,2));
    H = bsxfun(@times,H,1./hlen);

    for g=1:size(H,2)
        H(:,g)=H(:,g)./sum(H(:,g));
    end
    
elseif norm_nmf==0
    for g=1:size(H,2)
        H(:,g)=H(:,g)./sum(H(:,g));
    end
end
