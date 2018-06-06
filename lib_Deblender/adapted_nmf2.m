function H =adapted_nmf2(marker_mixed_data,num_iter,tol,repl,norm_nmf)
% % the function is based on the nnmf function from mathworks (https://se.mathworks.com/help/stats/nnmf.html)

a=marker_mixed_data;
[n,m]=size(a);
sqrteps = sqrt(eps);

dnorm_repl=zeros(repl,1);% gather results from replicates to find best H
dnorm_repl_h=cell(repl,2);
for i=1:repl
    w0=rand(n,1);
    h0=rand(1,m);

    for j=1:num_iter
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
        elseif j==num_iter
            break
        end
    end
    %     initialize for next iteration      
        dnorm0 = dnorm;
        w0 = w;  
        h0 = h;   
    end
    %     find best result and return it
    dnorm_repl(i,1)=dnorm;
    dnorm_repl_h{i,1}=h;
    dnorm_repl_h{i,2}=w;
end

%     find best result from all replicates and return it

[k1, ~]=find(dnorm_repl==min(dnorm_repl));
H = dnorm_repl_h{(k1(1)),1};
W = dnorm_repl_h{(k1(1)),2};

%     normalize nmf H solution based on the l2-norm of W
if norm_nmf==1
    wlen=norm(W,2);
    % wlen = norm(W,1); % use this if you want l1-norm
    H=H.*wlen;
end




