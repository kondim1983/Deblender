 function [standard_error, explainedVariances, degrees_freedom,residuals] = calc_standard_error(mixed_data,S_data, A_proportions)
  
% % mixed_data[genes x samples]
% % S_data [genes X cell types], 
% % A_proportions [samples x cell types]

% % Compute 'goodness of fit' metrics 
  residuals = mixed_data-S_data*A_proportions';
  meanSquaredResiduals = sum(residuals.^2,2)./(size(mixed_data,2)-size(S_data,2));
  explainedVariances = 1-(sum(residuals.^2,2)./sum(mixed_data.^2,2));
 
% % Estimate standard errors for each gene
  m = inv(A_proportions'*A_proportions);
  vars = zeros(size(mixed_data,1),size(S_data,2));
  for i=1:size(mixed_data,1)
    vars(i,:) = meanSquaredResiduals(i,:)*diag(m);
  end
  
  standard_error = sqrt(vars);
  degrees_freedom=size(mixed_data,2)-size(S_data,2);