function xopt  = nmf_upso(nmf_input, marker_mixed_data,minimal_val, maximal_val,swarm_initial,upso_iter,max_velocity,nmf_objective_function)
%% IMPLEMENTATION OF UNIFIED PSO (UPSO) %
%===================
% PARAMETER SETTING
%===================
% PROBLEM-RELATED PARAMETERS
if nmf_objective_function==1
    dim = size(marker_mixed_data,2); % Problem dimension
else
    dim = size(marker_mixed_data,1);
end

swarm=swarm_initial;
Xmax=maximal_val;
Xmin=minimal_val;
vmax = max_velocity; % Maximum velocity   
sqrteps=sqrt(eps);

%  PSO-RELATED PARAMETERS
u = 0.5;
tol = 1e-6; % tolerance
SS = size(swarm_initial,2); % Swarm size
MaxIt = upso_iter;% Maximum number of iterations
chi = 0.729; % Constriction coefficient
c1 = 2.05; % Cognitive parameter
c2 = 2.05; % Social parameter
NR =1; % Neighborhood radius
 
% SWARM, VELOCITIES & BEST POSITION INITIALIZATION
vel = rand(dim,SS)*0.0000000001;
bestpos = swarm;
fmin_counter=zeros(upso_iter,1);
xopt_counter=cell(upso_iter,1);
iter = 1; % Iteration counter
STOP = 0; % Experiment stopping flag
    % EVALUATE SWARM
if  nmf_objective_function==1
    for i=1:SS,
        fswarm(i)=sqrt(sum(sum((marker_mixed_data-nmf_input*swarm(:,i)').^2))/(size(marker_mixed_data,1)*size(marker_mixed_data,2)));
        fbestpos(i) = fswarm(i);
    end
else
    for i=1:SS,
        fswarm(i)=sqrt(sum(sum((marker_mixed_data-swarm(:,i)*nmf_input').^2))/(size(marker_mixed_data,1)*size(marker_mixed_data,2)));
        fbestpos(i) = fswarm(i);
    end
end
fmin_counter(1,1)=min(fbestpos);
    % UPDATE OVERALL BEST
[fxopt, g_over] = min(fbestpos);% objective function value and index of best particle
xopt = bestpos(:,g_over);% best particle
xopt_counter{1,1}=xopt;
    % UPDATE LOCAL BEST PARTICLE INDICES PER NEIGHBORHOOD
g_neig = 1:SS;
for i=1:SS,
    for j=i-NR:i+NR,
        if (j<=0)
           jc = SS+j;
        elseif (j>SS)
           jc = j-SS;
        else
           jc = j;
        end
        if (fbestpos(jc)<fbestpos(g_neig(i)))
           g_neig(i) = jc;
        end
    end
end
    %======================
    % SWARM EVOLUTION LOOP
    %======================
while (STOP == 0)
        % UPDATE ITERATION COUNTER
      iter = iter+1;
        % UPDATE VELOCITIES
      for i=1:SS,
            % Global PSO component
          R1 = rand(dim,1);
          R2 = rand(dim,1);
          G = chi .* (vel(:,i) + c1.*R1.*(bestpos(:,i)-swarm(:,i)) + c2.*R2.*(bestpos(:,g_over)-swarm(:,i)));
            % Local PSO component
          R1 = rand(dim,1);
          R2 = rand(dim,1);
          L = chi .* (vel(:,i) + c1.*R1.*(bestpos(:,i)-swarm(:,i)) + c2.*R2.*(bestpos(:,g_neig(i))-swarm(:,i)));
            % Standard UPSO velocity update
          vel(:,i) = u.*G + (1-u).*L;
      end
        % CONSTRAIN VELOCITIES
      vel(vel>vmax) = vmax;
      vel(vel<-vmax) = -vmax;
        % UPDATE SWARM
      swarm = swarm + vel;
        % CONSTRAIN SWARM
      swarm(swarm>Xmax) = Xmax;
      swarm(swarm<Xmin) = Xmin;
        % EVALUATE SWARM
      if  nmf_objective_function==1
          for i=1:SS,
              fswarm(i)=sqrt(sum(sum((marker_mixed_data-nmf_input*swarm(:,i)').^2))/(size(marker_mixed_data,1)*size(marker_mixed_data,2)));
          end
      else
          for i=1:SS,
              fswarm(i)=sqrt(sum(sum((marker_mixed_data-swarm(:,i)*nmf_input').^2))/(size(marker_mixed_data,1)*size(marker_mixed_data,2)));
          end
      end
        % UPDATE BEST POSITIONS
      bestpos(:,fswarm<fbestpos) = swarm(:,fswarm<fbestpos);
      fbestpos(fswarm<fbestpos) = fswarm(fswarm<fbestpos);
      fmin_counter(iter)=min(fbestpos);
        % UPDATE OVERALL BEST
      [fxopt, g_over] = min(fbestpos);
      xopt = bestpos(:,g_over);
      xopt_counter{iter,1}=xopt;
        % UPDATE LOCAL BEST PARTICLE INDICES PER NEIGHBORHOOD
      g_neig = 1:SS;
      for i=1:SS,
          for j=i-NR:i+NR,
              if (j<=0)
                 jc = SS+j;
              elseif (j>SS)
                 jc = j-SS;
              else
                 jc = j;
              end
              if (fbestpos(jc)<fbestpos(g_neig(i)))
                  g_neig(i) = jc;
              end
          end
      end
        % CHECK STOPPING CRITERION
      if iter>=MaxIt
         STOP = 1;
         % %     check every 300 iterations if the objective function and solution converged 
      elseif mod(iter,300)==1 && ((fmin_counter(iter-10,1)-min(fbestpos))<=tol)&&(max(max(abs(xopt-xopt_counter{iter-10,:}) / (sqrteps+max(max(abs(xopt_counter{iter-10,:})))))))<=tol
         STOP = 1; 
      end   
        
end % FINISH SWARM EVOLUTION
%===========================
end





