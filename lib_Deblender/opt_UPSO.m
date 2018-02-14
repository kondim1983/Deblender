function fS = opt_UPSO(C,d,swarm_initialization,minimal_val,maximal_val,velocity_max)
% % IMPLEMENTATION OF UNIFIED PSO (UPSO) %%
%===================
% PARAMETER SETTING
%===================
% PROBLEM-RELATED PARAMETERS
swarm=swarm_initialization;% initialized swarm (each column represents a particle solution)
Xmax=maximal_val;% Upper bound of particles per coordinate direction
Xmin = minimal_val; % Lower bound of particles per coordinate direction
vmax=velocity_max;

%  PSO-RELATED PARAMETERS
u=0.5; % unification factor
tol = 1e-12; % Desired accuracy
dim = size(swarm,1); % Problem dimension
SS = size(swarm,2); % Swarm size
MaxIt = 100000; % Maximum number of iterations
chi = 0.729; % Constriction coefficient
c1 = 2.05; % Cognitive parameter
c2 = 2.05; % Social parameter
NR = 1; % Neighborhood radius
%===================================
    % PARAMETER INITIALIZATION
iter = 1; % Iteration counter
STOP = 0; % Stopping flag
sqrteps=sqrt(eps);
    % SWARM, VELOCITIES & BEST POSITION INITIALIZATION
vel = rand(dim, SS)*0.0000001 - Xmin;
bestpos = swarm;
fmin_counter=zeros(MaxIt,1);
xopt_counter=cell(MaxIt,1);

    % EVALUATE SWARM
for i=1:SS,
    fswarm(i) = norm(C*swarm(:,i)-d,2)^2;
    fbestpos(i) = fswarm(i);
end
fmin_counter(1,1)=min(fbestpos);
    % UPDATE OVERALL BEST PARTICLE INDICES
[fxopt, g_over] = min(fbestpos);
fS = bestpos(:,g_over);
xopt_counter{1,1}=fS;
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
      for i=1:SS,
          fswarm(i) = norm(C*swarm(:,i)-d,2)^2;
      end
        % UPDATE BEST POSITIONS
      bestpos(:,fswarm<fbestpos) = swarm(:,fswarm<fbestpos);
      fbestpos(fswarm<fbestpos) = fswarm(fswarm<fbestpos);
      fmin_counter(iter,1)=min(fbestpos);
        % UPDATE OVERALL BEST PARTICLE INDICES
      [fxopt, g_over] = min(fbestpos);
      fS = bestpos(:,g_over);     
      xopt_counter{iter,1}=fS;
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
      elseif mod(iter,300)==1 && ((fmin_counter(iter-10,1)-min(fbestpos))<=tol)&&(max(max(abs(fS-xopt_counter{iter-10,:}) / (sqrteps+max(max(abs(xopt_counter{iter-10,:})))))))<=tol
         STOP = 1; 
      end  
end % FINISH SWARM EVOLUTION
%==========================
end