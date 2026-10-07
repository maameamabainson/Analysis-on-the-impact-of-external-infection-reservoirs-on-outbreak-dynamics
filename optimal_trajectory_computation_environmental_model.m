function [data, extinct_path, runtime] = optimal_trajectory_computation_environmental_model(algorithm)
% This code computes the optimal trajectory (solution of the equations of motion) for the SIS model with
% environmental reservoir using either time truncation t in [-\infty, \infty] ->
% t' in [-T,T] or time transformation where t' = tanh(At).
%
% Run with 'optimal_trajectory_plotting.m'
% 
%
%   Usage:
%       [data, extinct_path, runtime] = optimal_trajectory_computation_environmental_model('CTs')
%       [data, extinct_path, runtime] = optimal_trajectory_computation_environmental_model('CTr')
%  
%
%   Returns:
%       data:                  [varying parameter, R0, action, max_Hmt, endemic_distance, df_distance, path_length]
%       
%
%   Required:
%       algorithm:             Algorithm key: 'CTs' for time transformation;
%                                             'CTr' for time truncation.
%                             
%
%   References:
%   Clancy, D. and Stewart, J.J., 2025. Computing the extinction path for epidemic models. Mathematical Biosciences, 386, p.109454.
%
%
%
%   Created: 05-02-2026
%   Author: Maame Ama Bainson
%   Programme version: MATLAB 2024b



switch algorithm 
    case 'CTr'
        [data, extinct_path, runtime] = sisv_coll_trunc;
    case 'CTs'
        [data, extinct_path, runtime] = sisv_coll_transf;
    otherwise
        error('Invalid Algorithm type!')
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
end




%% BVP5c WITH TIME TRANSFORMATION
function [data, extinct_path, runtime] = sisv_coll_transf

tic 

%% Setting solver options 
options = bvpset('Vectorized','on','AbsTol',1e-04);

PP = 0; QQ = 0; %policy parameters; PP: percentage of strategy given to preventive measures; QQ: percentage of strategy given to curative measures

%% Parameter values

lambda = 0.55;
alpha = 0.35;
gam = 1*(1+(0.2*QQ));
phi = 5;

set = 1.1*(1-(0.2*PP));      


%Grid for solving bvp
gridnum = 100; 


%Initializing data to receive [Parameter, R0, action, max_Hmt, endemic_distance, df_distance, path_length]
data = zeros(numel(set),7);   

%For each parameter set
for iter_no = 1:numel(set)
    
    if iter_no == 1
       
        beta = set(iter_no);
        [R0, end_eq, df_eq, sisv_ode, Hmt] = sisv_bvp(beta, phi, lambda, alpha, gam);


        %R0 less than 1, warning! else proceed with function
        assert(R0>1,'R0 must be greater than 1.')
   
        %Boundary conditions for solving bvp to obtain extinction path
        sisv_bc = @(ya,yb) [ ya(1:2) - end_eq(1:2)
                             yb(3:4) - df_eq(3:4)];

        T = 1; %End point of time transformed interval
        
        t_prime = linspace(-T,T,gridnum);
        
        a_tilde = 0.045; %selection from various trials
        distance = sqrt(transpose(end_eq-df_eq)*(end_eq-df_eq));
        a = a_tilde*(distance);
        
        trans_sisv_ode  = trans_t(sisv_ode,a);

        %init_guess = @(t_tilde) ivp_check(t_tilde,beta,gam, lambda, alpha, df_eq, end_eq,T);
        init_guess = @(t) (df_eq.*(t+abs(T)) - end_eq*(t-abs(T))) ./ (2*abs(T));

        solinit = bvpinit(t_prime,init_guess);
       
        extinct_path = bvp5c(trans_sisv_ode,sisv_bc,solinit,options);
     
        %Compute action integral and convergence diagnostics
        [action, path_length] = path_integrals(extinct_path.y);
        extinct_path_grid = deval(extinct_path,linspace(-T,T,gridnum));
        max_Hmt = max(abs(Hmt(extinct_path_grid)));
        endemic_distance = sqrt((extinct_path.y(:,1)-end_eq)'*(extinct_path.y(:,1)-end_eq));
        df_distance = sqrt((extinct_path.y(:,end)-df_eq)'*(extinct_path.y(:,end)-df_eq));

        %Display results
        fprintf("a = %.8g\n",a)
        fprintf("beta = %.8g\n",beta)
        fprintf("Basic reproduction number R0 = %.6g\n",R0)
        fprintf("bvp5c action integral = %.8g\n",action)
        fprintf("Maximal absolute value of Hamiltonian along path = %.3g\n",max_Hmt)
        fprintf("Distance from endemic point = %.3g\n",endemic_distance)
        fprintf("Distance from disease free point = %.3g\n",df_distance)
        fprintf("Path length = %.3g\n",path_length)
        fprintf("\n")
        

        %Store results
        data(iter_no,:)=[beta R0 action max_Hmt endemic_distance df_distance path_length];

    else 
        
        beta = set(iter_no);

        [R0, end_eq, df_eq, sisv_ode, Hmt] = sisv_bvp(beta, phi, lambda, alpha, gam);


        %R0 less than 1, warning! else proceed with function
        assert(R0>1,'R0 must be greater than 1.')
   
        %Boundary conditions for solving bvp to obtain extinction path
        sisv_bc = @(ya,yb) [ ya(1:2) - end_eq(1:2)
                             yb(3:4) - df_eq(3:4)];

        distance = sqrt(transpose(end_eq-df_eq)*(end_eq-df_eq));
        a_new = (a_tilde+0.005)*(distance);
        
        t_prime_old = extinct_path.x;

        t_prime_new = tanh(a.*atanh(t_prime_old)./a_new);  

        extinct_path.x = t_prime_new;

        a = a_new;

        trans_sisv_ode  = trans_t(sisv_ode,a);
        
        solinit = extinct_path;

        extinct_path = bvp5c(trans_sisv_ode,sisv_bc,solinit,options);
        
        %Compute action integral and convergence diagnostics
        [action, path_length] = path_integrals(extinct_path.y);
        extinct_path_grid = deval(extinct_path,linspace(-T,T,gridnum));
        max_Hmt = max(abs(Hmt(extinct_path_grid)));
        endemic_distance = sqrt((extinct_path.y(:,1)-end_eq)'*(extinct_path.y(:,1)-end_eq));
        df_distance = sqrt((extinct_path.y(:,end)-df_eq)'*(extinct_path.y(:,end)-df_eq));

        %% Display results
        fprintf("a = %.8g\n",a)
        fprintf("beta = %.8g\n",beta)
        fprintf("Basic reproduction number R0 = %.6g\n",R0)
        fprintf("bvp5c action integral = %.8g\n",action)
        fprintf("Maximal absolute value of Hamiltonian along path = %.3g\n",max_Hmt)
        fprintf("Distance from endemic point = %.3g\n",endemic_distance)
        fprintf("Distance from disease free point = %.3g\n",df_distance)
        fprintf("Path length = %.3g\n",path_length)
        fprintf("\n")
       

        %Store results
        data(iter_no,:)=[beta R0 action max_Hmt endemic_distance df_distance path_length];
    end
end
    
        optimal_trajectory_plotting(extinct_path, end_eq, df_eq,T);  
        
        runtime = toc;

end





%% BVP5c WITH TIME TRUNCATION 
function [data, extinct_path, runtime] = sisv_coll_trunc

tic 

%setting solver options 
options = bvpset('Vectorized','on','AbsTol',1e-04,'RelTol',1e-05,Nmax=50000);


%% Parameter values (Figure 6)
beta = 60;
lambda = 2755.75;
alpha = 43.07;
gam = 52.195;

set = [0.00001,0.001];

%{
%% Parameter values (Figure 7)
beta = 36.7;
lambda = 50;
alpha = 43;
gam = 35;

set = 0.001:3700;
%}

%Grid for solving bvp
gridnum = 2; 

%Initializing data to receive [Parameter,R0, action, max_Hmt,endemic_distance, df_distance, path_length]
data = zeros(numel(set),7);


%For each parameter set
for iter_no = 1:numel(set)

    phi = set(iter_no);
    %beta = set(iter_no);

   [R0, end_eq, df_eq, sisv_ode, Hmt] = sisv_bvp(beta, phi, lambda, alpha, gam);


    %R0 less than 1, warning! else proceed with function
    assert(R0>1,'R0 must be greater than 1.')

    %Boundary conditions for solving bvp to obtain extinction path
    sisv_bc = @(ya,yb) [ ya(1:2) - end_eq(1:2)
                         yb(3:4) - df_eq(3:4)];

    if iter_no == 1 
    
    %% Obtaining extinction path for first set of parameters with time extension  
    T_set = 1:0.01:5;
            
    for t_num = 1:numel(T_set)
        T = T_set(t_num);

        if t_num == 1 %setting up the extinction path for the first time truncation. 

            %init_guess = @(t) (df_eq.*(t+abs(T)) - end_eq*(t-abs(T))) ./ (2*abs(T)); 
            init_guess = @(t) ivp_check(t,beta,gam, lambda, alpha, df_eq, end_eq,T);

            solinit = bvpinit(linspace(-T,T,gridnum),init_guess);
            
        else 

            %After first truncation, we proceed to time extension at each end and use it as the initial guess for the subsequent truncation ranges. 

            ext_left = bvpxtend(extinct_path,-T,end_eq);
            solinit= bvpxtend(ext_left,T,df_eq);
        
        end

        extinct_path = bvp5c(sisv_ode,sisv_bc,solinit,options);
    end
    
    else 
    
    %Parameter continuation
    solinit = extinct_path;
    extinct_path = bvp5c(sisv_ode,sisv_bc,solinit,options);

    %Making sure computed path is defined over the whole interval [-T,T]
    new_ext_left = bvpxtend(extinct_path, -T, end_eq);
    extinct_path = bvpxtend(new_ext_left, T, df_eq);
    end
    
    %Action integral and convergence diagnostics
    [action, path_length] = path_integrals(extinct_path.y);
    extinct_path_grid = deval(extinct_path,linspace(-T,T,gridnum));
    max_Hmt = max(abs(Hmt(extinct_path_grid)));
    endemic_distance = sqrt((extinct_path.y(:,1)-end_eq)'*(extinct_path.y(:,1)-end_eq));
    df_distance = sqrt((extinct_path.y(:,end)-df_eq)'*(extinct_path.y(:,end)-df_eq));

    %% Display results
    fprintf("beta = %.8g\n", beta)
    fprintf("phi = %.8g\n", phi)
    fprintf("Basic reproduction number R0 = %.6g\n",R0)
    fprintf("bvp5c action integral = %.8g\n",action)
    fprintf("Maximal absolute value of Hamiltonian along path = %.3g\n",max_Hmt)
    fprintf("Distance from endemic point = %.3g\n",endemic_distance)
    fprintf("Distance from disease free point = %.3g\n",df_distance)
    fprintf("Path length = %.3g\n",path_length)
    fprintf("\n")
  

    %Store the results
    data(iter_no,:)=[beta R0 action max_Hmt endemic_distance df_distance path_length];
   

end


  
  %% Trajectory plots 
  optimal_trajectory_plotting(extinct_path, end_eq, df_eq,T)
       
  runtime = toc;

end






%% SIS_phi model essentials: R0, equilibrium points, equations of motion, Hamiltonian expression
function [R0, end_eq, df_eq, sisv_ode, Hmt] = sisv_bvp(beta, phi, lambda, alpha, gam)
R0 = (beta/gam) + lambda*phi/(gam*phi+gam*alpha);

xi_star = 1-(1/R0);
xv_star = (lambda/(phi+alpha))*xi_star;

end_eq = [xi_star
          xv_star
          0
          0];

L = gam*(alpha + phi)/(alpha*beta + beta*phi + lambda*phi);
M = (beta*alpha^2 + alpha*beta*phi + alpha*lambda*phi + alpha*gam*phi + gam*phi^2)/((alpha + phi)*(alpha*beta + beta*phi + lambda*phi));

theta_i = log(L); 
theta_v = log(M);

df_eq = [0 
         0
         theta_i
         theta_v];



sisv_ode = @(t,y)[beta.*y(1,:).*(1-y(1,:)).*exp(y(3,:)) + phi.*(1-y(1,:)).*y(2,:).*exp(y(3,:)-y(4,:)) - gam.*y(1,:).*exp(-y(3,:));
                  
                 lambda.*y(1,:).*exp(y(4,:)) - phi.*(1-y(1,:)).*y(2,:).*exp(y(3,:)-y(4,:)) - (phi.*y(1,:) + alpha).*y(2,:).*exp(-y(4,:));

                 (2.*beta.*y(1,:) - beta).*(exp(y(3,:))-1) + phi.*y(2,:).*(exp(y(3,:)-y(4,:))-1) - gam.*(exp(-y(3,:))-1) - lambda.*(exp(y(4,:))-1) - phi.*y(2,:).*(exp(-y(4,:))-1);

                 (phi.*y(1,:) - phi).*(exp(y(3,:)-y(4,:))-1) - (phi.*y(1,:) + alpha).*(exp(-y(4,:))-1)];


Hmt = @(y) beta.*y(1,:).*(1-y(1,:)).*(exp(y(3,:))- 1) + phi.*(1-y(1,:)).*y(2,:).*(exp(y(3,:)-y(4,:)) - 1) + gam.*y(1,:).*(exp(-y(3,:))-1) + lambda.*y(1,:).*(exp(y(4,:))-1) + (alpha + phi.*y(1,:)).*y(2,:).*(exp(-y(4,:))-1); 


end





%% Only for 'CTs' option: transform the time variable usinf tanh
function odefuntrans = trans_t(odefun,A)
    % Inputs:
    %	odefun	        - characteristic equation derivatives
    %	A		        - scaling factor
    %
    % Outputs:
    %	odefuntrans     - characteristic equations in transformed time
    
    odefuntrans = @(t,x) ode_trans(t,x);
    
    function dxdt = ode_trans(tau,x)
    % Function representing the characteristic equations after the transformation t' = tanh(At) 
    
    dxdt = odefun(atanh(tau)./A,x) ./ A ./ (1-tau.^2);
    
    % Avoid division by zero errors at boundaries by replacing with zero derivatives.
    % Note that time-transformed derivatives at boundaries may not be zero;
    % check diagnostic values to ensure computed solution is sensible.
    indices = find(abs(tau)==1); 
    dxdt(:,indices)=zeros(size(x,1),length(indices));
    
    end

end





%% Compute action integral and path length

function [A, L] = path_integrals(ext_path_y)
    A=0;

    for i=1:size(ext_path_y,1)/2
        A=A+trapz(ext_path_y(i,:),ext_path_y(end/2+i,:));
    end
    
    L = sum(sqrt(sum((ext_path_y(:,2:end)-ext_path_y(:,1:end-1)).^2,1)));
end





%% Compute the initial guess for solver using partial deterministic solution
function init = ivp_check(t_tilde,beta,gamma, lambda, alpha, df_eq, end_eq,T)
        t = -t_tilde + T;

        init = NaN(4,numel(t));
        r = beta - gamma;
        k = r/beta;
        x0 = 1e-5; xv0 = 0;

        init(1,:) = (k.*x0.*exp(r.*t))./(k-x0+x0.*exp(r.*t));

        init(2,:) = (xv0 - x0.*(lambda/alpha)).*exp(-alpha.*t) + init(1,:).*(lambda/alpha);

        init(3,:) = log(gamma/(beta*(1-init(1,:))));

        init(4,:) = (df_eq(4).*(t_tilde+abs(T)) - end_eq(4)*(t_tilde-abs(T))) ./ (2*abs(T));

end







