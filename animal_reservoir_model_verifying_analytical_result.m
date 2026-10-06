function R0 = animal_reservoir_model_verifying_analytical_result(kappa)
% This code verifies the analytical results for the mean return time and
% the mean fadeout time with the eigenvalue approach to obtaining the same
% quantities.
%
% It produces Figures 2,3 and 4 in the paper. 
%   Usage:
%       R0 = animal_reservoir_model_verifying_analytical_result(kappa)
%  
%   Returns:
%       R0:                     Basic reproduction number
%       
%       Plots:                  Figure 2 or 3, Figure 4
%
%   Required:
%       kappa:                  Choice of kappa parameter, N*force of external infection.
%   
%   Created: 18-04-2025
%   Author: Maame Ama Bainson
%   Programme version: MATLAB 2024b



    %% Parameters
    beta = 130;
    gam = 100;
    N_set = 100:50:900;
    R0 = beta/gam;
    
    MFT_set1 = NaN(2,numel(N_set));
    MFT_set2 = NaN(2,numel(N_set));
    
    for idx = 1:numel(N_set)
    
        N = N_set(idx);
    
        k = kappa; %{1e-3 or 2}
    
        Q = sisk1D_Q(beta,k,gam,N);
        Q_trunc = Q(2:end,2:end);
    
        %% Computing stationary distribution & mean return time by the eigenvalue method
        [u,~] = eigs(transpose(Q),1,'sm');
        u_3 = u./sum(u);
        mean_return = 1/(u_3(1)*k);
  
    
        %% Computing the quasistationary distribution & mean fadeout time from endemicity by the eigenvalue method
        [~,uevl_qsd] = eigs(transpose(Q_trunc),1,'sm');
        alpha = -uevl_qsd;
        MFT_eigen = 1/alpha;
    
        %% Analytical tau Eq(3.15) 
        nu = k/beta;
        numer = sqrt(2*pi); 
        denom = (N^(0.5-nu))*gamma(nu + 1)*beta; 
        tau_analytical = (numer/denom)*((R0/(R0-1))^(1-nu))*exp(N*(log(R0) - 1 + (1/R0))); 
       
        MFT_set1(:,idx) = [tau_analytical;mean_return];
    
        %% Analytical tau_z Eq(3.17)
        tau_z_analytical = (R0/(R0-1))*tau_analytical;
        MFT_set2(:,idx) = [tau_z_analytical;MFT_eigen];
    
     
    end 
   


    
    %% Plots 
    figure()
    subplot(2,1,1)
    semilogy(N_set,MFT_set1(1,:),'k.',N_set,MFT_set1(2,:),'b*','MarkerSize',15)
    xlabel('Population size, N','Interpreter','latex')
    ylabel('Mean return time','Interpreter','latex')
    legend('$\tau$','Eigenvalue result','Interpreter','latex','Location','southeast')
    legend('boxoff')
    
    subplot(2,1,2)
    semilogy(N_set,MFT_set1(1,:)./MFT_set1(2,:),'k.','MarkerSize',12)
    ylim([1e-4 1e1])
    xlabel('Population size, N','Interpreter','latex')
    ylabel('Ratio','Interpreter','latex')
    
    figure()
    subplot(2,1,1)
    semilogy(N_set,MFT_set2(1,:),'r.',N_set,MFT_set2(2,:),'k*','MarkerSize',15)
    legend('$\tau_z$','Eigenvalue result','Interpreter','latex','Location','southeast')
    legend('boxoff')
    xlabel('Population size, N', 'Interpreter','latex')
    ylabel('Mean fade out time', 'Interpreter','latex')
    
    subplot(2,1,2)
    semilogy(N_set,MFT_set2(1,:)./MFT_set2(2,:),'r.','MarkerSize',15)
    xlabel('Population size, N','Interpreter','latex')
    ylabel('Ratio','Interpreter','latex')
    


    %% Creating the transition rate matrix for the animal reservoir model
    function Q = sisk1D_Q(beta,k,gam,N)
    
        rows = repmat(2:N,[3,1]);
        
        cols = NaN(size(rows));
        
        for j = rows(1,1):rows(1,end)
            cols(:,j-1) = [j-1,j,j+1]';
        end
        
        vals = NaN(size(rows));
        
        for j = rows(1,1):rows(1,end)
            vals(:,j-1) = [(j-1)*(gam), -(j-1)*(gam)-(beta/N)*(j-1)*(N-j+1)-(k/N)*(N-j+1), (beta/N)*(j-1)*(N-j+1)+(k/N)*(N-j+1)]';
        end
        
        %% Special cases : i= 0, i = N
        i = 0;
        rows_1 = [i+1,i+1]; cols_1 = [i+1,i+2]; vals_1 = [-k,k];
        
        i = N;
        rows_N = [i+1,i+1]; cols_N = [N,N+1]; vals_N = [N*(gam),-N*(gam)];
        
        
        Q_main_r = [rows_1,rows(:)',rows_N];
        Q_main_c = [cols_1,cols(:)',cols_N];
        Q_main_v = [vals_1,vals(:)',vals_N];
        
        %% Building the matrix
        Q = sparse(Q_main_r, Q_main_c, Q_main_v);
    end
    
   

end


