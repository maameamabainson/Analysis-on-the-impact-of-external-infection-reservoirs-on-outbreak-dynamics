function tok = optimal_trajectory_plotting(extinct_path, end_eq, df_eq,T)
% Underlying code used to plot optimal trajectories in Figures 6 in the paper. 
%
% Run with 'optimal_trajectory_computation_environmental_model.m'
%
%   Usage:
%       tok = optimal_trajectory_plotting(extinct_path, end_eq, df_eq,T) 
%
%   Returns:
%       tok:                    Computational time
%       
%
%   Required:
%       extinct_path:           Computed extinction optimal trajectory.
%       end_eq:                 Corresponding endemic equilibrium (start point)
%       df_eq:                  Corresponding disease-free equilibrium (end point)
%       T:                      Absolute value of truncated/transformed time endpoint
%
%   
%   Created: 26-02-2026
%   Author: Maame Ama Bainson
%   Programme version: MATLAB 2024b




try
    morepoints = linspace(-T,T,10000); 
    extinct_path_grid = deval(extinct_path,morepoints);
catch
    morepoints = extinct_path.x; 
    extinct_path_grid = extinct_path.y;
end


tic


        figure
        subplot(1,2,1)
        hold on
        plot(extinct_path_grid(1,:),extinct_path_grid(2,:),'k',LineWidth=2)
        xlabel('x_i')
        ylabel('x_w')
        plot(end_eq(1),end_eq(2),'r.','MarkerSize',8)
        plot(df_eq(1),df_eq(2),'g.','MarkerSize',8)
        legend('','Endemic','Disease-free', 'Location','southeast')
        legend('boxoff')
        
        subplot(1,2,2)
        hold on
        plot(extinct_path_grid(3,:),extinct_path_grid(4,:),'k',LineWidth=2)
        xlabel('\theta_i')
        ylabel('\theta_w')
        plot(end_eq(3),end_eq(4),'r.','MarkerSize',8)
        plot(df_eq(3),df_eq(4),'g.','MarkerSize',8)
        

        figure
        subplot(2,2,1)
        hold on
        plot(morepoints,extinct_path_grid(1,:),'k',LineWidth=2)
        xlabel('Time')
        ylabel('x_i')
        plot(extinct_path.x(1),end_eq(1),'r.','MarkerSize',15)
        plot(extinct_path.x(end),df_eq(1),'r*','MarkerSize',12)
        legend('','x_i*','x_i-df')

        subplot(2,2,2)
        hold on
        plot(morepoints,extinct_path_grid(2,:),'k',LineWidth=2)
        xlabel('Time')
        ylabel('x_w')
        plot(extinct_path.x(1),end_eq(2),'g.','MarkerSize',15)
        plot(extinct_path.x(end),df_eq(2),'g*','MarkerSize',12)
        legend('','x_w*','x_w-df')

        subplot(2,2,3)
        hold on
        plot(morepoints,extinct_path_grid(3,:),'k',LineWidth=2)
        xlabel('Time')
        ylabel('\theta_i')
        plot(extinct_path.x(1),end_eq(3),'k.','MarkerSize',15,'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980])
        plot(extinct_path.x(end),df_eq(3),'k*','MarkerSize',12,'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980])
        legend('','\theta_i-end','\theta_i*')


        subplot(2,2,4)
        hold on
        plot(morepoints,extinct_path_grid(4,:),'k',LineWidth=2)
        xlabel('Time')
        ylabel('\theta_w')
        plot(extinct_path.x(1),end_eq(4),'b.','MarkerSize',15)
        plot(extinct_path.x(end),df_eq(4),'b*','MarkerSize',12)
        legend('','\theta_w-end','\theta_w*')

tok = toc;

end