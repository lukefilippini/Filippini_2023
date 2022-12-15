%% Filippini 2023 - Test Cases %%
clc, clear, close all

% Selection of test case
% 'A' - disc/sphere with absorbing outer boundary
% 'B' - disc/sphere with semi-absorbing outer boundary
% 'C' - annulus/spherical shell with absorbing outer boundary
% 'D' - annulus/spherical shell with semi-absorbing outer boundary
% 'E' - annulus/spherical shell with two absorbing boundaries
% 'F' - annulus/spherical shell with two semi-absorbing boundaries
% 'G' - disc/sphere with absorbing outer boundary (two layers)
% 'H' - disc/sphere with semi-absorbing outer boundary (two layers)

Case = 'A'; % case type

% Model parameters
Nr = 501; % no. of spatial nodes (continuum model)
Nt = 1e4; % no. of time steps (continuum model)
Np1 = 50; Np2 = 500; % no. of particles (stochastic model)
Ns = 100; % no. of simulations (stochastic model)
delta = 1; tau = 1; % step size and duration (stochastic model)
c0 = 1; % initial condition (continuum model)

% Test cases for two and three dimensions (radially-symmetric geometries)
for d = 2:3

    % Test case parameters
    [P,IB,OB,IF] = get_case_parameters(Case,d,delta);
    
    % Generate final time from single-exponential rate parameter
    T = 34000;

    % Continuum model
    D = P*delta^2/(2*d*tau); % mass diffusivity
    if IF == 0
        x = linspace(IB.L0,OB.L1,Nr); % spatial nodes
        c = homogeneous_continuum_model(d,D,IB,OB,x,Nt,T,c0);
        if (IB.L0 == 0) && (d > 1) 
            c = [c(1,:);c];
        end
        dt = T/Nt; tspan = logspace(-3.5,0,8)*T; % step size and time points for plotting
        indx = floor(tspan/dt + 1); % columns in solution matrix
        indigo = [20/255,47/255,103/255]; % plotting colour
        
        figure
        hold on
        for n = indx % plot numerical solution
            plot(x,c(:,n),'LineWidth',2.5,'Color',indigo)
        end
        xlabel('$x$','FontSize',20,'Interpreter','latex')
        ylabel('$c(x,t)$','FontSize',20,'Interpreter','latex')
        xlim([IB.L0,OB.L1]), ylim([0,1])
    else
        L1 = IF; % interface point
    end

    % Stochastic model
end

