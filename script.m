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

    % Continuum model
    D = P*delta^2/(2*d*tau); % mass diffusivity
    lambda = d*(d+2)*D/(OB.L1^2 + 2*OB.b1/OB.a1); % single-parameter model
    T = 3*log(10)/lambda; % final time (particle concentration 10^(-3))
    if IF == 0
        x = linspace(IB.L0,OB.L1,Nr); % spatial nodes
        [c_avg,~] = homogeneous_continuum_model(d,D,IB,OB,x,Nt,T,c0); % numerical solution
    else
        L1 = IF; % interface point
        x1 = IB.L0:0.1:L1; % spatial nodes (L0 <= x <= L1)
        x2 = L1:0.1:OB.L2; % spatial nodes (L1 <= x <= L2)
        [c_avg,~] = heterogeneous_continuum_model(d,D,IB,OB,x1,x2,Nt,T,c0); % numerical solution
    end


    % Stochastic model
    if IF == 0 % homogeneous problem
        propConfintN1 = compute_min_max_intervals(Ns,d,IB,OB,Np1,T,tau,delta,P);
        propConfintN2 = compute_min_max_intervals(Ns,d,IB,OB,Np2,T,tau,delta,P);
    else % heterogeneous problem
        if d == 2 % choose the number of partitions
            parts = 36;
        else
            parts = 12;
        end
        propConfint = compute_min_max_intervals(Ns,d,IB,OB,Np,T,tau,delta,P,L1,parts); % compute the min-max intervals (two layer)
    end
    
    ts = 0:tau:T; % time points
    tConf = [ts ts(end:-1:1)]; % time points setup for plotting
    confSetupN1 = [propConfintN1(1,:) propConfintN1(2,end:-1:1)]; % confidence intervals (N1 = 50)
    confSetupN2 = [propConfintN2(1,:) propConfintN2(2,end:-1:1)]; % confidence intervals (N2 = 500)
    
    figure
    hold on
    
    % Confidence interval for N = 50 particles
    confPlotN1 = fill(tConf,confSetupN1,'','HandleVisibility','off');
    confPlotN1.FaceColor = "#758da3";
    confPlotN1.FaceAlpha = 0.2;
    confPlotN1.EdgeColor = 'none';
    
    % Confidence interval for N = 500 particles
    confPlotN2 = fill(tConf,confSetupN2,'','HandleVisibility','off');
    confPlotN2.FaceColor = "#758da3";
    confPlotN2.FaceAlpha = 0.4;
    confPlotN2.EdgeColor = 'none';
    
    % Spatial average
    plot(0:T/Nt:T,c_avg,'LineWidth',4,'Color','#05386B')
    xlim([0,T]), ylim([0,1])
end

