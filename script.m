%% Filippini 2023 - Test Cases %%
clc, clear, close all

% Selection of test case
% 'A' - disc/sphere with absorbing outer boundary
% 'B' - disc/sphere with semi-absorbing outer boundary
% 'C' - annulus/spherical shell with absorbing outer boundary
% 'D' - annulus/spherical shell with semi-absorbing outer boundary
% 'E' - annulus/spherical shell with two absorbing boundaries
% 'F' - disc/sphere with absorbing outer boundary (two layers)
% 'G' - disc/sphere with semi-absorbing outer boundary (two layers)

Case = 'C'; % case type

% Model parameters
Nr = 501; % no. of spatial nodes (continuum model)
Nt = 1e4; % no. of time steps (continuum model)
Np1 = 50; Np2 = 500; % no. of particles (stochastic model)
Ns = 100; % no. of simulations (stochastic model)
delta = 1; tau = 1; % step size and duration (stochastic model)
c0 = 1; % initial condition (continuum model)
models = {'single','two','weight'}; % surrogate model types
c_exp = cell(3,1); % array for surrogate models

% Test cases for two and three dimensions (radially-symmetric geometries)
for d = 2:3 
    % Test case parameters
    [P,IB,OB,IF] = get_case_parameters(Case,d,delta);

    % Continuum model
    D = P*delta^2/(2*d*tau); % mass diffusivity
    [~,lambda] = compute_surrogate_model(d,D,IB,OB,'single',IF); % compute rate parameter for one term exponential model
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
        minMaxNp1 = compute_min_max_intervals(Ns,d,IB,OB,Np1,T,tau,delta,P);
        minMaxNp2 = compute_min_max_intervals(Ns,d,IB,OB,Np2,T,tau,delta,P);
    else % heterogeneous problem
        if d == 2 % choose the number of partitions
            parts = 36;
        else
            parts = 12;
        end
        minMaxNp1 = compute_min_max_intervals(Ns,d,IB,OB,Np1,T,tau,delta,P,L1,parts); % compute the min-max intervals (two layer)
        minMaxNp2 = compute_min_max_intervals(Ns,d,IB,OB,Np2,T,tau,delta,P,L1,parts); % compute the min-max intervals (two layer)
    end

    ts = 0:tau:T; % stochastic model time points
    tSetup = [ts ts(end:-1:1)]; % time points setup for plotting
    minMaxSetupNp1 = [minMaxNp1(1,:) minMaxNp1(2,end:-1:1)]; % confidence intervals (N1 = 50)
    minMaxSetupNp2 = [minMaxNp2(1,:) minMaxNp2(2,end:-1:1)]; % confidence intervals (N2 = 500)
    
    % Plotting %
    figure
    hold on
    
    % Confidence interval for N = 50 particles
    minMaxPlotNp1 = fill(tSetup,minMaxSetupNp1,'','HandleVisibility','off');
    minMaxPlotNp1.FaceColor = "#758da3";
    minMaxPlotNp1.FaceAlpha = 0.2;
    minMaxPlotNp1.EdgeColor = 'none';
    
    % Confidence interval for N = 500 particles
    minMaxPlotNp2 = fill(tSetup,minMaxSetupNp2,'','HandleVisibility','off');
    minMaxPlotNp2.FaceColor = "#758da3";
    minMaxPlotNp2.FaceAlpha = 0.4;
    minMaxPlotNp2.EdgeColor = 'none';

    % Spatial average
    plot(tc,c_avg,'LineWidth',4,'Color','#05386B')

    % Surrogate models
    tc = 0:T/Nt:T; % continuum model time points
    if isscalar(D) % three surrogate models for homogeneous media
        % Single-parameter exponential model
        [c_exp{1},lambda,~] = compute_surrogate_model(d,D,IB,OB,'single',IF,tc);
        err = mean(abs(c_avg - c_exp{1})); % mean absolute error

        % Scientific notation for plotting %
        lambda_strTemp = sprintf('%0.2e',lambda);
        lambda_str = strcat(lambda_strTemp(1:4),' \times 10^{',...
            num2str(str2double(lambda_strTemp(6:end))),'}');

        err1_strTemp = sprintf('%0.2e',err);
        err1_str = strcat(err1_strTemp(1:4),' \times 10^{',...
            num2str(str2double(err1_strTemp(6:end))),'}');
        
        % Two-parameter exponential model
        [c_exp{2},lambda,~] = compute_surrogate_model(d,D,IB,OB,'two',IF,tc);
        err = mean(abs(c_avg - c_exp{2})); % mean absolute error

        % Scientific notation for plotting
        lambda1_strTemp = sprintf('%0.2e',lambda(1));
        lambda1_str = strcat(lambda1_strTemp(1:4),' \times 10^{',...
            num2str(str2double(lambda1_strTemp(6:end))),'}');
        
        lambda2_strTemp = sprintf('%0.2e',lambda(2));
        lambda2_str = strcat(lambda2_strTemp(1:4),' \times 10^{',...
            num2str(str2double(lambda2_strTemp(6:end))),'}');

        err2_strTemp = sprintf('%0.2e',err);
        err2_str = strcat(err2_strTemp(1:4),' \times 10^{',...
            num2str(str2double(err2_strTemp(6:end))),'}');
        
         % Three-parameter exponential model
        [c_exp{3},~,theta] = compute_surrogate_model(d,D,IB,OB,'two',IF,tc);
        err = mean(abs(c_avg - c_exp{3})); % mean absolute error
    
        % Scientific notation for plotting
        err3_strTemp = sprintf('%0.2e',err);
        err3_str = strcat(err3_strTemp(1:4),' \times 10^{',...
            num2str(str2double(err3_strTemp(6:end))),'}');

    else % Two surrogate models for heterogeneous media
        % Single-parameter exponential model
        [c_exp{1},lambda,~] = compute_surrogate_model(d,D,IB,OB,'single',IF,tc);
        err = mean(abs(c_avg - c_exp{1})); % mean absolute error

        % Scientific notation for plotting
        lambda_strTemp = sprintf('%0.2e',lambda);
        lambda_str = strcat(lambda_strTemp(1:4),' \times 10^{',...
            num2str(str2double(lambda_strTemp(6:end))),'}');

        err1_strTemp = sprintf('%0.2e',err);
        err1_str = strcat(err1_strTemp(1:4),' \times 10^{',...
            num2str(str2double(err1_strTemp(6:end))),'}');
        
        % Two-parameter exponential model
        [c_exp{2},lambda,~] = compute_surrogate_model(d,D,IB,OB,'two',IF,tc);
        err = mean(abs(c_avg - c_exp{2})); % mean absolute error

        % Scientific notation for plotting
        lambda1_strTemp = sprintf('%0.2e',lambda(1));
        lambda1_str = strcat(lambda1_strTemp(1:4),' \times 10^{',...
            num2str(str2double(lambda1_strTemp(6:end))),'}');
        
        lambda2_strTemp = sprintf('%0.2e',lambda(2));
        lambda2_str = strcat(lambda2_strTemp(1:4),' \times 10^{',...
            num2str(str2double(lambda2_strTemp(6:end))),'}');

        err2_strTemp = sprintf('%0.2e',err);
        err2_str = strcat(err2_strTemp(1:4),' \times 10^{',...
            num2str(str2double(err2_strTemp(6:end))),'}');
    end
   
    % Final time (scientific notation)
    T_strTemp = sprintf('%0.2e',T);
    T_str = strcat(T_strTemp(1:4),' \times 10^{',num2str(str2double(T_strTemp(6:end))),'}');

    % Plot surrogate models with parameters
    plot(tc,c_exp{1},'LineWidth',4,'Color','#7b4288')
    plot(tc,c_exp{2},'LineWidth',4,'LineStyle','--','Color','#ff764a')
    if isscalar(D)
        plot(tc,c_exp{3},'LineWidth',4,'LineStyle','--','Color','#f5bf03')
    else
    end
    xlim([0,T]), ylim([0,1])
end

