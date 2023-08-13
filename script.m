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
Nr = 501; % no. of spatial nodes (homogeneous continuum model)
% Nr = 1001; % no. of spatial nodes (heterogeneous continuum model)
Nt = 1e4; % no. of time steps (continuum model)
Np1 = 50; Np2 = 500; % no. of particles (stochastic model)
Ns = 100; % no. of simulations (stochastic model)
delta = 1; tau = 1; % step size and duration (stochastic model)
c0 = 1; % initial condition (continuum model)
models = {'one-term','two-term','weighted'}; % surrogate model types
s_exp = cell(3,1); % array for surrogate models

% Test cases for two and three dimensions (radially-symmetric geometries)
for d = 2
    % Test case parameters
    [P,IB,OB,IF] = get_case_parameters(Case,delta);

    % Continuum model
    tic
    D = P*delta^2/(2*d*tau); % mass diffusivity
    [~,lambda] = compute_surrogate_model(d,D,IB,OB,'one-term',IF); % compute rate parameter for one term exponential model
    T = 2*log(10)/lambda; % final time (particle concentration 10^(-3))
    if IF == 0
        x = linspace(IB.L0,OB.L1,Nr); % spatial nodes
        [Pc,~] = homogeneous_continuum_model(d,D,IB,OB,x,Nt,T,c0); % numerical solution
    else
        L1 = IF; % interface point
        x1 = linspace(IB.L0,L1,ceil(Nr*L1/OB.L2)); % spatial nodes (L0 <= x <= L1)
        x2 = linspace(L1,OB.L2,ceil(Nr*(OB.L2-L1)/OB.L2)); % spatial nodes (L1 <= x <= L2)
        [Pc,~] = heterogeneous_continuum_model(d,D,IB,OB,x1,x2,Nt,T,c0); % numerical solution
    end
    fprintf("Continuum model: ")
    toc

    % Stochastic model
    tic
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
    fprintf("Stochastic model: ")
    toc

    ts = 0:tau:T; % stochastic model time points
    tSetup = [ts ts(end:-1:1)]; % time points setup for plotting
    minMaxSetupNp1 = [minMaxNp1(1,:) minMaxNp1(2,end:-1:1)]; % confidence intervals (N1 = 50)
    minMaxSetupNp2 = [minMaxNp2(1,:) minMaxNp2(2,end:-1:1)]; % confidence intervals (N2 = 500)

    % Surrogate models
    tic
    tc = 0:T/Nt:T; % continuum model time points
    if isscalar(D) % three surrogate models for homogeneous media
        % Single-parameter exponential model
        [s_exp{1},lambda,~] = compute_surrogate_model(d,D,IB,OB,'one-term',IF,tc);
        err = mean(abs(Pc - s_exp{1})); % mean absolute error

        % Scientific notation for plotting %
        err1_strTemp = sprintf('%0.2e',err);
        err1_str = strcat(err1_strTemp(1:4),' \times 10^{',...
            num2str(str2double(err1_strTemp(6:end))),'}');
        
        % Two-parameter exponential model
        [s_exp{2},lambda12,~] = compute_surrogate_model(d,D,IB,OB,'two-term',IF,tc);
        err = mean(abs(Pc - s_exp{2})); % mean absolute error

        % Scientific notation for plotting
        err2_strTemp = sprintf('%0.2e',err);
        err2_str = strcat(err2_strTemp(1:4),' \times 10^{',...
            num2str(str2double(err2_strTemp(6:end))),'}');
        
         % Three-parameter exponential model
        [s_exp{3},lambdaw,theta] = compute_surrogate_model(d,D,IB,OB,'weighted',IF,tc);
        err = mean(abs(Pc - s_exp{3})); % mean absolute error
    
        % Scientific notation for plotting
        err3_strTemp = sprintf('%0.2e',err);
        err3_str = strcat(err3_strTemp(1:4),' \times 10^{',...
            num2str(str2double(err3_strTemp(6:end))),'}');

    else % Two surrogate models for heterogeneous media
        % Single-parameter exponential model
        [s_exp{1},lambda,~] = compute_surrogate_model(d,D,IB,OB,'one-term',IF,tc);
        err = mean(abs(Pc - s_exp{1})); % mean absolute error

        % Scientific notation for plotting
        err1_strTemp = sprintf('%0.2e',err);
        err1_str = strcat(err1_strTemp(1:4),' \times 10^{',...
            num2str(str2double(err1_strTemp(6:end))),'}');
        
        % Two-parameter exponential model
        [s_exp{2},lambda12,~] = compute_surrogate_model(d,D,IB,OB,'two-term',IF,tc);
        err = mean(abs(Pc - s_exp{2})); % mean absolute error

        % Scientific notation for plotting
        err2_strTemp = sprintf('%0.2e',err);
        err2_str = strcat(err2_strTemp(1:4),' \times 10^{',...
            num2str(str2double(err2_strTemp(6:end))),'}');
    end
    fprintf("Surrogate models: ")
    toc
   
    % Final time (scientific notation)
    T_strTemp = sprintf('%0.2e',T);
    T_str = strcat(T_strTemp(1:4),' \times 10^{',num2str(str2double(T_strTemp(6:end))),'}');
    
    % Plotting %
    figure
    set(gcf,'Renderer','Painters');
    hold on
    box on
    

    % Confidence interval for N = 50 particles
    minMaxPlotNp1 = fill(tSetup,minMaxSetupNp1,'','HandleVisibility','off');
    minMaxPlotNp1.FaceColor = "#758da3";
    minMaxPlotNp1.FaceAlpha = 0.3;
    minMaxPlotNp1.EdgeColor = 'none';
    
    % Confidence interval for N = 500 particles
    minMaxPlotNp2 = fill(tSetup,minMaxSetupNp2,'','HandleVisibility','off');
    minMaxPlotNp2.FaceColor = "#758da3";
    minMaxPlotNp2.FaceAlpha = 0.8;
    minMaxPlotNp2.EdgeColor = 'none';

    % Spatial average
    plot(tc,Pc,'LineWidth',6.5,'Color','#05386B')

    % Plot surrogate models with parameters
    plot(tc,s_exp{1},'LineWidth',6.5,'Color','#7B4288')
    plot(tc,s_exp{2},'LineWidth',6.5,'LineStyle','--','Color','#FF764A')
    if isscalar(D)
        plot(tc,s_exp{3},'LineWidth',6.5,'LineStyle','--','Color','#F5BF03')
        text(0.4*T,0.9,['Case ',Case],'FontSize',50,'Interpreter','latex')
        text(0.4*T,0.78,strcat('$\varepsilon_1 = ',err1_str,'$'),...
            'FontSize',40,'Interpreter','latex')
        text(0.4*T,0.66,strcat('$\varepsilon_2 = ',err2_str,'$'),...
            'FontSize',40,'Interpreter','latex')
        text(0.4*T,0.54,strcat('$\varepsilon_3 = ',err3_str,'$'),...
            'FontSize',40,'Interpreter','latex')
        text(0.4*T,0.42,strcat('$T = ',T_str,'$'),'FontSize',40,'Interpreter','latex')
        text(0.05*T,0.1,strcat('$\rm (',lower(Case),')$'),'FontSize',50,'Interpreter','latex')
    else
        text(0.4*T,0.9,['Case ',Case],'FontSize',50,'Interpreter','latex')
        text(0.4*T,0.78,strcat('$\varepsilon_1 = ',err1_str,'$'),...
            'FontSize',40,'Interpreter','latex')
        text(0.4*T,0.66,strcat('$\varepsilon_2 = ',err2_str,'$'),...
            'FontSize',40,'Interpreter','latex')
        text(0.4*T,0.54,strcat('$T = ',T_str,'$'),'FontSize',40,'Interpreter','latex')
        
        if isequal(Case,'F')
            text(0.05*T,0.1,'$\rm (a)$','FontSize',50,'Interpreter','latex')
        else
            text(0.05*T,0.1,'$\rm (b)$','FontSize',50,'Interpreter','latex')
        end
    end

    ylabel('$\mathcal{P}(t)$','FontSize',20,'Interpreter','latex')
    xlim([0,T]), ylim([0,1])

    % Axis
    axis square
    set(gca,'FontSize',50,'XTick',[0,T/2,T],'YTick',[0,1],'TickLabelInterpreter','latex')
    ax = gca; 
    ax.XTickLabel{1} = "0"; ax.XTickLabel{2} = "$t$"; ax.XTickLabel{3} = "$T$"; % axis ticks
    %axLab = get(gca,'XLabel'); axLab.Position = axLab.Position - 0.2;
    set(gcf,'Position',[1,49,1100,878])
end
