%% Thesis - Comparison of analytical and numerical moments %%
clc, clear all
%% Zeroth moment M0(r)
L0 = 50; L1 = 100; % inner and outer boundaries
a0 = 0; b0 = 1; % inner boundary condition
a1 = 1; b1 = 2; % outer boundary condition
P = 1; delta = 1; tau = 1; % stochastic parameters
r = linspace(0.001,L1,5000); % radial coordinates
beta1 = b1/a1;

%% Zeroth moment
% Initial guess
Minit = bvpinit(r,[1;1]);

% Numerical and analytical solution
for d = 1:3
    D = P*delta^2/(2*d*tau); % diffusivity
    M0 = bvp4c(@(r,M) zeroth_ode(d,D,r,M),@bound_moments,Minit);
    M0 = deval(M0,r);
    M0c = (L1^2 - r.^2 + 2*beta1*L1)/(2*d*D);
   
    figure
    hold on
    box on
    plot(r,M0(1,:),'LineWidth',6.5,'Color','#05386B')
    plot(r,M0c,'LineStyle','--','LineWidth',6.5,'Color','#FF764A')
    y = ylabel('$\vspace{1cm}M_0(r)$','FontSize',20,'Interpreter','latex');
    xlim([0,L1]), ylim([0,12000])
    legend('Numerical','Analytical','Interpreter','Latex')

    % Axis
    axis square
    set(gca,'FontSize',50,'XTick',[0,L1/2,L1],'YTick',[0,12000],'TickLabelInterpreter','latex')
    ax = gca;
    ax.XTickLabel{1} = "0"; ax.XTickLabel{2} = "$r$"; ax.XTickLabel{3} = "$L$"; % axis ticks
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    set(gcf,'Position',[1,49,1100,878])
    set(gcf,'Renderer','Painters')
end

% Spatial average
h = r(2) - r(1); % spacing
avgM0 = h*d/L1^d*trapz(r.^(d-1).*M0(1,:)); % trapezoidal rule

% Validate lambda
lambda_a = d*(d+2)*D/(L1^2 + beta1*L1*(d+2));
lambda_n = 1/avgM0;
fprintf('One-term absolute error: ')
disp(abs(lambda_a - lambda_n))

%% First moment
fig_part = {'a','b','c'};
for d = 1:3
    D = P*delta^2/(2*d*tau); % diffusivity
    M1 = bvp4c(@(r,M) first_ode(d,D,L1,beta1,r,M),@bound_moments,Minit);
    M1 = deval(M1,r);
    M1c = ((L1^2 + 2*beta1*L1)*((d+4)*(L1^2 + 2*beta1*L1) - 2*(d+2)*r.^2) + ...
        d*r.^4 + 4*d*beta1^2*L1^2)/(8*d^2*(d+2)*D^2);

    figure
    hold on
    box on
    plot(r,M1(1,:),'LineWidth',6.5,'Color','#05386B')
    plot(r,M1c,'LineStyle','--','LineWidth',6.5,'Color','#FF764A')
    y = ylabel('$\vspace{1cm}M_1(r)$','FontSize',20,'Interpreter','latex');
    xlim([0,L1]), ylim([0,1.1*max(M1(1,:))])
    if d == 1
        legend('Numerical','Analytical','Interpreter','Latex')
    end

    % Axis
    axis square
    format bank
    set(gca,'FontSize',50,'XTick',[0,L1/2,L1],'YTick',[0,1.1*max(M1(1,:))],...
        'TickLabelInterpreter','latex')
    %format default
    ax = gca;
    ax.XTickLabel{1} = "0"; ax.XTickLabel{2} = "$r$"; ax.XTickLabel{3} = "$L$"; % axis ticks
    set(y, 'Units', 'Normalized', 'Position', [-0.05, 0.5, 0]);
    set(gcf,'Position',[1,49,1100,878])
    set(gcf,'Renderer','Painters')
    text(0.05*L1,0.13*max(M1(1,:)),strcat('$\rm (',fig_part(d),')$'),'FontSize',75,'Interpreter','latex')
end

% Spatial average
avgM1 = h*d/L1^d*trapz(r.^(d-1).*M1(1,:));

lambda12_a = [d*(d+2)*D/(L1^2*(1 + sqrt(d/(d+4))) + beta1*L1*(d+2));
    d*(d+2)*D/(L1^2*(1 - sqrt(d/(d+4))) + beta1*L1*(d+2))];
lambda12_n = [1/(avgM0 + sqrt(avgM1 - avgM0^2));
    1/(avgM0 - sqrt(avgM1 - avgM0^2))];
fprintf('Two-term absolute error: ')
disp(abs(lambda12_a - lambda12_n))

%% Functions
% Zeroth moment ODE
function dMdr = zeroth_ode(d,D,r,M)
    dMdr = [M(2); -(1/D + (d-1)./r*M(2))];
end

% First moment ODE
function dMdr = first_ode(d,D,L1,beta1,r,M)
    M0 = @(r) (L1^2 - r.^2 + 2*beta1*L1)/(2*d*D);
    dMdr = [M(2); -(M0(r)/D + (d-1)./r*M(2))];
end

% Boundary conditions
function BC = bound_moments(Ma,Mb)
    BC = [Ma(2); Mb(1) + 2*Mb(2)]; % Ma = M_i(0) = 0 and Mb = M_i(L) = a
end