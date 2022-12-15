function propConfint = compute_min_max_intervals(Ns,d,IB,OB,Np,T,tau,delta,P,L1,parts)
% COMPUTE_MIN_MAX_INTERVALS computes minimum and maximum particle 
% proportions at each time step for a radial random walk simulated Ns 
% times with Np particles,
% Inputs:
%   d - no. of dimensions (scalar)
%   IB - struct with fields L0, P0, innerBound pertaining to the inner
%   boundary location, condition ('absorb', 'reflect' or 'semi-absorb') and 
%   probability of absorption for a semi-absorbing condition (struct)
%   OB - struct with fields L1 (or L2), P1, outerBound pertaining to the outer
%   boundary location, condition ('absorb', 'reflect' or 'semi-absorb') and 
%   probability of absorption for a semi-absorbing condition (struct)
%   Np - no. of non-interacting particles (scalar)
%   T - total time (scalar)
%   tau - time step size (scalar)
%   delta - spatial step size (scalar)
%   P - probability of particle moving (scalar)
%   L1 - interface point (scalar, optional)
%   parts - number of angular partitions (scalar, optional)
% Outputs:
%   propConfit - confidence intervals for prop. remaining at each time step
%   (2 x no. of time steps) (matrix)
%   propRecord - prop. of particles remaining at each time step for Ns
%   simulations (Ns x no. of time steps) (matrix)

if nargin == 9 % if interface and number of partitions are not included
    L1 = 0; parts = 0;
end

% Setup
t = 0:tau:T; % time steps
propRecord = zeros(Ns,length(t)); % storage for prop. of remaining particles (Ns simulations)
propConfint = zeros(2,length(t)); % storage for min-max intervals

% Radial random walk (Ns simulations)
parfor i = 1:Ns
    if isscalar(P)
        Ps = homogeneous_stochastic_model(d,IB,OB,Np,T,tau,delta,P); % random walk simulation
    else
        Ps = heterogeneous_stochastic_model(d,IB,OB,L1,Np,T,tau,delta,P,parts);
    end
    propRecord(i,:) = mean(Ps); % proportion of particles remaining at each time step
end


% 95% Confidence Intervals
propConfint(1,:) = min(propRecord);
propConfint(2,:) = max(propRecord);

end