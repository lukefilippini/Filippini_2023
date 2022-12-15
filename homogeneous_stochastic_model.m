function Ps = homogeneous_stochastic_model(d,IB,OB,Np,T,tau,delta,P)
% HOMOGENEOUS_STOCHASTIC_MODEL conducts a random walk simulation for a
% radially-symmetric object of d dimensions, defined over the radial domain
% [L0,L1]. For Np particles, each particle takes a step of size delta after 
% every tau "seconds" with probability P until the final time T is elapsed. 
% Inputs:
%   d - no. of dimensions (scalar)
%   IB - struct with fields L0, P0, innerBound pertaining to the inner
%   boundary location, condition ('absorb', 'reflect' or 'semi-absorb') and 
%   probability of absorption for a semi-absorbing condition (struct)
%   OB - struct with fields L1, P1, outerBound pertaining to the outer
%   boundary location, condition ('absorb', 'reflect' or 'semi-absorb') and 
%   probability of absorption for a semi-absorbing condition (struct)
%   Np - no. of non-interacting particles (scalar or vector)
%   T - total time (scalar)
%   tau - time step size (scalar)
%   delta - spatial step size (scalar)
%   P - probability of particle moving (scalar)
% Outputs:
%   Ps - non-absorbed particles at each time point (matrix)

% Setup
NSteps = floor(T/tau); % no. of time steps
t = 0:tau:T; % time points
Ps = zeros(Np,NSteps+1); % particles in system
Ps(:,1) = 1; % particles are all initially within domain
pos = zeros(Np,d); % particle positions
L0 = IB.L0; L1 = OB.L1; % inner and outer boundaries
P0 = IB.P0; P1 = OB.P1; % probabilities of absorption (semi-absorbing condition)
innerBound = IB.innerBound; % inner boundary condition
outerBound = OB.outerBound; % outer boundary condition


% Initial distribution
r = (rand(Np,1)*(L1^d - L0^d) + L0^d).^(1/d); % initial radial position
theta = 2*pi*rand; % initialize all angular positions
if d == 3
    phi = acos(1 - 2*rand); % initialize phi
    pos(:,1) = r .* cos(theta) .* sin(phi); % initial x-coordinates
    pos(:,2) = r .* sin(theta) .* sin(phi); % initial y-coordinates
    pos(:,3) = r .* cos(phi); % initial z-coordinates
else
    pos(:,1) = r .* cos(theta); pos(:,2) = r .* sin(theta); % initialize particle positions
end

% Random walk simulation 
for i = 1:Np % for each particle
    for j = 1:NSteps % for each time step
        particleMoves = P > rand; % particle moves (1) or stays still (0)
        if particleMoves
            theta = 2*pi*rand;
            if d == 2
                pos(i,1) = pos(i,1) + delta*cos(theta); % new x-coordinate (d = 2)
                pos(i,2) = pos(i,2) + delta*sin(theta); % new y-coordinate (d = 2)
            else
                phi = acos(1 - 2*rand);
                pos(i,1) = pos(i,1) + delta*cos(theta)*sin(phi); % new x-coordinate (d == 3)
                pos(i,2) = pos(i,2) + delta*sin(theta)*sin(phi); % new y-coordinate (d == 3)
                pos(i,3) = pos(i,3) + delta*cos(phi); % new z-coordinate (d == 3)
            end
        end

        % Calculate norm
        r = sqrt(sum(pos(i,:).^2)); % radial position
        
        % Boundary conditions
        if r >= L1 % particle has crossed outer boundary
            Ps(i,j+1) = atBoundary(outerBound,P1); % determine whether the particle is absorbed or not
            
            if isequal(outerBound,'reflect')
                if d == 2
                    pos(i,1) = pos(i,1) - delta*cos(theta); % keep original x-coordinate
                    pos(i,2) = pos(i,2) - delta*sin(theta); % keep original y-coordinate
                else
                    pos(i,1) = pos(i,1) - delta*cos(theta)*sin(phi); % keep original x-coordinate
                    pos(i,2) = pos(i,2) - delta*sin(theta)*sin(phi); % keep original y-coordinate
                    pos(i,3) = pos(i,3) - delta*cos(phi); % keep original z-coordinate
                end
            end
        elseif r <= L0 % particle has crossed inner boundary
            Ps(i,j+1) = atBoundary(innerBound,P0); % determine whether the particle is absorbed or not

            if isequal(innerBound,'reflect')
                if d == 1
                    pos(i) = pos(i) + delta; % shift particle back to the right 
                elseif d == 2
                    pos(i,1) = pos(i,1) - delta*cos(theta); % keep original x-coordinate
                    pos(i,2) = pos(i,2) - delta*sin(theta); % keep original y-coordinate
                else
                    pos(i,1) = pos(i,1) - delta*cos(theta)*sin(phi); % keep original x-coordinate
                    pos(i,2) = pos(i,2) - delta*sin(theta)*sin(phi); % keep original y-coordinate
                    pos(i,3) = pos(i,3) - delta*cos(phi); % keep original z-coordinate
                end
            end
        else
            Ps(i,j+1) = 1; % if particle has not been absorbed, it remains in motion
        end

        if Ps(i,j+1) == 0 % if particle has exited the domain, move to the next particle
            break
        end
    end
end