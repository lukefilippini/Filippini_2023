function Ps = heterogeneous_stochastic_model(d,IB,OB,L1,Np,T,tau,delta,P,parts)
% HETEROGENEOUS_STOCHASTIC_MODEL simulates a stochastic random walk
% process for radially-symmetric geometries of dimension 'd' with two
% concentric layers of different compositions (diffusivities). The
% simulation is run for Np particles, for a total time of T, where each
% particle takes a step of size delta every tau seconds. 
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
if (d <= 0) || (d > 3)
    error('The number of dimensions must be 1, 2 or 3.')
end

NSteps = floor(T/tau); % no. of time steps
Ps = zeros(Np,NSteps+1); % particles in system
Ps(:,1) = 1; % particles are all initially within domain
pos = zeros(Np,d); % particle positions
L0 = IB.L0; L2 = OB.L2; % inner and outer boundaries
P0 = IB.P0; P1 = OB.P1; % probabilities of absorption (semi-absorbing condition)
innerBound = IB.innerBound; % inner boundary condition
outerBound = OB.outerBound; % outer boundary condition

% Initial distribution
r = (rand(Np,1)*(L2^d - L0^d) + L0^d).^(1/d); % initial radial position
if d > 1
    theta = 2*pi*rand; % initialize all angular positions
    thetaAll = (2*pi*((1:parts) - 1) ./ parts)'; % all possible theta choices for movement near interface
    if d == 3
        phi = acos(1 - 2*rand); % initialize phi
        phiAll = acos(1 - 2*((1:parts) - 1) / parts)'; % all possible phi choices for movement near interface
        pos(:,1) = r .* cos(theta) .* sin(phi); % initialize x-coordinates
        pos(:,2) = r .* sin(theta) .* sin(phi); % initialize y-coordinates
        pos(:,3) = r .* cos(phi); % initialize z-coordinates
    else
        pos(:,1) = r .* cos(theta); pos(:,2) = r .* sin(theta); % initialize particle positions
    end
else
    pos(:) = r; % initialize particle positions (d = 1)
end

%posStore = zeros(Np,NSteps+1,d);
%posStore(:,1,1) = pos(:,1);

% Heterogeneous random walk
for i = 1:Np
    for j = 1:NSteps
        % Calculate distance between center of circles
        if d > 1
            dist = sqrt(sum(pos(i,:).^2));
        else
            dist = pos(i);
        end
        
        % Determine if circles intersect
        intersect = (dist <= L1 + delta) && (dist >= L1 - delta);
      
        if intersect % Determine whether the particle near the interface should move
            if d == 1
                [particleMoves,shiftDir,~,] = move_near_interface(d,pos(i),P,L1,delta,parts);
            elseif d == 2
                [particleMoves,~,thetaChoice,~,] = move_near_interface(d,pos(i,:),P,L1,delta,parts,thetaAll);
            else
                [particleMoves,~,thetaChoice,phiChoice] = move_near_interface(d,pos(i,:),P,L1,delta,parts,thetaAll,phiAll);
            end
            
            if particleMoves % move the particle
                if d == 1
                    pos(i) = pos(i) + delta*shiftDir; % new x-coordinate
                elseif d == 2
                    pos(i,1) = pos(i,1) + delta*cos(thetaChoice); % new x-coordinate
                    pos(i,2) = pos(i,2) + delta*sin(thetaChoice); % new y-coordinate
                else
                    pos(i,1) = pos(i,1) + delta*cos(thetaChoice)*sin(phiChoice); % new x-coordinate
                    pos(i,2) = pos(i,2) + delta*sin(thetaChoice)*sin(phiChoice); % new y-coordinate
                    pos(i,3) = pos(i,3) + delta*cos(phiChoice); % new z-coordinate
                end
            end

            Ps(i,j+1) = 1;
        else % Determine whether the particle away from the interface should move
            if d > 1 
                r = sqrt(sum(pos(i,:).^2)); % distance from origin
            else
                r = pos(i);
            end

            if r < L1
                particleMoves = P(1) > rand; % layer L0 < x <= l2
            else
                particleMoves = P(2) > rand; % layer l2 < x < L2
            end

            if particleMoves % shift the particle by a step of size delta
                if d > 1
                    theta = 2*pi*rand;
                    if d == 2
                        pos(i,1) = pos(i,1) + delta*cos(theta); % new x-coordinate
                        pos(i,2) = pos(i,2) + delta*sin(theta); % new y-coordinate
                    else
                        phi = acos(1 - 2*rand);
                        pos(i,1) = pos(i,1) + delta*cos(theta)*sin(phi); % new x-coordinate
                        pos(i,2) = pos(i,2) + delta*sin(theta)*sin(phi); % new y-coordinate
                        pos(i,3) = pos(i,3) + delta*cos(phi); % new z-coordinate
                    end
                else
                    shift = sign(rand - 0.5);
                    pos(i) = pos(i) + delta*shift; % new x-coordinate
                end
            end

            % Calculate new radial position
            if d > 1
                r = sqrt(sum(pos(i,:).^2));
            else
                r = pos(i);
            end
        
            % Boundary conditions
            if r >= L2 % particle has crossed outer boundary
                Ps(i,j+1) = at_boundary(outerBound,P1); % determine whether the particle is absorbed or not
                
                if isequal(outerBound,'reflect')
                    if d == 1
                        pos(i) = pos(i) - delta; % shift position back to the left
                    elseif d == 2
                        pos(i,1) = pos(i,1) - delta*cos(theta); % keep the same x-coordinate
                        pos(i,2) = pos(i,2) - delta*sin(theta); % keep the same y-coordinate
                    else
                        pos(i,1) = pos(i,1) - delta*cos(theta)*sin(phi); % keep the same x-coordinate
                        pos(i,2) = pos(i,2) - delta*sin(theta)*sin(phi); % keep the same y-coordinate
                        pos(i,3) = pos(i,3) - delta*cos(phi); % keep the same z-coordinate
                    end
                end
            elseif r <= L0 % particle has crossed inner boundary
                Ps(i,j+1) = at_boundary(innerBound,P0); % determine whether the particle is absorbed or not
    
                if isequal(innerBound,'reflect')
                    if d == 1
                        pos(i) = pos(i) + delta; % shif the position back to the right
                    elseif d == 2
                        pos(i,1) = pos(i,1) - delta*cos(theta); % keep the same x-coordinate
                        pos(i,2) = pos(i,2) - delta*sin(theta); % keep the same y-coordinate
                    else
                        pos(i,1) = pos(i,1) - delta*cos(theta)*sin(phi); % keep the same x-coordiante
                        pos(i,2) = pos(i,2) - delta*sin(theta)*sin(phi); % keep the same y-coordiante
                        pos(i,3) = pos(i,3) - delta*cos(phi); % keep the same z-coordinate
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
end