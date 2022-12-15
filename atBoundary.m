function [notAbsorbed] = atBoundary(boundCond,Pa)
% ATBOUNDARY applies one of three possible boundary conditions to a single
% particle that has crossed a boundary in a radial random walk simulation.
% Inputs:
%   boundCond - boundary condition to be applied
%   Pa - probability of absorption for a semi-absorbing condition (scalar)
% Outputs:
%   notAbsorbed - particle is still in domain (scalar, 1 or 0)

% Setup
notAbsorbed = 0; % initialize output

if isequal(boundCond,'reflect') % reflecting boundary
    notAbsorbed = 1; % particle remains in domain
else
    if isequal(boundCond,'absorb') % guaranteed absorption
        particleIsAbsorbed = true;
    else % semi-absorbing boundary
        particleIsAbsorbed = Pa > rand;
    end

    if ~particleIsAbsorbed
        notAbsorbed = 1;
    end
end

end