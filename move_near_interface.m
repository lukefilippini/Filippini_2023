function [particleMoves,shiftDir,thetaChoice,phiChoice] = move_near_interface(d,pos,P,L1,delta,parts,thetaAll,phiAll)
% MOVE_NEAR_INTERFACE moves the particles near the interface (i.e. within a
% step size of delta of the interface boundary).
% Inputs:
%   pos - coordinates of each particle at current time step (N x 1 x d)
%   not (vector)
%   P - probaility of movement associated with each layer (vector)
%   L1 - interface boundary (scalar)
%   delta - step size (scalar)
%   parts - no. of partitions for the circle (or sphere) of radius delta (scalar)
%   thetaAll - angular choices for given partitions (vector)
%   phiAll - angular choices for given partitions (d = 3 only) (vector)
% Output:
%   particleMoves - indicator of particle movement (boolean)
%   shiftDir - shift direction of particle (d = 1) (scalar)
%   thetaChoice - choice of angle (d = 2,3) (scalar)
%   phiChoice - choice of angle (d = 3) (scalar)

if nargin == 6
    thetaAll = 0; phiAll = 0;
    thetaChoice = 0; phiChoice = 0;
elseif nargin == 7
    shiftDir = 0;
    phiAll = 0; phiChoice = 0;
else
    shiftDir = 0;
end

% Setup
meanProb = 0; % initialize mean probability
particleMoves = false; % initially assume particle doesn't move

% Calculate probaility of taking a step for each partition
if d == 1 % one dimension (line)
    randShiftProb = rand;
    shiftDir = 0; shift = [-1 1];
    for k = 1:2
        posHalf = pos + delta/2*shift(k);

        if posHalf < L1
            meanProb = meanProb + P(1)/2;
        else
            meanProb = meanProb + P(2)/2;
        end

        if randShiftProb < meanProb
            shiftDir = shift(k);
            particleMoves = true;
            break
        end
    end
elseif d == 2 % two dimensions (disc, annulus)
    randThetaProb = rand; % initialize probability for assessment against theta choice
    thetaChoice = 0; % initialize angle choice
    for k = 1:parts
        xPosHalf = pos(1,1) + delta/2*cos(thetaAll(k)); % x-pos for half step
        yPosHalf = pos(1,2) + delta/2*sin(thetaAll(k)); % y-pos for half step
        dist = sqrt(xPosHalf.^2 + yPosHalf.^2); % distance from origin
    
        if dist < L1
            meanProb = meanProb + P(1)/parts; % within L0 < x < L1
        else
            meanProb = meanProb + P(2)/parts; % within L1 < x < L2
        end

        if randThetaProb < meanProb % if angle is chosen, then move the particle
            thetaChoice = thetaAll(k);
            particleMoves = true;
            break
        end
    end
else % three dimensions (sphere, spherical shell)
    randProb = rand; % initialize probability for assessment against theta choice
    thetaChoice = 0; phiChoice = 0; % initialize choices for angles
    for k = 1:parts
        for j = 1:parts
            xPosHalf = pos(1,1) + delta/2*cos(thetaAll(k))*sin(phiAll(j)); % x-pos for half step
            yPosHalf = pos(1,2) + delta/2*sin(thetaAll(k))*sin(phiAll(j)); % y-pos for half step
            zPosHalf = pos(1,3) + delta/2*cos(phiAll(j)); % z-pos for half step
            dist = sqrt(xPosHalf.^2 + yPosHalf.^2 + zPosHalf.^2); % distance from origin
        
            if dist < L1
                meanProb = meanProb + P(1)/parts^2; % within L0 < x < L1
            else
                meanProb = meanProb + P(2)/parts^2; % within L1 < x < L2
            end

            if randProb < meanProb % if angle is chosen, then move the particle
                thetaChoice = thetaAll(k);
                phiChoice = phiAll(j);
                particleMoves = true;
                break
            end
        end

        if particleMoves
            break
        end
    end
end

end


