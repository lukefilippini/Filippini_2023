function [P,IB,OB,IF] = get_case_parameters(Case,d,delta)
% GET_CASE_PARAMETERS returns the parameters associated with case 'A', 'B',
% 'C', 'D', 'E', 'F' or 'G'.
% Inputs:
%   Case - test case (string)
%   d - no. of dimensions (scalar)
%   delta - step size (scalar)
% Outputs:
%   P - probability of movement (scalar or array)
%   IB - parameters associated with the inner boundary: L0, a0, b0 (vector)
%   OB - parameters associated with the outer boundary: L1, a1, b1 (vector)
%   IF - interface point (scalar, applicable to two layer cases)

% Set interface to zero for homogeneous domains
IF = 0;

% Test case parameters (P0, P1 are set to zero for non-semi-absorbing
% conditions)
if isequal(Case,'A')
    P = 1; % probability of movement
    IB = struct('L0',0,'a0',0,'b0',1,'P0',0,'P1',0,'innerBound','reflect'); % inner bound parameters (disc/sphere must have a reflecting inner "boundary")
    OB = struct('L1',100,'a1',1,'b1',0,'P0',0,'P1',0,'outerBound','absorb'); % outer bound parameters (absorbing outer boundary)
elseif isequal(Case,'B')
    P = 1; % probability of movement
    P1 = 0.5; % probability of particle absorption at outer boundary
    IB = struct('L0',0,'a0',0,'b0',1,'innerBound','reflect'); % inner bound parameters (disc/sphere must have a reflecting inner "boundary")

    if d == 2
        b1 = delta / (P1 + pi*delta); % contribution of outer boundary flux term for disc
        OB = struct('L1',100,'a1',1,'b1',b1,'P0',0,'P1',P1,'outerBound','semi-absorb'); % outer bound parameters (semi-absorbing outer boundary with P1 = 0.5)
    else
        b1 = delta / (P1 + 2*pi*delta^2); % contribution of outer boundary flux term for sphere
        OB = struct('L1',100,'a1',1,'b1',b1,'P0',0,'P1',P1,'outerBound','semi-absorb'); % outer bound parameters (semi-absorbing outer boundary with P1 = 0.5)
    end
elseif isequal(Case,'C')
    P = 1; % probability of movement
    IB = struct('L0',50,'a0',0,'b0',1,'P0',0,'P1',0,'innerBound','reflect'); % inner bound parameters (reflecting inner boundary)
    OB = struct('L1',100,'a1',1,'b1',0,'P0',0,'P1',0,'outerBound','absorb'); % outer bound parameters (absorbing outer boundary)
elseif isequal(Case,'D')
    P = 1; % probability of movement
    P1 = 0.5; % probability of particle absorption at outer boundary
    IB = struct('L0',50,'a0',0,'b0',1,'innerBound','reflect'); % lower bound parameters (reflecting inner boundary)

    if d == 2
        b1 = delta / (P1 + pi*delta); % contribution of outer boundary flux term for disc
        OB = struct('L1',100,'a1',1,'b1',b1,'P0',0,'P1',P1,'outerBound','semi-absorb'); % outer bound parameters (semi-absorbing outer boundary with P1 = 0.5)
    else
        b1 = delta / (P1 + 2*pi*delta^2); % contribution of outer boundary flux term for sphere
        OB = struct('L1',100,'a1',1,'b1',b1,'P0',0,'P1',P1,'outerBound','semi-absorb'); % outer bound parameters (semi-absorbing outer boundary with P1 = 0.5)
    end
elseif isequal(Case,'E')
    P = 1; % probability of movement
    IB = struct('L0',50,'a0',1,'b0',0,'P0',0,'P1',0,'innerBound','absorb'); % inner bound parameters (absorbing inner boundary)
    OB = struct('L1',100,'a1',1,'b1',0,'P0',0,'P1',0,'outerBound','absorb'); % outer bound parameters (absorbing outer boundary)
elseif isequal(Case,'F')
    P = 1; % probability of movement
    P0 = 0.2; % probability of particle absorption at inner boundary
    P1 = 0.5; % probability of particle absorption at outer boundary
    
    if d == 2
        b0 = delta / (P0 + pi*delta); % contribution of outer boundary flux term for disc
        IB = struct('L0',0,'a0',1,'b0',b0,'P0',P0,'P1',0,'innerBound','semi-absorb'); % outer bound parameters (semi-absorbing outer boundary with P1 = 0.5)
    else
        b0 = delta / (P0 + 2*pi*delta^2); % contribution of outer boundary flux term for sphere
        IB = struct('L0',0,'a0',1,'b0',b0,'P0',P0,'P1',0,'innerBound','semi-absorb'); % outer bound parameters (semi-absorbing outer boundary with P1 = 0.5)
    end

    if d == 2
        b1 = delta / (P1 + pi*delta); % contribution of outer boundary flux term for disc
        OB = struct('L1',100,'a1',1,'b1',b1,'P0',0,'P1',P1,'outerBound','semi-absorb'); % outer bound parameters (semi-absorbing outer boundary with P1 = 0.5)
    else
        b1 = delta / (P1 + 2*pi*delta^2); % contribution of outer boundary flux term for sphere
        OB = struct('L1',100,'a1',1,'b1',b1,'P0',0,'P1',P1,'outerBound','semi-absorb'); % outer bound parameters (semi-absorbing outer boundary with P1 = 0.5)
    end
elseif isequal(Case,'G')
    P = [0.4,1]; % probability of movement in each layer
    IB = struct('L0',0,'a0',0,'b0',1,'P0',0,'P1',0,'innerBound','reflect'); % inner bound parameters (reflecting inner boundary)
    IF = 40; % interface point
    OB = struct('L2',100,'a1',1,'b1',0,'P0',0,'P1',0,'outerBound','absorb'); % outer bound parameters (absorbing outer boundary)
elseif isequal(Case,'H')
    P = [0.4,1]; % probability of movement in each layer
    P1 = 0.5; % probability of particle absorption at the outer boundary
    IB = struct('L0',0,'a0',0,'b0',1,'innerBound','reflect'); % inner bound parameters (reflecting inner boundary)
    IF = 40; % interface point

    if d == 2
        b1 = delta / (P1 + pi*delta); % contribution of outer boundary flux term for disc
        OB = struct('L1',100,'a1',1,'b1',b1,'P0',0,'P1',P1,'outerBound','semi-absorb'); % outer bound parameters (semi-absorbing outer boundary with P1 = 0.5)
    else
        b1 = delta / (P1 + 2*pi*delta^2); % contribution of outer boundary flux term for sphere
        OB = struct('L1',100,'a1',1,'b1',b1,'P0',0,'P1',P1,'outerBound','semi-absorb'); % outer bound parameters (semi-absorbing outer boundary with P1 = 0.5)
    end
end