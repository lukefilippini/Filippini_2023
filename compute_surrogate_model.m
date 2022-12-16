function [c_exp,lambda] = compute_surrogate_model(d,D,IB,OB,model,IF,tc)
% COMPUTE_LAMBDA computes the rate parameter of a surrogate model for
% diffusion-controlled transport in a d-dimensional radially-symmetric 
% geometry with inner and boundary parameters specified in IB and OB, 
% respectively.
% Inputs:
%   d - number of dimensions (scalar)
%   D - mass diffusivity (scalar or vector)
%   IB - struct with fields L0, a0, b0 pertaining to the inner
%   boundary and corresponding boundary condition constants (struct)
%   OB - struct with fields L1, a1, b1 pertaining to the outer
%   boundary and corresponding boundary condition coefficients (struct)
%   model - surrogate model type (string)
%       Types: 'single', 'two' or 'weight'
%   IF - interface point (scalar)
%   tc - time points for continuum model (vector)
% Outputs:
%   c_exp - surrogate model (vector)
%   lambda - surrogate model rate parameter (scalar)

if nargin == 6
    tc = 0;
end

% Setup
L0 = IB.L0; % inner boundary
a1 = OB.a1; b1 = OB.b1; % outer boundary coefficients
beta = b1/a1; % measure of flux term dominance
innerBound = IB.innerBound; % inner boundary condition
outerBound = OB.outerBound; % outer boundary condition

% Compute rate parameter(s)
if isscalar(D) % rate parameters for homogeneous media
    L1 = OB.L1; % outer boundary
    if L0 == 0 % rate parameters for disc/sphere
        if isequal(model,'single')
            lambda = d*(d+2)*D/(L1^2 + (d+2)*beta*L1);
        elseif isequal(model,'two')
            lambda = zeros(2,1); % two parameters
            lambda(1) = d*(d+2)*D/(L1^2 + (d+2)*beta*L1 + L1^2*sqrt(d/(d+4)));
            lambda(2) = d*(d+2)*D/(L1^2 + (d+2)*beta*L1 - L1^2*sqrt(d/(d+4)));
        else
            lambda = zeros(2,1); % two parameters
            tau = (d+4)/(d*L1^4)*(((6-d)*L1^2 + (d+2)*(d+6)*beta*L1)/(d+6))^2;
            theta = 1/2 + 1/2*sqrt(tau/(tau+4)); % weighting parameter
            lambda(1) = d*(d+2)*D/(L1^2 + (d+2)*beta*L1 + L1^2*sqrt((1-theta)/theta*d/(d+4)));
            lambda(2) = d*(d+2)*D/(L1^2 + (d+2)*beta*L1 - L1^2*sqrt(theta/(1-theta)*d/(d+4)));
        end
    else % rate parameters for annulus/spherical shell
        if isequal(model,'single') % single parameter exponential model
            if isequal(innerBound,'reflect') % absorbing or semi-absorbing outer boundary
                if d == 2
                    lambda = 8*D*(L1^2 - L0^2)/(L0^4 + 4*L0^4*log(L1/L0) - L1^4 + ...
                        4*beta/L1*(L1^2 - L0^2)^2);
                else
                    lambda = 15*D*(L1^3 - L0^3)/(L1^5 - 5*L0^3*L1^2 + 9*L0^5 - ...
                        5*L0^6/L1 + 10*beta/L1^2*(L1^3 - L0^3)^2);
                end
            elseif isequal(innerBound,'absorb') && isequal(outerBound,'absorb') % two absorbing boundaries
                if d == 2
                    lambda = 8*D*log(L1/L0)/((L0^2 + L1^2)*log(L1/L0) - (L1^2 - L0^2));
                else
                    lambda = 60*D*(L1^3 - L0^3)/((L1 - L0)^3 * (4*L0^2 + 7*L0*L1 + 4*L1^2));
                end
            else % two semi-absorbing boundaries
            end
        elseif isequal(model,'two') % two parameter exponential model
            lambda = zeros(2,1); % two parameters
            if isequal(innerBound,'reflect') % absorbing or semi-absorbing outer boundary
                if d == 2
                    lambda(1) = 8*D*(L1^2 - L0^2)/(L1^4 - 4*L0^2*L1^2 + L0^4*(4*log(L1/L0) + 3) + ...
                        4*beta/L1*(L1^2 - L0^2)^2 + 1/sqrt(3)*sqrt(L1^8 - 10*L0^2*L1^6 + 7*L0^8 - ...
                        2*L0^6*L1^2*(24*log(L1/L0)^2 + 12*log(L1/L0) + 11) + 24*L0^4*L1^4*(log(L1/L0) + 1)));
                    lambda(1) = 8*D*(L1^2 - L0^2)/(L1^4 - 4*L0^2*L1^2 + L0^4*(4*log(L1/L0) + 3) + ...
                        4*beta/L1*(L1^2 - L0^2)^2 - 1/sqrt(3)*sqrt(L1^8 - 10*L0^2*L1^6 + 7*L0^8 - ...
                        2*L0^6*L1^2*(24*log(L1/L0)^2 + 12*log(L1/L0) + 11) + 24*L0^4*L1^4*(log(L1/L0) + 1)));
                else
                    lambda(1) = 15*D*(L1^3 - L0^3)/((L1 - L0)^3*(L1^2 + 3*L0*L1 + 6*L0^2 + 5*L0^3/L1) + ...
                        5*beta/L1^2*(L1^3 - L0^3)^2 + (L1 - L0)^3*sqrt(3/7)*sqrt(L1^4 + 6*L0*L1^3 + ...
                        21*L0^2*L1^2 + 41*L0^3*L1 + 36*L0^4));
                    lambda(2) = 15*D*(L1^3 - L0^3)/((L1 - L0)^3*(L1^2 + 3*L0*L1 + 6*L0^2 + 5*L0^3/L1) + ...
                        5*beta/L1^2*(L1^3 - L0^3)^2 - (L1 - L0)^3*sqrt(3/7)*sqrt(L1^4 + 6*L0*L1^3 + ...
                        21*L0^2*L1^2 + 41*L0^3*L1 + 36*L0^4));
                end
            elseif isequal(innerBound,'absorb') && isequal(outerBound,'absorb') % two absorbing boundaries
                if d == 2
                    lambda(1) = 8*D*log(L1/L0)/((L0^2 + L1^2)*log(L1/L0) - (L1^2 - L0^2) + ...
                        1/sqrt(3)*sqrt((L1^2 - L0^2)*(3*(L1^2 - L0^2)*(log(L1/L0)+1) - (L0^2 + L1^2)*log(L1/L0)^2)));
                    lambda(2) = 8*D*log(L1/L0)/((L0^2 + L1^2)*log(L1/L0) - (L1^2 - L0^2) - ...
                        1/sqrt(3)*sqrt((L1^2 - L0^2)*(3*(L1^2 - L0^2)*(log(L1/L0)+1) - (L0^2 + L1^2)*log(L1/L0)^2)));
                else
                    lambda(1) = 15*D*(L1^3 - L0^3)/((L1 - L0)^3*((L0 + L1)^2 + sqrt(3/4)*sqrt(16*L1^4 + ...
                        26*L0*L1^3 + 21*L0^2*L1^2 + 26*L0^3*L1 + 16*L0^4)));
                    lambda(1) = 15*D*(L1^3 - L0^3)/((L1 - L0)^3*((L0 + L1)^2 - sqrt(3/4)*sqrt(16*L1^4 + ...
                        26*L0*L1^3 + 21*L0^2*L1^2 + 26*L0^3*L1 + 16*L0^4)));
                end
            else % two semi-absorbing boundaries
            end
        else % three parameter exponential model
        end
    end
else % rate parameters for heterogeneous media with two layers (disc/sphere only)
    L1 = IF; % interface point
    L2 = OB.L2; % outer boundary
    if isequal(model,'single') % single exponential model
        lambda = d*(d+2)*D(1)*D(2)/((L2^2 + (d+2)*beta*L2)*D(1) + L1^(d+2)/L2^d*(D(2)-D(1)));
    else % two exponential model
        lambda = zeros(2,1); % two parameters
        lambda(1) = d*(d+2)*D(1)*D(2)/((L2^2 + (d+2)*beta*L2)*D(1) + L1^(d+2)/L2^d*(D(2)-D(1)) + ...
            1/L2^d + sqrt((d*D(1)^2*L2^(2*d+4) + d*(d+4)*(D(2)-D(1))*D(1)*L2^(d+2)*L1^(d+2) + ...
            (d+2)*((d+2)*D(1)^2 - (d+4)*D(1)*D(2) + 2*D(2)^2)*L2^d*L1^(d+4) - ...
            (d+4)*(D(1)-D(2))^2*L1^(2*d+4))/(d+4)));
        lambda(2) = d*(d+2)*D(1)*D(2)/((L2^2 + (d+2)*beta*L2)*D(1) + L1^(d+2)/L2^d*(D(2)-D(1)) + ...
            1/L2^d - sqrt((d*D(1)^2*L2^(2*d+4) + d*(d+4)*(D(2)-D(1))*D(1)*L2^(d+2)*L1^(d+2) + ...
            (d+2)*((d+2)*D(1)^2 - (d+4)*D(1)*D(2) + 2*D(2)^2)*L2^d*L1^(d+4) - ...
            (d+4)*(D(1)-D(2))^2*L1^(2*d+4))/(d+4)));
    end
end

% Surrogate model
if isequal(model,'single')
    c_exp = exp(-lambda*tc);
elseif isequal(model,'two')
    c_exp = (exp(-lambda(1)*tc) + exp(-lambda(2)*tc))/2;
else
    c_exp = theta*exp(-lambda(1)*tc) + (1-theta)*exp(-lambda(2)*tc);
end