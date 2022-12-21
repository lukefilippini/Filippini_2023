function [c_exp,lambda,theta] = compute_surrogate_model(d,D,IB,OB,model,IF,tc)
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
%   theta - weighting parameter (scalar)

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
                    lambda = 8*D*(L1^2 - L0^2)/(L1^4 - 4*L0^2*L1^2 + L0^4*(4*log(L1/L0) + 3) + ...
                        4*beta/L1*(L1^2 - L0^2)^2);
                else
                    lambda = 15*D*(L1^3 - L0^3)/((L1^3 - L0^3)*(L1^2 + 5*L0^3/L1) - ...
                        9*L0^3*(L1^2 - L0^2) + 10*beta/L1^2*(L1^3 - L0^3)^2);
                end
            elseif isequal(innerBound,'absorb') && isequal(outerBound,'absorb') % two absorbing boundaries
                if d == 2
                    lambda = 8*D*log(L1/L0)/((L0^2 + L1^2)*log(L1/L0) - (L1^2 - L0^2));
                else
                    lambda = 60*D*(L1^3 - L0^3)/((L1 - L0)^3 * (4*L0^2 + 7*L0*L1 + 4*L1^2));
                end
            end
        else % two parameter or three-parameter exponential model
            lambda = zeros(2,1); % two parameters
            if isequal(model,'two')
                theta = 0.5;
            end

            if isequal(innerBound,'reflect') % absorbing or semi-absorbing outer boundary
                if d == 2
                    kappa1 = L0^4*L1^2*(L1^2 - L0^2); 
                    kappa2 = (L1^2 - L0^2)^3*(7*L0^2 - L1^2);

                    if isequal(model,'weight') % three-parameter model
                        kappa3 = L0^6*L1^2*(L1^2 - L0^2)*(L1*(5*L0^2 + L1^2) - 4*beta*(L1^2 - L0^2));
                        kappa4 = L0^8*L1^3*(L0^2 + L1^2);
                        kappa5 = L0^6*L1*(L1^2 - L0^2)^2*(12*L1^2 - 7*L0^2) - ...
                            24*beta*L0^4*L1^2*(L1^6 - 3*L0^2*L1^4 + 3*L0^4*L1^2 - L0^6);
                        kappa6 = (L1^2 - L0^2)^3*(L1*(145*L0^6 - 83*L0^4*L1^2 + 25*L0^2*L1^4 - 3*L1^6) - ...
                            24*beta*(L1^2 - L0^2)^2*(L1^2 - 7*L0^2));
                        tau = (288*kappa3*log(L1/L0)^2 - 1152*kappa4*log(L1/L0)^3 - ...
                            24*kappa5*log(L1/L0) - kappa6)^2/(12*L1^2*...
                            (24*kappa1*log(L1/L0) - 48*L0^6*L1^2*log(L1/L0)^2 - kappa2)^3);
                        theta = 1/2 + 1/2*sqrt(tau/(tau+4));
                    end

                    lambda(1) = 8*D*(L1^2 - L0^2)/(L1^4 - 4*L0^2*L1^2 + L0^4*(4*log(L1/L0) + 3) + ...
                        4*beta/L1*(L1^2 - L0^2)^2 + 1/sqrt(3)*sqrt((1-theta)/theta*(24*kappa1*log(L1/L0) - ...
                        48*L0^6*L1^2*log(L1/L0)^2 - kappa2)));
                    lambda(2) = 8*D*(L1^2 - L0^2)/(L1^4 - 4*L0^2*L1^2 + L0^4*(4*log(L1/L0) + 3) + ...
                        4*beta/L1*(L1^2 - L0^2)^2 - 1/sqrt(3)*sqrt(theta/(1-theta)*(24*kappa1*log(L1/L0) - ...
                        48*L0^6*L1^2*log(L1/L0)^2 - kappa2)));
                else
                    if isequal(model,'weight')
                        kappa1 = L1*(L1-L0)*(10*L0^2 + 4*L0*L1 + L1^2)*(54*L0^5 + 100*L0^4*L1 + 50*L0^3*L1^2 + ...
                            15*L0^2*L1^3 + 5*L0*L1^4 + L1^5);
                        kappa2 = (L0^2 + L0*L1 + L1^2)^2*(36*L0^4 + 41*L0^3*L1 + 21*L0^2*L1^2 + ...
                            6*L0*L1^3 + L1^4);
                        tau = 7*(kappa1 + 15*beta*kappa2)^2/(27*L1^4*(L1-L0)^2*(36*L0^4 + 41*L0^3*L1 + ...
                            21*L0^2*L1^2 + 6*L0*L1^3 + L1^4)^3);
                        theta = 1/2 + 1/2*sqrt(tau/(tau+4));
                    end
                    lambda(1) = 15*D*(L1^3 - L0^3)/((L1 - L0)^3*(L1^2 + 3*L0*L1 + 6*L0^2 + 5*L0^3/L1) + ...
                        5*beta/L1^2*(L1^3 - L0^3)^2 + (L1 - L0)^3*sqrt(3/7)*sqrt((1-theta)/theta*(L1^4 + 6*L0*L1^3 + ...
                        21*L0^2*L1^2 + 41*L0^3*L1 + 36*L0^4)));
                    lambda(2) = 15*D*(L1^3 - L0^3)/((L1 - L0)^3*(L1^2 + 3*L0*L1 + 6*L0^2 + 5*L0^3/L1) + ...
                        5*beta/L1^2*(L1^3 - L0^3)^2 - (L1 - L0)^3*sqrt(3/7)*sqrt(theta/(1-theta)*(L1^4 + 6*L0*L1^3 + ...
                        21*L0^2*L1^2 + 41*L0^3*L1 + 36*L0^4)));
                end
            elseif isequal(innerBound,'absorb') && isequal(outerBound,'absorb') % two absorbing boundaries
                if d == 2
                    if isequal(model,'weight')
                        kappa1 = (L1^2 - L0^2)^2*(19*L0^4 + 46*L0^2*L1^2 + 19*L1^4);
                        kappa2 = (L1^2 - L0^2)^3*(L0^2 + L1^2);
                        tau = (6*(L1^8 - L0^8)*log(L1/L0)^3 - kappa1*log(L1/L0)^2 + ...
                            18*kappa2*log(L1/L0))^2/(48*(L1^2 - L0^2)^5*((L1^2 - L0^2)*...
                            (log(L1/L0)^2 + 3) - 3*(L0^2 + L1^2)*log(L1/L0))^3);
                        theta = 1/2 + 1/2*sqrt(tau/(tau+4));
                    end

                    lambda(1) = 8*D*log(L1/L0)/((L0^2 + L1^2)*log(L1/L0) - (L1^2 - L0^2) + ...
                        (L1^2 - L0^2)/sqrt(3)*sqrt((1-theta)/theta*(log(L1/L0)^2 - ...
                        3*(L0^2 + L1^2)/(L1^2 - L0^2)*log(L1/L0) + 3)));
                    lambda(2) = 8*D*log(L1/L0)/((L0^2 + L1^2)*log(L1/L0) - (L1^2 - L0^2) - ...
                        (L1^2 - L0^2)/sqrt(3)*sqrt(theta/(1-theta)*(log(L1/L0)^2 - ...
                        3*(L0^2 + L1^2)/(L1^2 - L0^2)*log(L1/L0) + 3)));
                else
                    eta = 16*L0^4 + 26*L0^3*L1 + 21*L0^2*L1^2 + 26*L0*L1^3 + 16*L1^4;
                    if isequal(model,'weight')
                        tau = 7*(64*L0^6 + 471*L0^5*L1 + 780*L0^4*L1^2 + 745*L0^3*L1^3 + ...
                            780*L0^2*L1^4 + 471*L0*L1^5 + 64*L1^6)^2/(27*eta^3);
                        theta = 1/2 + 1/2*sqrt(tau/(tau+4));
                    end
                    lambda(1) = 60*D*(L1^3 - L0^3)/((L1 - L0)^3*(4*L0^2 + 7*L0*L1 + ...
                        4*L1^2 + sqrt(3/7)*sqrt((1-theta)/theta*eta)));
                    lambda(2) = 60*D*(L1^3 - L0^3)/((L1 - L0)^3*(4*L0^2 + 7*L0*L1 + ...
                        4*L1^2 - sqrt(3/7)*sqrt(theta/(1-theta)*eta)));
                end
            end
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
             1/L2^d*sqrt((d*D(1)^2*L2^(2*d+4) + d*(d+4)*(D(2)-D(1))*D(1)*L2^(d+2)*L1^(d+2) + ...
            (d+2)*((d+2)*D(1)^2 - (d+4)*D(1)*D(2) + 2*D(2)^2)*L2^d*L1^(d+4) - ...
            (d+4)*(D(1)-D(2))^2*L1^(2*d+4))/(d+4)));
        lambda(2) = d*(d+2)*D(1)*D(2)/((L2^2 + (d+2)*beta*L2)*D(1) + L1^(d+2)/L2^d*(D(2)-D(1)) - ...
             1/L2^d*sqrt((d*D(1)^2*L2^(2*d+4) + d*(d+4)*(D(2)-D(1))*D(1)*L2^(d+2)*L1^(d+2) + ...
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