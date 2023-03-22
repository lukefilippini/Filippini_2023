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
theta = 0; % initialize theta
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
            omega = (d+4)/(d*L1^4)*(((6-d)*L1^2 + (d+2)*(d+6)*beta*L1)/(d+6))^2;
            theta = 1/2 + 1/2*sqrt(omega/(omega+4)); % weighting parameter
            lambda(1) = d*(d+2)*D/(L1^2 + (d+2)*beta*L1 + L1^2*sqrt((1-theta)/theta*d/(d+4)));
            lambda(2) = d*(d+2)*D/(L1^2 + (d+2)*beta*L1 - L1^2*sqrt(theta/(1-theta)*d/(d+4)));
        end
    else % rate parameters for annulus/spherical shell
        if isequal(model,'single') % single parameter exponential model
            if isequal(innerBound,'reflect') % absorbing or semi-absorbing outer boundary
                lambda = d*(d+2)*D*(L1^d - L0^d)/(L1^(d+2) - (d+2)*L0^d*(L1^2-L0^2) + ...
                    (d+2)*L0^(2*d)*integral(@(r) 1./r.^(d-1),L0,L1) - L0^(d+2) + ...
                    (d+2)*beta*(L1^d-L0^d)^2/L1^(d-1));
            elseif isequal(innerBound,'absorb') && isequal(outerBound,'absorb') % two absorbing boundaries
                lambda = 4*d*(d+2)*D*integral(@(r) 1./r.^(d-1),L0,L1)*...
                    (L1^d - L0^d)/(4*(L1^(d+2) - L0^(d+2))*integral(@(r) 1./r.^(d-1),L0,L1) - ...
                    (d+2)*(L1^2-L0^2)^2);
            end
        else % two parameter or three-parameter exponential model
            lambda = zeros(2,1); % two parameters
            if isequal(model,'two')
                theta = 0.5;
            end

            if isequal(innerBound,'reflect') % absorbing or semi-absorbing outer boundary
                if d == 1
                    if isequal(model,'weight')
                        omega = 5*(5*(L1-L0) + 21*beta)^2/(49*(L1-L0)^2);
                        theta = 1/2 + 1/2*sqrt(omega/(omega+4));
                    end
                    lambda(1) = 3*D/((L1-L0)*(L1-L0+3*beta) + (L1-L0)^2*sqrt((1-theta)/(5*theta)));
                    lambda(2) = 3*D/((L1-L0)*(L1-L0+3*beta) - (L1-L0)^2*sqrt(theta/(5*(1-theta))));
                elseif d == 2
                    kappa21 = -(L1^2 - L0^2)^3*(7*L0^2 - L1^2) - ...
                        24*L0^4*L1^2*log(L1/L0)*(2*L0^2*log(L1/L0) + L0^2 - L1^2);

                    if isequal(model,'weight') % three-parameter model
                        kappa22 = L0^6*L1^2*(L1^2-L0^2)*(L1*(5*L0^2+L1^2)-4*beta*(L1^2-L0^2));
                        kappa23 = L0^4*L1*(L1^2-L0^2)^2*(7*L0^4-12*L0^2*L1^2+24*beta*L1*(L1^2-L0^2));
                        kappa24 = (L1^2-L0^2)^3*(L1*(3*L1^6-25*L0^2*L1^4+83*L0^4*L1^2-145*L0^6) - ...
                            24*beta*(7*L0^2-L1^2)*(L1^2-L0^2)^2);
                        omega = (288*kappa22*log(L1/L0)^2 - 1152*L0^8*L1^3*(L0^2+L1^2)*log(L1/L0)^3 + ...
                            24*kappa23*log(L1/L0) + kappa24)^2/(12*L1^2*kappa21^3);
                        theta = 1/2 + 1/2*sqrt(omega/(omega+4));
                    end

                    lambda(1) = 8*D*(L1^2 - L0^2)/((L1^2 - L0^2)*(L1^2 - 3*L0^2) + 4*L0^4*log(L1/L0) + ...
                        4*beta/L1*(L1^2 - L0^2)^2 + 1/sqrt(3)*sqrt((1-theta)/theta*kappa21));
                    lambda(2) = 8*D*(L1^2 - L0^2)/((L1^2 - L0^2)*(L1^2 - 3*L0^2) + 4*L0^4*log(L1/L0) + ...
                        4*beta/L1*(L1^2 - L0^2)^2 - 1/sqrt(3)*sqrt(theta/(1-theta)*kappa21));
                else
                    kappa31 = L1^4 + 6*L0*L1^3 + 21*L0^2*L1^2 + 41*L0^3*L1 + 36*L0^4;
                    if isequal(model,'weight')
                        kappa32 = L1^5 + 5*L0*L1^4 + 15*L0^2*L1^3 + 50*L0^3*L1^2 + 100*L0^4*L1 + 54*L0^5;
                        omega = 7*(L1*(L1-L0)^3*(L1^2+4*L0*L1+10*L0^2)*kappa32 + ...
                            15*beta*kappa31*(L1^3-L0^3)^2)^2/(27*L1^4*(L1-L0)^6*kappa31^3);
                        theta = 1/2 + 1/2*sqrt(omega/(omega+4));
                    end
                    lambda(1) = 15*D*(L1^3 - L0^3)/((L1 - L0)^3*(L1^2 + 3*L0*L1 + 6*L0^2 + 5*L0^3/L1) + ...
                        5*beta/L1^2*(L1^3 - L0^3)^2 + (L1 - L0)^3*sqrt(3/7)*sqrt((1-theta)/theta*kappa31));
                    lambda(2) = 15*D*(L1^3 - L0^3)/((L1 - L0)^3*(L1^2 + 3*L0*L1 + 6*L0^2 + 5*L0^3/L1) + ...
                        5*beta/L1^2*(L1^3 - L0^3)^2 - (L1 - L0)^3*sqrt(3/7)*sqrt(theta/(1-theta)*kappa31));
                end
            elseif isequal(innerBound,'absorb') && isequal(outerBound,'absorb') % two absorbing boundaries
                if d == 1
                    if isequal(model,'weight')
                        theta = 1/2 + 1/2*sqrt(125/321);
                    end
                    lambda(1) = 12*D/((L1-L0)^2 * (1 + sqrt((1-theta)/(5*theta))));
                    lambda(2) = 12*D/((L1-L0)^2 * (1 - sqrt(theta/(5*(1-theta)))));
                elseif d == 2
                    xi21 = 3*(L1^2 - L0^2)^2 - 3*(L1^4 - L0^4)*log(L1/L0) + ...
                        (L1^2 - L0^2)^2*log(L1/L0)^2;
                    if isequal(model,'weight')
                        xi22 = (L1^2-L0^2)*(19*L0^4 + 46*L0^2*L1^2 + 19*L1^4);
                        omega = log(L1/L0)^2*(18*(L0^2+L1^2)*(L1^2-L0^2)^2 + ...
                            6*(L0^2+L1^2)*(L0^4+L1^4)*log(L1/L0)^2 - xi22*log(L1/L0))^2/(48*xi21^3);
                        theta = 1/2 + 1/2*sqrt(omega/(omega+4));
                    end

                    lambda(1) = 8*D*log(L1/L0)/((L0^2 + L1^2)*log(L1/L0) - (L1^2 - L0^2) + ...
                        1/sqrt(3)*sqrt((1-theta)/theta*xi21));
                    lambda(2) = 8*D*log(L1/L0)/((L0^2 + L1^2)*log(L1/L0) - (L1^2 - L0^2) - ...
                        1/sqrt(3)*sqrt(theta/(1-theta)*xi21));
                else
                    xi31 = 16*L0^4 + 26*L0^3*L1 + 21*L0^2*L1^2 + 26*L0*L1^3 + 16*L1^4;
                    if isequal(model,'weight')
                        omega = 7*(64*L0^6 + 471*L0^5*L1 + 780*L0^4*L1^2 + 745*L0^3*L1^3 + ...
                            780*L0^2*L1^4 + 471*L0*L1^5 + 64*L1^6)^2/(27*xi31^3);
                        theta = 1/2 + 1/2*sqrt(omega/(omega+4));
                    end
                    lambda(1) = 60*D*(L1^3 - L0^3)/((L1 - L0)^3*(4*L0^2 + 7*L0*L1 + ...
                        4*L1^2 + sqrt(3/7)*sqrt((1-theta)/theta*xi31)));
                    lambda(2) = 60*D*(L1^3 - L0^3)/((L1 - L0)^3*(4*L0^2 + 7*L0*L1 + ...
                        4*L1^2 - sqrt(3/7)*sqrt(theta/(1-theta)*xi31)));
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