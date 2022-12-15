function c = heterogeneous_continuum_model(d,D,IB,OB,x1,x2,Nt,T,c0)
% HETEROGENEOUS_CONTINUUM_MODEL computes the numerical solution for diffusion out of a
% d-dimensional, radially-symmetric object with two concentric layers using a 
% finite volume scheme with Crank Nicolson time-stepping. 
% Inputs:
%   d - number of dimensions (scalar)
%   D - mass diffusivity (scalar)
%   IB - struct with fields L0, a0, b0 (in order) pertaining to the inner
%   boundary and corresponding boundary condition constants (struct)
%   OB - struct with fields L1, a1, b1 (in order) pertaining to the outer
%   boundary and corresponding boundary condition coefficients (struct)
%   x1 - spatial nodes for L0 <= x <= L1 (vector)
%   x2 - spatial nodes for L1 <= x <= L2 (vector)
%   Nt - number of time steps (scalar)
%   T - end time (scalar)
%   c0 - initial condition (scalar or vector)
% Outputs:
%   c - numerical solution, where the columns correspond to the solution 
%       at each time point (matrix).

% Setup
dt = T/Nt; % time step size
Nr1 = length(x1); Nr2 = length(x2); Nr = Nr1 + Nr2 - 1; % no. of nodes
h1 = diff(x1); h2 = diff(x2); h = [h1 h2]; % node spacing
x = [x1 x2]; % spatial nodes
a0 = IB.a0; b0 = IB.b0; % inner boundary coefficients
a1 = OB.a1; b1 = OB.b1; % outer boundary coefficients


% Control volumes
V = zeros(Nr,1); % no. of control volumes
V(1) = h(1)/2; V(2:Nr-1) = (h(1:Nr-2) + h(2:Nr-1))/2; V(Nr) = h(Nr-1)/2; % control volumes

% West control volume boundaries
w = zeros(Nr,1); % no. of west boundaries
w(1) = x(1); w(2:Nr) = (x(1:Nr-1) + x(2:Nr))/2; % west boundaries

% East control volume boundaries
e = zeros(Nr,1); % no. of east boundaries
e(1:Nr-1) = w(2:Nr); e(Nr) = x(Nr); % east boundaries

% Construct matrix of node coefficients %
A = zeros(Nr,Nr);

% First row
A(1,1) = D(1)*(-e(1)^(d-1)/h(1) + x(1)^(d-1)*a0/b0) / (V(1)*x(1)^(d-1));
A(1,2) = D(1)*e(1)^(d-1)/(V(1)*x(1)^(d-1)*h(1));

% Inner layer (D = D1)
for i = 2:(Nr1-1)
    A(i,i-1) = D(1)*w(i)^(d-1) / (V(i)*x(i)^(d-1)*h(i-1));

    if (i == 2) && ((x(1) == 0) && (d > 1))
        A(i,i) = -D(1)*e(i)^(d-1)/(V(i)*x(i)^(d-1)*h(i));
    else
        A(i,i) = -D(1)*(w(i)^(d-1)/h(i-1) + e(i)^(d-1)/h(i)) / (V(i)*x(i)^(d-1));
    end 

    A(i,i+1) = D(1)*e(i)^(d-1) / (V(i)*x(i)^(d-1)*h(i));
end

% Interface point
A(Nr1,Nr1-1) = D(1)*w(Nr1)^(d-1) / (V(Nr1)*x(Nr1)^(d-1)*h(Nr1-1));
A(Nr1,Nr1) = -(D(1)*w(Nr1)^(d-1)/h(Nr1-1) + D(2)*e(Nr1)^(d-1)/h(Nr1)) / (V(Nr1)*x(Nr1)^(d-1));
A(Nr1,Nr1+1) = D(2)*e(Nr1)^(d-1) / (V(Nr1)*x(Nr1)^(d-1)*h(Nr1));

% Outer layer (D = D2)
for i = (Nr1+1):(Nr-1)
    A(i,i-1) = D(2)*w(i)^(d-1) / (V(i)*x(i)^(d-1)*h(i-1));
    A(i,i) = -D(2)*(w(i)^(d-1)/h(i-1) + e(i)^(d-1)/h(i)) / (V(i)*x(i)^(d-1));
    A(i,i+1) = D(2)*e(i)^(d-1) / (V(i)*x(i)^(d-1)*h(i));
end

% Final row
A(Nr,Nr-1) = D(2)*w(Nr)^(d-1) / (V(Nr)*x(Nr)^(d-1)*h(Nr-1));
A(Nr,Nr) = -D(2)*(w(Nr)^(d-1)/h(Nr-1) + x(Nr)^(d-1)*a1/b1) / (V(Nr)*x(Nr)^(d-1));

% Solve for c (Crank-Nicolson Time Stepping)
if (x(1) == 0) && (d > 1) % for a disc/sphere solve only for nodes 2,...,Nr (node 1 is equal to node 2)
    A = A(2:Nr,2:Nr);
    c = zeros(Nr-1,Nt+1); % solution matrix
    I = eye(Nr-1); % identity matrix
else % in all other cases, solve for all nodes
    c = zeros(Nr,Nt+1); % solution matrix
    I = eye(Nr); % identity matrix
end

c(:,1) = c0; % initial condition
Atilde = I - dt/2*A;

for n = 1:Nt % Solve at each time step
    btilde = (I + dt/2*A)*c(:,n);

    % Incorporate Dirichlet conditions if necessary
    if b0 == 0
        Atilde(1,:) = [1 zeros(1,Nr-1)];
        btilde(1) = 0;
    end

    if b1 == 0
        % Only solving for c2,...,cN if l0 = 0 and d > 1
        if (x(1) == 0) && (d > 1)
            Atilde(Nr-1,:) = [zeros(1,Nr-2) 1];
            btilde(Nr-1) = 0;
        else
            Atilde(Nr,:) = [zeros(1,Nr-1) 1];
            btilde(Nr) = 0;
        end
    end

    % Solve
    c(:,n+1) = Atilde \ btilde;
end

end