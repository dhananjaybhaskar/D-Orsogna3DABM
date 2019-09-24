function [t,output] = sim_swarm_ode_3d_rhs_iterative(tspan,deltat,Z,alpha,beta,cA,cR,lA,lR,sigma)
%{
sim_swarm_ode_3d_rhs.m written 3-7-18 by JTN as an update of 
sim_swarm_ode_rhs.m

This file sets the rhs of the D'Orsogna model in 3d. It is called by

- sim_swarm_main_3d.m

It takes as input

- t .. time
- Z .. a vector of size 6 * no. of particles containing the vector of locations (x) and velocities (v)
    Z = (x_1 , x_2 , ... , x_n , y_1, y_2 , ... , y_n, z_1 , ..., z_N, 
        vx_1 , ... , vx_n, vy_1, ... , vy_n, vz_1 , ... , vz_N)

- model parameter values alpha, beta, cR, lR, cA, lA

%}

%[t,z] = ode45(@(t,z) sim_swarm_ode_3d_rhs(t,z,alpha,beta,cA,cR,lA,lR),tspan,z0, Opt);
% ensure that Z is a column
Z = Z(:);

% find N, the number of particles
N = length(Z)/6;
N = int16(N);

% find x, y, z and velocity terms
x = Z(1:N) + (0 + (sigma)*randn(N,1));
y = Z(N+1:2*N) + (0 + (sigma)*randn(N,1));
z = Z(2*N+1:3*N) + (0 + (sigma)*randn(N,1));
vx = Z(3*N+1:4*N);
vy = Z(4*N+1:5*N);
vz = Z(5*N+1:6*N);

% normal distributions
%pdx = fitdist(x, 'Normal');
%sigmax = pdx.sigma * sigma;
%pdy = fitdist(y, 'Normal');
%sigmay = pdy.sigma * sigma;
%pdz = fitdist(z, 'Normal');
%sigmaz = pdz.sigma * sigma;

% new positions with noise
%x = x + sigmax;
%y = y + sigmay;
%z = z + sigmaz;

%set dt
dt = deltat;
tfinal = tspan(length(tspan));

% counter
c = tspan(1);
t = [];
output = [];

for i = 1:length(tspan)

    % provide ODE of the locations
    dxdt = vx;
    dydt = vy;
    dzdt = vz;
    
    % compute distance matrixes in x, y, and z (note that this requires a
    % relatively recent Matlab version. Otherwise use repmat to produce
    % matrices of the same dimension!
    Xdiff = x - x';
    Ydiff = y - y';
    Zdiff = z - z';

    % compute a total distance matrix in R^2
    D = sqrt(Xdiff.^2 + Ydiff.^2 + Zdiff.^2);

    % Get a vector for the norm^2 of the velocity
    v_normSq = vx.^2 + vy.^2 + vz.^2;

    % define the u'(r) function
    u_prime = cA/lA*exp(-D/lA) - cR/lR*exp(-D/lR);

    % compute the ODEs for the velocities
    dvxdt = (alpha - beta*v_normSq).*vx - sum(u_prime.*Xdiff./(D+eps),2);
    dvydt = (alpha - beta*v_normSq).*vy - sum(u_prime.*Ydiff./(D+eps),2);
    dvzdt = (alpha - beta*v_normSq).*vz - sum(u_prime.*Zdiff./(D+eps),2);
    
    % put into output vector
    if (tspan(i) == c)
        t = [t; c];
        output = [output; [x' , y' , z' ,  vx' , vy' , vz']];
        if (tspan(i) ~= tfinal)
            c = tspan(i + 2);
        end
    end
    
    % calculate new positions
    x = (x + dxdt*dt) + ((sigma)*randn(N,1));
    y = (y + dydt*dt) + ((sigma)*randn(N,1));
    z = (z + dzdt*dt) + ((sigma)*randn(N,1));
    
    % calculate new velocities
    vx = vx + dvxdt*dt;
    vy = vy + dvydt*dt;
    vz = vz + dvzdt*dt;
    
end

end
