function [] = sim_3D_model(idx, sigma)

% Simulate 3-D D'Orsogna model with params iCr, iLr
% Last modified: Isaac Nathoo (Sep 20, 2019, use iterative solver to add noise)
% Last modified: Subhanik Purkayastha (Sep 19, 2019, reformatted for CCV)
% Last modified: Katie Storey (Aug, 1, 2018, added order paramter code)
% Last modified: John T. Nardini (Jul, 3, 2018, 3D simulation code)

% set index values
iCr = ceil(idx/5);
iLr = rem(idx,5)+1;

% set RNG seed based on param value
seed = 17 + 23*(iCr) + iLr;
rng(seed)

% time-step for iterative solver
deltat = 0.01;

% create directories to store results
mkdir(['data_3d_iCr_' num2str(iCr) '_iLr_' num2str(iLr)])
mkdir(['order_3d_data_iCr_' num2str(iCr) '_iLr_' num2str(iLr)])

% simulate 100 realizations
for i = 1:100

    % number of particles
    N = 200;

    % setting for odesolver
    Opt = odeset('AbsTol', 1e-2);

    % parameter values except cR and lR (compare D'Orsognal Model)
    % note that since we're setting cA=lA=1, cR and lR can be interpreted as the
    % ratios C=cR/cA and l=lR/lA
    alpha = 1.5;
    beta = 0.5;
    cA = 1.0;
    lA = 1.0;

    % time range
    tfinal = 2000; 		% final time in seconds
    tspan = 0:1:tfinal; 	% which time points are to be saved
    
    % initial positions - uniform on [-1,1]
    x = -2*rand(N,1) + 1;
    y = -2*rand(N,1) + 1;
    z = -2*rand(N,1) + 1;

    % initial velocities - uniform on [-1,1]
    vx = -2*rand(N,1) + 1;
    vy = -2*rand(N,1) + 1;
    vz = -2*rand(N,1) + 1;
    
    % select cR, lR
    cRvals=[0.1 0.5 0.9 2 3]; %create the cR range .. replusion strength
    lRvals=[0.1 0.5 0.9 2 3]; %create the lR range .. repulsion range
    cR = cRvals(iCr);
    lR = lRvals(iLr);
     
    % put it all in one big vector of ICs
    z0 = [x; y; z; vx; vy; vz];

    % simulate
    [t, q] = sim_swarm_ode_3d_rhs_iterative(tspan, deltat, z0, alpha, beta, cA, cR, lA, lR, sigma);
    
    % save results
    save(strcat('./data_3d_iCr_', num2str(iCr), '_iLr_', num2str(iLr),'/data_3d_iCr_', num2str(iCr), '_ilR_',num2str(iLr),...
    '_iR_',num2str(i),'.mat'));
             
    % number of time steps
    max_t=length(t);

    % rename pos,vel matrix to capital Z (so we can use z for third component)
    Z = q;
 
    % find x,y,z position and velocity vectors
    x = Z(:,1:N);
    y = Z(:,N+1:2*N);
    z = Z(:,2*N+1:3*N);
    vx = Z(:,3*N+1:4*N);
    vy = Z(:,4*N+1:5*N);
    vz = Z(:,5*N+1:6*N);

    % norm of the velocity
    v_norm = sqrt(vx.^2 + vy.^2 + vz.^2);

    % creates col vector of center of mass (position) coordinates for each t
    cx = sum(x,2)./N;
    cy = sum(y,2)./N;
    cz = sum(z,2)./N;
            
    % creates col vector of center of mass (velocity) for each t
    c_vx = sum(vx,2)./N;
    c_vy = sum(vy,2)./N;
    c_vz = sum(vz,2)./N;

    % r = distance between particle position and c.o.m. (rx_i= x_i-cx)
    rx = zeros(max_t,N);
    ry = zeros(max_t,N);
    rz = zeros(max_t,N);

    % cross-product of r and v
    rvcross = zeros(max_t,N,3);
    for s=1:N
        rx(:,s) = x(:,s)-cx;
        ry(:,s) = y(:,s)-cy;
        rz(:,s) = z(:,s)-cz;
        r_i=[rx(:,s) ry(:,s) rz(:,s)];
        v_i=[vx(:,s) vy(:,s) vz(:,s)];
        rvcross(:,s,:)=cross(r_i,v_i);

    end
    r_norm = sqrt(rx.^2 + ry.^2 + rz.^2);

    % polarization vector 
    P = sqrt((sum(vx,2)./sum(v_norm,2)).^2.+(sum(vy,2)./sum(v_norm,2)).^2.+(sum(vz,2)./sum(v_norm,2)).^2);

    % angular momentum 
    Mang_mat = sum(rvcross,2)./sum(v_norm.*r_norm,2);
    Mang = sqrt(Mang_mat(:,1,1).^2 + Mang_mat(:,1,2).^2 + Mang_mat(:,1,3).^2);

    % absolute angular momentum
    rvcross_norm = sum(sqrt(rvcross(:,:,1).^2+rvcross(:,:,2).^2+rvcross(:,:,3).^2),2);
    Mabs = sum(abs(rvcross_norm./sum(v_norm.*r_norm,2)),2);

    % average nearest neighbor distance
    nnd = zeros(max_t,1);
            
    % 3-D order parameter I_s = I_flock - I_mill
            
    % calculating I_mill and nnd at each t
    I_mill = zeros(max_t,1);

    for v = 1:max_t
        
	% x,y vectors at time i
        x_vec=x(v,:);
        y_vec=y(v,:);
        z_vec=z(v,:);

        % make column vectors
        x_vec = x_vec(:);
        y_vec = y_vec(:);
        z_vec = z_vec(:);

        % distance martrices in x and y
        Xdiff = x_vec - x_vec';
        Ydiff = y_vec - y_vec';
        Zdiff = z_vec - z_vec';

        % total distance matrix in R^2
        D = sqrt(Xdiff.^2 + Ydiff.^2+ Zdiff.^2);
        D(1:size(D,1)+1:end) = Inf;
        nnd(v) = mean(min(D,[],2));
                
        % Morse potential (eqn 4 in 3D paper)
        mp=(cR*exp(-D/lR)-cA*exp(-D/lA))/2;
        
        % derivative of Morse potential 
        dmp=cA/lA*exp(-D/lA) - cR/lR*exp(-D/lR);
        dmp_x= -sum(dmp.*Xdiff./(D+eps),2);
        dmp_y= -sum(dmp.*Ydiff./(D+eps),2);
        dmp_z= -sum(dmp.*Zdiff./(D+eps),2);
                
        % particle-particle force interaction (x,y,z vector for each
        % particle)
        f_mat = [dmp_x dmp_y dmp_z];
                   
        % velocity vector (at time i) for each particle 
        vx_i = vx(v,:);
        vy_i = vy(v,:);
        vz_i = vz(v,:);
        v_mat = [vx_i(:) vy_i(:) vz_i(:)];
               
        % rotational axis (eqn D3 of 3-D paper)
        omega = cross(v_mat,f_mat)./(sqrt(v_mat(:,1).^2+v_mat(:,2).^2+v_mat(:,3).^2).*sqrt(f_mat(:,1).^2+f_mat(:,2).^2+f_mat(:,3).^2));
        
	% calculate the degree of alignment between all omega
        for j = 1:N
            align_j = sum((sum(omega(j,:).*omega,2))) - sum(omega(j,:).*omega(j,:),2);
            I_mill(v) = I_mill(v) + align_j;
        end
                
        I_mill(v) = I_mill(v)/(N*(N-1));
	
     end
            
     % calculating I_flock
     % v_diff = difference between velocity and vel. center of mass
     vx_diff = zeros(max_t,N);
     vy_diff = zeros(max_t,N);
     vz_diff = zeros(max_t,N);
     
     for x = 1:N
        vx_diff(:,x) = vx(:,x)-c_vx;
        vy_diff(:,x) = vy(:,x)-c_vy;
        vz_diff(:,x) = vz(:,x)-c_vz;
     end
            
     I_flock = 1-(sum(sqrt(vx_diff.^2+vy_diff.^2+vz_diff.^2),2))/(N*sqrt(alpha/beta));
            
     I_s = I_flock - I_mill;
           
     % matrix with all 5 order params for each t
     order_par_mat = [P Mang Mabs nnd I_s];

     % create table of order parameter values for each time unit
     prelim_table = [t(1:10:end) P(1:10:end) Mang(1:10:end) Mabs(1:10:end) nnd(1:10:end) I_s(1:10:end)];
     
     % set Nan values to zero
     prelim_table(isnan(prelim_table)) = 0;
            
     % fill in with repeated values at final time (if simulation stopped early)
     prelim_table_height = size(prelim_table,1);
     if prelim_table_height < 2001
        prelim_table = [prelim_table ; zeros(2001-prelim_table_height,6)];
        prelim_table(prelim_table_height+1:end,2:6) = repmat([P(end) Mang(end) Mabs(end) nnd(end) I_s(end)],2001-prelim_table_height,1);
        prelim_table(:,1) = 0:2000;
     end
            
     final_table = prelim_table;
            
     % save data as a csv file
     write_file = ['./order_3d_data_iCr_' num2str(iCr) '_iLr_' num2str(iLr) '/order_3d_data_iCr_' num2str(iCr) '_iLr_' num2str(iLr) '_iR_' num2str(i) '_order_params.csv'];          
     csvwrite(write_file,final_table, 0, 0)
     
end

end
