function K = PSO_alg(sys,n_iter,n_part,nQ,nR,c1,c2,xmax_val,Q_max,R_max,Qmax_vel,Rmax_vel,R_mode,rang_coef,verbose)
    % Function that initialises the particle swarm optimisation (PSO)
    % algorithm.
    % Parameters:
    % sys: Plant that LQR will be used for
    % n_iter: Iterations number
    % n_part: Particles number
    % nQ: States weight number
    % nR: Input weight number
    % xmax_val: Maximum state value for absorbing walls
    % c1: PSO coefficient 1
    % c2: PSO coefficient 2
    % Q_max: Maximum state weight
    % R_max: Maximum input weight
    % Qmax_vel: Maximum speed for the particle associated with state weight
    % Rmax_vel: Maximum speed for the particle associated with input weight
    % R_mode: 0 if the R matrix is the identity, or 1 if PSO should search
    %         input weights
    % rang_coef: Range in which the particles will be searched
    % verbose: 1 if the function should give verbose instructions
    % The PSO algorithm requires the following toolboxes:
    %   - Statistics
    %   - Parallel Computing
    
    [Ad,Bd,~,~] = ssdata(sys);
    
    [vel_clamp,Kap,swarm,~] = Swarm_Init(n_part,R_mode,nQ,nR,c1,c2,Q_max,R_max,rang_coef,Qmax_vel,Rmax_vel);
    xmin = repmat(-xmax_val,dim,1); xmax = repmat(xmax_val,dim,1);
    
    rng(1,'twister'); if isempty(gcp('nocreate')), parpool; end
    fitness = inf(1,n_part);    % fitness value should be infinite
    for i = 1:n_iter            % Initialize particles
        parfor n = 1:n_part
            try
                sw = swarm(n,:,:);
                sw = squeeze(sw(1,1,:));
                fitness(n) = LQR_Search(Ts,n_h,w_vec,R_mode,sw,Ad_aug,Bd_aug, ...
                    Ad,Bd,Ard,Brd,Ha,beta_c,Vdc,w,A_ref,Tsim,T_settling);
            catch
                fitness(n) = 1e6;
                disp(['Evaluation for particle no. ' num2str(n) ' was aborted']);
            end
        end
        for n = 1:n_part
            if fitness(n) < swarm(n,4,1)
                swarm(n,3,:) = swarm(n,1,:);
                swarm(n,4,1) = fitness(n);
            end
        end
        b_swarm = PSO_update(i,vel_clamp,Kap,swarm,n_iter,n_part,xmin,xmax,verbose);
    end

    K = Final_Value(b_swarm,w,w_vec,Ad_aug,Bd_aug,nx,nu,nxr,n_h,R_mode);
    delete(gcp('nocreate'));
end