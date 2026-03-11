function b_swarm = PSO_update(i,vel_clamp,Kap,swarm,n_iter,n_part,xmin,xmax,verbose)
    % Function that updates the particles accross the iterations.
    % Parameters:
    % i: iteration
    % vel_clamp: Speed clamping
    % Kap: Value of the PSO gain
    % swarm: swarm that is evaluated in the iteration
    % n_iter: Number of iterations
    % n_part: Number of particles
    % xmin: Minimum value for the states
    % xmax: Maximum value for the states
    % verbose: 1 if the function should give verbose instructions

    if ~exist('b_fitness','var')
        b_fitness = zeros(n_iter,1); fitness = inf(1,n_part);
    end
    % Index of the global best particle
    [~, gbest] = min(swarm(:,4,1));
    gbest_val = swarm(gbest,4,1);
    % Updating velocity vectors and positions
    for n = 1:n_part
        for d = 1 : dim
            r1 = rand; r2 = rand;
            x = swarm(n,1,d);
            v = swarm(n,2,d);
            p = swarm(n,3,d); %pbest
            g = swarm(gbest,3,d); %gbest
            v = Kap*(v + c1*r1*(p - x) + c2*r2*(g - x));
            v = min(max(v, -vel_clamp(d)), vel_clamp(d));
            x = x + v;
            % Absorbing walls: When a particle hits the boundary of the
            % solution space, the velocity is zeroed in that dimension.
            if x < xmin(d)
                x = xmin(d);
                v = 0; % absorb: kill velocity in this dim
            elseif x > xmax(d)
                x = xmax(d);
                v = 0; % absorb: kill velocity in this dim
            end
            % % Update
            swarm(n,2,d) = v;
            swarm(n,1,d) = x;
        end
    end
    
    %swarms{i} = swarm;
    
    if i > 1
        b_fitness(i) = min(gbest_val, b_fitness(i - 1));
    
    else
        b_fitness(i) = gbest_val;
    end
    
    b_swarm = swarm(gbest,3,:);

    if verbose == 1
        display(['Iteration: ' num2str(i) '  Fitness: ' num2str(fitness(gbest)) ' Fitness(best): ' num2str(gbest_val)])
    end
end