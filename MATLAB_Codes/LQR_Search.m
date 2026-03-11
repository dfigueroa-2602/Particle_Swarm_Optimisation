function fitness = LQR_Search(Params,res_mode,Sim_sys,t_settling,t_step,Ts,R_mode,sw,beta_c)
    % Function that injects the given swarm into the lqr function, then call
    % the simulation function.
    % Parameters:
    % Params: Structure. It should contain:
    %   - Ad
    %   - Bd
    %   - Ard if resonant terms are used
    %   - Brd if resonant terms are used
    %   - Hx
    %   - w
    %   - Vdc
    
    Ad = Params.Ad; Bd = Params.Bd; Hx = Params.Hx; w = Params.w; u_sat = Params.Vdc;
    if ~isempty(Params.Ard)
        Ard = Params.Ard;
        Brd = Params.Brd;
    end

    nx  = size(Bd,1); 
    nu  = size(Bd,2);

    idx_x_end = nx/2;
    sw_x  = sw(1:idx_x_end);
    Q = diag(10.^repelem(sw_x,2));

    if R_mode == 1
        Ru = repelem(sw(end - nu/2 + 1:end),2);
        R = diag((10.^Ru)');
    else
        R = diag(repelem(1,nu));
    end

    if any(R(:) < 0,1) || any(Q(:) < 0,1)
        disp('One of the elements of the Q and R matrixes is negative!')
        fitness = 1e6; return;
    end
    
    try
        [K,~,~] = dlqr(A_aug, B_aug, Q, R);
    catch
        disp('Error in the LQR')
        fitness = 1e6; return;
    end

    Acl = A_aug - B_aug*K;
    spec = max(abs(eig(Acl)));
    if spec >= 1
        fitness = 1e6;
        disp('One of the particles is too close of the instable area!')
        return;
    end

    try
        fitness = Sim_ResSystem(Tsim,Ts,VDC,w,A_ref,T_settling,beta_c,Ad,Bd,Ard,Brd,Ha,Kx,Kr,Ku);
    catch ME
        fprintf('\nSimulation failed: %s\n', ME.message)
        disp(getReport(ME, 'extended'))
        fitness = 1e6;
        return;
    end
end