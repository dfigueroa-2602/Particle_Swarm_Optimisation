function cost = Sim_ab_sys(Ad,Bd,Ard,Brd,Hx,ref,K,w,u_sat,Tsim,Ts,t_settling,t_step,beta_c)
    % Function that simulates the closed-loop of a alpha-beta system 
    % controlled by a gain, which can be designed by LQR.
    % It gives the cost associated with reference tracking and input
    % control effort.
    % Parameters:
    % Ad: State matrix
    % Bd: Input matrix
    % Ard: Resonant states matrix
    % Brd: Resonant inputs matrix
    % Hx: Selector matrix for states with references
    % ref: Amplitude reference vector nx X 1
    % K: K gain designed by LQR
    % w: Sinusoidal frequency in rad/s
    % u_sat: Saturation value for the control input (for numerical
    % purposes)
    % Tsim: Simulation time
    % t_settling: Settling time for references
    % t_step: Step time
    % beta_c: Weight associated with the control effort cost
    
    nx = size(Bd,1); nu = size(Bd,2); n_ref = length(ref); nr = size(Brd,1);
    if n_ref ~= nu
        error('The number of references must be equal to the number of system inputs')
    end
    if mod(n_ref,2) ~= 0
        error('n_ref must be even because references are handled as alpha-beta pairs.');
    end

    % Time vector
    N = round(Tsim / Ts);               % Samples Calculator
    t = (0:N-1).' * Ts;
    % Reference Generator
    tau = t_settling / 4;
    env = 1 - exp(-t/tau);
    alpha_ref = env.*sin(w*t);
    beta_ref  = env.*cos(w*t);
    
    % Preallocate
    refs = zeros(N, n_ref);

    for ref_i = 1:n_ref
        if mod(ref_i,2) == 1
            refs(:,ref_i) = alpha_ref.*ref(ref_i);
        else
            refs(:,ref_i) = beta_ref.*ref(ref_i);
        end
    end
    
    refs(t < t_step,:) = 0; refs = refs.';

    % Preallocate states and inputs
    x_k  = zeros(nx,1);
    rho_k  = zeros(nr,1);
    xu_k = zeros(nu,1);
    u_prev  = zeros(nu,1);
    % Preallocate fitness functions
    e_log  = zeros(N,n_ref);
    du_log = zeros(N,nu);
    
    for k = 1:N
        % To take the states with references
        z_k = Hx*x_k;
        % References at this instant
        r_k = refs(:,k);
        % Error
        e_k = r_k - z_k;
        e_log(k,:) = e_k.';
        % Resonant System
        rho_k_1 = Ard*rho_k + Brd*e_k;
        % Control law
        u_k = -K*[x_k;rho_k];
        % Saturation for numerical purposes
        u_k_sat = max(min(u_k, u_sat), -u_sat);
        % Plant with delayed input
        x_k_1 = Ad*x_k + Bd*xu_k;
        % Difference of actuation calculation
        du_k = u_k - u_prev;
        du_log(k,:) = du_k.';
        u_prev = u_k_sat;
        
        % States update
        x_k = x_k_1;
        rho_k = rho_k_1;
        xu_k = u_k_sat;
    end
    
    term_e  = sum(e_log.^2, 2);                    % [N x 1]
    term_du = sum(du_log.^2, 2);                   % [N x 1]
    cost = (1/N) * sum(term_e + beta_c*term_du);
end