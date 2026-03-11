Variables;

ss_VSI = ss(A,B,eye(nx),[]); [Ad,Bd,~,~] = ssdata(c2d(ss_VSI,Ts,'zoh'));

Ad_delay = [Ad Bd; zeros(nu,nx) zeros(nu)]; Bd_delay = [zeros(nx,nu);eye(nu)];

Ar1 = [0 w; -w 0]; Br1 = [0;1];
[Ard1,Brd1,~,~] = ssdata(c2d(ss(Ar1,Br1,eye(2),[]),Ts,'tustin'));
Ard1 = blkdiag(Ard1,Ard1); Brd1 = blkdiag(Brd1,Brd1);

w_vec = w*[1]; n_h = length(w_vec);

Ard = blkdiag(Ard1); Brd = [Brd1];
nxr = size(Brd,1); nur = size(Brd,2);

C_delay = [C zeros(nx,nu)];

% currentFolder = pwd;
% C_Folder = extractBefore(currentFolder,'MATLAB');
% C_folder = append(C_Folder,'PLECS/');
% cd(C_folder)
% save 3L_Resonant.mat Ard Brd
% cd(currentFolder)

Ha = [0 0 1 0 0 0 ;
      0 0 0 1 0 0];

Hx = C_delay(3:4,:);

Ad_aug = [Ad_delay  zeros((nx + nu),nxr) ;
          -Brd*Hx   Ard                 ];
Bd_aug = [Bd_delay; zeros(nxr,nu)];

Tsim = 0.1; u_sat = Vdc;
ref = [200;200]; K = dlqr(Ad_aug,Bd_aug,randi(100)*eye(size(Bd_aug,1)),randi(100)*eye(size(Bd_aug,2)));
t_settling = 20e-3; t_step = 20e-3; beta_c = 2e-7;
K
Sim_ab_sys(Ad_delay,Bd_delay,Ard,Brd,Ha,ref,K,w,u_sat,Tsim,Ts,t_settling,t_step,beta_c)

%%
search = 1;

beta_c = 3e-5; n_part = 100; iter = 50;
c1 = 2.05; c2 = 2.05; R_mode = 0;
Qmax = 4; Rmax = 2; Qmax_vel = 1; Rmax_vel = 0.1; rang_coef = 0.6;
nQ = nx/2 + nu/2 + n_h; nR = size(Bd_aug,2)/2; dim = nQ + (R_mode==1)*(nu/2);



%currentFolder = pwd;
%C_Folder = extractBefore(currentFolder,'MATLAB');
%C_folder = append(C_Folder,'PLECS/');
%cd(C_folder)
%save 3L_Gains.mat Kx Ku Kr
%cd(currentFolder);

%Kgain2Ccode({Kx, Ku, Kr, Ard, Brd, Ad_aug, Bd_aug},{'Kx_val', 'Ku_val', 'Kr_val', 'Ard_val', 'Brd_val', 'Ad_aug_val', 'Bd_aug_val'},'Matrices')