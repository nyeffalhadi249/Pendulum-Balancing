clear; clc;
% - Pendulum Parameters
load("C:\Users\nyeff\OneDrive\Desktop\Pendulum_Cart\EstimatedParameters_SinglePendulum.mat")
mp = EstimatedParameters(2); % mass
xp = EstimatedParameters(1); % distance of cg from cart
Ip = EstimatedParameters(3); % moment of intertia about cg
kt = EstimatedParameters(4); % rotational damping
g = EstimatedParameters(5); % gravitational acceleration

% - Sampling Parameters
fs = 60; % sampling frequency (Hz)

% - Initial Conditions
theta_0 = 10 * pi / 180;
theta_d_0 = 0;
x_0 = 0;
x_d_0 = 0;

% - H_inf
Ibar = Ip + mp*xp^2;

A = 1 / (Ip + mp*xp^2) * [   0,    Ip + mp*xp^2;
                          mp*g*xp,    -kt       ];

A_aug = [0,              1,          0, 0;
         mp*g*xp / Ibar, -kt / Ibar, 0, 0;
         0,              0,          0, 1;
         0,              0,          0, 0];

B = 1 / (Ip + mp*xp^2) * [0; mp*xp];

B_u = [0; mp*xp / Ibar; 0; 1];
Bw = 0.0*eye(4);
B_hinf = [Bw, B_u];

C = [1, 0, 0, 0; 
     0, 0, 1, 0];

Wz = diag([100 5 5 1]);
Cz = [Wz;
      zeros(1,4)];

Cy = [1 0 0 0;
      0 0 1 0];

C_hinf = [Cz;
          Cy];

R = 0.01;
D11 = zeros(5,4);
D12 = [zeros(4,1);
       sqrt(R)];

Noise = [0.0001, 0.01];
D21 = diag(Noise)*eye(2,4);
D22 = zeros(2, 1);

D_hinf = [D11 D12;
          D21 D22];


ncon = 1;
nmeas = 2;
Pc = ss(A_aug, B_hinf, C_hinf, D_hinf);
P = c2d(Pc,1/fs);
[K, CL, gamma] = hinfsyn(P, nmeas, ncon);
[AK, BK, CK, DK] = ssdata(K);

