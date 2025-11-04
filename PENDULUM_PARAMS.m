% - Pendulum Parameters
mp = EstimatedParameters(2); % mass
xp = EstimatedParameters(1); % distance of cg from cart
Ip = EstimatedParameters(3); % moment of intertia about cg
kt = EstimatedParameters(4); % rotational damping
g = EstimatedParameters(5); % gravitational acceleration

% - Sampling Parameters
fs = 60; % angle sampling frequency (Hz)
measurement_variance = 1e-4; % noise

% - Initial Conditions
theta_d_0 = 0;
theta_0 = 30 * pi / 180;

% - LQR

A = 1 / (Ip + mp*xp^2) * [   0,    Ip + mp*xp^2;
                          mp*g*xp,    -kt       ];

B = 1 / (Ip + mp*xp^2) * [0; mp*xp];

Q = [100, 0;
    0, 1];
R = 10; 
N = zeros(2, 1);

K = lqr(A, B, Q, R, N);

% - Kalman Filter
W = [1e-3, 0;
      0,   1e-3];

V = measurement_variance;

P0 = zeros(2, 2);