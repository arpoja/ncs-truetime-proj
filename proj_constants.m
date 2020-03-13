%% System constants:
M = 0.5;    % mass cart [kg]
m = 0.4;    % mass pendulum [kg]
b = 0.1;    % friction coeff.[N/m/s]
l = 0.25;   % pendulum center of mass [m]
I = 0.025;  % pendulum inertia [kgm^2]
g = 9.81;   % gravity [m/s^2]

%% Continuous system model
den = I*(M+m)+M*m*l^2;
Ac = [0, 1,                 0,              0;
      0, -(I+m*l^2)*b/den,  m^2*g*l^2/den,  0;
      0, 0,                 0,              1;
      0, -m*l*b/den,        m*g*l*(M+m)/den,0];
Bc = [0;(I+m*l^2)/den; 0; m*l/den];
Cc = [1 0 0 0;
      0 0 1 0];
Dc = 0;
sysc = ss(Ac,Bc,Cc,Dc);

%% Discretizing constants
Q = diag([1 1 1000 1]); % state cost
R = 1;                  % input cost

%% Filtering constants
Pq = 0.0000001 * eye(4); % system noise covar
Pr = 0.0000001 * eye(2); % measurement noise covar

%% cleanup
clear den;