%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Networked Control Systems project
%%% 2020-03-05
%%% Group E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clearvars; clc;
%% functions
rmse = @(x,hat) sqrt(mean((x - hat).^2));

%% System model and constants:
M = 0.5;    % mass cart [kg]
m = 0.4;    % mass pendulum [kg]
b = 0.1;    % friction coeff.[N/m/s]
l = 0.25;   % pendulum center of mass [m]
I = 0.025;  % pendulum inertia [kgm^2]
g = 9.81;   % gravity [m/s^2]

den = I*(M+m)+M*m*l^2;
Ac = [0, 1,                 0,              0;
      0, -(I+m*l^2)*b/den,  m^2*g*l^2/den,  0;
      0, 0,                 0,              1;
      0, -m*l*b/den,        m*g*l*(M+m)/den,0];
Bc = [0;(I+m*l^2)/den; 0; m*l/den];
Cc = [1 0 0 0;
      0 0 1 0];
%Cc = eye(4);
Dc = 0;

sysc = ss(Ac,Bc,Cc,Dc);

clear den;

%% Discretizing
t_sample = 0.015;   % sampling time [s] TODO: bandwidth?
Q = diag([1 1 1000 1]); % state cost
R = 1;                  % input cost

% c2d to discretize system
sysd = c2d(sysc,t_sample,'zoh');
Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
Dd = sysd.D;

% transfer function for simulink
[Dn,Dd] = ss2tf(Ad,Bd,Cd,Dd);

% Controller
% Kd: LQR controller gain
% Sd: solution to LQR riccatti equation
% Pd: poles of LQR controller
[Kd, Sd, Pd]= lqr(sysd,Q,R,0);  

% Stationary Kalman filter
Qx = 0.0000001 * eye(4);  % system noise  TODO: check with 0 noise?
% Qx = zeros(4);
Ry = 0.0000001 * eye(2);    % measurement noise
% Ps: Stationary covariance
% Ks: Stationary Kalman gain
% Ls: Stationary poles
[Ps,Ks,Ls] = idare(Ad,Cd',Qx,Ry,0,eye(4));
Ks = Ks';

% Filter
X0 = [0;0;0;0]; % initial true state
%X_corr = [0;0;pi/2;0];  % correction for LQR stability

P0 = 0.00001 * eye(4);  % initial covariance
Xm0 = X0 + chol(P0) * randn(4,1);
Xm = Xm0; % for matlab simulation
%Xm = [0;0;chol(P0(3,3)) * randn(1,1);0]; % only angle is uncertain
%m_est = X0;% + chol(P0) * randn(4,1);  % initial guess for state
%Pk = P0; % initialize covariance
%X0 = [X0;chol(P0) * randn(4,1)];


endtime = 30;   % end time in seconds
steps = endtime/t_sample;  % stepcount for sim
% allocate memory for filter
%X = zeros(4,steps); % true state
%Y = zeros(2,steps); % measurements
%U = zeros(1,steps); % input
%X_est = X;                  % estimate by kalman filter
%P_est = zeros(4,4,steps);   % filter covariance stored

%Phi = [ Ad - Bd*Kd      , Bd * Kd;
%        zeros(size(Ad)) , Ad - Ks * Cd];
%C = [Cd, zeros(size(Cd))];

%% Simulation
U = 0;

for k = 1:steps
    if k == 1   % with initial state
        %X(:,k) = X0 + Bd * U(:,k) +  chol(Qx) * randn(4,1);
        X(:,k) = Ad * X0 + Bd*U(:,k) + chol(Qx) * randn(4,1);
        
    else        % with last state
        %X(:,k) = Phi * X(:,k-1) + [chol(Qx) * randn(4,1); zeros(4,1)];
        X(:,k) = Ad * X(:,k-1) + Bd*U(:,k) + chol(Qx) * randn(4,1);
    end
    Y(:,k) = Cd * (X(:,k)) + chol(Ry) * randn(2,1);
    Xm0 = (Ad - Ks * Cd * Ad) * Xm0 + Bd * U(:,k) + Ks*Y(:,k);
    X_est(:,k) = Xm0;
    % break simulation if fallen over
    if  abs(X(3,k)) >= pi/8
        fprintf(sprintf('fell over at step %d, time %.2f\n',k,(k-1) * t_sample))
        break;
    end
    % calculating input for next cycle
    if k <= steps -1
        U(:,k+1) = -1 * Kd * Xm0;% + Lc * ref;
    end
    
       
    
%     % prediction
%     m_est = Ad * (m_est) + Bd * U(:,k);
%     Pk = Ad * Pk * Ad' + Qx;
%     
%     % update
%     Sk = Cd * Pk * Cd' + Ry;  % innovation covariance
%     Kk = Pk * Cd' * inv(Sk); % Kalman gain
%     m_est = m_est + Kk * (Y(:,k) - Cd * (m_est)); % update mean
%     Pk = Pk - Kk * Sk * Kk';    % update covariance
%     % store kalman filte results
%     X_est(:,k) = m_est;
%     P_est(:,:,k) = Pk;
%    
end
RMSE_X = rmse(X(1,:),X_est(1,:))
RMSE_Theta = rmse(X(3,:),X_est(3,:))


%% Plotting
t = 0:t_sample:(k - 1)*t_sample;
figure(1)
% angle
subplot(4,1,[1 2])
plot(t,Y(2,:),'.','Color','#AAAAAA')
hold on
plot(t,X(3,:))  % true angle
plot(t,X_est(3,:)) % filter

title('Angle deviation from upright')
legend('Measurement','True','Estimate','Location','best')
ylabel('Angle (Rad)')
xlabel('Time (s)')

% position
subplot(4,1,3)
plot(t,Y(1,:),'.','Color','#AAAAAA')
hold on
plot(t,X(1,:))
plot(t,X_est(1,:))
title('Position')
legend('Measurement','True','Estimate','Location','best')
ylabel('Pos. (m)')
xlabel('Time (s)')

% input
subplot(4,1,4)
plot(t,U)
title('Input')
legend('Input force','Location','best')
ylabel('Force (N)')
xlabel('Time (s)')


