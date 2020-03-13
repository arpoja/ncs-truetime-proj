function baseline_kernel_init(samplingTime)

ttInitKernel('prioFP');
% run script for constants:
proj_constants;

starttime = 0.0;
period = samplingTime;

% Discretize system
sysd = c2d(sysc,samplingTime,'zoh');    % samplingTime comes from outside

% controller values
[K_lqr,S_lqr,P_lqr]= lqr(sysd,Q,R,0);

% Stationary Kalman filter
[P_kalman,K_kalman,L_kalman] = idare(sysd.A',sysd.C',Pq,Pr,0,eye(size(sysd.A)));
K_kalman = K_kalman'; % transpose to get correct orientation

% filtering initialization
X0 = [0 0 0 0]'; % initial state
Xm0 = gaussian_noise(X0,100*Pq); % initial estimate with "large" covariance
                                % "large" wrt. system nosie

% Data structure for controller/filter
data.K_lqr = K_lqr; % LQR gain for controller
data.K_kalman = K_kalman;   % Kalman gain for filter
data.A = sysd.A;    % system matrices
data.B = sysd.B;    % -
data.C = sysd.C;    % - 
data.D = sysd.D;    % -
data.Xm = Xm0;      % filter estimate
data.u = 0;         % controller output
data.ts = period;   % sampling time

% create periodic task that does filtering and control
ttCreatePeriodicTask('filt_and_ctrl_task',starttime,period,'filt_and_ctrl',data);