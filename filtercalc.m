function data = filtercalc(data,meas)
% state estimate with stationary Kalman
data.Xm = (data.A - data.K_kalman * data.C * data.A) * data.Xm +...
    data.B * data.u + data.K_kalman * meas;

