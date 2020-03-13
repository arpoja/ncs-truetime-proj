function data = controlcalc(data)
% calculates control signal
data.u = -1 * data.K_lqr * data.Xm;