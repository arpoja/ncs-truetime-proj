function [exectime, data] = filt_and_ctrl(seg,data)
switch seg
    case 1
        meas_X = ttAnalogIn(1);
        meas_Theta = ttAnalogIn(2);
        meas = [meas_X meas_Theta]';
        data = filtercalc(data,meas);
        data = controlcalc(data);
        exectime = data.ts;
    case 2
        ttAnalogOut(1, data.u)
        ttAnalogOut(2, data.Xm(1))
        ttAnalogOut(3, data.Xm(3))
        exectime = -1;
end