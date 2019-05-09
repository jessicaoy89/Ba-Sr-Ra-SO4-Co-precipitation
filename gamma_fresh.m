function [gamma] = gamma_fresh(IS)
A = 0.5085;
if IS < 0.005
    gamma = 10^(-1.0* A* 4* sqrt(IS));
elseif IS >= 0.005 && IS < 0.1
    gamma = 10^(-1.0* A* 4* sqrt(IS)/(1+ 3.281*0.5*sqrt(IS)));
elseif IS >= 0.1
    gamma = 10^(-1.0* A* 4* (sqrt(IS)/(1+ sqrt(IS)) - 0.3*IS));
end