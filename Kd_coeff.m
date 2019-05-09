function [Kd_RainBa,Kd_SrinBa,Kd_RainSr,Kd_BainSr] = Kd_coeff(gammaBa, gammaSr, gammaRa)
Ksp_celestite = 0.000000234;
Ksp_barite = 0.000000000107;
Ksp_RaSulfate = 0.000000000055;
W_SrinBa = 981.76; W_BainSr = 867.01; W_RainBa = 259.42; W_RainSr = 1891.39;
R = 1.987; T = 298;

Kd_RainBa = exp(log(Ksp_barite/Ksp_RaSulfate) + log(gammaRa/gammaBa) - W_RainBa/(R*T));
Kd_SrinBa = exp(log(Ksp_barite/Ksp_celestite) + log(gammaSr/gammaBa) - W_SrinBa/(R*T));
Kd_RainSr = exp(log(Ksp_celestite/Ksp_RaSulfate) + log(gammaRa/gammaSr) - W_RainSr/(R*T));
Kd_BainSr = exp(log(Ksp_celestite/Ksp_barite) + log(gammaBa/gammaSr) - W_BainSr/(R*T));

end