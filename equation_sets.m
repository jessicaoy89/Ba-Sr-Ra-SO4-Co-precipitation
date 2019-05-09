function F = equation_sets(x,iniVolume,perVolume,drips,gammaBa,gammaSr,gammaSO4,Ba,Sr,SO4,Kd_SrinBa,Kd_BainSr,Ksp_barite,Ksp_celestite)

% %%%%%% changed Ra input
% x(1) = Ba_tot, x(2) = Sr_tot, x(3) = SO4_tot (mmole at the end of precipitation), 
% x(4) = delta_SrinBa, x(5) = delta_BainSr, x(6) = delta_BaSO4, x(7) = delta_SrSO4

F(1) = x(1)*x(3)*gammaBa*gammaSO4 - ((iniVolume + perVolume*drips)^2)*Ksp_barite;
F(2) = x(2)*x(3)*gammaSr*gammaSO4 - ((iniVolume + perVolume*drips)^2)*Ksp_celestite;
F(3) = Ba - x(6) - x(5) + x(4) - x(1);
F(4) = Sr - x(7)- x(4) + x(5) - x(2);
F(5) = SO4 - x(7)- x(6) - x(3);  %12/2/2018 - change F(5), remove x(4) & x(5), because MSO4 is composed of both pure sulfate & M1inM2
F(6) = x(4)*x(1) - Kd_SrinBa*x(6)*x(2);
F(7) = x(5)*x(2) - Kd_BainSr*x(7)*x(1);

end