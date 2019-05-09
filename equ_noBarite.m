function F = equ_noBarite(x,iniVolume,perVolume,drips,gammaSr,gammaSO4,Ba,Sr,SO4,Kd_BainSr,Ksp_celestite)

% x(1) = Ba_tot, x(2) = Sr_tot, x(3) = SO4_tot (amounts), 
% x(4) = delta_BainSr, x(5) = delta_SrSO4

F(1) = x(2)*x(3)*gammaSr*gammaSO4 - ((iniVolume + perVolume*drips)^2)*Ksp_celestite;
F(2) = Ba - x(4) - x(1);
F(3) = Sr - x(5) + x(4) - x(2);
F(4) = SO4 - x(5) - x(3); % 12/2, remove x(4) because F(3) & F(4) need to match
F(5) = x(4)*x(2) - Kd_BainSr*x(5)*x(1);
end