function F = equ_noCeles(x,iniVolume,perVolume,drips,gammaBa,gammaSO4,Ba,Sr,SO4,Kd_SrinBa,Ksp_barite)

% x(1) = Ba_tot, x(2) = Sr_tot, x(3) = SO4_tot (amounts), 
% x(4) = delta_SrinBa, x(5) = delta_BaSO4

F(1) = x(1)*x(3)*gammaBa*gammaSO4 - ((iniVolume + perVolume*drips)^2)*Ksp_barite;
F(2) = Ba - x(5) + x(4) - x(1);
F(3) = Sr - x(4) - x(2);
F(4) = SO4 - x(5) - x(3); % 12/2, remove x(4) because F(3) & F(4) need to match
F(5) = x(4)*x(1) - Kd_SrinBa*x(5)*x(2);
end