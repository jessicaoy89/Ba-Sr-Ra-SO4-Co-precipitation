function [gammaBa,gammaSr,gammaRa,gammaSO4] = gamma_pitzer(I2,na,ca,cl,ba,sr,so4)
b0_NaCl = 0.07534; b1_NaCl = 0.2769; C_NaCl = 0.00074;
b0_KCl = 0.04808; b1_KCl = 0.2168; C_KCl = -0.000394;
b0_CaCl = 0.3159; b1_CaCl = 1.614; C_CaCl = -0.0000495;
b0_SrCl = 0.28575; b1_SrCl = 1.66725; C_SrCl = -0.0004596;
b0_BaCl = 0.5268; b1_BaCl = 0.687; C_BaCl = -0.0505581;
b0_RaCl = 0.050628; b1_RaCl = 0.294487; C_RaCl = 0.0;  
b0_NaSO4 = 0.0273; b1_NaSO4 = 0.956; C_NaSO4 = 0.0012085;
b0_CaSO4 = 0.0; b1_CaSO4 = 3.546; b2_CaSO4 = -59.3; C_CaSO4 = 0.0285;
theta_KCa = -0.00535; theta_KNa = -0.012; theta_ClSO4 = 0.03;
phi_NaCaSO4 = -0.055; phi_KCaCl = -0.025; phi_KNaCl = -0.0015; phi_NaCaCl = -0.0148;
phi_ClSO4Na = 0; phi_ClSO4Ca = -0.122;
A = 0.5085; A_phi = 0.3915; B = 3.281; a_i = 0.5; 
alpha = 2.0; b = 1.2; alpha1 = 1.4; alpha2 = 12.0;

gI = (2.0/(alpha*alpha*I2))*(1.0-(1.0+alpha*sqrt(I2))*exp(-1.0*alpha*sqrt(I2)));
ggI = (-2.0/(alpha*alpha*I2))*(1.0-(1.0+alpha*sqrt(I2) + 0.5*alpha*alpha*I2)*exp(-1.0*alpha*sqrt(I2)));
B_NaCl = b0_NaCl + gI*b1_NaCl;
Bb_NaCl = ggI*b1_NaCl/I2;
B_CaCl = b0_CaCl + gI*b1_CaCl;
Bb_CaCl = ggI*b1_CaCl/I2;
B_BaCl = b0_BaCl + gI*b1_BaCl;
B_SrCl = b0_SrCl + gI*b1_SrCl;
B_RaCl = b0_RaCl + gI*b1_RaCl;
B_NaSO4 = b0_NaSO4 + gI*b1_NaSO4;
Bb_NaSO4 = ggI*b1_NaSO4/I2;
B_KCl = b0_KCl + gI*b1_KCl;

g_alpha1 = (2.0/(alpha1*alpha1*I2))*(1.0-(1.0+alpha1*sqrt(I2))*exp(-1.0*alpha1*sqrt(I2)));
gg_alpha1 = (-2.0/(alpha1*alpha1*I2))*(1.0-(1.0+alpha1*sqrt(I2)+0.5*alpha1*alpha1*I2)*exp(-1.0*alpha1*sqrt(I2)));
g_alpha2 = (2.0/(alpha2*alpha2*I2))*(1.0-(1.0+alpha2*sqrt(I2))*exp(-1.0*alpha2*sqrt(I2)));
gg_alpha2 = (-2.0/(alpha2*alpha2*I2))*(1.0-(1.0+alpha2*sqrt(I2)+0.5*alpha2*alpha2*I2)*exp(-1.0*alpha2*sqrt(I2)));
B_CaSO4 = b0_CaSO4 + b1_CaSO4*g_alpha1 + b2_CaSO4*g_alpha2;
Bb_CaSO4 = b1_CaSO4*gg_alpha1+b2_CaSO4*gg_alpha2;
f_gamma = -1.0*A_phi*(sqrt(I2)/(1+b*sqrt(I2)) + (2.0/b)*log(1.0+b*sqrt(I2))); 
Z = na + cl + 2*ca + 2*ba + 2*sr + 2*so4;

% calculate KCl replacements:
F_KCl = (f_gamma + na*so4*Bb_NaSO4 + ca*so4*Bb_CaSO4);
F_MSO4 = 4.0*(f_gamma + na*cl*Bb_NaCl + ca*cl*Bb_CaCl);
F_KSO4 = 2.0*(f_gamma + na*cl*Bb_NaCl + ca*cl*Bb_CaCl);
same_MSO4 = (0.5)*(na*(2.0*B_NaSO4+Z*C_NaSO4) + ca*(2.0*B_CaSO4+Z*C_CaSO4)) + na*cl*(0.5)*(4*C_NaCl+phi_ClSO4Na) + ca*cl*(0.5)*(4*C_CaCl+phi_ClSO4Ca) + na*ca*(0.5)*phi_NaCaSO4; 
	
ln_KCl = F_KCl + 0.5*cl*(2*B_KCl + Z*C_KCl)+ 0.5*(na*(2.0*B_NaCl+Z*C_NaCl+2.0*theta_KNa) + ca*(2.0*B_CaCl+Z*C_CaCl+2.0*theta_KCa)) + 0.5*na*ca*phi_NaCaCl;
ln_KSO4 = F_KSO4 + (0.6667)*cl*(2*B_KCl+Z*C_KCl+2*theta_ClSO4) + (0.3333)*(na*(2*B_NaSO4+Z*C_NaSO4+2*2*theta_KNa) + ca*(2*B_CaSO4+Z*C_CaSO4+2*2*theta_KCa)) + na*cl*(0.3333)*(2*2*C_NaCl+2*phi_KNaCl+phi_ClSO4Na) + ca*cl*(0.3333)*(2*2*C_CaCl+2*phi_KCaCl+phi_ClSO4Ca) + na*ca*(0.3333)*phi_NaCaSO4;
	
temp_KCl = exp(ln_KCl);
temp_KSO4 = exp(ln_KSO4);

% calculate gammas:
ln_BaSO4 = F_MSO4 + same_MSO4 + (0.5)*cl*(2*B_BaCl+Z*C_BaCl+2*theta_ClSO4);
gamma_tempBa = exp(ln_BaSO4);
gammaBa = gamma_tempBa*gamma_tempBa*temp_KCl*temp_KCl/(temp_KSO4*temp_KSO4*temp_KSO4);

ln_SrSO4 = F_MSO4 + same_MSO4 + (0.5)*cl*(2*B_SrCl+Z*C_SrCl+2*theta_ClSO4);
gamma_tempSr = exp(ln_SrSO4);
gammaSr = gamma_tempSr*gamma_tempSr*temp_KCl*temp_KCl/(temp_KSO4*temp_KSO4*temp_KSO4);

ln_RaSO4 = F_MSO4 + same_MSO4 + (0.5)*cl*(2*B_RaCl+Z*C_RaCl+2*theta_ClSO4);
gamma_tempRa = exp(ln_RaSO4);
gammaRa = gamma_tempRa*gamma_tempRa*temp_KCl*temp_KCl/(temp_KSO4*temp_KSO4*temp_KSO4);

gamma_tempSO4 = exp(ln_KCl);
gammaSO4 = temp_KSO4*temp_KSO4*temp_KSO4/(temp_KCl*temp_KCl);

end