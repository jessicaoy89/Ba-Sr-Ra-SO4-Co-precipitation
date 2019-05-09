function [gammaBa,gammaSr,gammaRa,gammaSO4] = gamma_sit(I1,na,cl)
A = 0.5085; B = 3.281; a_i = 0.5; sqZ = 4.0;

eBa = 0.086+(0.012-0.086)/(I1+1);
gammaBa = 10^(-1*A*sqZ*sqrt(I1)/(1+1.5*sqrt(I1)) + eBa*cl);

eSr = 0.116+0.0046*I1 + (0.084-0.116-0.0046*I1)/(1+I1);
gammaSr = 10^(-1*A*sqZ*sqrt(I1)/(1+1.5*sqrt(I1)) + eSr*cl);

eRa = 0.057;
gammaRa = 10^(-1*A*sqZ*sqrt(I1)/(1+1.5*sqrt(I1)) + eRa*cl);

eSO4 = -0.09 + 0.0052*I1 + (-0.386 + 0.09 - 0.0052*I1)/(1+I1);
gammaSO4 = 10^(-1*A*sqZ*sqrt(I1)/(1+1.5*sqrt(I1)) + eSO4*na);

end