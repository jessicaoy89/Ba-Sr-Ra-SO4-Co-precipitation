%% Titration experiment:
% Solution A, 1000 mL with Na2SO4 as reactant;
% Solution B, 100 mL with BaCl2 and SrCl2 as reactant;
% Both A & B have the same "background ions" as chosen in line 16 "conditions";
% Titrate solution A into solution B to start co-precipitation.
%% Revisions:
% 12/03/2018: Change the lower bound of delta_BaSO4.
% 12/07/2018: Change the input X0 of equation sets
% 02/03/2019: remove all hardwire to make it easy to change at one time
% 05/09/2019: added a looping/QC process
%% Define input values
format compact;
fprintf('Conditons: 1. DI; 2. 2M NaCl; 3. 2M NaCl + 0.5M CaCl2\n') % choose background condition:
prompt = 'What is the condition: ';
cont = input(prompt);
fprintf('Method: 1. FreshWater; 2. Pitzer; 3. SIT\n') % choose activity coefficient method:
prompt = 'Which method to use: ';
med = input(prompt);
fprintf('Choose how to enter Sr and Ba: 1. Ratio; 2. Input values\n') % choose how to enter Sr and Ba
prompt = 'How to enter: ';
how = input(prompt);
if how == 1
    fprintf('Sr/Ba ratios:\n 1. 0.5:5; 2. 5:5; 3. 5:0.5; 4. 50:0.5\n') % choose Sr/Ba ratios:
    prompt2 = 'What is the Sr/Ba ratio: ';
    ratio = input(prompt2);
    % Sr and Ba values are set in units of mmole (total aqueous amounts).
    % because the concentration changes over different drips.
    if ratio == 1
        Sr = 0.5; Ba = 5;
    elseif ratio == 2
        Sr = 5; Ba = 5;
    elseif ratio == 3
        Sr = 5; Ba = 0.5;
    elseif ratio == 4
        Sr = 50; Ba = 0.5;
    end
elseif how == 2
    prompt3 = 'Sr (mmole): ';
    Sr = input(prompt3);
    prompt4 = 'Ba (mmole): ';
    Ba = input(prompt4);
end
% background conditions is set in units of M/L.
% because their concentrations don't change over different drips.
if cont == 1
    Na_conc0 = 0; Ca_conc0 = 0; Cl_conc0 = 0; % DI water
elseif cont == 2
    Na_conc0 = 2; Ca_conc0 = 0; Cl_conc0 = 2; % 2 M NaCl
elseif cont == 3
    Na_cont0 = 2; Ca_conc0 = 0.5; Cl_conc0 = 3; % 2 M NaCl & 0.5 M CaCl2
end
% set initial variables according to inputs of concentration & ratios:
Ra = 10.0*2.7027e-8/226; % in units of mmole
iniVolume = 100.0; % 02/03/2019: added initial volume, which starts at 100 mL
perVolume = 0.1; % 1203: change the volume per drip to test different step sizes (from 0.1 to 0.05 to 0.025 mL per drip)
SO4_drip = 0.005*perVolume; % 1024: change SO4 amount per drip, unit is mmole;
add_Na = 0; % initial addition of Na together with SO4
SO4 = 0; % initial SO4 amount in 100 mL BaCl2 & SrCl2 solution
% initialize six sum values
sum_BaSO4 = 0; sum_SrSO4 = 0; sum_RainBa = 0; sum_RainSr = 0; sum_BainSr = 0; sum_SrinBa = 0;
% initialize X0 values for solving equations:
x0_BaSO4 = 0.0004; x0_SrSO4 = 0; x0_BainSr = 0; x0_SrinBa = 0;
% x0_BaSO4 is set at 0.0002 because otherwise it might be too small for the
% 1st drip
    
Sr0 = Sr; Ba0 = Ba; Ra0 = Ra;
Ksp_celestite = 0.000000234;
Ksp_barite = 0.000000000107;
Ksp_RaSulfate = 0.000000000055;
K_check = 0;
K_last = 0;
step_BaSO4 = 0; step_SrSO4 = 0;
step_BainSr = 0; step_SrinBa = 0;
step_RainSr = 0; step_RainBa = 0;

% create a matrix of 10000 or 20000 or 40000 rows depending on steps of choice to store results 
number_drips = 1000/perVolume; % change the maximum number of drips according to perVolume size
results = zeros(number_drips,31); % set a matrix dimension to record output of each drip
% define the options for solving lsqnonlin equations
options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','TolFun',1e-18,'OptimalityTolerance',1e-21,'TolX',1e-18);
    
%% Start looping 10000 or 20000 or 40000 drips of liquid
% In each loop, the Ba, Sr, Ra, SO4 (in mmole) remain in solution is re-calculated;
% the changes in solid amounts (BaSO4, SrSO4, BainSr,SrinBa,RainBa,RainSr) of each drip is calculated
for drips = 1:number_drips
    % Perform quality control - determing whether the drips will continue or will repeat the last drip.
    if K_last == 0
        SO4 = SO4 + SO4_drip;
        add_Na = add_Na + 2*SO4_drip;
    else
        SO4 = SO4;
        add_Na = add_Na;
    end
    branch_unfinished = 0;
    Ba_conc = Ba/(iniVolume + perVolume*drips);
    Sr_conc = Sr/(iniVolume + perVolume*drips);
    Ra_conc = Ra/(iniVolume + perVolume*drips);
    SO4_conc = SO4/(iniVolume + perVolume*drips);
    Na_conc = Na_conc0 + add_Na/(iniVolume + perVolume*drips);
    Ca_conc = Ca_conc0;
    Cl_conc = Cl_conc0 + 2*(Ba0 + Sr0 + Ra0)/(iniVolume + perVolume*drips);
    IS = 0.5*(4*(Ba_conc+Sr_conc+Ra_conc+SO4_conc+Ca_conc) + Na_conc + Cl_conc); % calculate ionic strength
    
    % calculate gamma values based on different choice of methods:
    if med == 1
        gammaBa = gamma_fresh(IS);
        gammaSr = gammaBa; gammaRa = gammaBa; gammaSO4 = gammaBa;
    elseif med == 2
        [gammaBa,gammaSr,gammaRa,gammaSO4] = gamma_pitzer(IS,Na_conc,Ca_conc,Cl_conc,Ba_conc,Sr_conc,SO4_conc);
    elseif med == 3
        [gammaBa,gammaSr,gammaRa,gammaSO4] = gamma_sit(IS, Na_conc, Cl_conc);
    end
    
    % calculate Kd values based on calculated gamma values
    [Kd_RainBa,Kd_SrinBa,Kd_RainSr,Kd_BainSr] = Kd_coeff(gammaBa, gammaSr, gammaRa);
    
    % Set the change in solids as zero initially;
    % these values might be changed latter based on specific conditions;
    % if nothing happens to a particular value, they will be unchanged.
    delta_BaSO4 = 0; delta_SrSO4 = 0; 
    delta_RainBa = 0; delta_RainSr = 0; delta_BainSr = 0; delta_SrinBa = 0;
    %% Determine situations based on K values
    % calculate the solubility product of barite and celestite
    K_barite = Ba_conc*gammaBa*SO4_conc*gammaSO4;
    K_celestite = Sr_conc*gammaSr*SO4_conc*gammaSO4;
    if K_barite > Ksp_barite && K_celestite > Ksp_celestite
        % Both barite & celestite precipitate
        K_compare = 1;
    elseif K_barite > Ksp_barite && K_celestite <= Ksp_celestite
        % Only barite precipitate
        K_compare = 2;
    elseif K_barite <= Ksp_barite && K_celestite > Ksp_celestite
        % Only celestite precipitate
        K_compare = 3;
    elseif K_barite <= Ksp_barite && K_celestite <= Ksp_celestite
        % No precipitation
        K_compare = 4;
    end
    
    switch K_compare
        case 1
            lb = [0,0,0,0,0,0,0];
            % choose the smallest upper boundary based on the relative amount of Ba, Sr, and SO4.
            if Ba < SO4 && Sr < SO4
                ub = [Ba,Sr,SO4,Sr,Ba,Ba,Sr];
            elseif SO4 <= Ba && Sr < SO4
                ub = [Ba,Sr,SO4,Sr,Ba,SO4,Sr];
            elseif SO4 > Ba && Sr >= SO4
                ub = [Ba,Sr,SO4,Sr,Ba,Ba,SO4];
            elseif SO4 <= Ba && SO4 <=Sr
                ub = [Ba,Sr,SO4,Sr,Ba,SO4,SO4];
            end
            x0 = [Ba,Sr,SO4,x0_SrinBa,x0_SrinBa,x0_BaSO4,x0_BaSO4];
            rng default
            func = @(x) equation_sets(x,iniVolume,perVolume,drips,gammaBa,gammaSr,gammaSO4,Ba,Sr,SO4,Kd_SrinBa,Kd_BainSr,Ksp_barite,Ksp_celestite);
            x = lsqnonlin(func,x0,lb,ub,options);
            % x contains results for Ba, Sr, SO4, BaSO4, and SrSO4 after precipitation
            % where x(1) = Ba total in solution (mmole),
            % x(2) = Sr total in solution (mmole),
            % x(3) = total SO4 in solution (mmole),
            % x(4) = the newly co-precipitated delta_SrinBa in this drip (mmole)
            % x(5) = the newly co-precipitated delta_BainSr in this drip (mmole)
            % x(6) = the newly precipitated BaSO4 (delta_BaSO4) in this drip (mmole).
            % x(7) = the newly precipitated SrSO4 (delta_SrSO4) in this drip (mmole).
            delta_SrinBa = x(4); delta_BainSr = x(5);
            delta_BaSO4 = x(6); delta_SrSO4 = x(7);
            if delta_SrinBa >= delta_BaSO4
                %% Change condition to case 3 when Sr replaced all Ba in BaSO4 & exceeded the total solid-SO4
                lb = [0,0,0,0,0];
                if Sr<SO4
                    ub = [Ba,Sr,SO4,Ba,Sr];
                elseif SO4 <= Sr
                    ub = [Ba,Sr,SO4,Ba,SO4];
                end
                x0 = [Ba,Sr,SO4,x0_BainSr,x0_SrSO4];
                rng default
                func = @(x) equ_noBarite(x,iniVolume,perVolume,drips,gammaSr,gammaSO4,Ba,Sr,SO4,Kd_BainSr,Ksp_celestite);
                x = lsqnonlin(func,x0,lb,ub,options);
                % where x(1) = Ba total in solution (mmole),
                % x(2) = Sr total in solution (mmole),
                % x(3) = total SO4 in solution (mmole),
                % x(4) = the newly co-precipitated delta_BainSr in this drip (mmole)
                % x(5) = the newly precipitated SrSO4 (delta_SrSO4) in this drip (mmole).
                Ba = x(1); SO4 = x(3); delta_BainSr = x(4); delta_SrSO4 = x(5);
                % calculate the change of amounts in solid barite and celestite using Kd value
                delta_RainSr = Ra*delta_SrSO4*Kd_RainSr/Sr;
                if Ra > delta_RainSr && Ra > 0
                    Sr = x(2)+ delta_RainSr; 
                    Ra = Ra - delta_RainSr; 
                elseif Ra > 0 && Ra <= delta_RainSr
                    Sr = x(2)+ Ra; 
                    Ra = 0;
                else
                    Sr = x(2); Ra = 0;
                end
                delta_RainBa = 0; delta_SrinBa = 0; delta_BaSO4 = 0;
            else
                % calculate the change of amounts in solid barite and celestite using Kd_RainBa & Kd_RainSr
                % then determine the final value of Ba, Sr, BaSO4, SrSO4, and Ra in this drip
                SO4 = x(3); 
                delta_RainBa = Ra*delta_BaSO4*Kd_RainBa/Ba; delta_RainSr = Ra*delta_SrSO4*Kd_RainSr/Sr;
                if Ra > (delta_RainBa + delta_RainSr) && Ra > 0
                    Ba = x(1)+ delta_RainBa; Sr = x(2)+delta_RainSr; 
                    Ra = Ra - delta_RainBa - delta_RainSr;
                elseif Ra > 0 && Ra <= (delta_RainBa + delta_RainSr)
                    temp_RainBa = delta_RainBa; temp_RainSr = delta_RainSr;
                    delta_RainBa = Ra*temp_RainBa/(temp_RainBa + temp_RainSr);
                    delta_RainSr = Ra*temp_RainSr/(temp_RainBa + temp_RainSr);
                    Ba = x(1)+delta_RainBa; 
                    Sr = x(2)+delta_RainSr; 
                    Ra = 0;
                else
                    Ba = x(1); Sr = x(2); Ra = 0;
                    delta_RainBa = 0; delta_RainSr = 0;
                end
            end
            flag_BaSO4 = delta_BaSO4; flag_SrSO4 = delta_SrSO4;
            % In this condition, only precipitation happens, so flag as the change in barite and celestite in this step.
            % Note that later on, this flag may be changed due to dissolution
        case 2
            lb = [0,0,0,0,0]; % set lower boundary
            if Ba < SO4
                ub = [Ba,Sr,SO4,Sr,Ba]; % set upper boundary
            elseif SO4 <= Ba
                ub = [Ba,Sr,SO4,Sr,SO4];
            end
            x0 = [Ba,Sr,SO4,x0_SrinBa,x0_BaSO4];
            rng default
            func = @(x) equ_noCeles(x,iniVolume,perVolume,drips,gammaBa,gammaSO4,Ba,Sr,SO4,Kd_SrinBa,Ksp_barite);
            x = lsqnonlin(func,x0,lb,ub,options);
            % where x(1) = Ba total in solution (mmole),
            % x(2) = total Sr in solution (mmole),
            % x(3) = total SO4 in solution (mmole),
            % x(4) = the newly co-precipitated delta_SrinBa in this drip (mmole).
            % x(5) = the newly precipitated BaSO4 (delta_BaSO4) in this drip (mmole).
            Sr = x(2); SO4 = x(3); delta_SrinBa = x(4); delta_BaSO4 = x(5);
            % calculate the change of amounts in solid barite and celestite using Kd value
            delta_RainBa = Ra*x(5)*Kd_RainBa/Ba;
            if Ra > delta_RainBa && Ra > 0
                Ba = x(1) + delta_RainBa; 
                Ra = Ra - delta_RainBa;
            elseif Ra > 0 && Ra <= delta_RainBa
                Ba = x(1)+Ra ; 
                Ra = 0;
            else
                Ba = x(1); Ra = 0;
            end
            delta_RainSr = 0; delta_BainSr = 0; delta_SrSO4 = 0;
            flag_BaSO4 = delta_BaSO4;
        case 3
            % calculate the total celestite will precipitate
            lb = [0,0,0,0,0];
            if Sr<SO4
                ub = [Ba,Sr,SO4,Ba,Sr];
            elseif SO4 <= Sr
                ub = [Ba,Sr,SO4,Ba,SO4];
            end
            x0 = [Ba,Sr,SO4,x0_BainSr,x0_SrSO4];
            rng default
            func = @(x) equ_noBarite(x,iniVolume,perVolume,drips,gammaSr,gammaSO4,Ba,Sr,SO4,Kd_BainSr,Ksp_celestite);
            x = lsqnonlin(func,x0,lb,ub,options);
            % where x(1) = Ba total in solution (mmole),
            % x(2) = Sr total in solution (mmole),
            % x(3) = total SO4 in solution (mmole),
            % x(4) = the newly co-precipitated delta_BainSr in this drip (mmole)
            % x(5) = the newly precipitated SrSO4 (delta_SrSO4) in this drip (mmole).
            Ba = x(1); SO4 = x(3); delta_BainSr = x(4); delta_SrSO4 = x(5);
            % calculate the change of amounts in solid barite and celestite using Kd value
            delta_RainSr = Ra*delta_SrSO4*Kd_RainSr/Sr;
            if Ra > delta_RainSr && Ra > 0
                Sr = x(2)+ delta_RainSr; 
                Ra = Ra - delta_RainSr; 
            elseif Ra > 0 && Ra <= delta_RainSr
                Sr = x(2)+ Ra; 
                Ra = 0;
            else
                Sr = x(2); Ra = 0;
            end
            delta_RainBa = 0; delta_SrinBa = 0; delta_BaSO4 = 0;    
            flag_SrSO4 = delta_SrSO4;
        case 4
            delta_BaSO4 = 0; delta_SrSO4 = 0;
            delta_RainBa = 0; delta_SrinBa = 0; delta_RainSr = 0; delta_BainSr = 0; 
            flag_SrSO4 = 0; flag_BaSO4 = 0;
    end
    
    %% Quality control - check K value of celestite after precipitation
    Ba_conc = Ba/(iniVolume + perVolume*drips);
    Sr_conc = Sr/(iniVolume + perVolume*drips);
    Ra_conc = Ra/(iniVolume + perVolume*drips);
    SO4_conc = SO4/(iniVolume + perVolume*drips);
    IS = 0.5*(4*(Ba_conc+Sr_conc+SO4_conc+Ca_conc+Ra_conc) + Na_conc + Cl_conc); % calculate ionic strength
    if med == 1
        gammaBa = gamma_fresh(IS);
        gammaSr = gammaBa; gammaSO4 = gammaBa; gammaRa = gammaBa;
    elseif med == 2
        [gammaBa,gammaSr,gammaRa,gammaSO4] = gamma_pitzer(IS,Na_conc,Ca_conc,Cl_conc,Ba_conc,Sr_conc,SO4_conc);
    elseif med == 3
        [gammaBa,gammaSr,gammaRa,gammaSO4] = gamma_sit(IS, Na_conc, Cl_conc);
    end
    % calculate K value based on calculated gamma values
    K_celestite = Sr_conc*gammaSr*SO4_conc*gammaSO4;
    K_barite = Ba_conc*gammaBa*SO4_conc*gammaSO4;
    if K_compare == 1
        if K_celestite > (1.0001*Ksp_celestite) || K_barite > (1.0001*Ksp_barite)
            K_check = 1;
            branch_unfinished = 111;
        else
            K_check = 0;
        end
    elseif K_compare == 2
        if K_barite > (1.0001*Ksp_barite)
            K_check = 1;
            branch_unfinished = 222;
        else
            K_check = 0;
        end
    elseif K_compare == 3
        if K_celestite > (1.0001*Ksp_celestite)
            K_check = 1;
            branch_unfinished = 333;
        else
            K_check = 0;
        end
    elseif K_compare == 4
        K_check = 0;
    end
    
    %% Decide how to add results into the huge matrix
    sum_BaSO4 = sum_BaSO4 + delta_BaSO4; sum_SrSO4 = sum_SrSO4 + delta_SrSO4;
    sum_RainBa = sum_RainBa + delta_RainBa; sum_RainSr = sum_RainSr + delta_RainSr; 
    sum_BainSr = sum_BainSr + delta_BainSr; sum_SrinBa = sum_SrinBa + delta_SrinBa;
    % Delta values are updated depending on whether iteration happens 
    % if iteration happens, delta values will be the sum of this step and
    % previous step; if not, delta values will just be what is calculated.
    if K_last == 1
        delta_BaSO4 = delta_BaSO4 + step_BaSO4; delta_SrSO4 = delta_SrSO4 + step_SrSO4;
        delta_RainBa = delta_RainBa + step_RainBa; delta_RainSr = delta_RainSr + step_RainSr;
        delta_SrinBa = delta_SrinBa + step_SrinBa; delta_BainSr = delta_BainSr + step_BainSr;
    end
    % Record the changes in solids in this step, and K_check value, in case
    % they need to be used to compensate for quality control in the next step.
    step_BaSO4 = delta_BaSO4; step_SrSO4 = delta_SrSO4;
    step_BainSr = delta_BainSr; step_SrinBa = delta_SrinBa;
    step_RainSr = delta_RainSr; step_RainBa = delta_RainBa;
    K_last = K_check;
    % Update the initial values for solving equation sets.
    x0_BaSO4 = delta_BaSO4; x0_SrSO4 = delta_SrSO4; x0_BainSr = delta_BainSr; x0_SrinBa = delta_SrinBa;
    
    %% Update the matrix
    results(drips,1) = drips; results(drips,2) = Ba; results(drips,3) = Sr; results(drips,4) = Ra;
    results(drips,5) = SO4; results(drips,6) = Na_conc; results(drips,7) = Ca_conc; results(drips,8) = Cl_conc;
    results(drips,9) = gammaBa; results(drips,10) = gammaSr; results(drips,11) = gammaRa; results(drips,12) = gammaSO4;
    results(drips,13) = Kd_RainBa; results(drips,14) = Kd_RainSr; results(drips,15) = Kd_BainSr; results(drips,16) = Kd_SrinBa;
    results(drips,17) = delta_RainBa; results(drips,18) = delta_RainSr; results(drips,19) = delta_BainSr; results(drips,20) = delta_SrinBa;
    results(drips,21) = delta_BaSO4; results(drips,22) = delta_SrSO4;
    results(drips,23) = K_check; results(drips,24) = branch_unfinished; results(drips,25) = K_compare;
    results(drips,26) = sum_BaSO4; results(drips,27) = sum_SrSO4;
    results(drips,28) = sum_RainBa; results(drips,29) = sum_RainSr; 
    results(drips,30) = sum_BainSr; results(drips,31) = sum_SrinBa;
    
    % Decide if the calculation within the drip need to be iterated
    if K_last == 1
        drips = drips - 1;
    end
end

% merge into a table at last:
% now remember which column represents what:
colNames = {'drips','Ba_tot','Sr_tot','Ra_tot','SO4_tot','Na_conc','Ca_conc','Cl_conc','gamma_Ba','gamma_Sr','gamma_Ra','gamma_SO4','Kd_RaBa','Kd_RaSr','Kd_BaSr','Kd_SrBa','delta_RainBa','delta_RainSr','delta_BainSr','delta_SrinBa','delta_BaSO4','delta_SrSO4','K_check','branch_unfinished','Kd_compare','sum_BaSO4','sum_SrSO4','sum_RainBa','sum_RainSr','sum_BainSr','sum_SrinBa'};
sTable = array2table(results,'VariableNames',colNames);



% writetable(sTable,'titration_noDiss_1202.xls','Sheet',1)