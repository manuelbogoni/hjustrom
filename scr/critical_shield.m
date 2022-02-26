%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CRITICAL SHIELD NUMBER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dg : dimensional grain size (m)
% fl : switch for formula (1: Bronwlie 1981, 2: Van Rijn 1989)
% thetac : critical Shield number for incipient erosion / bedload
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [thetac,dstar] = critical_shield(dg, fl)

% kinematic vistosity of the water
nu = 1.01*10^(-6);

% sediment density to water density ratio
R = 1.65;

% gravity
g = 9.806;

% Reynolds particle number
Rp = (R*g)^0.5/(nu)*dg^1.5;

% dimensionless grain size
dstar = Rp^(2/3);

switch fl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Brownlie(1981)
case 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    thetac = 0.22*Rp^(-0.6) + 0.06*exp(-17.77*Rp^(-0.6));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Van Rijn (1989)
case 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if dstar <= 4
        thetac = 0.24 / dstar;
    elseif dstar > 4 && dstar <= 10
        thetac = 0.14 * dstar^(-0.64);
    elseif dstar > 10 && dstar <= 20
        thetac = 0.04 * dstar^(-0.10);
    elseif dstar > 20 && dstar <= 150
        thetac = 0.013 * dstar^0.29;
    elseif dstar > 150
        thetac = 0.055;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error
    otherwise
        error('Critical Shield number: wrong input switch')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end of function
return