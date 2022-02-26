%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTLING VELOCITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dg : dimensional grain size (m)
% fl : switch for formula (1: Grace 1986, 2: Parker 1976)
% wg : dimensional settling velocity (m/s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wg,dstar] = settling_velocity(dg, fl)

% kinematic vistosity of the water
nu = 1.01*10^(-6);

% sediment density to water density ratio
R = 1.65;

% gravity
g = 9.806;

% Reynolds particle number
Rp = (R*g)^0.5/(nu)*dg^1.5;

% dimensionless grain size
dstar = dg * Rp^(2/3);

switch fl
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grace (1986) and Wilson et al. (1992)
    case 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if dstar <= 3.8
            VTS = (dstar^2)/18 - 3.1234*10^(-4)*dstar^5+1.6415*10^(-6)*dstar^8;
        elseif dstar > 3.8 && dstar <= 7.58
            VTS = 10^(-1.5446 + 2.9162*log10(dstar) - 1.0432*(log10(dstar))^2);
        elseif dstar > 7.58 && dstar <= 227
            VTS = 10^(-1.64758 + 2.94786*log10(dstar) - ...
                1.0907*(log10(dstar))^2 + 0.17129*(log10(dstar))^3);
        elseif dstar > 227
            VTS = 10^(5.1837 - 4.51034*log10(dstar) + ...
                1.687*(log10(dstar))^2 - 0.189135*(log10(dstar))^3);
        end
        
        VFL = (R*g*nu)^(1/3);
        wg = VTS*VFL*0.26;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parker (1976)
    case 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        w = 10^(-1.181 + 0.966*log10(Rp) - 0.1804*(log10(Rp))^2 + ...
          0.003746*(log10(Rp))^3 + 0.0008782*(log10(Rp))^4);
        wg = w * (R*g*dg)^0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error
    otherwise
        error('Settling velocity: wrong input switch')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end of function
return