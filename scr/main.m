%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HJUSTROM CURVE FOR SEDIMENT DEPOSITION AND EROSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize code
clear all; close all; clc; tic

% set friction coefficient
Cf = [0.001 0.002 0.005 0.01 0.02 0.05 0.1];

%% COMPUTING THRESHOLD VELOCITIES

% dimensional sediment diameter
dgmin = 0.0001;
dgmax = 0.10;
dg = dgmin:dgmin:dgmax;
n = length(dg);

% initialize threshold velocity arrays
wg = zeros(n,1);
uc = zeros(n,1);

% loopo ver sediment diameters
for j = 1:length(Cf)
for i = 1:n
    
% lower threshold for deposition: settling velocity
    [wg(i),~] = settling_velocity(dg(i), 1);
    
% upper threshold for erosion: critical Shield number
    [thetac,dstar] = critical_shield(dg(i), 1);
    uc(i,j) = (8*9.806*thetac*dg(i)/Cf(j))^0.5;
    
% end loop over diameters
end
end

%% AREAS OF DEPOSITION, TRANSPORTATION, AND EROSION

% check array shape
if iscolumn(dg)==0
    dg=dg';
end
if iscolumn(wg)==0
    wg=wg';
end
if iscolumn(uc)==0
    uc=uc';
end

% initialize storage arrays for plotting
area_DP = zeros(2*n,2);
area_TR = zeros(2*n,2);
area_ER = zeros(2*n,2);
area_ERaux = zeros(2*n,2);

% deposition area
for i = 1:n
    area_DP(i,1) = dg(i);
    area_DP(i,2) = 10^(-12);
end
j = 0;
for i = n+1:2*n
    j = j + 1;
    area_DP(i,1) = dg(n-j+1);
    area_DP(i,2) = wg(n-j+1);
end

% transportation area
for i = 1:n
    area_TR(i,1) = dg(i);
    area_TR(i,2) = wg(i);
end
j = 0;
for i = n+1:2*n
    j = j + 1;
    area_TR(i,1) = dg(n-j+1);
    area_TR(i,2) = uc(length(Cf),n-j+1);
end

% erosion area
for i = 1:n
    area_ER(i,1) = dg(i);
    area_ER(i,2) = uc(1,i);
end
j = 0;
for i = n+1:2*n
    j = j + 1;
    area_ER(i,1) = dg(n-j+1);
    area_ER(i,2) = 10000;
end

% transitional erosion area
for i = 1:n
    area_ERaux(i,1) = dg(i);
    area_ERaux(i,2) = uc(length(Cf),i);
end
j = 0;
for i = n+1:2*n
    j = j + 1;
    area_ERaux(i,1) = dg(n-j+1);
    area_ERaux(i,2) = uc(1,n-j+1);
end

ymin = 10^(-5);
ymax = 10;

for i = 1:2*n
    if area_DP(i,2)<ymin 
        area_DP(i,2) = ymin;
    elseif area_DP(i,2)>ymax
        area_DP(i,2) = ymax;
    end
    if area_TR(i,2)<ymin 
        area_TR(i,2) = ymin;
    elseif area_TR(i,2)>ymax
        area_TR(i,2) = ymax;
    end
    if area_ER(i,2)<ymin 
        area_ER(i,2) = ymin;
    elseif area_ER(i,2)>ymax
        area_ER(i,2) = ymax;
    end
    if area_ERaux(i,2)<ymin 
        area_ERaux(i,2) = ymin;
    elseif area_ERaux(i,2)>ymax
        area_ERaux(i,2) = ymax;
    end
end

%% PLOTTING

font = 'times new roman';
tickfs = 20;
axisfs = 24;
titlefs = 30;

% new figure
figure

% plot threshold lines
thr1 = plot(dg*1000,wg,'k','linewidth',1.2);
hold on
thr2 = plot(dg*1000,uc,'k','linewidth',1.2);

% plot areas
dp = fill(area_DP(:,1)*1000,area_DP(:,2),[255 178 88]/255,'edgecolor','none');
tr = fill(area_TR(:,1)*1000,area_TR(:,2),[203 213 233]/255,'edgecolor','none');
er = fill(area_ER(:,1)*1000,area_ER(:,2),[239 255 123]/255,'edgecolor','none');
er_aux = fill(area_ERaux(:,1)*1000,area_ERaux(:,2),[221 234 178]/255,'edgecolor','none');
hold off

% stack order
uistack(thr1,'top')
uistack(thr1,'top')
uistack(dp,'bottom')
uistack(tr,'bottom')
uistack(er,'bottom')
uistack(er_aux,'bottom')

% axis scale and font
% xxt = 10.^[log10(dgmin*1000):1:log10(dgmax*1000)];
xxt = [0.1 0.2 0.5 1 2 5 10 20 50 100];

set(gca,'xscale','log','yscale','log','fontsize',tickfs,'fontname',font, ...
    'xtick',xxt, 'xticklabel', {'0.1','0.2','0.5','1','2','5','10','20','50','100'}, ...
    'yminortick','on',...
    'layer','top');
ylim([ymin ymax])

% area setting
grid on
box on

% labels
xlabel('d_{g} (mm)','fontsize',axisfs,'fontname',font)
ylabel('U (m/s)','fontsize',axisfs,'fontname',font)
title('Hjulström curves','fontsize',titlefs,'fontname',font)

%% end of code
toc