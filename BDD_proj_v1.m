%**************************************************************************
%                           Jason Hunter
%                    Blade design development
%                         main program
%  This program will determine target angle of attack and CL values for
% a proposed wind turbine based on selected airfoil data file, as well as
% compute and plot the distribution of chord and twist along the blade
% radius, and use a linear fit to fine tune the overall blade design, with a
% final pass completed by using a weighted blending of ideal and linear fits for
% the chord.
%
%   program based on examples provided by   Andrew J. Goupee
%                                           Prof. Mechanical Engineering
%                                           University of Maine
%
%   last modified: 6/15/23
%   by: Jason Hunter
%
%**************************************************************************
clear
clc
close all

%import airfoil data
fx6 = readtable("fx66196v_foil.txt");
fx6.Properties.VariableNames = {'alpha','CL','CD'};
fp = readtable("fx66196v.DAT"); %for profile plot
fp.Properties.VariableNames = {'xc','yc'};
fpsuc = table(fp.xc(1:44),fp.yc(1:44));
fpsuc.Properties.VariableNames = {'xc','yc'};
fppre = table(fp.xc(45:end),fp.yc(45:end));
fppre.Properties.VariableNames = {'xc','yc'}; 

%calculate lift/drag ratios
ii = find(fx6.alpha >= 0 & fx6.alpha <= 20);
alpha = fx6.alpha(ii);
CL = fx6.CL(ii);
CLCD = CL ./ fx6.CD(ii);

jj = find(CLCD == max(CLCD));


figure(1)
clf
box on
plot(alpha,CLCD)
xlabel('\alpha (deg)')
ylabel('C_L/C_D')
title('C_L/C_D as a function of \alpha')
text(10,160,['\alpha_t = ' num2str(alpha(jj),'%.2f') char(176) ])
text(10,140,['C_L @ \alpha_t = ' num2str(CL(jj),'%.3f')])

rhub = 3.5; %Hub radius in m
rtip = 110.5; %Tip radius in m
B = 3; %Number of blades
TSR = 8; %Tiip speed ratio

%Create radial stations
r = rhub:107/45:rtip;

%Compute local inflow angle phi
phi = atand(2./(3*TSR*(r/rtip)));

%Compute beta
beta = phi-alpha(jj);

%Compute chord
chord = rtip*(8*pi*sind(phi)/(3*B*CL(jj)*TSR));

%Plot result
figure(2)
clf
subplot(2,1,1), plot(r,beta), xlabel('Radius (m)'), ylabel('\beta (deg)')
subplot(2,1,2), plot(r,chord), xlabel('Radius (m)'), ylabel('chord (m)')

%Find indices to fit
kk = find(r >= rhub+0.5*(rtip-rhub));

%Find linear fits for outer half of blade
P1 = polyfit(r(kk),beta(kk),1);
P2 = polyfit(r(kk),chord(kk),1);

%Compute new linear distributions
betalin = polyval(P1,r);
chordlin = polyval(P2,r);

%Plot
figure(3)
clf
subplot(2,1,1), hold on, box on, plot(r,beta,r,betalin), ...
xlabel('Radius (m)'), ylabel('\beta (deg)'), ...
legend('Ideal','Linear Fit')
title('Ideal vs Linear Fit distributions')
subtitle('\beta top and chord length bottom')
subplot(2,1,2), hold on, box on, plot(r,chord,r,chordlin), ...
xlabel('Radius (m)'), ylabel('chord (m)')

%exported to excel for some hand manipulation
%chordfile = table(r,chord,chordlin);
%betafile = table(r,beta,betalin);

%design fit, messy need to fix comments
iii = length(r);
rd = 1:1:length(r);
dr = (rtip-rhub) / (iii);
ri = rhub + (rd-1)*dr + dr/2;
root = find(ri <= rtip*.05 + rhub);
blend = find(ri > rtip*.05 + rhub & ri <= rtip*.5 + rhub);
linfit = find(ri > rtip*.5);

chordD(root) = 6;
chordD(blend) = chord(blend).*(exp(-2.5)) + chordlin(blend);
chordD(linfit) = chordlin(linfit);
airfoil = zeros(1,length(chordD),'double')+2;
airfoil(root) = 1;
drtable = zeros(1,length(chordD),'double')+dr;

%setup and export blade input file for BEM code
Desfile = table(ri',betalin',drtable',chordD',airfoil','VariableNames',["rl","twist","dr","chordD","af"]);
writetable(Desfile,'bladeD.txt','Delimiter',' ');

%plot design vs ideal vs linear
figure(4)
clf
hold on
plot(r,chord,r,chordlin,r,chordD)
xlabel("Radius (m)")
ylabel("Chord (m)")
legend('Ideal','Linear fit','Design')
title('Chord length distribution')
subtitle('Ideal vs Linear Fit vs Proposed Design')

figure(5)
clf
hold on
axis equal
plot(fpsuc.xc,fpsuc.yc,fppre.xc,fppre.yc)
xlabel('x/c (m)')
ylabel('y/c (m)')
title('FX 66-S-196 v1 Airfoil Profile')
