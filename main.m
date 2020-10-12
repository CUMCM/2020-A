% CUMCM 2020 Problem A: The Furnace Temperature Curve
% zhou lvwen: zhou.lv.wen@gmail.com
% Wechat Official ID: MATHmodels 
% September 11, 2020

clear; clc;
%% ------------------------------------------------------------------------
figure; Tprofile(); drawnow

%% ------------------------------------------------------------------------

% length of front zone, back zone, small temperature zone and gap.
Lf = 25.0; Lb = 25.0; Lz = 30.5; Lg =  5.0;  % [cm]
nz = 11;                         % number of small temperature zones (STZs)
xz = Lb + (Lz+Lg)*[0:1:nz-1];    % Starting point of STZs
T0 = 25; 
v  = 70/60;                      % cm/s

%% -------------------------------------------------------------------------

% Tz = [175 175 175 175 175 195 235 255 255 25 25]; 
% fmin = @(hk) SimOven(T0, Tz, v, hk, 1);
% opts = optimoptions('lsqnonlin','Display','iter');
% hk0 = [0.01*ones(1,6), 5];
% [hk,resnorm,residual,~,~] = lsqnonlin(fmin,hk0,hk0*0.1,hk0*10,opts);
% h0 = hk(1); h = hk(2:end-1); k = hk(end);

hk = [0.0074  0.0196  0.0214  0.0311  0.0197  0.0109  4.6050];
figure; Hprofile();drawnow
figure; SimOven(); drawnow
%% -----------------------------Problem 1----------------------------------

Tz = [173*ones(1,5) 198 230 257*ones(1,2) 25 25];
[T, t] = SimOven(T0, Tz, v, hk);
x = v*t;

x3678 = [xz([3,6,7]) + Lz/2, xz(8)+Lz];
T3678 = interp1(x, T, x3678);
figure('position', [50,50,800,400])

reflowoven(Tz, [10, 280], 20); hold on
plot(x, T, '-b', x3678, T3678, 'ro', 'linewidth', 2)
text(x3678-30, T3678, string(T3678), 'HorizontalAlignment', 'center')


%% -----------------------------Problem 2----------------------------------  

Tz = [182*ones(1,5) 203 237 254*ones(1,2) 25 25];
figure('position', [50,450,1000,600])
reflowoven(Tz, [10, 280], 20); hold on
vmax = 100; % cm/min
[T, t] = SimOven(T0, Tz, vmax/60, hk);
[is, area, hplot] = prolim(vmax, t, T, []);
drawnow
pause(1)
while is<5
   vmax = vmax - 1;
   [T, t] = SimOven(T0, Tz, vmax/60, hk);
   [is, area, hplot] = prolim(vmax, t, T,hplot);
   drawnow
end
[is, hplot] = prolim(vmax, t, T,hplot);


%% -----------------------------Problem 3----------------------------------  

figure('position', [50,50,1000,600])
rng(0)
Tz0 = [175 195 235 255]+20*rand(1,4)-10;
v0 = (65+35*rand)/60;
SimOven(T0, [Tz0([1 1 1 1 1, 2, 3, 4 4]), 25,25], v0, hk, 2);
pause(1)
fmin = @(vTz) SimOven(T0, [vTz([1 1 1 1 1, 2, 3, 4 4]+1), 25,25], vTz(1), hk, 2);
vTz0 = [v0 Tz0];
opts = optimset('Display','iter');
vTz = fminsearch(fmin,vTz0, opts);
v = vTz(1); Tz = [vTz([1 1 1 1 1, 2, 3, 4 4]+1), 25,25];
[T, t] = SimOven(T0, Tz, v, hk);

figure('position', [1000,200,800,400])
reflowoven(Tz, [10, 280], 20); hold on
[is, area] = prolim(v*60, t, T,[]);