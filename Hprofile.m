function hi = Hprofile(h0, h, Tz, xi)
% CUMCM 2020 Problem A: The Furnace Temperature Curve
% zhou lvwen: zhou.lv.wen@gmail.com
% Wechat Official ID: MATHmodels 
% September 11, 2020

if nargin==0
   % Temperature of the workshop and small temperature zones
   Tz = [ 175  175  175  175  175  195  235  255  255   25   25];
   % Heat transfer coef. of front&back zone, zones 1-5, 6, 7, 8-9, 10-11.
   h0 = 0.0074; h = [0.0196 0.0214 0.0311 0.0197 0.0109 4.6050];  xi = [];  % [W/m^2-K]
end

% length of front zone, back zone, small temperature zone and gap.
Lf = 25.0; Lb = 25.0; Lz = 30.5; Lg =  5.0;  % [cm]
nz = length(Tz);                 % number of small temperature zones (STZs)

L  = Lf + Lz*nz + Lg*(nz-1) + Lb;% Total length of the Oven
xz = Lb + (Lz+Lg)*[0:1:nz-1];    % Starting point of STZs

% h (heat transfer coef.) profile
hz = [h(1) h(1) h(1) h(1) h(1) h(2) h(3) h(4) h(4) h(5) h(5)];% [W/m^2-K]
x = [0, Lb, Lz+Lg/2, repmat(Lz+Lg,1,nz-2), Lz+Lg/2, Lb];
x = cumsum(x); 
h = [h0 hz h0];

i = max(ceil(interp1(x, 0:length(h), xi)),1);
hi = h(i);

% -------------------------------------------------------------------------

% plot h (heat transfer coefficient)
if nargin~=0; return; end
reflowoven(Tz, [0.0025, 0.0375], 5e-3); hold on
xi = 0:0.1:435.5;
i = max(ceil(interp1(x, 0:length(h), xi, 'linear')),1);
plot(xi,h(i), 'linewidth', 2); 
xlabel('x (cm)'); ylabel('h')