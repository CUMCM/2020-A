function Ti = Tprofile(Tz, T0, xi)
% CUMCM 2020 Problem A: The Furnace Temperature Curve
% zhou lvwen: zhou.lv.wen@gmail.com
% Wechat Official ID: MATHmodels 
% September 11, 2020

if nargin<=1; T0 = 25; end
if nargin==0; Tz = [175 175 175 175 175 195 235 255 255 25 25]; xi =[]; end

% length of front zone, back zone, small temperature zone and gap.
[Lf, Lb, Lz, Lg, nz, L] = reflowoven();
xz = Lb + (Lz+Lg)*[0:1:nz-1];    % Starting point of STZs

% Temperature profile
x = [0, Lb, repmat([Lz,Lg],1,nz-3), Lz, Lg+Lz+Lg/2, Lg/2+Lz, Lb];
x = cumsum(x);
T = [T0,Tz(reshape([1:nz-1; 1:nz-1], 1,2*nz-2)), T0];
Ti = interp1(x, T, xi);

% -------------------------------------------------------------------------

% plot temperature profile
if nargin~=0; return; end
reflowoven(Tz, [10, 280], 20); hold on
plot(x,T,'+-', 'linewidth', 2, 'markersize', 6); 
xlabel('x (cm)'); ylabel('T (^\circ C)')

% load and plot experiment data
v = 70/60; 
dat = load('expt.dat'); texpt = dat(:,1); Texpt = dat(:,2);
plot(v*texpt, Texpt, 'b', 'linewidth',2)
