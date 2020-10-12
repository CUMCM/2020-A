function [R, t] = SimOven(T0, Tz, v, hk, isopt)
% CUMCM 2020 Problem A: The Furnace Temperature Curve
% zhou lvwen: zhou.lv.wen@gmail.com
% Wechat Official ID: MATHmodels 
% September 11, 2020

if nargin==0
   T0 = 25; Tz = [175 175 175 175 175 195 235 255 255 25 25];
   hk = [0.0074, 0.0196 0.0214 0.0311 0.0197 0.0109, 4.6081];
   v  = 70/60; 
end

if nargin<=4; isopt = 0; end

h0 = hk(1); h = hk(2:end-1); k = hk(end);

% length of front zone, back zone, small temperature zone and gap.
Lf = 25.0; Lb = 25.0; Lz = 30.5; Lg =  5.0;  % [cm]
nz = length(Tz);      % number of small temperature zones (STZs)
L  = Lf + Lz*nz + Lg*(nz-1) + Lb;% Total length of the Oven
nz = length(Tz);


if (nargin==0) | isopt==1
   dat = load('expt.dat'); texp = dat(:,1); Texp = dat(:,2);
end

dt = 5e-1;
t = 0:dt:L/v;
n = length(t);
x = v*t;
Ta = Tprofile(Tz, T0, x);
ha = Hprofile(h0, h, Tz, x);

nsx = 15; nst = 100;
Tg = T0*ones(nsx,1);
T  = T0*ones(1, n);

if nargin==0
    xz = Lb + (Lz+Lg)*[0:1:nz-1];    % Starting point of STZs
    
    H = 40;
    subplot(2,1,1); hold on; box on; axis image; xlim([0,L]);
    for i = [-1, 1]
        fill(xz+[0; 0; Lz; Lz], i*H+i*[0; 15; 15; 0].*ones(size(xz)), Tz)
        text(xz+Lz/2, i*H+i*15/2*ones(size(xz)), string(Tz), ...
                         'HorizontalAlignment', 'center', 'Color', 'red')
    end
    img1 = imagesc(x(1)+[-5,5], [-20,20], Tg); colorbar
    
    subplot(2,1,2); axis([0,L,0, max(Tz)+5]);hold on
    plot(v*texp, Texp, 'b', 'linewidth',2)
    ho = plot(x(1), T0, 'or', 'linewidth',2);
    xlabel('x (cm)'); ylabel('T (^\circ C)')
    img2 = imagesc(350+[-5,5], [50,130], Tg); colorbar; box on
end


for i = 2:n
    [T(i), Tg] = sensor(Ta(i), Tg, ha(i), k, dt, nsx, nst);
    if nargin==0
        set(img1, 'CData', Tg, 'XData', x(i)+[-5,5])
        set(img2, 'CData', Tg)
        set(ho, 'XData', x(i), 'YData', T(i)); colorbar
        drawnow
    end
end


if isopt==1
    Tint = interp1(t, T, texp);
    R =  Tint - Texp; 
elseif isopt==2
    clf
    dTz = sum(abs(Tz(5:8) - [175,195,235,255])<=10);
    reflowoven(Tz, [10, 280], 20);
    [is, area] = prolim(v*60, t, T, []);
    drawnow
    R = area + 2000 - 200*is + (is~=5)*2000 -250*dTz + (dTz~=4)*2000;
else
    R = T;
end

if nargin==0
    i = ceil(linspace(1,n,20));
    set(ho, 'XData', x(i), 'YData', T(i));
    legend('experiment','simulation','Location','NorthWest')
end
% -------------------------------------------------------------------------

function [Ts, T] = sensor(Ta, T, h, k, t, nx, nt)
L  =  15;
dt =  t/nt;
dx =  L/(nx-1); 
i = 2:nx-1;
for n=1:nt
    % Compute new temperature 
    T(i) = T(i) + k * (T(i+1)-2*T(i)+T(i-1))/dx^2*dt;
    % Set boundary conditions
    T(1)  = T(1)  + h*(Ta-T(1) )*dt;
    T(nx) = T(nx) + h*(Ta-T(nx))*dt;
end
Ts = T((nx+1)/2);
