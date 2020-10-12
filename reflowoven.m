function [Lf, Lb, Lz, Lg, nz, L] = reflowoven(Tz, y, dy)
% CUMCM 2020 Problem A: The Furnace Temperature Curve
% zhou lvwen: zhou.lv.wen@gmail.com
% Wechat Official ID: MATHmodels 
% September 11, 2020

% if nargin==0;
%     Tz = [175 175 175 175 175 195 235 255 255 25 25]; 
%     y = [10, 280]; dy = 20;
% end

% number of small temperature zones (STZs)
nz = 11;

% length of front zone, back zone, STZs and gap.
Lf = 25.0; Lb = 25.0; Lz = 30.5; Lg =  5.0;  % [cm]

L  = Lf + Lz*nz + Lg*(nz-1) + Lb;% Total length of the Oven
xz = Lb + (Lz+Lg)*[0:1:nz-1];    % Starting point of STZs

if nargin==0; return; end

ax(1) = axes('XAxisLocation','top', 'YAxisLocation', 'right');
axis([0,L, y+[-dy, dy]]); set(ax(1),'xtick',   xz); grid on; box off
ax(2) = axes('Color','none');
axis([0,L, y+[-dy, dy]]); set(ax(2),'xtick',xz+Lz); grid on; box off

hold on
for i = 1:2
    fill(xz+[0; 0; Lz; Lz], y(i)+dy/2*[-1; 1; 1; -1].*ones(size(xz)), Tz)
    text(xz+Lz/2, y(i)*ones(size(xz)), string(round(Tz*10)/10), ...
                         'HorizontalAlignment', 'center', 'Color', 'red')
end