function [is, area, hplot] = prolim(v, t, T, hplot)
% CUMCM 2020 Problem A: The Furnace Temperature Curve
% zhou lvwen: zhou.lv.wen@gmail.com
% Wechat Official ID: MATHmodels 
% September 11, 2020

if nargin<=3; hplot = []; end

x = v/60*t;
dT = [0, diff(T)./diff(t)];
id = find(dT>0, 1, 'last');

% Peak temperature
[Tmax, imax] = max(T);
isTmax = Tmax>=240 & Tmax<=250;

% Maximum slope of temperature rising & Minimum slope of temperature drop
[dTmax, idmax] = max(dT); [dTmin, idmin] = min(dT); 
isdT = (dTmax<= 3) & (dTmin>=-3) & all(dT(1:imax-5)>=0) & all(dT(imax+5:end)<=0);

% Time for temperature rising from 150 to 190
i150 = find(T(1:id)>150, 1, 'first'); i190 = find(T(1:id)<190, 1, 'last');
t150_190 = t(i190)-t(i150);
ist150_190 = (t150_190>=60) & (t150_190<=120);

% Time of temperature greater than 217
i217 = find(T>=217);
if isempty(i217)
    tg217 = 0;
    istg217 = false;
    area = 0;
else
    tg217 = t(i217(end))-t(i217(1));
    istg217 = (tg217>=40) & (tg217<=90);
    xpoly = t([i217(1):imax, imax,    i217(1)]);
    ypoly = T([i217(1):imax, i217(1), i217(1)]);
    area = polyarea(xpoly, ypoly);
end
is = isdT + ist150_190 + istg217 + isTmax + (v>=65&v<=100);

% -------------------------------------------------------------------------

if nargin<=3; return; end

if ishandle(hplot)
    hplot(1).XData = x([idmax, idmin]); hplot(1).YData = T([idmax, idmin]);
    hplot(2).String = sprintf('dT/dt_{max} = %5.2f, dT/dt_{min} = %5.2f', dTmax, dTmin);
    hplot(3).XData = x([i150:i190, i190]); hplot(3).YData = T([i150:i190, i150]);
    hplot(4).String = sprintf('t_{150-190} = %5.2f', t150_190);
    hplot(5).XData = x(i217); hplot(5).YData = T(i217);
    hplot(6).String = sprintf('t_{>217} = %5.2f', tg217);
    hplot(7).Position(1) = x(i217(1))+10; hplot(7).String = num2str(area);    
    hplot(8).String = sprintf('T_{max} = %5.2f', Tmax);
    hplot(9).String=  sprintf('v = %d cm/min', v);
    hplot(10).XData = xpoly*v/60; hplot(10).YData = ypoly;
    hplot(11).XData = x(imax); hplot(11).YData = Tmax;
    hplot(12).XData = x; hplot(12).YData = T;
else
    clear hplot; hold on
    hplot(1) = plot(x([idmax, idmin]), T([idmax, idmin]), 'ro', 'linewidth',2);
    hplot(2) = text(200, 130, sprintf('dT/dt_{max} = %5.2f, dT/dt_{min} = %5.2f', dTmax, dTmin));
    hplot(3) = fill(x([i150:i190, i190]), T([i150:i190, i150]), [1,0.8,0.8]);
    hplot(4) = text(200,110, sprintf('t_{150-190} = %5.2f', t150_190));
    hplot(5) = fill(x(i217), T(i217), [0.8,1,0.8]);
    hplot(6) = text(200,90, sprintf('t_{>217} = %5.2f', tg217));
    hplot(7)= text((x(i217(1)))+10, 200, num2str(area));
    hplot(7).Color = [0,0.5,0];
    hplot(8) = text(200, 70, sprintf('T_{max} = %5.2f', Tmax));
    hplot(9)= text(25,220, sprintf('v = %d cm/min', round(v)));
    hplot(10) = fill(xpoly*v/60, ypoly, 'g');
    hplot(11) = plot(x(imax), Tmax, 'ro', 'linewidth',2);
    hplot(12) = plot(x, T, '-b', 'linewidth', 2);
    set(gca, 'ytick', [0, 15, 150, 190, 217, 250]);
end

if ~isdT;       hplot(2).Color=[1,0,0]; else; hplot(2).Color = [0,0,0]; end
if ~ist150_190; hplot(4).Color=[1,0,0]; else; hplot(4).Color = [0,0,0]; end
if ~istg217;    hplot(6).Color=[1,0,0]; else; hplot(6).Color = [0,0,0]; end
if ~isTmax;     hplot(8).Color=[1,0,0]; else; hplot(8).Color = [0,0,0]; end
