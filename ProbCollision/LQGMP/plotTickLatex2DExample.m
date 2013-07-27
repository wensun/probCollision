% 
% Copyright (C) 2011 Alex Bikfalvi
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.

% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
%

%%

t = linspace(0,2,1000);
x = trimf(t,[0 1 2]);

%%

figure(1); clf;
set(gcf,'Color',[1 1 1]);

subplot(3,1,1);
plot(t,x);

n = length(t);
T = 2;
ts = T/n;
fs = 1/ts;

f = linspace(0,fs,n);
y = fft(x)./n;

xlabel('$t$ [s]');
ylabel('$x(t)$');

set(gca,'YTick',[0 0.5 1]);
set(gca,'YTickLabel',{'0';'Xm/2';'Xm'});

legend('Time-domain $x(t)$');
legend boxoff;

subplot(3,1,2);

stem(f, abs(y),'red', 'MarkerSize',4,'MarkerFaceColor','red');

xlim([0 10]);
ylim([0 0.5]);

xlabel('$f$ $[\mbox{s}^{-1}]$');
ylabel('$|\mathcal{F}\{x\}(f)|$');

set(gca,'YTick',[0 0.25 0.5]);
set(gca,'YTickLabel',{'0';'Xm/2';'Xm/4'});

legend('Fourier transform: magnitude');
legend boxoff;

subplot(3,1,3);

stem(f, angle(y),'green', 'MarkerSize',4,'MarkerFaceColor','green');

xlim([0 10]);
ylim([-pi pi]);

set(gca,'YTick',[-pi -pi/2 0 pi/2 pi]);
set(gca,'YTickLabel',{'-pi';'-pi/2';'0';'pi/2';'pi'});

xlabel('$f$ $[\mbox{s}^{-1}]$');
ylabel('$\arg \mathcal{F}\{x\}(f)$');

legend('Fourier transform: phase');
legend boxoff;

%%

figure(2); clf;
set(gcf,'Color',[1 1 1]);

subplot(3,1,1);
plot(t,x);

n = length(t);
T = 2;
ts = T/n;
fs = 1/ts;

f = linspace(0,fs,n);
y = fft(x)./n;

xlabel('$t$ [s]');
ylabel('$x(t)$');

set(gca,'XTick',0:0.2:2);
set(gca,'YTick',[0 0.5 1]);
set(gca,'YTickLabel',{'$0$';'$\displaystyle\frac{X_m}{2}$';'$X_m$'});

legend('Time-domain $x(t)$');
legend boxoff;

plotTickLatex2D('xlabeldy',0.02);

subplot(3,1,2);

stem(f, abs(y),'red', 'MarkerSize',4,'MarkerFaceColor','red');

xlim([0 10]);
ylim([0 0.5]);

xlabel('$f$ $[\mbox{s}^{-1}]$');
ylabel('$|\mathcal{F}\{x\}(f)|$');

set(gca,'YTick',[0 0.25 0.5]);
set(gca,'YTickLabel',{'$0$';'$\displaystyle\frac{X_m}{4}$';'$\displaystyle\frac{X_m}{2}$'});

legend('Fourier transform: magnitude');
legend boxoff;

plotTickLatex2D('xlabeldy',0.02);

subplot(3,1,3);

stem(f, angle(y),'green', 'MarkerSize',4,'MarkerFaceColor','green');

xlim([0 10]);
ylim([-pi pi]);

set(gca,'YTick',[-pi -pi/2 0 pi/2 pi]);
set(gca,'YTickLabel',{'$-\pi$';'$-\displaystyle\frac{\pi}{2}$';'$0$';'$\displaystyle\frac{\pi}{2}$';'$\pi$'});

xlabel('$f$ $[\mbox{s}^{-1}]$');
ylabel('$\arg \mathcal{F}\{x\}(f)$');

legend('Fourier transform: phase');
legend boxoff;

plotTickLatex2D('xlabeldy',0.02);
