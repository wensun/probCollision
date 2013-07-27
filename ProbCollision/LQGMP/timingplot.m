function timingplot()

x = [1,50,100,200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2250,2500,2750,3000];
t1 = 0.942*ones(1,size(x,2));
t2 = 0.356*x;

%FH = figure;
FH = fig('units','inches','width',6,'height',4);
axes1 = axes('Parent',FH,'YMinorTick','on','XMinorTick','on');

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 3025]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0 1000]);
hold(axes1,'all');

% Create multiple lines using matrix input to plot
plot(x,t1,'Parent',axes1,'LineWidth',3,'Color',[0 0 1],'DisplayName','Our Method');
plot(x,t2,'Parent',axes1,'LineWidth',3,'Color',[0 0.7 0],'DisplayName','Monte-Carlo');

% Create xlabel
xlabel('Number of Monte-Carlo samples','FontSize',20,'FontName','Helvetica');

% Create ylabel
ylabel('Computation Time (in ms)','FontSize',20,'FontName','Helvetica');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Position',[0.199 0.7522 0.3466 0.1875],'FontSize',20);

ti = get(gca,'TightInset');
set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

print -painters -dpdf -r600 mc-car-timing.pdf
