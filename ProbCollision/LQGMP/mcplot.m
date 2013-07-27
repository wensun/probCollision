function mcplot(p)

nr = size(p,1);
x = [50,100,200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2250,2500,2750,3000];
y = zeros(nr,1);
e = zeros(nr,1);
for i = 1:nr
    y(i) = mean(p(i,:));
    e(i) = (max(p(i,:)) - min(p(i,:)))*0.5;
    %e(i) = std(p(i,:));
end
pc = 0.6725*ones(nr,1);
pclqgmp = 0.6887*ones(nr,1);
%hold on; plot(x,pc,'b-','LineWidth',3.0); plot(x,pclqgmp,'g-','LineWidth',3.0); errorbar(x,y,e,'r-','LineWidth',3.0); hold off;

%FH = figure;
FH = fig('units','inches','width',8,'height',4);
%axes1 = axes('Parent',FH,'YMinorTick','on','XTick',[0,50,100,200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000,2250,2500,2750,3000],'XMinorTick','on');
axes1 = axes('Parent',FH,'YMinorTick','on','XMinorTick','on');

%set(gca,'XTickLabel',{'$0$';'$500t_{mc}$';'$1000t_{mc}$';'$1500t_{mc}$';'$2000t_{mc}$';'$2500t_{mc}$';'$3000t_{mc}$'});
%plotTickLatex2D();

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 3025]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0.46 0.89]);
hold(axes1,'all');

% Create multiple lines using matrix input to plot
plot(x,pc,'Parent',axes1,'LineWidth',3,'Color',[0 0 1],'DisplayName','Ground Truth');
plot(x,pclqgmp,'Parent',axes1,'LineWidth',3,'Color',[0 0.7 0],'DisplayName','Our Method (t = 0.94 ms)');

% Create errorbar
errorbar(x,y,e,'or','LineWidth',2,'MarkerSize',6,'DisplayName','Monte Carlo (t_{mc} = 0.36 ms)','Color',[1 0 0],'Parent',axes1);

% Create xlabel
%xlabel('Number of Monte-Carlo samples','FontSize',20,'FontName','Helvetica');
%xlabel('Computation Time (milliseconds)','FontSize',20,'FontName','Helvetica');
xlabel('Computation Time (milliseconds)','FontSize',20,'VerticalAlignment','Top');

% Create ylabel
%ylabel('Collision Probability','FontSize',20,'FontName','Helvetica');
ylabel('Collision Probability','FontSize',20);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Position',[0.4938 0.669 0.4863 0.3038],'FontSize',20);

ti = get(gca,'TightInset');
set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);

print -painters -dpdf -r600 mc-car.pdf
