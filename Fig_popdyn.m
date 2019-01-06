function Fig_popdyn(X1, YMatrix1,np)
%CREATEFIGURE(X1, YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 21-Aug-2016 19:14:15

% Create figure
figure1 = figure;
% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.146245059288538 0.129179331306991 0.80764163372859 0.795820668693009]);
hold(axes1,'on');
plot(0,0,'k-','LineWidth',2);  plot(0,0,'k-.','LineWidth',2);   plot(0,0,'k--','LineWidth',2);
%legend('Lysogens','CRISPR Bacteria','Bacteriophage');
plot(X1,YMatrix1(:,2*(2^np)+1),'b--','LineWidth',2);
plot(X1,YMatrix1(:,2*(2^np)+2),'r--','LineWidth',2);
plot(X1,YMatrix1(:,2*(2^np)+3),'k--','LineWidth',2);
% Create multiple lines using matrix input to semilogy
semilogy1 = semilogy(X1,YMatrix1(:,1:2*(2^np)),'LineWidth',2,'Parent',axes1);
set(semilogy1(1),'Color',[1 0 1],'DisplayName','L_{000}');
set(semilogy1(2),'Color',[0 0 1],'DisplayName','L_{001}');
set(semilogy1(3),'Color',[0 0.447058826684952 0.74117648601532],'DisplayName','L_{010}');
set(semilogy1(4),'Color',[1 0 0],'DisplayName','L_{011}');
set(semilogy1(5),'LineStyle','-',...
    'Color',[0.466666668653488 0.674509823322296 0.18823529779911],'DisplayName','L_{100}');
set(semilogy1(6),'Color',[0.749019622802734 0.749019622802734 0],'DisplayName','L_{101}');
set(semilogy1(7),...
    'Color',[0.929411768913269 0.694117665290833 0.125490203499794],'DisplayName','L_{110}');
set(semilogy1(8),'Color',[0 0 0],'DisplayName','L_{111}');
set(semilogy1(9),'LineStyle','-.','Color',[1 0 1],'DisplayName','C_{000}');
set(semilogy1(10),'LineStyle','-.','Color',[0 0 1],'DisplayName','C_{001}');
set(semilogy1(11),'LineStyle','-.',...
    'Color',[0 0.447058826684952 0.74117648601532],'DisplayName','C_{010}');
set(semilogy1(12),'LineStyle','-.','Color',[1 0 0],'DisplayName','C_{011}');
set(semilogy1(13),'LineStyle','-.',...
    'Color',[0.466666668653488 0.674509823322296 0.18823529779911],'DisplayName','C_{100}');
set(semilogy1(14),'LineStyle','-.',...
    'Color',[0.749019622802734 0.749019622802734 0],'DisplayName','C_{101}');
set(semilogy1(15),'LineStyle','-.',...
    'Color',[0.929411768913269 0.694117665290833 0.125490203499794],'DisplayName','C_{110}');
set(semilogy1(16),'LineStyle','-.','Color',[0 0 0],'DisplayName','C_{111}');

% Create xlabel
xlabel('Time t (hrs)');

% Create ylabel
ylabel('Population densities (virions or cells / cm^2)');

% Create title
title('Lysogen Insertion Sequence \{\}');

%% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0 1000]);
%% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[1e-07 395833388.749007]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',16,'LineStyleOrderIndex',3,'XMinorTick','on',...
    'YMinorTick','on','YScale','log');