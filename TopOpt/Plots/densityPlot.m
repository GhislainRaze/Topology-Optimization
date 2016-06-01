%% Density Plot
%
% Plots the material distribution in levels of gray. The nodes are also
% plotted in blue.

function  h = densityPlot(h,figureTitle)

GlobalConst
if nargin < 1
    h = figure;
end
if nargin < 2
    figureTitle = 'Material Distribution';
end

Nx = 4*mCon.nx;
Ny = 4*mCon.ny;
dx = pCon.Lx/Nx;
dy = pCon.Ly/Ny;


[XX,YY] = meshgrid(0:dx:pCon.Lx,-pCon.Ly/2:dy:pCon.Ly/2);

rhoA = zeros(size(XX));

xi = [mnodes.x];
thetai = [mnodes.theta];
dmi = [mnodes.l]/2;
mi = [mnodes.m];

for ii = 1 : size(XX,1)
    
    for jj = 1 : size(XX,2)
        
        xj = [ XX(ii,jj) , YY(ii,jj) ]; 
        
        rhoA(ii,jj) = asymptoticDensity(xj,xi,thetai,dmi,mi,mmCon.rhoMin,...
            mmCon.rhoMax,1);
    end
    
end

h = figure(h);
clf;
hold on
contourf(XX,YY,rhoA,50,'linestyle','none')
colormap gray
set(gca,'fontsize',20)
caxis([0 1])
colorbar('fontsize',20)
plot(xi(1,:),xi(2,:),'o','linewidth',2.5)
xlim([0 ; pCon.Lx])
ylim([-pCon.Ly/2;pCon.Ly/2])
set(gca,'XTick',[0 ; pCon.Lx])
set(gca,'YTick',[-pCon.Ly/2;0;pCon.Ly/2])
set(gca,'XTickLabel',{'' ; ''})
set(gca,'YTickLabel',{'';'';''})
format_ticks(gca,{'0','L_x'},...
    {'-L_y/2','0','L_y/2'});
hold off
view([0 0 1])
title(figureTitle)
end