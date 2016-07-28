%% Density plot
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Plots the material distribution in levels of gray. The mass
% nodes/elements centers are represented by blue dots.
function h = densityPlot(x,distrType,h,figureTitle)
    
    global pCon mnodes

    
    if nargin < 2
        distrType = 1;
    end
    if nargin < 1
        x = mnodesToVector(mnodes,distrType);
    end
    if nargin < 3
        h = figure;
    end
    if nargin < 4
        figureTitle = 'Material distribution';
    end
    
    
    [rho,X,Y,xi] = densityField(x,distrType);
    cla
    hold on
    contourf(X,Y,rho,50,'linestyle','none')
    colormap gray
    colormap(flipud(colormap))
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
    title(figureTitle)
    drawnow
end