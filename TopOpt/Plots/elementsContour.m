%% Elements contour plot
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Plots the mass nodes/elements contour. The blue points represent the
% center coordinates of the node/element, the dashed pale blue lines
% represent its domain of influence and the full purple lines represent
% where the density is half of its maximum value.

function h = elementsContour(x,distrType,h,figureTitle)
    global pCon

    if nargin < 2
        distrType = 1;
    end
    if nargin < 3
        h = figure;
    end
    if nargin < 4
        figureTitle = 'Elements Contour';
    end

    if distrType == 3
        nd = 5;
    else
        nd = distrType+1;
    end
    
    [xi,thetai,dmi] = mnodesData(nd,x);
    dmi = dmi/2;
    
    
    cla
    hold on
    for i = 1 : length(xi)
    xd = [xi(1,i)+dmi(1,i)*cos(thetai(i))-dmi(2,i)*sin(thetai(i)),...
        xi(1,i)-dmi(1,i)*cos(thetai(i))-dmi(2,i)*sin(thetai(i)),...
        xi(1,i)-dmi(1,i)*cos(thetai(i))+dmi(2,i)*sin(thetai(i)),...
        xi(1,i)+dmi(1,i)*cos(thetai(i))+dmi(2,i)*sin(thetai(i)),...
        xi(1,i)+dmi(1,i)*cos(thetai(i))-dmi(2,i)*sin(thetai(i));...
        xi(2,i)+dmi(1,i)*sin(thetai(i))+dmi(2,i)*cos(thetai(i)),...
        xi(2,i)-dmi(1,i)*sin(thetai(i))+dmi(2,i)*cos(thetai(i)),...
        xi(2,i)-dmi(1,i)*sin(thetai(i))-dmi(2,i)*cos(thetai(i)),...
        xi(2,i)+dmi(1,i)*sin(thetai(i))-dmi(2,i)*cos(thetai(i)),...
        xi(2,i)+dmi(1,i)*sin(thetai(i))+dmi(2,i)*cos(thetai(i))];
    
    xc = [xi(1,i)+0.37*dmi(1,i)*cos(thetai(i))-0.37*dmi(2,i)*sin(thetai(i)),...
        xi(1,i)-0.37*dmi(1,i)*cos(thetai(i))-0.37*dmi(2,i)*sin(thetai(i)),...
        xi(1,i)-0.37*dmi(1,i)*cos(thetai(i))+0.37*dmi(2,i)*sin(thetai(i)),...
        xi(1,i)+0.37*dmi(1,i)*cos(thetai(i))+0.37*dmi(2,i)*sin(thetai(i)),...
        xi(1,i)+0.37*dmi(1,i)*cos(thetai(i))-0.37*dmi(2,i)*sin(thetai(i));...
        xi(2,i)+0.37*dmi(1,i)*sin(thetai(i))+0.37*dmi(2,i)*cos(thetai(i)),...
        xi(2,i)-0.37*dmi(1,i)*sin(thetai(i))+0.37*dmi(2,i)*cos(thetai(i)),...
        xi(2,i)-0.37*dmi(1,i)*sin(thetai(i))-0.37*dmi(2,i)*cos(thetai(i)),...
        xi(2,i)+0.37*dmi(1,i)*sin(thetai(i))-0.37*dmi(2,i)*cos(thetai(i)),...
        xi(2,i)+0.37*dmi(1,i)*sin(thetai(i))+0.37*dmi(2,i)*cos(thetai(i))];
    
    plot(xi(1,i),xi(2,i),'ob','linewidth',2.5)
    line(xd(1,:),xd(2,:),'linewidth',1.5,'LineStyle','--','color',[0.75 0.75 1])
    line(xc(1,:),xc(2,:),'linewidth',2,'color',[0.375 0.375 0.5])
    end
    
    regionsPlot(pCon.holes,[0 0 0])
    regionsPlot(pCon.filledRegions,[0.125 0.125 0.25])
    
    
    hold off
    axis equal
    set(gca,'fontsize',20)
    xlim([0 ; pCon.Lx])
    ylim([-pCon.Ly/2;pCon.Ly/2])
    set(gca,'XTick',[0 ; pCon.Lx])
    set(gca,'YTick',[-pCon.Ly/2;0;pCon.Ly/2])
    set(gca,'XTickLabel',{'' ; ''})
    set(gca,'YTickLabel',{'';'';''})
    format_ticks(gca,{'0','L_x'},...
        {'-L_y/2','0','L_y/2'});
    title(figureTitle)
    drawnow
end

    