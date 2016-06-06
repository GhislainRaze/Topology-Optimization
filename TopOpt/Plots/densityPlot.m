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
    
    global pCon mmCon

    if nargin < 1
        x = mnodesToVector();
    end
    if nargin < 2
        distrType = 1;
    end
    if nargin < 3
        h = figure;
    end
    if nargin < 4
        figureTitle = 'Material distribution';
    end
    
    if distrType == 3
        nd = 5;
    else
        nd = distrType+1;
    end
    nmn = length(x)/nd;
    nodelabel = 1:nmn;
    
    
    dx = pCon.Lx/50;
    dy = pCon.Ly/50;
    
    globalConst;    
    [XX,YY] = meshgrid(0:dx:pCon.Lx,-pCon.Ly/2:dy:pCon.Ly/2);
    xi = zeros(2,nmn);
    thetai = zeros(1,nmn);
    dmi = xi;
    mi = thetai;
    
    
    rhoA = zeros(size(XX));
    for i = 1:nd:length(x)
        nn = ceil(i/nd);
        xi(1,nn) = x(i);
        xi(2,nn) = x(i+1);
        if distrType >= 2
            thetai(nn) = x(i+2);
        else
            thetai(nn) = 0;
        end
        if distrType ==3
            dmi(1,nn) = x(i+3);
            dmi(2,nn) = x(i+4);
            mi(nn) = x(i+3)*x(i+4)/(mmCon.d^2);
        else
            dmi(1:2,nn) = mmCon.dm;
            mi(nn) = mmCon.mi;
        end
    end   

    
    
    for ii = 1 : size(XX,1)

        for jj = 1 : size(XX,2)

            xj = [ XX(ii,jj) , YY(ii,jj) ]; 
            
            r1 = (xi(1,:)-xj(1)).*cos(thetai)+(xi(2,:)-xj(2)).*sin(thetai);
            r2 = -(xi(1,:)-xj(1)).*sin(thetai)+(xi(2,:)-xj(2)).*cos(thetai);
            
            nemn=nodelabel(and(abs(r1)<dmi(1,:),abs(r2)<dmi(2,:)));
            
            if ~isempty(nemn)
                xinn = xi(:,nemn);
                thetainn = thetai(nemn);
                dminn = dmi(:,nemn);
                minn = mi(nemn);

                rhoA(ii,jj) = asymptoticDensity(xj,xinn,thetainn,dminn,minn,mmCon.rhoMin,...
                    mmCon.rhoMax,distrType,false);
            else
                rhoA(ii,jj) = mmCon.rhoMin;
            end
        end

    end
    cla
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
    title(figureTitle)
    drawnow
end