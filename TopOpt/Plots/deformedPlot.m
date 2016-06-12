%% Deformed configuration plot
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Plots the field given by _data{n}_ over the deformed structure. The
% deformation is amplified by a factor _mag_.

function deformedPlot(data,n,mag)
    global pCon
    
    X = data{1};
    Y = data{2};
    uX = data{3};
    uY = data{4};
    field = data{n+2};
    
    cla
    %cm = [flipud(gray(64));jet(64)];
    hold on
    set(gca,'units','normalized','position',[1/8,1/4,3/4,1/2])
    pcol = pcolor(X+mag*uX,Y+mag*uY,field);
    shading interp
    colorbar
    surface(X,Y,field-max(max(field))+min(min(field)),...
        'FaceColor',[0.75 0.75 0.75],'FaceAlpha',1,'LineStyle','none')
    uistack(pcol,'top')
    caxis([min(min(field)),max(max(field))])
    hold off
    axis equal
    set(gca,'fontsize',20)
    xlim([-0.1*pCon.Lx ; 1.1*pCon.Lx])
    ylim([-1.1*pCon.Ly/2;1.1*pCon.Ly/2])
    set(gca,'XTick',[0 ; pCon.Lx])
    set(gca,'YTick',[-pCon.Ly/2;0;pCon.Ly/2])
    set(gca,'XTickLabel',{'' ; ''})
    set(gca,'YTickLabel',{'';'';''})
    format_ticks(gca,{'0','L_x'},...
        {'-L_y/2','0','L_y/2'});
    title('Deformed configuration')
    drawnow
    
end