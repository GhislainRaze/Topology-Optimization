%% Discretization plot
%
% Plots the discretization
%
% * The nodes are blue
% * The mass nodes are purple
% * The boundary nodes are red
% * The cells/elements contours are black
% * The Gauss points are green


%Plot the initial configuration
    figure, hold on
    %set(gcf,'Position',get(0,'ScreenSize'))
    plot([0 0 pCon.Lx pCon.Lx 0],[-pCon.Ly/2 pCon.Ly/2 pCon.Ly/2 -pCon.Ly/2 -pCon.Ly/2],'black','linewidth',2)
    for i=1:mCon.m
        x1=cells(i).x(1)-cells(i).dx(1)/2;
        x2=cells(i).x(1)+cells(i).dx(1)/2;
        y1=cells(i).x(2)-cells(i).dx(2)/2;
        y2=cells(i).x(2)+cells(i).dx(2)/2;
        plot([x1 x1 x2],[y1 y2 y2],'--','color','black')
    end
    x=[];
    for i=1:mCon.m
        for j=1:cells(i).ni
            x=[x; cells(i).int(j).x(1) cells(i).int(j).x(2)];
        end
    end
    plot(x(:,1),x(:,2),'.','color','g','linewidth',2,'markersize',16)
    x=[];
    for i=1:mCon.mb
        for j=1:bcells(i).ni
            x=[x; bcells(i).int(j).x(1) bcells(i).int(j).x(2)];
        end
    end
    plot(x(:,1),x(:,2),'.','color','r','linewidth',2,'markersize',16)
    x=[];
    if mCon.mp~=0
        for i=mCon.mb+1:mCon.mb+mCon.mp
            for j=1:bcells(i).ni
                x=[x; bcells(i).int(j).x(1) bcells(i).int(j).x(2)];
            end
        end
        plot(x(:,1),x(:,2),'x','color','r','linewidth',2,'markersize',16)
    end
    x=[nodes.x]';
    plot(x(:,1),x(:,2),'o','color','b','linewidth',2)
    x=[mnodes.x]';
    plot(x(:,1),x(:,2),'o','color','m','linewidth',2)
    set(gca,'XTick',[0;pCon.Lx])
    set(gca,'YTick',[-pCon.Ly/2;0;pCon.Ly/2])
    set(gca,'XTickLabel',{'' ; ''})
    set(gca,'YTickLabel',{'';'';''})
    [hx,hy] = format_ticks(gca,{'\fontsize{20}0','\fontsize{20}L_x'},...
    {'\fontsize{20}-L_y/2','\fontsize{20}0','\fontsize{20}L_y/2'});
    set(gca,'fontsize',20)
    %axis([-0.001*pCon.Lx +1.001*pCon.Lx -1.001*pCon.Ly/2 1.001*pCon.Ly/2])
    title('Discretization')