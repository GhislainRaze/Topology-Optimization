%% Deformed structure evolution plot
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Plots a given field over the deformed structure at a given iteration.
% _threshold_ (optional, default value = 0.25) is the threshold of density 
% under which the results are not plotted (it should therefore range from 0
% to 1). _np_ (optional, default value = 9) is the number of starting 
% points for the discretization of the density field.

function deformedEvolutionPlot(history,distrType,method,threshold,np)

    if nargin < 4
        threshold = 0.25;
    end
    if nargin < 5
        np = 9;
    end
    h = figure;
    mag = 0.01;
    iterMax = length(history.C)-1;
    data = deformedComputation(history.x(:,end),distrType,...
        history.u(:,end-1:end),method,threshold,np);

    txtMag = uicontrol(h,'style','text','Units','normalized','Position',...
        [1/16 1/16 3/16 1/24],'String','Magnifying factor','FontSize',10);
    editer =  uicontrol(h,'style','edit','String',mag,...
        'Fontsize',10,'Units','normalized','Position',[4/16 1/16 1/8 1/24]);
    popupDeformed2 = @(hObject,event) popupDeformed(hObject,event,data,editer);

    popupmenu = uicontrol(h,'style','popupmenu','String',{'Displacement along x',...
        'Displacement along y','Displacement','Normal x stress','Normal y stress',...
        'Shear xy stress','Von Mises equivalent stress'},...
        'Units','normalized','Position',[1/16 7/8 24/48 1/16],...
        'Fontsize',10,'Callback',popupDeformed2);
    popupDeformed2(popupmenu,1); 
    buttonDeformed2 = @(hObject,event) buttonDeformed(hObject,event,data,...
        popupmenu,editer);
    button =  uicontrol(h,'style','pushbutton','String','Update plot',...
        'Units','normalized','Position',[6/16 1/16 1/6 1/24],...
        'Fontsize',10,'Callback',buttonDeformed2);
    sliderDeformed2 = @(hObject,event) sliderDeformed(hObject,event);
    slider = uicontrol(h,'style','slider','Min',0,'Max',iterMax,'Value',iterMax,...
        'Units','normalized','Position',[15/16 1/8 1/48 6/8],...
        'SliderStep',[1/(iterMax-1) 10/(iterMax-1)],'Callback',sliderDeformed2);
    txtLimInf = uicontrol(h,'style','text','Units','normalized','Position',...
        [15/16 1/16 1/24 1/24],'String',...
        num2str(0));
    txtLimSup = uicontrol(h,'style','text','Units','normalized','Position',...
        [15/16 7/8 1/24 1/24],'String',...
        num2str(iterMax));

    function sliderDeformed(slider,event)

        n = get(popupmenu, 'Value');
        [mag,success] = str2num(get(editer, 'String'));
        if ~success
            warndlg('Please insert a number in the magnifying factor box.')
        end

        i = round(get(slider, 'Value'));
        x = history.x(:,i+1);
        u = history.u(:,2*i+1:2*i+2);
        data = deformedComputation(x,distrType,u,method,threshold,np);

        deformedPlot(data,n,mag)


        popupDeformed2 = @(hObject,event) popupDeformed(hObject,event,data,editer);
        buttonDeformed2 = @(hObject,event) buttonDeformed(hObject,event,data,...
            popupmenu,editer);
        set(popupmenu,'Callback',popupDeformed2);
        set(button,'Callback',buttonDeformed2);
    end

end