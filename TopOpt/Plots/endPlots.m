%% End plots
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Plots the miscellaneous figures asked by the user after the optimization
% algorithm.


if plotInitial  
    densityPlot(history.x(:,1),massChoice,figure,'Initial Configuration');
end

if plotFinal
    densityPlot(history.x(:,end),massChoice,figure,'Final Configuration');
end

if plotMesh
    discretizationPlot;
end

if plotDEvolution
    densityEvolutionPlot(history,massChoice)
end

if plotCEvolution
    contourEvolutionPlot(history,massChoice)
end


if plotCompliance
    figure
    set(gca,'fontsize',20)
    plot(0:length(history.C)-1,history.C,'linewidth',2);
    xlim([0 length(history.C)])
    grid on
    xlabel('Iteration')
    ylabel('Compliance')
    title('Compliance evolution')
end

if plotDeformed
end