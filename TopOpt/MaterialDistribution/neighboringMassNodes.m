%% Neighboring mass nodes
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Finds the integration points in the influence domain of the mass nodes.
% Updates the _cells_ structure.

function cells = neighboringMassNodes(mnodes,cells)

    if ~isempty(mnodes)
        x1=[mnodes.x];
        theta=[mnodes.theta];
        dm=[mnodes.l]/2;
        nodelabel=1:length(mnodes);
        
        %Find integraion points in influence domain of nodes
        for j=1:length(cells)
            for k=1:cells(j).ni
                r1 = (x1(1,:)-cells(j).int(k).x(1)).*cos(theta)+...
                    (x1(2,:)-cells(j).int(k).x(2)).*sin(theta);
                r2 = -(x1(1,:)-cells(j).int(k).x(1)).*sin(theta)+...
                    (x1(2,:)-cells(j).int(k).x(2)).*cos(theta);
                cells(j).int(k).nemn=nodelabel(and(abs(r1)<dm(1,:),abs(r2)<dm(2,:)));
            end
        end
    end

end