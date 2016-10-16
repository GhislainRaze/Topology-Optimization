%% Deformed structure computation
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Postprocesses the history data to give some fields in the data cell
% structure.
%
% * _data{1}_: x coordinate of all the fields
% * _data{2}_: y coordinate of all the fields
% * _data{3}_: displacement along x
% * _data{4}_: displacement along y
% * _data{5}_: displacement norm
% * _data{6}_: normal x stress
% * _data{7}_: normal y stress
% * _data{8}_: shear xy stress
% * _data{9}_: equivalent Von Mises stress

function data = deformedComputation(x,distrType,u,method,threshold)
    global pCon mCon cells nodes oCon;
    
    data = cell(11,1);
    if nargin < 5
        threshold = 0.25;
    end
    
    % Density evaluation
    [rho,X,Y] = densityField(x,distrType);
    
    % Fields initialization
    uX = nan(size(rho));
    uY = uX;
    epsXX = uX;
    epsYY = uX;
    epsXY = uX;
    sigXX = uX;
    sigYY = uX;
    sigXY = uX;
    sigVM = uX;
    
    if method == 1
        nodeLabel = 1:length(nodes);
        xN = [nodes.x];
    elseif method == 2
        xCells = [cells.x];
        dxCells = [cells.dx];
    end
    
    
    xc = [];
    yc = [];
    
    a = uX;
    
    for i = 1:size(rho,1)
        for j = 1:size(rho,2)
            
        	if rho(i,j) > threshold
                % Shape functions and strains evaluation
                if method == 1
                    nen=nodeLabel(and(abs(xN(1,:)-X(i,j))<mCon.dm(1),abs(xN(2,:)-Y(i,j))<mCon.dm(2)));
                    [phi,dphidx,dphidy]=MLSShape(xN(:,nen)',[X(i,j);Y(i,j)],mCon.dm,mCon.pn);
                    uX(i,j) = phi*u(nen,1);
                    uY(i,j) = phi*u(nen,2);
                    epsXX(i,j) = dphidx*u(nen,1);
                    epsYY(i,j) = dphidy*u(nen,2);
                elseif method == 2
                    tol = 1e-9;
                    nc = find(and(and(xCells(1,:)-dxCells(1,:)/2-X(i,j)<tol,xCells(2,:)-dxCells(2,:)/2-Y(i,j)<tol),...
                        and(xCells(1,:)+dxCells(1,:)/2-X(i,j)>-tol,xCells(2,:)+dxCells(2,:)/2-Y(i,j)>-tol)));
                    if ~isempty(nc)
                        nc = nc(1);
                        coord = 2*([X(i,j);Y(i,j)]-cells(nc).x)./cells(nc).dx;
                        [phi,dphidx,dphidy]=FEMShape(coord,length(cells(nc).nen));
                        uX(i,j) = u(cells(nc).nen,1)'*phi;
                        uY(i,j) = u(cells(nc).nen,2)'*phi;
                        epsXX(i,j) = u(cells(nc).nen,1)'*dphidx;
                        epsYY(i,j) = u(cells(nc).nen,2)'*dphidy;
                    end
                end
                epsXY(i,j) = (epsXX(i,j)+epsYY(i,j));           % Engineering strain
                
                % Stresses computation
                sig = pCon.E*rho(i,j)^oCon.p*pCon.D*...
                    [epsXX(i,j) ; epsYY(i,j); epsXY(i,j)];
                sigXX(i,j) = sig(1);
                sigYY(i,j) = sig(2);
                sigXY(i,j) = sig(3);
                sigVM(i,j) = sqrt(sigXX(i,j)^2+sigYY(i,j)^2-sigXX(i,j)*sigYY(i,j)+3*sigXY(i,j)^2);
            end
        end
    end
    
    % Displacement norm
    uM = sqrt(uX.^2+uY.^2);
    
    
    % Saving the result into data
    data{1} = X;
    data{2} = Y;
    data{3} = uX;
    data{4} = uY;
    data{5} = uM;
    data{6} = sigXX;
    data{7} = sigYY;
    data{8} = sigXY;
    data{9} = sigVM;
    
end