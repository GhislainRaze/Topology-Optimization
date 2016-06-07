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
% * _data{3}_: x coordinate the contour of the structure
% * _data{4}_: y coordinate the contour of the structure
% * _data{5}_: displacement along x
% * _data{6}_: displacement along y
% * _data{7}_: displacement norm
% * _data{8}_: normal x stress
% * _data{9}_: normal y stress
% * _data{10}_: shear xy stress
% * _data{11}_: equivalent Von Mises stress

function data = deformedComputation(x,distrType,u,method,threshold,n)
    global pCon mCon cells nodes oCon;
    
    data = cell(11,1);
    if nargin < 5
        threshold = 0.25;
    end
    if nargin < 6
        n = 9;
    end
    
    % Density evaluation
    [rho,X,Y] = densityField(x,distrType,n);
    
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
    end
    
    
    xc = [];
    yc = [];
    
    for i = 1:size(rho,1)
        Nl = false;
        for j = 1:size(rho,2)
            
        	if rho(i,j) > threshold
                % Contour identification
                if ~Nl
                    xc = [xc,X(i,j)];
                    yc = [yc,Y(i,j)];
                elseif i == 1 || i == size(rho,1) ||...
                        j == size(rho,2)
                    xc = [xc,X(i,j)];
                    yc = [yc,Y(i,j)];
                elseif i~= 1 && rho(i-1,j) < threshold
                    xc = [xc,X(i,j)];
                    yc = [yc,Y(i,j)];
                end
                Nl =  true;
                
                % Shape functions and strains evaluation
                if method == 1
                    nen=nodelabel(and(abs(xN(1,:)-X(i,j))<mCon.dm(1),abs(xN(2,:)-Y(i,j))<mCon.dm(2)));
                    [phi,dphidx,dphidy]=MLSShape(xN(:,nen),[X(i,j);Y(i,j)],dm,pn);
                    uX(i,j) = u(nen,1)'*phi;
                    uY(i,j) = u(nen,2)'*phi;
                    epsXX(i,j) = u(nen,1)'*dphidx;
                    epsYY(i,j) = u(nen,2)'*dphidy;
                elseif method == 2
                    nx = floor(X(i,j)/mCon.cdx);
                    if nx == mCon.mx
                        nx = nx-1;
                    end
                    ny = floor((Y(i,j)+pCon.Ly/2)/mCon.cdy);
                    if ny == mCon.my
                        ny = ny-1;
                    end
                    nc = 1 + nx*mCon.my + ny;
                    coord = 2*([X(i,j);Y(i,j)]-cells(nc).x)./cells(nc).dx;
                    [phi,dphidx,dphidy]=FEMShape(coord,length(cells(nc).nen));
                    uX(i,j) = u(cells(nc).nen,1)'*phi;
                    uY(i,j) = u(cells(nc).nen,2)'*phi;
                    epsXX(i,j) = u(cells(nc).nen,1)'*dphidx;
                    epsYY(i,j) = u(cells(nc).nen,2)'*dphidy;
                end
                epsXY(i,j) = (epsXX(i,j)+epsYY(i,j));           % Engineering strain
                
                % Stresses computation
                sig = rho(i,j)^oCon.p*pCon.D*...
                    [epsXX(i,j) ; epsYY(i,j); epsXY(i,j)];
                sigXX(i,j) = sig(1);
                sigYY(i,j) = sig(2);
                sigXY(i,j) = sig(3);
                sigVM(i,j) = sqrt(sigXX(i,j)^2+sigYY(i,j)^2-sigXX(i,j)*sigYY(i,j)+3*sigXY(i,j)^2);
            else
                % Contour identification
                if Nl
                   xc = [xc,X(i,j-1)];
                   yc = [yc,Y(i,j-1)];
                elseif i~= 1 && rho(i-1,j) > threshold 
                   xc = [xc,X(i-1,j)];
                   yc = [yc,Y(i-1,j)];
                end
                Nl = false;
            end
        end
    end
    
    % Displacement norm
    uM = sqrt(uX.^2+uY.^2);
    
    % Contour ordering
    xco = zeros(size(xc));
    yco = zeros(size(yc));
    xco(1) = xc(1);
    yco(1) = yc(1);
    xc(1) = [];
    yc(1) = [];
    for i = 1 : length(xco)-1
        [~,indMin] = min(sqrt(((xco(i)-xc)/pCon.Lx).^2+((yco(i)-yc)/pCon.Ly).^2));
        xco(i+1) = xc(indMin);
        yco(i+1) = yc(indMin);
        xc(indMin) = [];
        yc(indMin) = [];
    end
       
    % Saving the result into data
    data{1} = X;
    data{2} = Y;
    data{3} = xco;
    data{4} = yco;
    data{5} = uX;
    data{6} = uY;
    data{7} = uM;
    data{8} = sigXX;
    data{9} = sigYY;
    data{10} = sigXY;
    data{11} = sigVM;
    
end