

    merCon.merge = true;                    % Merging process after convergence
    merCon.tolDist = 0.1;                   % Relative merging tolerance on distances
    merCon.tolDim = 0.25;                   % Relative merging tolerance on distances
    merCon.tolTheta = 5*pi/180;             % Absolute merging tolerance on angles
    merCon.densityRadius = 0.37;            % Merging radius
    
    
    merCon.suppressZeroWidthNodes = true;
    merCon.suppressOutsideNodes = true;
    merCon.suppressIsolatedNodes = true;
    merCon.isolationRatio = 1;