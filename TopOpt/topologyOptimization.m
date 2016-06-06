%% Topology optimization
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>

clear all
close all
clc



addpath('Optimization/');
addpath('Optimization/LineSearch/');
addpath('Optimization/Optimizers/');
addpath('Optimization/Scripts/');
addpath('Constants/');
addpath('Discretization/');
addpath('MaterialDistribution/');
addpath('EFG/');
addpath('FEM/');
addpath('IIEFG/');
addpath('Plots/');
addpath('Display/');

GlobalConst;
optimizationConstants;
%% Method, material distribution and optimization algorithm
% The elastic problem can be discretized thanks to three different methods
% : the Element-Free Galerkin (EFG) method, the Finite Element Method (FEM)
% or the Improved Interpolating EFG (IIEFG) method.
%
% The material distribution 
methodChoice = 2;           % 1: EFG, 2: FEM, 3: IIEFG 
massChoice = 2;             % 1: Mass nodes, 2: Undeformable structural 
                            % members, 3: Deformable structural members 
optimChoice = 1;            % 1: Overvelde's algorithm, 2: fminunc or 
                            % fmincon, 3: Genetic algorithm
                            
                            
%% Plots and movies

plotInitial = true;         % Plots the initial mass distribution
plotFinal = true;           % Plots the final mass distribution
plotMesh = true;            % Plots the discretization mesh
plotDEvolution = true;      % Plots the density at given iteration
plotCEvolution = true;      % Plots the elements contour at given iteration
plotCompliance = true;      % Plots the compliance evolution
plotDeformed = true;        % Plots the final deformed configuration



%% Optimization process

switch optimChoice
    
    case 1
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                  Overvelde''s algorithm               ')
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        history = overvelde(massChoice,methodChoice);
    case 2
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                    Steepest descent                   ')  
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        history = steepestDescent(massChoice,methodChoice);
    case 3
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                  Conjugated gradients                 ')
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        history = conjugatedGradients(massChoice,methodChoice);
    case 4
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                    Quasi-Newton BFGS                  ')
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        history = quasiNewtonBFGS(massChoice,methodChoice);
    case 5
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                  Matlab fminunc/fmincon               ')  
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        history = matlabFminunc(massChoice,methodChoice);
    case 6
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                 Matlab genetic algorithm              ')
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        history = matlabGa(massChoice,methodChoice);    
end

endPlots;
