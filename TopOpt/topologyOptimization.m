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
addpath('Constants/');
addpath('Discretization/');
addpath('MaterialDistribution/');
addpath('EFG/');
addpath('FEM/');
addpath('IIEFG/');
addpath('Plots/');

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
optimChoice = 4;            % 1: Overvelde's algorithm, 2: fminunc or 
                            % fmincon, 3: Genetic algorithm
                            
                            
%% Plots and movies

plotInitial = true;         % Plots the initial mass distribution
plotFinal = true;           % Plots the final mass distribution
plotMesh = true;            % Plots the discretization mesh
plotIter = true;            % Plots the configuration at given iteration
plotCompliance = true;      % Plots the compliance evolution



%% Optimization process

switch optimChoice
    
    case 1
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                  Overvelde''s algorithm               ')
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        overvelde(massChoice,methodChoice);
    case 2
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                    Steepest descent                   ')  
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        steepestDescent(massChoice,methodChoice);
    case 3
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                  Conjugated gradients                 ')
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        conjugatedGradients(massChoice,methodChoice);
    case 4
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                    Quasi-Newton BFGS                  ')
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        quasiNewtonBFGS(massChoice,methodChoice);
    case 5
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                  Matlab fminunc/fmincon               ')  
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        matlabFminunc(massChoice,methodChoice);
    case 6
        disp('=======================================================')
        disp('                  TOPOLOGY OPTIMIZATION                ')
        disp('                 Matlab genetic algorithm              ')
        methodDisplay(methodChoice)
        massDisplay(massChoice)
        disp('=======================================================')
        matlabGa(massChoice,methodChoice);    
end


