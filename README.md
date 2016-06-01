# Topology-Optimization
Topology optimization based on the Moving Node Algorithm using EFG, FEM or IIEFG methods.

## Introduction
This set of Matlab files are used for topology optimization using the Moving Node Approach (MNA) for 2D plane strain problems. In this approach, the material distribution is decoupled from the discretization.

### Material distribution
The material distribution is based on mass nodes. The density at one point can be computed thanks to a kernel approximation using cubic spline shape functions.

## Main file
The main file topologyOptimization.m is a Matlab script that launches the optimizers. The user can change the script to set
* The discretization method (EFG or FEM)
* The optimization variables (mass nodes, undeformable structural members or deformable structural members)
* The optimization algorithm
