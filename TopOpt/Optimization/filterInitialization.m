%% Filter Initialization
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Initializes the filter convolution matrix. The inputs are
%
% *_cells_: the structure array containing the cells/elements definition
% *_nG_: the number of Gauss points in a cell along one direction
% *_rmin_: the radius under which variations are smoothed
%
% The outputs are
%
% * _H_: the filter convolution matrix
% * _Hs_: the sum of the filter convolution matrix lines


function [H,Hs] = filterInitialization(cells,nG,rmin)
m = length(cells);
H = zeros(m*nG^2);

for ic=1:m                                         
    for ip=1:cells(ic).ni
        i = (ic-1)*nG^2+ip;
        for icc=1:m
            for ipc = 1:cells(icc).ni 
                j = (icc-1)*nG^2+ipc;
                H(i,j) = max(0,rmin-norm(cells(icc).int(ipc).x-cells(ic).int(ip).x));
            end
        end
    end
end

H = sparse(H);
Hs = sum(H,2);

end