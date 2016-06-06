%% Continuation strategy
%
% Code developed by Ghislain Raze under the supervision of Prof. Joseph
% Morlier
%
% Initial code by Johannes T. B. Overvelde
%
% <http://www.overvelde.com>
%
% Changes the penalization factor _p_ according to the continuation step
% presented by Groenwold and Etman.


function p = continuation(p,iter,pmax)

    if iter <= 20
        p = 1;
    else
        p = min(pmax,1.02*p^(iter-1));
    end

end