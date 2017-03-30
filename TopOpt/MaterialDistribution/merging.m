
function [mnodes,nMerged] = merging(mnodes) 

    global merCon

    nMerged = 0;                                        % Total number of merged elements
    nMergedTemp = 1;                                    % Temporary number of merged element

    alpha = (1+merCon.tolDist)*merCon.densityRadius;    % Tolerance on distances


    while nMergedTemp ~= 0                              % Merging loop
        nMergedTemp = 0;
        i = 1;
        while i < length(mnodes)
            j = i+1;
            while j <= length(mnodes)

                % Angles criteria (trmif(mod(theta,pi/2),[0,pi/4,pi/2])
                % gives the gap between theta and the nearest multiple of
                % pi/2
                if pi/4*(abs(trimf(mod(mnodes(i).theta,pi/2),[0,pi/4,pi/2])+...
                        trimf(mod(mnodes(j).theta,pi/2),[0,pi/4,pi/2]))) < merCon.tolTheta
                    % Common basis comparison
                    theta = atan2(mnodes(j).x(2)-mnodes(i).x(2),mnodes(j).x(1)-mnodes(i).x(1));
                    d = norm(mnodes(j).x-mnodes(i).x);
                    [ll1,lt1] = directionalComponents(mnodes(i),theta);
                    [ll2,lt2] = directionalComponents(mnodes(j),theta);

                    % Distance between centers criterion
                    if d <= alpha*(ll1+ll2)

                        % Distance between edges criterion
                        if abs(lt1-lt2) <= merCon.tolDim*(lt1+lt2)

                            % Creating the new merging node and deleting the
                            % merged ones
                            newmnode.m = mnodes(i).m + mnodes(j).m;
                            newmnode.x = (mnodes(i).m*mnodes(i).x + mnodes(j).m*mnodes(j).x)/newmnode.m;
                            newmnode.theta = theta;
                            newmnode.l = [2*(ll1+ll2);(prod(mnodes(i).l)+prod(mnodes(j).l))/(2*(ll1+ll2))];
                            mnodes = [mnodes,newmnode];
                            mnodes([i,j]) = [];
                            nMergedTemp = nMergedTemp+1;
                        end

                    end
                end
                j = j+1;
            end
            i = i+1;
        end
        nMerged = nMerged + nMergedTemp;
    end
end
            