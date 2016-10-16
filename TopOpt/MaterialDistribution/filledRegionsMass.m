function m = filledRegionsMass(filledRegions)

    m = 0;
    for i = 1 : length(filledRegions)
       if filledRegions(i).type == 1
           m = m + prod(filledRegions(i).l);
       elseif filledRegions(i).type == 2
           m = m + pi * filledRegions(i).r^2;
       end
    end

end