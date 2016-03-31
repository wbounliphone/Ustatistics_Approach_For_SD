function cases = findCovCases()
% Author: Wacha Bounliphone - wacha.bounliphone@centralesupelec.fr
% Copyright (c) 2016
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
% 
% If you use this software in your research, please cite:%
% Bounliphone, W. &  Blaschko, M.B. (2016).  
% A U-statistic Approach to Hypothesis Testing for Structure Discovery in 
% Undirected Graphical Models
%
%-----------------------------------------------------------
%We encode each case 
%as a base 4 string of length 4 indexing a configuration of cov(S_ij, S_kl)
%where ijkl vary over a set of p variables.  We only need find the cases
%that indicate whether ijkl are equal or different, and we will only derive
%estimators for cases that cannot be computed by another case with a change
%of variables.
%-----------------------------------------------------------
cases = [];
for currcase=0:(4^3-1) % by symmetry i can always be assigned to variable 0, there are max 4^3 remainign configurations of the other variables
    c = currcase;
    l = mod(c,4);
    c = floor(c/4);
    k = mod(c,4);
    c = floor(c/4);
    j = mod(c,4);
    i=0;
    c = [i,j,k,l];
    newcase=true;
    for prevcase = 1:size(cases,1)
        if(isEquivalent(cases(prevcase,:),c))
            newcase=false;
        end
    end
    if(newcase)
        cases = [cases;c];
    end
end
end

function equiv = isEquivalent(c1,c2)
    equiv=false;
    assert(checkCanonical(c1))
    assert(checkCanonical(c2))
    if(isEqual(makeCanonical(c1),makeCanonical(c2)) || ... %cov(S_ij, S_kl)
       isEqual(makeCanonical(c1),makeCanonical([c2(3:4) c2(1:2)])) || ... %cov(S_kl,S_ij)
       isEqual(makeCanonical(c1),makeCanonical([c2(1:2) c2(4:-1:3)])) || ... %cov(S_ij,S_lk) 
       isEqual(makeCanonical(c1),makeCanonical([c2(2:-1:1) c2(4:-1:3)])) || ... %cov(S_ji,S_lk)
       isEqual(makeCanonical(c1),makeCanonical([c2(2:-1:1) c2(3:4)])) || ... %cov(S_ji,S_kl)
       isEqual(makeCanonical(c1),makeCanonical([c2(4:-1:3) c2(2:-1:1)])) || ... %cov(S_lk, S_ji)
       isEqual(makeCanonical(c1),makeCanonical([c2(4:-1:3) c2(1:2)])) || ... %cov(S_lk, S_ij)
        (length(unique(c1))==4 && length(unique(c2))==4) ...
       )
        equiv=true;
    end
end

function equiv = isEqual(c1,c2)
    equiv = (norm(c1-c2)==0);
end

function ok = checkCanonical(c)
    cc = makeCanonical(c);
    ok=true;
    for j=1:4
        for k=j+1:4
            if(isEqual(c(j),c(k))~=isEqual(cc(j),cc(k)))
                ok=false;
            end
        end
    end
    if(~ok)
        c
        cc
    end
end

function c = makeCanonical(c1)
    v = unique(c1);
    cnt = zeros(size(v));
    for i=1:length(v)
        cnt(i) = sum(c1==v(i));
    end
    [cnt,ind] = sort(cnt);
    cnt(ind) = 1:length(v);
    c = zeros(size(c1));
    for i=1:length(v)
        c(c1==v(i))=cnt(i);
    end
end
