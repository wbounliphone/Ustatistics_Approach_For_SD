function covCov = covCov_estimator(X)
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
% -------------------------------------------------------------------
% covCov  is the full covariance between the elements of \hat_{Sigma}
% covCov = cov[\hat_{Sigma}]
% -------------------------------------------------------------------

[~,ind_ijkl,ind_qr] = matching_indices(size(X,2));
[ind_c1,ind_c2,ind_c3,ind_c4,ind_c5,ind_c6,ind_c7]=indice_for_all_cases(ind_ijkl);


%%%%% CASE 1: 
for indice=1:length(ind_c1)
    Xi = X(:,ind_ijkl(ind_c1(indice),1));
    Xj = X(:,ind_ijkl(ind_c1(indice),2));
    Xk = X(:,ind_ijkl(ind_c1(indice),3));
    Xl = X(:,ind_ijkl(ind_c1(indice),4));
    covCovTheo(ind_c1(indice)) = case1(Xi,Xj,Xk,Xl);
end

%%%%% CASE 2: 
for indice=1:length(ind_c2)
    case2_a=unique(ind_ijkl(ind_c2(indice),:));
    assert(length(case2_a)==2)
    case2_l(1) = length(find(ind_ijkl(ind_c2(indice),:)==case2_a(1)));
    case2_l(2) = length(find(ind_ijkl(ind_c2(indice),:)==case2_a(2)));
    assert(max(case2_l)==2);
    case2_ind = zeros(2,1)*NaN;
    case2_ind(1) = case2_a(1);
    case2_ind(2) = case2_a(2);
    Xi = X(:,case2_ind(1));
    Xk = X(:,case2_ind(2));
    covCovTheo(ind_c2(indice)) = case2(Xi,Xk);
end


%%%%% CASE 3: 
for indice=1:length(ind_c3)
    case3_a=unique(ind_ijkl(ind_c3(indice),:));
    assert(length(case3_a)==3)
    case3_l(1) = length(find(ind_ijkl(ind_c3(indice),:)==case3_a(1)));
    case3_l(2) = length(find(ind_ijkl(ind_c3(indice),:)==case3_a(2)));
    case3_l(3) = length(find(ind_ijkl(ind_c3(indice),:)==case3_a(3)));
    assert(max(case3_l)==2);
    case3_ind = zeros(3,1)*NaN;
    case3_ind(1) = case3_a(case3_l==2);
    case3_ind(2:3) = case3_a(case3_l==1);
    Xi = X(:,case3_ind(1));
    Xk = X(:,case3_ind(2));
    Xl = X(:,case3_ind(3));
    covCovTheo(ind_c3(indice)) = case3(Xi,Xk,Xl);
end

%%%%% CASE 4: 
for indice=1:length(ind_c4)
    case4_a=unique(ind_ijkl(ind_c4(indice),:));
    assert(length(case4_a)==3)
    case4_l(1) = length(find(ind_ijkl(ind_c4(indice),:)==case4_a(1)));
    case4_l(2) = length(find(ind_ijkl(ind_c4(indice),:)==case4_a(2)));
    case4_l(3) = length(find(ind_ijkl(ind_c4(indice),:)==case4_a(3)));
    assert(max(case3_l)==2);
    case4_ind = zeros(3,1)*NaN;
    case4_ind(1) = case4_a(case4_l==2);
    case4_ind(2:3) = case4_a(case4_l==1);
    Xi = X(:,case4_ind(1));
    Xj = X(:,case4_ind(2));
    Xl = X(:,case4_ind(3));
    covCovTheo(ind_c4(indice)) = case4(Xi,Xj,Xl);
end



%%%%% CASE 5
for indice=1:length(ind_c5)
    Xi = X(:,ind_ijkl(ind_c5(indice),1));
    Xj = X(:,ind_ijkl(ind_c5(indice),2));
    covCovTheo(ind_c5(indice)) = case5(Xi,Xj);
end



%%%%% CASE 6
for indice=1:length(ind_c6)
    case6_a=unique(ind_ijkl(ind_c6(indice),:));
    assert(length(case6_a)==2)
    case6_l(1) = length(find(ind_ijkl(ind_c6(indice),:)==case6_a(1)));
    case6_l(2) = length(find(ind_ijkl(ind_c6(indice),:)==case6_a(2)));
    assert(max(case6_l)==3);
    case6_ind = zeros(2,1)*NaN;
    case6_ind(1) = case6_a(case6_l==3);
    case6_ind(2) = case6_a(case6_l==1);
    Xi = X(:,case6_ind(1));
    Xl = X(:,case6_ind(2));
    covCovTheo(ind_c6(indice)) = case6(Xi,Xl);
end

%%%%% CASE 7
for indice=1:length(ind_c7)
    Xi = X(:,ind_ijkl(ind_c7(indice),1));
    covCovTheo(ind_c7(indice)) = case7(Xi);
end



%transforme into a matrix
for k = 1:size(ind_qr,1)
    i = ind_qr(k,1);
    j = ind_qr(k,2);
    covCov(i,j) = covCovTheo(k);
    covCov(j,i) = covCovTheo(k);
end

% cannot have a negative variance
assert(min(diag(covCov))>0)

% should be a positive definite matrix
assert(min(eig(covCov))>=0)

% trace > max eigenvalue
thetrace = trace(covCov); 
max_l_covcov = max(eig(covCov));
assert(thetrace>max_l_covcov)

end

