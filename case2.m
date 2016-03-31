function mycovariance= case2(Xi,Xk)
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

% -------------------------------------------------------------------
% mycovariance is the empirical estimators for case 2
% -------------------------------------------------------------------

n = size(Xi,1);

t1_1 = mean(Xi.*Xi.*Xk.*Xk);
t1_2 = -2*mean(Xi.*Xk.*Xk)*mean(Xi);
t1_3 = -2*mean(Xi.*Xi.*Xk)*mean(Xk);
t1_4 = 4*mean(Xi.*Xk)*mean(Xi)*mean(Xk);

t2_1 = mean(Xi.*Xi)-2*mean(Xi)^2;
t2_2 = mean(Xk.*Xk)-2*mean(Xk)^2;

zeta1 = (t1_1+t1_2+t1_3+t1_4) - t2_1*t2_2;
zeta1 = (1/4)*zeta1;

%final_term
constante = 1/nchoosek(n,2);
mycovariance = (constante * (2*(n-2)*zeta1));

end