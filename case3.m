function mycovariance= case3(Xi,Xk,Xl)
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
% mycovariance is the empirical estimators for case 3
% -------------------------------------------------------------------
n = size(Xi,1);

%zeta1
zeta1_t1 = mean(Xi.*Xi.*Xk.*Xl);
zeta1_t2 = -2*mean(Xi)*mean(Xi.*Xk.*Xl);
zeta1_t3 = -mean(Xk)*mean(Xi.*Xi.*Xl);
zeta1_t4 = 2*mean(Xi.*Xl)*mean(Xi)*mean(Xk);
zeta1_t5 = -mean(Xi.*Xi.*Xk)*mean(Xl);
zeta1_t6 = 2*mean(Xi.*Xk)*mean(Xi)*mean(Xl);

zeta1_t7 = mean(Xi.*Xi) - 2*mean(Xi)^2;
zeta1_t8 = mean(Xk.*Xl) - 2*mean(Xk)*mean(Xl);

zeta1 = zeta1_t1+zeta1_t2+zeta1_t3+zeta1_t4+zeta1_t5+zeta1_t6 - (zeta1_t7*zeta1_t8);
zeta1 = (1/4)*zeta1;


%final_term
constante = 1/nchoosek(n,2);
mycovariance = (constante * (2*(n-2)*zeta1));

end