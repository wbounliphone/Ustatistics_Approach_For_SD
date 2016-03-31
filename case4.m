function mycovariance= case4(Xi,Xj,Xl)
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
% mycovariance is the empirical estimators for case 4
% -------------------------------------------------------------------

n = size(Xi,1);

%zeta1

t1_1 = mean(Xi.*Xi.*Xj.*Xl);
t1_2 = -mean(Xi)*mean(Xj.*Xi.*Xl);
t1_3 = -mean(Xj)*mean(Xi.*Xi.*Xl);
t1_4 = -mean(Xi)*mean(Xi.*Xj.*Xl);
t1_5 = mean(Xi)^2*mean(Xj.*Xl);
t1_6 = mean(Xj)*mean(Xi)*mean(Xi.*Xl);
t1_7 = -mean(Xi.*Xi.*Xj)*mean(Xl);
t1_8 = mean(Xi)*mean(Xl)*mean(Xj.*Xi);
t1_9 = mean(Xj)*mean(Xl)*mean(Xi.*Xi);

t2_1 = (mean(Xi.*Xj) - 2*mean(Xi)*mean(Xj));
t2_2 = (mean(Xi.*Xl) - 2*mean(Xi)*mean(Xl));

zeta1 = (t1_1+t1_2+t1_3+t1_4+t1_5+t1_6+t1_7+t1_8+t1_9) - t2_1*t2_2;
zeta1 = (1/4)*zeta1;

%final_term
constante = 1/nchoosek(n,2);
mycovariance = (constante * (2*(n-2)*zeta1));

end