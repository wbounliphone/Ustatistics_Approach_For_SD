function mycovariance = case7(Xi)
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
% mycovariance is the empirical estimators for case 7
% -------------------------------------------------------------------

n=size(Xi,1);

zeta1_t1 = mean(Xi.*Xi.*Xi.*Xi);
zeta1_t2 = -4*mean(Xi) * mean(Xi.*Xi.*Xi);
zeta1_t3 = 4*mean(Xi.*Xi)*mean(Xi)^2;
zeta1_t4 = mean(Xi.*Xi) - 2*mean(Xi)^2;

zeta1 = zeta1_t1+zeta1_t2+zeta1_t3 - (zeta1_t4)^2;
zeta1 = (1/4)*zeta1;


constante = 1/nchoosek(n,2);
mycovariance = constante*(2*(n-2)*zeta1);

end

