function mycovariance = case5(X,Y)
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
% mycovariance is the empirical estimators for case 5
% -------------------------------------------------------------------

n =length(X);
meanX = mean(X);
meanY = mean(Y);

%zeta1
% E[(h_1)^2]
zeta1_t1 = mean((X.*Y).^2);
zeta1_t2 = -2*meanX * mean(X.*Y.*Y);
zeta1_t3 = meanX^2*mean(Y.*Y);
zeta1_t4 = -2*meanY * mean(X.*X.*Y);
zeta1_t5 = 2*meanX*meanY*mean(X.*Y);
zeta1_t6 = meanY^2*mean(X.*X);


% E[(h_1)]^2
zeta1_t7 = mean(X.*Y);
zeta1_t8 = -2*meanY*meanX;

zeta1 = zeta1_t1+zeta1_t2+zeta1_t3+zeta1_t4+zeta1_t5+zeta1_t6 - (zeta1_t7+zeta1_t8)^2;
zeta1 = (1/4)*zeta1;

constante = 1/nchoosek(n,2);
mycovariance = constante*(2*(n-2)*zeta1);


end
