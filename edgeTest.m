function [thresh_trace,thresh_eig,rejectNulltrace,rejectNulleig]=edgeTest(X,delta,mu)
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
% ----------------------------------------------------------------
% rejectNull is a binary matrix saying whether we reject null hypothesis that the true precision matrix is zero at the ijth entry
% threshtrace is a scalar that is the threshold for the trace test
% ----------------------------------------------------------------

%unbiased estimator of Cov
sigma_hat = ustat_estimator_of_covariance(X); %etape 1 
theta_hat = inv(sigma_hat); % etape 2 
testStatistic = theta_hat;


%estimator of the covariance matrix Cov: 
covCov = covCov_estimator(X); %etape 3


%eigenvalue of empirical cov matrix
a_hat = eig(sigma_hat);
a_hat = sort(a_hat,'descend');


%threshold eig
max_l_covcov = max(eig(covCov)); % etape 4
errCeig = sqrt(2*max_l_covcov)*norminv(1-(delta/2),0,1); %etape 5 
%assert(errCeig<min(a_hat))
t_eig =  -errCeig ./ (a_hat .* (a_hat - errCeig)) ; %etape 6
t_eig = norm(t_eig,2);
thresh_eig = mu*t_eig;


%threshold trace
thetrace = trace(covCov); % etape 4
errCtrace = sqrt(2*thetrace)*norminv(1-(delta/2),0,1); %etape 5 
t_trace =  -errCtrace ./ (a_hat .* (a_hat - errCtrace)) ;        
t_trace = norm(t_trace,2);
thresh_trace = mu*t_trace;


assert(errCtrace>errCeig)

if min(a_hat)<errCtrace
    warning('insufficient samples to find any significant dependence using trace method')
    thresh_trace = Inf;
end


if min(a_hat)<errCeig
    warning('insufficient samples to find any significant dependence using eig method')
    thresh_eig= Inf;
end

assert(thresh_trace>thresh_eig || isinf(thresh_trace))

rejectNulltrace = abs(testStatistic)>thresh_trace; % reject H0
rejectNulleig = abs(testStatistic)>thresh_eig; % reject H0

end




function u_stat = ustat_estimator_of_covariance(X,Y);

m = size(X,1);
if(nargin==1)
    X = X-ones(m,1)*mean(X);
    u_stat=(X'*X)./(m-1);
    return
end
u_stat = dot(X-mean(X),Y-mean(Y))./(m-1);

end


