function [ind_uptri,ind_ijkl,ind_qr] = matching_indices(p)
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
%ind_uptri indice of the upper triang
%ind_ijkl is a matrix with ((p*(p+1))/2)*((p*(p+1))/2) rows and 4 columns of the indice of ijkl
%ind_qr is a matrix matrix ((p*(p+1))/2 rows and 2 colums of indice of the upper triangular
% -------------------------------------------------------------------


ind_uptri = [];
rows = [1:p]'*ones(1,p);
cols = ones(p,1)*[1:p];
for i=1:p
    ind_uptri = [ind_uptri;rows(i,i:end)' cols(i,i:end)'];
end

ind_ijkl=[]; ind_qr=[];
for q=1:size(ind_uptri,1)
    for r=q:size(ind_uptri,1)
        ind_ijkl = [ind_ijkl;ind_uptri(q,:) ind_uptri(r,:)];
        ind_qr = [ind_qr; q r];
    end
end

end