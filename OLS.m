% CODE: Function OLS
%
% GOAL: OLS regression
%
% AUTH: Monika Avila M?rquez
%
% Date: 28-06-16
%-------------------------------------------------

function [beta_ols]=OLS(X,Y)

beta_ols=(X'*X)\X'*Y; 

end 


