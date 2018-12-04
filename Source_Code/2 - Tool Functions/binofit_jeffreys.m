%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function [phat, ci]=binofit_jeffreys(x,alpha)
%this function returns a MLE of p for a binomial distribution represented
%by the vector x.
       %X is a 1x(n+1) vector with the number of trials with i successes in
       %element X(i+1);
       %alpha is the signficancelevel in the CI (ie, alpha=.05 is a 95% CI)
       %phat is the MLE for the success probability p
       %ci is the confidence interval computed using the equal-tailed
       %Jeffreys prior
       
       dims=size(x);
       if(min(dims)~=1)
           error('X must be a vector.')
       end
       n=length(x)-1;
       phat=sum(x.*(0:n))/(sum(x)*n);  %MLE of p
       
       xhat=phat*n; %expected successes in a single trial
       
       k=norminv(1-alpha/2);
       wl=(k*sqrt(4*phat*(1-phat))/n +(k^2-3)/(6*n^2))/(4*phat*(1-phat));
       wl=wl+(.5-phat)*(phat*(1-phat)*(k^2+2)-1/n)/(6*n*((phat*(1-phat))^2));
       wu=(-k*sqrt(4*phat*(1-phat))/n +(k^2-3)/(6*n^2))/(4*phat*(1-phat));
       wu=wu+(.5-phat)*(phat*(1-phat)*(k^2+2)-1/n)/(6*n*((phat*(1-phat))^2));
       
       ci_l=(xhat+.5)/(n+1+(n-xhat+.5)*(exp(2*wl)-1));
       ci_u=(xhat+.5)/(n+1+(n-xhat+.5)*(exp(2*wu)-1));
       
       ci=[ci_l,ci_u];
end