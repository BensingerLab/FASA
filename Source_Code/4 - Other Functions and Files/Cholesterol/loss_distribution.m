%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function loss_distr=loss_distribution(distr,lambda)

n=length(distr);
loss_distr=zeros(1,n+2);
temp=distr*lambda;

loss_distr(1:n)=distr*lambda;
loss_distr(3:(n+2))=loss_distr(3:(n+2))+distr*(1-lambda);

end