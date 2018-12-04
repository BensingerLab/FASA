%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function val=fittowindow(val,window)
%window is a vector [a,b], with a<=b.  this function changes val to be in
%the closed interval [a,b] such that if val<a, val=a, and if val>b val=b;

if(val<window(1))
    val=window(1);
elseif(val>window(2))
    val=window(2);
end