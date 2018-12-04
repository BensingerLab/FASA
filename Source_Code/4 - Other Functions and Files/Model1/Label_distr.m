%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function pdf=Label_distr(p_1,p_2,T)

%This function returns the pdf of the correlated binomial.
%It assumes you build a molecule of T 2-carbon building blocks.
%p_1 and p_2 are the probabilties of getting a block with 1 or 2 labeled
%carbons, respectively.


pdf=zeros(1,2*T+1);
p_0=1-p_1-p_2;
for x=0:2*T;
    lim=floor(x/2); %only sum over elements where you can use double labels.
    for i=0:lim
        temp=nchoosek(T,i)*(p_2^i); %probability of picking i double-labeled blocks, we sum up to total possible for a given xnd
        if((x-2*i)>(T-i))
            temp2=0;  %if there is no way to get the remaining x-2i labels in T-i slots using single-labeled carbons, set p=0;
        else
            temp2=nchoosek(T-i,x-2*i)*(p_1^(x-2*i))*(p_0^(T-x+i)); %otherwise, probability of choosing x-2i single labels
        end
        pdf(x+1)=pdf(x+1)+temp*temp2;  %sum up over all i P(using i double-lab)*P(using x-2i single labels)
    end
end
    
        