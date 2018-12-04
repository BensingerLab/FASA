%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function pdf=elongation_distribution(q,p,elong)
%This returns the pdf of labeling given the number of ACoA units T.
%q is the background labeling percentage in the media
%p is a 3x1 vector holding [p0,p1,p2], the probabilities of a zero, single,
%     a double labeled ACoA in the cellular pool. sum(pi)=1
%elong is a variable length vector of [en,e0,e1,......ei] sum(eI)=1.  e_n is the
%proportion that is made entirely in the cell. e0 is lipid taken out of the
%media.  e1-ei are the proportions that have 1,2,3,.... ACoA added within
%the cell

T=7+length(elong)-2;

pdf=p;

for i=2:T
    pdf=add_label(pdf,p);
end
pdf=pdf*elong(1);
for i=2:length(elong)
    to_add=i-2; %number of ACoA's to add;
    L=T-to_add; %number of ACoA's taken from the media
    temp=binopdf(0:(L*2),L*2,q);
    for j=1:to_add
        temp=add_label(temp,p);
    end
    pdf=pdf+temp*elong(i);
end

pdf=add_label(pdf,[1-q,q]);    
    

end



function added=add_label(pdf,distr)
%this adds length to the pdf distribution. distr is a variable length
%vector [d0,d1,.....dn] which is the probability of adding 0,1,2....label
%to pdf
    l=length(pdf);
    dummy=zeros(1,l+length(distr)-1);
    dummy(1:l)=distr(1)*pdf;
    for i=2:length(distr)
        dummy(i:l+i-1)=dummy(i:l+i-1)+pdf*distr(i);
    end
    added=dummy;
end