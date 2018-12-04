%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function pdf_mix=Cholest_distr(p,q,s,varargin)

%This function returns the pdf of the labeled cholesterol
%It assumes that Cholesterol is built from 10 A-CoA's and 7 single carbons
%from the A-CoA pool
%It then adds 3carbons and a silicon for the GC-MS tube.
%There are 32 possible heavy states (30 for all carbons, 2 for the two
%silicon species) and one state with no labels

if(isempty(varargin))
    enrich=.984;%%%enrichment fraction of glucose carbons
else
    enrich=varargin{1};
end
pdf=zeros(1,33);
q2=q; %%%%%%%for adding in new background carbons

s0=.9223; %%Silicon labeling frequencies
s1=.0468;
s2=.0309;

p_2=p*(enrich^2)+(1-p)*q^2; %Prob (A-CoA has 2 C13s)
p_1=2*q*(1-q)*(1-p)+2*p*enrich*(1-enrich); %Prob(A-CoA has 1 C13)
pdf_20=Label_distr(p_1,p_2,10); %This makes the pdf of the 10 A-CoA's (20 carbons)

p_c13=p*enrich+(1-p)*q;  %Probability that a random Carbon in A-CoA pool is heavy
pdf_7=binopdf(0:7,7,p_c13); %pdf of weights added in by single Carbons

%%%Mix pdf_7 and pdf_20 together%%%%%
pdf_27=zeros(1,28);
for i=1:21
    for j=1:8
        pdf_27(i+j-1)=pdf_27(i+j-1)+pdf_20(i)*pdf_7(j);
    end
end

pdf_5=pdf_sil(q2,s0,s1,s2);%get pdf of silicon/carbon unit added on for GC/MS

%%%Mix pdf_5 and pdf_27 together

for i=1:28
    for j=1:6
        pdf(i+j-1)=pdf(i+j-1)+pdf_27(i)*pdf_5(j);
    end
end

%%%%%%%%%%%%Make Bin_Distr%%%%%%%%%%%%
bin_distr=zeros(1,33);
temp_distr=binopdf(0:27,27,q);
for i=1:28
    for j=1:6
        bin_distr(i+j-1)=bin_distr(i+j-1)+temp_distr(i)*pdf_5(j);
    end
end

pdf_mix=s*pdf+(1-s)*bin_distr;

end

function pdf_5=pdf_sil(q2,s0,s1,s2)
    pdf_5=zeros(1,6);
    pdf_5(1)=s0*(1-q2)^3;
    pdf_5(2)=s1*(1-q2)^3 + 3*s0*q2*(1-q2)^2;
    pdf_5(3)=s2*(1-q2)^3 + 3*s1*q2*(1-q2)^2+ 3*s0*(1-q2)*q2^2;
    pdf_5(4)=3*s2*q2*(1-q2)^2+3*s1*(1-q2)*q2^2+s0*q2^3;
    pdf_5(5)=3*s2*(1-q2)*q2^2+s1*q2^3;
    pdf_5(6)=s2*q2^3;
end
    
    
    
    
