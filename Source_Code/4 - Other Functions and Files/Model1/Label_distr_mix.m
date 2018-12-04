%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function mix_distr=Label_distr_mix(p,q,s,T,enrich)
%% dylan: mix_distr=s*lab_distr+(1-s)*bin_distr;
%% lab_distr is the 
    tot=T*2;
    q2=q; %%%%%%%for adding in new background carbons
    p_2=p*(enrich^2)+(1-p)*q^2;
    p_1=2*q*(1-q)*(1-p)+2*p*enrich*(1-enrich);
    lab_distr=Label_distr(p_1,p_2,T);  %this gives distribution for labeled T units

    %%%%ADD IN 1 RANDOM BACKGROUND CARBON%%%%%%%
    temp_distr=lab_distr(1)*(1-q2); %
    for i=2:(tot+1)
        temp_distr(i)=lab_distr(i-1)*q2+lab_distr(i)*(1-q2);
    end
    temp_distr(tot+2)=lab_distr(tot+1)*q2;
    lab_distr=temp_distr;
    %%%%%%%%%%%%%%%%%%%%%%%%

    bin_distr=zeros(1,tot+1);
    for i=0:tot
        bin_distr(i+1)=(q^i)*((1-q)^(tot-i))*nchoosek(tot,i);
    end

    %%%%ADD IN 1 RANDOM BACKGROUND CARBON%%%%%%%
    temp_distr=bin_distr(1)*(1-q2); %
    for i=2:(tot+1)
        temp_distr(i)=bin_distr(i-1)*q2+bin_distr(i)*(1-q2);
    end
    temp_distr(tot+2)=bin_distr(tot+1)*q2;
    bin_distr=temp_distr;
    %%%%%%%%%%%%%%%%%%%%%%%%

    mix_distr=s*lab_distr+(1-s)*bin_distr;

end
