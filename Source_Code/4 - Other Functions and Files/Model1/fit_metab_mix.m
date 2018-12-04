%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function [p,s,log_lik]=fit_metab_mix(data,q,p_range,s_range,cutoff,enrich)
%This function gets a maximum liklihood estimate for an admixture data set
%of  metabolism of C-13 into fatty acids made of 2-carbon links(
%Palmitate T=8);
%data is a 1x(2T+1+1) vector of the measured amounts of heavy carbon in Mass
%Spec (first +1 is for 0 label, seconf +1 is for extra carbon added later)
%q is the pre-determined background level of C13 in the media
%p and s are the initial estimates of the probability of using glucose for
%a building source of the fatty acid (given that it is de novo synthesized
%as opposed to scavenged), and the probability of synthesis instead of scavenging, respectfully.

    to_fit=data>cutoff;
    data=data/sum(data); %turn data into relative probabilities.

    log_lik=zeros(length(p_range),length(s_range));
    for i=1:length(p_range)
        for j=1:length(s_range)
            %log_lik(i,j)=liklihood(data,p_range(i),q,s_range(j), to_fit,enrich);
            T_f=(length(data)-2)/2;
            mix_distr_f=Label_distr_mix(p_range(i),q,s_range(j),T_f,enrich);
            log_lik(i,j)=-sum(to_fit.*((data-mix_distr_f).^2));
            
        end
    end
    t=get(gca,'Title');
    t=get(t,'String');
    mesh(s_range,p_range,log_lik)
    title(t)
    xlabel('s')
    ylabel('p')
    zlabel('Normalized log(liklihood)')
    pause(.1)
    [I,J,val]=find(log_lik==max(max(log_lik)));
    p=median(p_range(I));
    s=median(s_range(J));
    log_lik=mean(mean(log_lik(I,J)));



    function lik=liklihood(data,p,q,s,to_fit,enrich)
        %this function gives log-liklihood of the data given parameters p,q,and s
        T=(length(data)-2)/2;
        mix_distr=Label_distr_mix(p,q,s,T,enrich);
        lik=-sum(to_fit.*((data-mix_distr).^2));
    end
end