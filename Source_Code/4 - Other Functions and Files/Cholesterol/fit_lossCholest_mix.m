%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function [p,s,log_lik]=fit_lossCholest_mix(data,q,lambda,p_range,s_range,cutoff,enrich)
%This function gets a maximum liklihood estimate for an admixture data set
%of  metabolism of C-13 into cholesterol;
%data is a 1x(27+1+5) vector of the measured amounts of heavy carbon in Mass
%Spec (first +1 is for 0 label, seconf +5 is for extra weight added later)
%q is the pre-determined background level of C13 in the media
%p and s are the initial estimates of the probability of using glucose for
%a building source of the fatty acid (given that it is de novo synthesized
%as opposed to scavenged), and the probability of synthesis instead of scavenging, respectfully.

    to_fit=data>cutoff;
    data=data/sum(data); %turn data into relative probabilities.

    log_lik=zeros(length(p_range),length(s_range));
    for i=1:length(p_range)
        for j=1:length(s_range)
            log_lik(i,j)=liklihood(data,p_range(i),q,lambda,s_range(j),to_fit,enrich);
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



    function lik=liklihood(data,p,q,lambda,s,to_fit,enrich)
        %this function gives log-liklihood of the data given parameters p,q,and s
        mix_distr=Cholest_distr(p,q,s,enrich);
        mix_distr=loss_distribution(mix_distr,lambda);
        lik=-sum(to_fit(3:length(to_fit)).*((data(3:length(data))-mix_distr(3:length(mix_distr))).^2)); %based on linear SD growth with AUC measures
       
    end

end