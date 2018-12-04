%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function [q,lambda,log_lik]=fit_Cholest_background(data,q_range,lambda_range,cutoff)
%This function gets a maximum liklihood estimate for an admixture data set
%of  metabolism of C-13 into cholesterol;
%data is a 1x(27+1+5+2) vector of the measured amounts of heavy carbon in Mass
%Spec (first +1 is for 0 label, seconf +5 is for extra weight added later,+2 is for loss)
%q is the pre-determined background level of C13 in the media
%p and s are the initial estimates of the probability of using glucose for
%a building source of the fatty acid (given that it is de novo synthesized
%as opposed to scavenged), and the probability of synthesis instead of scavenging, respectfully.

    to_fit=data>cutoff;
    data=data/sum(data); %turn data into relative probabilities.

    log_lik=zeros(length(q_range),length(lambda_range));
    for i=1:length(q_range)
        for j=1:length(lambda_range) 
           log_lik(i,j)=liklihood(data,q_range(i),lambda_range(j),to_fit);
        end
    end

    t=get(gca,'Title');
    t=get(t,'String');
    mesh(lambda_range,q_range,log_lik)
    title(t)
    xlabel('lambda')
    ylabel('q')
    zlabel('Normalized log(liklihood)')
    pause(.1)

    [I,J,val]=find(log_lik==max(max(log_lik)));

    q=mean(q_range(I));
    lambda=mean(lambda_range(J));
    log_lik=mean(log_lik(I,J));

end

    function lik=liklihood(data,q,lambda,to_fit)
        %this function gives log-liklihood of the data given parameterq
        mix_distr=Cholest_distr(0,q,0);
        mix_distr=loss_distribution(mix_distr,lambda);
       lik=-sum(to_fit.*((data-mix_distr).^2)); %based on linear SD growth with AUC measures
       
    end
