%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function [b,stats,lik,result_mix]=Model3_MonteCarlo(data,weights,p_s_e_user_sample, D_1_restriction, tol, maxit, q, ENRICH, v, trials, min_lik,min_lik_coefficient_LetJoeDecide)
%                                 Model3_MonteCarlo(dat,weights, p_s_e_user_sample,D_1_restriction,10^-6, maxit, q,ENRICH, .015, MC_reps, liks(j),min_lik_coefficient_LetJoeDecide);
%% this function can be called by BothModels.m to run FASA
% nomenclature: 
%   Model 3 means FASA, Model 1 mean classic ISA
%   p represents a vector with [D0, D1, D2]
% dependencies: 
%   NewModel_MLEFIT.m
%   elongation_distribution.m
% Inputs: dat,weights, p_s_e_user_sample,D_1_restriction,10^-6, maxit, q,ENRICH, .015, MC_reps, liks(j),min_lik_coefficient_LetJoeDecide
%   dat: AUC data for a sample
%   weights: 1 for data>cutoff, 0 for data<cutoff
%   p_s_e_user_sample: user restricted value
%   D_1_restriction: if =0, FASA iterates D0, D1, and D2=1-D0-D1;
%                    if =1, FASA iterates D0, with D1, D2 calculated from D0
%   tol: criteria for convergence
%   maxit: maximum times for iteration
%   q: the natural occuring 13C percentage, =0 the script will calclute q from control, if set =0.112, then this is the value that will be used in the fitting
%   ENRICH: enrichment of 13C in the isotopic tracer. 
%           e=0.99 for 13C-U-Glucose means 99% of the carbons are 13C 
%   v: a coefficient used in iteration to adjust gradient
%   trials: same as MC_reps in BothModels.m, same as MC_successes in Metabolism_Preprocess_UCLA.m
%                     number of trials for modeling. 
%   min_lik: used in iteration
%   min_lik_coefficient_LetJoeDecide, % when running FASA after ISA, only when the new 
%           likelihood > min_lik_coefficient_LetJoeDecide*likelihood from ISA is called a success. 
%           otherwise it's marked "fail", and only 5 times of fails are allowed.
% Outputs:
%   [b,stats,lik,result_mix]
%   b={p,elong}, 
%       where p is a 3x1 vector holding [D0,D1,D2], the probabilities of a 
%           zero, single, double labeled ACoA in the cellular pool. sum(pi)=1
%       where elong is a variable length vector of [S,I,IE1,......IEi] sum(eI)=1.
%           where S is the proportion that is made entirely in the cell. 
%                 I is lipid uptaked or preexisting  
%                 IE1-IEi are the proportions that 1,2,3,.... ACoA added
%                 within cells
%%%updates:
% original scripts from Moses
% 2016.02.26 Dylan: restrict p
% 2017.01.21 Dylan: add p_s_e_user_sample
% 2018.03.26 Dylan: add "ENRICH" for Joe's D and 1-D iteration
% 20181021 Dylan: add header
%%
    elongs=(length(data)-16)/2;  %number of possible elongations
    success=0;
    fails=0;
    while(success<trials && fails<5)
        %%%%%%%%% Get Seed p %%%%%%%%%%% 
        p0=zeros(1,3);
        p0(1)=rand;
        p0(2)=rand*(1-p0(1));
        p0(3)=1-sum(p0);
        p0=[.425,.025,.55]; % seed p is not random, it starts with a big p2
        %%%  Dylan fixing p0,p1,p2 for Joe
        fix_p=0;
        p_user_sample=p_s_e_user_sample(1:3); 
            p_user_sample(2)=1-p_user_sample(1)-p_user_sample(3);
        if ( sum(p_user_sample<=1 & p_user_sample>=0) ==3 )   % yes fix p
            fix_p=1;
            p0=p_user_sample;
        end
        
        %%%%%%%%% Get Seed s_e series %%%%%%%%%%% 
        e0=zeros(1,elongs+2); 
        for i=1:length(e0)-1
            e0(i)=rand()*(1-sum(e0));
        end
        e0(length(e0))=1-sum(e0);
        
        %%% Dylan: fixing s_e series for Joe
        s_e_user_sample=p_s_e_user_sample(4:10);  
        s_e_user_sample=s_e_user_sample(1:length(e0));
        s_e_fix_index=find(s_e_user_sample<=1 & s_e_user_sample>=0);
       
        %%%%% just for display
        parameters={'p012';'s';'e0';'e1';'e2';'e3';'e4';'e5'};
        fix_index=s_e_fix_index+1; 
            if (fix_p)
                fix_index=[1,fix_index];
            end
        if (~isempty(fix_index))
            fix_pa_display=[];
            for i=1:length(fix_index)
                fix_pa_display=[fix_pa_display,', ',parameters{fix_index}];
            end
            display(['fixed parameters:',fix_pa_display]);
        end
        %%%% end of display
            if ( sum (s_e_user_sample(s_e_fix_index))<0) 
                s_e_fix_index=[];
            end
        if (~isempty(s_e_fix_index))
            e0=zeros(1,elongs+2);             
            e0(s_e_fix_index)=s_e_user_sample(s_e_fix_index);
            
            s_e_unfix_index=setdiff(1:length(e0),s_e_fix_index);
             for i= 1:(length(s_e_unfix_index)-1)
                 e0(s_e_unfix_index(i))=rand()*(1-sum(e0));
             end
             e0(s_e_unfix_index(length(s_e_unfix_index)))=1-sum(e0);
          
        end
        %%%%%%%%% Combine p and  s_e series %%%%%%%%%%% 
        b0={p0,e0};
        [b,its,lik,mix]=NewModel_MLEFIT(data,weights,fix_p,s_e_fix_index,b0,D_1_restriction,tol,maxit,q,ENRICH,v);%v=0.015
        if(lik(its+1)>= (min_lik* min_lik_coefficient_LetJoeDecide))  % likelyhood: the bigger the better
            success=success+1;
            display(['Success=',num2str(success)]);
            fails=0;
            ps(success,:)=b{its+1,1};
            es(success,:)=b{its+1,2};
            liks(success)=lik(its+1);
        else
            fails=fails+1;
        end
        
    end
    if(success>0)
        best=find(liks==max(liks));
        if(length(best>1))
            best=best(1);
        end
        
        p=ps(best,:);
        elong=es(best,:);
        result_mix=elongation_distribution(q,p,elong);
        lik=-sum(weights.*((data/sum(data)-result_mix).^2)); %based on linear SD growth with AUC measures
        b={p,elong};
        if(success==1)
            stats=[lik,p,elong;liks,ps,es;zeros(1,6+elongs);liks',ps,es];
        else
            stats=[lik,p,elong;mean(liks),mean(ps),mean(es);var(liks),var(ps),var(es);liks',ps,es];
        end
    else
        p=zeros(1,3);
        elong=zeros(1,length(e0));
        result_mix=zeros(1,length(data));
        lik=-inf;
        b={p,elong};
        stats={[],[]};
    end
end