%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function [Model_results, Model_colnames]=Cholesterol_Automated_Batch(data,cutoff,num_controls,SaveFigs,q,lambda,user_p,ENRICH)
%This is an automated version of the Cholesterol Model for Mass-Spec Stastistical
%Fitting. . 

%%%%%%%%%%%%% revision history %%%%%%%%%%%%%%%%%
% dylan
% 2015.06, added SaveFigs, presented qs, 
% 2015.06, revise things for lambda and control
% 2015.08, to fix p=0.8 at local maxima, I added parameter user_p=[0,0.5]
% try to add try catch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%it should be called as:
% Cholesterol_Automated(data,limits,num_controls,SaveFile,SaveFigs, enrich)

%   The "data" variable should be a mX35 matrix of AUC data for the cholestorol, 
%   where m is the number samples measured(including controls).
%   The 35 refers to: (30 for M-2,M-1, to M+27, 2 for possible Si species, 
%   and 3 for possible extra carbon weight.
%   ****************************************************************
%   *****Control measurements MUST be the last rows of data*****
%*  ***************************************************************

%   "cutoff" is the AUC under which we ignore data.  If no extra knowledge
%   about the system is present, cutoff should be set as -1.

%   "num_controls" is the number of control runs that are present in the
%   "data" matrix

%   "SaveFigs" is a binary 0 or 1.  When set to 1, fitting results and
%   final grid search results are saved as .fig files in a subdirectory or
%   "saveFile" - not currently implemented

%   "enrich" is an optional parameter to set the carbon enrichment of the
%   labled glucose. Default enrichment is 98.4%

    %%%%%%%%%%Set Save directory%%%%%%%%%%%%%%%%%
% dylan, disable the figure path
    %if(SaveFigs)
    %    figdir=[SaveFile,' - Figures'];
     %   if(~exist(figdir,'dir'))
      %      mkdir(figdir);
       % end
    %end
    %%%%%%%%Bookkeeping%%%%%%%%%%%%%%%%%%%%%
    scrsz = get(0,'ScreenSize');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fitmix=zeros(size(data));
    [runs,isops]=size(data);
    runs=runs-num_controls;
    clear temp isops;
    
        %%%%%%%%%%%Set Enrichment Level%%%%%%%%%%%%%%%%
%{
    if isempty(varargin)
        enrich=.984;  %Default level
    else
        enrich=varargin{1};  %Specified level
    end
%}
        enrich=ENRICH;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %Fit background C13 levels
    if(isempty(lambda))
        lambda=0
    end
    
    for j=1:num_controls
        tofit=data(j+runs,:)>cutoff;
        if(q==0 | lambda==0 ) % if [q,lambda] are not both decided
            [qs(j),lambdas(j),liks(j+runs),fitmix(j+runs,:)]=fit_cholest_back_auto(data(j+runs,:),cutoff);
        end
        if (q>0)
            qs(j)=q;
        end
        if (lambda>0)
            lambdas(j)=lambda;
        end
        tempfit=Cholest_distr(0,qs(j),0); %function pdf_mix=Cholest_distr(p,q,s,varargin)
        fitmix(runs+j,:)=loss_distribution(tempfit,lambdas(j)); %function loss_distr=loss_distribution(distr,lambda)
    end
    
    if(num_controls==0)
        if (q>0)
            qs=q;
        else
            q=0.0116;
            qs=q;
        end 
        
        if (lambda>0)
            lambdas=lambda;
        else
            lambda=0.023;
            lambdas=lambda;
        end
        tempfit=Cholest_distr(0,q,0); 
        fitmix(runs+1,:)=loss_distribution(tempfit,lambda);
        num_controls=1;
    end
    
    q=mean(qs);
    lambda=mean(lambdas);
    display(['Background C13 levels for chol:']);
    display(qs); %dylan
    display(lambdas);
    fprintf('q=%.6f \n', q);
    fprintf('lambda=%.6f \n', lambda);
       
    h=gcf();
    close(h);  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    h=figure('name','Grid Search','Position',[scrsz(3)-scrsz(4)/2-10 scrsz(4)/2-80 scrsz(4)/2 scrsz(4)/2]);
    h2=figure('name','Mixture Fitting', 'Position',[10 scrsz(4)/2-80 scrsz(3)/2 scrsz(4)/2]);

    %%%%%%%%%%Fit all cell lines
    for j=1:(runs) 
        dat=data(j,:);
        
        figure(h)
        title(['Cholesterol - Plate',num2str(j)])       
        [ps(j),ss(j),liks(j),fitmix(j,:)]=fit_lossCholest_mix_auto(dat,q,lambda,cutoff,user_p,enrich);
       
        % display on the screen and save data
        display(['plate',num2str(j),': p=',num2str(ps(j)),'; s=',num2str(ss(j)),'; likelihood=',num2str(liks(j))]);
        Results_temp=[[q*ones(1,j)]',[lambda*ones(1,j)]',ps',ss',liks(1:j)',fitmix(1:j,:)];
        save('temp_result.mat','Results_temp')
        
        
        if(SaveFigs)
            saveas(h,[figdir,'\Cholesterol - Plate',num2str(j),' - Grid Search.fig']);
        end

        figure(h2)
        plotmix(fitmix(j,:),dat)
        title(['Cholesterol - Plate',num2str(j)])
        if(SaveFigs)
            saveas(h2,[figdir,'\Cholesterol - Plate',num2str(j),' - Mixture Fit.fig']);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%Add empty entries for excel buffer%%%%%%%%%%
    ps=[ps,zeros(1,num_controls)];
    ss=[ss,zeros(1,num_controls)];
    liks=[liks(1:runs),zeros(1,num_controls)];
    
    %%%%%%%%%collect data to write to excel%%%%%%%%%%%%
    Model_results=[[q*ones(1,runs),qs]',[lambda*ones(1,runs),lambdas]',ps',ss',liks',fitmix];
    Model_colnames={'q','lambda','p','s','log_lik','M-2','M-1','M+0'};
    for j=1:32;
        Model_colnames{length(Model_colnames)+1}=['M+',num2str(j)];
    end
    
    
end

function plotmix(fitmix,dat)
    bar([fitmix',dat'/sum(dat)]);
    ylabel('AUC Percentage')
    xlabel('Extra Mass Units')
    set(gca,'XTick',1:length(dat))
    labels={'0'};
    for i=1:length(dat)-1
        labels{i+1}=num2str(i);
    end
    set(gca,'XTickLabel',labels) 
end





