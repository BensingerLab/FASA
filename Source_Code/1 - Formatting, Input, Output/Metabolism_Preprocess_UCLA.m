%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function [str_Mod1_R, str_Mod1_C, str_Mod3_R, str_Mod3_C, str_MC_Stats]=Metabolism_Preprocess_UCLA(datafile,cutoff,q, MC_successes,sheets_2_compute, elongation, D_1_restriction,SaveFigs,SaveMCStats,lambda,user_p,min_lik_coefficient_LetJoeDecide, varargin)
%% this function can be called by FASA_Step2.m to run FASA
% dependency: BothModels.m
% nomenclature notes: 
%   Model 3 means FASA, Model 1 mean classic ISA
%   p represents a vector with [D0, D1, D2]
% Inputs: 
%	datafile, name of .xls file with FA AUC data
%	cutoff, threshold for AUC data used in FASA. 
%           Any AUC value < cutoff will not contribute to the cost function or the fitting.
%	q: the natural occuring 13C percentage, =0 the script will calclute q from control, if set =0.112, then this is the value that will be used in the fitting
%   MC_successes: number of trials for modeling. 
%                     the modeled result with SSE closest to 0 will be
%                     written in the final result.xls file
%	sheets_2_compute,same as sheets_to_model from FASA_Step2.m: 
%                   the sheetsyou want to run. e.g: =1:13 for all 13 sheets, [1,3,4,7] for selected sheets
%	elongation: elongation=0 only runs Model 1
%               elongation=1 runs both Model 3 and Model 1, 
%                     use the likelihood from Model 1 as a termination criteira for Model 3 iterations.
%               elongation=2 runs FASA only, 
%                           if there's no user setup p_s_e, it will set initial p =[.425,.025,.55], and randomize inital e and s
%	D_1_restriction: if =0, FASA iterates D0, D1, and D2=1-D0-D1;
%                    if =1, FASA iterates D0, with D1, D2 calculated from D0
%	SaveFigs, if =1, save figurs to the folder
%	SaveMCStats, if =1, stats of trials will be saved to the folder
%	lambda, only useful for modeling cholesterol. 
%           ratio of fragmentation for cholesterol data. if lambda=0, scripts will either calculate from unlabel controls or use default=0.023
%	user_p,
%	min_lik_coefficient_LetJoeDecide, % when running FASA after ISA, only when the new 
%           likelihood > min_lik_coefficient_LetJoeDecide*likelihood from ISA is called a success. 
%           otherwise it's marked "fail", and only 5 times of fails are allowed.
%	varargin, if there's a value, it'll be taken as enrichment of 13C carbons in the isotope tracer
% Outputs: 
%   [str_Mod1_R, str_Mod1_C, str_Mod3_R, str_Mod3_C, str_MC_Stats] 
%       same as [Model1_results, Model1_colnames, Model3_results,
%       Model3_colnames] in BothModels.m
%   =====files generated=====
%   Result-Model3..xls: 
%       every sheet is for a FAtty acid, s, ie, D series and new distribution generated
%   MCStats-x160... xls:
%       there will be a .xlsx per FAtty acid, where each sheet is for one
%       sample, stats summarizes caculated parameters for each MonteCarlo
%       iteration
%   Figures: comparison of origianl distribution and modeled distribution
%%%%updates
% original codes written by Moses 
% 2016.04.12 by dylan: adding p_source, so D0/1/2 can be restricted
% 2015.06.26 by dylan: adding a lot more parameters
% 2015.08.25 by dylan: adding user_p
% 20180809 by Joe (JA): update xls input parameter titles
% 20181021 Dylan: add header
% 20181118 by Joe (JA): update MCStats ouput column titles
%%
% we ignore data below cutoff, it is a nX1 vector for the n species we are
%looking at
    StartRow=1;
 %dylan   
%dylan    dir=uigetdir('','Select Folder for Save Data');
%dylan    cd(dir);
    SaveFile1=['Results-Model1-',datafile];
    SaveFile3=['Results-FASA-',datafile]; %181027 JA updated file name to reflect current nomenclature
    %SaveFile3=['Results-Model3-',datafile];
    %figdir1=[SaveFile1,' - Figures']; %181027 JA inactivated these lines to stop creation of unused and poorly-named folder
    %if(~exist(figdir1,'dir'))
    %        mkdir(figdir1);
    %end

    display('Reading in data...')
    all=importdata(datafile); 
    lipids=fieldnames(all.data);  %Read in names of all sheets
    
    if(length(cutoff)<length(lipids));  %create buffer of cutoff is not enough are passed to the function
        cutoff(length(cutoff):length(lipids))=cutoff(length(cutoff));   %dylan: they will be all the same!
    end
    
    if(length(MC_successes)<length(lipids)); %create buffer of MC_successes is not enough are passed to the function
        MC_successes(length(MC_successes):length(lipids))=MC_successes(length(MC_successes));
    end
    
    if isempty(varargin); % 181027 JA adds semicolons to these few lines so ENRICH is not displayed.
        ENRICH=.99;  %Default level % dylan they are replaced by enrich=0.99 outside, so never mind!
    else
        ENRICH=varargin{1};
    end

    
    
    %dylan: read in all the data files, identify the controls, and write 
    %the colnames and rownames to the outputfiles
    % this loop is for reading data and write blank files.
    for i= sheets_2_compute  %length(lipids) 
        eval(['temp=all.data.',lipids{i},';']); %dylan: smart way in matlab to execute sentences      
        
        %%%%% Dylan 2016.2.26 
        eval(['colnames=all.textdata.',lipids{i},'(1,:);']); %dylan: get colnames
        keep_col=find(~isnan(temp(1,:)))+1; 
        colnames=colnames(keep_col);                
        temp=temp(:,~isnan(temp(1,:))); % dylan: get rid of the columns that has no value
        
        code_col=strmatch('code',lower(colnames));
        code=temp(:,code_col);  % Dylan: get code according to header
        %%%%%%%%%%% 2016.4.12
        
        %code=temp(:,length(temp(1,:))); 
        
 
        numcontrols=sum(code==0); %
        numused=sum(code>=0);    
        %dylan: give the bothmodel function samplenames
        eval(['samplename=all.textdata.',lipids{i},';']);
        tempcode=[0;code];
        samplename=samplename(tempcode==1,1);
        %%dylan end of get samplenames
        
        temp_select=[temp(code==1,:); temp(code==0,:)];
        data=temp_select(:,1:(code_col-1));   % data is imported till the left of column of code, so you can do whatever beyond that haha
        %p_source=temp_select(:,strmatch('p_source',lower(colnames)));
        % read in p_s_e_user data %dylan 20170121
        % 180809 JA - Updated first input of strmatch to match the parameter nomenclature in the paper (used for xls input and output data).
        % 180809 JA - Because 'lower(colnames)' is used in the second input of strmatch, I used the lowercase version of the parameter nomenclature in the paper for the first input.
        % 180809 JA - Updated strmatch to contain a third input 'exact' to profent a possible issue with 'I' finding 'IEX' as well. 
        % 180809 JA - Note that internally, the original parameter nomenclature is still used.
        p0_user=temp_select(:,strmatch('d0',lower(colnames),'exact'));
        p1_user=temp_select(:,strmatch('d1',lower(colnames),'exact'));
        p2_user=temp_select(:,strmatch('d2',lower(colnames),'exact'));
        s_user=temp_select(:,strmatch('s',lower(colnames),'exact'));
        e0_user=temp_select(:,strmatch('i',lower(colnames),'exact'));
        e1_user=temp_select(:,strmatch('ie1',lower(colnames),'exact'));
        e2_user=temp_select(:,strmatch('ie2',lower(colnames),'exact'));
        e3_user=temp_select(:,strmatch('ie3',lower(colnames),'exact'));
        e4_user=temp_select(:,strmatch('ie4',lower(colnames),'exact'));
        e5_user=temp_select(:,strmatch('ie5',lower(colnames),'exact'));
        if (isempty(p0_user))
           p0_user=zeros(numused,1)-1;
        end 
        if (isempty(p1_user))
           p1_user=zeros(numused,1)-1;
        end
        if (isempty(p2_user))
           p2_user=zeros(numused,1)-1;
        end
        if (isempty(s_user))
            s_user=zeros(numused,1)-1;
        end
        if (isempty(e0_user))
            e0_user=zeros(numused,1)-1;
        end
        if (isempty(e1_user))
            e1_user=zeros(numused,1)-1;
        end
        if (isempty(e2_user))
            e2_user=zeros(numused,1)-1;
        end
        if (isempty(e3_user))
            e3_user=zeros(numused,1)-1;
        end
        if (isempty(e4_user))
            e4_user=zeros(numused,1)-1;
        end
        if (isempty(e5_user))
            e5_user=zeros(numused,1)-1;
        end 
        p_s_e_user{i}=[p0_user,p1_user,p2_user,s_user,e0_user,e1_user,e2_user,e3_user,e4_user,e5_user];
        
        %if(numcontrols==0)
         %   data(numused+1,:)=binopdf(0:cols(data)-1,cols(data)-1,q);  %dylan: create a control if there is no control, use a default q value=0.0116. this is just for fatty acid! act there's no point in create a new distr
          %  numcontrols=1;
        %end
        % ok, now numcontrols can be 0
        controls{i}=numcontrols;
        datas{i}=data;
        codes{i}=code;
        samplenames{i}=samplename; %dylan 
      
        eval(['text=all.textdata.',lipids{i},';']);
        text(:,2:cols(text))=[];    %dylan:change '' to []
        text{1}=['''',text{1}]; %dylan:the corner
        plates{i}=text;
        warning('off','all');
        if ( elongation == 0 || elongation ==1 )
            xlswrite(SaveFile1,plates{i},lipids{i},['A',num2str(StartRow)]);
        end
        if (elongation>0)
            xlswrite(SaveFile3,plates{i},lipids{i},['A',num2str(StartRow)]);
        end
        warning('on','all');
    end 

    
    for i=sheets_2_compute
        %%assignin(ws, 'lipid',lipid{i} ); %dylan:  return the which lipid is being processed to the workspace
        temp=lower(lipids{i});
        display(['Computing fits for lipid: ', temp, '...'])
        if(~strcmp(temp(1:4),'chol')) %dylan: if we are looking at Fatty Acids
         %chi-chi:[Model1_results, Model1_colnames, Model3_results, Model3_colnames,MC_stats ] = BothModels(samplename,data,cutoff,num_controls,lipid,elongation,,p_s_e_user, MC_reps,min_lik_coefficient_LetJoeDecide,SaveFigs,q,varargin) % dylan added user_p parameter
         [Mod1_R, Mod1_C, Mod3_R, Mod3_C, MC_Stats] = BothModels(samplenames{i},datas{i},cutoff(i),controls{i},lipids{i},elongation,p_s_e_user{i},D_1_restriction,MC_successes(i),min_lik_coefficient_LetJoeDecide,SaveFigs,q,ENRICH); % calculate for FA data with the control in it
        else %dylan here we calculate chol, where model is slightly different
            choldat=[datas{i},zeros(rows(datas{i}),35-cols(datas{i}))]; %dylan: shit I don't have to do that much work for cholesterol! this will create a big zero matrix for that!
            [Mod1_R, Mod1_C]=Cholesterol_Automated_Batch(choldat,cutoff(i),controls{i},0,q,lambda,user_p, ENRICH); %dylan : added q to modeling function 0.99 is for testing the difference in enrichment don't d 0.99 for the screen
           %[Model_results, Model_colnames]=Cholesterol_Automated_Batch(data,cutoff,num_controls,SaveFigs,q,lambda,user_p,ENRICH)
            
            Mod3_R=[];
            Mod3_C=[];
        end
        
        %dylan: return to workspace
        %%assignin(ws, 'Mod1_R', Mod1_R);
        %%assignin(ws, 'Mod1_C', Mod1_C);
        %%assignin(ws, 'Mod3_R', Mod3_R);
        %%assignin(ws, 'Mod3_C', Mod3_C);       
        %%assignin(ws,'MC_Stats',MC_Stats);
        
        
        %%dylan
        %display('freezing the result data, hands off the computer!')
        
            
        code=codes{i};
        count=0;
        count_ctrl=0;
        clear dat_2_write1 dat_2_write3;
        temp_plates=plates{i};
        for j=1:length(code); 
            if(code(j)==1)
                count=count+1;
                if ( elongation ==0 || elongation ==1)
                  dat_2_write1(j,:)=Mod1_R(count,:);
                end 
                plates_2_use{count}=temp_plates{j+1};
                if(~isempty(Mod3_R))
                    dat_2_write3(j,:)=Mod3_R(count,:);
                end
            elseif(code(j)==0)
                count_ctrl=count_ctrl+1;
                if ( elongation ==0 || elongation ==1)
                    dat_2_write1(j,:)=Mod1_R(rows(Mod1_R)-controls{i}+count_ctrl,:);
                end
                if(~isempty(Mod3_R))
                    dat_2_write3(j,:)=Mod3_R(rows(Mod3_R)-controls{i}+count_ctrl,:);
                end
            else
                if ( elongation ==0 || elongation ==1)
                  dat_2_write1(j,:)=-1*ones(1,cols(Mod1_R));
                end
                if(~isempty(Mod3_R))
                    dat_2_write3(j,:)=-1*ones(1,cols(Mod3_R));
                end
            end
        end
        display('Writing Results...')
        Sheet=lipids{i};
        warning('off','all');
        if ( elongation ==0 || elongation ==1)
           xlswrite(SaveFile1,dat_2_write1,Sheet,['B',num2str(StartRow+1)]);
           xlswrite(SaveFile1,Mod1_C,Sheet,['B',num2str(StartRow)]);     %Write column header
        end
        
        if(~isempty(Mod3_R))
            xlswrite(SaveFile3,dat_2_write3,Sheet,['B',num2str(StartRow+1)]);
            xlswrite(SaveFile3,Mod3_C,Sheet,['B',num2str(StartRow)]);     %Write column header
           	%dylan
            if (SaveMCStats)
                   write_MCStats_tofile(['MCStats-',lipids{i},'-',datafile],MC_Stats,plates_2_use);
            end
        end
        warning('on','all');
    end
    
    display('Analysis finished!')
    
end     
    
function write_MCStats_tofile(SaveFile, MC_Stats, plates_2_use)
    colhead={'-SSE','D0','D1','D2','S','I','IE1','IE2','IE3','IE4','IE5'}; %181118 JA update colhead values to match manuscript
    rowhead={'max';'mean';'var';'Results';'...'}; 

    for i=1:length(plates_2_use)
        Sheet=plates_2_use{i};
        xlswrite(SaveFile,MC_Stats{i},Sheet,'B2');
        xlswrite(SaveFile,colhead,Sheet,'B1');     %Write column header
        xlswrite(SaveFile,rowhead,Sheet,'A2');
    end
end


