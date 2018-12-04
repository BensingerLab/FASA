%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function [str_Mod1_R, str_Mod1_C, str_Mod3_R, str_Mod3_C, str_MC_Stats]=FASA_Step2(sheets_to_model,cutoff,q,e,assume_no_label_diffusion,MC_successes) % 181027 JA updated D_1_restriction naming.
% function [str_Mod1_R, str_Mod1_C, str_Mod3_R, str_Mod3_C, str_MC_Stats]=FASA_Step2(sheets_to_model,cutoff,q,e,D_1_restriction,MC_successes) % 181027 JA Old
%% this function runs FASA model to generate s, e0, e1, e2... and D0, D1, D2
% dependency: Metabolism_Preprocess_UCLA.m
% Inputs: 
%   sheets_to_model,cutoff,q,e,D_1_restriction,MC_successes

%   sheets_to_model: the sheetsyou want to run. e.g: =1:13 for all 13 sheets, [1,3,4,7] for selected sheets
%   cutoff: threshold for AUC data used in FASA. 
%           Any AUC value < cutoff will not contribute to the cost function or the fitting.
%           cutoff=[50,100,200], means cutoff=50 for the 1st FA, 100 for the 2nd, 200 for the rest              
%   q: the natural occuring 13C percentage, =0 the script will calclute q from control, if set =0.112, then this is the value that will be used in the fitting
%   e: enrichment of 13C in the isotopic tracer. 
%           e=0.99 for 13C-U-Glucose means 99% of the carbons are 13C 
%   assume_no_label_diffusion: if =0, FASA iterates D0, D1, and D2=1-D0-D1;
%                    if =1, FASA iterates D0, with D1, D2 calculated from D0
%   MC_successes: number of trials for modeling. 
%                     the modeled result with SSE closest to 0 will be
%                     written in the final result.xls file
%   (hidden) AUC data: users need to choose .xls file(s) with AUC data
%
% Outputs example:
% =====files generated=====
%   Result-Model3..xls: 
%       every sheet is for a FAtty acid, s, ie, D series and new distribution generated
%   MCStats-x160... xlsx:
%       there will be a .xlsx per FAtty acid, where each sheet is for one
%       sample, stats summarizes caculated parameters for each MonteCarlo
%       iteration
%   Figures: comparison of origianl distribution and modeled distribution
%% 
[datafiles,dir]=uigetfile('*.xls','select AUC data file','MultiSelect','on');
if (class(datafiles)=='char')
    temp={[]};
    temp{1}=datafiles;
    datafiles=temp;
end
%dir=uigetdir('','Select Folder for Save Data');
cd(dir);

% "Internal" User Parameters
sheets_2_compute=sheets_to_model;
ENRICH=e;
D_1_restriction=assume_no_label_diffusion; % 181027 JA updated D_1_restriction naming.
lambda=0.023; % percentage of fragmentation for cholesterol data. nonsenese for fatty acid. if lambda=0, scripts will either calculate from unlabel controls or use default=0.023
user_p=[0,1];  % this is the user handle for range of percentage contribution of glucose
elongation=2; % elongation=0 only runs model 1, 
               % elongation=1 runs both model 1 & 3, 
                     % use the likelihood from Model 1 as a termination criteira for model 3 iterations.
                     %and plot Model 1 and 3 results together
               % elongation=2 runs only Model 3, if there's no user setup p_s_e, it will set initial p =[.425,.025,.55],
                    %and randomize e and s
                    % since no likihood from Model1 is available, I set it to -1.
min_lik_coefficient_LetJoeDecide=1000000;  % 2016.6.14when running Model 3, only when the new likelihood > 1.02*likelihood from Model 1 is called a success. otherwise it's marked "fail", and only 5 times of fails are allowed.
SaveFigs=1; % if =1 will save all the fancy images in the folder
SaveMCStats=1; % ignore this


for i=1:length(datafiles)
   
    datafile=datafiles{i}
    
    starting_time=clock();
    display(['starting time: ',num2str(starting_time)])
    Metabolism_Preprocess_UCLA(datafile,cutoff,q,MC_successes,sheets_2_compute,elongation,D_1_restriction,SaveFigs,SaveMCStats,lambda,user_p,min_lik_coefficient_LetJoeDecide,ENRICH);

    finishing_time=clock();
    display(['finishing_time: ',num2str(finishing_time)])
end
