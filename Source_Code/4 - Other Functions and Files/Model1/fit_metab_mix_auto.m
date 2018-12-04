%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function [p,s,log_lik,mix]=fit_metab_mix_auto(data,q,cutoff,enrich)
%% this function can be called by BothModels.m to run Model 1 (classicISA)
% nomenclature: 
%   Model 3 means FASA, Model 1 mean classic ISA
%   p means contribution of tracer to Acetyle-coA pool.same as "D" in paper
% dependency: fit_metab_mix.m
%%
    steps=[.1,.01,.001,.0001];
    [p,s,log_lik]=fit_metab_mix(data,q,[0:.1:1],[0:.1:1],cutoff,enrich);
    h=gcf;
    i=2;
    while i<=length(steps)

        p_range=[max(0,p-(steps(i-1))):steps(i):min(1,p+(steps(i-1)))];
        if(max(p_range) < 1  && max(p_range)+steps(i)>1)
            p_range=[p_range,1];
        end
        s_range=[max(0,s-(steps(i-1))):steps(i):min(1,s+(steps(i-1)))];
        if(max(s_range) < 1  && max(s_range)+steps(i)>1)
            s_range=[s_range,1];
        end
        set(0, 'CurrentFigure', h)
        [p,s,log_lik]=fit_metab_mix(data,q,p_range,s_range,cutoff,enrich);
        
        if((p==max(p_range) || p==min(p_range)) && p~=1 && p~=0)
            p_badfit=1;
        else
            p_badfit=0;
        end
        if((s==max(s_range) || s==min(s_range)) && s~=1 && s~=0)
            s_badfit=1;
        else
            s_badfit=0;
        end
        
        
        if(p_badfit || s_badfit) %make sure we end at a local maxima 
            if(p_badfit)
                if(p==max(p_range))
                    p_range=p-steps(i-1):steps(i):min(1,p+(50*steps(i-1)));
                elseif(p==min(p_range))
                    p_range=max(0,p-50*steps(i-1)):steps(i):p+steps(i-1);
                end
            else
                p_range=[max(0,p-steps(i)*3):steps(i):min(1,p+3*steps(i))];
            end


            
            if(s_badfit)
                if(s==max(s_range))
                    s_range=s-steps(i-1):steps(i):min(1,s+(50*steps(i-1)));
                elseif(s==min(s_range))
                    s_range=max(0,s-50*steps(i-1)):steps(i):s+steps(i-1);
                end
            else
                s_range=[max(0,s-steps(i)*3):steps(i):min(1,s+3*steps(i))];
            end
            

            [p,s,log_lik]=fit_metab_mix(data,q,p_range,s_range,cutoff,enrich);
            i=i-1;
        end
        i=i+1;
    end
    
    mix=Label_distr_mix(p,q,s,(length(data)/2-1),enrich);