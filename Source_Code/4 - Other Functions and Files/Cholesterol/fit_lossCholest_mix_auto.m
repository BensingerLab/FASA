%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
% 2015.08.25 added parameter user_p to restrict the range of p
function [p,s,log_lik,mix_distr]=fit_lossCholest_mix_auto(data,q,lambda,cutoff,user_p,enrich) % Dylan added a new parameter: user_p for te range. if you label 50% glucose, it should be <0.5~ haven't changed bigger numbers
    
    steps=[.1,.01,.001,.0001];
    [p,s,log_lik]=fit_lossCholest_mix(data,q,lambda, user_p(1):.1:user_p(2) , 0:.1:1 ,cutoff,enrich); % dylan temp
    
    h=gcf;
    i=2;
    while i<=length(steps)

        p_range=max(user_p(1),p-(2*steps(i-1)))  :  steps(i)  :   min(user_p(2),p+(2*steps(i-1)));
        if(max(p_range) < user_p(2)  && max(p_range)+steps(i)>user_p(2)) % if range is close to the limit, than extend the range
            p_range=[p_range,user_p(2)];
        end
        
        s_range=max(0,s-(2*steps(i-1))):steps(i):min(1,s+(2*steps(i-1)));
        if(max(s_range) < 1  && max(s_range)+steps(i)>1)
            s_range=[s_range,1];
        end
        
        set(0, 'CurrentFigure', h)
        [p,s,log_lik]=fit_lossCholest_mix(data,q,lambda,p_range,s_range,cutoff,enrich);
        title(['p=',num2str(p),' s=',num2str(s)])
        if((p==max(p_range) || p==min(p_range)) && p~=user_p(2) && p~=user_p(1))
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
                    p_range=p-steps(i-1):steps(i):min(user_p(2),p+(50*steps(i-1))); %dylan
                    %p_range=p-steps(i-1):steps(i):min(1,p+(50*steps(i-1)));%moses
                    %dylan: adding the following one line
                    %user_p(2)=min(1,p+(4*steps(i-1))) % for tmp use
                elseif(p==min(p_range))
                    p_range=max(0,p-50*steps(i-1)):steps(i):p+steps(i-1);
                end
            else
                p_range=[max(user_p(1),p-steps(i)*3):steps(i):min(user_p(2),p+3*steps(i))];
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

            [p,s,log_lik]=fit_lossCholest_mix(data,q,lambda,p_range,s_range,cutoff,enrich);
            i=i-1;
        end
       
        i=i+1;
    end
    if ( p == user_p(2) )
        p=nan
    end

    mix_distr=Cholest_distr(p,q,s,enrich);
    mix_distr=loss_distribution(mix_distr,lambda);