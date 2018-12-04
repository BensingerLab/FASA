%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function [q,lambda,log_lik,mix]=fit_cholest_back_auto(data,cutoff)


    steps=[.01,.001,.0001];
    [q,lambda,log_lik]=fit_Cholest_background(data, .001:.01:.1, 0.001:.01:.2, cutoff);
    h=gcf;
    i=2;
    while i<=length(steps)

        q_range=max(0,q-(2*steps(i-1))):steps(i):min(1,q+(2*steps(i-1)));
        lambda_range=max(0,lambda-(2*steps(i-1))):steps(i):min(1,lambda+(2*steps(i-1)));
        set(0, 'CurrentFigure', h)
        [q,lambda,log_lik]=fit_Cholest_background(data,q_range,lambda_range,cutoff);

        if(q==max(q_range)||q==min(q_range)||lambda==max(lambda_range)||lambda==min(lambda_range)) %make sure we end at a local maxima
            if(q~=0 && q~=1 && lambda~=0 && lambda~=1)
                i=i-1;
            end
        end
        i=i+1;
    end
    mix=Cholest_distr(0,q,0);
    mix=loss_distribution(mix,lambda); 