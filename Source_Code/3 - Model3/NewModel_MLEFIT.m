%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).
function [b,its,liks,varargout]=NewModel_MLEFIT(y,weights,fix_p,s_e_fix_index,b0,D_1_restriction,thresh,max_it,q,ENRICH,v)
                                %NewModel_MLEFIT(data,weights,fix_p,s_e_fix_index,b0,D_1_restriction,tol,maxit,q,v);
                                %v=0.015 to adjust gradient
%% this function can be called by Model3_MonteCarlo to run interations for FASA
% nomenclature: 
%   Model 3 means FASA, Model 1 mean classic ISA
%   p represents a vector with [D0, D1, D2]
% Outputs:
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
% 2016.2.26 Dylan fixing [p0,p1,p2] when p_source is not 0 s
% 20180326: Dylan adding handle for only iterating D and 1-D rather
%   than p_2, p_1 and p_0
% 20181021: Dylan updated screen display, and add this header
%%
%b0={p,elong}
    p=b0{1};
    elong=b0{2};
    b={p,elong}; % b is a cell array. every row is like b0
    its=0; % interation#
    step=.001; 
    n=2+length(elong)-1;
    
    y=y/sum(y); % the real data
    inarow=0;
    badconverge=0;
    
    m=length(y);
    W=eye(m);
    for i=1:length(weights)
        W(i,i)=weights(i)^2;
    end
    converge=0;


    yhat=elongation_distribution(q,p,elong); % returns the distribution
    
    log_lik=-sum(weights.*((y-yhat).^2));
    liks(1)=log_lik;


    while(~converge && its<max_it)
        its=its+1;
        p=b{its,1};
        grad=zeros(1,n);  % grad is the vector with dLiklihood/db for each parameter

        %%%%%%%%%%%%%%%%%%% Set new p for next iteration %%%%%%%%%%%%%%%%
        if (fix_p)
            p=b0{1};
        else
            %%%%%%%%%%Calc d/dp
            if D_1_restriction==0 % generates new p_2, p_1, p_0 by varying p_0 and p_1 
                for i=1:2
                    stp=min([step,1-p(i),p(3)-0]);
                    p(i)=p(i)+stp;
                    p(3)=p(3)-stp;
                    y_test=elongation_distribution(q,p,elong);

                    lik=-sum(weights.*((y-y_test).^2));
                    grad(i)=(lik-liks(its))/stp;
                    p=b{its,1}; %back to its original value
                end        
                %%%% set new p according to Calc d/dp
                grad_p=grad(1:2);
                grad_p=grad_p*v;
                p(1)=fittowindow(p(1)+grad_p(1),[0,1]);
                p(2)=fittowindow(p(2)+grad_p(2),[0,1]);
                p(3)=fittowindow(1-p(1)-p(2),[.00001,1]);
                p=p/sum(p); 
            else
                % Dylan: now D_1_restriction=1, so varies p(1) or D_0 and restricts p(2), p(3) or D_1, D_2 accordingly
                %p is a 3x1 vector holding [p0,p1,p2], the probabilities of a zero, single,
                    %     a double labeled ACoA in the cellular pool. sum(pi)=1                
                %%%%%%%%% Dylan: formula from Joe:
                %% D = % of AcCoA pool carbon contributed by labeled metabolite(s)
                %% D_x = % of AcCoA pool that contains x 13Cs
                %% e = % 13C in labeled metabolite(s)
                %% q = % 13C in natural metabolites
                % then D_0 = (1-D)*(1-q)^2 + D * (1-e)^2; 
                % and D_1 = 2*(1-D) * q*(1-q) + 2*D * e* (1-e);
                % and D_2 = D * e^2 + (1-D) * q^2
                % here p(1)=D_0; p(2)=D_1 and p(3)=D_2
                % solving the equation we have D and D_2 below:
                
                stp=min([step,1-p(1),p(3)-0]);
                p(1)=p(1)+stp;
                
                e=ENRICH;
                D_0=p(1);                
                D= (D_0 - (1-q)^2 ) / ( (1-e)^2 -(1-q)^2 );
                D_1 = 2*(1-D) * q*(1-q) + 2*D * e* (1-e);
                
                p(2)=D_1;
                p(3)=1-p(1)-p(2);
                y_test=elongation_distribution(q,p,elong);
                lik=-sum(weights.*((y-y_test).^2));
                
                %grad(i)=(lik-liks(its))/stp;
                grad=(lik-liks(its))/stp;          
                
                p=b{its,1}; %back to its original value
                
                %%%% set new p according to Calc d/dp
                grad_p=grad;
                grad_p=grad_p*v;
                p(1)=fittowindow(p(1)+grad_p,[0,1]);
                
                D_0=p(1);                
                D= (D_0 - (1-q)^2 ) / ( (1-e)^2 -(1-q)^2 );
                D_1 = 2*(1-D) * q*(1-q) + 2*D * e* (1-e);
                
                p(2)=D_1;
                p(3)=fittowindow(1-p(1)-p(2),[.00001,1]);
                p=p/sum(p); 
            end
            

        end
      %%%%%%%%%%%%%%%%%%% End of Setting new p for next iteration %%%%%%%%%%%%%%%%
      
      %%%%%%%%%%%%%%%%%%% Set new s and e series for next iteration %%%%%%%%%%%%%%%%
      s_e_unfix_index=setdiff(1:length(elong),s_e_fix_index);
      elong_unfix=elong(s_e_unfix_index);  
      %%%%%%% Calc d/delong
      for i= 1:(length(elong_unfix)-1)
        stp=min([step,1-elong_unfix(i),elong_unfix(length(elong_unfix))-0]);
        elong_unfix(i)=elong_unfix(i)+stp;
        elong_unfix(length(elong_unfix))=elong_unfix(length(elong_unfix))-stp;
        %dylan
        elong(s_e_unfix_index)=elong_unfix;
        %
        y_test=elongation_distribution(q,p,elong);

        lik=-sum(weights.*((y-y_test).^2));
        grad(i+2)=(lik-liks(its))/stp;
        elong=b{its,2}; %elong back
        elong_unfix=elong(s_e_unfix_index);   %elong unfix back
      end
      
       %%%% set new s and e series according to Calc d/elong
       grad=grad*v;
       for i= 1:(length(elong_unfix)-1)
           elong_unfix(i)=fittowindow(elong_unfix(i)+grad(2+i),[0,1]);
       end
       elong_unfix(length(elong_unfix))=0;
       elong(s_e_unfix_index)=elong_unfix;
       elong_unfix(length(elong_unfix))=fittowindow(1-sum(elong),[0.0001,1]);
       elong(s_e_unfix_index)=elong_unfix;
       elong=elong/sum(elong);
       %%%%%%%%%%%%%%%%%%% End of Setting new s and slong for next iteration %%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%% End of Setting new p, s and e series for next iteration %%%%%%%%%%%%%%%%
        
        b(its+1,:)={p,elong};
        yhat=elongation_distribution(q,p,elong);
        
        log_lik=-sum(weights.*((y-yhat).^2));
        liks(its+1)=log_lik;

        if(log_lik<=liks(its))
            inarow=0;
            v=v/1.5;
        else
            inarow=inarow+1;
        end

        if(inarow==20)
            v=min(v*1.5,100);
            inarow=0;
        end

        if((its<=100 && mod(its,20)==0) || mod(its,100)==0)    
            if (numel(elong)>=4)
                 display(['Iteration:',num2str(its),' -SSE=',num2str(log_lik), ' D2=',num2str(p(3)),' S=',num2str(elong(1)), ' I=',num2str(elong(2)), ' IE1=',num2str(elong(3)), ' IE2=',num2str(elong(4))]);
            % 181021 Dylan updated to match nomenclature
            elseif  (numel(elong)>=3)
                display(['Iteration:',num2str(its),' -SSE=',num2str(log_lik), ' D2=',num2str(p(3)),' S=',num2str(elong(1)), ' I=',num2str(elong(2)), ' IE1=',num2str(elong(3))]);
            % 181021 Dylan updated to match nomenclature
            end 
        end

        if(sqrt(sum(grad.^2)/(2+length(elong)))<thresh)
            converge=1;
            display(['Convergence after', num2str(its), ' iterations.'])
        end
        if(isnan(log_lik))
            converge=1;
            badconverge=1;
        end
        bar([y',yhat'])
        %%%%%%%%%%%%%%%%%%%%
        pause(.01);
    end
    yhat=elongation_distribution(q,p,elong);
    
    varargout{1}=yhat;
    if its>=max_it
        display(['Tired Convergence after ',num2str(its),' iterations'])
    end
    if(badconverge)
        display('Broken Convergence')
    end
end
