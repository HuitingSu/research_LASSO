function [Model_NoCross,Psc_opt_currentmodel,Lambda_opt_currentmodel]=LassoSVTSFast_NoCross(X,Y,SigmaTrue)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PSC fast algorithm
%Input: X, Y,
% Model_SVTSFast_Cross=Model_SVTSFast_CrossFun;
% Model_SVTSFast_Interval=Model_SVTSFast_IntervalFun;
% Model_SVTSFast_Derivative=Model_SVTSFast_DerivativeFun;
% P_Model_SVTSFast_Intervalopt_All=P_Model_SVTSFast_Intervalopt_AllFun;
% Psc_Lambda_3=Psc_Lambda_3Max;
%Generate intervals from Lasso: Obtain the Lasso Path intervals; find Xa and Lambda_Current(Lambda_1=LambdaMin,Lambda_2=LambdaMax) for each interval
%For each interval:
%Obtain Beta_OLS_Current
%Obtain P2(Lambda_1), P1(Lambda_1),P2(Lambda_2), P1(Lambda_2)
%Obtain (Lambda_1, P2(Lambda_1)), (Lambda_2, P2(Lambda_2)) = (a2,-b2)
%Obtain (Lambda_1, P1(Lambda_1)), (Lambda_2, P1(Lambda_2)) = (a1,b1)

%Test 1: Lambda_0 \in (Lambda_1, Lambda_2), (0, Lambda_1), (Lambda_02, \infinity)

%Test 2: Lambda_opt=1/2(a2/b2-a1/b1), P_opt= (a1^2b2^2+b1^2a2^2-a1a2b1b2)/ b1b2


%Test 3: if: -b2/P2(Lambda_1)+b1/P2(Lambda_1)> 0, \exit Lambda_3 \in (Lambda_1,lambba_2) s.t. Psc(Lambda_3)>Psc(Lambda_1)
%        if:  b2/P2(Lambda_2)-b1/P2(Lambda_2)> 0, \exit Lambda_3 \in (Lambda_1,lambba_2) s.t. Psc(Lambda_3)>Psc(Lambda_2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
iik=1;
Design=X;
Pace=0.0001;
[n,p]=size(X);
[A_setA_Path,Beta_Path,Lambda_Path]=OriLasso(X,Y);
%Lambda_Path(size(Lambda_Path,2)+1)=0.0001;
%         A_setA_Path(1:6,1:6)
%Maxmodel=min(size(A_setA_Path,2)-6,n/2,8);
Maxmodel=min(min(floor(n/2)+2,floor(p/3)+2),size(Lambda_Path,2)-1);
Maxmodel=max(min(min(floor(n/2)+2,floor(p/3)+2),size(Lambda_Path,2)-1),size(A_setA_Path,2))-1;
Selectmodel=zeros(p,1);
Model_SVTSFast_CrossFun=zeros(p,1);
Model_SVTSFast_IntervalFun=zeros(p,1);
Model_SVTSFast_DerivativeFun=zeros(p,1);
Model_NoCross=zeros(p,1);

P_Model_SVTSFast_Intervalopt_AllFun=zeros(p,1);
Psc_Lambda_3Max=zeros(p,1);
Lambda_opt_k11i=zeros(p,1);
Psc_Lambdaopt_Max=zeros(p,1);
DistanceofOptLambda=100*ones(p,1);

Lambda_SVOE_OLS=0;
PCSCS_SVOE_OLS=0;
Lambda_3_AllFun=zeros(p,1);
Lambda_Model_SVTSFast_Intervalopt_AllFun=zeros(p,1);
P1_Lambda_All=zeros(p,1);
P2_Lambda_All=zeros(p,1);
Psc_Lambda_3Max=zeros(p,1);
Psc_Lambda_1Fun=zeros(p,1);
Psc_Lambda_2Fun=zeros(p,1);
P1_Lambda_1Fun=zeros(Maxmodel,Maxmodel+1);
P2_Lambda_1Fun=zeros(Maxmodel,Maxmodel+1);
Psc_Lambda_1Fun=zeros(Maxmodel,Maxmodel+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate the sigma, error  Cheng 2005, the largest model
% sigma^2=1/N Sum E_i^2
% Sigma=(1/n*(Y-X*Beta_Path(:,size(Beta_Path,2)))'*(Y-X*Beta_Path(:,size(Beta_Path,2))))^(1/2);
% for kki2=1:size(Beta_Path,2)
%     Beta_Path_Curent2=Beta_Path(:,kki2);
%     SigmaBeta2(kki2)=(1/n*(Y-X*Beta_Path_Curent2)'*(Y-X*Beta_Path_Curent2))^(1/2);
% end
% bb1=find(A_setA_Path(:,size(A_setA_Path,2))==0);
% X1=X(:,A_setA_Path(1:(bb1(1)-1),size(A_setA_Path,2)));
% Sigma1=(1/n*(Y-X1*inv(X1'*X1)*X1'*Y)'*(Y-X1*(X1'*X1)*X1'*Y))^(1/2);
% bb2=find(A_setA_Path(:,3)==0);
% X2=X(:,A_setA_Path(1:(bb2(1)-1),3));
% Sigma2=(1/n*(Y-X2*inv(X2'*X2)*X2'*Y)'*(Y-X2*(X2'*X2)*X2'*Y))^(1/2);

% for kki1=1:size(A_setA_Path,2)
%     bbkki1=find(A_setA_Path(:,kki1)==0);
%     Xkki1=X(:,A_setA_Path(1:(bbkki1(1)-1),kki1));
%     
%     Sigmakki1(kki1)=(1/n*(Y-Xkki1*inv(Xkki1'*Xkki1)*Xkki1'*Y)'*(Y-Xkki1*(Xkki1'*Xkki1)*Xkki1'*Y))^(1/2);
% end
%         Sigmakki1
%         SigmaBeta2
Sigma=SigmaTrue;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lambda_Path_Sort=sort(Lambda_Path);
Lambda_Current=Lambda_Path_Sort(size(Lambda_Path,2)-1);

xxplot_k11i=zeros(Maxmodel,71);
P1_Lambda_plot_k11i=zeros(Maxmodel,71);
P2_Lambda_plot_k11i=zeros(Maxmodel,71);%P2_Lambda_plot;
Psc_Lambda_plot_k11i=zeros(Maxmodel,71);%Psc_Lambda_plot;
xLambda1_k11i=zeros(Maxmodel,11);%linspace(Lambda_1_tem,Lambda_1_tem,11);
yLambda1_k11i=zeros(Maxmodel,11);%0:0.1:1;
xLambda2_k11i=zeros(Maxmodel,11);%linspace(Lambda_2_tem,Lambda_2_tem,11);
yLambda2_k11i=zeros(Maxmodel,11);%0:0.1:1;
xLambda3_k11i=zeros(Maxmodel,11);%linspace(Lambda_3,Lambda_3,11);
yLambda3_k11i=zeros(Maxmodel,11);%0:0.1:1;



for k11_i=2:Maxmodel
    
    A_setA_PathAll0=A_setA_Path(:,k11_i);
    if sum(A_setA_PathAll0==0)==0
        Index00=size(A_setA_PathAll0,1);
    else
        [a11,bbb]=find(A_setA_PathAll0==0);
        Index00=a11(1)-1;
    end


    A_setA_onecolumn0=sort(A_setA_PathAll0(1:Index00(1)));
    Position0=A_setA_onecolumn0;
    Position_Current=Position0;
    Xa=X(:,Position_Current);
    BetaOLS=inv(Xa'*Xa)*(Xa'*Y);
    BetaZeroOLS=zeros(p,1);
    BetaZeroOLS(A_setA_onecolumn0)=BetaOLS;
    Beta_Current_OLS=BetaZeroOLS;
    %                  Beta_Current_OLS_All(:,k11_i)=BetaZeroOLS;
    %             LambdaPath_All(iik,k11_i)=Lambda_Path(k11_i);
    %
    %             Lambda_2=Lambda_Path(iik,k11_i);
    %             Lambda_1=Lambda_Path(iik,k11_i+1);
    for k12_i=2:Maxmodel
        if p<15  && (size(find(Beta_Current_OLS~=0),1)<6)
            [P1_Lambda_1,P2_Lambda_1]=P1P2Lambda(Design,Beta_Current_OLS,Lambda_Path(k12_i),Sigma);
            P1_Lambda_1Fun(k11_i,k12_i)=P1_Lambda_1;
            P2_Lambda_1Fun(k11_i,k12_i)=P2_Lambda_1;
            Psc_Lambda_1Fun(k11_i,k12_i)=P1_Lambda_1*P2_Lambda_1;
        elseif p>15  && (size(find(Beta_Current_OLS~=0),1)<6)
            [P2_Lambda_1]=P2LambdaCDF(Design,Beta_Current_OLS,Lambda_Path(k12_i),Sigma);
            [P1_Lambda_1]=P1LambdaSimu(Design,Beta_Current_OLS,Lambda_Path(k12_i),Sigma);
            P1_Lambda_1Fun(k11_i,k12_i)=P1_Lambda_1;
            P2_Lambda_1Fun(k11_i,k12_i)=P2_Lambda_1;
            Psc_Lambda_1Fun(k11_i,k12_i)=P1_Lambda_1*P2_Lambda_1;            
        else
            [P2_Lambda_1]=P2LambdaSimu(Design,Beta_Current_OLS,Lambda_Path(k12_i),Sigma);
            [P1_Lambda_1]=P1LambdaSimu(Design,Beta_Current_OLS,Lambda_Path(k12_i),Sigma);
            P1_Lambda_1Fun(k11_i,k12_i)=P1_Lambda_1;
            P2_Lambda_1Fun(k11_i,k12_i)=P2_Lambda_1;
            Psc_Lambda_1Fun(k11_i,k12_i)=P1_Lambda_1*P2_Lambda_1;
        end
    end
    %%%%%%%%%%%Get the largest three of the adjecent Psc_Lambda_1Fun?
    %%%%%%%%%%%Not have to be adjecent.
    [alam,blam]=max(Psc_Lambda_1Fun(k11_i,1:Maxmodel));
    Lambda_max=Lambda_Path(blam);
    P1_Lambda_max=P1_Lambda_1Fun(k11_i,blam);
    P2_Lambda_max=P2_Lambda_1Fun(k11_i,blam);
    Psc_Lambda_max=P1_Lambda_max*P2_Lambda_max;

    Lambda_left=Lambda_Path(blam+1);
    P1_Lambda_left=P1_Lambda_1Fun(k11_i,blam+1);
    P2_Lambda_left=P2_Lambda_1Fun(k11_i,blam+1);
    Psc_Lambda_left=P1_Lambda_left*P2_Lambda_left;

    if blam>1
        Lambda_right=Lambda_Path(blam-1);
        P1_Lambda_right=P1_Lambda_1Fun(k11_i,blam-1);
        P2_Lambda_right=P2_Lambda_1Fun(k11_i,blam-1);
        Psc_Lambda_right=P1_Lambda_right*P2_Lambda_right;
    else
        Lambda_right=Lambda_Path(blam)+1/5*max((Lambda_Path(blam)-Lambda_Path(blam+1)),(Lambda_Path(blam)));
        if p<15  && (size(find(Beta_Current_OLS~=0),1)<6)
            [P1_Lambda_right,P2_Lambda_right]=P1P2Lambda(Design,Beta_Current_OLS,Lambda_right,Sigma);
            Psc_Lambda_right=P1_Lambda_right*P2_Lambda_right;
        elseif p>15  && (size(find(Beta_Current_OLS~=0),1)<6)
            [P2_Lambda_right]=P2LambdaCDF(Design,Beta_Current_OLS,Lambda_right,Sigma);
            [P1_Lambda_right]=P1LambdaSimu(Design,Beta_Current_OLS,Lambda_right,Sigma);
            Psc_Lambda_right=P1_Lambda_right*P2_Lambda_right;
        else
            [P2_Lambda_right]=P2LambdaSimu(Design,Beta_Current_OLS,Lambda_right,Sigma);
            [P1_Lambda_right]=P1LambdaSimu(Design,Beta_Current_OLS,Lambda_right,Sigma);
            Psc_Lambda_right=P1_Lambda_right*P2_Lambda_right;
        end
    end



    Lambda_maxrightInterval=(Psc_Lambda_max*Lambda_max+Psc_Lambda_right*Lambda_right)/(Psc_Lambda_max+Psc_Lambda_right);
    if p<15  && (size(find(Beta_Current_OLS~=0),1)<6)
        [P1_Lambda_maxrightInterval,P2_Lambda_maxrightInterval]=P1P2Lambda(Design,Beta_Current_OLS,Lambda_maxrightInterval,Sigma);
        Psc_Lambda_maxrightInterval=P1_Lambda_maxrightInterval*P2_Lambda_maxrightInterval;
    elseif p>15  && (size(find(Beta_Current_OLS~=0),1)<6)
        [P2_Lambda_maxrightInterval]=P2LambdaCDF(Design,Beta_Current_OLS,Lambda_maxrightInterval,Sigma);
        [P1_Lambda_maxrightInterval]=P1LambdaSimu(Design,Beta_Current_OLS,Lambda_maxrightInterval,Sigma);
        Psc_Lambda_maxrightInterval=P1_Lambda_maxrightInterval*P2_Lambda_maxrightInterval;
    else
        [P2_Lambda_maxrightInterval]=P2LambdaSimu(Design,Beta_Current_OLS,Lambda_maxrightInterval,Sigma);
        [P1_Lambda_maxrightInterval]=P1LambdaSimu(Design,Beta_Current_OLS,Lambda_maxrightInterval,Sigma);
        Psc_Lambda_maxrightInterval=P1_Lambda_maxrightInterval*P2_Lambda_maxrightInterval;
    end


    Lambda_leftmaxInterval=(Psc_Lambda_left*Lambda_left+Psc_Lambda_max*Lambda_max)/(Psc_Lambda_left+Psc_Lambda_max);
    if p<15  && (size(find(Beta_Current_OLS~=0),1)<6)
        [P1_Lambda_leftmaxInterval,P2_Lambda_leftmaxInterval]=P1P2Lambda(Design,Beta_Current_OLS,Lambda_leftmaxInterval,Sigma);
        Psc_Lambda_leftmaxInterval=P1_Lambda_leftmaxInterval*P2_Lambda_leftmaxInterval;
    elseif p>15  && (size(find(Beta_Current_OLS~=0),1)<6)
        [P2_Lambda_leftmaxInterval]=P2LambdaCDF(Design,Beta_Current_OLS,Lambda_leftmaxInterval,Sigma);
        [P1_Lambda_leftmaxInterval]=P1LambdaSimu(Design,Beta_Current_OLS,Lambda_leftmaxInterval,Sigma);
        Psc_Lambda_leftmaxInterval=P1_Lambda_leftmaxInterval*P2_Lambda_leftmaxInterval;
    else
        [P2_Lambda_leftmaxInterval]=P2LambdaSimu(Design,Beta_Current_OLS,Lambda_leftmaxInterval,Sigma);
        [P1_Lambda_leftmaxInterval]=P1LambdaSimu(Design,Beta_Current_OLS,Lambda_leftmaxInterval,Sigma);
        Psc_Lambda_leftmaxInterval=P1_Lambda_leftmaxInterval*P2_Lambda_leftmaxInterval;   
    end


    Psc_five=[Psc_Lambda_left,Psc_Lambda_leftmaxInterval,Psc_Lambda_max,Psc_Lambda_maxrightInterval,Psc_Lambda_right];
    Lambda_five=[Lambda_left,Lambda_leftmaxInterval,Lambda_max,Lambda_maxrightInterval,Lambda_right];
    [Psc_MaxNo1,MaxNo1]=max(Psc_five);
    Lambda_MaxNo1=Lambda_five(MaxNo1);
    Psc_MaxNo1;

    Psc_four=Psc_five;
    Lambda_four=Lambda_five;
    Psc_four(MaxNo1)=[];
    Lambda_four(MaxNo1)=[];
    [Psc_MaxNo2,MaxNo2]=max(Psc_four);
    Lambda_MaxNo2=Lambda_four(MaxNo2);
    Psc_MaxNo2;

    Pscfinaltwo=[Psc_MaxNo1,Psc_MaxNo2];
    Lambdafinaltwo=[Lambda_MaxNo1,Lambda_MaxNo2];

    [af2,bf2]=min(Pscfinaltwo);
    Psc_Max01=Pscfinaltwo(bf2);
    Lambda_Max01=Lambdafinaltwo(bf2);
    [af22,bf22]=max(Pscfinaltwo);
    Psc_Max02=Pscfinaltwo(bf22);
    Lambda_Max02=Lambdafinaltwo(bf22);

    Lambda_Max03=(Psc_Max01*Lambda_Max01+Psc_Max02*Lambda_Max02)/(Psc_Max01+Psc_Max02);
    if p<15  && (size(find(Beta_Current_OLS~=0),1)<6)
        [P1_Max03,P2_Max03]=P1P2Lambda(Design,Beta_Current_OLS,Lambda_Max03,Sigma);
        Psc_Max03=P1_Max03*P2_Max03;
    elseif p>15  && (size(find(Beta_Current_OLS~=0),1)<6)
        [P2_Max03]=P2LambdaCDF(Design,Beta_Current_OLS,Lambda_Max03,Sigma);
        [P1_Max03]=P1LambdaSimu(Design,Beta_Current_OLS,Lambda_Max03,Sigma);
        Psc_Max03=P1_Max03*P2_Max03;
    else
        [P2_Max03]=P2LambdaSimu(Design,Beta_Current_OLS,Lambda_Max03,Sigma);
        [P1_Max03]=P1LambdaSimu(Design,Beta_Current_OLS,Lambda_Max03,Sigma);
        Psc_Max03=P1_Max03*P2_Max03;
    end

%     Psc_Max03
%     Lambda_Max03
    Psc_opt_currentmodel(k11_i)=Psc_Max03;
    Lambda_opt_currentmodel(k11_i)=Lambda_Max03;

% Psc_opt_currentmodel,Lambda_opt_currentmodel
    Lambda_Path(k11_i);
    
    if Lambda_opt_currentmodel(k11_i)<Lambda_Path(k11_i+1) %|| Lambda_opt_currentmodel(k11_i)>Lambda_Path(k11_i)
        DistanceofOptLambda(k11_i,1)=abs(Lambda_opt_currentmodel(k11_i)-Lambda_Path(k11_i+1))/abs(Lambda_Path(k11_i)-Lambda_Path(k11_i+1));
    elseif Lambda_opt_currentmodel(k11_i)>Lambda_Path(k11_i)
        DistanceofOptLambda(k11_i,1)=abs(Lambda_opt_currentmodel(k11_i)-Lambda_Path(k11_i))/abs(Lambda_Path(k11_i)-Lambda_Path(k11_i+1));
    else
        DistanceofOptLambda(k11_i,1)=0;
    end
    
end

DistanceofOptLambda(1,1)=100;
[am1,bm1]=min(DistanceofOptLambda(1:Maxmodel));

Model_NoCross(bm1)=1;


Model_NoCross(bm1)=1;



    %
    % xxplot_k11i(k11_i,:)=xxplot;
    % P1_Lambda_plot_k11i(k11_i,:)=P1_Lambda_plot;
    % P2_Lambda_plot_k11i(k11_i,:)=P2_Lambda_plot;
    % Psc_Lambda_plot_k11i(k11_i,:)=Psc_Lambda_plot;
    %
    % [Psc_optcurrent,bla1]=max(Psc_Lambda_plot);
    % Psc_Lambdaopt_Max(k11_i,1)=Psc_optcurrent;
    % Lambda_optcurrent=xxplot(bla1);
    % Lambda_opt_k11i(k11_i,1)=Lambda_optcurrent;
    % if Lambda_optcurrent<Lambda_1 %|| Lambda_optcurrent>Lambda_2
    %     DistanceofOptLambda(k11_i,1)=abs(Lambda_optcurrent-Lambda_1)/abs(Lambda_2-Lambda_1);
    % elseif Lambda_optcurrent>Lambda_2
    %     DistanceofOptLambda(k11_i,1)=abs(Lambda_optcurrent-Lambda_2)/abs(Lambda_2-Lambda_1);
    % else
    %     DistanceofOptLambda(k11_i,1)=0;
    % end
    % DistanceofOptLambda(k11_i,1)
    %
    %
    % xLambda1_k11i(k11_i,:)=linspace(Lambda_1,Lambda_1,11);
    % yLambda1_k11i(k11_i,:)=0:0.1:1;
    % xLambda2_k11i(k11_i,:)=linspace(Lambda_2,Lambda_2,11);
    % yLambda2_k11i(k11_i,:)=0:0.1:1;
    % xLambda3_k11i(k11_i,:)=linspace(Lambda_3,Lambda_3,11);
    % yLambda3_k11i(k11_i,:)=0:0.1:1;
    %
    %
    %
    % xLambda1=linspace(Lambda_1,Lambda_1,11);
    % yLambda1=0:0.1:1;
    % xLambda2=linspace(Lambda_2,Lambda_2,11);
    % yLambda2=0:0.1:1;
    % xLambda3=linspace(Lambda_3,Lambda_3,11);
    % yLambda3=0:0.1:1;
    % plot(xxplot,P1_Lambda_plot,xxplot,P2_Lambda_plot,xxplot,Psc_Lambda_plot,xLambda1,yLambda1,'-.b',xLambda2,yLambda2,'-.b',xLambda3,yLambda3,'-.r')
    %
    %
    %
    % plot(xxplot_k11i(1,19:43),...
    %     Psc_Lambda_plot_k11i(1,19:43),xLambda1_k11i(1,:),yLambda1_k11i(1,:),'-b',xLambda2_k11i(1,:),yLambda2_k11i(1,:),...
    %     '-b',...
    %     xxplot_k11i(2,19:43),...
    %     Psc_Lambda_plot_k11i(2,19:43),xLambda1_k11i(2,:),yLambda1_k11i(2,:),'-b',xLambda2_k11i(2,:),yLambda2_k11i(2,:),...
    %     '-b',...
    %     xxplot_k11i(3,19:43),...
    %     Psc_Lambda_plot_k11i(3,19:43),xLambda1_k11i(3,:),yLambda1_k11i(3,:),'-b',xLambda2_k11i(3,:),yLambda2_k11i(3,:),...
    %     '-b',...
    %     xxplot_k11i(4,19:43),...
    %     Psc_Lambda_plot_k11i(4,19:43),xLambda1_k11i(4,:),yLambda1_k11i(4,:),'-b',xLambda2_k11i(4,:),yLambda2_k11i(4,:),...
    %     '-b',...
    %     xxplot_k11i(5,19:43),...
    %     Psc_Lambda_plot_k11i(5,19:43),xLambda1_k11i(5,:),yLambda1_k11i(5,:),'-b',xLambda2_k11i(5,:),yLambda2_k11i(5,:),...
    %     '-b',...
    %     xxplot_k11i(6,19:43),...
    %     Psc_Lambda_plot_k11i(6,19:43),xLambda1_k11i(6,:),yLambda1_k11i(6,:),'-b',xLambda2_k11i(6,:),yLambda2_k11i(6,:),...
    %     '-b')
    %
    %
    %
    % plot(xxplot_k11i(1,:),...
    %     Psc_Lambda_plot_k11i(1,:),xLambda1_k11i(1,:),yLambda1_k11i(1,:),'-b',xLambda2_k11i(1,:),yLambda2_k11i(1,:),...
    %     '-b',...
    %     xxplot_k11i(2,:),...
    %     Psc_Lambda_plot_k11i(2,:),xLambda1_k11i(2,:),yLambda1_k11i(2,:),'-b',xLambda2_k11i(2,:),yLambda2_k11i(2,:),...
    %     '-b',...
    %     xxplot_k11i(3,:),...
    %     Psc_Lambda_plot_k11i(3,:),xLambda1_k11i(3,:),yLambda1_k11i(3,:),'-b',xLambda2_k11i(3,:),yLambda2_k11i(3,:),...
    %     '-b',...
    %     xxplot_k11i(4,:),...
    %     Psc_Lambda_plot_k11i(4,:),xLambda1_k11i(4,:),yLambda1_k11i(4,:),'-b',xLambda2_k11i(4,:),yLambda2_k11i(4,:),...
    %     '-b',...
    %     xxplot_k11i(5,:),...
    %     Psc_Lambda_plot_k11i(5,:),xLambda1_k11i(5,:),yLambda1_k11i(5,:),'-b',xLambda2_k11i(5,:),yLambda2_k11i(5,:),...
    %     '-b',...
    %     xxplot_k11i(6,:),...
    %     Psc_Lambda_plot_k11i(6,:),xLambda1_k11i(6,:),yLambda1_k11i(6,:),'-b',xLambda2_k11i(6,:),yLambda2_k11i(6,:),...
    %     '-b')
    %
    %
    %
    % plot(xxplot_k11i(1,:),P1_Lambda_plot_k11i(1,:),xxplot_k11i(1,:),P2_Lambda_plot_k11i(1,:),xxplot_k11i(1,:),...
    %     Psc_Lambda_plot_k11i(1,:),xLambda1_k11i(1,:),yLambda1_k11i(1,:),'-.b',xLambda2_k11i(1,:),yLambda2_k11i(1,:),...
    %     '-.b',...
    %     xxplot_k11i(2,:),P1_Lambda_plot_k11i(2,:),xxplot_k11i(2,:),P2_Lambda_plot_k11i(2,:),xxplot_k11i(2,:),...
    %     Psc_Lambda_plot_k11i(2,:),xLambda1_k11i(2,:),yLambda1_k11i(2,:),'-.b',xLambda2_k11i(2,:),yLambda2_k11i(2,:),...
    %     '-.b',...
    %     xxplot_k11i(3,:),P1_Lambda_plot_k11i(3,:),xxplot_k11i(3,:),P2_Lambda_plot_k11i(3,:),xxplot_k11i(3,:),...
    %     Psc_Lambda_plot_k11i(3,:),xLambda1_k11i(3,:),yLambda1_k11i(3,:),'-.b',xLambda2_k11i(3,:),yLambda2_k11i(3,:),...
    %     '-.b',...
    %     xxplot_k11i(3,:),P1_Lambda_plot_k11i(4,:),xxplot_k11i(4,:),P2_Lambda_plot_k11i(4,:),xxplot_k11i(4,:),...
    %     Psc_Lambda_plot_k11i(4,:),xLambda1_k11i(4,:),yLambda1_k11i(4,:),'-.b',xLambda2_k11i(4,:),yLambda2_k11i(4,:),...
    %     '-.b')

    % plot(xxplot_k11i(1,:),P1_Lambda_plot_k11i(1,:),xxplot_k11i(1,:),P2_Lambda_plot_k11i(1,:),xxplot_k11i(1,:),...
    %     Psc_Lambda_plot_k11i(1,:),xLambda1_k11i(1,:),yLambda1_k11i(1,:),'-.b',xLambda2_k11i(1,:),yLambda2_k11i(1,:),...
    %     '-.b',...
    %     xxplot_k11i(2,:),P1_Lambda_plot_k11i(2,:),xxplot_k11i(2,:),P2_Lambda_plot_k11i(2,:),xxplot_k11i(2,:),...
    %     Psc_Lambda_plot_k11i(2,:),xLambda1_k11i(2,:),yLambda1_k11i(2,:),'-.b',xLambda2_k11i(2,:),yLambda2_k11i(2,:),...
    %     '-.b',...
    %     xxplot_k11i(3,:),P1_Lambda_plot_k11i(3,:),xxplot_k11i(3,:),P2_Lambda_plot_k11i(3,:),xxplot_k11i(3,:),...
    %     Psc_Lambda_plot_k11i(3,:),xLambda1_k11i(3,:),yLambda1_k11i(3,:),'-.b',xLambda2_k11i(3,:),yLambda2_k11i(3,:),...
    %     '-.b',...
    %     xxplot_k11i(3,:),P1_Lambda_plot_k11i(4,:),xxplot_k11i(4,:),P2_Lambda_plot_k11i(4,:),xxplot_k11i(4,:),...
    %     Psc_Lambda_plot_k11i(4,:),xLambda1_k11i(4,:),yLambda1_k11i(4,:),'-.b',xLambda2_k11i(4,:),yLambda2_k11i(4,:),...
    %     '-.b',...
    %     xxplot_k11i(5,:),P1_Lambda_plot_k11i(5,:),xxplot_k11i(5,:),P2_Lambda_plot_k11i(5,:),xxplot_k11i(5,:),...
    %     Psc_Lambda_plot_k11i(5,:),xLambda1_k11i(5,:),yLambda1_k11i(5,:),'-.b',xLambda2_k11i(5,:),yLambda2_k11i(5,:),...
    %     '-.b',...
    %     xxplot_k11i(6,:),P1_Lambda_plot_k11i(6,:),xxplot_k11i(6,:),P2_Lambda_plot_k11i(6,:),xxplot_k11i(6,:),...
    %     Psc_Lambda_plot_k11i(6,:),xLambda1_k11i(6,:),yLambda1_k11i(6,:),'-.b',xLambda2_k11i(6,:),yLambda2_k11i(6,:),...
    %     '-.b')


% [Model_SVTSFast_CrossFun,Model_SVTSFast_IntervalFun,Model_SVTSFast_DerivativeFun,P_Model_SVTSFast_Intervalopt_AllFun,Psc_Lambda_3Max]
% [Lambda_Model_SVTSFast_Intervalopt_AllFun,P_Model_SVTSFast_Intervalopt_AllFun,Lambda_3_AllFun,Psc_Lambda_1Fun,Psc_Lambda_3Max,Psc_Lambda_2Fun]


end