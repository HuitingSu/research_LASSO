function [Model_SVTSFast_CrossFun,Model_SVTSFast_DerivativeFun,Psc_lambda_3Max,A_setA_Path]=LassoSVTSFast_OLS_Sigma_New(X,Y,SigmaTrue)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PSC fast algorithm
%Input: X, Y,
% Model_SVTSFast_Cross=Model_SVTSFast_CrossFun;
% Model_SVTSFast_Interval=Model_SVTSFast_IntervalFun;
% Model_SVTSFast_Derivative=Model_SVTSFast_DerivativeFun;
% P_Model_SVTSFast_Intervalopt_All=P_Model_SVTSFast_Intervalopt_AllFun;
% Psc_lambda_3=Psc_lambda_3Max;
%Generate intervals from Lasso: Obtain the Lasso Path intervals; find Xa and Lambda_Current(Lambda_1=LambdaMin,Lambda_2=LambdaMax) for each interval
%For each interval:
%Obtain Beta_OLS_Current
%Obtain P2(lambda_1), P1(lambda_1),P2(lambda_2), P1(lambda_2)
%Obtain (lambda_1, P2(lambda_1)), (lambda_2, P2(lambda_2)) = (a2,-b2)
%Obtain (lambda_1, P1(lambda_1)), (lambda_2, P1(lambda_2)) = (a1,b1)
%Test 1: lambda_0 \in (lambda_1, lambda_2), (0, lambda_1), (lambda_02, \infinity)
%Test 2: lambda_opt=1/2(a2/b2-a1/b1), P_opt= (a1^2b2^2+b1^2a2^2-a1a2b1b2)/ b1b2
%Test 3: if: -b2/P2(lambda_1)+b1/P2(lambda_1)> 0, \exit lambda_3 \in (lambda_1,lambba_2) s.t. Psc(lambda_3)>Psc(lambda_1)
%        if:  b2/P2(lambda_2)-b1/P2(lambda_2)> 0, \exit lambda_3 \in (lambda_1,lambba_2) s.t. Psc(lambda_3)>Psc(lambda_2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%,%%%%%%%%%%55
iik=1;
Design=X;
Pace=0.0001;
[n,p]=size(X);
% [A_setA_Path,Beta_Path,Lambda_Path]=OriLasso(X,Y);     % Lambda_Path=[75.6656  71.2099 45.7350 37.6293 ...]

[Blasso,Slasso]=lasso(X,Y,'NumLambda',25);
Lambda_Path=fliplr(Slasso.Lambda).*50;
Beta_Path=fliplr(Blasso);
A_setA_Path=zeros(size(Beta_Path,1),size(Beta_Path,2));
tail=0;
for j=2:size(A_setA_Path,2)
    A_setA_Path(:,j)=A_setA_Path(:,j-1);
    for i=1:size(A_setA_Path,1)
        if Beta_Path(i,j)~=0 && Beta_Path(i,j-1)==0
            A_setA_Path(tail+1,j)=i;
            tail=tail+1;
        end
    end
end


%         A_setA_Path(1:6,1:6)
%Maxmodel=min(size(A_setA_Path,2)-6,n/2,8);
%Maxmodel=min(min(floor(n/2)+2,floor(p/3)+2),size(Lambda_Path,2)-2);
%Maxmodel=min(floor(n/2)+2,floor(p/3)+2)
Maxmodel=size(A_setA_Path,2)-1;
Selectmodel=zeros(p,1);
Model_SVTSFast_CrossFun=zeros(p,1);
Model_SVTSFast_IntervalFun=zeros(p,1);
Model_SVTSFast_DerivativeFun=zeros(p,1);
P_Model_SVTSFast_Intervalopt_AllFun=zeros(p,1);
Psc_lambda_3Max=zeros(p,1);
Lambda_opt_k11i=zeros(p,1);
Psc_lambdaopt_Max=zeros(p,1);
DistanceofOptLambda=zeros(p,1);

Lambda_SVOE_OLS=0;
PCSCS_SVOE_OLS=0;
Lambda_3_AllFun=zeros(p,1);
Lambda_Model_SVTSFast_Intervalopt_AllFun=zeros(p,1);
P1_lambda_All=zeros(p,1);
P2_lambda_All=zeros(p,1);
Psc_lambda_3Max=zeros(p,1);
Psc_lambda_1Fun=zeros(p,1);
Psc_lambda_2Fun=zeros(p,1);

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
%     Sigmakki1(kki1)=(1/n*(Y-Xkki1*inv(Xkki1'*Xkki1)*Xkki1'*Y)'*(Y-Xkki1*(Xkki1'*Xkki1)*Xkki1'*Y))^(1/2);
% end
%         Sigmakki1
%         SigmaBeta2
Sigma=SigmaTrue;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lambda_Path_Sort=sort(Lambda_Path);
Lambda_Current=Lambda_Path_Sort(size(Lambda_Path,2)-1);
for k11_i=2:Maxmodel

    Lambda_Current=Lambda_Current-Pace;         % This seems not used --Huiting
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
    LambdaPath_All(iik,k11_i)=Lambda_Path(k11_i);
    Lambda_2=Lambda_Path(iik,k11_i);   % Lambda_2 and Lambda_1 are both from the Lambda_Path
    if k11_i+1>size(Lambda_Path,2)
        Lambda_1=0.0001;
    else
        Lambda_1=Lambda_Path(iik,k11_i+1);
    end

    if  p<15  && (size(find(Beta_Current_OLS~=0),1)<6)
        [P1_lambda_1,P2_lambda_1]=P1P2Lambda(Design,Beta_Current_OLS,Lambda_1,Sigma);
        [P1_lambda_2,P2_lambda_2]=P1P2Lambda(Design,Beta_Current_OLS,Lambda_2,Sigma);
        Psc_lambda_1Fun(k11_i)=P1_lambda_1*P2_lambda_1;
        Psc_lambda_2Fun(k11_i)=P1_lambda_2*P2_lambda_2;
    elseif p>15  && (size(find(Beta_Current_OLS~=0),1)<6)
        [P1_lambda_1]=P1LambdaSimu(Design,Beta_Current_OLS,Lambda_1,Sigma);
        [P1_lambda_2]=P1LambdaSimu(Design,Beta_Current_OLS,Lambda_2,Sigma);
        [P2_lambda_1]=P2LambdaCDF(Design,Beta_Current_OLS,Lambda_1,Sigma);
        [P2_lambda_2]=P2LambdaCDF(Design,Beta_Current_OLS,Lambda_2,Sigma);
        Psc_lambda_1Fun(k11_i)=P1_lambda_1*P2_lambda_1;
        Psc_lambda_2Fun(k11_i)=P1_lambda_2*P2_lambda_2;
    else
        [P1_lambda_1]=P1LambdaSimu(Design,Beta_Current_OLS,Lambda_1,Sigma);
        [P1_lambda_2]=P1LambdaSimu(Design,Beta_Current_OLS,Lambda_2,Sigma);
        [P2_lambda_1]=P2LambdaSimu(Design,Beta_Current_OLS,Lambda_1,Sigma);
        [P2_lambda_2]=P2LambdaSimu(Design,Beta_Current_OLS,Lambda_2,Sigma);
        Psc_lambda_1Fun(k11_i)=P1_lambda_1*P2_lambda_1;
        Psc_lambda_2Fun(k11_i)=P1_lambda_2*P2_lambda_2;
    end
    b2=(P2_lambda_1-P2_lambda_2)/(Lambda_2-Lambda_1);
    a2=P2_lambda_2+b2*Lambda_2;
    b1=(P1_lambda_2-P1_lambda_1)/(Lambda_2-Lambda_1);
    a1=P1_lambda_2-b1*Lambda_2;
    [k11_i,P1_lambda_1,P2_lambda_1,P1_lambda_2,P2_lambda_2];
    [Lambda_1,P1_lambda_1,Lambda_2,P1_lambda_2,a1,b1;Lambda_1,P2_lambda_1,Lambda_2,P2_lambda_2,a2,b2];
    %Test 1: lambda_0 \in (lambda_1, lambda_2), (0, lambda_1), (lambda_02, \infinity)
    if (P2_lambda_1-P1_lambda_1)*(P2_lambda_2-P1_lambda_2)<0
        Model_SVTSFast_CrossFun(k11_i)=1;
    elseif (P2_lambda_2-P1_lambda_2)>0
        Model_SVTSFast_CrossFun(k11_i)=0;
    elseif (P2_lambda_2-P1_lambda_2)<0
        Model_SVTSFast_CrossFun(k11_i)=2;
    end
    %Test 4: if new Psc(lambda_3)>max(Psc(lambda_2),Psc(lambda_1)),then
    k3=0;
    P1_lambda_1_tem=P1_lambda_1;
    P2_lambda_1_tem=P2_lambda_1;
    Lambda_1_tem=Lambda_1;
    P1_lambda_2_tem=P1_lambda_2;
    P2_lambda_2_tem=P2_lambda_2;
    Lambda_2_tem=Lambda_2;
    Psc_lambda_1_tem=P1_lambda_1_tem*P2_lambda_1_tem;
    Psc_lambda_2_tem=P1_lambda_2_tem*P2_lambda_2_tem;
    Model_SVTSFast_DerivativeFun(k11_i)=0;
    stop=0;
    while stop==0
        k3=k3+1;
        Lambda_3=(P1_lambda_1_tem*P2_lambda_1_tem*Lambda_1_tem+P1_lambda_2_tem*P2_lambda_2_tem*Lambda_2_tem)/(P1_lambda_1_tem*P2_lambda_1_tem+P1_lambda_2_tem*P2_lambda_2_tem)
        if  p<15  && (size(find(Beta_Current_OLS~=0),1)<6)
            [P1_lambda_3,P2_lambda_3]=P1P2Lambda(Design,Beta_Current_OLS,Lambda_3,Sigma);
            Psc_lambda_3=P1_lambda_3*P2_lambda_3;
            Psc_lambda_3Max(k11_i)=max([Psc_lambda_3,Psc_lambda_1_tem,Psc_lambda_2_tem]);
        elseif p>15  && (size(find(Beta_Current_OLS~=0),1)<5)
            [P2_lambda_3]=P2LambdaCDF(Design,Beta_Current_OLS,Lambda_3,Sigma);
            [P1_lambda_3]=P1LambdaSimu(Design,Beta_Current_OLS,Lambda_3,Sigma);
            Psc_lambda_3=P1_lambda_3*P2_lambda_3;
            Psc_lambda_3Max(k11_i)=max([Psc_lambda_3,Psc_lambda_1_tem,Psc_lambda_2_tem]);
        else
            [P2_lambda_3]=P2LambdaSimu(Design,Beta_Current_OLS,Lambda_3,Sigma);
            [P1_lambda_3]=P1LambdaSimu(Design,Beta_Current_OLS,Lambda_3,Sigma);
            Psc_lambda_3=P1_lambda_3*P2_lambda_3;
            Psc_lambda_3Max(k11_i)=max([Psc_lambda_3,Psc_lambda_1_tem,Psc_lambda_2_tem]);
        end

        if Psc_lambda_3>=max(Psc_lambda_1_tem, Psc_lambda_2_tem) %|| (min([Psc_lambda_3,Psc_lambda_1_tem,Psc_lambda_2_tem])>0.90 && k11_i~=1)
            Model_SVTSFast_DerivativeFun(k11_i)=1;
            stop=1;
        elseif Psc_lambda_1_tem==1 && Psc_lambda_2_tem==1%|| (min([Psc_lambda_3,Psc_lambda_1_tem,Psc_lambda_2_tem])>0.90 && k11_i~=1)
            Model_SVTSFast_DerivativeFun(k11_i)=1;
            stop=1;
        else
            Lambda_triple=[Lambda_1_tem,Lambda_2_tem,Lambda_3];
            Psc_triple=[Psc_lambda_1_tem,Psc_lambda_2_tem,Psc_lambda_3];
            P1_triple=[P1_lambda_1_tem,P1_lambda_2_tem,P1_lambda_3];
            P2_triple=[P2_lambda_1_tem,P2_lambda_2_tem,P2_lambda_3];
            [aasort,bbsort]=sort(Psc_triple);
            Lambda_1_tem=Lambda_triple(bbsort(2));
            P1_lambda_1_tem=P1_triple(bbsort(2));
            P2_lambda_1_tem=P2_triple(bbsort(2));

            Lambda_2_tem=Lambda_triple(bbsort(3));
            P1_lambda_2_tem=P1_triple(bbsort(3));
            P2_lambda_2_tem=P2_triple(bbsort(3));

            Psc_lambda_1_tem=P1_lambda_1_tem*P2_lambda_1_tem;
            Psc_lambda_2_tem=P1_lambda_2_tem*P2_lambda_2_tem;

            %                         Lambda_length=abs(Lambda_2_tem-Lambda_1_tem);
            %                         Psc_length=abs(Psc_lambda_1_tem-Psc_lambda_2_tem);
            if k3>1 || max([Psc_lambda_3,Psc_lambda_1_tem,Psc_lambda_2_tem])<0.1
                stop=1;
            end
        end
    end
    Lambda_3_AllFun(k11_i)=Lambda_3;
    Lambda_3;   % add by Huiting 
    % k11_i
end

