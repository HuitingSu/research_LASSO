function [P1_lambda,P2_lambda]=P1P2Lambda(Design,BetaOri,Lambda,Sigma)
%[FixOptimalLambda,FixOptimalPCSCS]=OptimalLambdaforFixBetaCDF(Design,BetaOri,Sigma,Accuracy)
% Input Beta
% Output Lambdaoptimal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input  Beta_Current;Position_Current;Lambda_Current;A_setA_onecolumn_Current;k_Current;Lambda_Start;Lambda_End;Lambda_Pace
%get Lambda_Current update;
[n,p]=size(Design);
ppp=0;  %All Pacei
aaa=0;  %All Optimal
ccc=0;  %All computer cdfmv times;
Stop_LamOpItaration=0;
Beta_AllPath=zeros(p,1);
Lambda_AllPath=0;
Beta_Optimal_All=zeros(p,1);
Lambda_Optimal_All=0;
TooSmall=0;

k_select=sum(abs(BetaOri)>0);

Position0=find(abs(BetaOri)>0)';

Position_Current=Position0;
A_setA_onecolumn_Current=Position0;
Beta_Current=BetaOri;
k_Current=k_select;
SplitNumber=15;
Lambda_Start=1;
Y=Design*BetaOri;
Lambda_End=max(abs(2*Design'*Y));
Lambda_Pace=(Lambda_End-Lambda_Start)/SplitNumber;
Lambda_End=Lambda_End+1;

P_Epsilon1Ori_Aj=0;P_Epsilon2Ori_Aj=0;P_CSCSlb=0;P_mvncdf=0;Lambda_Current=0;Beta_Current_All=zeros(p,1);
Lambda_Current_All=0;
Stop_Paceless1=0;
%while Stop_Paceless1==0
Xa=Design(:,Position_Current);
Xc=Design;
Xc(:,Position_Current)=[];
SgnBetaA=zeros(size(Position_Current,2),1);
for i_Sign=1:size(Position_Current,2)
    SgnBetaA(i_Sign,:)=sign(Beta_Current(A_setA_onecolumn_Current(i_Sign)))';
end
P=Xa*inv(Xa'*Xa)*Xa';
R=Xc'*Xa*inv(Xa'*Xa)*SgnBetaA;
I=eye(n);

%get the Lambda_Optimal0
kkk=0;
Lambda_Moving=Lambda_Start-0.01;
P0De=0;P_Epsilon2De=0;
P_Epsilon1De=0;
P_Epsilon1=0;
P_Epsilon2=0;
P_CSCS=0;
P_Epsilon1Ori=0;
P_Epsilon2Ori=0;
P1mvncdf=0;
P2mvncdf=0;
PCSCSOri=0;
Lambda_XX=0;
%while Lambda_Moving<(Lambda_End+0.01)
ccc=ccc+1;
Lambda=Lambda;
kkk=kkk+1;
Lambda_XX(kkk)=Lambda_Moving;
Mu2i=1/2*inv(Xa'*Xa)*SgnBetaA;
Mu2iOri=Lambda/2*inv(Xa'*Xa)*SgnBetaA;
Sigma1jj=zeros(1,size(Xc,2));
Sigma1jjOri=zeros(1,size(Xc,2));
aj=zeros(1,size(Xc,2));
for j=1:size(Xc,2)
    Sigma1jj_matrix=(4*Xc'*(I-P)*Xc*Sigma^2);
    Sigma1jj_matrixLambda=(4*Xc'*(I-P)*Xc*Sigma^2)/Lambda^2;
    Sigma1jj(j)=Sigma1jj_matrix(j,j)^(1/2);
    Sigma1jjOri(j)=(Sigma1jj_matrix(j,j)/Lambda^2)^(1/2);
    aj(j)=R(j)/Sigma1jj(j);
end
BetaiA=zeros(1,size(Xa,2));
Sigma2ii=zeros(1,size(Xa,2));
ai=zeros(1,size(Xa,2));
bi=zeros(1,size(Xa,2));
for i=1:size(Xa,2)
    BetaiA(i)=Beta_Current(A_setA_onecolumn_Current(i));
    Sigma2ii_matrix=(inv(Xa'*Xa)*Sigma^2);
    Sigma2ii(i)=Sigma2ii_matrix(i,i)^(1/2);
    ai(i)=Mu2i(i)/Sigma2ii(i);
    bi(i)=BetaiA(i)/Sigma2ii(i);
end

%Original format  y = mvncdf([0,0],[100,100],0,[1,0;0,1]) (xl,xu,mu,SIGMA)
%For P(epsilon1)

xl1=-ones(p-k_Current,1); xu1=ones(p-k_Current,1);
mu1=Xc'*Xa*inv(Xa'*Xa)*SgnBetaA;
SIGMA1=(4*Xc'*(I-P)*Xc*Sigma^2)/Lambda^2+0.0000000001*eye(p-k_Current);
%SIGMA1=(4*Xc'*(I-P)*Xc*Sigma^2)/Lambda^2;
%P1mvncdf=zeros(1,kkk);
%[DTT,pTT]=chol(A),if A is positive definite,pTT=0.
[T,err] = cholcov(SIGMA1,0);
[DTT,p1TT]=chol(SIGMA1);
% if p1TT==0 && err==0 && size(Position0,2)<=p
% 
%     P1mvncdf(kkk)=mvncdf(xl1,xu1,mu1,SIGMA1);
% else
%     P1mvncdf(kkk)=0;
% end
if p1TT==0 && err==0 && size(Position0,2)<p
    P1mvncdf(kkk)=mvncdf(xl1,xu1,mu1,SIGMA1);
elseif p1TT==0 && err==0 && size(Position0,2)==p
    P1mvncdf(kkk)=1;
else
    P1mvncdf(kkk)=0;
end

%P1mvncdf(kkk)=0;
%For P(epsilon2) >= P(|N2|<|BetaiA|)
xl2=-abs(BetaiA)'; xu2=abs(BetaiA)';
mu2=Lambda/2*inv(Xa'*Xa)*SgnBetaA;
SIGMA2=(inv(Xa'*Xa)*Sigma^2);

[DTT,p2TT]=chol(SIGMA2);
[T,err] = cholcov(SIGMA2,0);
if p2TT==0 && err==0 && size(Position0,2)<=p
    P2mvncdf(kkk)=mvncdf(xl2,xu2,mu2,SIGMA2);
    if sum(mu2.*BetaiA'>0)~=size(mu2,1) && sum(mu2==0)==0
        %[xl2,xu2,mu2,BetaiA'];
        xl2New1=xl2;
        xu2New1=xu2;
        for i=1:size(mu2,1)
            if mu2(i)*BetaiA(i)<0
                kk=i;
                %           xl2New(i)=-5*max(abs(mu2(i)),BetaiA(i));
                %           xu2New(i)= 5*max(abs(mu2(i)),BetaiA(i));
                xl2New1(kk)=mu2(kk)-9*SIGMA2(kk,kk);
                xu2New1(kk)=mu2(kk)+9*SIGMA2(kk,kk);
            end
        end
        P2mvncdf(kkk)=mvncdf(xl2New1,xu2New1,mu2,SIGMA2);
    end
else
    P2mvncdf(kkk)=0;
end
P1_lambda=P1mvncdf(1:kkk);
P2_lambda=P2mvncdf(1:kkk);
%          [P1_lambda,P2_lambda]
% %For P2, we can also derive a dimension reduction methods similar with P1;
%1. BetaiA*mu<0 or |BetaiA|>|mu|+9*sigma; \phi~0.99;
% [mu2,BetaiA'];
% if sum(mu2.*BetaiA'>0)~=size(mu2,1) && sum(mu2==0)==0
%
% P2mvncdf01(kkk)=mvncdf(xl2,xu2,mu2,SIGMA2);
% [xl2,xu2,mu2,BetaiA'];
% xl2New1=xl2;
% xu2New1=xu2;
% xl2New2=xl2;
% xu2New2=xu2;
%    for i=1:size(mu2,1)
%        if mu2(i)*BetaiA(i)<0
%           kk=i;
%           xl2New(i)=-5*max(abs(mu2(i)),BetaiA(i));
%           xu2New(i)= 5*max(abs(mu2(i)),BetaiA(i));
%           xl2New1(kk)=mu2(kk)-9*SIGMA2(kk,kk);
%           xu2New1(kk)=mu2(kk)+9*SIGMA2(kk,kk);
%
%        end
%    end
%
%    P2mvncdf02(kkk)=mvncdf(xl2New1,xu2New1,mu2,SIGMA2);
%    P2mvncdf03(kkk)=mvncdf(xl2New2,xu2New2,mu2,SIGMA2);
%
%    [xl2,xu2,mu2,BetaiA']
%    [xl2New1,xu2New1,mu2,BetaiA']
%    [kk,P2mvncdf01,P2mvncdf02,P2mvncdf03]
%
%    [P2mvncdf01,P2mvncdf02,P2mvncdf03]
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%P2mvncdf(kkk)=0;
%
%                     PCSCSOri(kkk)=P1mvncdf(kkk)* P2mvncdf(kkk);
%
%                     %Original lower bound format
%                     Sum1Ori=0;
%                     for Sum_j=1:size(Xc,2)
%                         Sum1Ori=Sum1Ori+(Sigma1jjOri(Sum_j))*(1/(1+R(Sum_j))*exp(-(1+R(Sum_j))^2/(2*Sigma1jjOri(Sum_j)^2))+1/(1-R(Sum_j))*exp(-(1-R(Sum_j))^2/(2*Sigma1jjOri(Sum_j)^2)));
%                 %     1/(1+R(Sum_j))*exp(-(1+R(Sum_j))^2/(2*Sigma1jjOri(Sum_j)^2))
%                 %     1/(1-R(Sum_j))*exp(-(1-R(Sum_j))^2/(2*Sigma1jjOri(Sum_j)^2))
%                 %     Sum1Ori    6
%                     end
%                     Sum2Ori=0;
%                     for Sum_i=1:size(Xa,2)
%                         Sum2Ori=Sum2Ori+(Sigma2ii(Sum_i))*(BetaiA(Sum_i)/abs(BetaiA(Sum_i)))*(1/(BetaiA(Sum_i)-Mu2iOri(Sum_i)))*exp(-(BetaiA(Sum_i)-Mu2iOri(Sum_i))^2/(2*Sigma2ii(Sum_i)^2));
%                     end
%                     P_Epsilon1Ori(kkk)=1-1/((2*pi)^(1/2)) * Sum1Ori;
%                     P_Epsilon2Ori(kkk)=1-1/((2*pi)^(1/2)) * Sum2Ori;
%                     Sum1=0;
%                     for Sum_j=1:size(Xc,2)
%                         Sum1=Sum1+(1/Sigma1jj(Sum_j))*(2*Sigma1jj(Sum_j)^2/Lambda-Lambda+(3/4*aj(Sum_j)^2+1/(4*Sigma1jj(Sum_j)^3))*Lambda^3);
% %                       Sum3=Sum3+(1/Sigma1jj(Sum_j))*(-2*Sigma1jj(Sum_j)^2/Lambda^2-1+3*(3/4*aj(Sum_j)^2+1/(4*Sigma1jj(Sum_j)^3))*Lambda^2);
%                 %       Sum1=Sum1+(1/Sigma1jj(Sum_j))*(2*Sigma1jj(Sum_j)^2/Lambda-Lambda);
%                     end
%                     Sum2=0;
%                     for Sum_i=1:size(Xa,2)
%                         Sum2=Sum2+BetaiA(Sum_i)/abs(BetaiA(Sum_i))*(1/(bi(Sum_i)-Lambda*ai(Sum_i))-1/2*(bi(Sum_i)-Lambda*ai(Sum_i)) + 1/8*(bi(Sum_i)-Lambda*ai(Sum_i))^3);
%                      %   Sum4=Sum4+BetaiA(Sum_i)/abs(BetaiA(Sum_i))*(ai(Sum_i)/(bi(Sum_i)-Lambda*ai(Sum_i))^2+ai(Sum_i)/2-3*ai(Sum_i)/8*(bi(Sum_i)-Lambda*ai(Sum_i))^2);
%                     end
%                     P_Epsilon1(kkk)=1-1/((2*pi)^(1/2)) * Sum1;
%                     P_Epsilon2(kkk)=1-1/((2*pi)^(1/2)) * Sum2;
%                     P_CSCS(kkk)=P_Epsilon1(kkk)*P_Epsilon2(kkk);
%                     [P1mvncdf(kkk),P_Epsilon1Ori(kkk),P_Epsilon1(kkk)];
%                     [P2mvncdf(kkk),P_Epsilon2Ori(kkk),P_Epsilon2(kkk)];
%                     Sum3=0;
%                     for Sum_j=1:size(Xc,2)
%                         Sum3=Sum3+(1/Sigma1jj(Sum_j))*(-2*Sigma1jj(Sum_j)^2/Lambda^2-1+3*(3/4*aj(Sum_j)^2+1/(4*Sigma1jj(Sum_j)^3))*Lambda^2);
%                     end
%                     Sum4=0;
%                     for Sum_i=1:size(Xa,2)
%                         Sum4=Sum4+BetaiA(Sum_i)/abs(BetaiA(Sum_i))*(ai(Sum_i)/(bi(Sum_i)-Lambda*ai(Sum_i))^2+ai(Sum_i)/2-3*ai(Sum_i)/8*(bi(Sum_i)-Lambda*ai(Sum_i))^2);
%                     end
%                     P_Epsilon1De(kkk)= (-1)/((2*pi)^(1/2))*Sum3;
%                     P_Epsilon2De(kkk)= (-1)/((2*pi)^(1/2))*Sum4;
%                     P0De(kkk)=P_Epsilon1De(kkk)*P_Epsilon2(kkk) + P_Epsilon2De(kkk)*P_Epsilon1(kkk);
%
%             for i=1:kkk
%                 if   P_Epsilon1Ori(i)<0
%                      P_Epsilon1Ori_Aj(i)=0;
%                 else
%                     P_Epsilon1Ori_Aj(i)=P_Epsilon1Ori(i);
%                 end
%             end
%             Start=0;
%             for i=1:kkk
%                 if   P_Epsilon2Ori(i)<0 || Start==1
%                      P_Epsilon2Ori_Aj(i)=0;
%                      Start=1;
%                 elseif  Start==0;
%                     P_Epsilon2Ori_Aj(i)=P_Epsilon2Ori(i);
%                 end
%             end
%             for i=1:kkk
%                 P_CSCSlb(i)=P_Epsilon1Ori_Aj(i)*P_Epsilon2Ori_Aj(i);
%                 P_mvncdf(i)=P1mvncdf(i)*P2mvncdf(i);
%             end
%
%             [P1mvncdf(1:kkk);P_Epsilon1Ori_Aj(1:kkk);P2mvncdf(1:kkk);P_Epsilon2Ori_Aj(1:kkk)];
%             [PCSCS_CSCSlb,Lambda_1CSCSlb]=max(P_CSCSlb(1:kkk));
%             [PCSCS_mvncdf,Lambda_1mvncdf]=max(P_mvncdf(1:kkk));


end