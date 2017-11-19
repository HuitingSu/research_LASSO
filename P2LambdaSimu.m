function [P2_lambda]=P2LambdaSimu(Design,Beta_Current_OLS,Lambda,Sigma)
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
BetaOri=Beta_Current_OLS;
Total=100000;
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
P2mvnSimu=0;
PCSCSOri=0;
Lambda_XX=0;
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

xl2=-abs(BetaiA)'; xu2=abs(BetaiA)';
mu2=Lambda/2*inv(Xa'*Xa)*SgnBetaA;   SIGMA2=(inv(Xa'*Xa)*Sigma^2);

[DTT,p2TT]=chol(SIGMA2);
[T,err] = cholcov(SIGMA2,0);
if p2TT==0 && err==0 && size(Position0,2)<=p

    P2mvnSimu(kkk)=MVNMDSimu(xl2,xu2,mu2,SIGMA2,Total);
    if sum(mu2.*BetaiA'>0)~=size(mu2,1) && sum(mu2==0)==0
        xl2New1=xl2;
        xu2New1=xu2;
        for i=1:size(mu2,1)
            if mu2(i)*BetaiA(i)<0
                kk=i;
                xl2New1(kk)=mu2(kk)-9*SIGMA2(kk,kk);
                xu2New1(kk)=mu2(kk)+9*SIGMA2(kk,kk);
            end
        end
        P2mvnSimu(kkk)=MVNMDSimu(xl2New1,xu2New1,mu2,SIGMA2,Total);
    end

else
    P2mvnSimu(kkk)=0;
end
P2_lambda=P2mvnSimu(1:kkk);

end