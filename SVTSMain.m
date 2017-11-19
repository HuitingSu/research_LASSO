% The main function for Self-voting LASSO. 
% SigmaTrue is generally set to be 1.
function [Model_SVTSFast_Actual,fixed]=SVTSMain(X,Y,SigmaTrue)
TypeMulti=1;
[n,p]=size(X);
[A_setA_Path,Beta_Path,Lambda_Path]=OriLasso(X,Y);
MaxmodelNumber=min(min(floor(n/2)+2,floor(p/3)+2),size(Lambda_Path,2)-1);
if n<15
   MaxmodelNumber=(size(Lambda_Path,2)-2);
end
[Model_SVTSFast_Cross,Model_SVTSFast_Derivative,Psc_lambda_3,A_setA_Path]=LassoSVTSFast_OLS_Sigma_New(X,Y,SigmaTrue);

if sum(Model_SVTSFast_Derivative(2:size(Model_SVTSFast_Derivative,1)))==0
    fixed=1;

    [Model_NoCross,Psc_opt_currentmodel,Lambda_opt_currentmodel]=LassoSVTSFast_NoCross(X,Y,SigmaTrue);
    Model_SVTSFast_Derivative=Model_NoCross;
else
    fixed=0;
end

Test03Options=find(Model_SVTSFast_Derivative==1);
if TypeMulti==1
%max Psc
[Test03a2,Test03b]=max(Psc_lambda_3(Test03Options));
elseif TypeMulti==2
%largest model
Test03b=Test03Options(size(Test03Options,1));
elseif TypeMulti==3
%Psc*k
[Test03a2,Test03b]=max(Psc_lambda_3(Test03Options)*Test03Options);
end

k_SVTSHit=Test03Options(Test03b);
if size(k_SVTSHit,1)==0
    Model_SVTSFast_Actual=zeros(p,1);
else
    Model_SVTSFast_Actual=A_setA_Path(:,k_SVTSHit);
end



end