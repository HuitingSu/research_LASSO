function [PMVNMDSimu]=MVNMDSimu(xl,xu,mu,SIGMA,Total)
% %    PP1=0;
% PP2=0;
% %    Total=10000
% for i=1:Total
%     kk1=mvnrnd(mu,SIGMA,1)';
%     kk4=(kk1<=xl | kk1>=xu);
%     %kk2=abs(kk1)>=1;
%     %[kk1,kk2,kk4];
%     %PP1=(sum(kk2)>0)+PP1;
%     PP2=(sum(kk4)>0)+PP2;
% end
% P1MVNMDSimu=(Total-PP2)/Total;
% %    (Total-PP1)/Total
% %P1mvncdf=mvncdf(xl,xu,mu,SIGMA)
% %P1MVNMDSimu-0.6570

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PPTotal=zeros(1,Total);
xlTotal=xl*ones(1,Total);
xuTotal=xu*ones(1,Total);
uTotal=mvnrnd(mu,SIGMA,Total)';
if size(uTotal,1)==1
    kkTotalcompsum=   (uTotal<=xlTotal | uTotal>=xuTotal);
else
    kkTotalcompsum=sum(uTotal<=xlTotal | uTotal>=xuTotal);
end
kkTotalcomp=kkTotalcompsum>0;
kkTotal=-kkTotalcomp+ones(1,Total);
PPTotal=kkTotal;
PMVNMDSimu=sum(PPTotal)/Total;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555

end
