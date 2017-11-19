function [P1MVNMDSimu]=MVNMDSimuXc1Xc2(xl1r,xu1r,mu1r,mu12,G,SIGMA1r,Total)
%    PP1=0;
    PP2Total=zeros(1,Total);
    xl12Total=-ones(size(mu12,1),Total); 
    xu12Total=ones(size(mu12,1),Total);
 
%    xl1rTotal=-ones(size(mu1r,1),Total); 
    xl1rTotal=xl1r*ones(1,Total); 
%    xu1rTotal=ones(size(mu1r,1),Total);
    xu1rTotal=xu1r*ones(1,Total);
    mu1rTotal=mu1r*ones(1,Total); 
    mu12Total=mu12*ones(1,Total); 
%    Total=10000
%    for i=1:Total
        %Trianglelow=chol(SIGMA1r)';
        %Trianglelow*normrnd(0,1,r,10000)>xulr or <xl1r
        %sum product;
        u1Total=mvnrnd(mu1r,SIGMA1r,Total)';
        
        kk1Totalcompsum=sum(u1Total<=xl1rTotal | u1Total>=xu1rTotal);
        kk1Totalcomp=kk1Totalcompsum>0;
        kk1Total=-kk1Totalcomp+ones(1,Total);
        %if sum(kk1Total)==0
        u2Total=mu12Total+G*(u1Total-mu1rTotal);
        if size(u2Total,1)==1
           kk2Totalcompsum=(u2Total<=xl12Total | u2Total>=xu12Total);
        else
           kk2Totalcompsum=sum(u2Total<=xl12Total | u2Total>=xu12Total);
        end
        kk2Totalcomp=kk2Totalcompsum>0;
        kk2Total=-kk2Totalcomp+ones(1,Total);
           %if sum(kk2Total)==0
         PP2Total=kk1Total.*kk2Total;
%               u1
%               u2
%               
           %end
        %end

 %   end
    P1MVNMDSimu=sum(PP2Total)/Total;
%    (Total-PP1)/Total
    %P1mvncdf=mvncdf(xl1,xu1,mu1,SIGMA1)
    %P1MVNMDSimu-0.6570
    
    
    
    
    
end


% function [P1MVNMDSimu]=MVNMDSimu(xl1,xu1,mu1,SIGMA1,Total)
% %    PP1=0;
%     PP2=0;
% %    Total=10000
%     for i=1:Total
%         kk1=mvnrnd(mu1,SIGMA1,1)';
%         kk4=(kk1<=xl1 | kk1>=xu1);
%         %kk2=abs(kk1)>=1;
%         %[kk1,kk2,kk4];
%         %PP1=(sum(kk2)>0)+PP1;
%         PP2=(sum(kk4)>0)+PP2;
%     end
%     P1MVNMDSimu=(Total-PP2)/Total;
% %    (Total-PP1)/Total
%     %P1mvncdf=mvncdf(xl1,xu1,mu1,SIGMA1)
%     %P1MVNMDSimu-0.6570
% end
