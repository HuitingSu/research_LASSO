function [A_setA_Path,Beta_Path,Lambda_Path]=OriLasso(Design,Y)  

    X=Design;
    [n,p]=size(X);
    Matrix_0=zeros(p);Single_column0=Matrix_0(:,1);
    Single_row0=Matrix_0(1,:);
    [A_0,A_setAll]=sort(Single_row0);
    Beta=Matrix_0(:,1);
    A_setA_Path1=zeros(1,p);
    Max_index=find(abs(2*X'*Y)==max(abs(2*X'*Y)));   
    k_Lambda=1;
    LambdaD1D2(k_Lambda)=max(abs(2*X'*Y));
    A_setA=Max_index';
    A_setA_Path1(1,1:size(A_setA,2))=A_setA;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if size(Max_index,1)>1
       Max_index;
       'There are multiple ones in the first step. Max(abs(2Xt*Y)';
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A_setAll(A_setA)=[];
    A_setAc=A_setAll;
    [A_0,A_setAll]=sort(Single_row0);
    Del_L=-2*(X'*Y-X'*X*Beta);  
    V_setA=-sign(Del_L(A_setA));
    V_0=zeros(p-size(V_setA,1));
  
    V_setAc=V_0(:,1);  
    V_setAll=Single_column0;

    for ii_A_set=1:size(A_setA,2)
        V_setAll(A_setA(ii_A_set))=V_setA(ii_A_set);
    end
    Step_i=0;
    Beta_All(:,1)=Beta;
    while max(abs(Del_L))>0.0001 && Step_i<=min(n,p)          

% Step_i
          if size(A_setAc,2)==0
             break;
          end
        C_p=(X'*Y-X'*X*Beta);
        a_p=X'*X*V_setAll;
        Index_iid1=0;

        for ii_d1=1:size(A_setA,2)
            Index_iid1=Index_iid1+1;
            Index_jjd1=0;
            for jj_d1=1:size(A_setAc,2)
                Index_jjd1=Index_jjd1+1;
                if (C_p(A_setA(ii_d1))+C_p(A_setAc(jj_d1)))/(a_p(A_setA(ii_d1))+a_p(A_setAc(jj_d1)))>0
                    d_positive(Index_jjd1,Index_iid1)=(C_p(A_setA(ii_d1))+C_p(A_setAc(jj_d1)))/(a_p(A_setA(ii_d1))+a_p(A_setAc(jj_d1)));
                else
                    d_positive(Index_jjd1,Index_iid1)=1.0e+100;
                end
                if (C_p(A_setA(ii_d1))-C_p(A_setAc(jj_d1)))/(a_p(A_setA(ii_d1))-a_p(A_setAc(jj_d1)))>0
                    d_negative(Index_jjd1,Index_iid1)=(C_p(A_setA(ii_d1))-C_p(A_setAc(jj_d1)))/(a_p(A_setA(ii_d1))-a_p(A_setAc(jj_d1)));
                else
                    d_negative(Index_jjd1,Index_iid1)=1.0e+100;
                end
                d1_All(Index_jjd1,Index_iid1)=min(d_positive(Index_jjd1,Index_iid1),d_negative(Index_jjd1,Index_iid1));  %这里的min是两个中选一个，所以无所谓多个；
            end
        end

        if size(d1_All,1)==1
           Index_d1_ij=1;
        else
           d1_ij=min(d1_All);
           Index_d1_ij=zeros(size(A_setAc,2),size(A_setA,2));  
           for ii_Index_d1_ij=1:size(A_setA,2)
               AA1=find(d1_All(:,ii_Index_d1_ij)==d1_ij(ii_Index_d1_ij));
               for ii_AA1=1:size(AA1,1)
                   Index_d1_ij(ii_AA1,ii_Index_d1_ij)=AA1(ii_AA1);
               end
               clear AA1
           end
        end

        [d1,Index_d1_i]=min(min(d1_All));
        BB1=Index_d1_ij(:,Index_d1_i);
        CC1=BB1(find(BB1~=0));    
        clear BB1 
        Index_d1=A_setAc(CC1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
         %test3
        if size(Index_d1,2)>1
            'd1 has multiple ones.';
%             d1_All;
%             Index_d1;
%             Design_F ;
%             Y_F;
%             Beta_setting_vector;
%             ii_Beta_setting;
        end
        clear d_positive d_negative d1_All d1_ij Index_d1_ij 
            for ii_d2=1:size(A_setA,2)
                if -Beta(A_setA(ii_d2))/V_setAll(A_setA(ii_d2))>0
                   d2_All(ii_d2)=-Beta(A_setA(ii_d2))/V_setAll(A_setA(ii_d2));
                else
                   d2_All(ii_d2)=1.0e+100;
                end
            end
            d2=min(d2_All);
            Index_d2_ij=find(d2_All==d2);
            Index_d2=A_setA(Index_d2_ij);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
         %test4
        if size(Index_d2_ij,2)>1 && d2<1.0e+100
            'd2 has multiple ones';
            d2_All;
            Index_d2_ij;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

            clear d2_All 

            if d1==d2
                'd1 and d2 is equal'
%                 d1;
%                 Design_F;
%                 Y_F;
%                 Beta_setting_vector;
%                 ii_Beta_setting;
                break;
            end

            d=min(d1,d2);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %(b)        %take step:Beta=Beta+d*V_setAll
            Beta=Beta+d*V_setAll;
    %(c)        
            if d==d1          %；
               A_setA=[A_setA,Index_d1]; 
            end
            if d==d2            
               A_setA(Index_d2_ij)=[];%
            end
    %(d)    Calculate new direction   C=sum(X(A,i)*X(A,i)')
    Step_i=Step_i+1;
    Beta_Path1(:,Step_i)=Beta;
    A_setA_Path1(Step_i+1,1:size(A_setA,2))=A_setA;
    k_Lambda=k_Lambda+1;
    LambdaD1D2(k_Lambda)=LambdaD1D2(k_Lambda-1)-d;

          if Step_i==1
             LambdaVectorMax=-2*(-X'*Y);
             Lambda_Path(Step_i)=abs(LambdaVectorMax(A_setA(1)));
             
          end

          if Step_i>1
             LambdaVectorMax=abs(-2*(X'*X*Beta_Path1(:,(Step_i-1))-X'*Y));
             Lambda_Path(Step_i)=abs(LambdaVectorMax(A_setA(1)));             
          end    


            X_Ai=X(:,A_setA);
            C_newdir=X_Ai'*X_Ai;
            Beta_setA=Beta(A_setA,:);  
            
            if det(C_newdir)==0 || cond(C_newdir)>1.0e+013
                Skip=0;
                Exist_NaN_F=1;           
               'Mkk is not inversable';
               %A_setA=zeros(1,1);
               break;
            end 

             V_setA=inv(C_newdir)*(-sign(Del_L(A_setA)));                                                 
             A_setAll(A_setA)=[];%
            A_setAc=A_setAll;
            [A_0,A_setAll]=sort(Single_row0);

            Del_L=-2*(X'*Y-X'*X*Beta);  %p*1
% Del_L
% Beta
% max(abs(Del_L))
% max(abs(Del_L))>0.0001

            
            V_0=zeros(p-size(V_setA,1));  
            if size(V_setA,1)~=p
                V_setAc=V_0(:,1);
            else
                V_setAc=[];
            end
            V_setAll=Single_column0;
            for ii_A_set=1:size(A_setA,2)
                V_setAll(A_setA(ii_A_set))=V_setA(ii_A_set);
            end

            
%             max(abs(Del_L))>0.0001 && Step_i<=min(n,p)  
%             Beta_setA
% %             A_setAc
%             A_setA
%             Lambda_Path
%             size(Beta_Path1)
    end

         A_setA_Path=A_setA_Path1';
         LambdaVectorMax=-2*(X'*X*Beta_Path1(:,(Step_i))-X'*Y);        
         Lambda_Path(Step_i+1)=abs(LambdaVectorMax(A_setA(1)));    
         
         Beta_Path=[zeros(p,1),Beta_Path1];
         
end