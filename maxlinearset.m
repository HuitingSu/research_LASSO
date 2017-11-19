function [X1,X2,G,LinearSet]=maxlinearset(X)
stop=0;
%%%%%%%%%%%%%%%%%%%%%
% A1 = randint(8, 3, [1 3]);
% A=A1'*A1;
r = rank(X); 
[n,p] = size(X);
%repeat=10000;
%combos = combntns(1:n,r);

%%%%%%%%%%%%%%%%%%%%%
% if s1(1) <= s1(2)

%[Position_matrix]=RandomVariableSelect(n,r,repeat);
%       combos=unidrnd(n,1,r)

%NormMatrixLarge= normrnd(0,1,24,1000);
      i=0;
while stop==0
      i=i+1;
      

      %B(:, :, i) = X(combos(i, 1:r), :);
      %combos=RandomVariableSelect(n,r,1);
      combos=randsample(n,r);
      
      B(:, :, i) = X(combos, :);
      if rank(B(:, :, i)) == r
            C = B(:, :, i);
            LinearSet=combos;            
            stop=1;
      end
    
end
   % else
%     for i = 1 : s2(1)
%         B(:, :, i) = A(1 : s1(1), combos(i, 1:r));
%         if rank(B(:, :, i)) == r
%             C = B(:, :, i);
%         end
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%
% for i = 1 : size(C, 1)
%     a1 = max(C(i, :));
%     for j = 1 : size(C, 2)
%         a1 = gcd(a1, C(i, j));
%     end
%     C(i, :) = C(i, :) / a1;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X1=C;
X2=X;
X2(LinearSet,:)=[];
G=(X1'\X2')';


%%%%%%%%%%%%%%%%%%%%%%%%
% function [VariablesSelected]=RandomVariableSelect(p,s,repeat)
% %p is the variables in total and s is the number of the true variables
% 
% for kkP=1:repeat
%     A1=0;
%     while size(A1,2)~=s+1
%           x=unidrnd(p);
%           if x~=0
%              if size(find(A1==x),2)==0
%                 A1=[A1,x];
%              end
%           end
%     end
%     VariablesSelected(kkP,:)=sort(A1(2:s+1));
% end





end