%run Matlab lasso function for 100 column of results
load CaseInput_RSSdesign.mat
lasso_model=zeros(20,100);
lambda=23; % lambda is the column number of the chosen lambda

% Store the result in beta
for i=1:100
    beta(i).r=lasso(X,Y(:,i));
end

for i=1:100
    t=beta(i).r;
    lasso_r(:,i)=t(:,lambda);
end
    