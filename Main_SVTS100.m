%run Selfvoting for 100 column of results
Model_SVTSFast_Actual=zeros(20,100);
fixed=zeros(100);
load CaseInput_RSSdesign.mat

% For some Y(:,i), the selfvoting may fail. 
% In these cases, comment the first 2 lines, change i to the fail interation,
% then run again untial get the answer.
for i=1:100
    [Model_SVTSFast_Actual(:,i),fixed(i)]=SVTSMain(X,Y(:,i),1); 
end