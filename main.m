clear,clc;
addpath './OBDExtraction';
addpath './DataClass';
addpath './Feature'
addpath './Net'
videoName = {'ID002_T001','ID002_T002','ID002_T003', '118_07182017','023','028','112_07172017','106_07142017'};
data = [];
label = [];
frame = [];
total=[];
for i =1:length(videoName)
    Mdata = Maneuverdata(videoName{i});
    Mdata = Mdata.trainData();
    %save(sprintf('%s.mat', Mdata.dataID),'Mdata');
   % showDate(Mdata);
    ShowDemo(Mdata.dataID,Mdata.logField,Mdata.eventField,Mdata.trainingFrame,Mdata.trainingLabel )
    data = [data;Mdata.trainingData];
    label = [label;Mdata.trainingLabel];
    %frame = [label;Mdata.trainingLabel];
end
%%%%%%%%
% data = cellfun(@(x) Normalize(x), data,'UniformOutput',false);

[trainD,trainL,testD, testL] = SplitData(data,label,0.75);
net = ManeuversNet(trainD, trainL,1000);
y = classify(net,testD );
accuracy =  sum(y == testL)/numel(testL);
save('net.mat','net');





