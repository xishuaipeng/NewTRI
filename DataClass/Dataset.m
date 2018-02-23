classdef  Dataset
    
properties
logPath; %l 
labelPath;
videoPath;
segmentPath;
dataID;
logData;
sampleData;
logTable;
timeDelayforVideo;
logField;
eventField;
eventLabel;
eventFrame;
samplingProperty;
samplingStep;
frameRate;
logRate;
end

methods(Static)
    function sequence = table2sequence(table, window, step)
        [r,c] = size(table);
        numSequence = floor((r-window)/step)+1;
        sequence = cell(numSequence,1);     
        for i= 1:numSequence
            start = (i-1)*step;
            last = start + window;
            seq = table(start+1:last, :);
            sequence{i} = seq;
        end
    end
     
     
 end
methods
    function obj = Dataset(dataID,feature_field, event_field)
    % load the default path and data 
        obj.dataID = dataID;
        obj.logPath = sprintf('./input/raw_data/%s/%s_datalog.Csv',dataID,dataID);
        obj.labelPath = sprintf('./input/raw_data/%s/%s_labelresult.Csv',dataID,dataID);
        obj.videoPath = sprintf('./input/raw_data/%s/%s_video.avi',dataID,dataID);%avi
        obj.segmentPath = sprintf('./input/raw_data/%s/%s_labelresult.Csv',dataID,dataID);
        obj.logField = feature_field;
        obj.eventField = event_field;
        obj.timeDelayforVideo = 0;
        obj.frameRate = 29.9;
    end
    
    function obj = readLogdata(varargin)
        obj = varargin{1};
        obj.logTable  = TRI_GPS_extract_oneTrip(obj.logPath);
        obj.logTable = obj.logTable(:,obj.logField);
        %table 2 data
        obj.logData  = table();
        for i = 1:size(obj.logTable ,2)
            obj.logData (:,i) = table(str2num(char(table2array(obj.logTable (:,i)))));
         end
         obj.logData .Properties = obj.logTable .Properties;
         obj.logRate =  obj.logData.time(2) - obj.logData.time(1) ;
        %obj.logData = table2data(obj.logTable);   
        %obj.maxFrame = round(timeParser(obj.logData.time(end)) * obj.frameRate);
    end
    
    function obj = readEventlabel(varargin)
          obj = varargin{1};
          obj.eventLabel  = readtable(obj.labelPath );
          eventField =  obj.eventField;
          num_case = size(obj.eventLabel,1);
          obj.eventFrame = zeros(num_case, 2+length(eventField));
          for i=1: num_case
            startIndex = round(timeParser(obj.eventLabel.StartTime{i})* obj.frameRate);
            endIndex =  round(timeParser(obj.eventLabel.EndTime{i})* obj.frameRate);
            event= obj.eventLabel{i,eventField};
            if sum(event)>0
                %event_class = find(event>0);
               % frameLabel(startIndex:endIndex,:) =  repmat(event_class,endIndex - startIndex+1,1);
                obj.eventFrame(i,:) = [startIndex,endIndex, event ];
            end
          end
          %obj.sampleData.label = frameLabel(obj.sampleData.frame,:);
          obj.eventFrame(sum(obj.eventFrame,2)==0,:)=[];
    end
%     function obj = Frame2Logdata(varargin)
%             obj = varargin{1};
%             [num_case, num_feature] = size(obj.logData);
%             frameIndex = zeros(num_case,1);
%             for  i=1:num_case
%              frameIndex(i)=  round((obj.logData.time(i) - obj.timeDelayforVideo )* obj.frameRate )+1;
%             end   
%             obj.logData.frame = frameIndex;
%        
%     end     
    
    function obj = resampLogdata( varargin)
         % second: segWindow
         %third: segStep
         %forth: sampStep;
         %fifth: sampProperty
        obj = varargin{1};
        if isempty( obj.logData)
            obj = obj.readLogdata();
        end
        obj.samplingProperty = varargin{2};
        obj.samplingStep = varargin{3};
        [num_case, num_feature] = size(obj.logData);
        segVector = obj.logData.(obj.samplingProperty);
        max_sample = 1 + floor((segVector(end) - segVector(1)) / obj.samplingStep );
         %find the time field index
        timeIndex = cellfun(@(x) strcmp(x,'time'),obj.logData.Properties.VariableNames);
        timeIndex = find(timeIndex==1);
        %construct the array
        tempData = zeros(max_sample, num_feature);
        frameIndex = zeros(max_sample, 1);
        interValue = [segVector(1) : obj.samplingStep : segVector(end)];
        tempData(1,:) = obj.logData{1,:};
        frameIndex(1)=  round((tempData(1,timeIndex)  - obj.timeDelayforVideo )* obj.frameRate )+1;
        index = 1;  
        for  i=2:num_case
             current =interValue(index) + obj.samplingStep;%tempData(index, segIndex) 
             if obj.logData.( obj.samplingProperty )(i) >= current
                 x0 = obj.logData.( obj.samplingProperty )(i-1);
                 x1 = obj.logData.( obj.samplingProperty )(i);
                 x = current;
                 theta = (x-x0)/(x1-x0);
                 index = index + 1; 
                 tempData(index,1:end) = obj.logData{i-1,1:end}+ theta * (obj.logData{i,1:end}-obj.logData{i-1,1:end});  
                 %transfer time to frame index
                 frameIndex(index)=  round((tempData(index,timeIndex)  - obj.timeDelayforVideo )* obj.frameRate )+1;
             end   
        end  
        obj.sampleData = array2table(tempData);
        obj.sampleData .Properties = obj.logData.Properties;
        obj.sampleData .frame = frameIndex;
      end

    function obj = labelSampledata(varargin)
          obj = varargin{1};
          if isempty(obj.sampleData)
              'no sampleData please generate it first'   
          end
          if  isempty(obj.eventLabel)
              obj = obj.readEventlabel();
              
          end
          maxFrame = obj.sampleData.frame(end);
          tempFLabel = zeros(maxFrame,1);
          num_case = size(obj.eventFrame,1);
%           eventFrame = zeros(num_case, 2+length(eventField));
          for i=1: num_case
            startIndex = obj.eventFrame(i,1);
            endIndex =  obj.eventFrame(i,2);
            event= obj.eventFrame(i,3:end);
            event = find(event>0);
            tempFLabel(startIndex:endIndex,:) =  repmat(event,endIndex - startIndex+1,1);
          end
          obj.sampleData.label = tempFLabel(obj.sampleData.frame,:);
    end
    
     function obj = reSync(varargin)
        obj = varargin{1};
        if numel(varargin)==1
            obj.timeDelayforVideo = checkIfSync(obj.logData);
        elseif numel(varargin) == 4
            obj.timeDelayforVideo = checkIfSync(obj.logData,varargin{2},varargin{3},varargin{4});
        elseif numel(varargin) == 2
            obj.timeDelayforVideo = checkIfSync(obj.logData, obj.dataID, varargin{2});
        end
     end
    
    function obj = checkLabel(varargin)
      obj = varargin{1};
      x = obj.sampleData.GPS_long;
      y = obj.sampleData.GPS_lat;
      maxLength = length(x);
      label = obj.sampleData.label+1;
      %label =label - min(label);
      classes = max(label);
      radio = jet( classes);
      figure; hold on;
      for i = 1: classes
          index = find(label == i);
          p_x = x(index);
          p_y = y(index);
          color = radio(i,:);
          scatter(p_x,p_y, 'MarkerFaceColor',color,'MarkerEdgeColor',color);
          
      end
     legend(['No event',obj.eventField]);
      
      
      
%       maxLength = length(x);
%       event_start = round((obj.timeParser(obj.eventLabel.StartTime )+ obj.timeDelayforVideo ) * logRate)+1;
%       event_end =  round((obj.timeParser(obj.eventLabel.EndTime) + obj.timeDelayforVideo )* logRate)+1;
    %           event_time =[ event_start, event_end];
%       eventNum = length(event_start);
%       figure();hold on;
%       plot(x,y,'.k');
%       colorBar = ['r', 'g', 'y', 'b'];
%       event = table2array(obj.eventLabel(:,obj.eventField));
%       
%       
%       
%        for i =1 : eventNum
%            hold on;
%            label_index = event(i,:)
%            label_index = find(label_index==1);
%            if ~isempty(label_index)
%                obj.frameLabel(event_start(i):event_end(i))= obj.eventField{label_index};
%                plot(x(event_start(i):event_end(i)),y(event_start(i):event_end(i)),'.', ...
%         'Color', colorBar(label_index), 'MarkerSize', 6);
%            end
%        end               
    end 
     
     
%     function obj = FrameVgg19(varargin)
%            obj = varargin{1};
%            videoObj = VideoReader(obj.videoPath);
%            net = vgg19;
%            sz = net.Layers(1).InputSize;
%            index = 1;
%            obj.frameVgg19 = zeros(obj.maxFrame, 4096);
%            while hasFrame(videoObj)
%                 frame = readFrame(videoObj);
%                 frame = imresize(frame,[sz(1),sz(2)]);
%                 frame_feature = activations(net,frame,42);
%                 obj.frameVgg19(index,:) = frame_feature;
%                 index = index +1 ;
%            end
%     end 

%     
%      function obj = segtrip( varargin)
%          % second: segWindow
%          %third: segStep
%          %forth: sampStep;
%          %fifth: sampProperty
%         obj = varargin{1};
%         segType = varargin{5};
%         obj.samplingProperty = segType;
%         obj.samplingStep = varargin{4};
%         [num_case, num_feature] = size(obj.logData);
%         segVector = obj.logData.(segType);
%         max_sample = 1 + floor((segVector(end) - segVector(1)) / obj.samplingStep );
%         tempData = zeros(max_sample, num_feature);
%         interValue = [segVector(1) : obj.samplingStep : segVector(end)];
%         tempData(1,:) = obj.logData{1,:};
%         index = 1;  
%         for  i=2:num_case
%              current =interValue(index) + obj.samplingStep;%tempData(index, segIndex) 
%              if obj.logData.(segType)(i) >= current
%                  x0 = obj.logData.(segType)(i-1);
%                  x1 = obj.logData.(segType)(i);
%                  x = current;
%                  theta = (x-x0)/(x1-x0);
%                  index = index + 1; 
%                  tempData(index,1:end) = obj.logData{i-1,1:end}+ theta * (obj.logData{i,1:end}-obj.logData{i-1,1:end});     
%              end   
%         end  
%         %segment by distance
%         obj.segWindow = varargin{2};
%         obj.segWindow = floor(obj.segWindow / obj.samplingStep );
%         obj.segStep = varargin{3};
%         obj.segStep = floor(obj.segStep / obj.samplingStep );
%         fps = obj.frameRate ;
%         max_length = size(tempData,1);
%         index = 1;
%         begin = 1;
%         last = begin +  obj.segWindow;
%         while(last< max_length)
%             obj.segData(index).data = [ tempData(begin:last ,:)];
%             obj.segData(index).data = array2table(obj.segData(index).data);
%             obj.segData(index).data.Properties = obj.logData.Properties;
%     %             obj.segData(index).data.name = {};
%             frame =  round(  (obj.segData(index).data.time  - obj.timeDelayforVideo )* fps )  +1;
%             obj.segData(index).frameIndex = frame;
%             obj.segData(index).minFrame = min(frame(:));
%             obj.segData(index).maxFrame = max(frame(:));
%             begin = begin +  obj.segStep ;
%             last = begin +  obj.segWindow;
%             index = index + 1;
%         end
%       end
    
%%
%    function obj = checkLabel(varargin)
%           obj = varargin{1};
%           x = obj.logData.GPS_long;
%           y = obj.logData.GPS_lat;
%           maxLength = length(x);
%           event_start = round((obj.timeParser(obj.eventLabel.StartTime )+ obj.timeDelayforVideo ) * logRate)+1;
%           event_end =  round((obj.timeParser(obj.eventLabel.EndTime) + obj.timeDelayforVideo )* logRate)+1;
% %           event_time =[ event_start, event_end];
%           eventNum = length(event_start);
%           figure();hold on;
%           plot(x,y,'.k');
%           colorBar = ['r', 'g', 'y', 'b'];
%            for i =1 : eventNum
%                 hold on;
%                label_index = table2array(obj.eventLabel(i,obj.eventField));
%                label_index = find(label_index==1);
%                if ~isempty(label_index)
%                obj.frameLabel(event_start(i):event_end(i))= obj.eventField{label_index};
%                plot(x(event_start(i):event_end(i)),y(event_start(i):event_end(i)),'.', ...
%         'Color', colorBar(label_index), 'MarkerSize', 6);
%                end
%            end               
%    end 
%    function obj = appendLabel(varargin)
%           obj = varargin{1};
%           eventField =  obj.eventField;
%           segNum = size(obj.segData, 2);
%           data_time = [obj.segData.minFrame, obj.segData.maxFrame];
%           data_time = reshape(data_time,[segNum,2]);
%           event_start = obj.timeParser(obj.eventLabel.StartTime)* obj.frameRate;
%           event_end =  obj.timeParser(obj.eventLabel.EndTime)* obj.frameRate;
%           event_time =[ event_start, event_end];
%           eventNum = length(event_start);
%           for i= 1:segNum
%             for j = 1:eventNum
%                 if data_time(i,2) < event_time(j,1)
%                     break;
%                 end
%                 start = max(data_time(i,1), event_time(j,1));
%                 last = min(data_time(i,2), event_time(j,2));
%                 duration = last - start;
%                 radio_event = duration /(event_time(j,2) - event_time(j,1));
% 				radio_seg  =  duration /(data_time(i,2) - data_time(i,1));
%                 radio = max(radio_event,radio_seg);
%                 if radio > 0.7
%                     tLabel = obj.eventLabel(j,eventField);
%                     tLabel = table2array(tLabel);
%                     label_index = find(tLabel==1);
%                     if isempty(label_index)
%                         obj.segData(i).Label = obj.negativeField ;
%                     else
%                         obj.segData(i).Label = eventField{label_index};
%                      break;
%                     end
%                 else
%                    obj.segData(i).Label = obj.negativeField ;
%                 end
% 
%             end
%           end
%    end 
%    function obj = checkSeg(varargin)
%            obj = varargin{1};
%            color = '.k';
%            figure;hold on;
%            numSample = size(obj.segData,2);
%            x_event = [];
%            y_event = [];
%            event_color = color;
%            colorBar = ['.r'; '.g'; '.y'; '.b'];
%            for i = 1: numSample
%                x = obj.segData(i).data.GPS_long;
%                y = obj.segData(i).data.GPS_lat;
%                z = obj.segData(i).Label;  
%                if strcmp(color,'.k')
%                    color = '.w';
%                else
%                    color = '.k';
%                end
%                if any(strcmp(z,obj.eventField))
%                    event_index = find(strcmp(z,obj.eventField) == 1);
%                    x_event = x;
%                    y_event = y;
% %                    x_event = [x_event; x];
% %                    y_event = [y_event; y];
%                    event_color = colorBar(event_index,:);
%                end
%                plot(x,y,color);
%                plot(x_event,y_event,event_color);
%            end  
% %             plot(x_event,y_event,'.y');
% % plot_google_map  
%    end  
%    function obj = extractDrifting(varargin)
%            obj = varargin{1};
%            segNum = size(obj.segData,2);
%            videoObj = VideoReader(obj.videoPath);
%            maxFrame = videoObj.NumberOfFrames;
%            startIndex = 1;
%            lastIndex = maxFrame;
%            vstep = 10;
%            numImages = length([startIndex+1 : vstep : lastIndex]);
%            preFrame = read(videoObj, startIndex);
%            preFrame = rgb2gray(preFrame);
%            prePoints = detectSURFFeatures(preFrame);
%            [preFeatures,prePoints] = extractFeatures(preFrame,prePoints);
%            tforms(numImages) = projective2d(eye(3));
%            tIndex = 2;
%            for index = startIndex+vstep : vstep : lastIndex
%                curFrame = read(videoObj, index);
%                curFrame = rgb2gray(curFrame);
%                curPoints = detectSURFFeatures(curFrame); 
%                [curFeatures,curPoints] = extractFeatures(curFrame,curPoints);
%                matchedIndex = matchFeatures(curFeatures, preFeatures,'Unique', true);
%                curMatched = curPoints(matchedIndex(:,1),:);
%                preMatched = prePoints(matchedIndex(:,2),:);
%                    % Estimate the transformation between I(n) and I(n-1).
%                tforms(tIndex) = estimateGeometricTransform(curMatched, preMatched,...
%                     'projective', 'Confidence', 99.9, 'MaxNumTrials', 2000);
%                tforms(tIndex).T  =   tforms(tIndex).T * tforms(tIndex-1).T;
%                tIndex = tIndex+1; 
%                prePoints = curPoints;
%                preFeatures = curFeatures;
%            end
%            imageSize = size(preFrame);  % all the images are the same size
%             % Compute the output limits  for each transform
%             for i = 1:numel(tforms)
%                 [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(2)], [1 imageSize(1)]);
%             end
%            avgXLim = mean(xlim, 2);
%            avgYLim = mean(ylim, 2);  
% %            [~, idx] = sort(avgXLim);
% %            centerIdx = floor((numel(tforms)+1)/2);
% %            centerImageIdx = idx(centerIdx);
% %            Tinv = invert(tforms(centerImageIdx));
% %            for i = 1:numel(tforms)
% %                 tforms(i).T = tforms(i).T * Tinv.T;
% %            end
% %            for i = 1:numel(tforms)
% %                 [xlim(i,:), ylim(i,:)] = outputLimits(tforms(i), [1 imageSize(2)], [1 imageSize(1)]);
% %             end
% %             % Find the minimum and maximum output limits
% %             xMin = min([1; xlim(:)]);
% %             xMax = max([imageSize(2); xlim(:)]);
% % 
% %             yMin = min([1; ylim(:)]);
% %             yMax = max([imageSize(1); ylim(:)]);
% % 
% %             % Width and height of panorama.
% %             width  = round(xMax - xMin);
% %             height = round(yMax - yMin);
% % 
% %             % Initialize the "empty" panorama.
% %             Frame = read(videoObj, 1);
% %             panorama = zeros([height width 3], 'like', Frame);
% %             blender = vision.AlphaBlender('Operation', 'Binary mask', ...
% %             'MaskSource', 'Input port');
% %             % Create a 2-D spatial reference object defining the size of the panorama.
% %             xLimits = [xMin xMax];
% %             yLimits = [yMin yMax];
% %             panoramaView = imref2d([height width], xLimits, yLimits);
% %             % Create the panorama.
% %             tIndex = 1;
% %             for i = startIndex : vstep : lastIndex
% %                 I = read(videoObj, i);
% % %                 I = rgb2gray(I);
% %                 % Transform I into the panorama.
% %                 warpedImage = imwarp(I, tforms(tIndex), 'OutputView', panoramaView);
% %                 % Generate a binary mask.
% %                 mask = imwarp(true(size(I,1),size(I,2)), tforms(tIndex), 'OutputView', panoramaView);
% %                 % Overlay the warpedImage onto the panorama.
% %                 panorama = step(blender, panorama, warpedImage, mask);
% %                 tIndex = tIndex + 1;
% %             end
% %             figure
% %             imshow(panorama)
%    end
%    function obj = extractCurvature(varargin)
%            obj = varargin{1};
%            width = 64;
%            height = 64;
%            segNum = size(obj.segData,2);
%            stride = 5;
%            for i = 1:segNum
%                 x = obj.segData(i).data.GPS_lat;
%                 y = obj.segData(i).data.GPS_long;
%                 num = length(x);
%                 minx = min(x(:));
%                 maxx = max(x(:));
%                 miny =  min(y(:));
%                 maxy =  max(y(:));
%                 max_range = max((maxx -minx ),(maxy-miny));
%                 x = stride+(width-2*stride) *(x-minx)/max_range;%(maxx-minx);
%                 y = stride+(height-2*stride) *(y-miny)/max_range;%(maxy-miny);
%                 p = polyfit(x , y, 2);
%                 x_ = [1:width]';
%                 y_= polyval(p,x);
%                 %compute the curvature
%                 pp1 = 2 * p(1) * y_ + p(2);
%                 pp2 = 2 * p(1);
%                 k = abs(pp2) ./ sqrt((1+pp1.^2).^3);
%                 
%                 u = minx + ( x - stride ) * max_range /(width - 2 * stride) ;
%                 v = miny + ( y_ - stride) * max_range /(height - 2 * stride) ;  
%                 obj.segData(i).curvature = [u,v,k];
%                 % image plot
% %                 y_ = polyval(p,x_);
% %                 x_y = [x_, y_];
% %                 x_y(x_y(:,2)>(height-1),:)=[];
% %                 x_y(x_y(:,2)<0,:)=[];
% %                 img = zeros(height, width);
% %                 index = int64(floor(x_y(:,2)) * width + floor(x_y(:,1)) + 1);
% %                 img(index) = 255;
%            end    
%        end
%    function obj = checkCurvature(varargin)
%            obj = varargin{1};
%            color = '.r';
%            figure;hold on;
%            numSample = size(obj.segData,2);
%            x_event = [];
%            y_event = [];
%            for i = 1: numSample
%                x = obj.segData(i).curvature(:,2);
%                y = obj.segData(i).curvature(:,1);
%                z = max(obj.segData(i).curvature(:,3));  
%                if strcmp(color,'.r')
%                    color = '.b';
%                else
%                    color = '.r';
%                end
%                
%                if z > 0.03
%                    x_event = [x_event;x];
%                    y_event = [y_event;y];
%                end
%                plot(x,y,color);
%            end  
%             plot(x_event,y_event,'.g');
%            
%    end
%    function obj = eventAnalysis(varargin)
%           obj = varargin{1};
%           eventField = obj.eventField;
%           tLabel = obj.eventLabel(:,eventField);
%           event_start = obj.timeParser(obj.eventLabel.StartTime) + ...
%               obj.timeDelayforVideo;
%           event_end =  obj.timeParser(obj.eventLabel.EndTime) + ...
%               obj.timeDelayforVideo;
%           logRate = 1 / double(obj.logData.time(2));
%           event_start_idx = round(event_start * logRate)+1;
%           event_end_idx =  round(event_end * logRate)+1;
%           event_time = [event_start, event_end];
%           eventNum = length(event_start);
%           eventIndex = 0;
%           obj.eventAna = table();
%           for i = 1 : eventNum
%                label_index = table2array(tLabel(i,:));
%                label_index = find(label_index==1);
%                if ~isempty(label_index)
%                obj.eventAna.name(eventIndex + 1,1) = ...
%                    obj.eventField(label_index);               
%                obj.eventAna.duration(eventIndex + 1,1) = event_end(i) - ...
%                    event_start(i);
%                obj.eventAna.distance(eventIndex + 1,1) = ...
%                obj.logData.distance(event_end_idx(i)) - ...
%                obj.logData.distance(event_start_idx(i));
%                eventIndex = eventIndex + 1;
%                end
%           end          
%       end
%%
end

end