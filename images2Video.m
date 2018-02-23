function  images2Video(dirName,videoName)
list = dir(dirName);
[r,c] = size(list);
vWriter = VideoWriter(videoName, 'MPEG-4');
open(vWriter)

for i = 3 : r
    img = list(i).name;
    img = strcat(dirName,img);
    img = imread(img);
    img = imcrop(img,[300,50,2504,1363]);
    writeVideo(vWriter,img);
    
    
    
end
close(vWriter)
end

