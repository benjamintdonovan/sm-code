function [ Ind ] = MergeStack( fname, num )

% This function reads the input TIF image stack 'fname', merge the first
% 'num' number of frames into a new single frame TIF file called 'fname_m.TIF'
% fname must be the full filename of a stack tiff image

sname = [strtok(fname,'.'),'_m.TIF'];

info = imfinfo(fname);
img_height = info(1,1).Height;
img_width = info(1,1).Width;

% pre-allocate space for image stack
clear('imgstk', 'imgmrg');
imgstk = zeros(img_height, img_width, num, 'uint16');
imgmrg = zeros(img_height, round(img_width/2), 'uint16');

% load first few frames of the tiff stack, number of frames indicated by
% num
for k = 1:num
    imgstk(:,:,k) = imread(fname, k, 'Info', info);
end

% create the merged image from the average of the image stack
for i = 1:img_height
    for j = 1:round(img_width/2)
        imgmrg(i,j) = mean(imgstk(i,j,:));
    end
end

% plot the merged image
% figure
% image(imgmrg,'CDataMapping','scaled')
% colormap(gray)

imwrite(imgmrg, sname, 'tif');
Ind = j;

end

