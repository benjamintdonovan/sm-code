function [ N ] = MakeTracesDual( fname )

% this function reads a Dual-View TIF file and a .txt file that is output from
% 'find maxima' in ImageJ. These two input file should have the same
% filename. 
% The input filename 'fname' should not contain any extensions
% It then get the average intensities at each maxima position, save the
% donor and acceptor traces into struct array 'traces'
% 
% *******NOTE************************************************************
% in ImageJ, coordinate the upper left corner of an image is (0, 0)
% in matlab, it's (1, 1)

% ********888************************************************************
Q16 = [1/16,2/16,1/16;2/16,4/16,2/16;1/16,2/16,1/16];    % convolution factor
dlmname = [fname,'.txt'];
imgname = [fname,'.TIF'];
sname = [fname,'.mat'];

traces = struct('donr',{},'acptr',{},'IndT',0,'imName',{},'position',{});

P = dlmread(dlmname,'',1,0);  % positions, P(:,2) x coordinates, P(:,3) y coordinates

img_info = imfinfo(imgname);

N = size(P,1);                        % total number of 'possible' fluorophores (intensity maxima)

img_num = size(img_info,1);              % number of frames, or the length of the trace

for k = 1:N
traces(k).donr = zeros(img_num,1); 
traces(k).acptr = zeros(img_num,1); 
traces(k).IndT = 0;
traces(k).position = P(k,2:3);
traces(k).imName = imgname;
end

tic
for j = 1:img_num
    
%    tic
    img = imread(imgname, j, 'Info', img_info);
    
    for i = 1:N
        
        acptrM = double(img(P(i,3):P(i,3)+2, P(i,2):P(i,2)+2));               % coordinate index different in ImageJ and matlab!
        donrM = double(img(P(i,3):P(i,3)+2, P(i,2)+128:P(i,2)+130));
        traces(i).donr(j) = sum(sum(donrM.*Q16));
        traces(i).acptr(j) = sum(sum(acptrM.*Q16));
    end
%     toc
    
end
toc
save(sname, 'traces');

end

