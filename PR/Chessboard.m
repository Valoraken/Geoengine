function [res] = Chessboard(gt,fv,r,c)
%<Chessboard> function applies image segmentation with Chessboard algorithm
%   Detailed explanation goes here
%   INPUT
%       gt: groud truth
%
%       r: number of row in segment
%
%       c: number of column in segment
%
%       fv: feature vector
%
%
%   OUTPUT
%       res: segmentation result, contains 7 elements 
%
%           1.mask:      for each segment the correspounding value
%                        is the majorty vote of the ground truth
%           2.list_fv:   each row is the feature vector of a segment
%           3.list_label:each row is the label of a segment 
%           4.img_fv:    image with its segments filled with feature vector 
%           5.img_label: image with its segments filled with labels
%           6.accuracy:  overall segmentation accuracy 
%           7.num_seg:   number of segmentations

tic;
[mr,mc] = size(gt);

% number of segments
num_seg = ceil(mr/r) * ceil(mc/c);
segc = c;
segr = r;


% initialize the output
mask = zeros(mr,mc,'like',gt);
list_fv = zeros(num_seg, size(fv,3));
list_label = zeros(num_seg,1,'like',gt);
img_fv = zeros(size(fv));
img_label = zeros(mr,mc,'like',gt);
index = 1;

% calculate for each segment
for i=1:r:mr
    for j=1:c:mc
        
       % calculate index for current grid cur_r cur_c
       if i+r<= mr
           cur_r = [i:1:i+r-1];
       else
           cur_r = [i:1:mr];
           segr = mr-i+1;
           
       end
       
       if j+c<= mc
           cur_c = [j:1:j+c-1];
       else
           cur_c = [j:1:mc];
           segc = mc-j+1;
       end 
       % find the majority vote
       intM = mode(gt(cur_r,cur_c),1);
       M = mode(intM,2);
       
       % claculate the output
       mask(cur_r,cur_c) = M; 
       c_fv = sum(sum(fv(cur_r,cur_c,:)))/(segr*segc);
       list_fv(index,:) = c_fv;
       list_label(index) = M;
       img_fv(cur_r,cur_c,:)=ones(size(img_fv(cur_r,cur_c,:))).*c_fv;
       index = index + 1;
    
    end
end 
    
img_label = mask;

% calculate overall accuracy
accuracy = length( find((img_label - gt)==0)) / (mr*mc);

res = struct('mask',mask,'list_fv',list_fv,'list_label',list_label,...
             'img_fv', img_fv,'img_label',img_label,'accuracy',accuracy,...
             'num_segments',num_seg);
runtime=toc;
end