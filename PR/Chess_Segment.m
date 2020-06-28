function [Output] = Chess_Segment(gt,fv,r,c)
%	This function is to apply Chessboard Segmentation to images accroding
%	to the 
%   Tianqi Xiao
%   3371477
tic;
%% Generate Mask
% number of clusters
k=size(fv,1)/r*size(fv,2)/c;
Kernel=ones(r,c);
% Segment mask
mask =cell2mat(reshape(arrayfun(@(p) Kernel*p, 1:k, 'UniformOutput', false),size(fv,1)/r,size(fv,2)/c)');
% image feature vector
Img_fv=cat(3,fv,mask);
Seg_fv=blockproc(Img_fv, [r,c],@(p) mean(mean(p.data)));
% respective feature
Img_res=repelem(Seg_fv,r,c);
Seg_label=blockproc(gt, [r,c],@(p) mode(mode(p.data)));
Img_label=repelem(Seg_label,r,c);
% sort features and label from mat to list
List_fv=sortrows(reshape(Seg_fv,k,[]),6);
List_label=reshape(Seg_label',k,[]);
% outputImage = zeros(size(fv(:,:,1:3)),'like',fv);
% idx = label2idx(mask);
% numRows = size(fv,1);
% numCols = size(fv,2);
% for labelVal = 1:k
%     [m_r,m_c]=find(mask==labelVal);
%     segfv(labelVal,:)=mean(mean(Img_fv(m_r,m_c,:)));
% %     redIdx = idx{labelVal};
% %     greenIdx = idx{labelVal}+numRows*numCols;
% %     blueIdx = idx{labelVal}+2*numRows*numCols;
% %     outputImage(redIdx) = mean(fv(redIdx));
% %     outputImage(greenIdx) = mean(fv(greenIdx));
% %     outputImage(blueIdx) = mean(fv(blueIdx));
%    Img_rgb(m_r,m_c,1)=segfv(1);
% end  
% get boundaries of segment
Seg_bw = boundarymask(mask);
% imshow(imoverlay(fv(:,:,1:3),Seg_bw,'cyan'),'InitialMagnification',67)
accuracy = length(find((Img_label - gt)==0)) / (size(gt,1)*size(gt,2));

Output = struct('Mask',mask,'Feature_List',List_fv,'Label_List',List_label,'Seg_RGB', Img_res,'Seg_Label',Img_label,'Accuracy',accuracy,'Seg_k',k,'Seg_Border',Seg_bw);
runtime=toc;        


end

