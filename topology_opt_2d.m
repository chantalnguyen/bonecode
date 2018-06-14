% increase size of image by 4x, then erode (or dilate) and resize to
% original dimensions
immat = imread('~/Dropbox/2d-opt/309555.bmp');
immat4 = imresize(immat,4);
for i = 2:14
    se = strel('square',i);
    imerd = imerode(immat4,se);
    immat4_1 = imresize(imerd,0.25);
    immat4_1(immat4_1 >= 130) = 255;
    immat4_1(immat4_1 < 130) = 0;
    imwrite(immat4_1,['~/Dropbox/2d-opt/309555-e' num2str(i) '.bmp']);
end
for i = 2:14
    se = strel('square',i);
    imdil = imdilate(immat4,se);
    immat4_1 = imresize(imdil,0.25);
    immat4_1(immat4_1 >= 130) = 255;
    immat4_1(immat4_1 < 130) = 0;
    imwrite(immat4_1,['~/Dropbox/2d-opt/309555-d' num2str(i) '.bmp']);
end

% have to manually compute trabecular thicknesses in bonej and save results

tbresults = csvread('~/Documents/bone-networks/2Dnets/309555Results.csv',0,1);
slice=imread(['~/Dropbox/2d-opt/309555.bmp']);
skel = Skeleton3D((slice./255));
[node,link]=convertSkelToGraph(skel,2);
thickmap = imread(['~/Documents/bone-networks/2Dnets/309555_Tb.bmp']);
link = getFIJITbTh(thickmap,link,tbresults(tbresults(:,1)==0,4));
save(['~/Documents/bone-networks/2Dnets/309555.mat'],'node','link','skel');
generate_inp_from_2Dnet_topopt_failure(node,link,['309555']);
%%
% afrac = zeros(14,1);
for i = 2:14
%     slice=imread(['~/Dropbox/2d-opt/309555-e' num2str(i) '.bmp']);
%     skel = Skeleton3D((slice./255));
%     [node,link]=convertSkelToGraph(skel,2);
%     thickmap = imread(['~/Documents/bone-networks/2Dnets/309555-e' num2str(i) '_Tb.bmp']);
%     link = getFIJITbTh(thickmap,link,tbresults(tbresults(:,1)==i,4));
%     save(['~/Documents/bone-networks/2Dnets/309555-e' num2str(i) '.mat'],'node','link','skel');
    load(['~/Documents/bone-networks/2Dnets/309555-e' num2str(i) '.mat'])
%     afrac(i-1) = calculate_areafrac_2d(node,link);
%     generate_inp_from_2Dnet_topopt_failure(node,link,['309555-d' num2str(i)]);
    generate_inp_from_2Dnet_topopt_tension(node,link,['309555-e' num2str(i)]);

    
end
%%
for i = 1:13  
    modify_2Dnet_topopt_failure(node,link,'309555',afrac(i),['e' num2str(i+1)]);
end
%%
for i = 1:13  
    modify_2Dnet_topopt_failure(node,link,'309555',0.55,['e3-' num2str(i+1)]);
end

