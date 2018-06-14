% parse unsorted txt file of element statuses for each frame, horizontally
% concatenated
% 1st column is element #
% 2nd column is element status
clc;
df =9;%h2=figure();
load(['~/Documents/bone-networks/2Dnets/newimg1-MLd' num2str(df) '.mat'])
% [node,link]=generate_2D_slice_mirr(node,link);
numElem = length(link);
allstats = zeros(numElem,34,1);%17076,17,1);%16262 for d1
for i = 1:size(allstats,3)
%     path = ['~/Dropbox/abaqus/slicetopologyOpt_001_' num2str(i) '_status.txt'];
    path = (['~/Dropbox/abaqus/newimg1-MLd' num2str(df) '-mises_status.txt']);
    statuses = csvread(path,0,0,[0 0 size(allstats,1)*size(allstats,2)-1 1]);
    first1s = find(statuses(:,1)==1,2,'first');
    numElem = first1s(2) - first1s(1);
    numFrames = size(statuses,1)/numElem;
    if mod(numFrames,1)~=0
        disp('something is wonky');
        numFrames = floor(numFrames);
        statuses(numFrames*numElem:end,:) = [];
    end

    % reshape statuses so that each frame is a slice of a 3D array
    statuses = permute(reshape(statuses',[2,numElem,numFrames]),[2,1,3]);

    % sort each slice 
    [m,n,p] = size(statuses);
    [~, row_ind] = sort(statuses(:,1,:), 1);
    lin_ind = bsxfun(@plus, bsxfun(@plus, row_ind, (0:n-1)*m), reshape((0:p-1)*m*n, 1, 1, p));
    statuses = statuses(lin_ind);

    % reshape into matrix containing #elem rows and #frames columns
    statuses(:,1,:) = [];
    statuses = squeeze(statuses);
    allstats(:,:,i) = statuses;
%     figure();
%     hold on
%     for i = 1:size(statuses,1)
%         plot(statuses(i,:))
%     end
end
% %%
% laststat = squeeze(allstats(:,end,:));
%     figure();
% hold on;
% for i = 1:size(laststat,2)
%     plot(1:length(laststat),(laststat(:,i)==0))
% end
% %%
% [node,link]=generate_2D_slice(node,link);
centroids = zeros(length(link),2);
for i = 1:length(link)
    centroids(i,1) = (node(link(i).n1).comx + node(link(i).n2).comx)/2;
    centroids(i,2) = (node(link(i).n1).comy + node(link(i).n2).comy)/2;
end 
centroids_x = num2cell(centroids(:,1));
centroids_y = num2cell(centroids(:,2));
[link.centroids_x] = centroids_x{:};
[link.centroids_y] = centroids_y{:};

figure()
% scatter3(centroids(:,1),centroids(:,2),zeros(size(centroids(:,2))),'markeredgecolor','g');
t = 30;
crackpts1 = [centroids(find(allstats(:,t)==0),1),centroids(find(allstats(:,t)==0),2)];

scatter3(crackpts1(:,1),crackpts1(:,2),t*ones(size(crackpts1(:,2))),'markeredgecolor',[ 0    0.4470    0.7410],'markerfacecolor',[ 0    0.4470    0.7410]);
hold on;
t = t+1;
crackpts2 = [centroids(find(allstats(:,t)==0),1),centroids(find(allstats(:,t)==0),2)];
scatter3(crackpts2(:,1),crackpts2(:,2),t*ones(size(crackpts2(:,2))),'markeredgecolor',[ 0.8500    0.3250    0.0980],'markerfacecolor',[ 0.8500    0.3250    0.0980]);
commonpts = intersect(crackpts1,crackpts2,'rows');
for i = 1:length(commonpts)
    plot3([commonpts(i,1);commonpts(i,1)],[commonpts(i,2);commonpts(i,2)],[t-1;t],'color','k')
end
t = t+1;
crackpts3 = [centroids(find(allstats(:,t)==0),1),centroids(find(allstats(:,t)==0),2)];
scatter3(crackpts3(:,1),crackpts3(:,2),t*ones(size(crackpts3(:,2))),'markeredgecolor',[0.4660    0.6740    0.1880],'markerfacecolor',[0.4660    0.6740    0.1880]);
commonpts = intersect(crackpts2,crackpts3,'rows');
for i = 1:length(commonpts)
    plot3([commonpts(i,1);commonpts(i,1)],[commonpts(i,2);commonpts(i,2)],[t-1;t],'color','k')
end
t = t+1;
crackpts4 = [centroids(find(allstats(:,t)==0),1),centroids(find(allstats(:,t)==0),2)];
scatter3(crackpts4(:,1),crackpts4(:,2),t*ones(size(crackpts4(:,2))),'markeredgecolor',[ 0    0.4470    0.7410],'markerfacecolor',[ 0    0.4470    0.7410]);
commonpts = intersect(crackpts3,crackpts4,'rows');
for i = 1:length(commonpts)
    plot3([commonpts(i,1);commonpts(i,1)],[commonpts(i,2);commonpts(i,2)],[t-1;t],'color','k')
end
t = t+1;
crackpts5 = [centroids(find(allstats(:,t)==0),1),centroids(find(allstats(:,t)==0),2)];
scatter3(crackpts5(:,1),crackpts5(:,2),t*ones(size(crackpts5(:,2))),'markeredgecolor',[ 0.8500    0.3250    0.0980],'markerfacecolor',[ 0.8500    0.3250    0.0980]);
commonpts = intersect(crackpts4,crackpts5,'rows');
for i = 1:length(commonpts)
    plot3([commonpts(i,1);commonpts(i,1)],[commonpts(i,2);commonpts(i,2)],[t-1;t],'color','k')
end
% edges = 200:200:1400;
% [N,~,bin]=histcounts(crackpts(:,2),edges);
% ubin = unique(bin);
% crackpath = zeros(length(ubin),2);
% count=1;
% for i = 1:length(ubin)
%     crackpath(i,:) = mean(crackpts(count:(count+N(ubin(i))-1),:));
%     count=count+N(ubin(i));
%     
% end
% plot(crackpath(:,1),crackpath(:,2),'linewidth',3,'color','k')
title(['dilation factor = ' num2str(df)])
% title(['t = ' num2str(t) ', dilation factor = ' num2str(df)])
%%

% h = figure;
% axis tight manual % this ensures that getframe() returns a consistent size
filename = 'test_t13.gif';
      % Capture the plot as an image 
frame = getframe(h2); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256); 
% Write to the GIF File 
imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 

frame = getframe(h3); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256);  
imwrite(imind,cm,filename,'gif','WriteMode','append'); 

frame = getframe(h4); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256);  
imwrite(imind,cm,filename,'gif','WriteMode','append'); 

frame = getframe(h5); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256);  
imwrite(imind,cm,filename,'gif','WriteMode','append'); 

frame = getframe(h6); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256);  
imwrite(imind,cm,filename,'gif','WriteMode','append'); 

frame = getframe(h7); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256);  
imwrite(imind,cm,filename,'gif','WriteMode','append'); 

frame = getframe(h8); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256);  
imwrite(imind,cm,filename,'gif','WriteMode','append'); 

frame = getframe(h9); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256);  
imwrite(imind,cm,filename,'gif','WriteMode','append'); 

frame = getframe(h10); 
im = frame2im(frame); 
[imind,cm] = rgb2ind(im,256);  
imwrite(imind,cm,filename,'gif','WriteMode','append'); 

      
%% Determining crack width
crx = centroids(find(allstats(:,t)==0),1);
cry = centroids(find(allstats(:,t)==0),2);

[coeff,score] = pca([crx cry]);
figure()
scatter(score(:,1),score(:,2))
asp_ratio = (max(score(:,1)) - min(score(:,1)))/std(score(:,2))
width=max(score(:,1)) - min(score(:,1))
xlabel('PC1')
ylabel('PC2')
title(['t = ' num2str(t) ', dilation = ' num2str(df) ', asp ratio = ' num2str(asp_ratio) ', width = ' num2str(width)])
%%
load(['~/Documents/bone-networks/2Dnets/newimg1-MLd2.mat'])
numElem = 26772;%27082;
allstats = zeros(numElem,9,1);%17076,17,1);%16262 for d1
for i = 1:size(allstats,3)
%     path = ['~/Dropbox/abaqus/slicetopologyOpt_001_' num2str(i) '_status.txt'];
    path = (['~/Dropbox/abaqus/newimg1-MLd2-slice-mirr-pbc-mises_status.txt']);
    statuses = csvread(path,0,0,[0 0 size(allstats,1)*size(allstats,2)-1 1]);
    first1s = find(statuses(:,1)==1,2,'first');
    numElem = first1s(2) - first1s(1);
    numFrames = size(statuses,1)/numElem;
    if mod(numFrames,1)~=0
        disp('something is wonky');
        numFrames = floor(numFrames);
        statuses(numFrames*numElem:end,:) = [];
    end

    % reshape statuses so that each frame is a slice of a 3D array
    statuses = permute(reshape(statuses',[2,numElem,numFrames]),[2,1,3]);

    % sort each slice 
    [m,n,p] = size(statuses);
    [~, row_ind] = sort(statuses(:,1,:), 1);
    lin_ind = bsxfun(@plus, bsxfun(@plus, row_ind, (0:n-1)*m), reshape((0:p-1)*m*n, 1, 1, p));
    statuses = statuses(lin_ind);

    % reshape into matrix containing #elem rows and #frames columns
    statuses(:,1,:) = [];
    statuses = squeeze(statuses);
    allstats(:,:,i) = statuses;
%     figure();
%     hold on
%     for i = 1:size(statuses,1)
%         plot(statuses(i,:))
%     end
end
%%
% laststat = squeeze(allstats(:,end,:));
%     figure();
% hold on;
% for i = 1:size(laststat,2)
%     plot(1:length(laststat),(laststat(:,i)==0))
% end
% %%
[node,link]=generate_2D_slice_mirr(node,link);
centroids = zeros(length(link),2);
for i = 1:length(link)
    centroids(i,1) = (node(link(i).n1).comx + node(link(i).n2).comx)/2;
    centroids(i,2) = (node(link(i).n1).comy + node(link(i).n2).comy)/2;
end 
centroids_x = num2cell(centroids(:,1));
centroids_y = num2cell(centroids(:,2));
[link.centroids_x] = centroids_x{:};
[link.centroids_y] = centroids_y{:};
%%
figure();
scatter(centroids(:,1),centroids(:,2),'markeredgecolor','g');
hold on
t = 3;
scatter(centroids(find(allstats(:,t)==0),1),centroids(find(allstats(:,t)==0),2));
title(['t = ' num2str(t)])