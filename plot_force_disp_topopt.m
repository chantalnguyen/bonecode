% plot force-displacement curves

figure()
df = 2;
load(['~/Documents/bone-networks/2Dnets/309555-e' num2str(df) '.mat'])

res = 0.037;
nodet = struct2table(node);
linkt = struct2table(link);
nodemat = [(1:height(nodet))', res*nodet.comx, res*nodet.comy,zeros(size(nodet.comx))];
nodemat = sortrows(nodemat,2);
nodemat = [(1:length(nodemat))',nodemat];

% sort tables and remap nodes
linkmat = [(1:height(linkt))',linkt.n1,linkt.n2];
linkmat = double(linkmat);
linkmat = [linkmat,res*linkt.avgthickness];

map = containers.Map(nodemat(:,2)',nodemat(:,1)');
startNodes = arrayfun(@(x) map(x), linkmat(:,2));
endNodes = arrayfun(@(x) map(x), linkmat(:,3));
linkmat(:,2:3) = [startNodes endNodes];
nodemat(:,2) = [];

topNodes = find(nodemat(:,2) >= (max(nodemat(:,2))-0.5));

path = ['~/Dropbox/abaqus/309555-e' num2str(df) '_force.txt'];
data = csvread(path);
numSteps = 81;
numNode = length(data)/numSteps;
forces = zeros(numSteps,numNode);
for i = 1:numNode
    forces(:,i) = data(((i-1)*numSteps+1):(i*numSteps),3)';
end
forces = -1*forces;

path = ['~/Dropbox/abaqus/309555-e' num2str(df) '_disp.txt'];
data = csvread(path);
disps = zeros(numSteps,numNode);
for i = 1:numNode
    disps(:,i) = data(((i-1)*numSteps+1):(i*numSteps),3)';
end
disps = -1*disps;

filename = ['forcedisp-309e' num2str(df) '.txt'];
dlmwrite(filename,[mean(disps(:,topNodes),3),sum(forces(:,topNodes),3)]);

plot(mean(disps(:,topNodes),3),sum(forces(:,topNodes),3),'color',colors(count,:))
hold on

%% for modlinks:
figure()

df=2;

filename = ['forcedisp-modlinks-309e' num2str(df) '.txt'];
data = dlmread(filename);
plot(data(:,1),data(:,2));
hold on

for k = 2:9
    
    path = ['~/Dropbox/abaqus/309555-modlinks-e'  num2str(df) '-' num2str(k) '_force.txt'];
    data = csvread(path);
    numSteps = 81;
    numNode = length(data)/numSteps;
    forces = zeros(numSteps,numNode);
    for i = 1:numNode
        forces(:,i) = data(((i-1)*numSteps+1):(i*numSteps),3)';
    end
    forces = -1*forces;

    path = ['~/Dropbox/abaqus/309555-modlinks-e'  num2str(df) '-' num2str(k) '_disp.txt'];
    data = csvread(path);
    disps = zeros(numSteps,numNode);
    for i = 1:numNode
        disps(:,i) = data(((i-1)*numSteps+1):(i*numSteps),3)';
    end
    disps = -1*disps;
    
    fid = fopen(['~/Dropbox/abaqus/309555_mod_links_e' num2str(df) '-' num2str(k) '.inp']);
    inpfile = textscan(fid,'%s');
    inpfile = inpfile{1};
    idx = find(contains(inpfile,'nset=topnodes'),1,'first')+2;
    topNode = str2double(inpfile{idx});
    actualNumNode = str2double(inpfile{idx+1});
    fclose(fid);
    topData = topNode;
    if actualNumNode ~= numNode
        count=1;
        topData = [];
        while isempty(topData)
            topData = (find(data(:,1)==topNode,1,'first'));
            count=count+1;
        end
        topData = (topData-1)/numSteps + 1;
    end
    
    filename = ['forcedisp-modlinks-309e' num2str(df) '-' num2str(k) '.txt'];
    dlmwrite(filename,[mean(disps(:,topData:end),2),sum(forces(:,topData:end),2)]);
    plot(mean(disps(:,topData:end),2),sum(forces(:,topData:end),2));%,'color',colors(count,:))
    hold on

end

xlabel('Displacement (mm)','fontsize',18)
ylabel('Force (N)','fontsize',18)
set(gca,'fontsize',18)
