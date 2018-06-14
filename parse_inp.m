% parse inp file to obtain nodemat and linkmat
df = 2;
fid = fopen(['~/Dropbox/abaqus/309555_mod_links_e' num2str(df) '.inp']);
data = textscan(fid,'%s');
data = data{1};
nodeind1 = find(contains(data,'*NODE'),1,'first')+1;
nodeind2 = find(contains(data,'*Element'),1,'first')-1;
nodemat = data(nodeind1:nodeind2);
nodemat = cellfun(@str2num,nodemat,'UniformOutput',false);
nodemat = cell2mat(nodemat);
%%
linkind1 = nodeind2+3;
linkind2 = find(contains(data,'*Elset'),1,'first')-1;
linkmat = data(linkind1:linkind2);
linkmat = cellfun(@str2num,linkmat,'UniformOutput',false);
linkmat = cell2mat(linkmat);
%%
linkmat = [linkmat, zeros(length(linkmat),1)];
startind = linkind2+1;
for i = 1:length(linkmat)
    linkmat(i,4) = str2double(data{startind+18*i+2*(i-1)});
end

csvwrite(['309555-e' num2str(df) '-modlinks-node.txt'],nodemat);
csvwrite(['309555-e' num2str(df) '-modlinks-link.txt'],linkmat);