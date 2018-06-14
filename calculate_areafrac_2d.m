function areafrac = calculate_areafrac_2d(node,link)
res = 0.037;
nodet = struct2table(node);
linkt = struct2table(link);
nodemat = [(1:height(nodet))', res*nodet.comx, res*nodet.comy,zeros(size(nodet.comx))];
% x (2nd column) direction is z direction??
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

% calculate area fraction
vectors = horzcat(nodemat(linkmat(:,3),2) - nodemat(linkmat(:,2),2),...
    nodemat(linkmat(:,3),3) - nodemat(linkmat(:,2),3),...
    nodemat(linkmat(:,3),4) - nodemat(linkmat(:,2),4));
lengths = sqrt(vectors(:,1).^2 + vectors(:,2).^2 + vectors(:,3).^2);
areas = lengths.*linkmat(:,4)/2;
tot_area = (max(nodemat(:,2)) - min(nodemat(:,2)))*(max(nodemat(:,3)) - min(nodemat(:,3)));
areafrac = sum(areas)/tot_area;
