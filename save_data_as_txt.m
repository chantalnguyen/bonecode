% saves node and link data as text files
res = 0.037;

for df = 2:14
    load(['~/Documents/bone-networks/2Dnets/309555-e' num2str(df) '.mat'])

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
    
    csvwrite(['309555-e' num2str(df) '-node.txt'],nodemat);
    csvwrite(['309555-e' num2str(df) '-link.txt'],linkmat);
end