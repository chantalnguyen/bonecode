% Modifies existing inp file by randomly removing links to obtain specified
% area fraction
function modify_2Dnet_topopt_failure(node,link,filename,afrac,outname,radMult,res)
% inputs:
% node: struct containing node information
% link: struct containing link information
% filename: name of original input file
% afrac: desired area fraction
% outname: string to append to filename (optional)
% radMult: factor by which to multiply radii (optional, default=1)
% res: resolution in mm/pix (optional, default = 0.037)

if nargin < 6
    radMult = 1;
end
if nargin < 7
    res = 0.037; % resolution (mm)
end

if nargin < 5
    outname = [];
end

nodet = struct2table(node);
linkt = struct2table(link);
nodemat = [(1:height(nodet))', res*nodet.comx, res*nodet.comy,zeros(size(nodet.comx))];
% x (2nd column) direction is direction of loading
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

while areafrac > afrac
    removeind = randi(length(linkmat));
    linkmat(removeind,:) = [];
    lengths(removeind,:) = [];
    areas = lengths.*linkmat(:,4)/2;
    areafrac = sum(areas)/tot_area;
end
linkmat(:,1) = (1:length(linkmat))';

bottomNodes = find(nodemat(:,2) <= (min(nodemat(:,2))+0.5));
topNodes = find(nodemat(:,2) >= (max(nodemat(:,2))-0.5));
%% Create .inp file

fid = fopen(['~/Dropbox/abaqus/' filename '_mod_links_' outname '.inp'],'w');
fprintf(fid,['**Modified input file derived from ' filename '\n' ]);
fprintf(fid,['**with area fraction of ' num2str(areafrac) '\n' ]);
% first add the heading
fprintf(fid,'*Heading\nUnits: millimetres(mm)\n');

%% Part: Bone Part
% print main bone part
fprintf(fid,'**\n** PARTS\n**\n*Part, name=BONE-PART\n');

% then list the nodes using coordinates from the table
fprintf(fid,'*NODE\n');
fprintf(fid,'%i,%f,%f,%f\n',nodemat');

% then list the elements (links) using coordinates from the table
fprintf(fid,'*Element, type=B31\n');
fprintf(fid,'%i,%i,%i\n',linkmat(:,1:3)');

% then define one element set for each link, because each link has a
% different thickness
for k = 1:length(linkmat)
    elsetName = ['link' num2str(k)];
    fprintf(fid,['*Elset, elset=' elsetName ', internal, generate\n '...
        ' ' num2str(k) ', ' num2str(k) '\n']);
    fprintf(fid,['** Section: Section-' num2str(k) ' Profile: Profile-' num2str(k) '\n*Beam Section, elset=' ...
        elsetName ', material=Material-1, temperature=GRADIENTS, POISSON=0.16, section=CIRC\n' ...
        num2str(radMult*linkmat(k,4)/2) '\n0.,0.,-1.\n']);
end

fprintf(fid,'*Nset, nset=allnodes, generate\n');
fprintf(fid,['1, ' num2str(size(nodemat,1)) ', 1\n']);
fprintf(fid,'*Nset, nset=topnodes, generate\n');
fprintf(fid,[num2str(min(topNodes)) ', ' num2str(max(topNodes)) ', 1\n']);
fprintf(fid,'*Nset, nset=bottomnodes, generate\n');
fprintf(fid,[num2str(min(bottomNodes)) ', ' num2str(max(bottomNodes)) ', 1\n']);
fprintf(fid,'*Elset, elset=all-elem, generate\n');
fprintf(fid,['1, ' num2str(size(linkmat,1)) ', 1\n']);
fprintf(fid,'*End Part\n');

%% Assembly
fprintf(fid,['**\n** ASSEMBLY\n**\n*Assembly, name=Assembly\n**\n' ...
    '*Instance, name=BONE-PART-1, part=BONE-PART\n*End Instance\n**\n' ...
    '*End Instance\n**\n']);
fprintf(fid,'**\n*Nset, nset=SET-1, instance=BONE-PART-1, generate\n');
fprintf(fid,[num2str(min(topNodes)) ', ' num2str(max(topNodes)) ', 1\n']);
fprintf(fid,'*End Assembly\n');

%% Amplitude
% print amplitude for displacement
fprintf(fid,'*Amplitude, name=AMP-1\n 0., 0., 1., 1.\n');
fprintf(fid,'*Amplitude, name=smooth-amp, definition=SMOOTH STEP\n 0., 0., 0.0001, 1.\n');
%% Materials
% print materials
fprintf(fid,'**\n** MATERIALS\n**\n');
% material with young's modulus and poisson ratio
fprintf(fid,'*Material, name=MATERIAL-1\n');
fprintf(fid,'*Density\n');
fprintf(fid,'2e-06,\n');
fprintf(fid,'*Depvar, delete=1\n');
fprintf(fid,'1,\n');
fprintf(fid,'*Elastic\n10000., 0.16\n');
fprintf(fid,'*User Defined Field\n');
fprintf(fid,'**\n** INTERACTION PROPERTIES\n**\n*Surface Interaction, name=IntProp-1\n');
%% Boundary Conditions
% print boundary conditions
fprintf(fid,'**\n** BOUNDARY CONDITIONS\n**\n');
fprintf(fid,'** Name:Disp-BC-1 Type: Symmetry/Antisymmetry/Encastre\n');
fprintf(fid,'*Boundary\n');
fprintf(fid,'BONE-PART-1.BOTTOMNODES, ENCASTRE\n');

%% Step
fprintf(fid,'** -----------------------------------------------\n**\n');
fprintf(fid,'** STEP: displace-top\n**\n*Step, name=displace-top, nlgeom=YES\n');
fprintf(fid,'*Dynamic, Explicit\n, 0.02\n*Bulk Viscosity\n0.06,1.2\n');%, stabilize=0.0002, allsdtol=0.05, continue=NO\n');
fprintf(fid,'**\n** BOUNDARY CONDITIONS\n**\n');
fprintf(fid,'** Name:Disp-BC-2 Type: Velocity/Angular velocity\n');
fprintf(fid,'*Boundary, amplitude=smooth-amp, type=VELOCITY\n');
fprintf(fid,'BONE-PART-1.TOPNODES, 1, 1, -0.004\n');
fprintf(fid,'**\n** INTERACTIONS\n**\n');
fprintf(fid,'** Interaction: Int-1\n*Contact, op=NEW\n*Contact Inclusions, ALL EXTERIOR\n');
fprintf(fid,'*Contact Property Assignment\n , , IntProp-1\n*Surface Property Assignment, property=THICKNESS\n');
fprintf(fid,' , ORIGINAL, 1.\n');
fprintf(fid,'**\n** OUTPUT REQUESTS\n**\n');
fprintf(fid,'*Restart, write, number interval=1, time marks=NO\n');
fprintf(fid,'**\n** FIELD OUTPUT; F-Output-1\n**\n');
fprintf(fid,'*Output, field, time interval=0.00025\n');
fprintf(fid,'*Element Output, directions=YES\nE,FV,S,STATUS\n');
fprintf(fid,'**\n*Node Output\nRF, U\n');
fprintf(fid,'**\n** HISTORY OUTPUT: H-Output-1\n**\n');
fprintf(fid,'*Output, history, time interval=0.00025\n');
fprintf(fid,'*Node Output, nset=BONE-PART-1.ALLNODES\n');
fprintf(fid,'RF1, U1\n');
fprintf(fid,'*End Step');
fclose(fid);