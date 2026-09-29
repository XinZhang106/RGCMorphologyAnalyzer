function [nodes,edges,radii,nodeTypes,abort,soma_size] = readArborTrace(fileName,validNodeTypes)
abort = false; 
nodes = []; 
edges = []; 
nodeTypes = [];
soma_size = 0;
validNodeTypes = setdiff(validNodeTypes,1); % 1 is for soma %%setdiff(A,B) returns data in a that is not in b
% read the SWC file
[nodeID, nodeType, xPos, yPos, zPos, radii, parentNodeID] = textread(fileName, '%u%d%f%f%f%f%d','commentstyle', 'shell');
% every tree should start from a node of type 1 (soma)
nodeType(find(parentNodeID==-1))=1;

%calculate the size of the soma, doesn't need flattening
soma_end_idx = find(diff(parentNodeID) < 0, 1);
%somapixN = sum(find(parentNodeID==-1));
soma_path = zeros([soma_end_idx,3]);
soma_path(:, 1) = xPos(1:soma_end_idx);
soma_path(:, 2) = yPos(1:soma_end_idx);
soma_path(:, 3) = zPos(1:soma_end_idx);
%soma_path = [x(:) y(:)];
s = svd(soma_path - mean(soma_path,1));
if numel(s) < 2 || s(2) < 1e-6*s(1)          % collinear
    somavec = soma_path(end, 1) - soma_path(1, :);
    soma_size = norm(somavec);

else                                        % minimum enclosing circle (Welzl)
    H = soma_path(convhull(soma_path),:);  H = H(randperm(end),:);
    c = H(1,:);  r = 0;
    for i = 2:size(H,1), if norm(H(i,:)-c) > r*(1+1e-12)
        c = H(i,:);  r = 0;
        for j = 1:i-1, if norm(H(j,:)-c) > r*(1+1e-12)
            c = (H(i,:)+H(j,:))/2;  r = norm(H(i,:)-c);
            for k = 1:j-1, if norm(H(k,:)-c) > r*(1+1e-12)
                c = ((2*[H(j,:)-H(i,:); H(k,:)-H(i,:)]) \ ...
                     [sum(H(j,:).^2)-sum(H(i,:).^2); sum(H(k,:).^2)-sum(H(i,:).^2)])';
                r = norm(H(i,:)-c);
            end, end
        end, end
    end, end
    soma_size = 2*r;                                % diameter; c is the centre
end


% find the first soma node in the list (more than one node can be labeled as soma)
firstSomaNode = find(nodeType == 1 & parentNodeID == -1, 1);
% find the average position of all the soma nodes, and assign it as THE soma node position
somaNodes = find(nodeType == 1);
somaX = mean(xPos(somaNodes)); 
somaY = mean(yPos(somaNodes));
somaZ = mean(zPos(somaNodes));
somaRadius = mean(radii(somaNodes));
xPos(firstSomaNode) = somaX; yPos(firstSomaNode) = somaY; zPos(firstSomaNode) = somaZ;
radii(firstSomaNode) = somaRadius;
% change parenthood so that there is a single soma parent
parentNodeID(ismember(parentNodeID,somaNodes)) = firstSomaNode;
% delete all the soma nodes except for the firstSomaNode
nodesToDelete = setdiff(somaNodes,firstSomaNode);
nodeID(nodesToDelete)=[]; nodeType(nodesToDelete)=[];
xPos(nodesToDelete)=[]; yPos(nodesToDelete)=[]; zPos(nodesToDelete)=[];
radii(nodesToDelete)=[]; parentNodeID(nodesToDelete)=[];
% reassign node IDs due to deletions
for kk = 1:numel(nodeID)
  while ~any(nodeID==kk) %logical 1 is true
    nodeID(nodeID>kk) = nodeID(nodeID>kk)-1;
    parentNodeID(parentNodeID>kk) = parentNodeID(parentNodeID>kk)-1;
  end
end
% of all the nodes, retain the ones indicated in validNodeTypes
% ensure connectedness of the tree if a child is marked as valid but not some of its ancestors
validNodes = nodeID(ismember(nodeType,validNodeTypes));% outputs 1 (true) when a is in b (a,b)
additionalValidNodes = [];
for kk = 1:numel(validNodes)
  thisParentNodeID = parentNodeID(validNodes(kk)); thisParentNodeType = nodeType(thisParentNodeID);
  while ~ismember(thisParentNodeType,validNodeTypes)
    if thisParentNodeType == 1
      break;
    end
    additionalValidNodes = union(additionalValidNodes, thisParentNodeID); nodeType(thisParentNodeID) = validNodeTypes(1);
    thisParentNodeID = parentNodeID(thisParentNodeID); thisParentNodeType = nodeType(thisParentNodeID);
  end
end
% retain the valid nodes only
% the soma node is always a valid node to ensure connectedness of the tree
validNodes = [firstSomaNode; validNodes; additionalValidNodes']; validNodes = unique(validNodes);% unique function returns same data without repeats
nodeID = nodeID(validNodes); nodeType = nodeType(validNodes); parentNodeID = parentNodeID(validNodes);
xPos = xPos(validNodes); yPos = yPos(validNodes); zPos = zPos(validNodes); radii = radii(validNodes);
% reassign node IDs after deletions
for kk = 1:numel(nodeID)
  while ~any(nodeID==kk)
    nodeID(nodeID>kk) = nodeID(nodeID>kk)-1;
    parentNodeID(parentNodeID>kk) = parentNodeID(parentNodeID>kk)-1;
  end
end
% return the resulting tree data
nodes = [xPos yPos zPos];
edges = [nodeID parentNodeID];
edges(any(edges==-1,2),:) = [];
nodeTypes = unique(nodeType)';
