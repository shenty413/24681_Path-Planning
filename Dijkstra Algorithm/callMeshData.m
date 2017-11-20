function [vdata,fdata,numV] = callMeshData(step)
% [vdata,fdata,numV] = callMeshData returns a 15*15 meshgrid
%
% Output:
%   vdata: 225*3 vertice data matrix
% 	fdata: 392*3 face data matrix
%   numV: number of vertices
% 

%     [x,y] = meshgrid(1:0.1:15,1:0.1:15);

    [x,y] = meshgrid(1:step:15,1:step:15);
    tri = delaunay(x,y);
    e = floor(14/step+1);
    z = peaks(e);
%     trimesh(tri,x,y,z) % plot the mesh data

    numV = e^2;
    vdata = [reshape(x,numV,1) reshape(y,numV,1) reshape(z,numV,1)];
    fdata = tri;

end