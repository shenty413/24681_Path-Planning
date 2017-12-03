clear;
clc;

step = 0.5;
[x,y] = meshgrid(1:step:15,1:step:15);
tri = delaunay(x,y);
e = floor(14/step+1);
z = peaks(e);
%     trimesh(tri,x,y,z) % plot the mesh data

numV = e^2;
vdata = [reshape(x,numV,1) reshape(y,numV,1) reshape(z,numV,1)];
fdata = tri;

%%
fid=fopen('meshdata.obj','w');
for i = 1:size(vdata,1)
    fprintf(fid,'v %d %d %d\n',vdata(i,1),vdata(i,2),vdata(i,3));
end

for i = 1:size(fdata,1)
    fprintf(fid,'f %d %d %d\n',fdata(i,1),fdata(i,2),fdata(i,3));
end