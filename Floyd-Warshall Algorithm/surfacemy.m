clear all;%close all;
% [x,y] = meshgrid(1:0.5:15,1:0.5:15);
%     tri = delaunay(x,y);
% %     z = ones(15,15);
%     z = peaks(29);
% %     trimesh(tri,x,y,z) % plot the mesh data
% numV = 29*29;
%     v = [reshape(x,29*29,1) reshape(y,29*29,1) reshape(z,29*29,1)];
%     f = tri;
[v,f,numV] = callMeshData(0.3);

%extract edges according to faces   
edges = [];
for i = 1:size(f,1)
    for j=1:size(f(i,:),2)
        if j==size(f(i,:),2)
            edges = [edges;[f(i,j),f(i,1)]];
        else
            edges = [edges;[f(i,j),f(i,j+1)]];
        end
    end
end
for i = 1:size(f,1)
    j=size(f(i,:),2);
    while(j>=1)
        if j==1
            edges = [edges;[f(i,1),f(i,size(f(i,:),2))]];
        else
            edges = [edges;[f(i,j),f(i,j-1)]];
        end
        j=j-1;
    end
end
%delete the repeated edges
uedges = unique(edges,'rows');  
length=[];
for i =1:size(uedges,1)
    length=[length;norm(v(uedges(i,1),:)-v(uedges(i,2),:))];
end
uedges = [uedges,length];
for i =1:size(uedges,1)
    if(uedges(i,1)==uedges(i,2))
        uedges(i,:)=[0,0,0];
    end
end
uedges1 = unique(uedges,'rows');
if(sum(uedges1(1,:)==[0,0,0])==3)
uedges1(1,:)=[];
end

%construct the weight/distance matrix
startpoint = 29;
endpoint = 900;
spoints = unique(uedges1(:,1));
currentpoint = startpoint;
dist = ones(size(spoints,1),size(spoints,1))*1000000;
path =zeros(size(spoints,1),size(spoints,1));
for i =1:size(spoints,1)
    index = find(uedges1(:,1)==i);
    connectpoint = uedges1(index,2);
    for j = 1:size(spoints,1)
        if(i==j)
            dist(i,j)=0;
        elseif(ismember(j,connectpoint))
            dist(i,j)=norm(v(i,:)-v(j,:));
        end
          
    end
end
%floyd-warshall process get the full path matrix
fulldist = dist;
for i = 1:size(spoints,1)
    for j =1:size(spoints,1)
        for k =1:size(spoints,1)
            if((fulldist(i,j)>(fulldist(i,k)+fulldist(k,j)))&&fulldist(i,k)<1000000&&fulldist(i,k)<1000000)
                     fulldist(i,j)=fulldist(i,k)+fulldist(k,j);
                     path(i,j)= k;
                     
            end
        end
    end
end
%output the path according to the start and end points

spath = output(startpoint,endpoint,path);
spath = [startpoint,spath];
figure
hold on
for i = 1:size(f,1)
    v1=v(f(i,1),:);v2=v(f(i,2),:);v3=v(f(i,3),:);
    plot3([v1(1) v2(1) v3(1) v1(1)],[v1(2) v2(2) v3(2) v1(2)],[v1(3) v2(3) v3(3) v1(3)],'k-');
end
plot3(v(spath,1),v(spath,2),v(spath,3),'r');
scatter3(v(startpoint,1),v(startpoint,2),v(startpoint,3),'b');
scatter3(v(endpoint,1),v(endpoint,2),v(endpoint,3),'g');
hold off

routeLength = 0;
for i = 1:size(spath,2)-1
    ddd = pdist([v(spath(i),1),v(spath(i),2),v(spath(i),3);...
        v(spath(i+1),1),v(spath(i+1),2),v(spath(i+1),3)]);
    routeLength = routeLength+ddd;
end
routeLength