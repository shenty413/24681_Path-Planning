clear all;
% close all;
%this script is written by Disheng Hou, a
%read obj files
% [v,f]=obj__read("donut.obj");
% v=v';
% f=f';

figure();

edges = [];
%extract edges according to faces
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
scatter3(v(:,3),v(:,2),v(:,1));
plot3(v(spath,1),v(spath,2),v(spath,3));
hold off
figure
hold on
plot3(v(spath,1),v(spath,2),v(spath,3));
hold off
