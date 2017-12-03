[v,f]=obj__read("donut.obj");
v=v';
f=f';
scatter3(v(:,3),v(:,2),v(:,1));
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