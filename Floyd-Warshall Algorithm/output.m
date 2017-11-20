function [path] = output(s,e,paths)
path=[];
if(paths(s,e)==0)
    path = [path,e];
else
    path = [path,output(s,paths(s,e),paths)];
     path = [path,output(paths(s,e),e,paths)];
end
end

