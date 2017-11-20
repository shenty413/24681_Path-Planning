function idx = FindMinDistToStart(nState,distToStart)
% Function
%   idx = FindMinDistToStart(nState,distToStart)
%
% Input
%   nState: n*1 vertice state array
%   distToStart: 
%
% Output
%   idx: the index of the vertice
%
    idx = distToStart~=0 & nState==0;
    [value,i] = min(distToStart(idx));
    f = find(distToStart~=0 & nState==0);
    idx = f(i);
end
