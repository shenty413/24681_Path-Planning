function [d,dToS,rArray] = updateDistToStart(p1,vLast,vdata,dToS,rArray)
% Function
%   [d,distToStart,index] = updateDistToStart(p1,vLast,vdata,dToS,dToL,route)
%   updates the distance from p1 to start point
%
% Input
%   p1: given point index
%   vLast: the last chosen point index
%   vdata: n*3 vdata
%   dToS: n*1 distToStart array
%   dToL: n*1 distToLast array
%   rArray: n*1 route array that represents all the routes from start point
%           to all the points
% Output
%   d: distance from p1 to start
%   dToS: updated n*1 distToStart array
%   rArray: updated n*1 route array
%
% Example:
%   [dist(i),distToStart,route] = updateDistToStart(vAroundvLast(i),vLast,vdata,distToStart,route);
%

    x(1,:) = vdata(p1,:);
    x(2,:) = vdata(vLast,:);
    if pdist(x)+dToS(vLast)<dToS(p1) || dToS(p1)==0
        % if the distance from last point to this point
        % plus the last point to start point is smaller
        % or the point has not been reached before,
        % update the distance to start point
        dToS(p1) = pdist(x)+dToS(vLast);
        rArray{p1} = [rArray{vLast} p1];
    end
    d = dToS(p1);
end