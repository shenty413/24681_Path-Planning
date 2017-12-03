function drawTriangle(v1,v2,v3)
% Function
%   drawTriangle(v1,v2,v3) draws the triangle
    plot3([v1(1) v2(1) v3(1) v1(1)],[v1(2) v2(2) v3(2) v1(2)],[v1(3) v2(3) v3(3) v1(3)],'k-');
    hold on;
end