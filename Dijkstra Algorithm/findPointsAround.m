function [vAroundvLast] = findPointsAround(vLast,vInFaceIdx,fdata,nState,vAroundvLast)
%     vAroundvLast = findPointsAround(vLast,vInFaceIdx,fdata,nState)
%     6 - - - 5
%     | \     |
%     |  \    |
%     |   3 - 4
%     |  /|  /
%     | / | /
%     1 - 2
%     [2,3,6] = findPointsAround(1,vInFaceIdx,fdata,nState)

    fAroundvLast = vInFaceIdx(vLast,:);
    [value,idx] = min(fAroundvLast);
    for i = 1:idx-1
        vv(i,:) = fdata(fAroundvLast(i),:);
    end
    % find all vertice around given v, including those have been chosen and
    % the last point
    vAroundvLastCopy = unique(vv);
    
    n = 0;
    length = max([size(vAroundvLastCopy,1) size(vAroundvLastCopy,2)]);
%     if nState(vAroundvLastCopy)==ones(length,1)
%         stopFlag = 1;
% %         vAroundvLast = vAroundvLast;
%     else
            for i = 1:length
                if isnan(vAroundvLastCopy(i))==1
                    break;
                end
                if nState(vAroundvLastCopy(i))==0 && vAroundvLastCopy(i)~=vLast
                    % find those that need to be 
                    n = n+1;
                    vAroundvLast(n)=vAroundvLastCopy(i);
        %         else
        %             stopFlag = 1;
        %             break;
                end
            end
%     end

end