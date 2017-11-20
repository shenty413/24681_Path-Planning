function [vdata,fdata,numV] = getDataFromFile(fileName,vStartFrom0)
% Function
%   [vdata,fdata,numV] = getDataFromFile(fileName,vStartFrom0)
%
% Input
%   fileName: '.txt'
%   vStartFrom0: boolen value, true if face data v starts from 0
% Output
%   vdata: n*3 vertice data matrix
%   fdata: m*3 face data matrix
%   numV: number of vertices
%
% Example
%   [vdata,fdata,numV] = getDataFromFile('triceratops.dat',true)
%

    %%% open the file 
    fid = fopen(char(fileName));
    C = textscan(fid, '%s%f%f%f');
    DataColumn = 3;
    fclose(fid);

    %%% read the data
    type=C{1};
    for n = 1:DataColumn
        % if we directly read data from the file, the last row of some data
        % may be omitted. So we need to verify the size of the matrix first

        if size(C{n+1})==size(C{1})
            % if the size is the same as the first column, store the data
            data(:,n) = C{n+1};
        else
            % if the size is different, add the blank array with NaN
            for j = 1:(size(C{1},1)-size(C{n+1},1))
                datax(:,1) = C{n+1};
                datax(size(datax,1)+j,1) = NaN;
            end
            % store the adjusted data
            data(:,n) = datax(:,1);
        end
        clear datax;
    end

    %%% for the vertices
    % initialize the variables
    numV = 0;
    length = size(type,1);

    for n = 1:length
        % get the number of vertices
        if strcmp(type(n),'v')
            numV = numV+1;
        end
    end

    vdata = data(1:numV,:);
    if vStartFrom0==true
        fdata = data(numV+1:end,:)+1;
    else
        fdata = data(numV+1:end,:);
    end
end