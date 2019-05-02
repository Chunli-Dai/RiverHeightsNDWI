function [Meanshift_Pts Window_Size_Pts]= MeanShiftCentralLine(dataPts,wMinRiverW);
%perform MeanShift Clustering of data using a flat kernel
%
% ---INPUT---
% dataPts           - input data, (numDim x numPts)
% wMinRiverW         - is the minimum river width
%---OUTPU---
% Meanshift_Pts     - the location of shifted points
% Window_Size_Pts   - the corresponding window size used at each point

% Written by C.Zeng on 3 Nov 2014, for River central line detection


%*** Check input ****
if nargin < 2
    error('no sWindow specified')
end

%**** Initialize stuff ***
[numDim,numPts] = size(dataPts);

%define the two window size: for point density and cluster
%the process of one cluster by given a random location.
    n=10;  %the widest river is 10 times of the minimum River
    k=1:n;
    window_point_density_size=(2*k+1).*wMinRiverW;
    window_cluster_size=k.*wMinRiverW;
  
    %the vector to store the width of each point
    Window_Size_Pts=zeros(numPts,1);
    
    progress= fix(numPts/10);
    h = waitbar(0,'Initializing waitbar...');
for index=1:numPts
    
    myMean=dataPts(:,index);
    %the shift of the data
    shiftData=dataPts-repmat(myMean,1,numPts);

    %search for the appropriate point density and given an appropriate
    %window size
    for idx=1:n  %find the best window size iteratively
        tempWindowSize=window_point_density_size(idx);

        %%used for uniform(square) kernal
        tempWindowKernal=tempWindowSize/2; %the width of the kernal is only half of the box width
        idx_candidates=find(shiftData(1,:)>= -1* tempWindowKernal & ...
            shiftData(1,:)<= tempWindowKernal & ...
            shiftData(2,:)>= -1* tempWindowKernal &...
            shiftData(2,:)<= tempWindowKernal);  %find points in the box
        nPts=length(idx_candidates);
        point_density= nPts/tempWindowSize^2;

        %decide whether it is the end of the iteration
        pt_Threshold=idx/(2*idx+1);
        if point_density < pt_Threshold  
            % Window_Size_Pts(index)=window_cluster_size(idx);
            break;
        end
    end
    sWindow= window_cluster_size(idx); 
    Window_Size_Pts(index)=sWindow;
    %find the points in the cluster window
    sWindow_kernal=sWindow/2; %the width of the kernal is only half of the box width
    idx_candidates= shiftData(1,:)>= -1* sWindow_kernal & ...
            shiftData(1,:)<= sWindow_kernal & ...
            shiftData(2,:)>= -1* sWindow_kernal &...
            shiftData(2,:)<= sWindow_kernal;
    
    %calcuate the new Mean in the neighbour area
    NewMean=mean(dataPts(:,idx_candidates),2);
    %%%%replace the curent point with its mean in the neighbours
    dataPts(:,index)=NewMean;
    
    if mod(index, progress)==0
        perc=index/numPts;
        waitbar(perc,h,sprintf('%d%% along...',ceil(perc*100)))
    end
end

    close(h)
    Meanshift_Pts=dataPts;
end
