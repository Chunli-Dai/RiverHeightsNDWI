function nedgelist=Improved_CentralLine(Mask_Path_filter,Max_River_width)
%
% this function is designed to extract the river centerline
% INput parameter: 
% Mask_Path_filter: 	the final water mask after river connection
% Max_River_width:		the  estimated width of the river, it is used to remove spurs.
%
% Output parameter:
% nedgelist: the organized edge list that represents the centrelines, in vector format.
%
	%%-----extract the central line of the connected river mask -----
	Img_thin= bwmorph(Mask_Path_filter,'thin',Inf);
	%link the edge
	edgelist = edgelink(Img_thin);
	%clean the edge
	nedgelist = cleanedgelist(edgelist, Max_River_width/2);
	Old_edge_num=0;
	while Old_edge_num~=length(nedgelist)  %iterate till stable edge numbers
		Old_edge_num=length(nedgelist);
		%nedgelist = cleanedgelist(edgelist, Max_River_width/2);
		nedgelist = cleanedgelist(nedgelist, Max_River_width/2);
	end
	clear edgelist

end