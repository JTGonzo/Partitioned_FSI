function [cnnBC] = getElemBC(BC, cnn)

	[tf,idx] = ismember(BC(:,1),cnn(:,1));
	[tf,idx(:,2)] = ismember(BC(:,1),cnn(:,2));
	[tf,idx(:,3)] = ismember(BC(:,1),cnn(:,3));
	[tf,idx(:,4)] = ismember(BC(:,1),cnn(:,4));

	elem = unique(idx(:), 'stable');
% 	elem = unique(idx(:));
	elem(elem==0)=[];

	[tf,idx2] = ismember(BC(:,2),cnn(:,1));
	[tf,idx2(:,2)] = ismember(BC(:,2),cnn(:,2));
	[tf,idx2(:,3)] = ismember(BC(:,2),cnn(:,3));
	[tf,idx2(:,4)] = ismember(BC(:,2),cnn(:,4));

	elem2 = unique(idx2(:), 'stable');
% 	elem2 = unique(idx2(:));
	elem2(elem2==0)=[];
	   
    elem = intersect(elem, elem2, 'stable');

	cnnBC = cnn(elem,:);

	clear elem elem2 idx idx2

end