function [cnnBCnew] = reorderElemBC(BC, cnnBC)
	for i=1:size(BC,1)
	   [tf, idx1(i)] = ismember(BC(i,1),cnnBC(i,:));
	   cnnBCnew(i,1) = cnnBC(i,idx1(i));
	   if (idx1(i) == 1)
		   cnnBCnew(i,2) = cnnBC(i,2);
		   cnnBCnew(i,3) = cnnBC(i,3);
		   cnnBCnew(i,4) = cnnBC(i,4);    
	   elseif (idx1(i) == 2)
		   cnnBCnew(i,2) = cnnBC(i,3);
		   cnnBCnew(i,3) = cnnBC(i,4);
		   cnnBCnew(i,4) = cnnBC(i,1);
	   elseif (idx1(i) == 3)
		   cnnBCnew(i,2) = cnnBC(i,4);
		   cnnBCnew(i,3) = cnnBC(i,1);
		   cnnBCnew(i,4) = cnnBC(i,2);        
	   elseif (idx1(i) == 4)
		   cnnBCnew(i,2) = cnnBC(i,1);
		   cnnBCnew(i,3) = cnnBC(i,2);
		   cnnBCnew(i,4) = cnnBC(i,3);
	   end
	   [tf, idx2(i)] = ismember(BC(i,2),cnnBCnew(i,:));
	   if (idx2(i) == 4)
		   temp = cnnBCnew(i,1);
		   cnnBCnew(i,1) = cnnBCnew(i,4);
		   cnnBCnew(i,4) = cnnBCnew(i,3);
		   cnnBCnew(i,3) = cnnBCnew(i,2);
		   cnnBCnew(i,2) = temp;
       end	   
	   clear tf idx1 idx2		   
	end	
end