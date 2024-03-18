function [] = plt_export(MESH, uF, t, problemString)

Flag3D = 0;  
elemType = '3Tri';

%% Output the *.plt files
    if (mod(timeStep,solver.outFreq)==0)

		elemelem = MESH.elements;
					
		data = [MESH.nodes uF(:,1:2) uF(:,3) ];
		
        filename = sprintf('Figures/%s.%d.plt',problemString,timeStep);
		fileId = fopen(filename,'w');
		
        FileTitle = 'velocity plot';

		fprintf(fileId,' TITLE = \"%s\"\n',FileTitle);
		if (strcmp(FileTitle,'velocity plot') & Flag3D == 1)
			fprintf(fileId,' VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\", \"p\", \"vapFrac\"\n');
		elseif (strcmp(FileTitle,'velocity plot') & Flag3D ~= 1)
			fprintf(fileId,' VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"p\",  \n');
		end
		
		if (strcmp(elemType,'3Tri'))
			fprintf(fileId,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',t,size(MESH.nodes,1),size(elemelem,1));
		elseif (strcmp(elemType,'6Prism'))
			fprintf(fileId,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK\n',t,size(MESH.nodes,1),size(elemelem,1));
        end
		
	    fprintf(fileId,'%12.10f %12.10f %12.10f %12.10f %12.10f \n',data');

		if (strcmp(elemType,'3Tri'))
			fprintf(fileId,'%d %d %d\n',elemelem(1:3,:)');
        elseif (strcmp(elemType,'6Prism'))
			fprintf(fileId,'%d %d %d %d %d %d\n',elemelem(1:3,:)');
        end
        
		fclose(fileId);
    end   
end