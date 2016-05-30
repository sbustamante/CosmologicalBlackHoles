#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include "hdf5.h"
#define filename "./snapshot"

struct Particle{
	
	double xp;
	double yp;
	double zp;
	double mp;

  }*parts;

struct param{

	int Lbox;
	int N_files;
	int N_part;

}PRM;  

int main(  ){

	int i,j,counter,rank=2;
	double sdata[1][3],MT[6];
	char box_filename[80];
	int dumb,dumb1[6],dumb2[6],dumb0;

	hsize_t dim[2],count[2],offset[2],stride[2],block[2];
	hsize_t memspace_id, dataspace_id, dset_id;
	hid_t att, root, file_id, dataset_id;
	herr_t status;

	dim[0]=1,dim[1]=3;
	//Determines how many blocks to select from the 
	//dataspace, in ea	ch dimension
    count[0]=1,count[1]=3;
    //Steps beteween two elements
    stride[0]=1,stride[1]=1;
    block[0]=1,block[1]=1;

	printf( "\n\n***********************************************************************\n" ); 
	printf("Reading HDF5 snapshot %s\n", filename );	

	//Open an existing file
	sprintf(box_filename, "%s.%d.hdf5", filename,0);
	file_id = H5Fopen(box_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	
	root = H5Gopen(file_id, "Header",H5P_DEFAULT);

	//Reading total particle number
	att = H5Aopen(root,"NumPart_Total", H5P_DEFAULT); 
	status = H5Aread( att,H5T_NATIVE_INT,dumb1 );
	PRM.N_part = dumb1[1];
	
	//Memory allocation for part variable
	parts = (struct Particle *) calloc( (size_t) PRM.N_part,sizeof(struct Particle ) );

	//Reading Lbox 
	att = H5Aopen_name(root, "BoxSize"); 
	status = H5Aread(att,H5T_NATIVE_INT, &dumb);
	PRM.Lbox = dumb;

	//Reading DM mass 
	att = H5Aopen_name(root, "MassTable"); 
	status = H5Aread(att,H5T_NATIVE_DOUBLE, MT);

	//Reading Num Files Per Snapshot
	att = H5Aopen_name(root, "NumFilesPerSnapshot"); 
	status = H5Aread(att,H5T_NATIVE_INT,&dumb0 );
	PRM.N_files = dumb0;

	j = 0;
	//for ( counter = 0; counter < PRM.N_files; counter++ ){		
	for ( counter = 0; counter < 1; counter++ ){		

		if (counter != 0){
		   sprintf(box_filename, "%s.%d.hdf5", filename,0);
		   file_id = H5Fopen(box_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
		   root = H5Gopen(file_id, "Header",H5P_DEFAULT);
		   }

		//Reading Particles in this file
		att = H5Aopen(root,"NumPart_ThisFile", H5P_DEFAULT); 
		status = H5Aread( att,H5T_NATIVE_INT,dumb2 );

		// Close the attribute
    	status = H5Aclose( att );
    	// Close the Group
		status = H5Gclose( root );	 

		//Open an existing dataset	
		dataset_id = H5Dopen2(file_id, "/PartType1/Coordinates", H5P_DEFAULT);
					   
		memspace_id = H5Screate_simple(rank, dim, NULL); 
	    dataspace_id = H5Dget_space (dataset_id);

		for ( i = 0; i < dumb2[1]; i++ ){	
			
			//Specifies the offset of the starting element 
			//of the specified hyperslab.
		    offset[0] = i;
	    	offset[1] = 0;

	    	status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, 
	    								 stride, dim, block);

	    	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, 
	    								 H5P_DEFAULT, sdata);
	    	parts[j].xp = sdata[0][0];
	    	parts[j].yp = sdata[0][1];
	    	parts[j].zp = sdata[0][2];
	    	parts[j].mp = MT[1];
	   		j++;
			
		}
		
	    status = H5Sclose (memspace_id);
    	status = H5Sclose (dataspace_id);
    	//Close dataset
		status = H5Dclose( dataset_id );
		// Close the file
		status = H5Fclose(file_id);

   } 

return 0;	
}
