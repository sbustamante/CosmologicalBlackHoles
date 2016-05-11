#include "hdf5.h"
#include <stdlib.h>
#include <math.h>
#define FILE "snapshot.hdf5"

#define NMAX 100000

int main()
{
    int i, ndims;
    char box_filename[80];
    hsize_t dims[1] = {NMAX};
    hid_t file_id, dataset_id,memtype, space; 
    herr_t status;
	
    double **dset_data = malloc(NMAX * sizeof(double *));
    dset_data[0] = malloc(NMAX * 3 * sizeof(double));
    
    for( i = 1; i < NMAX; i++ )
	dset_data[i] = dset_data[0] + i * 3;

    /* Open an existing file. */
    file_id = H5Fopen(FILE, H5F_ACC_RDWR, H5P_DEFAULT);

    /* Open an existing dataset. */
    dataset_id = H5Dopen2(file_id, "/PartType1/Coordinates", H5P_DEFAULT);

    /* Read the dataset. */
    status = H5Dread( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *dset_data ); 
    printf("%lf\n",dset_data[0][0] );

    /* Close the dataset. */
    status = H5Dclose( dataset_id );

    /* Close the file. */
    status = H5Fclose(file_id);
   
    free(*dset_data); 
    
return 0;	
}
