#pragma OPENCL EXTENSION cl_amd_printf : enable
#pragma OPENCL EXTENSION cl_khr_byte_addressable_store: enable
			
void kernel simple_add(global const int* A, global const int* B, global int* C){

	//printf("A %i B %i ", A[get_global_id(0)], B[get_global_id(0)]);
	
	C[get_global_id(0)]=A[get_global_id(0)]+B[get_global_id(0)];
	
}

void kernel print_string(global const char* A, global const char* B){

	printf("ID-A %i ID-B %i ", get_global_id(0), get_global_id(0));
	printf("A %s B %s ", A[get_global_id(0)], B[get_global_id(0)]);
		
}