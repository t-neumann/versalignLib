#include "ocl_testing.h"

std::string get_device_type(cl_int const & device_type) {
	return device_type == CL_DEVICE_TYPE_GPU ? "GPU" : "CPU";
}

void run_ocl_test() {

// Testing OpenCL

	//get all platforms (drivers)
	std::vector<cl::Platform> all_platforms;
	cl::Platform::get(&all_platforms);
	if(all_platforms.size()==0){
		std::cout<<" No platforms found. Check OpenCL installation!\n";
		exit(-1);
	}
	cl::Platform default_platform=all_platforms[0];
	std::cout << "Using platform: "<<default_platform.getInfo<CL_PLATFORM_NAME>()<<"\n";

	//get default device of the default platform
	std::vector<cl::Device> all_devices;
	default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
	if(all_devices.size()>0){
		std::cout<< "Available devices:" << std::endl;
		for (std::vector<cl::Device>::iterator i = all_devices.begin(); i != all_devices.end(); ++i) {
			std::cout << i->getInfo<CL_DEVICE_NAME>() << "\n\n";
			std::cout << "\tType: " << get_device_type(i->getInfo<CL_DEVICE_TYPE>()) << std::endl;
            std::cout << "\tVendor: " << i->getInfo<CL_DEVICE_VENDOR>() << std::endl;
            std::cout << "\tMax Compute Units: " << i->getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << std::endl;
            std::cout << "\tGlobal Memory: " << i->getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() << std::endl;
            std::cout << "\tMax Clock Frequency: " << i->getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << std::endl;
            std::cout << "\tMax Allocateable Memory: " << i->getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() << std::endl;
            std::cout << "\tLocal Memory: " << i->getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << std::endl;
            std::cout << "\tAvailable: " << i->getInfo< CL_DEVICE_AVAILABLE>() << std::endl;
		}
		std::cout << std::endl;

	} else {
		std::cout<< "No devices found. Check OpenCL installation!\n";
		exit(-1);
	}

	cl::Device default_device = all_devices[0];
	std::cout<< "Using device: " << default_device.getInfo<CL_DEVICE_NAME>() << std::endl;
	std::vector<size_t> work_item_sizes = default_device.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>();
	for (std::vector<size_t>::iterator itor = work_item_sizes.begin(); itor != work_item_sizes.end(); ++itor) {
		std::cout<< "Available work item sizes: " << *itor << std::endl;
	}
	std::cout<< "Build in kernels: " << default_device.getInfo<CL_DEVICE_BUILT_IN_KERNELS>() << std::endl;

	cl::Context context(default_device);

	cl::Program::Sources sources;

	// kernel calculates for each element C=A+B
	std::string kernel_code=
			"   #pragma OPENCL EXTENSION cl_amd_printf : enable\n"
			"   #pragma OPENCL EXTENSION cl_khr_byte_addressable_store: enable\n"
			"   void kernel simple_add(global const int* A, global const int* B, global int* C){\n"
			"   printf(\"hello world A %i B %i \", A[get_global_id(0)], B[get_global_id(0)]);\n"
			"       C[get_global_id(0)]=A[get_global_id(0)]+B[get_global_id(0)];\n"
			"   }                                                                               ";

	sources.push_back(std::make_pair(kernel_code.data(), kernel_code.length()));

	cl::Program program(context,sources);
	if(program.build()!=CL_SUCCESS){
		std::cout<<" Error building: "<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device)<<"\n";
		exit(1);
	}
	// create buffers on the device
	cl::Buffer buffer_A(context,CL_MEM_READ_WRITE,sizeof(int)*10);
	cl::Buffer buffer_B(context,CL_MEM_READ_WRITE,sizeof(int)*10);
	cl::Buffer buffer_C(context,CL_MEM_READ_WRITE,sizeof(int)*10);

	int A[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	int B[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0};

	//create queue to which we will push commands for the device.
	cl::CommandQueue queue(context,default_device);

	//write arrays A and B to the device
	queue.enqueueWriteBuffer(buffer_A,CL_TRUE,0,sizeof(int)*10,A);
	queue.enqueueWriteBuffer(buffer_B,CL_TRUE,0,sizeof(int)*10,B);

	//alternative way to run the kernel
	cl::Kernel kernel_add=cl::Kernel(program,"simple_add");
	kernel_add.setArg(0,buffer_A);
	kernel_add.setArg(1,buffer_B);
	kernel_add.setArg(2,buffer_C);
	queue.enqueueNDRangeKernel(kernel_add,cl::NullRange,cl::NDRange(10),cl::NullRange);
	queue.finish();

	int C[10];
	//read result C from the device to array C
	queue.enqueueReadBuffer(buffer_C,CL_TRUE,0,sizeof(int)*10,C);

	std::cout<<" result: \n";
	for(int i=0;i<10;i++){
		std::cout<<C[i]<<" ";
	}
}
