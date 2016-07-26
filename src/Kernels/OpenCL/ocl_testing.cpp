#include "ocl_testing.h"
#include "opencl_definitions.h"
#include "scoring_kernels.h"
#include "alignment_kernels.h"
#include "AlignmentKernel.h"

#include <sstream>
#include <cstdio>

using std::stringstream;

extern char const opencl_definitions[];
extern char const scoring_kernels[];
extern char const alignment_kernels[];

std::string get_device_type(cl_int const & device_type) {
	return device_type == CL_DEVICE_TYPE_GPU ? "GPU" : "CPU";
}

void run_ocl_test(char const * const * const reads,
		char const * const * const refs, int const & seqNumber,
		size_t const & max_read_length, size_t const & max_ref_length) {

// Testing OpenCL

//get all platforms (drivers)
	std::vector<cl::Platform> all_platforms;
	cl::Platform::get(&all_platforms);
	if (all_platforms.size() == 0) {
		std::cout << " No platforms found. Check OpenCL installation!\n";
		exit(-1);
	}
	cl::Platform default_platform = all_platforms[0];
	std::cout << "Using platform: "
			<< default_platform.getInfo<CL_PLATFORM_NAME>() << "\n";

	//get default device of the default platform
	std::vector<cl::Device> all_devices;
	default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
	if (all_devices.size() > 0) {
		std::cout << "Available devices:" << std::endl;
		for (std::vector<cl::Device>::iterator i = all_devices.begin();
				i != all_devices.end(); ++i) {
			std::cout << i->getInfo<CL_DEVICE_NAME>() << "\n\n";
			std::cout << "\tType: "
					<< get_device_type(i->getInfo<CL_DEVICE_TYPE>())
					<< std::endl;
			std::cout << "\tVendor: " << i->getInfo<CL_DEVICE_VENDOR>()
					<< std::endl;
			std::cout << "\tMax Compute Units: "
					<< i->getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << std::endl;
			std::cout << "\tGlobal Memory: "
					<< i->getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() << std::endl;
			std::cout << "\tMax Clock Frequency: "
					<< i->getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << std::endl;
			std::cout << "\tMax Allocateable Memory: "
					<< i->getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() << std::endl;
			std::cout << "\tLocal Memory: "
					<< i->getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << std::endl;
			std::cout << "\tAvailable: " << i->getInfo< CL_DEVICE_AVAILABLE>()
					<< std::endl;
		}
		std::cout << std::endl;

	} else {
		std::cout << "No devices found. Check OpenCL installation!\n";
		exit(-1);
	}

	cl::Device default_device = all_devices[0]; // Pick CPU

//	// Calculate number of used cores/hyperthreading threads
//	// Always only use 3 quarters for computation and leave 1 quarter for background processes
//	int max_devices = default_device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
//	int fission = max_devices / 4 * 3;
//
//	// Partition CPU
//
//	std::vector<cl::Device> subdevices;
//
//	cl_device_partition_property props[4];
//	props[0] = CL_DEVICE_PARTITION_BY_COUNTS;
//	props[1] = std::max(fission, 1);
//	props[2] = CL_DEVICE_PARTITION_BY_COUNTS_LIST_END;
//	props[3] = 0;
//
//	default_device.createSubDevices(props, &subdevices);
//
//	std::cout << "Created " << subdevices.size() << " subdevices.\n";
//
//	default_device = subdevices[0];
	std::cout << "Using device: " << default_device.getInfo<CL_DEVICE_NAME>()
			<< std::endl;
	std::vector<size_t> work_item_sizes = default_device.getInfo<
	CL_DEVICE_MAX_WORK_ITEM_SIZES>();
	for (std::vector<size_t>::iterator itor = work_item_sizes.begin();
			itor != work_item_sizes.end(); ++itor) {
		std::cout << "Available work item sizes: " << *itor << std::endl;
	}
	std::cout << "Max work group size: "
			<< default_device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>()
			<< std::endl;
	std::cout << "Build in kernels: "
			<< default_device.getInfo<CL_DEVICE_BUILT_IN_KERNELS>()
			<< std::endl;

	cl::Context context(default_device);

	cl::Program::Sources sources;

	// Load Kernel

	//std::stringstream kernel_code_hex;
	//kernel_code_hex << test_kernel;
	//std::string kernel_code_hex (test_kernel);

	std::stringstream test;
	test << opencl_definitions << scoring_kernels << alignment_kernels;

	std::string input(test.str());

	//sources.push_back(std::make_pair(kernel_code.data(), kernel_code.length()));
	sources.push_back(std::make_pair(input.data(), input.length()));

	stringstream compilerDefines;
	compilerDefines << "-D read_length=" << max_read_length << " -D ref_length="
			<< max_ref_length << " -D aln_length="
			<< max_ref_length + max_read_length << " -D score_gap_read=-3"
			<< " -D score_gap_ref=-3" << " -D score_match=2"
			<< " -D score_mismatch=-1";

	cl::Program program(context, sources);

	if (program.build(compilerDefines.str().c_str()) != CL_SUCCESS) {
		std::cout << " Error building: "
				<< program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device)
				<< "\n";
		exit(1);
	}
//	// create buffers on the device
//	cl::Buffer buffer_A(context,CL_MEM_READ_WRITE,sizeof(int)*10);
//	cl::Buffer buffer_B(context,CL_MEM_READ_WRITE,sizeof(int)*10);
//	cl::Buffer buffer_C(context,CL_MEM_READ_WRITE,sizeof(int)*10);
//
//	int A[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
//	int B[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0};
//
//	//create queue to which we will push commands for the device.
	cl::CommandQueue queue(context, default_device);
//
//	//write arrays A and B to the device
//	queue.enqueueWriteBuffer(buffer_A,CL_TRUE,0,sizeof(int)*10,A);
//	queue.enqueueWriteBuffer(buffer_B,CL_TRUE,0,sizeof(int)*10,B);
//
//	//alternative way to run the kernel
//	cl::Kernel kernel_add=cl::Kernel(program,"simple_add");
//	kernel_add.setArg(0,buffer_A);
//	kernel_add.setArg(1,buffer_B);
//	kernel_add.setArg(2,buffer_C);
//	queue.enqueueNDRangeKernel(kernel_add,cl::NullRange,cl::NDRange(10),cl::NullRange);
//	queue.finish();
//
//	int C[10];
//	//read result C from the device to array C
//	queue.enqueueReadBuffer(buffer_C,CL_TRUE,0,sizeof(int)*10,C);
//
//	//std::cout<<" result: \n";
//	for(int i=0;i<10;i++){
//		//std::cout<<C[i]<<" ";
//	}

	//alternative way to run the kernel
	cl::Kernel kernel_print = cl::Kernel(program,
			"calc_alignment_needleman_wunsch");

	const size_t max_work_items = kernel_print.getWorkGroupInfo<
			CL_KERNEL_WORK_GROUP_SIZE>(default_device);

	std::cout << "Kernel group size: "
			<< kernel_print.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(
					default_device) << std::endl;

	std::cout << "Number of reads: " << seqNumber << std::endl;
	std::cout << "GlobalMax mem:\t"
			<< default_device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() << std::endl;
	std::cout << "GlobalAll mem:\t"
			<< default_device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>()
			<< std::endl;

	const int vectors_per_workitem = 16;

	const int aln_length = max_ref_length + max_read_length;
	const int matrix_size = (max_read_length + 1) * (max_ref_length + 1);

	size_t const _500MB = 1024 * 1024;
	//size_t const _500MB = 1024 * 1024 * 500;
//	size_t const _one_alignment = sizeof(char) * max_read_length
//			+ sizeof(char) * max_ref_length + sizeof(short) * 2;
	size_t const _one_alignment = sizeof(char) * max_read_length
			+ sizeof(char) * max_ref_length + sizeof(char) * aln_length * 2
			+ sizeof(short) * 2 + sizeof(short) * matrix_size;

	std::cout << "500MB in bytes:\t" << _500MB << std::endl;
	std::cout << "one in bytes:\t" << _one_alignment << std::endl;

	size_t batch_size =
			((_500MB / _one_alignment) / vectors_per_workitem
					* vectors_per_workitem) == 0 ?
					vectors_per_workitem :
					(_500MB / _one_alignment) / vectors_per_workitem
							* vectors_per_workitem;

	batch_size = std::min(batch_size, max_work_items * vectors_per_workitem);

	std::cout << "Batch size:\t" << batch_size << std::endl;

	size_t const num_batches = seqNumber / batch_size;
	size_t const overhang = seqNumber % batch_size;

	std::cout << "Num batches:\t" << num_batches << std::endl;
	std::cout << "Overhang:\t" << overhang << std::endl;
	std::cout << "Aln length:\t" << aln_length << std::endl;
	std::cout << "Read length:\t" << max_read_length << std::endl;
	std::cout << "Ref length:\t" << max_ref_length << std::endl;

	char * host_reads = new char[batch_size * max_read_length];
	char * host_refs = new char[batch_size * max_ref_length];
	//		short * host_results = new short[batch_size]();
	char * host_results = new char[aln_length * 2 * batch_size];
	short * host_indices = new short[2 * batch_size];
	short * host_matrix = new short[matrix_size * batch_size];

	char a;
	std::cin >> a;

	for (int batch = 0; batch < num_batches; ++batch) {

		memset(host_reads, 0, sizeof(char) * batch_size * max_read_length);
		memset(host_refs, 0, sizeof(char) * batch_size * max_ref_length);
		memset(host_results, 0, sizeof(char) * batch_size * aln_length * 2);
		memset(host_indices, 0, sizeof(short) * batch_size * 2);
		memset(host_matrix, 0, sizeof(short) * matrix_size * batch_size);

		for (int i = 0; i < batch_size; ++i) {
			memcpy(&host_reads[i * max_read_length],
					reads[batch * batch_size + i],
					sizeof(char) * max_read_length);
			std::cout << "Read: " << reads[batch * batch_size + i] << std::endl;
			memcpy(&host_refs[i * max_ref_length], refs[batch * batch_size + i],
					sizeof(char) * max_ref_length);
			std::cout << "Ref: " << refs[batch * batch_size + i] << std::endl;
		}

		cl::Buffer read_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				sizeof(char) * batch_size * max_read_length, host_reads);
		cl::Buffer ref_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				sizeof(char) * batch_size * max_ref_length, host_refs);
		//		cl::Buffer result_buffer(context,
		//		CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(short) * batch_size,
		//				host_results);
		cl::Buffer result_buffer(context,
		CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
				sizeof(char) * aln_length * 2 * batch_size, host_results);
		cl::Buffer index_buffer(context,
		CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(short) * batch_size * 2,
				host_indices);
		cl::Buffer matrix_buffer(context,
		CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
				sizeof(short) * matrix_size * batch_size, host_matrix);

		//printf("Reads %.*s\n", max_read_length, host_reads);
		//printf("Refs %.*s\n", max_ref_length, host_refs);

		//printf("A %s B %s ", host_reads, host_refs);
		//std::cout << host_reads << std::endl << host_refs << std::endl;

		//cl::Buffer read_buffer(context,CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,sizeof(char)*seqNumber*max_read_length);
		//cl::Buffer ref_buffer(context,CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,sizeof(char)*seqNumber*max_ref_length);

		//queue.enqueueMapBuffer(read_buffer,CL_TRUE,CL_MAP_READ,0,sizeof(char)*seqNumber*max_read_length);
		//queue.enqueueMapBuffer(ref_buffer,CL_TRUE,CL_MAP_READ,0,sizeof(int)*seqNumber*max_ref_length);

		std::cout << "Running \"calc_alignment_smith_waterman\" Kernel...\n";

		kernel_print.setArg(0, read_buffer);
		kernel_print.setArg(1, ref_buffer);
		kernel_print.setArg(2, result_buffer);
		kernel_print.setArg(3, index_buffer);
		kernel_print.setArg(4, matrix_buffer);
		int stuff = batch_size / 16;
		std::cout << "Work groups: " << stuff << std::endl;
		queue.enqueueNDRangeKernel(kernel_print, cl::NullRange,
				cl::NDRange(stuff), cl::NullRange);
		queue.finish();
		std::cout << "Finished Kernel.\n";
		//std::cin >> a;

		Alignment * alignment = new Alignment[batch_size];

		for (int i = 0; i < batch_size; ++i) {
			char * read_back = new char[max_read_length];
			memcpy(read_back, host_reads + i * max_read_length,
					max_read_length);
			std::cout << "Read: " << read_back << std::endl;
			char * ref_back = new char[max_ref_length];
			memcpy(ref_back, host_refs + i * max_ref_length, max_ref_length);
			std::cout << "Ref: " << ref_back << std::endl;
			//std::cout << "Alignment score: " << host_results[i] << std::endl;
			//std::cin >> a;
			alignment[i].read = new char[aln_length];
			alignment[i].ref = new char[aln_length];

			memcpy(alignment[i].read, host_results + 2 * i * aln_length,
					aln_length * sizeof(char));
			memcpy(alignment[i].ref,
					host_results + 2 * i * aln_length + aln_length,
					aln_length * sizeof(char));

			alignment[i].readStart = host_indices[2 * i];
			alignment[i].refStart = host_indices[2 * i + 1];

			alignment[i].readEnd = aln_length - 1;
			alignment[i].refEnd = aln_length - 1;

			std::cout << "==================" << std::endl << "\"";
			std::cout << "Start:\t" << host_indices[2 * i] << "\t"
					<< alignment[i].read + alignment[i].readStart;
			std::cout << "\"" << std::endl << "\"";
			std::cout << "End:\t" << host_indices[2 * i + 1] << "\t"
					<< alignment[i].ref + alignment[i].refStart;
			std::cout << "\"" << std::endl << "==================" << std::endl;
			//std::cin >> a;

		}
	}

	if (overhang != 0) {

		std::cout << "Only overhang!\n";

		memset(host_reads, 0, sizeof(char) * batch_size * max_read_length);
		memset(host_refs, 0, sizeof(char) * batch_size * max_ref_length);
		memset(host_results, 0, sizeof(char) * batch_size * aln_length * 2);
		memset(host_indices, 0, sizeof(short) * batch_size * 2);
		memset(host_matrix, 0, sizeof(short) * matrix_size * batch_size);

		for (int remainder = 0; remainder < overhang; ++remainder) {

			memcpy(&host_reads[remainder * max_read_length],
					reads[num_batches * batch_size + remainder],
					sizeof(char) * max_read_length);
			std::cout << "Read: " << reads[num_batches * batch_size + remainder]
					<< std::endl;
			memcpy(&host_refs[remainder * max_ref_length],
					refs[num_batches * batch_size + remainder],
					sizeof(char) * max_ref_length);
			std::cout << "Ref: " << refs[num_batches * batch_size + remainder]
					<< std::endl;
		}

		cl::Buffer read_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				sizeof(char) * batch_size * max_read_length, host_reads);
		cl::Buffer ref_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				sizeof(char) * batch_size * max_ref_length, host_refs);
		//		cl::Buffer result_buffer(context,
		//		CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(short) * batch_size,
		//				host_results);
		cl::Buffer result_buffer(context,
		CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
				sizeof(char) * aln_length * 2 * batch_size, host_results);
		cl::Buffer index_buffer(context,
		CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(short) * batch_size * 2,
				host_indices);
		cl::Buffer matrix_buffer(context,
		CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
				sizeof(short) * matrix_size * batch_size, host_matrix);

		std::cout << "Running \"calc_alignment_smith_waterman\" Kernel...\n";
		//alternative way to run the kernel

		kernel_print.setArg(0, read_buffer);
		kernel_print.setArg(1, ref_buffer);
		kernel_print.setArg(2, result_buffer);
		kernel_print.setArg(3, index_buffer);
		kernel_print.setArg(4, matrix_buffer);
		int stuff = batch_size / 16;
		std::cout << "Running with " << stuff << " workers.\n";

		queue.enqueueNDRangeKernel(kernel_print, cl::NullRange,
				cl::NDRange(stuff), cl::NullRange);
		queue.finish();

		Alignment * alignment = new Alignment[overhang];

		for (int i = 0; i < overhang; ++i) {
			char * read_back = new char[max_read_length];
			memcpy(read_back, host_reads + i * max_read_length,
					max_read_length);
			std::cout << "Read: " << read_back << std::endl;
			char * ref_back = new char[max_ref_length];
			memcpy(ref_back, host_refs + i * max_ref_length, max_ref_length);
			std::cout << "Ref: " << ref_back << std::endl;
//			std::cout << "Alignment score: " << host_results[i] << std::endl;
			alignment[i].read = new char[aln_length];
			alignment[i].ref = new char[aln_length];

			memcpy(alignment[i].read, host_results + 2 * i * aln_length,
					aln_length * sizeof(char));
			memcpy(alignment[i].ref,
					host_results + 2 * i * aln_length + aln_length,
					aln_length * sizeof(char));

			alignment[i].readStart = host_indices[2 * i];
			alignment[i].refStart = host_indices[2 * i + 1];

			alignment[i].readEnd = aln_length - 1;
			alignment[i].refEnd = aln_length - 1;

			std::cout << "==================" << std::endl << "\"";
			std::cout << "Start:\t" << host_indices[2 * i] << "\t"
					<< alignment[i].read + alignment[i].readStart;
			std::cout << "\"" << std::endl << "\"";
			std::cout << "End:\t" << host_indices[2 * i + 1] << "\t"
					<< alignment[i].ref + alignment[i].refStart;
			std::cout << "\"" << std::endl << "==================" << std::endl;
			//std::cin >> a;

		}
	}
}
