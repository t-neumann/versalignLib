/*
 * OpenCLKernel.cpp
 *
 *  Created on: Jul 22, 2016
 *      Author: tobias.neumann
 */

#include "OpenCLKernel.h"

#include "opencl_definitions.h"
#include "scoring_kernels.h"
#include "alignment_kernels.h"

#include <sstream>
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>

using std::stringstream;

extern char const opencl_definitions[];
extern char const scoring_kernels[];
extern char const alignment_kernels[];

void OpenCLKernel::score_alignments(int const & opt, int const & aln_number,
		char const * const * const reads, char const * const * const refs,
		short * const scores) {

	int alignment_algorithm = opt & 0xF;

	char const * kernel_name = "";

	switch (alignment_algorithm) {
	case 0:
		kernel_name = score_smith_waterman_kernel;
		break;
	case 1:
		kernel_name = score_needleman_wunsch_kernel;
		break;
	default:
		// Unsupported mode
		break;
	}

	kernel = setup_kernel(program, kernel_name);

	size_t batch_size = calculate_batch_size_from_memory(kernel, device, true);

	size_t batch_num;
	size_t overhang;

	partition_load(aln_number, batch_size, batch_num, overhang);

	for (int batch = 0; batch < batch_num; ++batch) {

		init_host_memory(batch_size, true);

		for (int i = 0; i < batch_size; ++i) {
			memcpy(&host_reads[i * readLength], reads[batch * batch_size + i],
					sizeof(char) * readLength);
			std::cout << "Read: " << reads[batch * batch_size + i] << std::endl;
			memcpy(&host_refs[i * refLength], refs[batch * batch_size + i],
					sizeof(char) * refLength);
			std::cout << "Ref: " << refs[batch * batch_size + i] << std::endl;
		}

		cl::Buffer read_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				sizeof(char) * batch_size * readLength, host_reads);
		cl::Buffer ref_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				sizeof(char) * batch_size * refLength, host_refs);
		cl::Buffer result_buffer(context,
		CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(short) * batch_size,
				host_scores);

		std::cout << "Running \"" << kernel_name << "\" Kernel...\n";

		kernel.setArg(0, read_buffer);
		kernel.setArg(1, ref_buffer);
		kernel.setArg(2, result_buffer);

		int workers = batch_size / VECTORS_PER_WORKITEM;

		std::cout << "Work groups: " << workers << std::endl;

		queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(workers),
				cl::NullRange);
		queue.finish();

		std::cout << "Finished batch " << batch << ".\n";
	}

	if (overhang != 0) {

		init_host_memory(batch_size, true);

		for (int remainder = 0; remainder < overhang; ++remainder) {

			memcpy(&host_reads[remainder * readLength],
					reads[batch_num * batch_size + remainder],
					sizeof(char) * readLength);
			std::cout << "Read: " << reads[batch_num * batch_size + remainder]
					<< std::endl;
			memcpy(&host_refs[remainder * refLength],
					refs[batch_num * batch_size + remainder],
					sizeof(char) * refLength);
			std::cout << "Ref: " << refs[batch_num * batch_size + remainder]
					<< std::endl;
		}

		cl::Buffer read_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				sizeof(char) * batch_size * readLength, host_reads);
		cl::Buffer ref_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				sizeof(char) * batch_size * refLength, host_refs);
		cl::Buffer result_buffer(context,
		CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(short) * batch_size,
				host_scores);

		std::cout << "Running \"" << kernel_name << "\" Kernel...\n";

		kernel.setArg(0, read_buffer);
		kernel.setArg(1, ref_buffer);
		kernel.setArg(2, result_buffer);

		int workers = batch_size / VECTORS_PER_WORKITEM;

		std::cout << "Work groups: " << workers << std::endl;

		queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(workers),
				cl::NullRange);
		queue.finish();

		std::cout << "Finished overhang of " << overhang << ".\n";
	}
}

void OpenCLKernel::compute_alignments(int const & opt, int const & aln_number,
		char const * const * const reads, char const * const * const refs,
		Alignment * const alignments) {

	int alignment_algorithm = opt & 0xF;

	char const * kernel_name = "";

	switch (alignment_algorithm) {
	case 0:
		kernel_name = align_smith_waterman_kernel;
		break;
	case 1:
		kernel_name = align_needleman_wunsch_kernel;
		break;
	default:
		// Unsupported mode
		break;
	}

	kernel = setup_kernel(program, kernel_name);

	size_t batch_size = calculate_batch_size_from_memory(kernel, device, false);

	size_t batch_num;
	size_t overhang;

	partition_load(aln_number, batch_size, batch_num, overhang);

	for (int batch = 0; batch < batch_num; ++batch) {

		init_host_memory(batch_size, false);

		for (int i = 0; i < batch_size; ++i) {
			memcpy(&host_reads[i * readLength], reads[batch * batch_size + i],
					sizeof(char) * readLength);
			std::cout << "Read: " << reads[batch * batch_size + i] << std::endl;
			memcpy(&host_refs[i * refLength], refs[batch * batch_size + i],
					sizeof(char) * refLength);
			std::cout << "Ref: " << refs[batch * batch_size + i] << std::endl;
		}

		cl::Buffer read_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				sizeof(char) * batch_size * readLength, host_reads);
		cl::Buffer ref_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				sizeof(char) * batch_size * refLength, host_refs);
		cl::Buffer result_buffer(context,
				CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
				sizeof(char) * alnLength * 2 * batch_size, host_alignments);
		cl::Buffer index_buffer(context,
				CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
				sizeof(short) * batch_size * 2, host_indices);
		cl::Buffer matrix_buffer(context,
				CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
				sizeof(short) * matrix_size * batch_size, host_matrix);

		std::cout << "Running \"" << kernel_name << "\" Kernel...\n";

		kernel.setArg(0, read_buffer);
		kernel.setArg(1, ref_buffer);
		kernel.setArg(2, result_buffer);
		kernel.setArg(3, index_buffer);
		kernel.setArg(4, matrix_buffer);

		int workers = batch_size / VECTORS_PER_WORKITEM;

		std::cout << "Work groups: " << workers << std::endl;

		queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(workers),
				cl::NullRange);
		queue.finish();

		std::cout << "Finished batch " << batch << ".\n";
	}

	if (overhang != 0) {

		init_host_memory(batch_size, false);

		for (int remainder = 0; remainder < overhang; ++remainder) {

			memcpy(&host_reads[remainder * readLength],
					reads[batch_num * batch_size + remainder],
					sizeof(char) * readLength);
			std::cout << "Read: " << reads[batch_num * batch_size + remainder]
					<< std::endl;
			memcpy(&host_refs[remainder * refLength],
					refs[batch_num * batch_size + remainder],
					sizeof(char) * refLength);
			std::cout << "Ref: " << refs[batch_num * batch_size + remainder]
					<< std::endl;
		}

		cl::Buffer read_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				sizeof(char) * batch_size * readLength, host_reads);
		cl::Buffer ref_buffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				sizeof(char) * batch_size * refLength, host_refs);
		cl::Buffer result_buffer(context,
				CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
				sizeof(char) * alnLength * 2 * batch_size, host_alignments);
		cl::Buffer index_buffer(context,
				CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
				sizeof(short) * batch_size * 2, host_indices);
		cl::Buffer matrix_buffer(context,
				CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
				sizeof(short) * matrix_size * batch_size, host_matrix);

		std::cout << "Running \"" << kernel_name << "\" Kernel...\n";

		kernel.setArg(0, read_buffer);
		kernel.setArg(1, ref_buffer);
		kernel.setArg(2, result_buffer);
		kernel.setArg(3, index_buffer);
		kernel.setArg(4, matrix_buffer);

		int workers = batch_size / VECTORS_PER_WORKITEM;

		std::cout << "Work groups: " << workers << std::endl;

		queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(workers),
				cl::NullRange);
		queue.finish();

		std::cout << "Finished overhang of " << overhang << ".\n";
	}
}

void OpenCLKernel::initialize_opencl_environment() {

	cl::Device cpu_device = setup_opencl_device(CL_DEVICE_TYPE_CPU);

	std::vector<cl::Device> sub_devices = fission_opencl_device(cpu_device);

	if (sub_devices.size() == 0) {
		throw "Cannot partition OpenCL CPU device.";
	}

	device = sub_devices[0];
	context = setup_context(device);
	std::cout << cpu_device.getInfo<CL_DEVICE_NAME>() << "\n\n";
	char a;
	std::cin >> a;
	program = setup_program(context);

	queue = setup_queue(context, device);
}

std::string get_device_type(cl_int const & device_type) {
	return device_type == CL_DEVICE_TYPE_GPU ? "GPU" : "CPU";
}

void print_device_info(cl::Device const & device) {
	std::cout << device.getInfo<CL_DEVICE_NAME>() << "\n\n";
	std::cout << "\tType: " << get_device_type(device.getInfo<CL_DEVICE_TYPE>())
			<< std::endl;
	std::cout << "\tVendor: " << device.getInfo<CL_DEVICE_VENDOR>()
			<< std::endl;
	std::cout << "\tMax Compute Units: "
			<< device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << std::endl;
	std::cout << "\tGlobal Memory: "
			<< device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() << std::endl;
	std::cout << "\tMax Clock Frequency: "
			<< device.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << std::endl;
	std::cout << "\tMax Allocateable Memory: "
			<< device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() << std::endl;
	std::cout << "\tLocal Memory: "
			<< device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << std::endl;
	std::cout << "\tAvailable: " << device.getInfo< CL_DEVICE_AVAILABLE>()
			<< std::endl;
}

void present_devices(std::vector<cl::Device> & all_devices) {
	std::cout << "Available devices:" << std::endl;

	for (std::vector<cl::Device>::iterator i = all_devices.begin();
			i != all_devices.end(); ++i) {
		print_device_info(*i);
	}
}

cl::Context OpenCLKernel::setup_context(cl::Device const & device) {
	return cl::Context(device);
}

cl::Program OpenCLKernel::setup_program(cl::Context const & context) {

	cl_int cl_error_num = CL_SUCCESS;

	cl::Program::Sources sources;

	// Load Kernel
	std::stringstream source_loader;
	source_loader << opencl_definitions << scoring_kernels << alignment_kernels;

	std::string source(source_loader.str());

	sources.push_back(std::make_pair(source.data(), source.length()));

	stringstream compilerDefines;
	compilerDefines << "-D read_length=" << this->readLength
			<< " -D ref_length=" << this->refLength << " -D aln_length="
			<< this->readLength + this->refLength << " -D score_gap_read="
			<< this->scoreGapRead << " -D score_gap_ref=" << this->scoreGapRef
			<< " -D score_match=" << this->scoreMatch << " -D score_mismatch="
			<< this->scoreMismatch;

	cl::Program program(context, sources);

	std::cout << "Only thing left: Building the program." << std::endl;
	char a;
	std::cin >> a;

	cl_error_num = program.build(compilerDefines.str().c_str());
	check_opencl_success("Error building OpenCL program: ", cl_error_num);

	std::cout << "Done: Building the program." << std::endl;
	std::cin >> a;

	return program;
}

cl::CommandQueue OpenCLKernel::setup_queue(cl::Context const & context,
		cl::Device const & device) {
	return cl::CommandQueue(context, device);
}

cl::Device OpenCLKernel::setup_opencl_device(
		cl_device_type const & device_type) {

	cl_int cl_error_num = CL_SUCCESS;

	//get all platforms (drivers)
	std::vector<cl::Platform> all_platforms;
	cl_error_num = cl::Platform::get(&all_platforms);
	check_opencl_success("OpenCL platform query failed: ", cl_error_num);

	if (all_platforms.size() == 0) {
		std::cout << " No platforms found. Check OpenCL installation!\n";
		exit(-1);
	}

	cl::Platform default_platform = all_platforms[0];
	std::cout << "Using platform: "
			<< default_platform.getInfo<CL_PLATFORM_NAME>() << "\n";

	//get default device of the default platform
	std::vector<cl::Device> all_devices;
	cl_error_num = default_platform.getDevices(CL_DEVICE_TYPE_ALL,
			&all_devices);
	check_opencl_success("OpenCL device query failed: ", cl_error_num);

	present_devices(all_devices);

	std::vector<cl::Device> cpu_devices;
	cl_error_num = default_platform.getDevices(device_type, &cpu_devices);
	check_opencl_success("OpenCL device query failed: ", cl_error_num);

	if (cpu_devices.size() == 0) {
		std::cout << "No devices found. Check OpenCL installation!\n";
		exit(-1);
	}

	cl::Device cpu_device = cpu_devices[0];
	std::cout << "Using device: " << cpu_device.getInfo<CL_DEVICE_NAME>()
			<< std::endl;

	return cpu_device;

}

std::vector<cl::Device> OpenCLKernel::fission_opencl_device(
		cl::Device & device) {

	std::vector<cl::Device> subdevices;

	if (device.getInfo<CL_DEVICE_EXTENSIONS>().find("cl_ext_device_fission")
			== std::string::npos) {
		std::cout << "No device fission support! Returning entire device."
				<< std::endl;
		subdevices.push_back(device);
		return subdevices;
	}

	// Calculate number of used cores/hyperthreading threads
	// Always only use 3 quarters for computation and leave 1 quarter for background processes
	int max_devices = device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
	int fission = max_devices / 4 * 3;

	// Partition CPU
	cl_device_partition_property props[4];
	props[0] = CL_DEVICE_PARTITION_BY_COUNTS;
	//props[1] = std::max(fission, 1);
	props[1] = 1;
	props[2] = CL_DEVICE_PARTITION_BY_COUNTS_LIST_END;
	props[3] = 0;

	device.createSubDevices(props, &subdevices);

	std::cout << "Created " << subdevices.size() << " subdevices.\n";

	return subdevices;
}

cl::Kernel OpenCLKernel::setup_kernel(cl::Program const & program,
		char const * kernel_name) {
	return cl::Kernel(program, kernel_name);
}

size_t OpenCLKernel::calculate_batch_size_from_memory(cl::Kernel const & kernel,
		cl::Device const & device, bool const & score) {

	size_t max_alloc_mem = device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>();
	size_t max_work_items = kernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(
			device);
	std::cout << "Kernel group size: " << max_work_items << std::endl;

	size_t _one_alignment;

	if (score) {
		_one_alignment = sizeof(char) * readLength + sizeof(char) * refLength
				+ sizeof(short) * 2;
	} else {
		_one_alignment = sizeof(char) * readLength + sizeof(char) * refLength
				+ sizeof(char) * alnLength * 2 + sizeof(short) * 2
				+ sizeof(short) * matrix_size;
	}

	size_t const _1MB = 1024 * 1024;

	std::cout << "1MB in bytes:\t" << _1MB << std::endl;
	std::cout << "Sizeof Alignment in bytes:\t" << _one_alignment << std::endl;

	size_t batch_size =
			((_1MB / _one_alignment) / VECTORS_PER_WORKITEM
					* VECTORS_PER_WORKITEM) == 0 ?
					VECTORS_PER_WORKITEM :
					(_1MB / _one_alignment) / VECTORS_PER_WORKITEM
							* VECTORS_PER_WORKITEM;

	batch_size = std::min(batch_size, max_work_items * VECTORS_PER_WORKITEM);

	std::cout << "Batch size:\t" << batch_size << std::endl;

	return batch_size;
}

void OpenCLKernel::partition_load(int const & aln_number,
		size_t const & batch_size, size_t & batch_num, size_t & overhang) {
	batch_num = aln_number / batch_size;
	overhang = aln_number & batch_size;

	std::cout << "Num batches:\t" << batch_num << std::endl;
	std::cout << "Overhang:\t" << overhang << std::endl;
}

void OpenCLKernel::init_host_memory(size_t const & batch_size,
		bool const & score) {

	if (host_reads == 0)
		host_reads = new char[batch_size * readLength];
	if (host_refs == 0)
		host_refs = new char[batch_size * refLength];
	if (score && host_scores == 0)
		host_scores = new short[batch_size];
	if (!score && host_alignments == 0)
		host_alignments = new char[alnLength * 2 * batch_size];
	if (!score && host_indices == 0)
		host_indices = new short[2 * batch_size];
	if (!score && host_matrix == 0)
		host_matrix = new short[matrix_size * batch_size];

	memset(host_reads, 0, sizeof(char) * batch_size * readLength);
	memset(host_refs, 0, sizeof(char) * batch_size * refLength);

	if (score) {
		memset(host_scores, 0, sizeof(short) * batch_size);
	} else {
		memset(host_alignments, 0, sizeof(char) * batch_size * alnLength * 2);
		memset(host_indices, 0, sizeof(short) * batch_size * 2);
		memset(host_matrix, 0, sizeof(short) * matrix_size * batch_size);
	}
}

void OpenCLKernel::check_opencl_success(char const * msg, cl_int ci_error_num) {
	if (ci_error_num != CL_SUCCESS) {
		std::cout << cl_error_to_string(ci_error_num);
		throw;
	}
}

char const * OpenCLKernel::cl_error_to_string(cl_int ci_error_num) {
	switch (ci_error_num) {
	case CL_SUCCESS:
		return strdup("Success.");
	case CL_DEVICE_NOT_FOUND:
		return strdup("Device not found.");
	case CL_DEVICE_NOT_AVAILABLE:
		return strdup("Device unavailable");
	case CL_COMPILER_NOT_AVAILABLE:
		return strdup("Compiler unavailable");
	case CL_MEM_OBJECT_ALLOCATION_FAILURE:
		return strdup("Memory object allocation failure");
	case CL_OUT_OF_RESOURCES:
		return strdup("Out of resources");
	case CL_OUT_OF_HOST_MEMORY:
		return strdup("Out of host memory");
	case CL_PROFILING_INFO_NOT_AVAILABLE:
		return strdup("Profiling information not available");
	case CL_MEM_COPY_OVERLAP:
		return strdup("Memory copy overlap");
	case CL_IMAGE_FORMAT_MISMATCH:
		return strdup("Image format mismatch");
	case CL_IMAGE_FORMAT_NOT_SUPPORTED:
		return strdup("Image format not supported");
	case CL_BUILD_PROGRAM_FAILURE:
		return strdup("Program build failure");
	case CL_MAP_FAILURE:
		return strdup("Map failure");
	case CL_INVALID_VALUE:
		return strdup("Invalid value");
	case CL_INVALID_DEVICE_TYPE:
		return strdup("Invalid device type");
	case CL_INVALID_PLATFORM:
		return strdup("Invalid platform");
	case CL_INVALID_DEVICE:
		return strdup("Invalid device");
	case CL_INVALID_CONTEXT:
		return strdup("Invalid context");
	case CL_INVALID_QUEUE_PROPERTIES:
		return strdup("Invalid queue properties");
	case CL_INVALID_COMMAND_QUEUE:
		return strdup("Invalid command queue");
	case CL_INVALID_HOST_PTR:
		return strdup("Invalid host pointer");
	case CL_INVALID_MEM_OBJECT:
		return strdup("Invalid memory object");
	case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
		return strdup("Invalid image format descriptor");
	case CL_INVALID_IMAGE_SIZE:
		return strdup("Invalid image size");
	case CL_INVALID_SAMPLER:
		return strdup("Invalid sampler");
	case CL_INVALID_BINARY:
		return strdup("Invalid binary");
	case CL_INVALID_BUILD_OPTIONS:
		return strdup("Invalid build options");
	case CL_INVALID_PROGRAM:
		return strdup("Invalid program");
	case CL_INVALID_PROGRAM_EXECUTABLE:
		return strdup("Invalid program executable");
	case CL_INVALID_KERNEL_NAME:
		return strdup("Invalid kernel name");
	case CL_INVALID_KERNEL_DEFINITION:
		return strdup("Invalid kernel definition");
	case CL_INVALID_KERNEL:
		return strdup("Invalid kernel");
	case CL_INVALID_ARG_INDEX:
		return strdup("Invalid argument index");
	case CL_INVALID_ARG_VALUE:
		return strdup("Invalid argument value");
	case CL_INVALID_ARG_SIZE:
		return strdup("Invalid argument size");
	case CL_INVALID_KERNEL_ARGS:
		return strdup("Invalid kernel arguments");
	case CL_INVALID_WORK_DIMENSION:
		return strdup("Invalid work dimension");
	case CL_INVALID_WORK_GROUP_SIZE:
		return strdup("Invalid work group size");
	case CL_INVALID_WORK_ITEM_SIZE:
		return strdup("Invalid work item size");
	case CL_INVALID_GLOBAL_OFFSET:
		return strdup("Invalid global offset");
	case CL_INVALID_EVENT_WAIT_LIST:
		return strdup("Invalid event wait list");
	case CL_INVALID_EVENT:
		return strdup("Invalid event");
	case CL_INVALID_OPERATION:
		return strdup("Invalid operation");
	case CL_INVALID_GL_OBJECT:
		return strdup("Invalid OpenGL object");
	case CL_INVALID_BUFFER_SIZE:
		return strdup("Invalid buffer size");
	case CL_INVALID_MIP_LEVEL:
		return strdup("Invalid mip-map level");
	default:
		return strdup("Unknown");
	}
}
