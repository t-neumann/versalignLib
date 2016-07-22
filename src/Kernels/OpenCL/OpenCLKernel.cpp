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

void present_devices(std::vector<cl::Device> const & all_devices) {
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

	cl_error_num = program.build(compilerDefines.str().c_str());
	check_opencl_success("Error building OpenCL program: ", cl_error_num);

	return program;
}

cl::CommandQueue OpenCLKernel::setup_queue(cl::Context const & context, cl::Device const & device) {
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
	cl_error_num = default_platform.getDevices(device_type, &all_devices);
	check_opencl_success("OpenCL device query failed: ", cl_error_num);

	if (cpu_devices.size() == 0) {
		std::cout << "No devices found. Check OpenCL installation!\n";
		exit(-1);
	}

	cl::Device cpu_device = cpu_devices[0];
	std::cout << "Using device: " << cpu_device.getInfo<CL_DEVICE_NAME>()
			<< std::endl;

//	std::vector<size_t> work_item_sizes = default_device.getInfo<
//	CL_DEVICE_MAX_WORK_ITEM_SIZES>();
//	for (std::vector<size_t>::iterator itor = work_item_sizes.begin();
//			itor != work_item_sizes.end(); ++itor) {
//		std::cout << "Available work item sizes: " << *itor << std::endl;
//	}

	std::cout << "Build in kernels: "
			<< cpu_device.getInfo<CL_DEVICE_BUILT_IN_KERNELS>() << std::endl;

	return cpu_device;

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
