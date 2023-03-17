#include "xcl2.hpp"
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "wsat_hls.h"

void init(std::vector<int, aligned_allocator<int>>& ClausesVec, std::vector<int, aligned_allocator<int>>& VarsLocVec, 
    std::vector<short, aligned_allocator<short>>& clsK_v, std::vector<short, aligned_allocator<short>>& varND_v,
    std::string fileName, int& numVars, int& numClauses) {

    int VarsLocVec2[MAXNLIT * R];
    short varR_off[MAXNLIT];

    std::ifstream fileDIMACS(fileName);

	int maxk = 0;
	int maxr = 0;

	if(fileDIMACS.is_open()){
		std::string line;
		char pp;
		while ((pp = fileDIMACS.peek()) == 'c') {
			fileDIMACS.ignore(256, '\n');
		}

		// parsing the first line
		fileDIMACS >> line;	// p
		fileDIMACS >> line;	// cnf
		fileDIMACS >> numVars;	// nv
		fileDIMACS >> numClauses;	// nc

		int nextLit;
		int clauseCounted = 1;

		char p;
		int nvinc = 0;	// #var in a clause
		while (!fileDIMACS.eof()) {
			fileDIMACS.ignore(256, '\n');
			p = fileDIMACS.peek();
			if (p == EOF || fileDIMACS.eof()) {
				break;
			}
			if (p == 'c') {
				fileDIMACS.ignore(256, '\n');
				continue;
			}
			if ((p > '0' && p <= '9') || p == '-' || p == ' ') {
				if (p == ' ') {
					fileDIMACS.ignore(256, ' ');
				}

				// inside clause line
				while (true) {
					fileDIMACS >> nextLit;

					int idx = clauseCounted * K + nvinc ;
					ClausesVec[idx] = nextLit;

					if (nextLit == 0){	// end of line
						if (nvinc > maxk) {
							maxk = nvinc;
						}
						clsK_v[clauseCounted] = nvinc;
						nvinc = 0;
						clauseCounted++;
						break;
					}
					nvinc++;

					for (int i = 0; i < R; i++) {
						if (VarsLocVec2[GETPOS(nextLit) * R + i] == 0) {
							VarsLocVec2[GETPOS(nextLit)* R + i] = clauseCounted;

							varR_off[GETPOS(nextLit)] = i + 1;
							if (i + 1 > maxr) {
								maxr = i + 1;
							}
							break;
						}
						if (i == R - 1) {
							std::cout << "R excessed at " << nextLit << "\n";
						}
					}
				}
			}
		}

		if (clauseCounted - 1 != numClauses) {
			std::cout << "#clauses error || written: " << numClauses << " counted: " << clauseCounted - 1 << "\n";
		}

	} else {
		std::cerr << "Cannot open file: " << fileName << "\n";
	}


	for (int v = 0; v < numVars * 2; v++) {

		std::vector<int> waitlist[DSIZE];

		int endR = varR_off[v];
		for (int i = 0; i < endR; i++) {
			int loc = VarsLocVec2[v * R + i];
			waitlist[loc % DSIZE].push_back(loc);
		}

		for (int t = 0; t < R; t++) {	// all lines
			bool finished = true;

			for (int i = 0; i < DSIZE; i++) {	// waitlist one line
				// int idx = (v * DSIZE * R) + (DSIZE * t) + i;
				int idx = (v * R) + t;
				if (not waitlist[i].empty()) {
					finished = false;

					VarsLocVec[idx * DSIZE + i] = waitlist[i][0];
					waitlist[i].erase(waitlist[i].begin());
				}
				else {
					VarsLocVec[idx * DSIZE + i] = 0;
				}
			}
			if (finished) {
				// varND_off[v] = t;
				varND_v[v] = t;
				break;
			}
		}
	}
	std::cout << "numvars: " << numVars << " numcls: " << numClauses << " maxk: " << maxk << " maxr: " << maxr << "\n";
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " <XCLBIN File> <fileName>" << std::endl;
        return EXIT_FAILURE;
    }
    std::string binaryFile = argv[1];
    std::string fileName = "/home/centos/src/project_data/aws-fpga/Vitis/examples/xilinx_2021.2/mm2/src/k3-r4.26-v600-c2556-043.cnf";

    cl_int err;
    cl::Context context;
    cl::Kernel krnl_yalsat;
    cl::CommandQueue q;
    // Allocate Memory in Host Memory
    // When creating a buffer with user pointer (CL_MEM_USE_HOST_PTR), under the
    // hood user ptr
    // is used if it is properly aligned. when not aligned, runtime had no choice
    // but to create
    // its own host side buffer. So it is recommended to use this allocator if
    // user wish to
    // create buffer using CL_MEM_USE_HOST_PTR to align user buffer to page
    // boundary. It will
    // ensure that user buffer is used when user create Buffer/Mem object with
    // CL_MEM_USE_HOST_PTR

	// std::vector<DTYPE, aligned_allocator<DTYPE> > At(SIZE*SIZE); 
	// std::vector<DTYPE, aligned_allocator<DTYPE> > B(SIZE*SIZE); 
	// std::vector<DTYPE, aligned_allocator<DTYPE> > AB_sw(SIZE*SIZE); 
	// std::vector<DTYPE, aligned_allocator<DTYPE> > AB_hw(SIZE*SIZE); 
    
    std::vector<int, aligned_allocator<int>> ClausesVec(MAXNCLS * K);
    std::vector<int, aligned_allocator<int>> VarsLocVec(MAXNLIT * R * DSIZE);
    std::vector<short, aligned_allocator<short>> clsK_v(MAXNCLS);
    std::vector<short, aligned_allocator<short>> varND_v(MAXNLIT);

	std::vector<int, aligned_allocator<int>> retval(1);

    int numVars, numClauses;
    int seed = time(0);
	seed = 1678121930;
    srand(seed);

    std::cout << fileName << " seed: " << seed << "\n";
    init(ClausesVec, VarsLocVec, clsK_v, varND_v, fileName, numVars, numClauses);

	printf("Done initializing vectors\n");

    // OPENCL HOST CODE AREA START
    // get_xil_devices() is a utility API which will find the xilinx
    // platforms and will return list of devices connected to Xilinx platform
    auto devices = xcl::get_xil_devices();
    // read_binary_file() is a utility API which will load the binaryFile
    // and will return the pointer to file buffer.
    auto fileBuf = xcl::read_binary_file(binaryFile);
    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};
    bool valid_device = false;
    for (unsigned int i = 0; i < devices.size(); i++) {
        auto device = devices[i];
        // Creating Context and Command Queue for selected Device
        OCL_CHECK(err, context = cl::Context(device, nullptr, nullptr, nullptr, &err));
        OCL_CHECK(err, q = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err));
        std::cout << "Trying to program device[" << i << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
        cl::Program program(context, {device}, bins, nullptr, &err);
        if (err != CL_SUCCESS) {
            std::cout << "Failed to program device[" << i << "] with xclbin file!\n";
        } else {
            std::cout << "Device[" << i << "]: program successful!\n";
            OCL_CHECK(err, krnl_yalsat = cl::Kernel(program, "mm", &err));
            valid_device = true;
            break; // we break because we found a valid device
        }
    }
    if (!valid_device) {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

    // Allocate Buffer in Global Memory
    // Buffers are allocated using CL_MEM_USE_HOST_PTR for efficient memory and
    // Device-to-host communication
    OCL_CHECK(err, cl::Buffer b_varND_off(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(short)* MAXNLIT, varND_v.data(), &err));
    OCL_CHECK(err, cl::Buffer b_clsK_off(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(short)* MAXNCLS, clsK_v.data(), &err));
    OCL_CHECK(err, cl::Buffer b_CV(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(int)* MAXNCLS * K, ClausesVec.data(), &err));
    OCL_CHECK(err, cl::Buffer b_VLV(context, CL_MEM_USE_HOST_PTR | CL_MEM_READ_ONLY, sizeof(int)* MAXNLIT * R * DSIZE, VarsLocVec.data(), &err));
	
	OCL_CHECK(err, cl::Buffer b_ret(context, CL_MEM_USE_HOST_PTR | CL_MEM_WRITE_ONLY, sizeof(int)*1, retval.data(), &err));

    OCL_CHECK(err, err = krnl_yalsat.setArg(0, b_varND_off));
    OCL_CHECK(err, err = krnl_yalsat.setArg(1, b_clsK_off));
    OCL_CHECK(err, err = krnl_yalsat.setArg(2, b_CV));
    OCL_CHECK(err, err = krnl_yalsat.setArg(3, b_VLV));
    OCL_CHECK(err, err = krnl_yalsat.setArg(4, numVars));
    OCL_CHECK(err, err = krnl_yalsat.setArg(5, numClauses));
    OCL_CHECK(err, err = krnl_yalsat.setArg(6, b_ret));
    // OCL_CHECK(err, err = krnl_yalsat.setArg(7, &issolved));
    OCL_CHECK(err, err = krnl_yalsat.setArg(7, seed));

	OCL_CHECK(err, err = q.enqueueMigrateMemObjects({b_varND_off, b_clsK_off, b_CV, b_VLV}, 0 /* 0 means from host*/));
	q.finish();
	
	std::cout << "Running FPGA...\n";
	auto start = std::chrono::steady_clock::now();

	OCL_CHECK(err, err = q.enqueueTask(krnl_yalsat));
	q.finish();

	auto end = std::chrono::steady_clock::now();
	std::cout << "Done.\n";
	double exec_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
	

	OCL_CHECK(err, err = q.enqueueMigrateMemObjects({b_ret}, CL_MIGRATE_MEM_OBJECT_HOST));
	// OCL_CHECK(err, err = q.enqueueMigrateMemObjects({buffer_AB}, CL_MIGRATE_MEM_OBJECT_HOST));
	q.finish();

	if(retval[0] > 0){
        std::cout << "Solver found a solution. | ";
		std::cout << "Time: " << exec_time*1e-9 << ", flips: " << retval[0] << std::endl;
	}
	else{
		std::cout << "Solver could not find a solution. | ";
		std::cout << "Time: " << exec_time*1e-9 << ", flips: " << MAX_FLIPS << std::endl;
	}

	return EXIT_SUCCESS;
}
