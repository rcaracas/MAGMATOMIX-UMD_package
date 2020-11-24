#!/usr/bin/env python3
"""
Created on Fri Sep  4 17:14:42 2020

@author: Laurent Gilquin
"""

from os.path import dirname, abspath, join
import warnings
import numpy as np
import pyopencl as cl
from itertools import chain



def global_bits(shape, float_t, int_t):
    """
    Function returning the total size to allocate in the gpu global memory.
    """
    bits = np.prod(shape) * float_t # coordinate matrix
    bits += shape[0] * (shape[0] - 1) // 2 * float_t # distance matrix
    bits += 3 * float_t # coefficients
    bits += int(4 * float_t + 6 * int_t) # other variables
    return bits


def closest_two_pow(x, n=0):
    """
    Find the closest power of two higher than x.
    """
    if (2**n) > x:
        return n
    else:
        return closest_two_pow(x, n+1)



class gpu():
    
    def __init__(self, input_shape, custom_flags=None,
                 platforms=None, devices=None):
        # instantiate arguments
        self.build_flags = {}
        self.context = None
        self.device = devices
        self.input_shape = input_shape
        self.kernel = None
        self.memory_pools = None
        self.platforms = platforms
        self.program = None
        # run at start
        self.get()
        if self.device is not None:
            self.set_flags()
            if self.device is not None:
                self.set_custom_flags(custom_flags)
                self.build_program()
                self.create_memory_pools()
                self.get_kernel()


    def get(self):
        """
        Utility function that retrieves the "best" gpu, in terms of memory space,
        from the detected platforms.
        """
        # get GPUs from list of devices
        if self.platforms is None:
            self.platforms = cl.get_platforms()
        if self.device is None:
            self.device = list(chain.from_iterable([[device for device in pl.get_devices()]\
                                                     for pl in self.platforms]))
        gpus = [device for device in self.device if device.type == cl.device_type.GPU]
        
        if len(gpus) != 0:
            # select best GPU
            global_mems = [gpu.global_mem_size for gpu in gpus]
            best_global = global_mems.index(max(global_mems))
            self.device = gpus[best_global]
        else:
            warnings.warn("""
            Pyopencl found no GPU, the distance calculation will default to use the C extension.\n
            Please consult the installation page of pyopencl to make sure that all
            additional softwares have been installed:
            \thttps://documen.tician.de/pyopencl/misc.html\n
            Additionally, you should check that your GPU is compatible with OpenCL:
            \thttps://www.khronos.org/conformance/adopters/conformant-products/opencl\n
            """, ResourceWarning)



    def set_flags(self):
        """    
        Function that checks which version should be used to calculate the distance matrix
        according to the device global and local memory spaces as well as the input/output
        sizes.
        """
        # check the integer type to use
        if (self.input_shape[0] * (self.input_shape[0] - 1) //2) < 2**16:
            self.build_flags["FINT"] = "ushort"
            int_t = 16
        else:
            self.build_flags["FINT"] = "int"
            int_t = 32
        # check the floating point type
        self.build_flags["FHALF"] = False
        float_types = [64, 32, 16]
        sizes = np.array(list(map(lambda x: global_bits(self.input_shape, x, int_t),
                                  float_types)))
        choice = (self.device.global_mem_size - sizes) > 1000
        if choice[0] and "cl_khr_fp64" in self.device.extensions:
            self.build_flags["FREAL"] = "double"
        elif choice[1]:
            self.build_flags["FREAL"] = "float"
        elif choice[2] and "cl_khr_fp16" in self.device.extensions:
            self.build_flags["FHALF"] = True
            self.build_flags["FREAL"] = "float"
        else:
            warnings.warn("""
            Not enough space in the global memory, the distance calculation
            will default to use the C extension.\n
            """, ResourceWarning)
            self.device = None



    def set_custom_flags(self, custom_flags):
        """
        Check if user flags can be used.
        """
        if custom_flags is not None:
            for key, value in custom_flags.items():
                if key == "FINT":
                    if value in ["ushort", "int"]:
                        self.build_flags[key] = value
                if key == "FREAL":
                    if (value == "float") and self.build_flags[key] == "double":
                        self.build_flags[key] = value
                if key == "FHALF":
                    if value and not self.build_flags[key]:
                        self.build_flags[key] = value



    def build_program(self):
        """
        Function that builds the program from the kernels source file.
        """
        # create opencl context
        self.context = cl.Context([self.device])
        self.queue = cl.CommandQueue(self.context, self.device)
        # get the cl file
        path_cl = abspath(join(dirname(__file__), 'gpu_utils.cl'))
        kernelfile = open(path_cl, 'r')
        kernelsource = kernelfile.read()
        kernelfile.close()
        # build the program
        
        flags = ["FINT={}".format(self.build_flags["FINT"]),
                 "FREAL={}".format(self.build_flags["FREAL"]),
                 "K=3"]
        flags = ["-D"+ele for ele in flags]
        flags += ["-cl-fast-relaxed-math",
                  "-cl-denorms-are-zero",
                  "-cl-finite-math-only"]
        self.program = cl.Program(self.context, kernelsource).build(options=flags)



    def create_memory_pools(self):
        """
        Create memory pool for instantiating inputs and output buffers.
        """
        # create inputs memory pool with read only flag
        input_mem_pool = cl.tools.MemoryPool(\
                cl.tools.ImmediateAllocator(self.queue, cl.mem_flags.READ_ONLY))
        # create inputs memory pool with read only flag
        output_mem_pool = cl.tools.MemoryPool(\
                cl.tools.ImmediateAllocator(self.queue, cl.mem_flags.WRITE_ONLY))
        # store the two pools
        self.memory_pools = {"inputs": input_mem_pool, "output": output_mem_pool}
        # store the float type and the size for the buffers allocation
        cl_flags = {"double": cl.cltypes.double, "float": cl.cltypes.float}
        if self.build_flags["FHALF"]:
            self.float_t = cl.cltypes.half
        else:
            self.float_t = cl_flags[self.build_flags["FREAL"]]
        self.float_size = 8*self.float_t().itemsize



    def free_memory_pools(self):
        """
        Free memory hold by the pools.
        """
        for pool in self.memory_pools.values():
            pool.stop_holding()



    def get_kernel(self):
        """
        Select the kernel to use according to the build flags.
        """
        # select the kernel to use
        if self.build_flags["FHALF"]:
            self.kernel = self.program.half_pdist
        else:
            self.kernel = self.program.pdist
        # set kernel arguments
        cl_flags= {"int": cl.cltypes.int, "ushort": cl.cltypes.ushort}
        self.kernel.set_scalar_arg_dtypes([cl_flags[self.build_flags["FINT"]],
                                           None, None, None])
        self.global_grid = tuple([int(2**closest_two_pow(self.input_shape[0], 1))]*2)
        if self.global_grid[0] > 16:
            self.local_grid = tuple([int(np.sqrt(self.device.max_work_group_size))]*2)
        else:
            self.local_grid = None



    def compute_distance(self, X, coeffs):
        """
        Allocate buffers from memory pools, enqueue the values to device and
        compute the distance.
        """
        # allocate buffer for X (matrix of atoms coordinates) and coeffs
        X_dev = self.memory_pools["inputs"].allocate(self.float_size*X.size)
        C_dev = self.memory_pools["inputs"].allocate(self.float_size*coeffs.size)
        cl.enqueue_copy(self.queue, X_dev, X.astype(self.float_t))
        cl.enqueue_copy(self.queue, C_dev, coeffs.astype(self.float_t))
        # allocate buffer for output
        dim = X.shape[0] * (X.shape[0] - 1) // 2
        res = np.empty(dim, dtype=self.float_t)
        res_dev = self.memory_pools["output"].allocate(self.float_size*res.size)
        # call the kernel
        event = self.kernel(self.queue, self.global_grid, self.local_grid,
                            self.input_shape[0], X_dev, res_dev, C_dev)
        event.wait()
        cl.enqueue_copy(self.queue, res, res_dev)
        # release buffers memory into the pools
        X_dev.release()
        C_dev.release()
        res_dev.release()

        return res
