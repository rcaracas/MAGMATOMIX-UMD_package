"""
Created on Fri Sep  4 17:14:42 2020

@author: Laurent Gilquin
"""

from os.path import dirname, abspath, join
import warnings
import numpy as np
import pyopencl as cl
from itertools import chain


def closest_two_pow(n, x):
    """
    Find the closest power of two higher than x.
    """

    if (2**n) > x:
        return n
    else:
        return closest_two_pow(n+1,x)



def get_flags(shape, global_mem, local_mem, extensions):
    """    
    Function that checks which version should be used to calculate the distance matrix
    according to the device global and local memory spaces as well as the input/output
    sizes.
    """
    # initiliaze dict flags
    dict_flags = {}
    
    # function returning the total size to allocate in the gpu global memory
    def global_bits(shape, float_t, int_t):
        
        bits = np.prod(shape) * float_t # coordinate matrix
        bits += shape[0] * (shape[0] - 1) // 2 * float_t # distance matrix
        bits += int(4 * float_t + 6 * int_t) # other variables
        return bits

    # check the integer type
    if (shape[0] * (shape[0] - 1) //2) < 2**16:
        dict_flags["INT"] = "ushort"
        int_t = 16
    else:
        dict_flags["INT"] = "int"
        int_t = 32

    #TODO: correct ugly repetition of raise BufferError
    # check the floating point type
    float_types = [64, 32, 16]
    sizes = np.array(list(map(lambda x: global_bits(shape, x, int_t), float_types)))
    choice = np.where((global_mem - sizes) > 1000)[0]
    if len(choice) != 0:
        if choice[0] == 0:
            dict_flags["HALF"] = False
            if "cl_khr_fp64" in extensions:
                dict_flags["REAL"] = "double"
            else:
                dict_flags["REAL"] = "float"
        elif choice[0] == 1:
            dict_flags["HALF"] = False
            dict_flags["REAL"] = "float"
        elif choice[0] == 2:
            if "cl_khr_fp16" in extensions:
                dict_flags["HALF"] = True
                dict_flags["REAL"] = "float"
            else:
                raise BufferError(
                """
                Not enough space in the global memory, the distance calculation
                will default to use the C extension.\n
                """
                )
    else:
        raise BufferError(
        """
        Not enough space in the global memory, the distance calculation
        will default to use the C extension.\n
        """
        )

    # check if tiling is possible and which floating type
    sizes = np.array(list(map(lambda x: int(2*16*3*x), float_types[:2])))
    choice = np.where((local_mem/sizes) > 0.5)[0]
    if len(choice) == 0:
        dict_flags["TILED"] = False
        dict_flags["REAL_IN"] = "float" # needs to be default to build whole program
    else:
        dict_flags["TILED"] = True
        if (choice[0] == 0) and ("cl_khr_fp64" in extensions):
            dict_flags["REAL_IN"] = "double"
        else:
            dict_flags["REAL_IN"] = "float"

    return dict_flags
    


def get_gpu(plateforms=None, device=None):
    """
    Utility function that retrieves the "best" gpu, in terms of memory space,
    from the platform.
    """

    # get GPUs from list of devices
    platforms = cl.get_platforms()
    devices = list(chain.from_iterable([[device for device in pl.get_devices()]\
                                         for pl in platforms]))
    gpus = [device for device in devices if device.type == cl.device_type.GPU]
    
    if len(gpus) != 0:
        # select best GPU
        global_mems = [gpu.global_mem_size for gpu in gpus]
        best_global = global_mems.index(max(global_mems))
        gpu = gpus[best_global] # only one on my laptop
        return gpu
    else:
        warnings.warn("""
        Pyopencl found no GPU, the distance calculation will default to use the C extension.\n
        Please consult the installation page of pyopencl to make sure that all
        additional softwares have been installed:
        \thttps://documen.tician.de/pyopencl/misc.html\n
        Additionally, you should check that your GPU is compatible with OpenCL:
        \thttps://www.khronos.org/conformance/adopters/conformant-products/opencl\n
        """, ResourceWarning)
        return None



def gpumd_pdist(X, coeffs, gpu):
    """
    Utility function that calculates the umd distance on a gpu.
    """

    # create opencl context
    context = cl.Context([gpu])
    queue = cl.CommandQueue(context, gpu)
    # get the cl file
    path_cl = abspath(join(dirname(__file__), 'gpu_utils.cl'))
    kernelsource = open(path_cl).read()
    # build the program with the according flags
    dict_flags = get_flags(X.shape,gpu.global_mem_size, gpu.local_mem_size,
                           gpu.extensions.split(" "))
    TS = int(np.sqrt(gpu.max_work_group_size))
    flags = ["INT={}".format(dict_flags["INT"]),
             "REAL={}".format(dict_flags["REAL"]),
             "REAL_IN={}".format(dict_flags["REAL_IN"]),
             "TS={}".format(TS), "K=3"]
    flags = ["-D"+ele for ele in flags]
    flags += ["-cl-fast-relaxed-math"]
    program = cl.Program(context, kernelsource).build(options=flags)
    
    # initialize buffers
    dim = X.shape[0] * (X.shape[0] - 1) // 2
    if dict_flags["HALF"]:
        hR_out = np.empty(dim, dtype=cl.cltypes.half)
        dX_in = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                          hostbuf=X.astype(cl.cltypes.half))
        dC_in = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                          hostbuf=coeffs.astype(cl.cltypes.half))
    else:
        cl_flags = {"int": cl.cltypes.int,
                    "ushort": cl.cltypes.ushort,
                    "double": cl.cltypes.double,
                    "float": cl.cltypes.float}
        hR_out = np.empty(dim, dtype=cl_flags[dict_flags["REAL"]])
        dX_in = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                          hostbuf=X.astype(cl_flags[dict_flags["REAL"]]))
        dC_in = cl.Buffer(context, cl.mem_flags.READ_ONLY | cl.mem_flags.COPY_HOST_PTR,
                          hostbuf=coeffs.astype(cl_flags[dict_flags["REAL"]]))
    # create output array in device memory
    dR_out = cl.Buffer(context, cl.mem_flags.WRITE_ONLY, hR_out.nbytes)
    
    # select the kernel to use
    if dict_flags["HALF"]:
        if dict_flags["TILED"]:
            kernel = program.half_tiled_pdist
        else:
            kernel = program.half_pdist
    else:
        if dict_flags["TILED"]:
            kernel = program.tiled_pdist
        else:
            kernel = program.pdist

    # set kernel arguments
    kernel.set_scalar_arg_dtypes([cl_flags[dict_flags["INT"]], None, None, None])
    global_grid = tuple([int(2**closest_two_pow(1, X.shape[0]))]*2)
    if global_grid[0] > 16:
        local_grid = tuple([TS]*2)
    else:
        local_grid = None

    # call the kernel
    event = kernel(queue, global_grid, local_grid, X.shape[0], dX_in, dR_out, dC_in)
    event.wait()
    cl.enqueue_copy(queue, hR_out, dR_out)

    return hR_out
