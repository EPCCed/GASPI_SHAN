# Heat Benchmark

## Description
The Heat simulation uses an iterative Gauss-Seidel method to solve the heat equation,
which is a parabolic partial differential equation that describes the distribution of
heat (or variation in temperature) in a given region over time.

The heat equation is of fundamental importance in a wide range of science fields. In
mathematics, it is the parabolic partial differential equation par excellence. In statistics,
it is related to the study of the Brownian motion. Also, the diffusion equation is a generic
version of the heat equation, and it is related to the study of chemical diffusion processes.

## Requirements
The requirements of this application are shown in the following lists. The main requirements are:

  * **GNU Compiler Collection**.
  * **OmpSs-2**: OmpSs-2 is the second generation of the OmpSs programming model. It is a task-based
    programming model originated from the ideas of the OpenMP and StarSs programming models. The
    specification and user-guide are available at https://pm.bsc.es/ompss-2-docs/spec/ and
    https://pm.bsc.es/ompss-2-docs/user-guide/, respectively. OmpSs-2 requires both Mercurium and
    Nanos6 tools. Mercurium is a source-to-source compiler which provides the necessary support for
    transforming the high-level directives into a parallelized version of the application. The Nanos6
    runtime system library provides the services to manage all the parallelism in the application (e.g. task
    creation, synchronization, scheduling, etc). Both can be downloaded from https://github.com/bsc-pm.
  * **MPI**: This application requires an MPI library supporting the multi-threading mode. It mainly targets
    the MPICH implementation, however, it should work with other libraries by adding the needed implementation-specific
    flags in the Makefile.

In addition, there are optional tools which enable the building of other application versions:

  * **Task-Aware MPI (TAMPI)**: The Task-Aware MPI library provides the interoperability mechanism for MPI
    and OmpSs. Please contact <pm-tools@bsc.es> to get access to the TAMPI library.

## Versions

The heat application has several versions which are compiled in different 
binaries, by executing the `make` command. They are:

  * **01.heat_seq.bin**: Sequential version of this benchmark.
  * **02.heat_ompss.bin**: Parallel version with OmpSs-2. It divides the 2-D matrix into different 2-D blocks of consecutive elements, and
    a task is created to compute each block. Tasks declare fine-grained dependencies on the target block and its adjacent blocks.
  * **03.heat_mpi.bin**: Parallel version using MPI. It divides horizontally the 2-D matrix and each MPI process is responsible for computing
    its assigned rows.
  * **04.heat_mpi_ompss_forkjoin.bin**: Parallel version using MPI + OmpSs-2. It uses a fork-join parallelization strategy, where computation phases
    are parallelized with tasks and communication phases are sequential.
  * **05.heat_mpi_ompss_tasks.bin**: Parallel version using MPI + OmpSs-2 tasks. Both computation and communication are parallelized with tasks.
    However, communication tasks are serialized declaring a dependency on a sentinel variable. This is to prevent deadlocks between processes,
    since communications tasks perform blocking MPI calls.
  * **06.heat_mpi_ompss_tasks_interop.bin**: Parallel version using MPI + OmpSs-2 tasks + Interoperability library. This version disables the
    artificial dependency on the sentinel variable, so communication tasks can run totally in parallel. The interoperability library is in charge
    of intercepting the blocking MPI calls to avoid the blocking of the underlying hardware thread.
  * **07.heat_mpi_ompss_tasks_interop_async.bin**: Parallel version using MPI + OmpSs-2 tasks + Interoperability library. This version
    disables the artificial dependency on the sentinel variable, so communication tasks can run totally in parallel. These tasks do not
    call to blocking MPI procedures (MPI_Recv & MPI_Send), but their non-blocking counterparts (MPI_Irecv & MPI_Isend). The resulting
    MPI requests are linked to the release of the dependencies of the communication task by calling to the non-blocking TAMPI_Iwaitall
    function (or TAMPI_Iwait), and then, the task can finish its execution. Once the linked requests complete, the dependencies of the
    task will be released by the Interoperability library services.
  * **08.heat_mpi_nbuffer.bin**: Parallel version using MPI, but it is a more elaborate version than **03.heat_mpi.bin**. It exchanges the
    block boundaries as soon as possible. In addition, it tries to overlap computation and communication phases by using non-blocking MPI
    calls.
  * **09.heat_mpi_nbuffer.bin**: Parallel version using MPI, but it is a more elaborate version than **08.heat_mpi.bin**. Transposed
    row / columns for faster buildup of wavefront.
  * **010.solver_gaspi_shan.cpp**: GASPI_SHAN version, derived from the MPI version. Uses task model approach, where all potential steps
    (send/recv upper/lower rows, solve matrix) are modeled as mutual dependencies. Uses notified communication with the node (SHAN) 
    and across nodes (GASPI) for fine grain execution of tasks. Uses a simple greedy task execution across all rows.    

  The simplest way to compile this package is:

  1. Stay in Heat root directory to recursively build all the versions.
     The Heat MPI + OmpSs-2 tasks + Interoperability library version is
     compiled only if the environment variable `MPI_INTEROP_HOME`
     is set to the Task-Aware MPI (TAMPI) library's installation directory.

  2. Type `make` to compile the selected benchmark's version(s).
     Optionally, you can use a different block size in each dimension
     (BSX and BSY for vertical and horizontal dimensions, respectively)
     when building the benchmark (by default 1024). Type
     `make BSX=MY_BLOCK_SIZE_X BSY=MY_BLOCK_SIZE_Y` in order to change
     this value. If you want the same value in each dimension, type
     `make BSX=MY_BLOCK_SIZE`.

  3. In addition, you can type 'make check' to check the correctness
     of the built versions. By default, the pure MPI version runs with
     4 processes and the hybrid versions run with 2 MPI processes and 2
     hardware threads for each process. You can change these parameters
     when executing `scripts/run-tests.sh` (see the available options
     passing `-h`).

## Execution instructions

The binaries accept several options. The most relevant options are the size 
of the matrix in each dimension with `-s`, and the number of timesteps with 
`-t`. More options can be seen passing the `-h` option. An example of execution
could be:

```
$ mpiexec -n 4 -bind-to hwthread:16 heat_mpi.task.1024bs.exe -t 150 -s 8192
```

in which the application will perform 150 timesteps in 4 MPI processes with 16 
hardware threads in each process (used by the OmpSs-2 runtime). The size of the
matrix in each dimension will be 8192 (8192^2 elements in total), this means
that each process will have 2048 * 8192 elements (16 blocks per process).

