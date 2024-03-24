# DDT

This is the codebase of DDT, which is based on the libgrape-lite.

### Dependencies

- [CMake](https://cmake.org/) (>=2.8)
- A modern C++ compiler compliant with C++-11 standard. (g++ >= 4.8.1 or clang++ >= 3.3)
- [MPICH](https://www.mpich.org/) (>= 2.1.4) or [OpenMPI](https://www.open-mpi.org/) (>= 3.0.0)
- [glog](https://github.com/google/glog) (>= 0.3.4)

### Usage

- ifconfig: `docker0, eno1, lo, veth0a86a69`

- `mpirun -n 96 --hostfile HOSTFILE --mca btl_base_warn_component_unused 0 --mca btl_tcp_if_include eno1 ./run_app --application=ddt --vfile vertex_file_path --efile edge_file_path --out_prefix ./output_ddt`
