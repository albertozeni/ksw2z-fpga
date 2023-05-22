# KSW2z - Genomic Alignment Acceleration

## Project Description
<p align="justify">
Pairwise sequence alignment is a fundamental step for many genomics and molecular biology applications.
Given the quadratic time complexity of alignment algorithms, the community demands innovative, fast, and efficient techniques to perform this task.
Furthermore, general-purpose architectures lack the necessary performance to address the computational load of these algorithms.
In this context, we present the first open-source FPGA implementation of the popular KSW2z algorithm employed by minimap2. 
Our design also implements the $Z$-drop heuristic and banded alignment as the original software to further reduce the processing time if needed.
The proposed multi-core accelerator achieves up to 7.70x improvement in speedup and 20.07x in energy efficiency compared to the multi-threaded software implementation run on a Xeon Platinum 8167M processor. 
</p>

## Usage & Compilation

First clone the repository and pull its submodules as follows:
```
git clone https://github.com/albertozeni/ksw2z-fpga.git
cd ksw2z-fpga
git submodule update --init --recursive
```
Before executing the code in this repo be sure so source XRT and Vitis.
You are required to have Vitis 2020.2 and Xilinx Runtime installed together with all the necessary files to build for the Xilinx Alveo U250.

### Building the bitstream

To build the project simply type:
```
make all TARGET=hw
```
This will both generate the bitstream and the host executable.

To set up the architecture to align sequences up to a certain specific length (default is 256), please edit the file `src/datatypes.h` accordingly.
The proposed architecture will always deploy four cores, each with the indicated `CUs` in `src/datatypes.h`, remember to update the `NUM_CU` define to reflect the length of sequences you want to align.

### Executing the kernel

The command line inputs for our kernel are:
```
./host-ksw2 [bitstream] [npairs] [w] [zdrop]
```
Where `w` is the bandwidth of the target alignment, set it to `-1` if you want the bandwidth to be disabled, and `zdrop` is the value used for the $Z$-drop heuristic,
set this also to `-1` if you want the heuristic to be disabled.
The host will automatically generate random sequences that match the maximum length indicated in the `src/datatypes.h` and executed the same alignments using all the available threads on your host machines,
it will then output performance metrics for our FPGA implementation and the software implementation of KSW2z, also printing out a check for the correspondence of CPU and FPGA results.

## Notes

Our implementation is a full 1:1 replacement for the KSW2z scoring algorithm, as such it can be used to replace KSW2z in all the pipelines that exploit.
Furthermore, our kernel supports the alignment of protein sequences also, as it implements a generic scoring method, making it possible to align sequences with different scores depending on the comparing characters.
A link to the proceedings of our paper @ISCAS23 will be added shortly.
