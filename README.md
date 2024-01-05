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

## Performance Analysis

Comparison of our design against software baseline [KSW2z](https://github.com/lh3/ksw2/blob/master/ksw2_extz2_sse.c) executed using 52-threads on an Intel Xeon 8167M using different combinations of sequence length, Z, and bandwidth on 1M sequence pairs.

Alveo U250 Exec. Time [s] | Xeon 8167M Exec. Time [s] | Sequence Length | Bandwidth | Z-drop | Speedup wrt. Xeon | Energy Eff. wrt. Xeon |
|--------: |--------:	|--------: |--------:	|--------: |--------: |--------:
1.30  | 6.40  | 256  | -   | -   | 4.94 | 12.86|
3.64  | 23.60 | 512  | -   | -   | 6.49 | 16.91|
11.69 | 90.00 | 1024 | -   | -   | 7.70 | 20.07|
1.02  | 6.40  | 256  | 250 | -   | 6.29 | 16.39|
2.54  | 17.76 | 512  | 250 | -   | 6.99 | 18.23|
7.43  | 40.36 | 1024 | 250 | -   | 5.43 | 14.15|
1.30  | 6.40  | 256  | -   | 250 | 4.93 | 12.86|
3.27  | 14.21 | 512  | -   | 250 | 4.35 | 11.33|
7.92  | 14.77 | 1024 | -   | 250 | 1.87 | 4.86 |
1.02  | 6.39  | 256  | 250 | 250 | 6.25 | 16.28|
2.80  | 10.26 | 512  | 250 | 250 | 3.66 | 9.54 |
7.43  | 10.27 | 1024 | 250 | 250 | 1.38 | 3.61 |

## Notes

Our implementation is a full 1:1 replacement for the KSW2z scoring algorithm, as such it can be used to replace KSW2z in all the pipelines that exploit.
Furthermore, our kernel supports the alignment of protein sequences also, as it implements a generic scoring method, making it possible to align sequences with different scores depending on the comparing characters.

## Citation

To cite our work or to know more about our methods, please refer to:

```
@inproceedings{zeni2023genome,
  title={On the genome sequence alignment fpga acceleration via ksw2z},
  author={Zeni, Alberto and Di Donato, Guido Walter and Della Valle, Alessia and Carloni, Filippo and Santambrogio, Marco D},
  booktitle={2023 IEEE International Symposium on Circuits and Systems (ISCAS)},
  pages={1--5},
  year={2023},
  organization={IEEE}
}
```
