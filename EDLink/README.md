## Dependencies

- Linux, MacOS, FreeBSD:
  - GCC ≥ 5 [limited GCC-4.9 support on Linux]
  - Clang/LLVM ≥ 3.6 [limited Clang-3.5 support on Linux]
  - Intel Compiler ≥ 17.0.0 on Linux
- Windows:
  * Visual C++ ≥ 14.0 / Visual Studio ≥ 2015
  * Intel Compiler ≥ 17.0.0 / Visual Studio ≥ 2015u2
  * Clang/C2 ≥ 3.8.0 / Visual Studio ≥ 2015u3 [experimental, requires CMake ≥ 3.6]
- Boost
- OpenMP


## To compile

`./compile.sh`


## To build FM-Index

``` bash
./generator -b ${fasta_file}
${fasta_file}: reference fasta file
```


## To build tree

``` bash
./generator -l ${fasta_file} -b ${tree_level} -d ${tree_dir}
${fasta_file}: reference fasta file
${tree_level}: depth of the tree
${tree_dir}: dir to store the tree info
```

## To calculate neighbors of the tree

``` bash
./processor ${tree_dir} ${ED_threshold} ${max_depth} ${result_dir}
${tree_dir}: dir that stores tree
${ED_threshold}: threshold to store F
${max_depth}: the maximum depth of nodes to process in the tree
${result_dir}: directory to store all the results
```


## Example

``` bash
./generator -b ../ref.fasta
./generator -l ../ref.fasta -b 8 -d tree_dir
./processor tree_dir 2 4 ed_dir
```

