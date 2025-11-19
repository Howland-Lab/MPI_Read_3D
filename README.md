# MPI Fortran program for reading and diagnosing large PadeOps output files

If reading and analyzing PadeOps output files using the Python scripts in PadeOpsIO is too heavy, this Fortran program can be used. 

The program uses the 2decomp_fft library to read PadeOps output files, which is the same library used by PadeOps to export the output files. 

The program relies on MPI to divide reading the heavy PadeOps output files among the MPI processes, making it much faster than a numpy-based read in Python. 

The **build_*.sh** attached script should be modified bases on the HPC platform. Ideally, set **DECOMP_DIR** to the same 2decomp_fft path that PadeOps uses. 

The current version of the program performs horizontal averaging of different fields as specified by the input file **input.dat**. The program writes the generated horizontally averaged profiles to **csv** files. Further capabilities will be included, such as exporting slices through a three-dimensional field.

To use the program, please run the build shell script, after modification based on the used platform. A sample job script is included (**run**), which should also be modified accordingly.