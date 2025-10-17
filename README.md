# SSIPe
SSIPe is a method to calculate binding free-energy changes (ΔΔ<i>G</i><sub>bind</sub>) of protein-protein interactions (PPIs) upon mutations at protein-protein interface. Starting from a PPI complex structure, SSPIe first generates multple structure and sequence alignments from non-redundant interface library (NIL) and STRING databases separately. The structural and sequence profiles are then combined with the physical energy function EvoEF to predict the impact of the mutations on PPI binding free energies. SSIPe can be used to guide the designing and engineering of protein-protein interactions with enhanced binding affinity, and/or for understanding the roles of disease-related mutations associated with protein-protein interactions.

Due to the restriction of uploading large files (> 25 Mb per file), only the source code, binaries and training and test sets are uploaded to this repository. The NIL and STRING databases, which are 1.1 Gb and 2.4 Gb respectively, cannot be uploaded, but they can be found at https://zhanggroup.org/SSIPe/download/. The source code and training/test data are also provided in our website. After the two libraries (NIL and STRING) are downloaded from the website, they can be unzipped and put into the directory /bin/ssip. Alternatively, users can download the whole package (source code + databases) from https://zenodo.org/records/17364388.

# SSIPe server and standalone program
The SSIPe webserver is available at https://zhanggroup.org/SSIPe/. The SSIPe standalone program is available at https://zenodo.org/records/17364388. Please can contact Dr. Xiaoqiang Huang (xiaoqiah@umich.edu or xiaoqiah@outlook.com) to report any issues.

# Run the SSIPe standalone program
Note that SSIPe can only run in a Linux environment. Users can download the program and libraries which are mentioned above and directly run the run_SSIPe.pl script following the instructions in example/README. The results will be printed out in result.txt. Before running the program, the user should check if some of the external binaries are executable. Specifically, the user should check:<br>
bin/evoef/EvoEF/src/EvoEF<br>
bin/extint/extint<br>
bin/ssip/ialign/extint<br>
bin/ssip/ialign/IS-align<br>
bin/ssip/ncbi-blast-2.7.1+/bin/psiblast<br>
bin/ssip/ncbi-blast-2.7.1+/bin/makeblastdb<br>
bin/ssip/seqalign/build_cpx_alignment<br>
bin/ssip/seqalign/PDB2FAS<br>
bin/ssip/SSIPserver/SSIP<br>
Due to size restriction, the psibalst and makeblastad program are uploaded as zipped file (.zip), the user can unzip them into the same directory. For the other programs, the user may need to use 'chmod +x' to make the other binaries executable or if neccessary, the user can recompile and build the binaries. For recompiling EvoEF, the user can go to the bin/evoef/EvoEF/src/ directory and run 'g++ -O3 -o EvoEF \*.cpp' to rebuild EvoEF. For recompiling SSIP, the user can go to the bin/ssip/SSIPserver/ and run 'g++ -O3 -o SSIP src/\*.cpp' to rebuild SSIP.

# Reference

Huang X, Zheng W, Pearce R, Zhang Y. SSIPe: accurately estimating protein-protein binding affinity change upon mutations using evolutionary profiles in combination with an optimized physical energy function. Bioinformatics (2020) 36:2429-2437.
