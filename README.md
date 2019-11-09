# SSIPe
SSIPe is a method to calculate binding free-energy changes (ΔΔGbind) of protein-protein interactions (PPIs) upon mutations at protein-protein interface. Starting from a PPI complex structure, SSPIe first generates multple structure and sequence alignments from non-redundant interface library (NIL) and STRING databases separately. The structural and sequence profiles are then combined with the physical energy function EvoEF to predict the impact of the mutations on PPI binding free energies. SSIPe can be used to guide the designing and engineering of protein-protein interactions with enhanced binding affinity, and/or for understanding the roles of disease-related mutations associated with protein-protein interactions.
  Due to the restriction of uploading large files (> 25 Mb per file), only the source code, binaries and training and test sets are uploaded to this repository. The NIL and STRING databases, which are 1.1 Gb and 2.4 Gb respectively, cannot be uploaded, but they can be found at https://zhanglab.ccmb.med.umich.edu/SSIPe/download/. The source code and training/test data are also provided in our website. After the two libraries (NIL and STRING) are downloaded from the website, they can be unzipped and put into the directory /bin/ssip.
