# Divide-and-Conquer-method-in-TDA

This folder contains the code and data for Li, C. and Cisewski-Kehe, J., 2024. A Divide-and-Conquer Approach to Persistent Homology. arXiv preprint arXiv:2410.01839.

Full citation:  [To be added once finalized]



This is the divide-and-conquer method for computing a persistence diagram using a Vietoris-Rips filtration, implemented in R. 

There are two illustration implementation codes for the 1D and 2D cases, where test data is provided in the data repo. Wisconsin lake data implementation is in large_sample_implementation; this is a parallel implementation of this method for high-performance computing.

We provide detailed implementation code in each folder.


1.  2D-DaC-Example.R
- Illustrates the proposed DaC method on a 2D point cloud that includes a large and a small loop.
- Outputs a DaC H1 persistence diagram.


2.  3D-DaC-Example.R
- Illustrates the proposed DaC method on a 3D point cloud of a 2-sphere.
- Outputs a DaC H2 persistence diagram


3.  data
- Simulated data sets used in the paper, including a 1-dimensional circle, a 2-dimensional sphere, the Stanford Bunny example, and 2 disjoint circles.


4. Functions
- Support functions for the DaC method when m=1.


5.  Functions3Combine
- Support functions for the DaC method when m=2.


6.  large_sample_implementation
- Support functions for the DaC method are implemented using parallel computation.


7.  wi_lakes_analysis
- This folder is self-contained and includes all the data and files necessary to carry out the WI lakes analysis presented in the manuscript and supplementary material.




