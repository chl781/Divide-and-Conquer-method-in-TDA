# Divide-and-Conquer-method-in-TDA

This folder contains the code and data for the Wisconsin lake analysis for Li, C. and Cisewski-Kehe, J., 2024. A Divide-and-Conquer Approach to Persistent Homology. arXiv preprint arXiv:2410.01839.

This folder is self-contained, and includes the following:

1.  wi_lakes_analysis.R 
 - This R script takes you through the WI Lakes analysis.  Hopefully there are enough comments to follow the procedures, but more details are available in the paper.  Email Jessi Cisewski-Kehe (jjkehe@wisc.edu) with any additional questions or issues.

2. lakes_wi_updated.csv
- The data were obtained from https://apps.dnr.wi.gov/lakes/lakepages/Results.aspx
- This data file is the cleaned version as described in the manuscript.

3.  pd_north.rds and pd_south.rds
- The northern region and southern region persistence diagrams
- Both sets include 8 persistence diagrams

4.  partition_coordinates.rds
- The latitude and longitude coordinates of the northern and southern regions

5. GetDiagram.R
- This is a ggplot function for persistence diagrams
