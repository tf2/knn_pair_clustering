# knn_pair_clustering
Input data format is expected to be a csv file containing at least 4 columns with a header

xnm - x axis cooridinates

ynm - y axis cooridinates

znm - z axis cooridinates

locprecnm - error

NOTE: error recalulation method should be refined

Basic description: iterative k nearest neighbour approach, k1 neighbours paired and collapsed, followed by overlapping pairs being grouped together into clusters and so on until all clusters are closed. Then we search again across both the singles and clusters having adjusted the errors of each cluster based on the maximum union of the error spheres within each cluster to see if any additional points/clusters can be merged together. And so on and so on until we end when there are no more mergeable points/clusters. 
