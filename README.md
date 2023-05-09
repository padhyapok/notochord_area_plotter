# notochord_area_plotter

MATLAB code extracting notochord segment areas accompanying the Current Biology submission titled "Dynamic BMP signaling mediates notochord segmentation in zebrafish". The input for this code uses fluorescent images labeling segmentation in the zebrafish notochord. The code extracts clusters of expression after making masks of the notochord surface and thresholding gaussian filtered images. Clusters detected within a certain distance are merged and identified as part of a single segment using the accompanying function 'join_clusters_close'.
