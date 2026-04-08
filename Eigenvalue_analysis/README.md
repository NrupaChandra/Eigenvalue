# MATLAB Eigenvalue Analysis Script

JSON File: Make sure the JSON file (combined_data_standardized.json) containing the eigenvalue data is in the same directory as the script or provide the correct path in the configuration section.

**Place Files:**
Ensure that the script and the JSON data file (combined_data_standardized.json) reside in the same folder or update the filename path accordingly.


**Analysis Options**

The script offers the following analysis options:

Plot Eigenvalues for a Specific Re and ω
Visualize eigenvalues for a target Reynolds number and a designated complex frequency. The script calculates an "implied alpha" based on the provided Re for comparison.

Compare Ratios (ω/Eigenvalue) with Group Velocity
Compute the ratio of ω to each eigenvalue and plot the relationship with both real and imaginary components of the group velocity. An optional filter lets you restrict the analysis to a specific Reynolds number.

Plot Unstable Eigenvalues for a Specified Range of Re
Display unstable eigenvalues (for both 2D and 3D modes) within a given Reynolds number range. The script prints out each unstable eigenvalue’s associated Re and ω.

Plot Unstable Eigenvalues for a Specified ω Range (for a Given Re)
For a fixed Reynolds number, this option plots only those unstable eigenvalues where the real part of ω falls within a user-defined range.

Combined Plot (Target Re & ω with Unstable Modes over a Re Range)
Combines a plot of the eigenvalues for a specified (Re, ω) case with unstable eigenvalue data over a Reynolds number range. It also displays the calculated "implied alpha" marker.

Plot Unstable Eigenvalues for All Re for a Given ω Range
Visualize unstable eigenvalues for data entries matching a specified range of the ω real part. The plotted markers are color-coded based on the Reynolds number.

Combined Plot of Option 1 and Option 6 in a Single Plot
Merges the plots from Option 1 (specific eigenvalues and implied alpha) and Option 6 (unstable eigenvalues over all Re within a certain ω range) into one comprehensive visualization. The plot is enhanced with color coding based on Reynolds numbers.


