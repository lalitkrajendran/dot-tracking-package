This package implements a dot tracking methodology to estimate image displacements from Background Oriented Schlieren (BOS) experiments, as described in:

Rajendran, L. K., Bane, S. P., & Vlachos, P. P. (2019). Dot tracking methodology for background-oriented schlieren (BOS). Experiments in Fluids, 60(11), 162.

Please cite the above paper if you use the code package for your work.

Run the codes in the following order:
1) batch_run_calibration_example.m: runs a calibration code to obtain prior information regarding the dot pattern such as the dot spacing and arrangement, and generates a mapping function for the camera
2) batch_ptv_processing_example.m: runs the dot tracking method to estimate displacements for the sample data provided with the package
3) batch_calculate_uncertainty_ptv.m: calculates the displacement uncertainty for the tracking estimates 
