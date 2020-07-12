This code is the Fortran 77 version of the UMAT, FLOW, and SDVINI subroutines of the cartilage model, I firstly proposed in my Master's thesis. The model with some minor modifications was used in several publications.

If you use this code in your work, please cite the following article as the reference, where you can also find more information about this code at:

https://doi.org/10.1177/0954411919854011

Please note that for the DDSDDE array, I used two different equations, one for debugging of the code and the other is based on the above paper. The latter should work better if and only if it is used together with the right linearization method, otherwise use the former if you cannot find the best linearization method.
