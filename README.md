# BioUMAT: Abaqus Fortran subroutine for cartilage multiphasic modeling

This code is the Fortran 77 version of the subroutines of the cartilage model I firstly proposed in my MSc thesis, and I used it with some modifications in several publications. For more information, read our relevant [paper](https://shayansss.github.io/files/2019_09_preprint.pdf).

## Citation
If this research data is useful for your work, kindly please consider citing our work ([DOI](http://doi.org/10.1177/0954411919854011) | [PDF](https://shayansss.github.io/files/2019_09_preprint.pdf)):

```
@article{doi:10.1177/0954411919854011,
author = {Seyed Shayan Sajjadinia and Mohammad Haghpanahi and Mohammad Razi},
title ={Computational simulation of the multiphasic degeneration of the bone-cartilage unit during osteoarthritis via indentation and unconfined compression tests},
journal = {Proceedings of the Institution of Mechanical Engineers, Part H: Journal of Engineering in Medicine},
volume = {233},
number = {9},
pages = {871-882},
year = {2019},
doi = {10.1177/0954411919854011},
}
```

## Important notes
- In the paper, there are a few typing errors (that might be clear from the context), including: 1) in Table 1, the symbol of the initial solid volume fraction was mistakenly written by that of its updated value; 2) in Figure 1, the numbers of layers were mistakenly illustrated from 1 to 9, while they should be from 9 to 1. These errors did not exist in the code and therefore did not affect the results, and since they might be realized from the context, we published no correction.
- When the model is run, the cartilage part should first expand due to the pre-stressing effect. Accordingly, we manually made the model smaller to reach the correct geometry (i.e., the stress-free condition), as shown in Figure 1, and then after inducing the initial osmotic pre-stressing, the loading or boundary conditions were applied. For this reason, your results might differ slightly from ours (while not affecting the conclusions of this study), since we did not illustrate the manually developed stress-free geometries. We published another [paper](https://shayansss.github.io/files/2021_02.pdf), where the cartilage model is extended using a novel automatic pre-stressing algorithm.
- The DDSDDE array (i.e., the consistent Jacobian matrix) was implemented based on the derived equations in the paper. Alternatively for future studies, an approximation method was also implemented as the exact array we developed here may fail to work well if the code is changed or if some extensive local transition between the compressibility and incompressibility occurs (which was not handled in the exact DDSDDE code as it was not of importance for our simulation tests).
