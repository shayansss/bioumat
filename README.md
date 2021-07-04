# Fortran UMAT, FLOW, and SDVINI subroutines for a multiphasic and fibril-reinforced model of articular cartilage

This code is the Fortran 77 version of the UMAT, FLOW, and SDVINI subroutines of the cartilage model, I firstly proposed in my Master's thesis, and I used it with some modifications in several publications.

## Citation
If this research data is useful for your work, kindly please consider citing:

*Sajjadinia, Seyed Shayan, Mohammad Haghpanahi, and Mohammad Razi. "Computational simulation of the multiphasic degeneration of the bone-cartilage unit during osteoarthritis via indentation and unconfined compression tests." Proceedings of the Institution of Mechanical Engineers, Part H: Journal of Engineering in Medicine 233.9 (2019): 871-882.<br/>[DOI](http://doi.org/10.1177/0954411919854011){:target="_blank"} | [PDF](https://shayansss.github.io/files/2019_09_preprint.pdf){:target="_blank"}


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
    note ={PMID: 31232647},
URL = {https://doi.org/10.1177/0954411919854011}
}
```

## Further information
Note that for the DDSDDE array, I used two different equations, one for debugging of the code and the other is based on equations in the appendix of the [paper](https://shayansss.github.io/files/2019_09_preprint.pdf){:target="_blank"}. The latter should work better if it is used together with the right linearization method, otherwise use the former.

Please, read the [paper](https://shayansss.github.io/files/2019_09_preprint.pdf){:target="_blank"} to find more information about implementation of this code.
