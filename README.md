## Scripts and notebooks used in Ogane et al. 2022

This repository contains scripts and Jupyter notebooks used in "Hidden Markov Modeling of Molecular Orientations and Structures from High-Speed Atomic Force Microscopy Time-Series Images" by Tomonori Ogane, Daisuke Noshiro, Toshio Ando, Atsuko Yamashita, Yuji Sugita, and Yasuhiro Matsunaga (2022).

## Descriptions on files

The files are organied as the following:

- `src/01_generate_test_data.jl` 

 - Script to generate pseudo AFM images for twin experiment. 
 
- `src/02_analyze.jl` 

 - Calculate likelihoods of frame-by-frame analysis.
 
- `src/03_0_calc_rmsd_freq_origin.jl - 03_4_calc_rmsd_freq_rand_no_rotate.jl` 

 - Calculate RMSDs from the ground truth structures with various likelihoods.
 
- `src/fig01.ipynb` 

 - Auxiliary notebook for generating Fig. 1.
 
- `src/fig02.ipynb` 

 - Auxiliary notebook for generating Fig. 2.
 
- `src/fig03.ipynb` 

 - Notebook to generate Fig. 3.
 
- `src/fig04.ipynb` 

 - Notebook to generate Fig. 4.
 
- `src/fig05.ipynb` 

 - Notebook to generate Fig. 5.
 
- `src/fig06.ipynb` 

 - Notebook to generate Fig. 6.
 
## Contact

- Yasuhiro Matsunaga
- ymatsunaga@mail.saitama-u.ac.jp


