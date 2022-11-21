## Files used in Ogane et al.

This repository contains scripts and Jupyter notebooks used in "Development of Hidden Markov Modeling Method for Molecular Orientations and Structure Estimation from High-Speed Atomic Force Microscopy Time-Series Images" by Tomonori Ogane, Daisuke Noshiro, Toshio Ando, Atsuko Yamashita, Yuji Sugita, and Yasuhiro Matsunaga.

## Descriptions on files

The files are organized as follows:

- [src/afm.jl](https://github.com/matsunagalab/paper_ogane2022/blob/main/src/afm.jl)

  - Defines functions for hidden Makrov modeling (including the Viterbi and the Baum-Welch algorithms)
 
- [md/](https://github.com/matsunagalab/paper_ogane2022/tree/main/md)

  - Contains input files for the molecular dynamics simulations of a coarse-grained model for taste receptor
 
- [data/quaternion](https://github.com/matsunagalab/paper_ogane2022/tree/main/data/quaternion)

  - Quaternions for uniform rotations in the SO(3) group
 
- [data/t1r](https://github.com/matsunagalab/paper_ogane2022/tree/main/data/t1r)

  - Contains cluster structures (Markov state model structures) obtained by the clustering of MD trajectories
 
- [00_setup/](https://github.com/matsunagalab/paper_ogane2022/tree/main/00_setup)

  - Contains a script to generate pseudo-AFM images for likelihood calculations
 
- [01_psudo_test/01_generate_realization01.jl](https://github.com/matsunagalab/paper_ogane2022/blob/main/01_psudo_test/01_generate_realization01.jl)

  - Script to generate pseudo-AFM images for twin experiment
 
- [01_psudo_test/02_analyze_realization01.jl](https://github.com/matsunagalab/paper_ogane2022/blob/main/01_psudo_test/02_analyze_realization01.jl)

  - Calculates frame-by-frame likelihood for the pseudo AFM images generated by [paper_ogane2022/01_psudo_test/01_generate_realization01.jl](https://github.com/matsunagalab/paper_ogane2022/blob/main/01_psudo_test/01_generate_realization01.jl)
 
- [01_psudo_test/03_0_calc_rmsd_freq_origin.jl](https://github.com/matsunagalab/paper_ogane2022/blob/main/01_psudo_test/03_0_calc_rmsd_freq_origin.jl)

  - Frame-by-frame rigid-body fitting analysis
 
- [01_psudo_test/03_1_calc_rmsd_freq_MD.jl](https://github.com/matsunagalab/paper_ogane2022/blob/main/01_psudo_test/03_2_calc_rmsd_freq_baumwelch.jl)

  - Applies the Viterbi algorithm with the transition probabilities of the ground-truth
 
- [01_psudo_test/03_2_calc_rmsd_freq_baumwelch.jl](https://github.com/matsunagalab/paper_ogane2022/blob/main/01_psudo_test/03_2_calc_rmsd_freq_baumwelch.jl)

  - Applies the Baum-Welch algorithm to estimate the transition probabilities and then applies the Viterbi algorithm

- [01_psudo_test/03_3_calc_rmsd_freq_rand.jl](https://github.com/matsunagalab/paper_ogane2022/blob/main/01_psudo_test/03_3_calc_rmsd_freq_rand.jl)

  - Applies the Viterbi algorithm with a transition probabilities created by a random matrix

## Large files on OneDrive

Pseudo-AFM images (in the file format of [BSON](https://github.com/JuliaIO/BSON.jl)) for twin experiments generated by [01_psudo_test/02_analyze_realization01.jl](https://github.com/matsunagalab/paper_ogane2022/blob/main/01_psudo_test/02_analyze_realization01.jl) are uploaded on OneDrive, which can be downloaded from [this link (15 GB)](https://suitc-my.sharepoint.com/:u:/g/personal/ymatsunaga_mail_saitama-u_ac_jp/EeXYLbP6Zl5Mr4_V2mAB7qEB4Oj3zuIa1Gwi205QZ0idGA?e=wVo6RJ). 

## Required packages

Codes in thie repository are writtein in Julia programming language. 
You need to install julia to run the codes. 
Also, some codes depend on several packages. 
The packages can be installed as follows:

```julia
$ julia
julia> 
# enter the package mode by pressing ]
pkg> add IJulia Plots Statistics StatsBase ProgressMeter Revise
pkg> add https://github.com/matsunagalab/MDToolbox.jl.git
# return to the REPL mode by pressing BACKSPACE or DELETE
julia> using IJulia
julia> exit()
```
 
## License

This repository is licensed under the under the terms of GNU General Public License v3.0. 

Quaternion data contained in `data/quaternion/` directory were taken from the repository of the BioEM program written by Cossio et al. https://github.com/bio-phys/BioEM. These are separately licensed under the terms of the GNU General Public License. Please check the license file `data/quaternion/LICENSE`. 


## Contact

Yasuhiro Matsunaga

ymatsunaga@mail.saitama-u.ac.jp

