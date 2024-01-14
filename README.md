# LAMMPS_intermol_heatflux
Paper: [Atomistic Mechanism of Anisotropic Heat Conduction in the Liquid Crystal 4-Heptyl-4′-cyanobiphenyl: All-Atom Molecular Dynamics](https://pubs.acs.org/doi/10.1021/acs.jpcb.9b08158)  

Decomposing nonbonded intra/inter-molecular contribution in "compute stress/atom" command is implemented for thermal conductivity decomposition analysis.  
The sample inputs are located in "sample inputs" directory.  
1) Calculate thermal conductivity by s[Num]e[Num]Ne[Num].in
2) Thermal conductivity decomposition analysis as the post analysis by post_ssk.in

The current input files of the decomposition is based on the stress tensor obtained by "compute stress/atom" command.
**However, the "compute stress/atom" command has been reported to inadequate to calculate the contribution from angle, dihedral, improper and constraint interactions (See https://docs.lammps.org/compute_heat_flux.html).**  
The version of LAMMPS in this repository is 16Mar2018, so the centroid/stress/atom commands officially implemented in 3Mar2020 cannot be used.  

For accurate thermal conductivity decomposition, you are strongly advised to use "compute centroid/stress/atom" command for angle, dihedral, improper terms by using LAMMPS after version 3Mar2020.  
It should be noted that "compute stress/atom" command can be used for two-body intra/inter-molecular interactions (ex. Bonding, Lennard-Jones, Coulomb, etc...).   

