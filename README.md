[![Build Status](https://travis-ci.org/Crespo-Otero-group/fromage.svg?branch=master)](https://travis-ci.org/Crespo-Otero-group/fromage) [![Docs](
https://readthedocs.org/projects/fromage/badge/?version=latest&style=plastic)](https://fromage.readthedocs.io/en/latest/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/Crespo-Otero-group/fromage.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/Crespo-Otero-group/fromage/context:python)

<p align="left">
  <img height="150" src="doc/logo.png">
</p>

**fromage** (FRamewOrk for Molecular AGgregate Excitations) is a Python framework designed to facilitate the study of molecular aggregates in the excited state. It contains utilities for geometry manipulation going from periodic to finite models, exciton analysis and ONIOM calculations. Here we present our new version 2.0

**fromage** is developed at Queen Mary University of London by the [Crespo-Otero group](https://crespootero.wordpress.com/). We acknowledge the Leverhulme Trust (RPG-2019-122).

The documentation can be found [here](https://fromage.readthedocs.io/).

To cite the use of the program, please use:

Rivera, M., Dommett, M., Sidat, A., Rahim, W., Crespo-Otero, R. fromage: A library for the study of molecular crystal excited states at the aggregate scale. *J Comput Chem* 2020; 1-14. https://doi.org/10.1002/jcc.26144

And if you are using one of the ONIOM implementations:

Rivera, M., Dommett, M., Crespo-Otero, R. ONIOM(QM:QM') Electrostatic Embedding Schemes for Photochemistry in Molecular Crystals. *J. Chem. Theory Comput.* 2019; 15, 4, 2504-2516 https://doi.org/10.1021/acs.jctc.8b01180

## Features

   - Nonadiabatic molecular dynamics with Surface Hopping (GSH, FSSH):
     Generalized FSSH and ZNSH with nonadibatic coupling and spin-orbit coupling
   - Constrained DFT (charge and/or spin) optimizations 
   - Conical intersections and S/T crossing search in the crystal, gas phase and solution (implicit solvent)
   - Supports: Gaussian, Turbomole, Q-Chem, NWChem, Orca, Molcas, xTB, DFTB+ and MOPAC (including FOMO-CI)
    
## 1 Installation
### 1.1 Easy installation using pip
1. Make sure that you have the following installed:

  - Python 3.7 and above
  - swig
  - [Modified version of Ewald](https://github.com/Crespo-Otero-group/Ewald) (only necessary for Ewald embedding calculations)
  - Cython. Necessary for the Surface Hopping module.

2. Clone this repository to wherever you want to install it:

  ```bash
   cd /path/to/dir/
   git clone https://github.com/Crespo-Otero-group/fromage.git
   cd fromage/
  ```

3. Install

  ```bash
  sudo pip install .
  ```
  N.B. this will install `numpy`, `scipy` and `pytest`  on your system.

4. Set environment variables. In your `.bashrc`, write
  ```bash
  export FRO_GAUSS=g16
  export FRO_EWALD=Ewald
  ```
  If you are using different binaries for Gaussian or Ewald, change accordingly.
  
  Voila!
  
### 1.2 Common pitfall

When installing, you might find the error:
```
Python.h: No such file or directory
```
In order to install Python packages which contain C or C++ on Linux, you need `python-dev` or `python3-dev` which provides the header `Python.h`.

## 2 Preparing the calculation

`fro_prep_run.py` requires:

- A `.xyz` file of a unit cell
- A `config` file
- A file containing the high level embedding charges
- A file containting the low level embedding charges
- A set of `*.template` files for `mh`, `ml` and `rl` (`mg` too for MECI search)

And optionally:

- A `.xyz` target shell file
- A self consistent calculation template file

### 2.1 Unit cell file

This is the kind of file which is typically produced from a periodic DFT calculation. Avoid double counted atoms outside the unit cell which might be snuck in by visualization programs such as VESTA. For this kind of problem, the `fro_uc_tools.py` utility can come in useful. Just make a file containing the unit cell vectors and run `fro_uc_tools.py cell.xyz vectors -d`. If you deem that geometry optimisation of your cell is unnecessary you will need to convert a `.cif` file into a unit cell `.xyz` file. For this we recommend [Open Babel](http://openbabel.org/wiki/Main_Page) which would use the syntax:

```bash
    babel -i -cif cell.cif -o -xyz cell.xyz --filluc
```

The `--filluc` keyword is crucial otherwise you will end up with an asymmetric unit of the cell.

### 2.2 Configuration file

The `config` file is a list of keywords followed by their values which the user should input. A set of reasonable default values is coded in `parse_config_file` for every keyword except for the lattice vectors which can not be assumed.

If you use the `target_shell` keyword you will need to supply the program with an additional target shell `.xyz` file to specify the shell that you want in your real system.

### 2.3 Population analysis files

These files contain the starting (and in some cases also ending) values for the charges that you intend to use in your embedding. The high level charges will eventually be used in the embedding of the `mh` calculation and can be fitted to the Ewald potentially directly or self consitently. They need to be extracted from a Gaussian output file or in the special case of direct Ewald embedding a CP2K >= 4.1 file will do as well. For the low level embedding, only Gaussian is available.

It is crucial that you match the population analysis from the low level charges to the low level of theory in order to properly cancel out the electrostatic intreactions from `rl`. The choice of high level charges is more subtle but consistency would have you use the same level of theory as in `mh`.

### 2.4 Primitive template files

These files are named `mh.template`, `ml.template`, `rl.template` and optionally `mg.template`. They serve as model template files with a blank name for the checkpoint file, blank calculation name, blank atomic positions and blank point charges. Prepare calculations populates all of these fields except for the atomic positions which will later be repeatedly updated by the optimisation procedure. It is important to include the following keywords:

- `charge` : allows for the use of point charge embedding (not actually used in `rl`)
- `symmetry=none` : conserves the input geometry throughout the calculation, making the position of the charges correct
- `force` : to calculate the energy gradients necessary for the optimisation

### 2.5 Target shell file

In certain cases, the cluster of molecules generated radially will be impractical due to the packing of the crystal. For example it may need to include a large number of distant molecules to also include a certain nearest neighbour molecule. In those cases a shell file can be manually edited in the user's favourite chemistry visualisation software to remove extra molecules. In that case a target shell file can be supplied which will be used to generate `rl` and `ml`. Be extra careful that the central molecule has the correct orientation with respect to your generated cluster. It is recommended to use the shell file from a large radius calculation and manually strip it down to only the nearest neighbour molecules.

### 2.6 Self consistent template file

This is the Gaussian template file used if the Ewald charge background is computed self consistently. It has the same form as the primitive template files. The level of theory here can be chosen to be in excited state for a fully excited crystal. If this is your intention, don't forget to use `density=current` to make sure that your population analysis is in the excited state and not the ground state. Of course ground state self-constistent calculations are also possible.

### 2.7 Outputs of the preparation

Once you have finished running `fro_prep_run.py`, you will end up with a few files:

- `prep.out` which gives you information about how your preparation went. If the last line gives you an ending time, that is good news
- `mol.init.xyz` will be the initial position of your molecule for the optimisation
- `shell.xyz` is the cluster of molecules without the central one
- `mh`, `ml`, `rl` and `mg` directories containing `.temp` files corresponding to all of the parallel calculations
- `ewald` directory where the ewald calculation is run

### 2.8 Prepare the input for nonadiabatic dynamics calculation

 - This step comes after the fromage calculation has been prepared
 - `dynamixsampling.py -i freq_file` must be invoked to prepare the initial conditions (initconds file). 
   Currently, only Orca, G16, Bagel and Molden formats are supported
 - `setup_dyn.py` must be invoked to prepare all the necessary files to run the trajectories.

### 2.9 fromage-dynamics (Also check the README.md in `fro_dyn_test/`)

    List of fromage.in file parameters:

          Name |   Format   | Description
    -----------|------------|-------------------------------------------------------
       types   |  List<str> | List of atomic symbols for each atom in high layer
               |            |
         R     |  np.array  | XYZ coordinates for each atom, units = Angstrom
               |            |
         M     |  np.array  | Atomic masses for each atom, units = AMU
               |            |
         V     |  np.array  | XYZ velocities for each atom, units = Bohr/a.u.
               |            |
       State   |    <str>   | Root of initial electronic state (ground state = "1")
               |            |
      natoms   |    <int>   | Number of atoms in high layer
               |            |
     stepsize  |    <str>   | Time step for trajectory, units = femtoseconds
               |            |
         T     |    <str>   | Total simulation time, units = femtoseconds
               |            |
      HopType  |    <str>   | Surface-hopping method to use, one of "FSSH", "GSH, "NOSH"
               |            |
      nactype  |    <str>   | Nonadiabatic coupling to compute, "nac" or "ktdc"
               |            |
         spin  |    <int>   | Spin numbers
               |            |
       states  |    <str>   | Number of states in each spin multiplicity
               |            |    
    singlestate|    <int>   | Compute gradient for current state only (1) or all state (0)
               |            |
       temp    |    <str>   | (optional) Temperature for NVT simulation, units = Kelvin
               |            |
       check   |    <str>   | (optional) Name of checkpoint file to write
               |            |
      chk_stp  |    <int>   | (optional) Number of steps every which the system's state should be written to restart file (default 25).
               |            |
    hop_thresh |    <str>   | (optional) Threshold for computing Zhu-Nakamura internal conversion surface-hopping
               |            |
    isc_thresh |    <str>   | (optional) Threshold for computing Zhu-Nakamura intersystem crossing surface-hopping
               |            |
    e_cons     |    <str>   | (optional) Threshold for total energy conservation (default 0.2 eV)
               |            |
    dyn_restart|    <str>   | (optional) Option to restart the dynamics from the last step
               |            |
    -----------|------------|-------------------------------------------------------


        List of possible init_dict parameters:

         Name  |   Format   | Description
       --------|------------|-------------------------------------------------------
        types  |   <list>   | List of atomic symbols for each atom
               |            |
          R    | <np.array> | XYZ coordinates for each atom, units = Angstrom
               |            |
          M    | <np.array> | Atomic masses for each atom, units = a.u. of mass (not AMU)
               |            |
          V    | <np.array> | XYZ velocities for each atom, units = Bohr/a.u.
               |            |
        state  |    <int>   | Root number of initial electronic state (ground state = 1)
               |            |
       natoms  |    <int>   | Number of atoms in the high layer
               |            |
         spin  |    <int>   | Spin numbers
               |            |
       states  |    <int>   | Number of states in each spin multiplicity
               |            |
    singlestate|    <int>   | Compute gradient for current state only (1) or all state (0)
               |            |
     stepSize  |   <float>  | Trajectory time step, units = femtoseconds
               |            |
          T    |   <float>  | Total simulation time, units = femtoseconds
               |            |
        temp   |   <float>  | (optional) Temperature of NVT simulation, units = Kelvin
               |            |     no default, NVE ensemble is used if no temp is provided
               |            |
      HopType  |    <str>   | (optional) Type of surface-hopping probability,
               |            |     one of "FSSH", "GSH", "NOSH", default = "NOSH"
               |            |
      nactype  |    <str>   | (optional) Type of nonadiabatic coupling,
               |            |     "nac" for exact computed NACs, "ktdc" for approximated NACs
               |            |
        check  |    <str>   | (optional) Name of checkpoint file to use for writing data
               |            |     default = "trajectory.chk"
               |            |
         gap   |   <float>  | (optional) Energy-difference threshold for computing 
               |            |     Zhu-Nakamura internal conversion surface-
               |            |     hopping probability, units = eV, default = 0.5
               |            |
      gapsoc   |   <float>  | (optional) Energy-difference threshold for computing 
               |            |     Zhu-Nakamura intersystem crossing surface-
               |            |     hopping probability, units = eV, default = 0.5
               |            |
       --------|------------|-------------------------------------------------------
        
        
        
     List of internal Trajectory parameters

       Name |   Format   | Description
     -------|------------|-------------------------------------------------------
        Vs  | <np.array> | Thermostat array, set and used by NoseHoover()
            |            |
        E   | <np.array> | Array of state-specific energies, units = Eh
            |            |
        G   | <np.array> | Array of state-specific gradients, units = Eh/Bohr
            |            |
        N   | <np.array> | Array of Nonadiabatic coupling, units = 1/Bohr
            |            |
       iter |    <int>   | Current step number of the trajectory
            |            |
        t   |   <float>  | Current total simulation time, units = femtoseconds
            |            |
     Thermo |   <bool>   | Switch to turn on NVT ensemble and Nose-Hoover thermostat
            |            |
       Maxh |    <int>   | Maximum number of hopping between states
            |            |
       Delt |   <float>  | Probability integration time step, unit = a.u.
            |            |
      Aprev | <np.array> | Density matrix for previous step
            |            |
      Hprev | <np.array> | Energy matrix for previous step
            |            |
      Dprev | <np.array> | Nonadiabatic coupling matrix for previous step
            |            |
      Acurr | <np.array> | Density matrix for current step
            |            |
      Hcurr | <np.array> | Energy matrix for current step
            |            |
      Dcurr | <np.array> | Nonadiabatic coupling matrix for current step
            |            |
       Ekin |   <float>  | Current kinetic energy, units = Eh
            |            |
       Deco |   <float>  | Decoherence energy, units = Eh
            |            |
     Hopped |    <int>   | Switch for hopping at the current time step
            |            |
     -------|------------|-------------------------------------------------------

## 3 Running the calculation

### 3.1 Last few steps before running

To run fromage, all you need is:


- A `fromage.in` file
- `mol.init.xyz`
- `shell.xyz`
- `mh`, `ml` and `rl` directories containing `.temp` files (`mg` is only needed for MECI serach)
- `velocity` (only needed for a SH dynamics calculation)

The `fromage.in` file has a similar structure to the `config` file but is much simpler and is not even necessary for geometry optimisation in Gaussian. If you want to change a program used in a specific level of theory from Gaussian ot something else, simply add `high_level [program]` or `low_level [program]`. For MECI search, add `bool_ci 1`. For a SH dynamics add `dynamics` (and `dyn_restart` to restart a SH dynamics). Also, when Q-Chem, Molcas, NWChem or MOPAC are used, the variable `nprocs  <int>` has to be added. 

For Turbomole RI-CC2, run a define and then write in all of the point charges from `mh.temp` under the block `$point_charges` after scaling.

**IMPORTANT**: Turbomole uses Bohr units in its control file and as such the x, y and z columns should be scaled accordingly

For Molcas RASCF, prepare an input file in the directory called `molcas.input` with geom.xyz as the coordinate. To add point charges, use in `&GATEWAY`:


```
  xfield
[number of charges] Angstrom
 -9.237176  -2.137048   3.557237   0.432218
-10.014996  -1.455739   0.568597  -0.168284
-10.112382  -2.173251   1.384633   0.146427
                .
                .
                .
```

You should be all set now. Run `fro_run.py` to begin the calculation.

#### 3.1.1 MECI search in gas phase or solution (implicit solvent)

   - You only need `fromage.in`, `mol.init.xyz` `mh` and `mg` directories containing `.temp` files. 
   - Run `ci_search.py` to begin the calculation 

## 4 Additional utilities

A couple of useful utilities are included here for manipulation of molecular crystal clusters. They may come in handy when making sure that your calculation is doing what you want.

### 4.1 fro_pick_mol.py

This program selects molecules out of a cluster of molecules and writes them to another file. This is particulary useful when you are dealing with a large molecular cluster and want to extract something like a dimer.

The syntax is intuitive:

```bash
  fro_pick_mol.py cluster_file.xyz 1 56 22
```

In this case we have selected a trimer of the molecules containing the atoms labeled 1, 56 and 22 in the input .xyz file. This program will get angry if you feed it more than one atom label of one same molecule. For additional options use `fro_pick_mol.py --help`

### 4.2 fro_assign_charges.py

This is more of a debugging tool for checking to see that you are using a sensible bond length in your definition of your molecules. It reads charges and positions from a Gaussian output file for one molecule and assigns those charges to a cluster of atoms made up of those same molecules in a .`xyz`-like file.

It could come in handy as a standalone utility if you want to assign charges in a forcefiled calculation.

The syntax is:

```bash
  fro_assign_charges.py population_analysis.log cluster_of_molecules.xyz
```

As usual, use `assign_charges.py --help` for more options related to bond length, charge kind etc.

### 4.3 fro_uc_tools.py

This utility does operations on xyz files paired up with unit cell vector files. The vector file is of the form:
```
8.9638004303         0.0000000000         0.0000000000
0.0000000000        10.5200004578         0.0000000000
-3.8748910079         0.0000000000        10.7924653741
```
Options include extracting the nonequivalent monomers from the cell, generating a tessalating cell but with all complete molecules (therefore spilling out of the bounding box of the unit cell vectors), confining a cell to the bounding box and creating supercells.

### 4.4 Analyse dimers in aggregate geoemtries

The script `fro_dimer_tools.py` can identify the unique dimers in the supplied geometry, taking into account periodicity if relevant. The dimers can further be characterised by the angles between their principal, secondary and normal axes, as well as their centroid-to-centroid distance. As before, many parameters can be altered so using `fro_dimer_tools.py -h` is encouraged.


### 4.5 Voronoi volume evaluation

It can be useful to determine the available volume of a molecule in an aggregate environment. To do this, one could use the union of the van der Waals volumes of each atom, or the Voronoi volume of the molecule, scaled by van der Waals radii.

`fro_volumetrics.py clust.xyz -l 13`

This will produce cube files of the available volume of the molecule containing atom 13 (-l 13) within the cluster of molecules (clust.xyz). The ouputs are the Voronoi volume (voro.cube), the van der Waals volume (vdw.cube) and the union of the two (add.cube). A file called volumes prints the integrated volume of each of the three.

### 4.6 Exciton Coupling Evaluation

Exciton coupling evaluation from Gaussian output files can also be carried out, using fro_coupling.py. A diabatisation of the Hamiltonian is employed which relies on the calculation of excited state properties such as population analysis or transition dipole moments.[4] More options are also available.

As an example of use, the line:

`fro_coupling.py -m DIA -p TDM -mf a.log b.log -of dim_ab.log -os 2`

Will use the diabatisation method (-m DIA) and use the transition dipole moment property (-p TDM) to read the Gaussian log files of monomer S1 calculations (-mf a.log b.log) and the dimer S2 calculation (-of dim_ab.log) with state of interest S2 (-os 2). The output will show the diabatic Hamiltonian, whose off-diagonal elements are the exciton coupling values.

### 4.7 And more!!!


## 5 Some parting words

If you find yourself with a bunch of error files during your optimisation, ask yourself where some calculations might have failed in the geometry optimisation. Maybe some SCF did not converge or your central molecule escaped the cluster to be with its one true love: infinitely attractive point charges. To combat this, try adding more molecules to your cluster.

If you are happily preparing a calculation, see no error, run `fro_run.py` with sensible geometries and find that your quantum chemistry program is very upset about something, it might be that the point charges you are feeding it are unreasonable. Indeed when you start fiddling with large numbers of constrained point charges in Ewald, the system of linear equations which fits them becomes highly linearly dependent and you end up with point charges with values in the thousands. If this happens to you just try a smaller number of constrained point charges while still containing your central molecule.

Hacking this program for your own personal needs is a perfectly good idea and you may be able to get a lot from just importing the I/O modules and the `Atom` and `Mol` classes.

If all you want to do is integrate your favourite quantum chemistry package into fromage, all you need to do is a) add new io routines in `read_file` and `edit_file` b) make a new `Calc` object modeled after one of the existing ones in the `calc` module c) Add the class and corresponding keyword to the `calc_types` at the top of the `calc` module

The Ewald program is often the source of all of your problems when tinkering with the embedding methods, even as a regular user pushing the program to its limits. It uses a deprecated lapack function and needs to be modified very specifically to be used with `fro_prep_run.py`.

--------------------------------------------------------------------------------

More detailed instructions can be found in the [documentation](https://fromage.readthedocs.io/).
For any questions about usage, citing or contributing, please email our group at f.hernandez@qmul.ac.uk or r.crespo-otero@qmul.ac.uk

- Miguel Rivera and Federico Hernandez
