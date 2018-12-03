
<p align="left">
  <img height="150" src="doc/logo.png">
</p>

**fromage** (FRamewOrk for Molecular AGgregate Excitations) is a Python framework designed to facilitate the study of molecular aggregates in the excited state. 

**fromage** is developed at Queen Mary University of London by the [Crespo-Otero group](https://crespootero.wordpress.com/).

The documentation can be found [here](https://fromage.readthedocs.io/).


## 1 Installation

1. Make sure that you have all the required dependencies:

  - Python 2.7 and above (Python 3 recommended)
  - numpy
  - scipy
  - swig
  - [Modified version of Ewald](https://github.com/Crespo-Otero-group/Ewald) (only necessary for Ewald embedding calculations)

2. Clone this repository to wherever you want to install it:

  ```bash
   cd /path/to/dir/
   git clone https://github.research.its.qmul.ac.uk/btx156/fromage.git
   cd fromage/
  ```

3. Compile

  ```bash
  cd fromage/fdist/
  swig -c++ -python fdist.i
  cd ../../
  sudo python setup.py build_ext install
  ```
4. Set environment variables. In your `.bashrc`, write
  ```bash
  export FRO_GAUSS=g16
  export FRO_EWALD=Ewald
  ```
  If you are using different binaries for Gaussian or Ewald, change accordingly.
  
  Voilà!

The main two modules in fromage are `prepare_calculation.py` and `fromage.py`. The former produces template files and geometry files to be used in the latter for geometry optimisation or minimal energy conical intersection (MECI) search. Following the standard ONIOM nomenclature, the central system which is only treated at the high level of theory is called the 'model' system. The whole system is called the 'real' system. As such the three parallel calculations which are carried out are called `mh` for `model high`, `ml` for `model low`, `rl` for `real low`. For MECI search, an additional calculation of the high level region for the ground state gradients is necessary and is labeled `mg`.

## 2 Preparing the calculation

`prepare_calculation.py` requires:

- A `.xyz` file of a unit cell
- A `config`
- A file containing the high level embedding charges
- A file containting the low level embedding charges
- A set of `*.template` files for `mh`, `ml` and `rl` (`mg` too for MECI search)

And optionally:

- A `.xyz` target shell file
- A self consistent calculation template file

### 2.1 Unit cell file

This is the kind of file which is typically produced from a periodic DFT calculation. Avoid double counted atoms outside the unit cell which might be snuk in by visualization programs such as VESTA. For this kind of problem, the `uc_tools.py` utility can come in useful. Just make a file containing the unit cell vectors and run `uc_tools.py cell.xyz vectors -d`. If you deem that geometry optimisation of your cell is unnecessary (publish this result because I am interested) you will need to convert a `.cif` file into a unit cell `.xyz` file. For this I recommend [Open Babel](http://openbabel.org/wiki/Main_Page) which would use the syntax:

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

Once you have finished running `prepare_calculation.py`, you will end up with a few files:

- `prep.out` which gives you information about how your preparation went. If the last line gives you an ending time, that is good news
- `mol.init.xyz` will be the initial position of your molecule for the optimisation
- `fixed_cell.xyz` is the unit cell after any necessary transformations have been applied to complete important molecules
- `clust.xyz` is the real system i.e. a cluster of molecules with `mol.init.xyz` in the middle
- `shell.xyz` is the cluster of molecules without the central one
- `mh`, `ml`, `rl` and `mg` directories containing `.temp` files corresponding to all of the parallel calculations

## 3 Running the calculation

### 3.1 Last few steps before running

To run fromage, all you need is:

- A `fromage.in` file
- `mol.init.xyz`
- `shell.xyz`
- `mh`, `ml` and `rl` directories containing `.temp` files (`mg` is only needed for MECI serach)

The `fromage.in` file has a similar structure to the `config` file but is much simpler and is not even necessary for geometry optimisation in Gaussian. If you want to change a program used in a specific level of theory from Gaussian ot something else, simply add `high_level [program]` or `low_level [program]`. For MECI search, add `bool_ci 1`

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

You should be all set now. Run `run_fromage.py` to begin the calculation.

### 3.2 Outputs

The program only has three main outputs:

- `fromage.out` which gives updates on all of the individual energies being calculated, the total gradient norm and the energy gap
- `geom_mol.xyz` which keeps a record of the optimising geometry
- `geom_cluster.xyz` which combines `geom_mol.xyz` and `shell.xyz` for a better view of intermolecular interactions

If the minimisation ran smoothly, the last line of `fromage.out` should be the ending time.

## 4 Additional utilities

A couple of useful utilities are included here for manipulation of molecular crystal clusters. They may come in handy when making sure that your calculation is doing what you want.

### 4.1 pick_mol.py

This program selects molecules out of a cluster of molecules and writes them to another file. This is particulary useful when you are dealing with a large molecular cluster and want to extract something like a dimer.

The syntax is intuitive:

```bash
  pick_mol.py cluster_file.xyz 1 56 22
```

In this case we have selected a trimer of the molecules containing the atoms labeled 1, 56 and 22 in the input .xyz file. This program will get angry if you feed it more than one atom label of one same molecule. For additional options use `pick_mol.py --help`

### 4.2 assign_charges.py

This is more of a debugging tool for checking to see that you are using a sensible bond length in your definition of your molecules. It reads charges and positions from a Gaussian output file for one molecule and assigns those charges to a cluster of atoms made up of those same molecules in a .`xyz`-like file.

It could come in handy as a standalone utility if you want to assign charges in a forcefiled calculation.

The syntax is:

```bash
  assign_charges.py population_analysis.log cluster_of_molecules.xyz
```

As usual, use `assign_charges.py --help` for more options related to bond length, charge kind etc.

### 4.3 uc_tools.py

This utility does operations on xyz files paired up with unit cell vector files. The vector file is of the form:
```
8.9638004303         0.0000000000         0.0000000000
0.0000000000        10.5200004578         0.0000000000
-3.8748910079         0.0000000000        10.7924653741
```
Options include extracting the nonequivalent monomers from the cell, generating a tessalating cell but with all complete molecules (therefore spilling out of the bounding box of the unit cell vectors), confining a cell to the bounding box and creating supercells.

## 5 Some parting words

If you find yourself with a bunch of `core.` files during your optimisation, ask yourself where some calculations might have failed in the geometry optimisation. Maybe some SCF did not converge or your central molecule escaped the cluster to be with its one true love: infinitely attractive point charges. To combat this, try adding more molecules to your cluster.

If you are happily preparing a calculation, see no error, run `fromage.py` with sensible geometries and find that your quantum chemistry program is very upset about something, it might be that the point charges you are feeding it are unreasonable. Indeed when you start fiddling with large numbers of constrained point charges in Ewald, the system of linear equations which fits them becomes highly linearly dependent and you end up with point charges with values in the thousands. If this happens to you just try a smaller number of constrained point charges while still containing your central molecule.

Hacking this program for your own personal needs is a perfectly good idea and you may be able to get a lot from just importing the I/O modules and the `Atom` and `Mol` classes.

If all you want to do is integrate your favourite quantum chemistry package into fromage, all you need to do is a) add new io routines in `read_file` and `edit_file` b) make a new `Calc` object modeled after one of the existing ones in the `calc` module c) Add the class and corresponding keyword to the `calc_types` at the top of the `calc` module

The Ewald program is often the source of all of your problems when tinkering with the embedding methods, even as a regular user pushing the program to its limits. It uses a deprecated lapack function and needs to be modified very specifically to be used with `prepare_calculation`.

--------------------------------------------------------------------------------

More detailed instructions can be found in the [documentation](https://fromage.readthedocs.io/).
For any questions about usage, citing or contributing, please email our group at r.crespo-otero@qmul.ac.uk

-Miguel Rivera
