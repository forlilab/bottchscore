# bottchscore (Improved Accuracy)
Calculate Böttcher score on small molecules as described in [Demoret et al. (ChemRxiv 2020)](https://chemrxiv.org/articles/Synthesis_and_Mechanistic_Interrogation_of_Ginkgo_biloba_Chemical_Space_en_route_to_-Bilobalide/12132939) according to the definition from [Böttcher, J.Chem.Inf.Mod. 2016](https://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00723). This repository contains a modified script originally developed by ForliLab (https://github.com/forlilab/bottchscore). The modified script accounts for E/Z double bond isomers and changes the s_i term accordingly for the appropriate atom. An option to change the maximum memory size for automorphic calculations was also added, allowing to aviod the memory limit errors for highly symmetrical molecules. The complexity score can also be calculated by directly calling the appropriate function from another script. 

```
usage: bottchscore4.py [-h] -i filename.ext [-m] [-p] [-c] [-v] [-x MaxMemory]

Tool to calculate Böttcher score on small molecules.

optional arguments:
  -h, --help       show this help message and exit
  -i filename.ext  input structure; support all input formats supported by OB,
                   including multi-structure formats
  -m               disable mesomeric effect estimate
  -p               generate PNG image of the structure
  -c               add a progressive counter to the list of results shown
  -v               verbose mode; print the full table of the terms used to
                   estimate the score, as described in the paper
  -x MaxMemory     specify the maximum memory that will be available for the automorphism/symmetry calculations; the
                   default value is set to 3000000
```
Require OpenBabel v.3.0 or newer.

All reading formats supported by [OpenBabel](http://openbabel.org/wiki/Main_Page) are supported. Although, there might be errors with molecule labels when using the ChemDraw format (just delete molecule labels).

This version of the program supports calculating scores for optical stereoisomers as well as E/Z double bond isomers as described in [Böttcher, J.Chem.Inf.Mod. 2016](https://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00723); however, it does not support axial isomers (axial chirality).

Python function:
```
calculate_bottchscore_from_smiles(smiles: str, verbose_response=False, debug_arg=False, disable_mesomer=False, automorp_memory_maxsize=3000000) -> float
```
can be called to calculate a Böttcher score for a molecule directly from SMILES passed to the function as a string. 

The modified script returns the same scores as the original script (see the 'Calculated Examples' folder), however the accuracy improvement is seen for such molecules as:
```
SMILES                                            Original Score               New Score
CC(/C=C/C1=CC=CC=C1)=O                            73.43                        80.60 (the same as in SI figure S6 of the Bottcher paper)
Cl/C=C\C=C\Br                                     54.25                        68.59
CC/C(C1=CC=CC=C1)=C(C2=CC=CC=C2)/CC               60.26                        66.26
C/C(=C(/C=C/C)\CCC)/CC                            64.34                        83.51
CC/C=C(C)/[2H]                                    24.34                        31.51
CC/C=C1CCC[C@H](Br)C/1                            99.80                        106.97
```
Compounds that are not E/Z isomer but have a double bond also work as expected:
```
SMILES                                Original Score              New Score
C/C=C(C)/C                            19.17                       19.17
CC/C=C1CCCCC/1                        38.17                       38.17 
```
