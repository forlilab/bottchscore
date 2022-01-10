# bottchscore
Calculate Böttcher score on small molecules as described in [Demoret et al. (ChemRxiv 2020)](https://chemrxiv.org/articles/Synthesis_and_Mechanistic_Interrogation_of_Ginkgo_biloba_Chemical_Space_en_route_to_-Bilobalide/12132939) according to the definition from [Böttcher, J.Chem.Inf.Mod. 2016](https://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00723)

```
usage: bottchscore3.py [-h] -i filename.ext [-m] [-p] [-c] [-v]

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
```
Require OpenBabel v.3.0 or newer.

All reading formats supported by [OpenBabel](http://openbabel.org/wiki/Main_Page) are supported. Although, there might be errors with molecule labels when using the ChemDraw format (just delete molecule labels).

## Support
For help and support, please subscribe to the [CCSB mailing list](https://autodocksuite.scripps.edu/support/).
