# bottchscore
Calculate Böttcher score on small molecules as described in [Demoret et al. (ChemRxiv 2020)](https://chemrxiv.org/articles/Synthesis_and_Mechanistic_Interrogation_of_Ginkgo_biloba_Chemical_Space_en_route_to_-Bilobalide/12132939) according to the definition from [Böttcher, J.Chem.Inf.Mod. 2016](https://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00723)

```
Usage: ./bottchscore3.py molfile.ext   [ all supported OB types are acceptable ] [-v]
Output: Bottch score value

If verbose [-v], the full table with all terms is reported
```
Require OpenBabel v.3.0 or newer.

All reading formats supported by [OpenBabel](http://openbabel.org/wiki/Main_Page) are supported. Although, there might be errors with molecule labels when using the ChemDraw format (just delete molecule labels).
