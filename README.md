# bottchscore
Calculate Böttcher score on small molecules as described by [Böttcher (J.Chem.Inf.Mod. 2016)](https://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00723)

```Usage: ./bottchinator2000_v3.py molfile.ext   [ all supported OB types are acceptable ] [-v]
Output: Bottch score value

If verbose [-v], the full table with all terms is reported
```

All reading formats supported by [OpenBabel](http://openbabel.org/wiki/Main_Page) are supported. Although, there might be errors with molecule labels when using the ChemDraw format (just delete molecule labels).
