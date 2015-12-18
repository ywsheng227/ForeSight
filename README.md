# ForeSight
A Comprehensive Mass Spectrometry Deconvolution Tool for Analyzing Glycosaminoglycan (GAG)</br>
</br>
This software takes the numbers of monosaccharides and adduct ions as input parameters, and generates a database of GAG species within certain m/z and charge ranges.

## How to Use 
### Input 
<b>min_charge</b>: minimum charge</br>
<b>max_charge</b>: maximum charge</br>
<b>dHexA</b>: number of unsaturated hexenuronic acid</br>
<b>HexA</b>: number of hexenuronic acid</br>
<b>HexN</b>: number of hexosamine</br>
<b>HexNAc</b>: number of N-acetylhexosamine</br>
<b>Mann</b>: number of 2,5-anhydromannose</br>
<b>Ac</b>: maximum number of acetate</br>
<b>SO3</b>: maximum number of sulfate</br>
<b>NH4</b>: maximum number of ammonium</br>
<b>HOHloss</b>: water loss (maximum value: 1)</br>
<b>Na</b>: maximum number of sodium ion</br>
<b>K</b>: maximum number of potassium ion</br>
<b>Ca</b>: maximum number of calcium ion</br>
<b>Li</b>: maximum number of lithium ion</br>

### How to generate a GAG database
    exp = GagProspector(4, 8, 200, 2000, 1, 6, 7, 0, 0, 2, 20, 1, 1, 0)
    db, _ = exp.build_database()
    
