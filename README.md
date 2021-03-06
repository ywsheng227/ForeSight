# ForeSight
A Comprehensive Mass Spectrometry Deconvolution Tool for Analyzing Glycosaminoglycan (GAG)</br>
</br>
This software identifies GAG species from a mass spectrum. 

## How to Use 
### Input parameters
<b>spectrum</b>: a spectrum file in plain text format</br>
<b>gag_type</b>: type of GAG ('hs': heparan sulfate, 'cs': chondroitin sulfate, 'ds': dermatan sulfate, 'ha': hyaluronic acid)</br>
<b>accuracy</b>: mass accuracy in ppm</br>
<b>num_charge</b>: minimum number of charge states</br>
<b>length_type</b>: 'Fixed' or 'Variable'</br>
<b>range_low</b>: lower end of length range</br>
<b>range_high</b>: higher end of length range</br>
<b>dHexA</b>: number of unsaturated hexeruronic acid</br>
<b>Mann</b>: number of 2,5-anhydromannose</br>
<b>NH4</b>: maximum number of ammonium</br>
<b>Na</b>: maximum number of sodium ion</br>
<b>K</b>: maximum number of potassium ion</br>
<b>Ca</b>: maximum number of calcium ion</br>
<b>Li</b>: maximum number of lithium ion</br>

### How to deconvolute a mass spectrum
    from gaginterpreter import GagInterpreter
    
    spectrum = #some spectrum file
    gag_type = 'hs'
    accuracy = 3
    num_charge = 3
    length_type = 'Variable'
    range_low = 4
    range_high = 10
    dHexA = 1
    exp = GagIdentifier(spectrum, gag_type, accuracy, num_charge, length_type, range_low, range_high, dHexA)
    matches = exp.identify_species()
    
    for mass, mz, charge, gag, _ in matches:
        print "%.4f\t%.4f\t%d\t%d %d %d %d %d %d %d %d %d %d %d %d %d" % (mass, mz, charge, gag.dHexA, gag.HexA, gag.HexN, gag.HexNAc,
            gag.Mann, gag.Ac, gag.SO3, gag.NH4, gag.HOHloss, gag.Na, gag.K, gag.Ca, gag.Li)
    
