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
<b>length_type</b>: 'Fixed' or 'Flexible'</br>
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
    spectrum = #some spectrum file
    gag_type = 'hs'
    accuracy = 3
    num_charge = 3
    length_type = 'Flexible'
    range_low = 4
    range_high = 10
    dHexA = 1
    exp = GagIdentifier(spectrum, gag_type, accuracy, num_charge, length_type, range_low, range_high, dHexA)
    matches = exp.identify_species()
    
    for mass, mz, charge, gag, _ in matches:
            
    
