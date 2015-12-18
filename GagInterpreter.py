from __future__ import print_function
import logging
import re
from collections import namedtuple
from GagDBdualmode import GagProspector, GAG, gag_to_formula
import operator
from brainpy import isotopic_variants

WINDOW = 7
MIN_CHARGE = 2
MAX_CHARGE = 12
CUTOFF_MZ = 370
CHALLENGE_MATCH_THRESHOLD = 1
NUM_ISOTOPIC_PEAKS = 10
AVG_MASS_DISACCHARIDE = {'hs': 520, 'ha': 400, 'ds': 520, 'cs': 520}

ATOM_MASS = {'C':12, 'H':1.007825, 'N':14.003074, 'O':15.9949146,
             'S':31.972071, 'Na+':22.989222, 'K+':38.963158, 'Ca2+':39.961494,
             'H+':1.007276, '13C':13.00335483, 'Li+':7.015455} 

def read_spec_file(spec_file):
    #read a spectrum file
    regexp = r'[0-9]+\.?[0-9]*'
    mz, intensity = [], []
    while True:
        line = spec_file.readline()
        if not line: break
        if re.match(regexp, line):
            single_mz, single_inten = re.findall(regexp, line)
            mz.append(float(single_mz))
            intensity.append(float(single_inten))
    return mz, intensity

def read_spectrum(spectrum):
    #read a spectrum list
    regexp = r'[0-9]+\.?[0-9]*'
    i, L = 0, len(spectrum)
    while i < L and not re.match('Mass', spectrum[i:]) and not re.match('m/z', spectrum[i:]):
        i += 1
    nums = re.findall(regexp, spectrum[i:])
    mz = [float(nums[i]) for i in range(0, len(nums), 2)]
    intensity = [float(nums[i]) for i in range(1, len(nums), 2)]
    return mz, intensity

def generate_isotopic_clusters_brainpy(gag, charge):
    f = gag_to_formula(gag)
    formula = dict(C=f.C, H=f.H, O=f.O, N=f.N, S=f.S)
    return isotopic_variants(formula, n_peaks=NUM_ISOTOPIC_PEAKS, charge=-charge)

def gagLength(gag):
    return gag.dHexA + gag.HexA + gag.HexN + gag.HexNAc + gag.Mann

class GagIdentifier:
    def __init__(self, spectrum, gag_type, num_charge, accuracy, length_type='Fixed', range_low=0, range_high=0, 
                dHexA=0, Mann=0, NH4=0, Na=0, K=0, Ca=0, Li=0):
        self.spectrum = spectrum
        self.gag_type = gag_type
        self.accuracy = float(accuracy) / 1000000 
        self.num_charge = num_charge
        self.length_type = length_type
        self.range_low = range_low
        self.range_high = range_high
        self.dHexA = dHexA
        self.Mann = Mann
        self.NH4 = NH4
        self.Na = Na
        self.K = K
        self.Ca = Ca
        self.Li = Li
        self.maxi = 0

    def peak_picking(self, mz, intensity):
        #find m/z peaks from spectrum
        i, peak_list = 0, []
        for i in range(len(mz)):
            if mz[i] < CUTOFF_MZ: continue
            self.maxi = max(self.maxi, intensity[i])
        min_intensity = self.maxi / 10000
        i = 0
        while mz[i] < CUTOFF_MZ:
            i += 1
        for j in range(i, len(mz) - WINDOW):
            if max(intensity[j:j+WINDOW]) == intensity[j + WINDOW / 2] and \
               intensity[j + WINDOW / 2] >= min_intensity:
                peak_list += [(mz[j + WINDOW / 2], intensity[j + WINDOW / 2])]
        peak_list.sort()
        return peak_list

    def find_monoisotopic_charge(self, peak_list):
        #find monoisotopic mass and charge of each peak
        L, i, found, mass_list, charge_list = len(peak_list), 0, False, [], []
        while i <= L - 3:
            adj_diff = abs((peak_list[i+1][0] - peak_list[i][0]) - (peak_list[i+2][0] - peak_list[i+1][0]))
            if not adj_diff <= self.accuracy * peak_list[i][0]:
                i += 1
                continue
            charge = int(round(1 / (peak_list[i+1][0] - peak_list[i][0])))
            if not (charge <= MAX_CHARGE and charge >= MIN_CHARGE):
                i += 1
                continue
            mass = (peak_list[i][0] + ATOM_MASS['H+']) * charge
            mass_list.append(mass)
            charge_list.append(charge)
            while i < L - 1 and abs(peak_list[i+1][0] - peak_list[i][0] - (ATOM_MASS['13C'] - ATOM_MASS['C']) / charge) \
                  <= self.accuracy * peak_list[i][0]:
                i += 1
            i += 1
        mass_list.sort()
        return mass_list, charge_list

    def match_peaks_normal(self, mass_list, db):
        #under normal mode, match mass_list to a database 
        matches = []
        for mass in mass_list:
            match = self.binary_search(mass, db)
            if match:
                matched_gag = match[1]
                matches.append([mass, matched_gag])
        matches.sort()
        return matches
    
    def match_peaks_challenge(self, peak_list, db):
        #under challenge mode, match a database to mass_list
        matches = []
        for mz, charge, gag in db:
            isotopic_peaks = generate_isotopic_clusters_brainpy(gag, charge)
            count, matched_peaks, score = 0, [], 0
            for p in isotopic_peaks:
                match = self.binary_search(-p.mz, peak_list)
                if match:
                    count += 1
                    score += match[1] ** 2
                    matched_peaks.append(match) #add matched peak to a list
            if count >= self.matched_isotopic_peaks:
                matches += [((mz + ATOM_MASS['H+']) * charge, mz, charge, gag, score)]
                for peak in matched_peaks:
                    peak_list.remove(peak) #remove matched peaks from peak list
        return matches

    def binary_search(self, target, db):
        #match a target to a database
        lo, hi = 0, len(db) - 1
        while lo <= hi:
            mid = lo + (hi - lo) / 2
            if abs(target - db[mid][0]) <= self.accuracy * db[mid][0]:
                return db[mid]
            elif target < db[mid][0]:
                hi = mid - 1
            elif target > db[mid][0]:
                lo = mid + 1
    
    def identify_species(self):
        #main program
        mz, intensity = read_spectrum(self.spectrum) if isinstance(self.spectrum, basestring) \
                        else read_spec_file(self.spectrum)
        peak_list = self.peak_picking(mz, intensity)
        mass_list, charge_list = self.find_monoisotopic_charge(peak_list)
        
        charge_count = dict((charge, charge_list.count(charge))
                            for charge in charge_list)
        charge_count = sorted(charge_count.items(), key = operator.itemgetter(1),
                              reverse = True)
        charge_high_count = [charge for charge, count in charge_count[:self.num_charge * 2]]

        avg_disaccharide = AVG_MASS_DISACCHARIDE[self.gag_type]
        
        matches = []
        if self.length_type == 'Fixed':
            length = int(mass_list[int(len(mass_list)*0.6)] / float(avg_disaccharide) * 2)
            length = length if length % 2 == 0 else length + 1
            lengths = [length - 2, length] 
        else:
            lengths = range(self.range_low, self.range_high + 1)

        decon_type = 'challenge'
        total_num_gag = 0
        for length in lengths:
            self.matched_isotopic_peaks = min(6, max(4, length / 3))
            if self.gag_type == 'hs':
                HexNAc = 0
                HexN = length / 2 - self.Mann if length % 2 == 0 else (length + 1) / 2 - self.Mann
            else:
                HexN = 0
                HexNAc = length / 2 - self.Mann if length % 2 == 0 else (length + 1) / 2 - self.Mann

            if self.dHexA:
                HexA = length / 2 - self.dHexA if length % 2 == 0 else (length + 1) / 2 - self.dHexA
            else:
                HexA = length / 2
            Ac = min(3, HexN) 
            SO3 = HexN * 2 + HexNAc * 2 + self.dHexA + HexA + self.Mann * 2
            HOHloss = True if self.dHexA else False

            gag = GagProspector(min(charge_high_count), max(charge_high_count), CUTOFF_MZ, max(intensity),
                                           self.dHexA, HexA, HexN, HexNAc, self.Mann,
                                           Ac, SO3, self.NH4, HOHloss, self.Na, self.K, 
                                           self.Ca, self.Li, decon_type)
            db, num_gag = gag.build_database()
            total_num_gag += num_gag
        
            if decon_type == 'normal':
                matches += self.match_peaks_normal(mass_list, db)
            else:
                matches += self.match_peaks_challenge(peak_list, db)

            if self.dHexA:
                gag = GagProspector(min(charge_high_count), max(charge_high_count), CUTOFF_MZ, max(intensity),
                                           0, self.dHexA + HexA, HexN, HexNAc, self.Mann,
                                           Ac, SO3, self.NH4, 0, self.Na, self.K, 
                                           self.Ca, self.Li, decon_type)
                db, num_gag = gag.build_database()
                total_num_gag += num_gag
        
                if decon_type == 'normal':
                    matches += self.match_peaks_normal(mass_list, db)
                else:
                    matches += self.match_peaks_challenge(peak_list, db)

        return matches
                                                                          
