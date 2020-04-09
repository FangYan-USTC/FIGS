# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import pymzml
from bidict import bidict
from scipy import sparse
import struct
import zlib
import sqlite3
import sparse_nnls

def LibraryNormalization(SpectraLibrary):
    for key in SpectraLibrary.keys():
        lib = SpectraLibrary[key]
        spectrum = np.array(lib['Spectrum'])
        spectrum[:,1] = spectrum[:,1]/sum(spectrum[:,1])
        spectrum = spectrum[spectrum[:,0].argsort()]
        lib['Spectrum'] = spectrum
    return SpectraLibrary

def LoadMS2(path):
    DIArun = pymzml.run.Reader(path)
    E = enumerate(DIArun)
    MS2 = [[spectrum.peaks('raw'),(spectrum['MS:1000827']-spectrum['MS:1000828'], spectrum['MS:1000827']+spectrum['MS:1000829']),
            spectrum['MS:1000016'],i] for i,spectrum in E if spectrum['ms level']==2.0]
    return MS2

def DivideLibraryByWindows(SpectraLibrary, windows):
    DividedLibraries = {}
    for window in windows:
        DividedLibraries[window] = {key:value for key, value in SpectraLibrary.items() if window[0] <= value['PrecursorMZ'] < window[1]}
    return DividedLibraries

def GenerateDecoyLibrary(SpectraLibrary, distance):
    DecoyLibrary = {}
    charges = set(np.array(list(SpectraLibrary.keys()))[:,1].astype(int))
    for charge in charges:
        LibraryDividedByCharge = {key:value for key, value in SpectraLibrary.items() if key[1]==charge}
        DecoyLibraryDividedByCharge = {}
        swapped_keys = bidict()

        key_premz = []
        for key, value in LibraryDividedByCharge.items():
            key_premz.append([key, value['PrecursorMZ']])
        key_premz = np.array(key_premz)
        key_premz = key_premz[np.argsort(key_premz[:,1])]

        i, j = 0, 1
        while i < len(key_premz) and j < len(key_premz):
            i_key, i_premz, j_key, j_premz = key_premz[i][0], key_premz[i][1], key_premz[j][0], key_premz[j][1]
            if abs(j_premz-i_premz) >= distance:
                decoy_i_key = ('DECOY-'+i_key[0], i_key[1])
                decoy_j_key = ('DECOY-'+j_key[0], j_key[1])
                
                DecoyLibraryDividedByCharge[decoy_i_key] = SpectraLibrary[i_key].copy()
                DecoyLibraryDividedByCharge[decoy_i_key]['Spectrum'] = SpectraLibrary[j_key]['Spectrum'].copy()

                DecoyLibraryDividedByCharge[decoy_j_key] = SpectraLibrary[j_key].copy()
                DecoyLibraryDividedByCharge[decoy_j_key]['Spectrum'] = SpectraLibrary[i_key]['Spectrum'].copy()
                
                delta = i_premz - j_premz
                DecoyLibraryDividedByCharge[decoy_i_key]['Spectrum'][:,0] += delta
                DecoyLibraryDividedByCharge[decoy_j_key]['Spectrum'][:,0] -= delta

                swapped_keys[i_key] = j_key
            else:
                j += 1

            while (i < len(key_premz)) and (key_premz[i][0] in (set(swapped_keys) | set(swapped_keys.inv))):
                i += 1
            while (j < len(key_premz)) and (key_premz[j][0] in (set(swapped_keys) | set(swapped_keys.inv))):
                j += 1

        if len(LibraryDividedByCharge) != len(DecoyLibraryDividedByCharge):
            unswapped_keys = set(LibraryDividedByCharge.keys()) - (set(swapped_keys) | set(swapped_keys.inv))
            unswapped_keys = list(unswapped_keys)
            unswapped_keys.sort()

            for unswapped_key in unswapped_keys:
                unswapped_premz = SpectraLibrary[unswapped_key]['PrecursorMZ']
                for swapped_key in sorted(swapped_keys):
                    swapped_key2 = swapped_keys[swapped_key]
                    if abs(unswapped_premz - SpectraLibrary[swapped_key]['PrecursorMZ'])>= distance and \
                       abs(unswapped_premz - SpectraLibrary[swapped_key2]['PrecursorMZ'])>= distance:
                        DecoyLibraryDividedByCharge[('DECOY-'+unswapped_key[0], unswapped_key[1])] = SpectraLibrary[unswapped_key].copy()
                        DecoyLibraryDividedByCharge[('DECOY-'+unswapped_key[0], unswapped_key[1])]['Spectrum'] = SpectraLibrary[swapped_key2]['Spectrum'].copy()

                        DecoyLibraryDividedByCharge[('DECOY-'+swapped_key[0], swapped_key[1])]['Spectrum'] = SpectraLibrary[unswapped_key]['Spectrum'].copy()
                        
                        DecoyLibraryDividedByCharge[('DECOY-'+swapped_key[0], swapped_key[1])]['Spectrum'][:,0] += (SpectraLibrary[swapped_key]['PrecursorMZ']-unswapped_premz)
                        DecoyLibraryDividedByCharge[('DECOY-'+unswapped_key[0], unswapped_key[1])]['Spectrum'][:,0] += (SpectraLibrary[unswapped_key]['PrecursorMZ']-SpectraLibrary[swapped_key2]['PrecursorMZ'])
                        
                        break
                if swapped_key in swapped_keys: swapped_keys.pop(swapped_key)
        DecoyLibrary.update(DecoyLibraryDividedByCharge)
    return DecoyLibrary

def cal_nnls(LibIntensity, MS2Intensity, penalty):
    RowIndex = list(range(len(LibIntensity) + 1))
    ColIndex = [0] * (len(LibIntensity) + 1)
    LibIntensity.append(penalty)

    MS2Intensity.append(0)
    MS2Intensity = np.array(MS2Intensity)

    LibraryVector = sparse.coo_matrix((LibIntensity, (RowIndex, ColIndex)))
    LibraryCoeffs = sparse_nnls.lsqnonneg(LibraryVector, MS2Intensity, {'show_progress': False})
    LibraryCoeffs = LibraryCoeffs['x']
    LibraryCoeffs = LibraryCoeffs[0]

    return LibraryCoeffs

def LoadBlib(libPath):
    Lib = sqlite3.connect(libPath)
    LibPrecursorInfo = pd.read_sql("SELECT * FROM RefSpectra", Lib)

    SpectraLibrary = {}

    for i in range(len(LibPrecursorInfo)):
        precursorID = str(LibPrecursorInfo['id'][i])
        precursorKey = (LibPrecursorInfo['peptideModSeq'][i], LibPrecursorInfo['precursorCharge'][i])
        NumPeaks = pd.read_sql("SELECT numPeaks FROM RefSpectra WHERE id = " + precursorID, Lib)['numPeaks'][0]

        SpectrumMZ = pd.read_sql("SELECT peakMZ FROM RefSpectraPeaks WHERE RefSpectraID = " + precursorID, Lib)['peakMZ'][0]
        SpectrumIntensities = pd.read_sql("SELECT peakIntensity FROM RefSpectraPeaks WHERE RefSpectraID = " + precursorID, Lib)['peakIntensity'][0]

        if len(SpectrumMZ) == 8 * NumPeaks and len(SpectrumIntensities) == 4 * NumPeaks:
            SpectraLibrary.setdefault(precursorKey, {})
            SpectrumMZ = struct.unpack('d' * NumPeaks, SpectrumMZ)
            SpectrumIntensities = struct.unpack('f' * NumPeaks, SpectrumIntensities)
            SpectraLibrary[precursorKey]['Spectrum'] = np.array((SpectrumMZ, SpectrumIntensities)).T
            SpectraLibrary[precursorKey]['PrecursorMZ'] = LibPrecursorInfo['precursorMZ'][i]
            SpectraLibrary[precursorKey]['PrecursorRT'] = LibPrecursorInfo['retentionTime'][i]
        elif len(SpectrumIntensities) == 4 * NumPeaks:
            SpectraLibrary.setdefault(precursorKey, {})
            SpectrumMZ = struct.unpack('d' * NumPeaks, zlib.decompress(SpectrumMZ))
            SpectrumIntensities = struct.unpack('f' * NumPeaks, SpectrumIntensities)
            SpectraLibrary[precursorKey]['Spectrum'] = np.array((SpectrumMZ, SpectrumIntensities)).T
            SpectraLibrary[precursorKey]['PrecursorMZ'] = LibPrecursorInfo['precursorMZ'][i]
            SpectraLibrary[precursorKey]['PrecursorRT'] = LibPrecursorInfo['retentionTime'][i]
        elif len(SpectrumMZ) == 8 * NumPeaks:
            SpectraLibrary.setdefault(precursorKey, {})
            SpectrumMZ = struct.unpack('d' * NumPeaks, SpectrumMZ)
            SpectrumIntensities = struct.unpack('f' * NumPeaks, zlib.decompress(SpectrumIntensities))
            SpectraLibrary[precursorKey]['Spectrum'] = np.array((SpectrumMZ, SpectrumIntensities)).T
            SpectraLibrary[precursorKey]['PrecursorMZ'] = LibPrecursorInfo['precursorMZ'][i]
            SpectraLibrary[precursorKey]['PrecursorRT'] = LibPrecursorInfo['retentionTime'][i]
        elif len(zlib.decompress(SpectrumMZ)) == 8 * NumPeaks and len(
                zlib.decompress(SpectrumIntensities)) == 4 * NumPeaks:
            SpectraLibrary.setdefault(precursorKey, {})
            SpectrumMZ = struct.unpack('d' * NumPeaks, zlib.decompress(SpectrumMZ))
            SpectrumIntensities = struct.unpack('f' * NumPeaks, zlib.decompress(SpectrumIntensities))
            SpectraLibrary[precursorKey]['Spectrum'] = np.array((SpectrumMZ, SpectrumIntensities)).T
            SpectraLibrary[precursorKey]['PrecursorMZ'] = LibPrecursorInfo['precursorMZ'][i]
            SpectraLibrary[precursorKey]['PrecursorRT'] = LibPrecursorInfo['retentionTime'][i]

    return SpectraLibrary
