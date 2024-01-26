## Periodic table for PyRAIMD
## Jingbai Li Feb 10 2020

class Element:
    ## This class is periodic table
    ## This class read atom name in various format
    ## This class return atomic properties

    def __init__(self,name):

        Periodic_Table = {
              "HYDROGEN"     :  1,   "H":  1,   "H":  1,   "1":  1,     1:  1,
              "HELIUM"       :  2,  "He":  2,  "HE":  2,   "2":  2,     2:  2,
              "LITHIUM"      :  3,  "Li":  3,  "LI":  3,   "3":  3,     3:  3,
              "BERYLLIUM"    :  4,  "Be":  4,  "BE":  4,   "4":  4,     4:  4,
              "BORON"        :  5,   "B":  5,   "B":  5,   "5":  5,     5:  5,
              "CARBON"       :  6,   "C":  6,   "C":  6,   "6":  6,     6:  6,
              "NITROGEN"     :  7,   "N":  7,   "N":  7,   "7":  7,     7:  7,
              "OXYGEN"       :  8,   "O":  8,   "O":  8,   "8":  8,     8:  8,
              "FLUORINE"     :  9,   "F":  9,   "F":  9,   "9":  9,     9:  9,
              "NEON"         : 10,  "Ne": 10,  "NE": 10,  "10": 10,    10: 10,
              "SODIUM"       : 11,  "Na": 11,  "NA": 11,  "11": 11,    11: 11,
              "MAGNESIUM"    : 12,  "Mg": 12,  "MG": 12,  "12": 12,    12: 12,
              "ALUMINUM"     : 13,  "Al": 13,  "AL": 13,  "13": 13,    13: 13,
              "SILICON"      : 14,  "Si": 14,  "SI": 14,  "14": 14,    14: 14,
              "PHOSPHORUS"   : 15,   "P": 15,   "P": 15,  "15": 15,    15: 15,
              "SULFUR"       : 16,   "S": 16,   "S": 16,  "16": 16,    16: 16,
              "CHLORINE"     : 17,  "Cl": 17,  "CL": 17,  "17": 17,    17: 17,
              "ARGON"        : 18,  "Ar": 18,  "AR": 18,  "18": 18,    18: 18,
              "POTASSIUM"    : 19,   "K": 19,   "K": 19,  "19": 19,    19: 19,
              "CALCIUM"      : 20,  "Ca": 20,  "CA": 20,  "20": 20,    20: 20,
              "SCANDIUM"     : 21,  "Sc": 21,  "SC": 21,  "21": 21,    21: 21,
              "TITANIUM"     : 22,  "Ti": 22,  "TI": 22,  "22": 22,    22: 22,
              "VANADIUM"     : 23,   "V": 23,   "V": 23,  "23": 23,    23: 23,
              "CHROMIUM"     : 24,  "Cr": 24,  "CR": 24,  "24": 24,    24: 24,
              "MANGANESE"    : 25,  "Mn": 25,  "MN": 25,  "25": 25,    25: 25,
              "IRON"         : 26,  "Fe": 26,  "FE": 26,  "26": 26,    26: 26,
              "COBALT"       : 27,  "Co": 27,  "CO": 27,  "27": 27,    27: 27,
              "NICKEL"       : 28,  "Ni": 28,  "NI": 28,  "28": 28,    28: 28,
              "COPPER"       : 29,  "Cu": 29,  "CU": 29,  "29": 29,    29: 29,
              "ZINC"         : 30,  "Zn": 30,  "ZN": 30,  "30": 30,    30: 30,
              "GALLIUM"      : 31,  "Ga": 31,  "GA": 31,  "31": 31,    31: 31,
              "GERMANIUM"    : 32,  "Ge": 32,  "GE": 32,  "32": 32,    32: 32,
              "ARSENIC"      : 33,  "As": 33,  "AS": 33,  "33": 33,    33: 33,
              "SELENIUM"     : 34,  "Se": 34,  "SE": 34,  "34": 34,    34: 34,
              "BROMINE"      : 35,  "Br": 35,  "BR": 35,  "35": 35,    35: 35,
              "KRYPTON"      : 36,  "Kr": 36,  "KR": 36,  "36": 36,    36: 36,
              "RUBIDIUM"     : 37,  "Rb": 37,  "RB": 37,  "37": 37,    37: 37,
              "STRONTIUM"    : 38,  "Sr": 38,  "SR": 38,  "38": 38,    38: 38,
              "YTTRIUM"      : 39,   "Y": 39,   "Y": 39,  "39": 39,    39: 39,
              "ZIRCONIUM"    : 40,  "Zr": 40,  "ZR": 40,  "40": 40,    40: 40,
              "NIOBIUM"      : 41,  "Nb": 41,  "NB": 41,  "41": 41,    41: 41,
              "MOLYBDENUM"   : 42,  "Mo": 42,  "MO": 42,  "42": 42,    42: 42,
              "TECHNETIUM"   : 43,  "Tc": 43,  "TC": 43,  "43": 43,    43: 43,
              "RUTHENIUM"    : 44,  "Ru": 44,  "RU": 44,  "44": 44,    44: 44,
              "RHODIUM"      : 45,  "Rh": 45,  "RH": 45,  "45": 45,    45: 45,
              "PALLADIUM"    : 46,  "Pd": 46,  "PD": 46,  "46": 46,    46: 46,
              "SILVER"       : 47,  "Ag": 47,  "AG": 47,  "47": 47,    47: 47,
              "CADMIUM"      : 48,  "Cd": 48,  "CD": 48,  "48": 48,    48: 48,
              "INDIUM"       : 49,  "In": 49,  "IN": 49,  "49": 49,    49: 49,
              "TIN"          : 50,  "Sn": 50,  "SN": 50,  "50": 50,    50: 50,
              "ANTIMONY"     : 51,  "Sb": 51,  "SB": 51,  "51": 51,    51: 51,
              "TELLURIUM"    : 52,  "Te": 52,  "TE": 52,  "52": 52,    52: 52,
              "IODINE"       : 53,   "I": 53,   "I": 53,  "53": 53,    53: 53,
              "XENON"        : 54,  "Xe": 54,  "XE": 54,  "54": 54,    54: 54,
              "CESIUM"       : 55,  "Cs": 55,  "CS": 55,  "55": 55,    55: 55,
              "BARIUM"       : 56,  "Ba": 56,  "BA": 56,  "56": 56,    56: 56,
              "LANTHANUM"    : 57,  "La": 57,  "LA": 57,  "57": 57,    57: 57,
              "CERIUM"       : 58,  "Ce": 58,  "CE": 58,  "58": 58,    58: 58,
              "PRASEODYMIUM" : 59,  "Pr": 59,  "PR": 59,  "59": 59,    59: 59,
              "NEODYMIUM"    : 60,  "Nd": 60,  "ND": 60,  "60": 60,    60: 60,
              "PROMETHIUM"   : 61,  "Pm": 61,  "PM": 61,  "61": 61,    61: 61,
              "SAMARIUM"     : 62,  "Sm": 62,  "SM": 62,  "62": 62,    62: 62,
              "EUROPIUM"     : 63,  "Eu": 63,  "EU": 63,  "63": 63,    63: 63,
              "GADOLINIUM"   : 64,  "Gd": 64,  "GD": 64,  "64": 64,    64: 64,
              "TERBIUM"      : 65,  "Tb": 65,  "TB": 65,  "65": 65,    65: 65,
              "DYSPROSIUM"   : 66,  "Dy": 66,  "DY": 66,  "66": 66,    66: 66,
              "HOLMIUM"      : 67,  "Ho": 67,  "HO": 67,  "67": 67,    67: 67,
              "ERBIUM"       : 68,  "Er": 68,  "ER": 68,  "68": 68,    68: 68,
              "THULIUM"      : 69,  "TM": 69,  "TM": 69,  "69": 69,    69: 69,
              "YTTERBIUM"    : 70,  "Yb": 70,  "YB": 70,  "70": 70,    70: 70,
              "LUTETIUM"     : 71,  "Lu": 71,  "LU": 71,  "71": 71,    71: 71,
              "HAFNIUM"      : 72,  "Hf": 72,  "HF": 72,  "72": 72,    72: 72,
              "TANTALUM"     : 73,  "Ta": 73,  "TA": 73,  "73": 73,    73: 73,
              "TUNGSTEN"     : 74,   "W": 74,   "W": 74,  "74": 74,    74: 74,
              "RHENIUM"      : 75,  "Re": 75,  "RE": 75,  "75": 75,    75: 75,
              "OSMIUM"       : 76,  "Os": 76,  "OS": 76,  "76": 76,    76: 76,
              "IRIDIUM"      : 77,  "Ir": 77,  "IR": 77,  "77": 77,    77: 77,
              "PLATINUM"     : 78,  "Pt": 78,  "PT": 78,  "78": 78,    78: 78,
              "GOLD"         : 79,  "Au": 79,  "AU": 79,  "79": 79,    79: 79,
              "MERCURY"      : 80,  "Hg": 80,  "HG": 80,  "80": 80,    80: 80,
              "THALLIUM"     : 81,  "Tl": 81,  "TL": 81,  "81": 81,    81: 81,
              "LEAD"         : 82,  "Pb": 82,  "PB": 82,  "82": 82,    82: 82,
              "BISMUTH"      : 83,  "Bi": 83,  "BI": 83,  "83": 83,    83: 83,
              "POLONIUM"     : 84,  "Po": 84,  "PO": 84,  "84": 84,    84: 84,
              "ASTATINE"     : 85,  "At": 85,  "AT": 85,  "85": 85,    85: 85,
              "RADON"        : 86,  "Rn": 86,  "RN": 86,  "86": 86,    86: 86,
              }

        FullName=["HYDROGEN", "HELIUM", "LITHIUM", "BERYLLIUM", "BORON", "CARBON", "NITROGEN", "OXYGEN", "FLUORINE", "NEON", 
              "SODIUM", "MAGNESIUM", "ALUMINUM", "SILICON", "PHOSPHORUS", "SULFUR", "CHLORINE", "ARGON", "POTASSIUM", "CALCIUM", 
              "SCANDIUM", "TITANIUM", "VANADIUM", "CHROMIUM", "MANGANESE", "IRON", "COBALT", "NICKEL", "COPPER", "ZINC", 
              "GALLIUM", "GERMANIUM", "ARSENIC", "SELENIUM", "BROMINE", "KRYPTON", "RUBIDIUM", "STRONTIUM", "YTTRIUM", "ZIRCONIUM", 
              "NIOBIUM", "MOLYBDENUM", "TECHNETIUM", "RUTHENIUM", "RHODIUM", "PALLADIUM", "SILVER", "CADMIUM", "INDIUM", "TIN", 
              "ANTIMONY", "TELLURIUM", "IODINE", "XENON", "CESIUM", "BARIUM", "LANTHANUM", "CERIUM", "PRASEODYMIUM", "NEODYMIUM", 
              "PROMETHIUM", "SAMARIUM", "EUROPIUM", "GADOLINIUM", "TERBIUM", "DYSPROSIUM", "HOLMIUM", "ERBIUM", "THULIUM", "YTTERBIUM", 
              "LUTETIUM", "HAFNIUM", "TANTALUM", "TUNGSTEN", "RHENIUM", "OSMIUM", "IRIDIUM", "PLATINUM", "GOLD", "MERCURY", 
              "THALLIUM", "LEAD", "BISMUTH", "POLONIUM", "ASTATINE", "RADON"]

        Symbol=[ "H","He","Li","Be","B","C","N","O","F","Ne",
                "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
                "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
                "Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
                "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
                "Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
                "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","TM","Yb",
                "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
                "Tl","Pb","Bi","Po","At","Rn"]

        Mass=[1.008,4.003,6.941,9.012,10.811,12.011,14.007,15.999,18.998,20.180,
              22.990,24.305,26.982,28.086,30.974,32.065,35.453,39.948,39.098,40.078,
              44.956,47.867,50.942,51.996,54.938,55.845,58.933,58.693,63.546,65.390,
              69.723,72.640,74.922,78.960,79.904,83.800,85.468,87.620,88.906,91.224,
              92.906,95.940,98.000,101.070,102.906,106.420,107.868,112.411,114.818,118.710,
              121.760,127.600,126.905,131.293,132.906,137.327,138.906,140.116,140.908,144.240,
              145.000,150.360,151.964,157.250,158.925,162.500,164.930,167.259,168.934,173.040,
              174.967,178.490,180.948,183.840,186.207,190.230,192.217,195.078,196.967,200.590,
              204.383,207.200,208.980,209.000,210.000,222.000]

        # Van der Waals Radius, missing data replaced by 2.00
        Radii=[1.20,1.40,1.82,1.53,1.92,1.70,1.55,1.52,1.47,1.54,
               2.27,1.73,1.84,2.10,1.80,1.80,1.75,1.88,2.75,2.31,
               2.11,2.00,2.00,2.00,2.00,2.00,2.00,1.63,1.40,1.39,
               1.87,2.11,1.85,1.90,1.85,2.02,3.03,2.49,2.00,2.00,
               2.00,2.00,2.00,2.00,2.00,1.63,1.72,1.58,1.93,2.17,
               2.00,2.06,1.98,2.16,3.43,2.68,2.00,2.00,2.00,2.00,
               2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,
               2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.75,1.66,1.55,
               1.96,2.02,2.07,1.97,2.02,2.20]

        self.__name = Periodic_Table[name]
        self.__FullName = FullName[self.__name-1]
        self.__Symbol = Symbol[self.__name-1]
        self.__Mass = Mass[self.__name-1]
        self.__Radii = Radii[self.__name-1]

    def getFullName(self):
        return self.__FullName
    def getSymbol(self):
        return self.__Symbol
    def getUpperSymbol(self):
        return self.__Symbol.upper()
    def getMass(self):
        return self.__Mass
    def getNuc(self):
        return self.__name
    def getNelectron(self):
        return self.__name
    def getRadii(self):
        return self.__Radii

