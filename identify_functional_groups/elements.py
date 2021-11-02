from collections import defaultdict


# Set default value for the defaultdict used for functional group pattern lookup. If no matching element name is found,
# return 'unknown element'.
def def_value():
    return "unknown element!"


# Only heteroatoms are relevant for lookup, so H and C are omitted from the dictionary. Lanthanoids and Actinoids not
# included.
def get_element_names():
    element_names = defaultdict(def_value, {
        'He': 'HELIUM',
        'Li': 'LITHIUM',
        'Be': 'BERYLLIUM',
        'B': 'BORON',
        'N': 'NITROGEN',
        'O': 'OXYGEN',
        'F': 'FLUORINE',
        'Ne': 'NEON',
        'Na': 'SODIUM',
        'Mg': 'MAGNESIUM',
        'Al': 'ALUMINIUM',
        'Si': 'SILICON',
        'P': 'PHOSPHORUS',
        'S': 'SULFUR',
        'Cl': 'CHLORINE',
        'Ar': 'ARGON',
        'K': 'POTASSIUM',
        'Ca': 'CALCIUM',
        'Sc': 'SCANDIUM',
        'Ti': 'TITANIUM',
        'V': 'VANADIUM',
        'Cr': 'CHROMIUM',
        'Mn': 'MANGANESE',
        'Fe': 'IRON',
        'Co': 'COBALT',
        'Ni': 'NICKEL',
        'Cu': 'COPPER',
        'Zn': 'ZINC',
        'Ga': 'GALLIUM',
        'Ge': 'GERMANIUM',
        'As': 'ARSENIC',
        'Se': 'SELENIUM',
        'Br': 'BROMIDE',
        'Kr': 'KRYPTON',
        'Rb': 'RUBIDIUM',
        'Sr': 'STRONTIUM',
        'Y': 'YTTRIUM',
        'Zr': 'ZIRCONIUM',
        'Nb': 'NIOBIUM',
        'Mo': 'MOLYBDENUM',
        'Tc': 'TECHNETIUM',
        'Ru': 'Ruthenium',
        'Rh': 'RHODIUM',
        'Pd': 'PALLADIUM',
        'Ag': 'SILVER',
        'Cd': 'CADMIUM',
        'In': 'INDIUM',
        'Sn': 'TIN',
        'Sb': 'ANTIMONY',
        'Te': 'TELLURIUM',
        'I': 'IODINE',
        'Xe': 'XENON',
        'Cs': 'CESIUM',
        'Ba': 'BARIUM',
        'Hf': 'HAFNIUM',
        'Ta': 'TANTALUM',
        'W': 'TUNGSTEN',
        'Re': 'RHENIUM',
        'Os': 'OSMIUM',
        'Ir': 'IRIDIUM',
        'Pt': 'PLATINUM',
        'Au': 'GOLD',
        'Hg': 'MERCURY',
        'Tl': 'THALLIUM',
        'Pb': 'LEAD',
        'Bi': 'BISMUTH',
        'Po': 'POLONIUM',
        'At': 'ASTATINE',
        'Rn': 'RADON',
        'Fr': 'FRANCIUM',
        'Ra': 'RADIUM',
        'Rf': 'RUTHERFORDIUM',
        'Db': 'DUBNIUM',
        'Sg': 'SEABORGIUM',
        'Bh': 'BOHRIUM',
        'Hs': 'HASSIUM',
        'Mt': 'MEITNERIUM',
        'Ds': 'DARMSTADTIUM',
        'Rg': 'ROENTGENIUM',
        'Cn': 'COPERNICIUM',
        'Nh': 'NIHONIUM',
        'Fl': 'FLEROVIUM',
        'Mc': 'MOSCOVIUM',
        'Lv': 'LIVERMORIUM',
        'Ts': 'TENNESSINE',
        'Og': 'OGANESSON'
    })
    return element_names
