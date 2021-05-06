import numpy as np

from sciava.data import *


def checkAtom(atom):
    return True if (type(atom) is Atom or (type(atom) is list and len(atom) == 4)) else False

def getAtom(atom):
    return Atom(atom[0], atom[1], atom[2], atom[3]) if type(atom) is list else atom

def defaultAtom():
    return None

def niceAtom(atom):
    return str(atom)


def checkAtoms(atoms):
    if type(atoms) is not list:
        return False

    for atom in atoms:
        if not checkAtom(atom):
            return False

    return True

def getAtoms(atoms):
    return [Atom(atom[0], atom[1], atom[2], atom[3]) if type(atom) is list else atom for atom in atoms]

def defaultAtoms():
    return list()

def niceAtoms(atoms):
    return ', '.join([atom.element for atom in atoms])


def checkNumAtoms(numAtoms):
    return True if type(numAtoms) is int and numAtoms >= 0 else False

def getNumAtoms(numAtoms):
    return numAtoms

def defaultNumAtoms():
    return int(0)

def niceNumAtoms(numAtoms):
    return '{:<3d}'.format(numAtoms)







def cellAllowed(cell):
    assert type(cell) is str, '{} cell must be specified by a string.'.format(cell)

    c = cell.upper()

    return True if c in cellsKnown else False

def valueAllowed(cell, value):
    assert type(cell) is str, '{} cell must be specified by a string.'.format(cell)

    c = cell.upper()
    val = value.upper() if type(value) is str else value

    if not cellAllowed(c):
        return False

    funcToCheckVal = cellCheckValueFuncs.get(c)

    return funcToCheckVal(val)

def getShortCell(cell):
    assert type(cell) is str, '{} cell must be specified by a string.'.format(cell)

    c = cell.upper()

    assert cellAllowed(c), '{} parameter not known.'.format(cell)

    return cellsManyToOne.get(c)

def getShortValue(cell, value):
    assert type(cell) is str, '{} cell must be specified by a string.'.format(cell)
    assert value is not None, '{} cell value must be specified.'.format(cell)

    c = cell.upper()
    val = value.upper() if type(value) is str else value

    assert cellAllowed(c), '{} cell not known.'.format(cell)
    assert valueAllowed(c, val), '{} not an acceptable value for {}'.format(value, cell)

    funcToGetVal = cellGetValueFuncs.get(c)

    return funcToGetVal(val)

def getDefaultCellValue(cell):
    assert type(cell) is str, '{} cell must be specified by a string.'.format(cell)

    c = cell.upper()

    assert cellAllowed(c), '{} cell not known.'.format(cell)

    funcToGetDefault = cellGetDefaultFuncs.get(c)

    return funcToGetDefault()

def getNiceCellName(cell):
    assert type(cell) is str, '{} cell must be specified by a string.'.format(cell)

    c = cell.upper()

    assert cellAllowed(c), '{} cell not known.'.format(cell)

    return cellsNiceNames.get(c)

def getNiceValue(cell, value):
    assert type(cell) is str, '{} cell must be specified by a string.'.format(cell)
    assert value is not None, '{} cell value must be specified.'.format(cell)

    c = cell.upper()
    val = value.upper() if type(value) is str else value

    assert cellAllowed(c), '{} cell not known.'.format(cell)
    assert valueAllowed(c, val), '{} not an acceptable value for {}'.format(value, cell)

    funcToGetNiceNames = valuesNice.get(c)

    return funcToGetNiceNames(val)





''' This is just the basic list of cells. '''
cellsKnown = ['ATOM', 'ATOMCOORD', 'ATOMICCOORD', 'ATOMCOORDINATE', 'ATOMICCOORDINATE',
              'ATOM_COORD', 'ATOMIC_COORD', 'ATOM_COORDINATE', 'ATOMIC_COORDINATE',

              'ATOMS', 'ATOMCOORDS', 'ATOMICCOORDS', 'ATOMCOORDINATES', 'ATOMICCOORDINATES',
              'ATOM_COORDS', 'ATOMIC_COORDS', 'ATOM_COORDINATES', 'ATOMIC_COORDINATES',

              'NUMATOMS', 'NUMBERATOMS', 'NUMOFATOMS', 'NUMBEROFATOMS',
              'NUM_ATOMS', 'NUMBER_ATOMS', 'NUM_OF_ATOMS', 'NUMBER_OF_ATOMS'
              ]

''' This dictionary takes many cells and gives the version which is used in the code. '''
cellsManyToOne = { 'ATOM'               : 'ATOM',
                   'ATOMCOORD'          : 'ATOM',
                   'ATOMICCOORD'        : 'ATOM',
                   'ATOMCOORDINATE'     : 'ATOM',
                   'ATOMICCOORDINATE'   : 'ATOM',
                   'ATOM_COORD'         : 'ATOM',
                   'ATOMIC_COORD'       : 'ATOM',
                   'ATOM_COORDINATE'    : 'ATOM',
                   'ATOMIC_COORDINATE'  : 'ATOM',

                   'ATOMS'              : 'ATOMS',
                   'ATOMCOORDS'         : 'ATOMS',
                   'ATOMICCOORDS'       : 'ATOMS',
                   'ATOMCOORDINATES'    : 'ATOMS',
                   'ATOMICCOORDINATES'  : 'ATOMS',
                   'ATOM_COORDS'        : 'ATOMS',
                   'ATOMIC_COORDS'      : 'ATOMS',
                   'ATOM_COORDINATES'   : 'ATOMS',
                   'ATOMIC_COORDINATES' : 'ATOMS',

                   'NUMATOMS'        : 'NUMATOMS',
                   'NUMBERATOMS'     : 'NUMATOMS',
                   'NUMOFATOMS'      : 'NUMATOMS',
                   'NUMBEROFATOMS'   : 'NUMATOMS',
                   'NUM_ATOMS'       : 'NUMATOMS',
                   'NUMBER_ATOMS'    : 'NUMATOMS',
                   'NUM_OF_ATOMS'    : 'NUMATOMS',
                   'NUMBER_OF_ATOMS' : 'NUMATOMS'
                   }

''' This dictionary holds the functions that check the value of the passed cell. '''
cellCheckValueFuncs = { 'ATOM'     : checkAtom,
                        'ATOMS'    : checkAtoms,
                        'NUMATOMS' : checkNumAtoms
                        }

''' This dictionary holds the functions that get the value of the passed cell. '''
cellGetValueFuncs = { 'ATOM'     : getAtom,
                      'ATOMS'    : getAtoms,
                      'NUMATOMS' : getNumAtoms
                      }

''' This dictionary holds the functions that get the default values for cells. '''
cellGetDefaultFuncs = { 'ATOM'     : defaultAtom,
                        'ATOMS'    : defaultAtoms,
                        'NUMATOMS' : defaultNumAtoms
                        }

''' This dictionary holds the nice presentable names for printing cells to screen. '''
cellsNiceNames = { 'ATOM'     : 'Atom',
                   'ATOMS'    : 'Atoms',
                   'NUMATOMS' : 'Num. atoms'
                   }

'''
'''''' This dictionary takes many cell values and gives the version which is used in the code. ''''''
valuesManyToOne = { 'TASK' : { 'SP'          : 'SP',
                               'SINGLEPOINT' : 'SP'
                               }
                    }
'''

''' This dictionary holds the functions that give the nice presentable forms for printing cell values to screen. '''
valuesNice = { 'ATOM'     : niceAtom,
               'ATOMS'    : niceAtoms,
               'NUMATOMS' : niceNumAtoms,
               }






class System:
    def __init__(self, currentParams=None, **kwargs):
        self.atoms    = []
        self.numAtoms = 0

        # Record what cells have been specified by the user so they are not updated when deciding default values.
        self.userSpecified = []

        self.update(currentParams, **kwargs)

    def update(self, currentParams=None, **kwargs):
        """ This function updates n things about the system with associated values. """

        for cell, value in kwargs.items():
            if value is None:
                self.remove(cell, currentParams)
                continue

            c = getShortCell(cell)
            val = getShortValue(c, value)

            self.userSpecified.append(c)

            if c == getShortCell('ATOM'):
                self.atoms.append(val)
                self.numAtoms = len(self.atoms)

                # Because ATOM and ATOMS correspond to the same list (self.atoms) but are treated separately:
                #    let's just sort out the self.userSpecified manually to be sure we have no issues
                if c in self.userSpecified:
                    self.userSpecified.remove(c)
                self.userSpecified.append(getShortCell('ATOMS'))

            elif c == getShortCell('ATOMS'):
                self.atoms += val
                self.numAtoms = len(self.atoms)

                # Because ATOM and ATOMS correspond to the same list (self.atoms) but are treated separately:
                #    let's just sort out the self.userSpecified manually to be sure we have no issues
                if c in self.userSpecified:
                    self.userSpecified.remove(c)
                self.userSpecified.append(getShortCell('ATOMS'))

        # Get defaults of system first.
        self.getDefaults(currentParams)
        if currentParams is not None:
            currentParams.getDefaults(self)

    def remove(self, cell=None, currentParams=None):
        """ This function removes a single cell from the configuration. """

        assert cell is not None, 'Need to supply cell to remove.'

        c = getShortCell(cell)

        if c in self.userSpecified:
            self.userSpecified.remove(c)

        if c == getShortCell('ATOM'):
            raise NotImplementedError('Removing single atoms currently not implemented.')
        elif c == getShortCell('ATOMS'):
            self.atoms = None
            self.numAtoms = None

        # Get defaults of system first.
        self.getDefaults(currentParams)
        if currentParams is not None:
            currentParams.getDefaults(self)

    def getDefaults(self, currentParams=None):
        """ This function gets default values for cells based off the current configuration. """

        if self.atoms is None or getShortCell('ATOMS') not in self.userSpecified:
            self.atoms = getDefaultCellValue('ATOMS')
            self.numAtoms = len(self.atoms)

    def getCurrentSystem(self, currentParams=None):
        """ The function returns a dict of system properties that are currently set. """
        dct = dict()

        if currentParams is None:
            dct[getNiceCellName('ATOMS')] = getNiceValue('ATOMS', self.atoms)
            dct[getNiceCellName('NUMATOMS')] = getNiceValue('NUMATOMS', self.numAtoms)

        # !! *** Be careful to ensure the if statements here use the shortParam and shortValue in parameters. *** !! #
        elif currentParams.task == 'SP':
            if self.numAtoms > 0:
                dct[getNiceCellName('ATOMS')] = getNiceValue('ATOMS', self.atoms)
                dct[getNiceCellName('NUMATOMS')] = getNiceValue('NUMATOMS', self.numAtoms)

        return dct


class Atom:
    def __init__(self, element=None, x=None, y=None, z=None):
        if element is None:
            raise ValueError('Atom requires an element.')

        if x is None or y is None or z is None:
            raise ValueError('Atom requires a 3D coordinate.')

        assert type(element) is str, 'Element {} must be a string.'.format(element)
        assert element[0].upper() + element[1:].lower() in periodicTable, 'Element {} not recognised.'.format(element)

        for coord in [x, y, z]:
            if type(coord) not in [float, int]:
                raise ValueError('{} not an acceptable positional value'.format(coord))

        self.element = element[0].upper() + element[1:].lower()
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.r = np.array((self.x, self.y, self.z))
        self.Z = atomicNumbers.get(self.element)

    def __repr__(self):
        return str(self)

    def __str__(self):
        return '{:>4}  {:>9.5f}  {:>9.5f}  {:>9.5f}'.format(self.element, self.x, self.y, self.z)

