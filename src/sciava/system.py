import numpy as np

from sciava.data import *


knownCells = ['ATOMCOORD', 'ATOMICCOORD', 'ATOMCOORDINATE', 'ATOMICCOORDINATE', 'ATOM_COORD', 'ATOMIC_COORD',
              'ATOM_COORDINATE', 'ATOMIC_COORDINATE', 'ATOMCOORDS', 'ATOMICCOORDS', 'ATOMCOORDINATES', 'ATOMICCOORDINATES',
              'ATOM_COORDS', 'ATOMIC_COORDS', 'ATOM_COORDINATES', 'ATOMIC_COORDINATES']

cellValues = { }

cellDefaults = { 'NUMATOMS'   : 0
                 }

class System:
    def __init__(self, currentParams=None, **kwargs):
        self.atoms    = []
        self.numAtoms = 0

        # Record what cells have been specified by the user so they are not updated when deciding default values.
        self.userSpecified = []

        self.update(currentParams, **kwargs)
        self.check()

    def update(self, currentParams=None, **kwargs):
        """ This function updates n things about the system with associated values. """

        for cell, value in kwargs.items():
            c = cell.upper()
            assert c in knownCells, '{} cell not known.'.format(cell)

            if value is None:
                if c in self.userSpecified:
                    self.userSpecified.remove(c)

            # We do value checking manually in system as there are so many possibilities.

            if c in ['ATOMCOORD', 'ATOMICCOORD', 'ATOMCOORDINATE', 'ATOMICCOORDINATE', 'ATOM_COORD', 'ATOMIC_COORD',
              'ATOM_COORDINATE', 'ATOMIC_COORDINATE' 'ATOMCOORDS', 'ATOMICCOORDS', 'ATOMCOORDINATES', 'ATOMICCOORDINATES',
              'ATOM_COORDS', 'ATOMIC_COORDS', 'ATOM_COORDINATES', 'ATOMIC_COORDINATES']:
                assert type(value) in [Atom, list], '{} value not acceptable for {}'.format(value, cell)

                if type(value) is list:
                    assert len(value) == 4, 'List entry of atomic coordinate must be [element, x, y, z] format.'
                    atom = Atom(value[0], value[1], value[2], value[3])
                else:
                    atom = value

                self.atoms.append(atom)
                self.numAtoms = len(self.atoms)
                self.userSpecified.append('ATOMS')

        # Get defaults of system first.
        self.getDefaults(currentParams)
        if currentParams is not None:
            currentParams.getDefaults(self)

    def remove(self, cell=None, currentParams=None):
        """ This function removes a single cell from the configuration. """

        if cell is None:
            raise ValueError('Need to supply cell to remove.')

        assert type(cell) is str, 'Cell to remove must be specified by a string.'

        c = cell.upper()
        assert c in knownCells, '{} cell not known.'.format(cell)

        if c in ['ATOMCOORD', 'ATOMICCOORD', 'ATOMCOORDINATE', 'ATOMICCOORDINATE', 'ATOM_COORD', 'ATOMIC_COORD',
              'ATOM_COORDINATE', 'ATOMIC_COORDINATE', 'ATOMCOORDS', 'ATOMICCOORDS', 'ATOMCOORDINATES', 'ATOMICCOORDINATES',
              'ATOM_COORDS', 'ATOMIC_COORDS', 'ATOM_COORDINATES', 'ATOMIC_COORDINATES']:
            self.atoms = []
            self.numAtoms = 0
            if 'ATOMS' in self.userSpecified:
                self.userSpecified.remove('ATOMS')

        # Get defaults of system first.
        self.getDefaults(currentParams)
        if currentParams is not None:
            currentParams.getDefaults(self)

    def getDefaults(self, currentParams=None):
        """ This function gets default values for cells based off the current configuration. """
        pass

    def check(self):
        """ This function checks that the current system is logical and the model will run. """
        pass

    def getCurrentSystem(self, currentParams=None):
        """ The function returns a dict of system properties that are currently set. """
        dct = dict()

        if currentParams is None:
            dct['Atoms']        = ', '.join([e.element for e in self.atoms])
            dct['Num of atoms'] = self.numAtoms

        elif currentParams.task == 'SINGLEPOINT':
            dct['Atoms']        = ', '.join([e.element for e in self.atoms])
            dct['Num of atoms'] = self.numAtoms

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