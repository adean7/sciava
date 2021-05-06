from time import time

''' This is just the basic list of parameters. '''
parametersKnown = ['TASK',
                   'TASKMETHOD',
                   'THEORY', 'XCFUNCTIONAL',
                   'ATOMICSOLVER']

''' This dictionary takes many parameters and gives the version which is used in the code. '''
parametersManyToOne = { 'TASK' : 'TASK',

                        'TASKMETHOD' : 'TASKMETHOD',

                        'THEORY'       : 'THEORY',
                        'XCFUNCTIONAL' : 'THEORY',

                        'ATOMICSOLVER' : 'ATOMICSOLVER'
                        }

''' This dictionary holds all the possible values for each parameter. '''
parameterValues = { 'TASK'         : ['SP', 'SINGLEPOINT'],
                    'TASKMETHOD'   : ['ATOMISTIC', 'HF', 'HARTREE-FOCK', 'HARTREEFOCK', 'HARTREE_FOCK'],

                    'THEORY' : ['LDA', 'PBE', 'BLYP'],

                    'ATOMICSOLVER' : ['SH', 'SCHRODINGER', 'SCHROEDINGER', 'SCHRO', 'SCHROE', 'SCHROD', 'SCHROED']
                    }

''' This dictionary holds all the default values for parameters. There may be specified configurations that have separate 
defaults and they are coded manually in the setDefaults function of the Parameters class. '''
parameterDefaults = { 'TASK'         : 'SP',
                      'TASKMETHOD'   : 'HF',

                      'THEORY' : 'LDA',

                      'ATOMICSOLVER' : 'SH'
                      }

''' This dictionary holds the nice presentable names for printing parameters to screen. '''
parametersNiceNames = { 'TASK'       : 'Task',
                        'TASKMETHOD' : 'Task method',

                        'THEORY' : 'Theory',

                        'ATOMICSOLVER' : 'Atomic solver'
                        }

''' This dictionary takes many parameter values and gives the version which is used in the code. '''
valuesManyToOne = { 'TASK' : { 'SP'          : 'SP',
                               'SINGLEPOINT' : 'SP'
                               },

                    'TASKMETHOD' : { 'ATOMISTIC'    : 'ATOMISTIC',
                                     'HF'           : 'HF',
                                     'HARTREEFOCK'  : 'HF',
                                     'HARTREE-FOCK' : 'HF',
                                     'HARTREE_FOCK' : 'HF'
                                     },

                    'THEORY' : { 'LDA'  : 'LDA',
                                 'PBE'  : 'PBE',
                                 'BLYP' : 'BLYP'
                                 },

                    'ATOMICSOLVER' : { 'SH'           : 'SH',
                                       'SCHRODINGER'  : 'SH',
                                       'SCHROEDINGER' : 'SH',
                                       'SCHRO'        : 'SH',
                                       'SCHROE'       : 'SH',
                                       'SCHROD'       : 'SH',
                                       'SCHROED'      : 'SH'
                                       }
                    }

''' This dictionary holds the nice presentable names for printing parameter values to screen. '''
valuesNiceNames = { 'TASK' : { 'SP' : 'Singlepoint'
                               },

                    'TASKMETHOD' : { 'ATOMISTIC' : 'Atomistic',
                                     'HF'        : 'Hartree-Fock'
                                     },

                    'THEORY' : { 'LDA'  : 'LDA',
                                 'PBE'  : 'PBE',
                                 'BLYP' : 'BLYP'
                                 },

                    'ATOMICSOLVER' : { 'SH' : 'Schroedinger'
                                       }
                    }





def parameterAllowed(parameter):
    assert type(parameter) is str, '{} parameter must be specified by a string.'.format(parameter)

    param = parameter.upper()

    return True if param in parametersKnown else False

def valueAllowed(parameter, value):
    assert type(parameter) is str, '{} parameter must be specified by a string.'.format(parameter)

    param = parameter.upper()
    val = value.upper() if type(value) is str else value

    return True if parameterAllowed(param) and val in parameterValues.get(param) else False

def getShortParam(parameter):
    assert type(parameter) is str, '{} parameter must be specified by a string.'.format(parameter)

    param = parameter.upper()

    assert parameterAllowed(param), '{} parameter not known.'.format(parameter)

    return parametersManyToOne.get(param)

def getShortValue(parameter, value):
    assert type(parameter) is str, '{} parameter must be specified by a string.'.format(parameter)
    assert value is not None, '{} parameter value must be specified.'.format(parameter)

    param = parameter.upper()
    val = value.upper() if type(value) is str else value

    assert parameterAllowed(param), '{} parameter not known.'.format(parameter)
    assert valueAllowed(param, val), '{} not an acceptable value for {}'.format(value, parameter)

    return valuesManyToOne.get(param).get(val)

def getDefaultParamValue(parameter):
    assert type(parameter) is str, '{} parameter must be specified by a string.'.format(parameter)

    param = parameter.upper()

    assert parameterAllowed(param), '{} parameter not known.'.format(parameter)

    return parameterDefaults.get(param)

def getNiceParamName(parameter):
    assert type(parameter) is str, '{} parameter must be specified by a string.'.format(parameter)

    param = parameter.upper()

    assert parameterAllowed(param), '{} parameter not known.'.format(parameter)

    return parametersNiceNames.get(param)

def getNiceValueName(parameter, value):
    assert type(parameter) is str, '{} parameter must be specified by a string.'.format(parameter)
    assert value is not None, '{} parameter value must be specified.'.format(parameter)

    param = parameter.upper()
    val = value.upper() if type(value) is str else value

    assert parameterAllowed(param), '{} parameter not known.'.format(parameter)
    assert valueAllowed(param, val), '{} not an acceptable value for {}'.format(value, parameter)

    return valuesNiceNames.get(param).get(val)





class Parameters:
    def __init__(self, currentSystem=None, **kwargs):
        self.task       = None
        self.taskMethod = None

        self.theory = None

        self.atomicSolver = None

        # Record what parameters have been specified by the user so they are not updated when deciding default values.
        self.userSpecified = []

        self.update(currentSystem, **kwargs)

        self.startTime = None
        self.stopTime  = None

    def update(self, currentSystem=None, **kwargs):
        """ This function updates n parameters with associated values. """

        for parameter, value in kwargs.items():
            if value is None:
                self.remove(parameter, currentSystem)
                continue

            param = getShortParam(parameter)
            val = getShortValue(param, value)

            self.userSpecified.append(param)

            if param == getShortParam('TASK'):
                self.task = val

            elif param == getShortParam('TASKMETHOD'):
                self.taskMethod = val

            elif param == getShortParam('THEORY'):
                self.theory = val

            elif param == getShortParam('ATOMICSOLVER'):
                self.atomicSolver = val

        # Get defaults of system first.
        if currentSystem is not None:
            currentSystem.getDefaults(self)
        self.getDefaults(currentSystem)

    def remove(self, parameter=None, currentSystem=None):
        """ This function removes a single parameter from the configuration. """

        assert parameter is not None, 'Need to supply parameter to remove.'

        param = getShortParam(parameter)

        if param in self.userSpecified:
            self.userSpecified.remove(param)

        if param == getShortParam('TASK'):
            self.task = None

        elif param == getShortParam('TASKMETHOD'):
            self.taskMethod = None

        elif param == getShortParam('THEORY'):
            self.theory = None

        elif param == getShortParam('ATOMICSOLVER'):
            self.atomicSolver = None

        # Get defaults of system first.
        if currentSystem is not None:
            currentSystem.getDefaults(self)
        self.getDefaults(currentSystem)

    def getDefaults(self, currentSystem=None):
        """ This function gets default values for parameters based off the current configuration. """

        # If task hasn't been set or if it's set but not by the user then we can get a default.
        if self.task is None or getShortParam('TASK') not in self.userSpecified:
            self.task = getDefaultParamValue('TASK')

        if self.task == getShortValue('TASK', 'SINGLEPOINT'):

            # If task method hasn't been set or if it's set but not by the user we can give a default.
            if self.taskMethod is None or getShortParam('TASKMETHOD') not in self.userSpecified:

                # Special case when we have one atom (therefore spherically symmetric) we default to a full atomistic run.
                if currentSystem is not None and currentSystem.numAtoms == 1:
                    self.taskMethod = getShortValue('TASKMETHOD', 'ATOMISTIC')

                else:   # If no special case then let's just get the usual default for task method.
                    self.taskMethod = getDefaultParamValue('TASKMETHOD')

            # If theory hasn't been set or if it's set but not by the user we can give a default.
            if self.theory is None or getShortParam('THEORY') not in self.userSpecified:

                # We only need the DFT theory if we're doing an atomistic calculation, otherwise it's redundant.
                if self.taskMethod == 'ATOMISTIC':
                    self.theory = getDefaultParamValue('THEORY')

            # If atomic solver hasn't been set or if it's set but not by the user we can give a default.
            if self.atomicSolver is None or getShortParam('ATOMICSOLVER') not in self.userSpecified:

                # We only need the atomic solver if we're doing an atomistic calculation, otherwise it's redundant.
                if self.taskMethod == 'ATOMISTIC':
                    self.atomicSolver = getDefaultParamValue('ATOMICSOLVER')

    def getCurrentParams(self, currentSystem=None):
        """ The function returns a dict of parameters that are currently set. """
        dct = dict()

        dct[getNiceParamName('TASK')] = getNiceValueName('TASK', self.task)

        # Only output certain information for singlepoint tasks.
        if self.task == getShortValue('TASK', 'SINGLEPOINT'):
            dct[getNiceParamName('TASKMETHOD')] = getNiceValueName('TASKMETHOD', self.taskMethod)

            # Output some extra information if we're doing an atomistic calculation.
            if self.taskMethod == getShortValue('TASKMETHOD', 'ATOMISTIC'):
                dct[getNiceParamName('THEORY')] = getNiceValueName('THEORY', self.theory)
                dct[getNiceParamName('ATOMICSOLVER')] = getNiceValueName('ATOMICSOLVER', self.atomicSolver)

        return dct

    def startTimer(self):
        self.startTime = time()

    def stopTimer(self):
        self.stopTime = time()

    def getRunTime(self):
        return self.stopTime - self.startTime
