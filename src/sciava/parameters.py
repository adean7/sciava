
knownParameters = ['TASK', 'TASKMETHOD']

parameterValues = { 'TASK'       : ['SP', 'SINGLEPOINT'],
                    'TASKMETHOD' : ['ATOMISTIC', 'HF', 'HARTREE-FOCK', 'HARTREEFOCK', 'HARTREE_FOCK']
                    }

parameterDefaults = { 'TASK'       : 'SINGLEPOINT',
                      'TASKMETHOD' : 'HARTREE-FOCK'
                      }

class Parameters:
    def __init__(self, currentSystem=None, **kwargs):
        self.task       = None
        self.taskMethod = None

        # Record what parameters have been specified by the user so they are not updated when deciding default values.
        self.userSpecified = []

        self.update(currentSystem, **kwargs)
        self.check()

    def update(self, currentSystem=None, **kwargs):
        """ This function updates n parameters with associated values. """

        for parameter, value in kwargs.items():
            param = parameter.upper()
            assert param in knownParameters, '{} parameter not known.'.format(parameter)

            if value is None:
                val = None
                if param in self.userSpecified:
                    self.userSpecified.remove(param)
            else:
                val = value.upper()
                assert val in parameterValues.get(param), '{} not an acceptable value for {}'.format(value, parameter)

            if param == 'TASK':
                self.userSpecified.append('TASK')
                if val in ['SP', 'SINGLEPOINT']:
                    self.task = 'SINGLEPOINT'

            elif param == 'TASKMETHOD':
                self.userSpecified.append('TASKMETHOD')
                if val == 'ATOMISTIC':
                    self.taskMethod = 'ATOMISTIC'
                elif val in ['HF', 'HARTREE-FOCK', 'HARTREEFOCK', 'HARTREE_FOCK']:
                    self.taskMethod = 'HARTREE-FOCK'

        # Get defaults of system first.
        if currentSystem is not None:
            currentSystem.getDefaults(self)
        self.getDefaults(currentSystem)

    def remove(self, parameter=None, currentSystem=None):
        """ This function removes a single parameter from the configuration. """

        if parameter is None:
            raise ValueError('Need to supply parameter to remove.')

        assert type(parameter) is str, 'Parameter to remove must be specified by a string.'

        param = parameter.upper()
        assert param in knownParameters, '{} parameter not known.'.format(parameter)

        if param == 'TASK':
            self.task = None
            if 'TASK' in self.userSpecified:
                self.userSpecified.remove('TASK')

        elif param == 'TASKMETHOD':
            self.taskMethod = None
            if 'TASKMETHOD' in self.userSpecified:
                self.userSpecified.remove('TASKMETHOD')

        # Get defaults of system first.
        if currentSystem is not None:
            currentSystem.getDefaults(self)
        self.getDefaults(currentSystem)

    def getDefaults(self, currentSystem=None):
        """ This function gets default values for parameters based off the current configuration. """

        # Task.
        if self.task is None:
            self.task = parameterDefaults.get('TASK')

        # Task method.
        if self.task == 'SINGLEPOINT':
            if self.taskMethod is None:
                if currentSystem is not None and currentSystem.numAtoms == 1:
                    self.taskMethod = 'ATOMISTIC'
                else:
                    self.taskMethod = parameterDefaults.get('TASKMETHOD')
            elif 'TASKMETHOD' not in self.userSpecified:
                if currentSystem is not None and currentSystem.numAtoms == 1:
                    self.taskMethod = 'ATOMISTIC'
                else:
                    self.taskMethod = parameterDefaults.get('TASKMETHOD')

    def check(self):
        """ This function checks that the current parameters of the system are logical and the model will run. """
        pass

    def getCurrentParams(self, currentSystem=None):
        """ The function returns a dict of parameters that are currently set. """
        dct = dict()

        if self.task == 'SINGLEPOINT':
            dct['Task']        = 'Singlepoint'
            dct['Task method'] = self.taskMethod[0].upper() + self.taskMethod[1:].lower()

        return dct