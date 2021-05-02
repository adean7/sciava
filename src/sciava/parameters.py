from sciava.data import *


class Parameters:
    def __init__(self, **kwargs):
        self.task = None

        self.update(**kwargs)
        self.getDefaults()
        self.check()

    def update(self, **kwargs):
        """ This function updates n parameters with associated values. """

        for parameter, value in kwargs.items():
            param = parameter.lower()
            val   = value.lower()

            assert param in knownParameters, '{} parameter not known.'.format(parameter)

            assert val in parameterValues.get(param), '{} not an acceptable value for {}'.format(value, parameter)

            if param == 'task':
                self.task = val

    def getDefaults(self):
        """ This function gets default values for parameters based off the current configuration. """

        self.task = self.task if self.task is not None else parameterDefaults.get('task')

    def check(self):
        """ This function checks that the current parameters of the system are logical and the model will run. """
        pass

    def getSetParams(self):
        """ The function returns a dict of parameters that are currently set. """
        return { 'task' : self.task }