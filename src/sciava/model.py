from sciava.parameters import Parameters
from sciava.system import System


class Model:
    """ Contains all the information for the calculation.
        This includes the physical system itself and all parameters. """

    def __init__(self, name=None, params=None, system=None):
        self.name   = name   if name   is not None else 'Untitled'
        self.params = params if params is not None else Parameters()
        self.system = system if system is not None else System()

    def updateParams(self, **kwargs):
        self.params.update(**kwargs)

    def updateSystem(self, **kwargs):
        self.system.update(**kwargs)

    def __repr__(self):
        return str(self)

    def __str__(self):
        string = '-- {} --\n\n'.format(self.name)

        string += 'Parameters -->\n'
        for param, value in self.params.getSetParams().items():
            string += '  {:>18} : {:<18}\n'.format(param, value)

        # TODO: system.

        # TODO: results of runs.

        return string