from sciava.parameters import Parameters
from sciava.system import System

from sciava.qm import run as qmRun


class Model:
    """ Contains all the information for the calculation.
        This includes the physical system itself and all parameters. """

    def __init__(self, name=None, params=None, system=None):
        self.name   = name   if name   is not None else '<untitled-model>'
        self.params = params if params is not None else Parameters()
        self.system = system if system is not None else System()

    def updateParams(self, **kwargs):
        self.params.update(currentSystem=self.system, **kwargs)

    def updateSystem(self, **kwargs):
        self.system.update(currentParams=self.params, **kwargs)

    def removeParam(self, parameter):
        self.params.remove(parameter, currentSystem=self.system)

    def removeSystem(self, cell):
        self.system.remove(cell, currentParams=self.params)

    def run(self):
        self.params.startTimer()

        if self.params.task == 'SP':
            qmRun(self)

        self.params.stopTimer()

        print('\nTime of run: {:>9.3f} s.'.format(self.params.getRunTime()))


    def __repr__(self):
        return str(self)

    def __str__(self):
        string = '-- {} --\n'.format(self.name)

        currentParams = self.params.getCurrentParams(currentSystem=self.system).items()
        if len(currentParams) > 0:
            string += '\nParameters -->\n'
            for param, value in currentParams:
                if value is not None:
                    string += '  {:>18} : {:<18}\n'.format(param, value)

        currentSystem = self.system.getCurrentSystem(currentParams=self.params).items()
        if any([True for _, value in currentSystem if value is not None]):
            string += '\nSystem -->\n'
            for cell, value in currentSystem:
                if value is not None:
                    string += '  {:>18} : {:<18}\n'.format(cell, value)

        # TODO: results of runs.

        return string