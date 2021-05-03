from sciava.basis import AllElectronBasis
from sciava.scf import AllElectronAtom, atomisticCycleSchrod


def run(model):
    print('Beginning quantum calculation...')

    if model.params.taskMethod == 'ATOMISTIC':
        atomisticRun(model.params, model.system)

    print('Cycle complete!')


def atomisticRun(params, system):
    atom = system.atoms[0]

    # Generate the basis for the calculation.
    ab = AllElectronBasis(atom.Z)

    # Initialise the all electron atom.
    aeat = AllElectronAtom(params, atom.element, atom.Z, ab)

    # Do the all electron calculation for the atom.
    if params.atomicSolver == 'SH':
        atomisticCycleSchrod(ab, aeat)
