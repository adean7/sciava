
import sciava


model = sciava.Model()

atom = sciava.Atom(element='Mt',
                   x=0.0, y=0.0, z=0.0)

model.updateSystem(atomcoord=atom)
model.updateParams(theory='blyp')

print(model)

model.run()