
import sciava


model = sciava.Model()

model.updateSystem(atomcoord=['h', 1.0, 3.0, 2.0])
model.updateParams(atomicsolver='schrod', xcfunctional='lda')

print(model)

model.run()