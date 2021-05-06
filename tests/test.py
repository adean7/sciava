
import sciava


model = sciava.Model()

model.updateSystem(atomcoord=['He', 1.0, 3.0, 2.0])
model.updateParams(taskmethod='atomistic', atomicsolver='schrod', xcfunctional='lda')

print(model)

model.run()
