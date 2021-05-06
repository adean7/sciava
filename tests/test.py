import sciava

model = sciava.Model()

flourine = sciava.Atom('F', 1.0, 2.0, 3.0)

model.updateSystem(atomic_coord=flourine)
model.updateParams(taskmethod='hf', atomicsolver='schrod', xcfunctional='lda')

print(model)

#model.run()
