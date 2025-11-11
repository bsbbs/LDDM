from neuron import n
n.nrnmpi_init()       # initialize MPI
pc = n.ParallelContext()
print('I am {} of {}'.format(pc.id(), pc.nhost()))
n.quit()              # necessary to avoid a warning message on parallel exit on some systems
