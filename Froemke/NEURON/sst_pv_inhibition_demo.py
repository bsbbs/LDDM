from neuron import h
import matplotlib.pyplot as plt

h.load_file('stdrun.hoc')
h.dt = 0.025
h.tstop = 100

### Create pyramidal cell with compartments ###
class PyramidalCell:
    def __init__(self):
        self.soma = h.Section(name='soma')
        self.soma.L = self.soma.diam = 20
        self.soma.insert('hh')

        self.dend = h.Section(name='dend')
        self.dend.L = 300
        self.dend.diam = 2
        self.dend.insert('pas')
        self.dend.connect(self.soma(1))

        self.v_soma = h.Vector().record(self.soma(0.5)._ref_v)
        self.v_dend = h.Vector().record(self.dend(1.0)._ref_v)

cell = PyramidalCell()

### Excitatory input to dendrite ###
exc_syn = h.ExpSyn(cell.dend(1.0))
exc_syn.e = 0
exc_syn.tau = 2
exc_stim = h.NetStim()
exc_stim.number = 1
exc_stim.start = 20
exc_stim.noise = 0
nc_exc = h.NetCon(exc_stim, exc_syn)
nc_exc.weight[0] = 0.01

### Inhibitory synapse for PV (soma) ###
pv_syn = h.ExpSyn(cell.soma(0.5))
pv_syn.e = -75
pv_syn.tau = 5
pv_stim = h.NetStim()
pv_stim.number = 1
pv_stim.start = 25  # shortly after excitation
pv_stim.noise = 0
nc_pv = h.NetCon(pv_stim, pv_syn)
nc_pv.weight[0] = 0.01

### Inhibitory synapse for SST (distal dendrite) ###
sst_syn = h.ExpSyn(cell.dend(1.0))
sst_syn.e = -75
sst_syn.tau = 5
sst_stim = h.NetStim()
sst_stim.number = 1
sst_stim.start = 25  # same timing as PV
sst_stim.noise = 0
nc_sst = h.NetCon(sst_stim, sst_syn)
nc_sst.weight[0] = 0.01

### Record time ###
t_vec = h.Vector().record(h._ref_t)

### Helper to run and record ###
def run(condition='none'):
    # Turn on or off inhibition types
    nc_pv.weight[0] = 0.01 if condition == 'pv' else 0
    nc_sst.weight[0] = 0.01 if condition == 'sst' else 0

    h.finitialize(-65)
    h.run()
    return list(t_vec), list(cell.v_soma), list(cell.v_dend)

### Run simulations ###
t_base, v_soma_base, _ = run('none')
t_pv, v_soma_pv, _ = run('pv')
t_sst, v_soma_sst, _ = run('sst')

### Plot results ###
plt.figure(figsize=(10, 6))
plt.plot(t_base, v_soma_base, label='No Inhibition')
plt.plot(t_pv, v_soma_pv, label='PV Inhibition (soma)')
plt.plot(t_sst, v_soma_sst, label='SST Inhibition (distal dend)')
plt.xlabel('Time (ms)')
plt.ylabel('Somatic Membrane Potential (mV)')
plt.title('Effect of SST vs PV Inhibition on Pyramidal Soma')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()