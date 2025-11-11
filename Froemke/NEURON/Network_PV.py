from neuron import h
import matplotlib.pyplot as plt

h.load_file('stdrun.hoc')
h.tstop = 200  # short trial for comparison
h.dt = 0.025

### Multi-compartment E neuron ###
class MultiCompECell:
    def __init__(self):
        self.soma = h.Section(name='soma')
        self.soma.L = self.soma.diam = 20
        self.soma.insert('hh')

        self.dend = h.Section(name='apical')
        self.dend.L = 300
        self.dend.diam = 2
        self.dend.insert('pas')
        self.dend.connect(self.soma(1))

        self.v_soma = h.Vector().record(self.soma(0.5)._ref_v)
        self.v_dend = h.Vector().record(self.dend(1.0)._ref_v)
        self.spike_times = h.Vector()
        self.nc = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        self.nc.threshold = -20
        self.nc.record(self.spike_times)

### PV-like I neuron ###
class PVNeuron:
    def __init__(self):
        self.soma = h.Section()
        self.soma.L = self.soma.diam = 10
        self.soma.insert('hh')
        self.spike_times = h.Vector()
        self.nc = h.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        self.nc.threshold = -20
        self.nc.record(self.spike_times)

### Build network ###
E = MultiCompECell()
I = PVNeuron()

### Input sources ###
def create_poisson_input(target_syn, rate, start=10, dur=150):
    stim = h.NetStim()
    stim.number = 1e9
    stim.start = start
    stim.noise = 1
    stim.interval = 1000.0 / rate  # ms
    nc = h.NetCon(stim, target_syn)
    nc.weight[0] = 0.002
    return stim, nc

# Proximal input (peri-soma)
prox_syn = h.ExpSyn(E.soma(0.5))
prox_syn.e = 0
prox_syn.tau = 2
_, prox_nc = create_poisson_input(prox_syn, rate=20)

# Distal input (apical dendrite)
dist_syn = h.ExpSyn(E.dend(1.0))
dist_syn.e = 0
dist_syn.tau = 2
_, dist_nc = create_poisson_input(dist_syn, rate=20)

# PV neuron input (Poisson as well)
inh_drive_syn = h.ExpSyn(I.soma(0.5))
inh_drive_syn.e = 0
inh_drive_syn.tau = 2
_, _ = create_poisson_input(inh_drive_syn, rate=20)

# PV inhibitory synapse to E soma
inh_syn = h.ExpSyn(E.soma(0.5))
inh_syn.e = -75
inh_syn.tau = 5
inh_nc = h.NetCon(I.soma(0.5)._ref_v, inh_syn)
inh_nc.weight[0] = 0.01  # strong inhibitory synapse

### Time vector ###
t_vec = h.Vector().record(h._ref_t)

### Run simulations ###
def run_with_inhibition(on=True):
    inh_nc.weight[0] = 0.01 if on else 0.0
    h.finitialize(-65)
    h.run()
    return list(t_vec), list(E.v_soma), list(E.v_dend)

### Run with and without inhibition ###
t1, soma_v_on, dend_v_on = run_with_inhibition(True)
t2, soma_v_off, dend_v_off = run_with_inhibition(False)

### Plot ###
plt.figure(figsize=(12, 6))
plt.subplot(2, 1, 1)
plt.plot(t1, soma_v_on, label='Soma Vm (with inhibition)')
plt.plot(t2, soma_v_off, '--', label='Soma Vm (no inhibition)')
plt.ylabel('mV')
plt.legend()
plt.title('Somatic potential response')

plt.subplot(2, 1, 2)
plt.plot(t1, dend_v_on, label='Dend Vm (with inhibition)')
plt.plot(t2, dend_v_off, '--', label='Dend Vm (no inhibition)')
plt.xlabel('Time (ms)')
plt.ylabel('mV')
plt.legend()
plt.title('Dendritic potential response')

plt.tight_layout()
plt.show()