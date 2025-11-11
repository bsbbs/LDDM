from neuron import h, gui
import matplotlib.pyplot as plt

# Create cell sections
soma = h.Section(name='soma')
soma.L = soma.diam = 20  # spherical soma
soma.insert('hh')  # Hodgkin-Huxley channels

dend = h.Section(name='dend')
dend.L = 200
dend.diam = 2
dend.insert('pas')  # passive properties

# Connect dend to soma
dend.connect(soma(1))

# Insert recording vectors
t_vec = h.Vector().record(h._ref_t)
v_soma = h.Vector().record(soma(0.5)._ref_v)
v_dend = h.Vector().record(dend(0.5)._ref_v)

# Interneuron spike generators
sst_spike = h.NetStim()
sst_spike.number = 1
sst_spike.start = 20  # ms

pv_spike = h.NetStim()
pv_spike.number = 1
pv_spike.start = 20  # ms

# Inhibitory synapses on pyramidal cell
sst_syn = h.ExpSyn(dend(0.8))  # distal dendrite
sst_syn.e = -75  # inhibitory reversal
sst_syn.tau = 5

pv_syn = h.ExpSyn(soma(0.5))  # soma
pv_syn.e = -75
pv_syn.tau = 1

# NetCons from spike generators to synapses
nc_sst = h.NetCon(sst_spike, sst_syn)
nc_sst.weight[0] = 0.005

nc_pv = h.NetCon(pv_spike, pv_syn)
nc_pv.weight[0] = 0.005

# Inject current into soma to test inhibition
stim = h.IClamp(soma(0.5))
stim.delay = 5
stim.dur = 40
stim.amp = 0.15  # Enough to spike if not inhibited

# Run the simulation
h.tstop = 60
h.run()

# Plotting
plt.figure(figsize=(10, 5))
plt.plot(t_vec, v_soma, label='Soma Vm')
plt.plot(t_vec, v_dend, label='Dendrite Vm')
plt.title('Pyramidal neuron with PV and SST inhibition')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane Potential (mV)')
plt.legend()
plt.grid(True)
plt.show()