from neuron import h, gui
import matplotlib.pyplot as plt

# Create soma and two dendrites
soma = h.Section(name='soma')
dend1 = h.Section(name='dend1')
dend2 = h.Section(name='dend2')

dend1.connect(soma(1))
dend2.connect(soma(0))

soma.L = soma.diam = 20  # spherical soma
for dend in [dend1, dend2]:
    dend.L = 200
    dend.diam = 1

# Insert active properties
for sec in [soma, dend1, dend2]:
    sec.insert('hh')

# Add excitatory synapse on dendrite 1
exc_syn = h.ExpSyn(dend1(0.9))
exc_syn.e = 0  # reversal for excitatory

# Add SST inhibitory synapse on dendrite 1
sst_syn = h.ExpSyn(dend1(0.9))
sst_syn.e = -75  # reversal for GABA_A

# Add PV inhibitory synapse near soma
pv_syn = h.ExpSyn(soma(0.1))
pv_syn.e = -75  # reversal for GABA_A

# Create stimulators
exc_stim = h.NetStim()
exc_stim.number = 1
exc_stim.start = 5  # ms

sst_stim = h.NetStim()
sst_stim.number = 1
sst_stim.start = 5  # ms

pv_stim = h.NetStim()
pv_stim.number = 1
pv_stim.start = 5  # ms

# Connect NetCons
exc_nc = h.NetCon(exc_stim, exc_syn)
exc_nc.weight[0] = 0.05  # excitatory weight

sst_nc = h.NetCon(sst_stim, sst_syn)
sst_nc.weight[0] = 0.1  # inhibitory weight

pv_nc = h.NetCon(pv_stim, pv_syn)
pv_nc.weight[0] = 0.1  # inhibitory weight

# Record membrane potentials
t_vec = h.Vector().record(h._ref_t)
v_soma = h.Vector().record(soma(0.5)._ref_v)

# Run different conditions
results = {}

for condition in ['exc_only', 'exc_sst', 'exc_pv']:
    # Reset
    h.finitialize(-65)
    exc_nc.weight[0] = 0.05
    sst_nc.weight[0] = 0.0
    pv_nc.weight[0] = 0.0

    if condition == 'exc_sst':
        sst_nc.weight[0] = 0.1  # turn on SST inhibition
    elif condition == 'exc_pv':
        pv_nc.weight[0] = 0.1  # turn on PV inhibition

    h.continuerun(40)

    # Save trace
    results[condition] = list(v_soma)

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(t_vec, results['exc_only'], label='Excitation only')
plt.plot(t_vec, results['exc_sst'], label='Excitation + SST inhibition')
plt.plot(t_vec, results['exc_pv'], label='Excitation + PV inhibition')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane potential (mV)')
plt.title('Differential Effects of SST vs PV Inhibition on AP Firing')
plt.legend()
plt.show()