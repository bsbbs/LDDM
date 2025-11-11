import random
import numpy as np
from neuron import n
from neuron.units import ms, mV, µm

# --- 2. Neuronal Morphology Parameters ---
AxonRange_EtoE = 130 # um, standard deviation of exponential decay of connection probability over somatic distance from E to E
AxonRange_EtoI = 100 # um, standard deviation of exponential decay of connection probability over somatic distance from E to I
AxonRange_ItoE = 97 # um, standard deviation of exponential decay of connection probability over somatic distance from I to E
CnnctProb_EtoE = .5 # the maximum connection probability from E to E
CnnctProb_EtoI = .5 # the maximum connection probability from E to I
CnnctProb_ItoE = .5 # the maximum connection probability from I to E

# --- 3. Neuron Biophysiology Parameters ---
Cm = 1 # micro Farads / cm^2, Membrane capacitance
Ra = 35.4 # Ohm * cm, Axial resistance
# Hodgkin-Huxley kinetics parameters
g_Nabar = 0.12  # S/cm2, Na channel conductance 
g_Kbar = 0.036 # S/cm2, K channel conductance
g_L = 0.0003   # S/cm2, Leak conductance
E_L = -54.3 * mV # Reversal potential for leaky channel

# --- 4. Neuron and Synapse Parameters ---
# Synaptic parameters
Delay_EtoE = 1.7 * ms # synaptic delay from E to E
Delay_ItoE = 1 * ms # synaptic delay from I to E
Delay_EtoI = 1.3 * ms # synaptic delay from E to I
# Campagnola et al., 2022 Recurrent excitatory connections in human cortex (E→E)
# have longer latency than those with a pre- or postsynaptic inhibitory cell (E→E median 1.73 versus I→E 1.04 ms, E→I 1.34 ms)
AMPA_TAU = 2 * ms # the decay time constant of EPSC,  AMPA 
GABAA_TAU = 6 * ms # the decay time constant of IPSC, like GABAa
AMPA_RvrslP = 0 * mV # reversal potential for AMPA receptor
GABA_RvrslP = -80 * mV # reversal potential for GABA receptor
gbar_AMPA = 1.4e-3 # µS -> 1.4 nS, maximum conductance of a single point AMPA receptor
gbar_GABA = 3.5e-3 # µS -> 3.5 nS, maximum conductance of a single point GABA receptor
NMDA_TAU1 = 2 * ms # rising time constant for NMDA, Feldmeyer et al., 2002
NMDA_TAU2 = 26 * ms # decaying time constant for NMDA, Feldmeyer et al., 2002
NMDA_a = 0.062; # mV^-1
NMDA_b = 3.57; # mM
NMDA_Mg2 = 1; # nM

# --- 5. E and I STDP Parameters ---
init_weight = .5
max_weight = 1
# STDP parameters (for E -> I connections)
STDP_TAU_PLUS = 20 * ms
STDP_TAU_MINUS = 20 * ms
STDP_A_PLUS = 0.01  # LTP strength
STDP_A_MINUS = 0.012 # LTD strength
W_MAX = 0.01        # Maximum synaptic weight

class Cell:
    # Allow to define the physical location and rotation of a cell
    def __init__(self, gid, x, y, z, theta):
        self._gid = gid
        self._setup_morphology()
        self._setup_biophysics()
        self.x = self.y = self.z = 0
        n.define_shape()
        self._rotate_z(theta)
        self._set_position(x, y, z)
        self._spike_detector = n.NetCon(self.soma(0.5)._ref_v, None, sec=self.soma)
        self.spike_times = n.Vector()
        self._spike_detector.record(self.spike_times)

        self._ncs = [] # connections

    def __repr__(self):
        return "{}[{}]".format(self.name, self._gid)

    def _set_position(self, x, y, z):
        for sec in self.all:
            for i in range(sec.n3d()):
                sec.pt3dchange(
                    i,
                    x - self.x + sec.x3d(i),
                    y - self.y + sec.y3d(i),
                    z - self.z + sec.z3d(i),
                    sec.diam3d(i),
                )
        self.x, self.y, self.z = x, y, z

    def _rotate_z(self, theta):
        """Rotate the cell about the Z axis."""
        for sec in self.all:
            for i in range(sec.n3d()):
                x = sec.x3d(i)
                y = sec.y3d(i)
                c = n.cos(theta)
                s = n.sin(theta)
                xprime = x * c - y * s
                yprime = x * s + y * c
                sec.pt3dchange(i, xprime, yprime, sec.z3d(i), sec.diam3d(i))


class ENeuron(Cell):
    name = "E"

    def _setup_morphology(self):
        # define a pyramidal neuron as two compartments, a soma body and a dendrite
        self.soma = n.Section("soma", self)
        self.soma.L = self.soma.diam = 20 # 10
        self.dend1 = n.Section("dend1", self)
        self.dend2 = n.Section("dend2", self)
        for dend in [self.dend1, self.dend2]:
            dend.connect(self.soma)
            dend.L = 200
            dend.diam = 1
            dend.nseg = 3 # 3 segments for distal/proximal targeting
        self.all = self.soma.wholetree()
 
    def _setup_biophysics(self):
        for sec in self.all:
            sec.Ra = Ra  # Axial resistance in Ohm * cm
            sec.cm = Cm  # Membrane capacitance in micro Farads / cm^2
        # Insert hh channel to all segments of the cell
        self.soma.insert(n.hh)
        for seg in self.soma:
            seg.hh.gnabar = g_Nabar  # Sodium conductance in S/cm2
            seg.hh.gkbar = g_Kbar  # Potassium conductance in S/cm2
            seg.hh.gl = g_L  # Leak conductance in S/cm2
            seg.hh.el = E_L  # Reversal potential for leaky channel in mV
        for dend in [self.dend1, self.dend2]:
            dend.insert(n.hh)
            for seg in dend:
                seg.hh.gnabar = g_Nabar  # Sodium conductance in S/cm2
                seg.hh.gkbar = g_Kbar  # Potassium conductance in S/cm2
                seg.hh.gl = g_L  # Leak conductance in S/cm2
                seg.hh.el = E_L  # Reversal potential for leaky channel in mV
            
        # Setup the synapses
        # AMPA 
        self.AMPA_soma = n.ExpSyn(self.soma(0.5)) # receiving excitatory input from the perisoma area
        self.AMPA_soma.tau = AMPA_TAU
        self.AMPA_soma.e = AMPA_RvrslP
        self.AMPA_dend = n.ExpSyn(self.dend1(1)) # receiving excitatory input from the distal dendrites
        self.AMPA_dend.tau = AMPA_TAU
        self.AMPA_dend.e = AMPA_RvrslP
        self.AMPA_dend2 = n.ExpSyn(self.dend2(1)) # receiving excitatory input from the distal dendrites
        self.AMPA_dend2.tau = AMPA_TAU
        self.AMPA_dend2.e = AMPA_RvrslP
        # GABAa
        self.GABA_soma = n.ExpSyn(self.soma(0.5)) # receiving inhibitory input from the perisoma area
        self.GABA_soma.tau = GABAA_TAU
        self.GABA_soma.e = GABA_RvrslP
        self.GABA_dend = n.ExpSyn(self.dend1(1)) # receiving inhibitory input from the distal dendrites
        self.GABA_dend.tau = GABAA_TAU
        self.GABA_dend.e = GABA_RvrslP
        self.GABA_perisoma = n.ExpSyn(self.dend1(0)) # receiving inhibitory input from the distal dendrites
        self.GABA_perisoma.tau = GABAA_TAU
        self.GABA_perisoma.e = GABA_RvrslP
        self.GABA_dend2 = n.ExpSyn(self.dend2(1)) # receiving inhibitory input from the distal dendrites
        self.GABA_dend2.tau = GABAA_TAU
        self.GABA_dend2.e = GABA_RvrslP


class INeuron(Cell):
    name = "I"

    def _setup_morphology(self):
        # simplified neuromorphology for I enuron as a point model
        self.soma = n.Section("soma", self)
        self.soma.L = self.soma.diam = 10 # 5
        self.all = self.soma.wholetree()

    def _setup_biophysics(self):
        for sec in self.all:
            sec.Ra = Ra  # Axial resistance in Ohm * cm
            sec.cm = Cm  # Membrane capacitance in micro Farads / cm^2
        # Insert hh channel to point model of the soma
        self.soma.insert(n.hh)
        for seg in self.soma:
            seg.hh.gnabar = g_Nabar  # Sodium conductance in S/cm2
            seg.hh.gkbar = g_Kbar  # Potassium conductance in S/cm2
            seg.hh.gl = g_L  # Leak conductance in S/cm2
            seg.hh.el = E_L  # Reversal potential for leaky channel in mV
     
        # Setup the synapses
        # AMPA only, not consider disinhibition
        self.AMPA = n.ExpSyn(self.soma(0.5))
        self.AMPA.tau = AMPA_TAU
        self.AMPA.e = AMPA_RvrslP # the reversal potential of AMPA receptor


class Belt:
    """A network of E cells and I cells
    with their connections defined as E to E,
    E to I, and I to E. The I to E connection 
    can be defined like SST, targeting the distal dendrite of the E cell
    or like PV, targeting the perisomatic area of the E cell
    """
    def __init__(self, NE=100, NI = 25, SPACE_EXTENT_X = 300, SPACE_EXTENT_Y = 100, Itype = "SST"):
        """
        :param NE: Number of E cells
        :param NI: Number of I cells
        :param SPACE_EXTENT_X: physical spreading over the x axis
        :param SPACE_EXTENT_Y: physical spreading over the y axis
        :param Itype: type of the I cells
        """
        self._create_cells(NE, NI, SPACE_EXTENT_X, SPACE_EXTENT_Y)
        # self._connect_EtoE()
        self._connect_EtoI()
        if Itype == "SST":
            self._connect_ItoE_SST()
        if Itype == "PV":
            self._connect_ItoE_PV()

    def _create_cells(self, NE, NI, SPACE_EXTENT_X, SPACE_EXTENT_Y):
        self.Ecells = []
        for i in range(NE):
            self.Ecells.append(
                ENeuron(i, (NE-i)/NE * SPACE_EXTENT_X, 0.5 * SPACE_EXTENT_Y, 0, n.PI/2)
            )
        self.Icells = []
        for i in range(NI):
            self.Icells.append(
                INeuron(i, (NI-i)/NI * SPACE_EXTENT_X,  -0.5 * SPACE_EXTENT_Y, 0, n.PI/2)
            )
            
    # From E to E
    def _connect_EtoE(self):
        for source in self.Ecells:
            for target in self.Ecells:
                dist = ((source.x - target.x)**2 + (source.y - target.y)**2 + (source.z - target.z)**2)**0.5
                conn_prob = CnnctProb_EtoE*np.exp(-.5*(dist/AxonRange_EtoE)**2)
                if random.random() < conn_prob:
                    nc = n.NetCon(source.soma(0)._ref_v, target.AMPA_soma, sec=source.soma)
                    nc.weight[0] = gbar_AMPA * init_weight * random.random()
                    nc.delay = Delay_EtoE
                    source._ncs.append(nc)

    # From E to I
    def _connect_EtoI(self):
        for source in self.Ecells:
            for target in self.Icells:
                dist = ((source.x - target.x)**2 + (source.y - target.y)**2 + (source.z - target.z)**2)**0.5
                conn_prob = CnnctProb_EtoI#*np.exp(-.5*(dist/AxonRange_EtoI)**2)
                if random.random() < conn_prob:
                    nc = n.NetCon(source.soma(0)._ref_v, target.AMPA, sec=source.soma)
                    nc.weight[0] = gbar_AMPA * init_weight * 0 * random.uniform(0.3,1)
                    nc.delay = Delay_EtoI
                    target._ncs.append(nc)
                
    # From I to E, like SST
    def _connect_ItoE_SST(self):
        for source in self.Icells:
            for target in self.Ecells:
                dist = ((source.x - target.x)**2 + (source.y - target.y)**2 + (source.z - target.z)**2)**0.5
                conn_prob = CnnctProb_ItoE#*np.exp(-.5*(dist/AxonRange_ItoE)**2)
                if random.random() < conn_prob:
                    nc = n.NetCon(source.soma(0)._ref_v, target.GABA_dend, sec=source.soma)
                    nc.weight[0] = gbar_GABA * init_weight * random.uniform(0.3,1) #random.random()
                    nc.delay = Delay_ItoE
                    target._ncs.append(nc)
                
    # From I to E, like PV
    def _connect_ItoE_PV(self):
        for source in self.Icells:
            for target in self.Ecells:
                dist = ((source.x - target.x)**2 + (source.y - target.y)**2 + (source.z - target.z)**2)**0.5
                conn_prob = CnnctProb_ItoE#*np.exp(-.5*(dist/AxonRange_ItoE)**2)
                if random.random() < conn_prob:
                    nc = n.NetCon(source.soma(0)._ref_v, target.GABA_perisoma, sec=source.soma)
                    nc.weight[0] = gbar_GABA * init_weight * random.random()
                    nc.delay = Delay_ItoE
                    source._ncs.append(nc)