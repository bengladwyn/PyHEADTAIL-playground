import numpy as np
np.random.seed(42)

import matplotlib.pyplot as plt

import seaborn as sns
sns.set_context('talk')

import copy

from scipy.constants import c, e, m_p

import PyHEADTAIL

from FERMImachines import RR

from PyHEADTAIL.particles.slicing import UniformBinSlicer
from PyHEADTAIL.impedances.wakes import CircularResistiveWall, CircularResonator, WakeField, WakeTable


tune_freqs = []
n_turns = 2048

    
machine = RR(n_segments=1, machine_configuration='53MHz', 
            optics='smooth', printer=PyHEADTAIL.general.printers.SilentPrinter())

C = machine.circumference

epsn_x = epsn_y = 2.5e-6*np.pi # in [m.rad]

bunch = machine.generate_6D_Gaussian_bunch_matched(
    n_macroparticles=int(5e4), intensity=5e11, 
    epsn_x=epsn_x, epsn_y=epsn_y, sigma_z=0.57)

n_sigma_z = 2
n_slices = 50
uniform_bin_slicer = UniformBinSlicer(n_slices=n_slices, n_sigma_z=n_sigma_z)

#wake_table = CircularResonator(R_shunt=1e6, frequency=1e9, Q=1)

wakefile = 'RR_fullWake_resistiveWall_PyHTConvention.dat'
wake_table = WakeTable(wakefile, ['time', 'dipole_x', 'dipole_y', 'quadrupole_x', 'quadrupole_y','nonsense'],n_turns_wake=5)   # Follow the order of the columns in the wake file here. If you want to exclude a component, change the name e.g. for the dipole y to 'no_dipole_y'.
wake_field = WakeField(uniform_bin_slicer, wake_table) #, wake_table_k, wake_table_k)

#wake_table = CircularResistiveWall(pipe_radius=5e-2, resistive_wall_length=C, conductivity=3e9, dt_min=1e-3/c, beta=bunch.beta, n_turns_wake=100)
#wake_field = WakeField(uniform_bin_slicer, wake_table)
machine.one_turn_map.append(wake_field)

m_x = []
m_y = []

for i in range(n_turns):
    machine.track(bunch)
    m_x.append(bunch.mean_x())
    m_y.append(bunch.mean_y())
    if i % 100 == 0:
        print(i)

m_x = np.array(m_x)
m_y = np.array(m_y)

max_index = np.argmax(np.abs(np.fft.rfft(m_x)))
freqs = np.fft.rfftfreq(n_turns)
max_freq = freqs[max_index]

print(max_freq)

# 0.5e11: 0.4599609375
# 5e11: 0.46044921875