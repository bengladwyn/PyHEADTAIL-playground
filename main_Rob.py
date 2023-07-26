from __future__ import division, print_function

import numpy as np
np.random.seed(42)

from scipy.constants import c, epsilon_0, e, m_p

import matplotlib as mpl
import matplotlib.pyplot as plt

import h5py
import sys, os, shutil
sys.path.append('/Users/rainswor/Installs/PyHEADTAIL/')
from PyHEADTAIL.particles.slicing import UniformBinSlicer
from PyHEADTAIL.spacecharge.spacecharge import TransverseGaussianSpaceCharge
#from PyHEADTAIL.spacecharge.pypic_factory import create_mesh, create_3dmesh_from_beam
#from PyHEADTAIL.spacecharge.pypic_spacecharge import SpaceChargePIC
from PyHEADTAIL.monitors.monitors import ParticleMonitor

from PyHEADTAIL.general.printers import SilentPrinter
from PyHEADTAIL.general.contextmanager import CPU, GPU
contextmanager = CPU
from FERMImachines import RR

from PyHEADTAIL.impedances.wakes import WakeField, WakeTable, Resonator, CircularResonator, ParallelPlatesResonator
from PyHEADTAIL.impedances.wakes import ResistiveWall, CircularResistiveWall, ParallelPlatesResistiveWall

n_macroparticles = int(5e4)
n_slices_sc = 32

intensity = 5e11
epsn_x = epsn_y = 2.5e-6*np.pi # in [m.rad]
sigma_z = 0.57 # in [m]

def make_machine(n_segments=1):
    return RR(n_segments=n_segments,
               machine_configuration='53MHz',
               optics='smooth',
               printer=SilentPrinter(),
              )

def make_beam(machine=make_machine(), n_macroparticles=n_macroparticles):
    return machine.generate_6D_Gaussian_bunch_matched(
        n_macroparticles, intensity, epsn_x, epsn_y, sigma_z)

m = make_machine()
beam = make_beam(m)

sig_x = beam.sigma_x()
sig_y = beam.sigma_y()

slicing_interval = m.longitudinal_map.get_bucket(beam).interval
slicer_sc = UniformBinSlicer(n_slices_sc, z_cuts=slicing_interval)

lmbda = intensity * e / (np.sqrt(2*np.pi) * sigma_z)
Ksc = e / (beam.gamma**3 * m_p * (beam.beta * c)**2) * lmbda / (2*np.pi*epsilon_0)
R = m.circumference / (2*np.pi)

def dQ_inc(thissize, theothersize, thistune, Ksc=Ksc):
    'incoherent KV tune shift'
    return Ksc * R**2 / (4 * thistune * thissize * (thissize+theothersize))

print ('dQ_x = {0:.3f} and dQ_y = {1:.3f}'.format(
    dQ_inc(beam.sigma_x(), beam.sigma_y(), m.Q_x),
    dQ_inc(beam.sigma_y(), beam.sigma_x(), m.Q_y)))

##add space charge
#assert (m.optics == 'smooth')
#sc_integration_length = m.circumference / len(m.transverse_map)
#sc_node = TransverseGaussianSpaceCharge(slicer=slicer_sc, length=sc_integration_length)
#m.install_after_each_transverse_segment(sc_node)

##add wake
n_sigma_z = 2
n_slices = 50
uniform_bin_slicer = UniformBinSlicer(n_slices=n_slices, n_sigma_z=n_sigma_z)
#resis_para = ParallelPlatesResistiveWall(pipe_radius=2.2e-2, resistive_wall_length=m.circumference,
                                    #conductivity=1.35e6, dt_min=1e-3/c)
#wake_field = WakeField(uniform_bin_slicer, resis_para)
wakefile = 'RR_fullWake_resistiveWall_PyHTConvention.dat'
#wakefile_k = 'RR_fullWake_kicker2_PyHTConvention.dat'
wake_table = WakeTable(wakefile, ['time', 'dipole_x', 'dipole_y', 'quadrupole_x', 'quadrupole_y','nonsense'],n_turns_wake=5)   # Follow the order of the columns in the wake file here. If you want to exclude a component, change the name e.g. for the dipole y to 'no_dipole_y'.
#wake_table_k = WakeTable(wakefile_k, ['time', 'dipole_x', 'dipole_y', 'quadrupole_x', 'quadrupole_y','nonsense'])   # Follow the order of the columns in the wake file here. If you want to exclude a component, change the name e.g. for the dipole y to 'no_dipole_y'.
wake_field = WakeField(uniform_bin_slicer, wake_table) #, wake_table_k, wake_table_k)
#wake_field_k = WakeField(uniform_bin_slicer, wake_table_k)
m.one_turn_map.append(wake_field)
#m.one_turn_map.append(wake_field_k)

##Do the tracking
#m.one_turn_map.append(wake_field)
m_x=[]
m_y=[]
turns=2048
beam.y+=0e-3
print(m)
for i in range(turns):
    m.track(beam)
    m_x.append(beam.mean_x())
    m_y.append(beam.mean_y())
    if i%100==0:
        print (i)
m_x = np.array(m_x)
m_y = np.array(m_y)
dict={'m_x':m_x,'m_y':m_x,'intensity':intensity}
import pickle
from decimal import Decimal
name = '%.2E' % Decimal(dict['intensity'])
fname = 'int_'+name+'_IW2D.dat'
f=open(fname,'wb')
pickle.dump(dict,f)
f.close()
ytune= (np.argmax(np.abs(np.fft.fft(m_y)))/float(turns))
xtune= (np.argmax(np.abs(np.fft.fft(m_x)))/float(turns))

if ytune > 0.5:
    print (1-ytune)
else:
    print (ytune)
if xtune > 0.5:
    print (1-xtune)
else:
    print (xtune)

# 0.5e11
# ================================
# 0.41162109375
# 0.4599609375

# 5e11
# ================================
# 0.46044921875