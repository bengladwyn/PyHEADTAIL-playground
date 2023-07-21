from __future__ import division

import numpy as np
from scipy.constants import c, e, m_p

from PyCERNmachines.machines import Synchrotron
import PyCERNmachines.SPS.SPSOctupoles as SPSOctupoles


class RR(Synchrotron):

    def __init__(self, *args, **kwargs):

        if 'n_segments' not in kwargs.keys():
            raise ValueError('Number of segments must be specified')

        if 'machine_configuration' not in kwargs.keys():
            raise ValueError('machine_configuration must be specified')

        self.n_segments = kwargs['n_segments']
        self.machine_configuration = kwargs['machine_configuration']

        self.circumference  = 3319.418 
        self.s = (np.arange(0, self.n_segments + 1)
                  * self.circumference / self.n_segments)

        if self.machine_configuration == '53MHz':
            self.charge = e
            self.mass = m_p

            self.gamma = 8e9*e/(self.mass*c**2) + 1

            self.Q_x     = 25.4601
            self.Q_y     = 24.412

            self.Qp_x    = [0]
            self.Qp_y    = [0]

            self.app_x   = 0.0000e-9
            self.app_y   = 0.0000e-9
            self.app_xy  = 0

            self.alpha_x = 0 * np.ones(self.n_segments + 1)
            self.beta_x  = self.circumference/(2*np.pi*self.Q_x) * np.ones(self.n_segments + 1)
            self.D_x     = 0 * np.ones(self.n_segments + 1)
            self.alpha_y = 0 * np.ones(self.n_segments + 1)
            self.beta_y  = self.circumference/(2*np.pi*self.Q_y) * np.ones(self.n_segments + 1)
            self.D_y     = 0 * np.ones(self.n_segments + 1)

            self.alpha       = 0.00243665009164
            self.h1          = 588
            self.h2          = 588
            self.V1          = 80e3
            self.V2          = 0
            self.dphi1       = np.pi
            self.dphi2       = np.pi
            self.p_increment = 0 * e/c * self.circumference/(self.beta*c)

            self.longitudinal_focusing = 'non-linear'

        elif self.machine_configuration == '2.5MHz':
            self.charge = e
            self.mass = m_p

            self.gamma = 8e9*e/(self.mass*c**2) + 1

            self.Q_x     = 25.4601
            self.Q_y     = 24.412

            self.Qp_x    = [0]
            self.Qp_y    = [0]

            self.app_x   = 0.0000e-9
            self.app_y   = 0.0000e-9
            self.app_xy  = 0

            self.alpha_x = 0 * np.ones(self.n_segments + 1)
            self.beta_x  = self.circumference/(2*np.pi*self.Q_x) * np.ones(self.n_segments + 1)
            self.D_x     = 0 * np.ones(self.n_segments + 1)
            self.alpha_y = 0 * np.ones(self.n_segments + 1)
            self.beta_y  = self.circumference/(2*np.pi*self.Q_y) * np.ones(self.n_segments + 1)
            self.D_y     = 0 * np.ones(self.n_segments + 1)

            self.alpha       = 0.00243665009164
            self.h1          = 28 
            self.h2          = 28
            self.V1          = 90e3
            self.V2          = 0
            self.dphi1       = np.pi
            self.dphi2       = np.pi
            self.p_increment = 0 * e/c * self.circumference/(self.beta*c)

            self.longitudinal_focusing = 'non-linear'

        else:
            raise ValueError('ERROR: unknown machine configuration ' +
                             self.machine_configuration)

        super(RR, self).__init__(*args, **kwargs)



