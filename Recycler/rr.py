from cpymad.madx import Madx
madx = Madx()
madx.call('recycler.madx')
madx.survey(sequence="RING605_FODO",x0 = 29607.699942, y0 = 219.575544,z0 = 30941.487902,theta0=2.596497)
