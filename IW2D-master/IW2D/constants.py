from cppyy.gbl import amp, precision, C, mu0, eps0, Z0

fMP = amp.ampf[precision]
cMP = amp.campf[precision]
oneMP = fMP(1)
twoMP = fMP(2)
jimagMP = cMP(0)
jimagMP.y = 1
piMP = amp.pi[precision]()
twopiMP = amp.twopi[precision]()
