import numpy


def fluid_prop(h, T):
    # Densidad [kg/m3] Tref=15 C
    rho0 = 1.225
    #Exponential aprox
    #rho = rho0*numpy.exp(-h/10400)
    #Interp ISA Data
    href = [0, 305, 914, 1524, 1829, 2438, 3048, 3353, 3962, 4572, 4877, 5406, 6096, 6401, 7010, 7620, 7925, 9144, 10058, 12192]
    r = [1, 0.971, 0.915, 0.862, 0.836, 0.786, 0.738, 0.716, 0.671, 0.63, 0.61, 0.57, 0.53, 0.515, 0.48, 0.45, 0.43, 0.374, 0.335, 0.246]
    rho = rho0*numpy.interp(h, href, r)

    # V iscosidad dinamica [kg/m s]
    mu = 1.802e-5
    # c Velocidad del sonido
    #c = 340.2
    hrefc = [0, 1524, 3048, 4572, 6096, 7620, 9144, 10686, 12192]
    cref = [340.3, 334.4, 328.4, 322.2, 316, 309.6, 303.1, 295.4, 294.9]
    c = numpy.interp(h, hrefc, cref)
    return [rho, mu, c]
