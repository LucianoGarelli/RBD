def fluid_prop(h, T):
    # Densidad [kg/m3] Tref=15 C
    rho = 1.225
    # V iscosidad dinamica [kg/m s]
    mu = 1.802e-5
    # c Velocidad del sonido
    c = 343.2
    return [rho, mu, c]
