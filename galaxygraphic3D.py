import matplotlib.pyplot as plt
import numpy as np
import commonFunctions as cf

def generate3Dgraphic(i, minvarphiX, minvarphi, minrho, galaxdata):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ssurf = np.linspace(0.0015, 0.001, 50)
    minorssurf = np.linspace(0.0015, 0.001, 8)
    rhosurf = np.linspace(0.01, 0.03, 50)
    # csurf = cf.chiquad(rhosurf, ssurf, galaxdata) # Scatter data
    # surf = ax.scatter3D(ssurf, rhosurf, csurf)    # Scatter data
    Xsurf, Ysurf = np.meshgrid(ssurf, rhosurf)
    chisurf = np.zeros((len(Xsurf), len(Ysurf)))
    for k in range(len(Xsurf)):
        for j in range(len(Ysurf)):
            chielem = cf.chiquad(np.array([Ysurf[k,j]]), np.array([Xsurf[k,j]]), galaxdata)
            chisurf[k][j] = chielem
    surf = ax.plot_surface(Xsurf, Ysurf, chisurf, cmap="winter", alpha=0.4)     # Superficie
    # Función rho_0(s) para s in (0, infty) // varphi(s) es la proyección de esta rho_0(s) en el plano rho_0 = 0 #
    ax.plot(ssurf, cf.rho(ssurf, galaxdata), cf.chiquad(cf.rho(ssurf, galaxdata), ssurf, galaxdata), c='r')
    # Función chi^2(·, s)  para un s dado #
    ax.plot([minvarphiX] * 50, rhosurf, cf.chiquad(rhosurf, np.asarray([minvarphiX] * 50), galaxdata), c='green')
    ax.scatter3D(minvarphiX, minrho, minvarphi, c='black', linewidths=2)
    for x in minorssurf:
        ax.plot([x] * 50, rhosurf, cf.chiquad(rhosurf, np.asarray([x] * 50), galaxdata), c='green', alpha=0.4)
        convexmin = min(cf.chiquad(rhosurf, np.asarray([x] * 50), galaxdata))
        ax.scatter3D(x, cf.rho(np.array([x]), galaxdata), convexmin, c='red', linewidths=0.1)       # Funciones convexas
    ax.view_init(0, 90)
    ax.yaxis.labelpad = 5
    ax.xaxis.labelpad = 10
    ax.set_xlabel('s')
    ax.set_ylabel(r'$\rho_0$')
    ax.set_zlabel(r'$\chi^2(\rho_0,s)$')
    plt.title("Galaxia " + i + " con perfil " + galaxdata["profile"])
    plt.show()
    #plt.savefig("galaxies/graphics/" + galaxdata["profile"] + '-' + i + '-3D_convexity0190.png')
