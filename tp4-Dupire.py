import numpy as np
import math as  mt
from numpy import linalg as LA
import matplotlib.pyplot as plt 
import math as  mt
from pylab import *

Tmax = 0.5
Kmax=20
B1 = 1
B2 =  1
r = 0.1
S0 = 10
N = 199
M = 49
K = np.linspace(0,Kmax,N+2)
dK = Kmax/(N+1)
T = np.linspace(0,Tmax,M+2)
dT = Tmax/(M+1)

def prix_dupire(sigma):
   
    V = [[0 for j in range(N+2)] for i in range(M+2)]
    for i in range(N+2):
        V[0][i] = max(S0 - K[i], 0)

    for n in range(M+2):
        V[n][0] = S0
        V[n][N+1] = 0
    
    #Crank-Nicolson
    A = [0 for i in range(N+1)]
    B = [0 for i in range(N+1)]
    D = [0 for i in range(N+1)]
    C = [[0 for j in range(N+1)] for i in range(M+1)]

    for i in range(1,N+1): 
        A[i] = (dT/4)*(K[i]*r/dK - ((sigma[i]*K[i]/dK)**2))
        B[i] = -(dT/4)*(K[i]*r/dK + ((sigma[i]*K[i])/dK)**2)
        D[i] = 1 + (dT/2)*((sigma[i]*K[i]/dK)**2)
    
    Cetoile = [[0 for j in range(N+1)] for i in range(M+1)]
    Detoile = np.zeros(N+1)
    
    for n in range(M+1):
        for i in range(1,N+1):
            if i==1:
                C[n][i]=-(dT/4)*(K[i]*r/dK - ((sigma[i]**2) * (K[i]**2))/(dK**2))*V[n][i+1] + (1-dT/2*(sigma[i]**2)*(K[i]**2)/(dK**2))*V[n][i] + dT/4*(r*K[i]/dK + (sigma[i]**2) * (K[i]**2) /(dK**2))*V[n][i-1] - S0*B[2]
            else:
                C[n][i]=-(dT/4)*(K[i]*r/dK - ((sigma[i]**2) * (K[i]**2))/(dK**2))*V[n][i+1] + (1-dT/2*(sigma[i]**2)*(K[i]**2)/(dK**2))*V[n][i] + dT/4*(r*K[i]/dK + (sigma[i]**2) * (K[i]**2) /(dK**2))*V[n][i-1] 
        
        Cetoile[n][1] =  C[n][1] 
        Detoile[1] = D[1] 
        
        for i in range(2,N+1):
            Detoile[i] = D[i] - (B[i] * A[i-1] / Detoile[i-1])
            Cetoile[n][i] = C[n][i]- (B[i] * Cetoile[n][i-1] / Detoile[i-1])
        
        V[n+1][N] = Cetoile[n][N] / Detoile[N]
        
        for i in range(N-1,0,-1):
            V[n+1][i] = (Cetoile[n][i] -A[i] * V[n+1][i+1])/ Detoile[i]
    return V

def grapheVega(V, title): 
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(title, fontsize=16) 
    
    # Graphique t=0
    axs[0, 0].plot(K, V[0], label='t=0')
    axs[0, 0].set_title('Vega Dupire t=0')
    axs[0, 0].set_xlabel('Strike K')
    axs[0, 0].set_ylabel('Volatilité')
    axs[0, 0].legend()

    # Graphique t=T/2
    axs[0, 1].plot(K[3:], V[int((M)/2)][3:], label='t=T/2')
    axs[0, 1].set_title('Vega Dupire t=T/2')
    axs[0, 1].set_xlabel('Strike K')
    axs[0, 1].set_ylabel('Volatilité')
    axs[0, 1].legend()

    # Graphique t=T
    axs[1, 0].plot(K[3:], V[M][3:], label='t=T')
    axs[1, 0].set_title('Vega Dupire t=T')
    axs[1, 0].set_xlabel('Strike K')
    axs[1, 0].set_ylabel('Volatilité')
    axs[1, 0].legend()

    # Surface 3D dans le dernier subplot
    X1, Y1 = np.meshgrid(K, T)
    Z1 = np.array(V)
    ax3d = fig.add_subplot(2, 2, 4, projection="3d")
    ax3d.plot_surface(X1, Y1, Z1, cmap=plt.cm.coolwarm, alpha=0.8)
    ax3d.set_title('Surface Vega Dupire')
    ax3d.set_xlabel('Strike K')
    ax3d.set_ylabel('Maturité T')
    ax3d.set_zlabel('Volatilité')

    # Ajustement des espacements
    plt.tight_layout(rect=[0, 0, 1, 0.95])  # Laisser de la place pour le titre global
    plt.show()


def grapheDupire(V, title): 
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(title, fontsize=16)  # Titre global pour distinguer les cas
    
    # Graphique t=0
    axs[0, 0].plot(K, V[0], label='t=0')
    axs[0, 0].set_title('Pricing Call Dupire t=0')
    axs[0, 0].set_xlabel('Strike K')
    axs[0, 0].set_ylabel('Prix V(t,K)')
    axs[0, 0].legend()

    # Graphique t=T/2
    axs[0, 1].plot(K[3:], V[int((M)/2)][3:], label='t=T/2')
    axs[0, 1].set_title('Pricing Call Dupire t=T/2')
    axs[0, 1].set_xlabel('Strike K')
    axs[0, 1].set_ylabel('Prix V(t,K)')
    axs[0, 1].legend()

    # Graphique t=T
    axs[1, 0].plot(K[3:], V[M][3:], label='t=T')
    axs[1, 0].set_title('Pricing Call Dupire t=T')
    axs[1, 0].set_xlabel('Strike K')
    axs[1, 0].set_ylabel('Prix V(t,K)')
    axs[1, 0].legend()

    # Surface 3D dans le dernier subplot
    X1, Y1 = np.meshgrid(K, T)
    Z1 = np.array(V)
    ax3d = fig.add_subplot(2, 2, 4, projection="3d")
    ax3d.plot_surface(X1, Y1, Z1, cmap=plt.cm.coolwarm, alpha=0.8)
    ax3d.set_title('Surface Pricing Call Dupire')
    ax3d.set_xlabel('Strike K')
    ax3d.set_ylabel('Maturité T')
    ax3d.set_zlabel('Prix V(K,T)')

    # Ajustement des espacements
    plt.tight_layout(rect=[0, 0, 1, 0.95])  # Laisser de la place pour le titre global
    plt.show()


def vega_dupire(sigma, h, title=None):
    V = prix_dupire(sigma)
    sigmah = [sigma[i] + h for i in range(N+2)]
    V_dilate = prix_dupire(sigmah)
    Vega = [np.zeros(N+2) for i in range(M+2)]
    for n in range(M+2):
        for i in range(N+2):
            Vega[n][i] = (V_dilate[n][i] - V[n][i]) / h
    
    if title is None:
        title = "Vega Dupire"
    grapheVega(Vega, title)


h = 0.01

# Cas 1 : Pricing Dupire avec sigma constant
sig = [0.3 for i in range(N+2)]
#V = prix_dupire(sig)
#grapheDupire(V, title="Dupire pour σ = 0.3 (constante)")

# Cas 2 : Pricing Dupire avec une loi inverse
sig = [1/(K[i]) for i in range(N+2)]
#V = prix_dupire(sig)
#grapheDupire(V, title="Dupire pour σ = 1/K")

# Cas 1 : Vega Dupire avec sigma constant
sig = [0.3 for i in range(N+2)]
#vega_dupire(sig, h, title="Vega Dupire pour σ = 0.3 (constante)")

# Cas 2 : Vega Dupire avec une loi inverse
sig = [1/(K[i]) for i in range(N+2)]
#vega_dupire(sig, h, title="Vega Dupire pour σ = 1/K")


def prix_dupire_utiles(Kp):
    ip=int((Kp/dK))
    return ip

def Partie2(Kp, Vp):
    a = 5
    m = 5
    beta = [1, 1]
    d = [1, 1]
    eps = 10**(-5)
    b = 0.05
    ro = 0.1
    lambdas = 0.0001
    k = 0  
    Vdupire = np.zeros(len(Vp))
    Vdupire_dilate = np.zeros(len(Vp))
    vega = np.zeros(len(Vp))
    residus = np.zeros(len(Vp))
    J = np.array([np.zeros(2) for i in range(len(Vp))])

    while LA.norm(d) > eps: 
        for p in range(len(Vp)):
            sigma = [(beta[0] / (K[i]**beta[1])) for i in range(N+2)]
            sigmah = [h + sigma[i] for i in range(N+2)]
            V = prix_dupire(sigma)
            V_dilate = prix_dupire(sigmah)

            j = prix_dupire_utiles(Kp[p])
            Vdupire[p] = V[M+1][j]
            Vdupire_dilate[p] = V_dilate[M+1][j]

            vega[p] = -(Vdupire[p] - Vdupire_dilate[p]) / h

            residus[p] = Vp[p] - Vdupire[p]

            J[p][0] = -vega[p] / (Kp[p]**beta[1])
            J[p][1] = vega[p] * mt.log(Kp[p]) * beta[0] / (Kp[p]**beta[1])
        
        # Mise à jour des paramètres
        d = np.dot(-LA.inv(np.dot(J.transpose(), J) + lambdas * np.identity(2)), np.dot(J.transpose(), residus))
        beta = beta + d.transpose()    
        k += 1  # Incrémentation du compteur
        
        # Affichage des informations de progression
        #print(f"{k:<10}{LA.norm(d):<20.6e}{beta[0]:<15.6f}{beta[1]:<15.6f}")
    
    print('Beta1 =', beta[0])
    print('Beta2 =', beta[1])
    print('Nombre d iterations :', k)

    # Graphiques
    plt.plot(Kp, Vp, 'o', label='V_marche')
    plt.plot(Kp, Vdupire, label='V_Dupire')
    plt.xlabel('K_marche')
    plt.ylabel('V')
    plt.title('Prix de marché - prix de dupire à T=T_max post calibration')
    plt.legend()
    plt.show()

    X1, Y1 = np.meshgrid(K, T)
    Z1 = np.array(V)
    plt.rcParams["figure.figsize"] = [16, 9]
    plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot_surface(X1, Y1, Z1, cmap=plt.cm.coolwarm, alpha=0.8)
    plt.title('Pricing Call EU post calibration')
    ax.set_xlabel('K')
    ax.set_ylabel('T')
    ax.set_zlabel('V(K,T)')
    plt.show()

    ZZ = np.subtract(np.array(V_dilate), np.array(V)) / h
    Z1 = np.array(ZZ)
    plt.rcParams["figure.figsize"] = [16, 9]
    plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot_surface(X1, Y1, Z1, cmap=plt.cm.coolwarm, alpha=0.8)
    plt.title('Pricing Vega Call 3D post calibration ')
    ax.set_xlabel('K')
    ax.set_ylabel('T')
    ax.set_zlabel('Vega')
    plt.show()

Kp = [7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14]
Vp=[3.3634,2.9092,2.4703,2.0536,1.6666,1.3167,1.0100,0.7504,0.5389,0.3733,0.2491,0.1599,0.0986,0.0584,0.0332]
Partie2(Kp,Vp)

def Partie3(Kp,Vp):
    a=5
    m=5
    beta=[1,1]
    d=[1,1]
    eps=10**(-5)
    b=0.05
    ro=0.1
    lambdas=0.001
    k=0
    Vdupire = np.zeros(len(Vp))
    Vdupire_dilate = np.zeros(len(Vp))
    vega = np.zeros(len(Vp))
    residus = np.zeros(len(Vp))
    J = np.array([np.zeros(2) for i in range(len(Vp))])
    

    sigma = [b*(ro*(K[i]-m) + np.sqrt((K[i] - m) ** 2 + a ** 2 )) for i in range(N+2)]
    sigmah = [h + sigma[i]  for i in range(N+2)]
    V=prix_dupire (sigma)
    V_dilate= prix_dupire(sigmah)
    X1, Y1 = np.meshgrid(K, T)
    Z1=np.array(V)
    plt.rcParams["figure.figsize"]=[16,9]
    plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot_surface(X1, Y1, Z1, cmap=plt.cm.coolwarm, alpha=0.8)
    plt.title('Surface de Dupire pour a = 5 et m = 5')
    ax.set_xlabel('K')
    ax.set_ylabel('T')
    ax.set_zlabel('V(K,T)')
    plt.show()

    
    ZZ= np.subtract(np.array(V_dilate), np.array(V))/h
    Z1=np.array(ZZ)
    plt.rcParams["figure.figsize"]=[16,9]
    plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot_surface(X1, Y1, Z1, cmap=plt.cm.coolwarm, alpha=0.8)
    plt.title('Surface de Vega pour a = 5 et m = 5 ')
    ax.set_xlabel('K')
    ax.set_ylabel('T')
    ax.set_zlabel('V(K,T)')
    plt.show()

    while LA.norm(d)>eps: 
        for p in range(len(Vp)):
                
            sigma = [b*(ro*(K[i]-beta[1]) + np.sqrt((K[i] - beta[1] ) ** 2 + beta[0] ** 2 )) for i in range(N+2)]
           
            sigmah = [h + sigma[i]  for i in range(N+2)]
            V=prix_dupire (sigma)
            V_dilate= prix_dupire(sigmah)

            j=prix_dupire_utiles(Kp[p])
            Vdupire[p]=V[M+1][j]
            Vdupire_dilate[p]=V_dilate[M+1][j]
            
            vega[p]=-(Vdupire[p]-Vdupire_dilate[p])/h
            
            residus[p]=Vp[p]-Vdupire[p]

            J[p][0] = -vega[p] * b * (beta[0] / mt.sqrt(((Kp[p] - beta[1])**2) + beta[0]**2))
            J[p][1] = vega[p] * b *( -ro + ((Kp[p] - beta[1]) / mt.sqrt( (Kp[p] - beta[1])**2 + (beta[0]**2))))
    
        d = np.dot(-LA.inv(np.dot(J.transpose(), J) + lambdas*np.identity(2)), np.dot(J.transpose(), residus))
        beta = beta + d.transpose()    
        k+=1

    print('beta1 =',beta[0])
    print('beta2 =',beta[1])
    print(k)
 
     
    sigma = [b*(ro*(K[i]-beta[1]) + np.sqrt((K[i] - beta[1]) ** 2 + beta[0] ** 2 )) for i in range(len(K))]
    plt.plot(K,sigma)
    plt.xlabel('K')
    plt.ylabel('sigma')
    plt.title('Volatilité locale calibrée')
    plt.show()

    plt.plot(Kp,Vp, 'o', label = 'V_marche')
    plt.plot(Kp,Vdupire, label = 'V_Dupire')
    plt.xlabel('K_marche')
    plt.ylabel('V')
    plt.title('prix de marché - le prix de dupire à T=T_max post calibration')
    plt.legend()
    plt.show()

    X1, Y1 = np.meshgrid(K, T)
    Z1=np.array(V)
    plt.rcParams["figure.figsize"]=[16,9]
    plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot_surface(X1, Y1, Z1, cmap=plt.cm.coolwarm, alpha=0.8)
    plt.title('Pricing Dupire 3D calibrée')
    ax.set_xlabel('K')
    ax.set_ylabel('T')
    ax.set_zlabel('V(K,T)')
    plt.show()

    ZZ= np.subtract(np.array(V_dilate), np.array(V))/h
    Z1=np.array(ZZ)
    plt.rcParams["figure.figsize"]=[16,9]
    plt.figure()
    ax = plt.axes(projection="3d")
    ax.plot_surface(X1, Y1, Z1, cmap=plt.cm.coolwarm, alpha=0.8)
    plt.title('Pricing Vega 3D calibrée')
    ax.set_xlabel('K')
    ax.set_ylabel('T')
    ax.set_zlabel('V(K,T)')
    plt.show()


Kp = [5,6,7,8,9,10,11,12,13,14,15,16,17,18]
Vp=[5.2705, 4.3783, 3.5510, 2.8138, 2.1833,1.6651, 1.2541, 0.9374, 0.6983, 0.5195, 0.3851, 0.2817,0.1987,0.1277]
#Partie3(Kp,Vp)