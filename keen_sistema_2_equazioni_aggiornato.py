# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 01:36:45 2020

@author: Gabriele
"""

#%%
# VERSIONE INGEGNERIZZATA PER FARE QUANTE PIU' COSE AUTOMATICAMENTE ED
# EVITARE IL PROBLEMA DI NON TROVARE DEI PUNTI FISSI COME ACCADUTO NELLA RELAZIONE CON LA SADDLE

#%%
import numpy as np
import matplotlib.pyplot as plt
import sympy as sm
#from scipy.integrate import odeint
#from mpl_toolkits.mplot3d import Axes3D
#%%
from numpy import linalg as LA
#%%
from scipy import sqrt # secondo me si poteva anche non mettere e usare invece np.sqrt() però se fai così il centro che dovrebbe avere la parte reale degli autovalori nulla viene "nan" invece ora viene 10^{-17} ciò zero ma per motivi numerici non è proprio 0
#%%
# DEFINISCO LA FORMA FUNZIONALE CHE SERVE PER 4 FUNZIONI DIVERSE NEL PAPER

def functionalform(c1,c2,c3,c4,l):
    funza = (c1/((c2 - (c3 * l))**2)) - c4
    return funza
#%%
def runge_kutta(funzione,Xn,Dt):
    k1 = funzione(Xn) * Dt
    k2 = funzione(Xn + (1/2)*k1) * Dt
    k3 = funzione(Xn + (1/2)*k2) * Dt
    k4 = funzione(Xn + k3) * Dt
    Xnuovo = Xn + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    #Xnuovo = Xnuovo[:,None] # lo faccio diventare un vero vettore (m,1) invece che (m,)
    
    return Xnuovo
#%%
def flusso_traiettoria_singola(x_iniz,y_iniz,dt,punti_di_tempo,FUNZA):
    X = np.array([x_iniz,y_iniz])
    risultati_integrazione = X
    for r in range(1,punti_di_tempo):
        X = runge_kutta(FUNZA,X,dt)
        risultati_integrazione = np.vstack((risultati_integrazione,X))
    x_flusso = risultati_integrazione[:,0]
    y_flusso = risultati_integrazione[:,1]
    return x_flusso, y_flusso

    #np.vstack((A,b)) questo comando aggiunge il vettore b, inteso come vettore riga, in fondo alla matrice A aggiungendone una riga. Quindi se A è una matrice (m,n), np.vstack((A,b)) restituisce una matrice (m+1,n). Esiste anche il comando analogo per le colonne: np.hstack((A,b)) restituisce una matrice (m,n+1) con l'aggiunta di b inteso come vettore colonna
    
    #risultati_integrazione=np.array([risultati_integrazione,X])
    # QUESTO COMANDO NON FUNZIONAVA poichè non hai capito come funzionano gli array in Python!
    # Di fatto creavi una matrice con 2 sole righe invece di aggiungere il nuovo X come ultima riga della matrice precedente. La prima maxi riga era risultati_integrazione e la seconda era X. Cioè stavi creando un array composto da 2 array risultati_integrazione e X, invece di aggiungere l'array X agli ALTRI precedenti array!!!! risultati_integrazione viene considerato come un array nel suo complesso invece che come una collezione di tanti array
#%%
# SISTEMA PRINCIPALE DOVE CI SONO TUTTE LE ISTRUZIONI
# caso Goodwin (Legge di Say: k[pi/ni] = 1 - omega)

def sistema_2_eq_mio(X):
    #spacchettamento del vettore    
    omega = X[0]
    lanbda = X[1]
    
    pi = 1 - omega
    
    domegadt = omega * (functionalform(A,B,C,D,lanbda) - alfa)
    dlanbdadt = lanbda * ((pi/ni) - alfa - beta - gamma)
    vettore_derivate = np.array([domegadt,dlanbdadt])
    return vettore_derivate
#%%
# PARAMETRI e COSTANTI DEL MODELLO 

A = 0.0000641
B = 1
C = 1
D = 0.0400641
E = 0.0175
F = 0.53 
G = 6
H = 0.065

a = 1
N = 100
alfa = 0.015
beta = 0.035
gamma = 0.02
ni = 3

#omega0 = 0.96  # uniche mostrate da Keen
#lanbda0 = 0.9  # uniche mostrate da Keen

#omega0 = 1 - ni*(alfa + beta + gamma) # punto fisso equilibrio stabile
#lanbda0 = (B/C) - ( np.sqrt(A/(alfa + D)) )/C # punto fisso equilibrio stabile per cui functionalform(A,B,C,D,lanbda0)=alfa anche se dopo tanti zeri differisce per cui potrebbe non andare bene
# ovvero questi sono:
omega0eq1 = 0.79
lanbda0eq1 = 0.965881121906885
#%%
# PUNTI FISSI
# Per trovare automaticamente gli zeri del sistema di equazioni differenziali

# define the system in this way (asuming a predator-prey-system with no negative values)
# to avoid interference omega = u (predator) and lambda = v (prey)
# Qua devi riscrivere le eq del sistema con i nuovi nomi delle variabili per evitare di usare le stesse lettere già usate per non creare confusioni. Visto che queste saranno simboli
# Se ci fosse bisogno di una terza variabile per una terza eq. differenziale usa "w" e "W"

u, v = sm.symbols('u, v', negative=False)
U = u * (functionalform(A,B,C,D,v) - alfa)
V = v * (((1-u)/ni) - alfa - beta - gamma) # non funziona bene, dà un sacco di zeri
#V = -v*(1-u)*((1-u-u)/2) # questo funziona bene, ma non so perchè ???NO???
#V = -(u**2)*v + (3/2)*u*v -(1/2)*v  # non funziona

# use sympy's way of setting equations to zero
UEqual = sm.Eq(U, 0)
VEqual = sm.Eq(V, 0)

# compute fixed points
equilibria = sm.solve( (UEqual, VEqual), u, v )
print(equilibria)
#%%    
# INTEGRAZIONE

# definisco le condizioni iniziali di ogni flusso ma è ridondante perchè potrei direttamente farlo nella chiamata di flusso_traiettoria_singola(x_0,y_0,..,..) mettendo direttamente i numeri
# cerchio piccolo
x0_1 = 0.96
y0_1 = 0.9
# asse x da destra verso l'origine poi va verso l'alto lungo asse y
x0_2 = 0.9    
y0_2 = 0
# cerchio più vicino agli assi
x0_3 = 0.05
y0_3 = 0.05
# cerchio esterno  PROVA
x0_4 = 1.4
y0_4 = 0.3
# altro cerchio
x0_5 = 0.25
y0_5 = .25
# altro cerchio  OK
x0_6 = 0.3
y0_6 = 0.3 
# altro cerchio   PROVA
x0_7 = 0.8
y0_7 = 0.4
# altro cerchio  OK
x0_8 = 0.6
y0_8 = 0.6
# altro cerchio  OK
x0_9 = 0.65
y0_9 = 0.65
# altro cerchio  OK
x0_10 = 0.75
y0_10 = 0.75  
# alto sinistra verso il basso che va all'origine
x0_11 = 1
y0_11 = 0.2
# esterno sotto
#x0_12 = 1  
#y0_12 = 0.5

#dt = 0.01
#numero_di_punti_tempo_voluti = 5000  # deve essere un numero tale che tempo_finale venga un numero intero
dt=0.005
#numero_di_punti_tempo_voluti = 50000
#numero_di_punti_tempo_voluti = 25000
numero_di_punti_tempo_voluti = 15000
tempo_finale = numero_di_punti_tempo_voluti * dt
t = np.linspace(0,tempo_finale,numero_di_punti_tempo_voluti) # per fare il grafico dopo, l'ultimo è proprio tempo_finale

# ovviamente tutto ciò si poteva saltare mettendo direttamente i numeri al posto di x0_i e y0_i definendoli qui e non sopra che è ridondante
x_1_flusso, y_1_flusso = flusso_traiettoria_singola(x0_1,y0_1,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio)
x_2_flusso, y_2_flusso = flusso_traiettoria_singola(x0_2,y0_2,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio)
x_3_flusso, y_3_flusso = flusso_traiettoria_singola(x0_3,y0_3,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio)
x_4_flusso, y_4_flusso = flusso_traiettoria_singola(x0_4,y0_4,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio)
x_5_flusso, y_5_flusso = flusso_traiettoria_singola(x0_5,y0_5,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio)
x_6_flusso, y_6_flusso = flusso_traiettoria_singola(x0_6,y0_6,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio)
x_7_flusso, y_7_flusso = flusso_traiettoria_singola(x0_7,y0_7,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio)
x_8_flusso, y_8_flusso = flusso_traiettoria_singola(x0_8,y0_8,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio)
x_9_flusso, y_9_flusso = flusso_traiettoria_singola(x0_9,y0_9,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio)
x_10_flusso, y_10_flusso = flusso_traiettoria_singola(x0_10,y0_10,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio)
x_11_flusso, y_11_flusso = flusso_traiettoria_singola(x0_11,y0_11,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio)
#x_12_flusso, y_12_flusso = flusso_traiettoria_singola(x0_12,y0_12,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio)
#%%
# GRAFICO SPAZIO DELLE FASI
font = {'family':'serif','color':'darkred','weight':'normal','size': 10.7,}
# plt.title('Damped exponential decay', fontdict=font)
# plt.text(2, 0.65, r'$\cos(2 \pi t) \exp(-t)$', fontdict=font)
# plt.xlabel('time (s)', fontdict=font)
# plt.ylabel('voltage (mV)', fontdict=font)


# FRECCE DEL CAMPO VETTORIALE

# creo la griglia, in realtà saranno 2 matrici una con coordinata x e l'altra y così prese elemento per elemento danno il punto (x_i,y_i)
x_griglia = np.linspace(0, 2, 20) # buono anche solo 14 però poi le frecce non sono tutte ortogonali come dovrebbe invece essere
y_griglia = np.linspace(0, 2, 20)

# il comando sotto è interessante perchè assegna automaticamente i termini in ordine: es. gaia=np.array([3,4]) poi con hippie1,hippie2=gaia si ha hippie1=gaia[0]=3 e hippie2=gaia[1]=4
X1 , Y1  = np.meshgrid(x_griglia, y_griglia)    # create a grid
# ritorna una prima matrice che ha come righe il vettore x ripetuto tante volte quante la lunghezza del vettore y, cioè len(y) volte, e una seconda matrice che ha il vettore y come colonne ripetuto tante volte quanta la lunghezza del vettore x, len(x)
DX1, DY1 = sistema_2_eq_mio([X1, Y1])    # compute growth rate on the grid
# la moltiplicazione * usata in "sistema_II_semplificato" è quella normale e se X come in questo caso è composto da 2 matrici, quando spacchetta X assegna la prima matrice (cioè le coordinate x su tutta la griglia) a x e la seconda (cioè la griglia per le coordinate y) a y, quando poi si trova a fare x*y moltiplica elemento per elemento entrambe le matrici senza fare il prodotto matriciale! Perciò le due matrici devono avere entrambe la stessa dimensione, tutte e due devono essere nxn
M = (np.hypot(DX1, DY1))                        # norm growth rate 
M[ M == 0] = 1.                                 # avoid zero division errors 
DX1 /= M                                        # normalize each arrows
DY1 /= M
# se fai M.shape dà (20,20) è infatti una matrice 20x20 siccome hai usato due volte x = np.linspace(0, 2, 20) o x e y, ma sono uguali x=y=np.linspace(0, 2, 20)

#plt.quiver(X1, Y1, DX1, DY1, color='deepskyblue',pivot='mid')
#plt.quiver(X1, Y1, DX1, DY1, M, pivot='mid') # M gli dà una scala di colori che sono i numeri


plt.plot(x_1_flusso, y_1_flusso,'b')
#plt.plot(x_2_flusso, y_2_flusso,'b')
#plt.plot(x_3_flusso, y_3_flusso,'b')
plt.plot(x_4_flusso[1:6300], y_4_flusso[1:6300],'b')
#plt.plot(x_5_flusso, y_5_flusso,'b')
#plt.plot(x_6_flusso, y_6_flusso,'b')
plt.plot(x_7_flusso, y_7_flusso,'b')
plt.plot(x_8_flusso, y_8_flusso,'b')
plt.plot(x_9_flusso, y_9_flusso,'b')
plt.plot(x_10_flusso, y_10_flusso,'b')
#plt.plot(x_11_flusso, y_11_flusso,'b')
#plt.plot(x_12_flusso, y_12_flusso,'b')

# Nel grafico che sto creando metto anche tutti punti fissi. Ci sono 2 modi possibili per farlo
# versione "automatica" però non puoi variare i colori:
#for point in equilibria:
#    plt.plot(point[0],point[1],"royalblue", marker = "o", markersize = 7.0)

# versione "manuale":
plt.plot(0,0,'ob')
plt.plot(omega0eq1,lanbda0eq1,'o',color='royalblue')   # condizione iniziale o meglio punto del centro dei cicli. Basta allontanarsi poco da questo punto con la condizione iniziale che si formano i cycle!

ip=500
ip1=100
ip2=300
ip3=1820
ip4=2900
ip5=4250
plt.arrow(x_1_flusso[ip], y_1_flusso[ip],x_1_flusso[ip+1]-x_1_flusso[ip], y_1_flusso[ip+1]-y_1_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
#plt.arrow(x_2_flusso[ip], y_2_flusso[ip],x_2_flusso[ip+1]-x_2_flusso[ip], y_2_flusso[ip+1]-y_2_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
#plt.arrow(x_3_flusso[ip1], y_3_flusso[ip1],x_3_flusso[ip1+1]-x_3_flusso[ip1], y_3_flusso[ip1+1]-y_3_flusso[ip1],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_4_flusso[ip], y_4_flusso[ip],x_4_flusso[ip+1]-x_4_flusso[ip], y_4_flusso[ip+1]-y_4_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
#plt.arrow(x_5_flusso[ip], y_5_flusso[ip],x_5_flusso[ip+1]-x_5_flusso[ip], y_5_flusso[ip+1]-y_5_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
#plt.arrow(x_6_flusso[ip], y_6_flusso[ip],x_6_flusso[ip+1]-x_6_flusso[ip], y_6_flusso[ip+1]-y_6_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_7_flusso[ip], y_7_flusso[ip],x_7_flusso[ip+1]-x_7_flusso[ip], y_7_flusso[ip+1]-y_7_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_8_flusso[ip], y_8_flusso[ip],x_8_flusso[ip+1]-x_8_flusso[ip], y_8_flusso[ip+1]-y_8_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_9_flusso[ip], y_9_flusso[ip],x_9_flusso[ip+1]-x_9_flusso[ip], y_9_flusso[ip+1]-y_9_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_10_flusso[ip3], y_10_flusso[ip3],x_10_flusso[ip3+1]-x_10_flusso[ip3], y_10_flusso[ip3+1]-y_10_flusso[ip3],shape='full', lw=0, length_includes_head=True, head_width=.04)
#plt.arrow(x_11_flusso[ip2], y_11_flusso[ip2],x_11_flusso[ip2+1]-x_11_flusso[ip2], y_11_flusso[ip2+1]-y_11_flusso[ip2],shape='full', lw=0, length_includes_head=True, head_width=.04)
#plt.arrow(x_12_flusso[ip1], y_12_flusso[ip1],x_12_flusso[ip1+1]-x_12_flusso[ip1], y_12_flusso[ip1+1]-y_12_flusso[ip1],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(0.5, 0, -0.01, 0,shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(1.1, 0, -0.01, 0,shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(0, 0.8, 0, 0.01,shape='full', lw=0, length_includes_head=True, head_width=.05)
plt.arrow(0, 0.4, 0, 0.01,shape='full', lw=0, length_includes_head=True, head_width=.05)

#plt.legend(loc='best', prop={'size': 7})
#plt.title('spazio delle fasi Goodwin con Legge di Say',fontdict = font)
plt.xlabel('quota salari $\omega$')
plt.ylabel('occupazione $\lambda$')
plt.axis([0, 1.4, 0, 1.1]) # X=1.4 buono
plt.grid()
#plt.savefig('relazione corretta Goodwin Say spazio delle fasi.jpg', dpi=850, transparent=False)
plt.show()
#%%
# ANDAMENTO TEMPORALE PARAGONE DIVERSE DINAMICHE COME PARE VOLERE TURCHETTI
# Forse aggiungere legenda col punto iniziale o informazioni ad esempio "dinamica sella"

# grafico temporale con una bella dinamica ma forse non è quello che vuole Turchetti, oppure ne vuole tante a confronto 
plt.plot(t, x_1_flusso,'r',label='quota salari $\omega(t)$')
plt.plot(t, y_1_flusso,'black',label='occupazione $\lambda(t)$')
plt.legend(loc='lower right', prop={'size': 6})
#plt.title('andamento temporale',fontdict = font)
plt.xlabel('tempo t')
plt.ylabel('$\omega, \lambda$')
plt.grid()
plt.xlim(0,60)
#plt.savefig('relazione corretta Goodwin Say andamento temporale 1.jpg', dpi=650, transparent=False)
plt.show()

plt.plot(t, x_4_flusso,'r',label='quota salari $\omega(t)$')
plt.plot(t, y_4_flusso,'k',label='occupazione $\lambda(t)$')
plt.legend(loc='best',prop={'size': 6})
#plt.legend(loc='lower right', prop={'size': 6})
#plt.title('andamento temporale',fontdict = font)
plt.xlabel('tempo t')
plt.ylabel('$\omega, \lambda$')
plt.grid()
plt.xlim(0,60)
#plt.savefig('relazione corretta Goodwin Say andamento temporale 2.jpg', dpi=650, transparent=False)
plt.show()

plt.plot(t, x_7_flusso,'r',label='quota salari $\omega(t)$')
plt.plot(t, y_7_flusso,'k',label='occupazione $\lambda(t)$')
plt.legend(loc='best',prop={'size': 6})
#plt.title('andamento temporale',fontdict = font)
plt.xlabel('tempo t')
plt.ylabel('$\omega, \lambda$')
plt.grid()
plt.xlim(0,60)
#plt.savefig('relazione corretta Goodwin Say andamento temporale 3.jpg', dpi=650, transparent=False)
plt.show()

plt.plot(t, x_8_flusso,'r',label='quota salari $\omega(t)$')
plt.plot(t, y_8_flusso,'k',label='occupazione $\lambda(t)$')
plt.legend(loc='best',prop={'size': 6})
#plt.title('andamento temporale',fontdict = font)
plt.xlabel('tempo t')
plt.ylabel('$\omega, \lambda$')
plt.grid()
plt.xlim(0,60)
#plt.savefig('relazione corretta Goodwin Say andamento temporale 4.jpg', dpi=650, transparent=False)
plt.show()

plt.plot(t, x_10_flusso,'r',label='quota salari $\omega(t)$')
plt.plot(t, y_10_flusso,'k',label='occupazione $\lambda(t)$')
plt.legend(loc='best',prop={'size': 6})
#plt.title('andamento temporale',fontdict = font)
plt.xlabel('tempo t')
plt.ylabel('$\omega, \lambda$')
plt.grid()
plt.xlim(0,60)
#plt.savefig('relazione corretta Goodwin Say andamento temporale 5.jpg', dpi=650, transparent=False)
plt.show()
#%%
#%%
# DISEGNO NULLCLINE e GRAFICO GLI ZERI DELLE FUNZIONI 

# nullcline e zeri della funzione f(omega) che compare in dlambda/dt come vuole Turchetti
assex=np.linspace(-1,2,1000)

#    domegadt = omega * (functionalform(A,B,C,D,lanbda) - alfa)
#    dlanbdadt = lanbda * (((1-omega)/ni) - alfa - beta - gamma)

# Grafico omega - f(omega)
plt.plot(assex, (1-assex)/ni - alfa - beta - gamma,'blue')
plt.plot([0,0],[-1,2],'k')
plt.plot([-1,2],[0,0],'k')
plt.plot(omega0eq1,0,'ob')
plt.legend(loc='best',prop={'size': 6})
#plt.title('grafico omega f(omega) NON è la nullcline di $\dot{\lambda}$',fontdict = font)
plt.xlabel('$\omega$')
plt.ylabel('$f(\omega)$')
plt.grid()
#plt.xlim(-1,2)
#plt.savefig('relazione corretta Goodwin Say grafico omega f(omega).jpg', dpi=650, transparent=False)
plt.show()

# Grafico lambda - g(lambda)
plt.plot(assex, functionalform(A,B,C,D,assex) - alfa,'red')
plt.plot([0,0],[-1,3],'k')
plt.plot([-1,3],[0,0],'k')
plt.plot([1,1],[-1,3],'k')
plt.plot(lanbda0eq1,0,'or')
plt.legend(loc='best',prop={'size': 6})
#plt.title('grafico lambda g(lambda) NON è la nullcline di $\dot{\omega}$',fontdict = font)
plt.xlabel('$\lambda$')
plt.ylabel('$g(\lambda)$')
plt.grid()
plt.xlim(-0.1,1.2)
plt.ylim(-0.2,1.2)
#plt.savefig('relazione corretta Goodwin Say grafico lambda g(lambda).jpg', dpi=650, transparent=False)
plt.show()

# Grafico lambda - w[lambda]
plt.plot(assex, functionalform(A,B,C,D,assex),'blue')
plt.plot([0,0],[-1,3],'k')
plt.plot([-1,3],[0,0],'k')
plt.plot([1,1],[-1,3],'k')
plt.plot(1-(A/D)**0.5,0,'ob')
#plt.plot(1+(A/D)**0.5,0,'ob')
plt.legend(loc='best',prop={'size': 6})
#plt.title('funzione risposta lavoratori grafico lambda w[lambda]',fontdict = font)
plt.xlabel('$\lambda$')
plt.ylabel('$w[\lambda]$')
plt.grid()
plt.xlim(-0.1,1.2)
plt.ylim(-0.2,1.2)
#plt.savefig('relazione corretta Goodwin Say grafico lambda w[lambda].jpg', dpi=650, transparent=False)
plt.show()

# Faccio vedere a scopo dimostrativo che il campo vettoriale è perpendicolare alle nullcline vicino ad esse
# Draw Nullclines and Quiver plot
# plot nullclines
plt.plot([-5,5],[(B/C) - ( np.sqrt(A/(alfa + D)) )/C ,(B/C) - ( np.sqrt(A/(alfa + D)) )/C ], 'b-', lw=2, label='$\dot{\omega}$ nullcline')
plt.plot([-5,5],[(B/C) + ( np.sqrt(A/(alfa + D)) )/C ,(B/C) + ( np.sqrt(A/(alfa + D)) )/C ], 'b-', lw=2)
plt.plot([1-ni*(alfa+beta+gamma),omega0eq1],[-5,5], 'r-', lw=2, label='$\dot{\lambda}$ nullcline')
plt.plot(functionalform(A,B,C,D,assex) - alfa,assex,'green') # questa in realtà andrebbe tolta perchè N.B.: g(lambda) NON è la nullcline, essa sono i due valori di lambda* infatti la nullcline di lambda non dipende da omega, è un retta orizzontale quindi
plt.plot(assex, (1-assex)/ni - alfa - beta - gamma,'green') #  questa in realtà andrebbe tolta perchè N.B.: f(omega) NON è la nullcline, essa è il valore di omega* infatti la nullcline di omega non dipende da lambda, è un retta verticale quindi
plt.plot([0,0],[-5,5],'k',lw=2)
plt.plot([-5,5],[0,0],'k',lw=2)
plt.plot([-5,5],[1,1],'--k')
plt.legend(loc='best',prop={'size': 9})
plt.xlabel('$\omega$')
plt.ylabel('$\lambda$')
plt.xlim(-0.1,1.2)
plt.ylim(-0.1,1.2)
# plot fixed points
for point in equilibria:
    plt.plot(point[0],point[1],"blue", marker = "o", markersize = 6.0)
plt.title("Quiverplot with nullclines",fontdict = font)
plt.quiver(X1, Y1, DX1, DY1, color='deepskyblue',pivot='mid')
plt.grid()
# carina da tenere che si vede meglio anche se NON LA DEVI METTERE NELLA RELAZIONE
#plt.savefig('relazione corretta Goodwin Say nullcline spazio delle fasi.jpg', dpi=650, transparent=False)
plt.show()

# Grafico omega,lambda - f(omega),g(lambda)
plt.plot(assex,(((1-assex)/ni) - alfa - beta - gamma),'blue',label='$f(\omega)$')
plt.plot(assex,functionalform(A,B,C,D,assex)-alfa,'red', label='$g(\lambda)$')
plt.plot([-1,2],[0,0],'k', [0,0],[-1,2],'k',[omega0eq1,lanbda0eq1],[0,0],'or')
# oppure tutto insieme ma non funziona il comando label non so perchè!
#plt.plot(assex,functionalform(A,B,C,D,assex)-alfa,'red', assex,(((1-assex)/ni) - alfa - beta - gamma),'blue', [-1,2],[0,0],'k', [0,0],[-1,2],'k',[omega0eq1,lanbda0eq1],[0,0],'or')
plt.legend(loc='best',prop={'size': 8})
plt.axis([-0.2, 1.2, -0.2, 1.3])
plt.xlabel('$\omega, \lambda$')
plt.ylabel('$f(\omega), g(\lambda)$')
#plt.xlim(-1,2)
#plt.ylim(-1,1)
plt.grid()
plt.show()
#%%
#%%
# CALCOLO AUTOVALORI DEI PUNTI FISSI PER STUDIARE STABILITA'. N.B.: NON FIDARSI DI QUELLO CHE DICE CHE SONO. Inoltre devi inserire manualmente i punti fissi trovati sopra e la matrice Jacobiana generica
# non funziona molto bene e per i punti "centri" non funziona proprio! Siamo messi bene...

def eigenvalues(x,y):
#    domegadt = omega * (functionalform(A,B,C,D,lanbda) - alfa)
    a11 = functionalform(A,B,C,D,y)-alfa         # differentiated with respect to x
    a12 = x* ( (2*A*C)/(B-C*y)**3 )              # differentiated with respect to y
    
#    dlanbdadt = lanbda * (((1-omega)/ni) - alfa - beta - gamma)
    a21 = -(y/ni)                                # differentiated with respect to x
    a22 = ((1-x)/ni)-alfa-beta-gamma             # differentiated with respect to y

    tr = a11 + a22
    det = a11*a22 - a12*a21
    lambda1 = (tr - sqrt(tr**2 - 4*det))/2
    lambda2 = (tr + sqrt(tr**2 - 4*det))/2
    print('Check the fixed point  %s, %s' % (x,y)) 
    print('The real part of the first eigenvalue is %s' %lambda1.real)
    print('The real part of the second eigenvalue is %s' % lambda2.real)    
    
    if (lambda1.real < 0 and lambda2.real < 0):
        print('The fixed point in %s, %s is a sink. It is stable' % (x,y))
    if (lambda1.real > 0 and lambda2.real > 0):
        print('The fixed point in %s, %s is a source. It is unstable' % (x,y))
    if (lambda1.real > 0 and lambda2.real < 0):
        print('The fixed point in %s, %s is a saddle. It is unstable' % (x,y))
    if (lambda1.real < 0 and lambda2.real > 0):
        print('The fixed point in %s, %s is a saddle. It is unstable' % (x,y))
#    if (lambda1.real < 0 and lambda2.real > 0):
#        print('The fixed point in %s, %s is unstable' % (x,y))
    print('----------------------------')
    return lambda1 , lambda2

equilibria_manuale=np.array([[0,0],[omega0eq1,lanbda0eq1],[omega0eq1,(B/C) + ( np.sqrt(A/(alfa + D)) )/C ]]) # perchè con equilibria non funziona siccome ci sono tutti i decimali, dice che non può fare sqrt() di un float...non so perchè! 
# iterate through list of fixed points
for x,y in equilibria_manuale:
    eigenvalues(x,y)    

#for item in equilibria:
#    eigenvalues(item[0],item[1])
#%%
# AUTOVALORI JACOBIANA SISTEMA 2 EQUAZIONI SENZA FUNZIONE INVESTIMENTI: k[pi/ni] = 1-omega
def jacobiana_sistema_2_eq_basic_goodwin(omegaeq,lanbdaeq):
    primo = functionalform(A,B,C,D,lanbdaeq) - alfa
    secondo = omegaeq * ( (2*A*C)/(B-C*lanbdaeq)**3 )
    terzo = (-(lanbdaeq/ni)) +0.0 # sennò viene -0.0 per lanbdaeq=0 
    quarto = ((1-omegaeq)/ni) - alfa - beta - gamma
    primo=round(primo,9) # sennò non viene zero perchè tiene 17 cifre decimali!
    secondo=round(secondo,9)
    terzo=round(terzo,9)
    quarto=round(quarto,9)
    #quarto=round(quarto) # perchè sennò viene -0.0 quando omegaeq è quello che annulla "quarto". N.B.: TTENZIONE CHE IN TUTTI GLI ALTRI CASI TI FA VENIRE RISULTATI SBAGLIATI!!! Meglio lasciare sempre commentato, non usare
    return np.array([[primo, secondo],[terzo, quarto]])

omegaeq=0 # è un punto sella!
lanbdaeq=0 # è un punto sella!
#omegaeq = 1 - ni*(alfa + beta + gamma) # è un centro
#lanbdaeq = (B/C) - ( np.sqrt(A/(alfa + D)) )/C # è un centro, per cui functionalform(A,B,C,D,lanbda0)=alfa anche se dopo tanti zeri differisce per cui potrebbe non andare bene
# ovvero questi sono:
#omegaeq = 0.79
#lanbdaeq = 0.965881121906885


jacobiana_sistema_basic_goodwin = jacobiana_sistema_2_eq_basic_goodwin(omegaeq,lanbdaeq)
traccia_sistema_basic_goodwin = jacobiana_sistema_basic_goodwin[0,0] + jacobiana_sistema_basic_goodwin[1,1]
determinante_sistema_basic_goodwin = jacobiana_sistema_basic_goodwin[0,0]*jacobiana_sistema_basic_goodwin[1,1] - jacobiana_sistema_basic_goodwin[0,1]*jacobiana_sistema_basic_goodwin[1,0]

autovalori_sistema_basic_goodwin, autovettori_sistema_basic_goodwin = LA.eig(jacobiana_sistema_basic_goodwin)

print('\nLa matrice Jacobiana del sistema a 2 equazioni SENZA funzione investimenti k[pi/ni] è:\n',jacobiana_sistema_basic_goodwin,'\n')
print('Il determinante det(J) è:\n', determinante_sistema_basic_goodwin)
print('La traccia tr(J) è:\n', traccia_sistema_basic_goodwin)
print('Gli autovalori della matrice J sono:\n', autovalori_sistema_basic_goodwin)
print('Una possibile scelta per i 2 autovettori di J è:\n',autovettori_sistema_basic_goodwin)

plt.plot(determinante_sistema_basic_goodwin, traccia_sistema_basic_goodwin, marker='o', linestyle='--', color='r')
plt.title('Grafico determinante e traccia della matrice Jacobiana \n dalla quale si ricavano gli autovalori per quel punto fisso')
plt.xlabel('det(J)')
plt.ylabel('tr(J)')
plt.xlim(-10,10)
plt.ylim(-10,10)
plt.grid()
plt.show()
#%%
#%%
# STREAMPLOT per disegnare il ritratto di fase ma NON USARE perchè è meglio quiver!!!! 

plt.streamplot(X1, Y1, DX1, DY1, linewidth=1, density=2, arrowstyle='->', arrowsize=1.5)
#plt.streamplot(X1, Y1, DX1, DY1,color=M, linewidth=2, cmap=plt.cm.autumn)
# alla fine puoi aggiungere anche il comando density = 1.4
plt.xlim(0,2)
plt.ylim(0,2)
plt.grid()
plt.show()
#%%
#%%
#%%
#%%
#%%
#                                      ALTRO PROGRAMMA
#%%
#%%
#%%
#%%
#%%
# SISTEMA PRINCIPALE DOVE CI SONO TUTTE LE ISTRUZIONI
# caso Goodwin (Funzione Investimenti: k[pi/ni])

def sistema_2_eq_mio_con_funzione_investimenti(X):
    #spacchettamento del vettore    
    omega = X[0]
    lanbda = X[1]
    
    pi = 1 - omega
    
    domegadt = omega * (functionalform(A,B,C,D,lanbda) - alfa)
    dlanbdadt = lanbda * (((functionalform(E,F,G,H,(pi/ni)))/ni) - alfa - beta - gamma)
    vettore_derivate = np.array([domegadt,dlanbdadt])
    return vettore_derivate
#%%
# INFO DEL MODELLO. Abbastanza inutile, forse si può cancellare!

#omega0 = 0.85  # anche se queste condizioni iniziali sono molto vicine al centro omega0eq e lanbda0eq, il cycle si forma lo stesso essendo appunto un punto fisso centro. Il cerchio più le condizioni iniziali si allontanano dal centro diventa una elisse, per condizioni troppo vicine a (0,0) che il nodo stabile il sistema collassa, ciò accade già per (0.6,0.6) condzione iniziale
#lanbda0 = 0.97 # anche se queste condizioni iniziali sono molto vicine al centro omega0eq e lanbda0eq, il cycle si forma lo stesso essendo appunto un punto fisso centro. Il cerchio più le condizioni iniziali si allontanano dal centro diventa una elisse, per condizioni troppo vicine a (0,0) che il nodo stabile il sistema collassa, ciò accade già per (0.6,0.6) condzione iniziale
#omega0 = 0.77  # Il punto (0,0) è un nodo stabile quindi attrae ad esso tutte le condizioni iniziali vicine come (0.76,0.76) ad esempio. Per (0.77,0.77) come condizione iniziale si ha un comportamento interessante e il cycle comincia ad emergere
#lanbda0 = 0.77  # Il punto (0,0) è un nodo stabile quindi attrae ad esso tutte le condizioni iniziali vicine come (0.76,0.76) ad esempio. Per (0.77,0.77) come condizione iniziale si ha un comportamento interessante e il cycle comincia ad emergere
omega0 = 0.96
lanbda0 = 0.9
#omega0 = 1 - ((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ) # punto fisso equilibrio stabile
#lanbda0 = (B/C) - ( np.sqrt(A/(alfa + D)) )/C # punto fisso equilibrio stabile per cui functionalform(A,B,C,D,lanbda0)=alfa anche se dopo tanti zeri differisce per cui potrebbe non andare bene
# ovvero questi sono:
omega0eq = 0.8611312447773782
lanbda0eq = 0.965881121906885
#%%
# PUNTI FISSI
# Per trovare automaticamente gli zeri del sistema di equazioni differenziali

# define the system in this way (asuming a predator-prey-system with no negative values)
# to avoid interference omega = u (predator) and lambda = v (prey)
# Qua devi riscrivere le eq del sistema con i nuovi nomi delle variabili per evitare di usare le stesse lettere già usate per non creare confusioni. Visto che queste saranno simboli
# Se ci fosse bisogno di una terza variabile per una terza eq. differenziale usa "w" e "W"

u, v = sm.symbols('u, v', negative=False)
U = u * (functionalform(A,B,C,D,v) - alfa)
V = v * (((functionalform(E,F,G,H,((1-u)/ni)))/ni) - alfa - beta - gamma) # non funziona bene, dà un sacco di zeri
#V = -v*(1-u)*((1-u-u)/2) # questo funziona bene, ma non so perchè ???NO???
#V = -(u**2)*v + (3/2)*u*v -(1/2)*v  # non funziona

# use sympy's way of setting equations to zero
UEqual = sm.Eq(U, 0)
VEqual = sm.Eq(V, 0)

# compute fixed points
equilibria = sm.solve( (UEqual, VEqual), u, v )
print(equilibria)
# TI DA' LA BELLEZZA DI 5 PUNTI DI CUI PERO' 2 HANNO LA STESSA lambda>1 !!!
#%%    
# INTEGRAZIONE

# definisco le condizioni iniziali di ogni flusso ma è ridondante perchè potrei direttamente farlo nella chiamata di flusso_traiettoria_singola(x_0,y_0,..,..) mettendo direttamente i numeri

# I QUADRANTE: tutte verso l'origine stabile + sella disegnare le varietà
# verso l'origine dall'alto
x0_1 = 0.3
y0_1 = 0.99
# cerchi di beccare il braccio che sopra il punto scende giù a sinistra verso il punto sella
x0_20 = 0.471 # 0.47 buonino e 0.475 decisamente troppo. 0.471 è lui!
y0_20 = 0.996
# sella varietà instabile della sella che credo faccia parta dal centro in basso, faccia il giro del pesce lungo lambda=1 e credo ritorni di nuovo in basso dove era partita dall'infinito, è come se facesse una D larghissima
x0_2 = 0.72774 # con 0.74 ha uno strano comportamento: becca l'omega asintoto della funzione k[] cioè 1-((F*ni)/G)=0.735 e va a lambda all'infinito e poi torna giù verso lambda=1
y0_2 = 0.01    # funziona tenere!
# parte in basso al centro fa una gobba e va all'origine
x0_3 = 0.726
y0_3 = 0.01
#x0_3 = 0.61
#y0_3 = 0.97
# braccio instabile della sella: quello che non va all'origine. Il 4 va bene!
x0_4 = 0.6
y0_4 = 0.96

# II QUADRANTE: i cerchi del centro
# cerchio piccolo. QUELLO SCELTO DA STEVE KEEN
x0_5 = 0.96
y0_5 = 0.9
# da in basso a destra fa tutto il cerchio lungo
x0_6 = 1.4
y0_6 = 0.1
# cerchio interno successivo a quello di Keen
x0_8 = 1
y0_8 = 0.7
# attento che andando omega all'infinito lungo lambda=1 poi fa 2 giri perciò 2 cerchi che vedi nell'immagine e che escono da omega=1.4 in realtà sono lo stesso flusso che non si ricongiunge con le condizioni iniziali come dovrebbe perchè va all'infinito
x0_9 = 1.5
y0_9 = 0.5

# III QUADRANTE: parte non nel modello, i cerchi del centro in alto
# il cerchio pi piccolo
x0_10 = 0.6
y0_10 = 1.2
# cerchio più esterno, scende dall'alto a sinistra e va a lanbda=1
x0_11 = 0.05
y0_11 = 2.3  # prima era 2
# cerchio molto bello uno dei più esterni ed è completo
x0_12 = 0.6
y0_12 = 2
# cerchio esterno attenzione che parte ma dove inizia l'immagine ma poi fa il giro e ci torna perciò rischia di venire più spesso dalla parte da lanbda=2 dove passa 2 volte e meno sottile a lanbda>2 dove passa una sola 
x0_13 = 0.3
y0_13 = 2.2
# altro cerchio a metà
x0_14 = 0.6
y0_14 = 1.5


# IV QUADRANTE: parte non nel modello, la sella in alto difficile diseganre varietà
# lo facciamo scendere dall'alto per cercare di beccare la sella     ???
x0_15 = 1.5193 # 1.5193 miracolosamente funziona! 1.52 non va bene è troppo, 1.51 era poco
y0_15 = 2.295  # 2.295 miracolosamente funziona! 2.3 non andava bene come anche 2.29 ma buonino
x0_16 = 1.2
y0_16 = 2.3
#x0_16 = 0.85 # PARE CHE NON SERVA ALLORA USARE PER ALTRO!
#y0_16 = 1.04
# prova per beccare un braccio della sella partendo vicinissimo a destra del punto. QUESTO FUNZIONA MA DI FATTO NON E' VISIBILE NEL GRAFICO, TROPPO PICCOLO
x0_17 = 0.87
y0_17 = 1.03
# ottimi tenere
x0_18 = 0.8
y0_18 = 1.01
x0_19 = 1.5
y0_19 = 1.8
# cerco l'ultimo braccio della sella la varietà stabile che ci va da sotto a sinistra
x0_7 = 0.8138 # 0.81 è poco e 0.82 è troppo
y0_7 = 1.01  # 1.01 benino, ma sarebbe buono partire un po' più in basso però è impossibile


#dt = 0.01
#numero_di_punti_tempo_voluti = 5000  # deve essere un numero tale che tempo_finale venga un numero intero
#dt=0.005
dt=0.004 # o forse meglio 0.002
#numero_di_punti_tempo_voluti = 50000
#numero_di_punti_tempo_voluti = 25000
numero_di_punti_tempo_voluti = 15000
tempo_finale = numero_di_punti_tempo_voluti * dt
t = np.linspace(0,tempo_finale,numero_di_punti_tempo_voluti) # per fare il grafico dopo, l'ultimo è proprio tempo_finale

# ovviamente tutto ciò si poteva saltare mettendo direttamente i numeri al posto di x0_i e y0_i definendoli qui e non sopra che è ridondante
x_1_flusso, y_1_flusso = flusso_traiettoria_singola(x0_1,y0_1,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_2_flusso, y_2_flusso = flusso_traiettoria_singola(x0_2,y0_2,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_3_flusso, y_3_flusso = flusso_traiettoria_singola(x0_3,y0_3,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_4_flusso, y_4_flusso = flusso_traiettoria_singola(x0_4,y0_4,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_5_flusso, y_5_flusso = flusso_traiettoria_singola(x0_5,y0_5,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_6_flusso, y_6_flusso = flusso_traiettoria_singola(x0_6,y0_6,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_7_flusso, y_7_flusso = flusso_traiettoria_singola(x0_7,y0_7,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_8_flusso, y_8_flusso = flusso_traiettoria_singola(x0_8,y0_8,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_9_flusso, y_9_flusso = flusso_traiettoria_singola(x0_9,y0_9,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_10_flusso, y_10_flusso = flusso_traiettoria_singola(x0_10,y0_10,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_11_flusso, y_11_flusso = flusso_traiettoria_singola(x0_11,y0_11,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_12_flusso, y_12_flusso = flusso_traiettoria_singola(x0_12,y0_12,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_13_flusso, y_13_flusso = flusso_traiettoria_singola(x0_13,y0_13,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_14_flusso, y_14_flusso = flusso_traiettoria_singola(x0_14,y0_14,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_15_flusso, y_15_flusso = flusso_traiettoria_singola(x0_15,y0_15,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_16_flusso, y_16_flusso = flusso_traiettoria_singola(x0_16,y0_16,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_17_flusso, y_17_flusso = flusso_traiettoria_singola(x0_17,y0_17,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_18_flusso, y_18_flusso = flusso_traiettoria_singola(x0_18,y0_18,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_19_flusso, y_19_flusso = flusso_traiettoria_singola(x0_19,y0_19,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
x_20_flusso, y_20_flusso = flusso_traiettoria_singola(x0_20,y0_20,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
#%%
# GRAFICO SPAZIO DELLE FASI
font = {'family':'serif','color':'darkred','weight':'normal','size': 10.7,}
# plt.title('Damped exponential decay', fontdict=font)
# plt.text(2, 0.65, r'$\cos(2 \pi t) \exp(-t)$', fontdict=font)
# plt.xlabel('time (s)', fontdict=font)
# plt.ylabel('voltage (mV)', fontdict=font)


# FRECCE DEL CAMPO VETTORIALE

# creo la griglia, in realtà saranno 2 matrici una con coordinata x e l'altra y così prese elemento per elemento danno il punto (x_i,y_i)
x_griglia = np.linspace(0, 2, 20) # buono anche solo 14 però poi le frecce non sono tutte ortogonali come dovrebbe invece essere
y_griglia = np.linspace(0, 2, 20)

# il comando sotto è interessante perchè assegna automaticamente i termini in ordine: es. gaia=np.array([3,4]) poi con hippie1,hippie2=gaia si ha hippie1=gaia[0]=3 e hippie2=gaia[1]=4
X1 , Y1  = np.meshgrid(x_griglia, y_griglia)    # create a grid
# ritorna una prima matrice che ha come righe il vettore x ripetuto tante volte quante la lunghezza del vettore y, cioè len(y) volte, e una seconda matrice che ha il vettore y come colonne ripetuto tante volte quanta la lunghezza del vettore x, len(x)
DX1, DY1 = sistema_2_eq_mio_con_funzione_investimenti([X1, Y1])    # compute growth rate on the grid
# la moltiplicazione * usata in "sistema_II_semplificato" è quella normale e se X come in questo caso è composto da 2 matrici, quando spacchetta X assegna la prima matrice (cioè le coordinate x su tutta la griglia) a x e la seconda (cioè la griglia per le coordinate y) a y, quando poi si trova a fare x*y moltiplica elemento per elemento entrambe le matrici senza fare il prodotto matriciale! Perciò le due matrici devono avere entrambe la stessa dimensione, tutte e due devono essere nxn
M = (np.hypot(DX1, DY1))                        # norm growth rate 
M[ M == 0] = 1.                                 # avoid zero division errors 
DX1 /= M                                        # normalize each arrows
DY1 /= M
# se fai M.shape dà (20,20) è infatti una matrice 20x20 siccome hai usato due volte x = np.linspace(0, 2, 20) o x e y, ma sono uguali x=y=np.linspace(0, 2, 20)

# ATTIVARE SE VUOI VEDERE LE FRECCE DEL CAMPO VETTORIALE ma non penso che Turchetti le voglia
#plt.quiver(X1, Y1, DX1, DY1, color='deepskyblue',pivot='mid')
#plt.quiver(X1, Y1, DX1, DY1, M, pivot='mid') # M gli dà una scala di colori che sono i numeri


plt.plot(x_1_flusso, y_1_flusso,'b')
plt.plot(x_2_flusso, y_2_flusso,'b')
plt.plot(x_3_flusso, y_3_flusso,'b')
plt.plot(x_4_flusso, y_4_flusso,'b')
plt.plot(x_5_flusso, y_5_flusso,'b')
plt.plot(x_6_flusso, y_6_flusso,'b')
plt.plot(x_7_flusso, y_7_flusso,'b') # questo è un braccio nel IV quadrante
plt.plot(x_8_flusso, y_8_flusso,'b')
#plt.plot(x_9_flusso, y_9_flusso,'b')
plt.plot(x_9_flusso[0:6500], y_9_flusso[0:6500],'b')
plt.plot(x_10_flusso, y_10_flusso,'b')
plt.plot(x_11_flusso, y_11_flusso,'b')
plt.plot(x_12_flusso, y_12_flusso,'b')
#plt.plot(x_13_flusso, y_13_flusso,'b')
plt.plot(x_13_flusso[0:2297], y_13_flusso[0:2297],'b')
plt.plot(x_14_flusso, y_14_flusso,'b')
plt.plot(x_15_flusso, y_15_flusso,'b')
plt.plot(x_16_flusso, y_16_flusso,'b')
plt.plot(x_17_flusso, y_17_flusso,'b')
plt.plot(x_18_flusso, y_18_flusso,'b')
plt.plot(x_19_flusso, y_19_flusso,'b')
plt.plot(x_20_flusso, y_20_flusso,'b')

# Nel grafico che sto creando metto anche tutti punti fissi. Ci sono 2 modi possibili per farlo
# versione "automatica" però non puoi variare i colori:
#for point in equilibria:
#    plt.plot(point[0],point[1],"royalblue", marker = "o", markersize = 7.0)

# versione "manuale":
plt.plot(0,0,'ob')
plt.plot(1 - ((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),(B/C) - ( np.sqrt(A/(alfa + D)) )/C,'o',color='royalblue')   # condizione iniziale o meglio punto del centro dei cicli. Basta allontanarsi poco da questo punto con la condizione iniziale che si formano i cycle!
plt.plot(1 - ((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),(B/C) - ( np.sqrt(A/(alfa + D)) )/C,'o',color='royalblue')
plt.plot(1 - ((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),(B/C) + ( np.sqrt(A/(alfa + D)) )/C,'o',color='royalblue') 
plt.plot(1 - ((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),(B/C) + ( np.sqrt(A/(alfa + D)) )/C,'o',color='royalblue')

ip=500
ip1=100
ip2=300
ip3=1820
ip4=3840
ip5=200
plt.arrow(x_1_flusso[ip], y_1_flusso[ip],x_1_flusso[ip+1]-x_1_flusso[ip], y_1_flusso[ip+1]-y_1_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_2_flusso[ip], y_2_flusso[ip],x_2_flusso[ip+1]-x_2_flusso[ip], y_2_flusso[ip+1]-y_2_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_3_flusso[ip3], y_3_flusso[ip3],x_3_flusso[ip3+1]-x_3_flusso[ip3], y_3_flusso[ip3+1]-y_3_flusso[ip3],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_4_flusso[ip3], y_4_flusso[ip3],x_4_flusso[ip3+1]-x_4_flusso[ip3], y_4_flusso[ip3+1]-y_4_flusso[ip3],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_5_flusso[ip], y_5_flusso[ip],x_5_flusso[ip+1]-x_5_flusso[ip], y_5_flusso[ip+1]-y_5_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_6_flusso[ip], y_6_flusso[ip],x_6_flusso[ip+1]-x_6_flusso[ip], y_6_flusso[ip+1]-y_6_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_7_flusso[ip3], y_7_flusso[ip3],x_7_flusso[ip3+1]-x_7_flusso[ip3], y_7_flusso[ip3+1]-y_7_flusso[ip3],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_8_flusso[ip], y_8_flusso[ip],x_8_flusso[ip+1]-x_8_flusso[ip], y_8_flusso[ip+1]-y_8_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_9_flusso[ip], y_9_flusso[ip],x_9_flusso[ip+1]-x_9_flusso[ip], y_9_flusso[ip+1]-y_9_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_9_flusso[ip4], y_9_flusso[ip4],x_9_flusso[ip4+1]-x_9_flusso[ip4], y_9_flusso[ip4+1]-y_9_flusso[ip4],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_10_flusso[ip3], y_10_flusso[ip3],x_10_flusso[ip3+1]-x_10_flusso[ip3], y_10_flusso[ip3+1]-y_10_flusso[ip3],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_11_flusso[ip3], y_11_flusso[ip3],x_11_flusso[ip3+1]-x_11_flusso[ip3], y_11_flusso[ip3+1]-y_11_flusso[ip3],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_12_flusso[ip3], y_12_flusso[ip3],x_12_flusso[ip3+1]-x_12_flusso[ip3], y_12_flusso[ip3+1]-y_12_flusso[ip3],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_13_flusso[ip], y_13_flusso[ip],x_13_flusso[ip+1]-x_13_flusso[ip], y_13_flusso[ip+1]-y_13_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_14_flusso[ip], y_14_flusso[ip],x_14_flusso[ip+1]-x_14_flusso[ip], y_14_flusso[ip+1]-y_14_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_15_flusso[ip], y_15_flusso[ip],x_15_flusso[ip+1]-x_15_flusso[ip], y_15_flusso[ip+1]-y_15_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_15_flusso[ip4], y_15_flusso[ip4],x_15_flusso[ip4+1]-x_15_flusso[ip4], y_15_flusso[ip4+1]-y_15_flusso[ip4],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_16_flusso[ip], y_16_flusso[ip],x_16_flusso[ip+1]-x_16_flusso[ip], y_16_flusso[ip+1]-y_16_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_17_flusso[ip], y_17_flusso[ip],x_17_flusso[ip+1]-x_17_flusso[ip], y_17_flusso[ip+1]-y_17_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_18_flusso[ip1], y_18_flusso[ip1],x_18_flusso[ip1+1]-x_18_flusso[ip1], y_18_flusso[ip1+1]-y_18_flusso[ip1],shape='full', lw=0, length_includes_head=True, head_width=.02)
plt.arrow(x_19_flusso[ip], y_19_flusso[ip],x_19_flusso[ip+1]-x_19_flusso[ip], y_19_flusso[ip+1]-y_19_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_20_flusso[ip2], y_20_flusso[ip2],x_20_flusso[ip2+1]-x_20_flusso[ip2], y_20_flusso[ip2+1]-y_20_flusso[ip2],shape='full', lw=0, length_includes_head=True, head_width=.01)

plt.arrow(0.5, 0, -0.01, 0,shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(1.1, 0, -0.01, 0,shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(0, 0.9, 0, -0.01,shape='full', lw=0, length_includes_head=True, head_width=.05)
plt.arrow(0, 0.4, 0, -0.01,shape='full', lw=0, length_includes_head=True, head_width=.05)

#plt.legend(loc='best', prop={'size': 7})
#plt.title('spazio delle fasi Goodwin con Funzione Investimenti',fontdict = font)
plt.xlabel('quota salari $\omega$')
plt.ylabel('occupazione $\lambda$')
#plt.axis([0, 1.4, 0, 2.1]) # X=1.4 buono # PER FARE ANCHE PARTE SOPRA CON OCCUPATI>POPOLAZIONE cioì lambda>1
#plt.axis([0, 1.2, .9, 1.1]) # zoom sulla parte incasinata
plt.axis([0, 1.4, 0, 1.0]) # PER FARE GRAFICO SOLO DELLA PARTE SENSATA cioè con lambda NON superiore 1
plt.grid()
#plt.savefig('relazione corretta Goodwin Funzione Investimenti spazio delle fasi.jpg', dpi=850, transparent=False)
plt.show()
#%%
#%% # PROVA PER CAPIRE DOVE VANNO LE SINGOLE ORBITE. POI SI POTREBBE ANCHE ELIMINARE  ???
x0_19=0.72774
y0_19=0.01
#x0_19=0.8138
#y0_19=1.01
# la sella è a x0_19=0.608868755222622, y0_19=0.965881121906885
x_19_flusso, y_19_flusso = flusso_traiettoria_singola(x0_19,y0_19,dt,numero_di_punti_tempo_voluti,sistema_2_eq_mio_con_funzione_investimenti)
plt.plot(x_19_flusso, y_19_flusso,'b')
plt.axis([0, 1.4, 0, 2.1]) # X=1.4 buono
#plt.axis([0, 1.2, .9, 1.1]) # zoom sulla parte incasinata
plt.grid()
plt.show
#%%
# ANDAMENTO TEMPORALE PARAGONE DIVERSE DINAMICHE COME PARE VOLERE TURCHETTI
# Forse aggiungere legenda col punto iniziale o informazioni ad esempio "dinamica sella"

# grafico temporale con una bella dinamica ma forse non è quello che vuole Turchetti, oppure ne vuole tante a confronto 
plt.plot(t, x_5_flusso,'r',label='quota salari $\omega(t)$')
plt.plot(t, y_5_flusso,'black',label='occupazione $\lambda(t)$')
plt.legend(loc='lower right', prop={'size': 6})
#plt.title('andamento temporale, quello scelto da STEVE KEEN',fontdict = font)
plt.xlabel('tempo t')
plt.ylabel('$\omega, \lambda$')
plt.grid()
plt.xlim(0,60)
#plt.savefig('relazione corretta Goodwin Funzione Investimenti andamento temporale 1.jpg', dpi=650, transparent=False)
plt.show()

plt.plot(t, x_3_flusso,'r',label='quota salari $\omega(t)$')
plt.plot(t, y_3_flusso,'k',label='occupazione $\lambda(t)$')
plt.legend(loc='best',prop={'size': 6})
#plt.legend(loc='lower right', prop={'size': 6})
#plt.title('andamento temporale, va a zero passando sotto la sella',fontdict = font)
plt.xlabel('tempo t')
plt.ylabel('$\omega, \lambda$')
plt.grid()
plt.xlim(0,60)
#plt.savefig('relazione corretta Goodwin Funzione Investimenti andamento temporale 2.jpg', dpi=650, transparent=False)
plt.show()

plt.plot(t, x_2_flusso,'r',label='quota salari $\omega(t)$')
plt.plot(t, y_2_flusso,'k',label='occupazione $\lambda(t)$')
plt.legend(loc='best',prop={'size': 6})
#plt.title('andamento temporale, sulla sella va all'infinito omega',fontdict = font)
plt.xlabel('tempo t')
plt.ylabel('$\omega, \lambda$')
plt.grid()
plt.xlim(0,60)
#plt.ylim(0,100)
#plt.savefig('relazione corretta Goodwin Funzione Investimenti andamento temporale 3.jpg', dpi=650, transparent=False)
plt.show()

plt.plot(t, x_8_flusso,'r',label='quota salari $\omega(t)$')
plt.plot(t, y_8_flusso,'k',label='occupazione $\lambda(t)$')
plt.legend(loc='best',prop={'size': 6})
#plt.title('andamento temporale',fontdict = font)
plt.xlabel('tempo t')
plt.ylabel('$\omega, \lambda$')
plt.grid()
plt.xlim(0,60)
#plt.savefig('relazione corretta Goodwin Funzione Investimenti andamento temporale 4.jpg', dpi=650, transparent=False)
plt.show()
# IN TOTALE QUINDI NELLA RELAZIONE NE METTI 4 DI GRAFICI TEMPORALI

plt.plot(t, x_9_flusso,'r',label='quota salari $\omega(t)$')
plt.plot(t, y_9_flusso,'k',label='occupazione $\lambda(t)$')
plt.legend(loc='best',prop={'size': 6})
#plt.title('andamento temporale, è quello maledetto che fa 2 giri non usare',fontdict = font)
plt.xlabel('tempo t')
plt.ylabel('$\omega, \lambda$')
plt.grid()
plt.xlim(0,60)
#plt.savefig('relazione corretta Goodwin Funzione Investimenti andamento temporale 5.jpg', dpi=650, transparent=False)
plt.show()
#%%
#%%
# DISEGNO NULLCLINE e GRAFICO GLI ZERI DELLE FUNZIONI 

# nullcline e zeri della funzione f(omega) che compare in dlambda/dt come vuole Turchetti
assex=np.linspace(-1,2,1000)

#    domegadt = omega * (functionalform(A,B,C,D,lanbda) - alfa)
#    dlanbdadt = lanbda * (((functionalform(E,F,G,H,((1-omega)/ni)))/ni) - alfa - beta - gamma)

#omegaeq = 1 - ((ni*F)/G) +o- (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) )   # è un centro
#lanbdaeq = (B/C) -o+ ( np.sqrt(A/(alfa + D)) )/C   # è un centro, per cui functionalform(A,B,C,D,lanbda0)=alfa anche se dopo tanti zeri differisce per cui potrebbe non andare bene

# Grafico omega - f(omega)
plt.plot(assex, ((functionalform(E,F,G,H,((1-assex)/ni)))/ni) - alfa - beta - gamma,'blue')
plt.plot([0,0],[-10,20],'k')
plt.plot([-10,20],[0,0],'k')
plt.plot([1-(ni*F)/G,1-(ni*F)/G],[-10,20],'k')
plt.plot(1 - ((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),0,'ob')
plt.plot(1 - ((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),0,'ob')
plt.legend(loc='best',prop={'size': 6})
#plt.title('grafico omega f(omega) NON è la nullcline di $\dot{\lambda}$',fontdict = font)
plt.xlabel('$\omega$')
plt.ylabel('$f(\omega)$')
plt.grid()
plt.xlim(-0.3,1.2)
plt.ylim(-0.4,3.2)
#plt.savefig('relazione corretta Goodwin Funzione Investimenti grafico omega f(omega).jpg', dpi=650, transparent=False)
plt.show()

# Grafico lambda - g(lambda)
plt.plot(assex, functionalform(A,B,C,D,assex) - alfa,'red')
plt.plot([0,0],[-1,3],'k')
plt.plot([-1,3],[0,0],'k')
plt.plot([1,1],[-1,3],'k')
plt.plot(lanbda0eq1,0,'or')
plt.legend(loc='best',prop={'size': 6})
#plt.title('grafico lambda g(lambda) NON è la nullcline di $\dot{\omega}$',fontdict = font)
plt.xlabel('$\lambda$')
plt.ylabel('$g(\lambda)$')
plt.grid()
plt.xlim(-0.1,1.2)
plt.ylim(-0.2,1.2)
#plt.savefig('relazione corretta Goodwin Funzione Investimenti grafico lambda g(lambda).jpg', dpi=650, transparent=False)
plt.show()

# Grafico lambda - w[lambda]
plt.plot(assex, functionalform(A,B,C,D,assex),'blue')
plt.plot([0,0],[-1,3],'k')
plt.plot([-1,3],[0,0],'k')
plt.plot([1,1],[-1,3],'k')
plt.plot(1-(A/D)**0.5,0,'ob')
#plt.plot(1+(A/D)**0.5,0,'ob')
plt.legend(loc='best',prop={'size': 6})
#plt.title('funzione risposta lavoratori grafico lambda w[lambda]',fontdict = font)
plt.xlabel('$\lambda$')
plt.ylabel('$w[\lambda]$')
plt.grid()
plt.xlim(-0.1,1.2)
plt.ylim(-0.2,1.2)
#plt.savefig('relazione corretta Goodwin Funzione Investimenti grafico lambda w[lambda].jpg', dpi=650, transparent=False)
plt.show()

# Grafico omega - k[(1-omega)/ni]
plt.plot(assex, functionalform(E,F,G,H,((1-assex)/ni)),'red')
plt.plot([0,0],[-1,3],'k')
plt.plot([-1,3],[0,0],'k')
plt.plot([1-(ni*F)/G,1-(ni*F)/G],[-10,20],'k')
plt.plot(1 - ((ni*F)/G) + (ni/G)*( np.sqrt(E/H) ),0,'or')
plt.plot(1 - ((ni*F)/G) - (ni/G)*( np.sqrt(E/H) ),0,'or')
plt.legend(loc='best',prop={'size': 6})
#plt.title('funzione risposta capitalisti grafico omega k[(1-omega)/ni]',fontdict = font)
plt.xlabel('$\omega$')
plt.ylabel('$k[(1-\omega)/ \\nu]$')
plt.grid()
plt.xlim(-0.1,1.2)
plt.ylim(-0.2,1.2)
#plt.savefig('relazione corretta Goodwin Funzione Investimenti grafico omega k[1-omega/nu].jpg', dpi=650, transparent=False)
plt.show()

# Faccio vedere a scopo dimostrativo che il campo vettoriale è perpendicolare alle nullcline vicino ad esse
# Draw Nullclines and Quiver plot
# plot nullclines
plt.plot([-5,5],[(B/C) - ( np.sqrt(A/(alfa + D)) )/C ,(B/C) - ( np.sqrt(A/(alfa + D)) )/C ], 'b-', lw=2, label='$\dot{\omega}$ nullcline')
plt.plot([-5,5],[(B/C) + ( np.sqrt(A/(alfa + D)) )/C ,(B/C) + ( np.sqrt(A/(alfa + D)) )/C ], 'b-', lw=2)
plt.plot([1 - ((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),1 - ((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) )],[-5,5], 'r-', lw=2, label='$\dot{\lambda}$ nullcline')
plt.plot([1 - ((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),1 - ((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) )],[-5,5], 'r-', lw=2)
plt.plot(functionalform(A,B,C,D,assex) - alfa,assex,'green') # questa in realtà andrebbe tolta perchè N.B.: g(lambda) NON è la nullcline, essa sono i due valori di lambda* infatti la nullcline di lambda non dipende da omega, è un retta orizzontale quindi
plt.plot(assex, ((functionalform(E,F,G,H,((1-assex)/ni)))/ni) - alfa - beta - gamma,'green') #  questa in realtà andrebbe tolta perchè N.B.: f(omega) NON è la nullcline, essa è il valore di omega* infatti la nullcline di omega non dipende da lambda, è un retta verticale quindi
plt.plot([0,0],[-5,5],'k',lw=2)
plt.plot([-5,5],[0,0],'k',lw=2)
plt.plot([-5,5],[1,1],'--k')
plt.plot([1-(ni*F)/G,1-(ni*F)/G],[-5,5],'--k')
plt.legend(loc='best',prop={'size': 9})
plt.xlabel('$\omega$')
plt.ylabel('$\lambda$')
plt.xlim(-0.1,1.2)
plt.ylim(-0.1,1.2)
# plot fixed points
for point in equilibria:
    plt.plot(point[0],point[1],"blue", marker = "o", markersize = 6.0)
plt.title("Quiverplot with nullclines",fontdict = font)
plt.quiver(X1, Y1, DX1, DY1, color='deepskyblue',pivot='mid')
plt.grid()
# carina da tenere che si vede meglio anche se NON LA DEVI METTERE NELLA RELAZIONE
#plt.savefig('relazione corretta Goodwin Funzione Investimenti nullcline spazio delle fasi.jpg', dpi=650, transparent=False)
plt.show()

# Grafico omega,lambda - f(omega),g(lambda)
plt.plot(assex,((functionalform(E,F,G,H,((1-assex)/ni)))/ni) - alfa - beta - gamma,'blue',label='$f(\omega)$')
plt.plot(assex,functionalform(A,B,C,D,assex)-alfa,'red', label='$g(\lambda)$')
plt.plot([-1,2],[0,0],'k', [0,0],[-1,2],'k',[1 - ((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),1 - ((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),lanbda0eq1],[0,0,0],'or')
# oppure tutto insieme ma non funziona il comando label non so perchè!
#plt.plot(assex,functionalform(A,B,C,D,assex)-alfa,'red', assex,(((1-assex)/ni) - alfa - beta - gamma),'blue', [-1,2],[0,0],'k', [0,0],[-1,2],'k',[omega0eq1,lanbda0eq1],[0,0],'or')
plt.legend(loc='upper center',prop={'size': 8})
plt.plot([1-(ni*F)/G,1-(ni*F)/G],[-5,5],'--k')
plt.plot([1,1],[-5,5],'--k')
plt.axis([-0.2, 1.2, -0.2, 1.3])
plt.xlabel('$\omega, \lambda$')
plt.ylabel('$f(\omega), g(\lambda)$')
#plt.xlim(-1,2)
#plt.ylim(-1,1)
plt.grid()
plt.show()
#%%
#%%
# CALCOLO AUTOVALORI DEI PUNTI FISSI PER STUDIARE STABILITA'. N.B.: NON FIDARSI DI QUELLO CHE DICE CHE SONO. Inoltre devi inserire manualmente i punti fissi trovati sopra e la matrice Jacobiana generica
# non funziona molto bene e per i punti "centri" non funziona proprio! Siamo messi bene...

def eigenvalues(x,y):
#    domegadt = omega * (functionalform(A,B,C,D,lanbda) - alfa)
    a11 = functionalform(A,B,C,D,y)-alfa         # differentiated with respect to x
    a12 = x* ( (2*A*C)/(B-C*y)**3 )              # differentiated with respect to y
    
#    dlanbdadt = lanbda * (((k[pi/ni])/ni) - alfa - beta - gamma)
    a21 = (-(y/ni)) * ( ((2*E*G)/ni)/(F-G*((1-x)/ni))**3 )                # differentiated with respect to x
    a22 = (functionalform(E,F,G,H,((1-x)/ni))/ni) - alfa - beta - gamma   # differentiated with respect to y

    tr = a11 + a22
    det = a11*a22 - a12*a21
    lambda1 = (tr - sqrt(tr**2 - 4*det))/2
    lambda2 = (tr + sqrt(tr**2 - 4*det))/2
    print('Check the fixed point  %s, %s' % (x,y)) 
    print('The real part of the first eigenvalue is %s' %lambda1.real)
    print('The real part of the second eigenvalue is %s' % lambda2.real)    
    
    if (lambda1.real < 0 and lambda2.real < 0):
        print('The fixed point in %s, %s is a sink. It is stable' % (x,y))
    if (lambda1.real > 0 and lambda2.real > 0):
        print('The fixed point in %s, %s is a source. It is unstable' % (x,y))
    if (lambda1.real > 0 and lambda2.real < 0):
        print('The fixed point in %s, %s is a saddle. It is unstable' % (x,y))
    if (lambda1.real < 0 and lambda2.real > 0):
        print('The fixed point in %s, %s is a saddle. It is unstable' % (x,y))
#    if (lambda1.real < 0 and lambda2.real > 0):
#        print('The fixed point in %s, %s is unstable' % (x,y))
    print('----------------------------')
    return lambda1 , lambda2

equilibria_manuale=np.array([[0,0], [1 - ((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),(B/C) - ( np.sqrt(A/(alfa + D)) )/C], [1 - ((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),(B/C) - ( np.sqrt(A/(alfa + D)) )/C], [1 - ((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),(B/C) + ( np.sqrt(A/(alfa + D)) )/C], [1 - ((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),(B/C) + ( np.sqrt(A/(alfa + D)) )/C]])
# perchè con equilibria non funziona siccome ci sono tutti i decimali, dice che non può fare sqrt() di un float...non so perchè! 
# iterate through list of fixed points
for x,y in equilibria_manuale:
    eigenvalues(x,y)    

#for item in equilibria:
#    eigenvalues(item[0],item[1])
#%%
# AUTOVALORI JACOBIANA SISTEMA 2 EQUAZIONI CON FUNZIONE INVESTIMENTI k[pi/ni]
def jacobiana_sistema_2_eq_goodwin_con_funzione_investimenti(omegaeq,lanbdaeq):
    primo = functionalform(A,B,C,D,lanbdaeq) - alfa
    secondo = omegaeq * ( (2*A*C)/(B-C*lanbdaeq)**3 )
    terzo = (-(lanbdaeq/ni)) * ( ((2*E*G)/ni)/(F-G*((1-omegaeq)/ni))**3 )
    quarto = (functionalform(E,F,G,H,((1-omegaeq)/ni))/ni) - alfa - beta - gamma
    primo=round(primo,9) # sennò non viene zero perchè tiene 17 cifre decimali!
    secondo=round(secondo,9)
    terzo=round(terzo,9)
    quarto=round(quarto,9)
    #quarto=round(quarto) # perchè sennò viene -0.0 quando omegaeq è quello che annulla "quarto". N.B.: TTENZIONE CHE IN TUTTI GLI ALTRI CASI TI FA VENIRE RISULTATI SBAGLIATI!!! Meglio lasciare sempre commentato, non usare
    return np.array([[primo, secondo],[terzo, quarto]])


#omegaeq=0  # è un punto fisso stabile
#lanbdaeq=0

#omegaeq = 1 - ((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ) # è un centro
#lanbdaeq = (B/C) - ( np.sqrt(A/(alfa + D)) )/C # è un centro, per cui functionalform(A,B,C,D,lanbda0)=alfa anche se dopo tanti zeri differisce per cui potrebbe non andare bene
# ovvero questi sono:
#omegaeq = 0.8611312447773782
#lanbdaeq = 0.965881121906885

omegaeq = 1 - ((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ) # è un punto sella quello che non trovavi!
lanbdaeq = (B/C) - ( np.sqrt(A/(alfa + D)) )/C # è una sella, per cui functionalform(A,B,C,D,lanbda0)=alfa anche se dopo tanti zeri differisce per cui potrebbe non andare bene
# ovvero questi sono:
#omegaeq = 0.608868755223  ?forse servono altre cifre?
#lanbdaeq = 0.965881121906885

#omegaeq = 1 - ((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ) # è un punto sella
#lanbdaeq = (B/C) + ( np.sqrt(A/(alfa + D)) )/C # è una sella, per cui functionalform(A,B,C,D,lanbda0)=alfa anche se dopo tanti zeri differisce per cui potrebbe non andare bene
# ovvero questi sono:
#omegaeq = 0.8611312447773782
#lanbdaeq = 1.03411887809 >1  ?forse servono altre cifre?

#omegaeq = 1 - ((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ) # è un centro
#lanbdaeq = (B/C) + ( np.sqrt(A/(alfa + D)) )/C # è un centro, per cui functionalform(A,B,C,D,lanbda0)=alfa anche se dopo tanti zeri differisce per cui potrebbe non andare bene
# ovvero questi sono:
#omegaeq = 0.608868755223  ?forse servono altre cifre?
#lanbdaeq = 1.03411887809 >1  ?forse servono altre cifre?


jacobiana_sistema_goodwin_con_funzione_investimenti = jacobiana_sistema_2_eq_goodwin_con_funzione_investimenti(omegaeq,lanbdaeq)
traccia_sistema_goodwin_con_funzione_investimenti = jacobiana_sistema_goodwin_con_funzione_investimenti[0,0] + jacobiana_sistema_goodwin_con_funzione_investimenti[1,1]
determinante_sistema_goodwin_con_funzione_investimenti = jacobiana_sistema_goodwin_con_funzione_investimenti[0,0]*jacobiana_sistema_goodwin_con_funzione_investimenti[1,1] - jacobiana_sistema_goodwin_con_funzione_investimenti[0,1]*jacobiana_sistema_goodwin_con_funzione_investimenti[1,0]

autovalori_sistema_goodwin_con_funzione_investimenti, autovettori_sistema_goodwin_con_funzione_investimenti = LA.eig(jacobiana_sistema_goodwin_con_funzione_investimenti)

print('\nLa matrice Jacobiana del sistema a 2 equazioni con funzione investimenti k[pi/ni] è:\n',jacobiana_sistema_goodwin_con_funzione_investimenti,'\n')
print('Il determinante det(J) è:\n', determinante_sistema_goodwin_con_funzione_investimenti)
print('La traccia tr(J) è:\n', traccia_sistema_goodwin_con_funzione_investimenti)
print('Gli autovalori della matrice J sono:\n', autovalori_sistema_goodwin_con_funzione_investimenti)
print('Una possibile scelta per i 2 autovettori di J è:\n',autovettori_sistema_goodwin_con_funzione_investimenti)

plt.plot(determinante_sistema_goodwin_con_funzione_investimenti, traccia_sistema_goodwin_con_funzione_investimenti, marker='o', linestyle='--', color='r')
plt.title('Grafico determinante e traccia della matrice Jacobiana \n dalla quale si ricavano gli autovalori per quel punto fisso')
plt.xlabel('det(J)')
plt.ylabel('tr(J)')
plt.xlim(-10,10)
plt.ylim(-10,10)
plt.grid()
plt.show()
#%%
#%%
# STREAMPLOT per disegnare il ritratto di fase ma NON USARE perchè è meglio quiver!!!! 

plt.streamplot(X1, Y1, DX1, DY1, linewidth=1, density=2, arrowstyle='->', arrowsize=1.5)
#plt.streamplot(X1, Y1, DX1, DY1,color=M, linewidth=2, cmap=plt.cm.autumn)
# alla fine puoi aggiungere anche il comando density = 1.4
plt.xlim(0,2)
plt.ylim(0,2)
plt.grid()
plt.show()