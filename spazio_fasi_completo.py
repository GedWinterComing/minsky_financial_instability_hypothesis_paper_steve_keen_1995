# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 13:28:05 2020

@author: Gabriele
"""

#%%
# VERSIONE INGEGNERIZZATA PER FARE QUANTE PIU' COSE AUTOMATICAMENTE ED
# EVITARE IL PROBLEMA DI NON TROVARE DEI PUNTI FISSI COME ACCADUTO NELLA RELAZIONE CON LA SADDLE

#%%
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import sympy as sm
from mpl_toolkits.mplot3d import Axes3D
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
# SISTEMA PRINCIPALE DOVE CI SONO TUTTE LE ISTRUZIONI caso II modello semplificato di Turchetti nella controrelazione

def sistema_II_semplificato(X):
    #spacchettamento del vettore
    x = X[0]
    y = X[1]
    
    dxdt = -x*(1-y)
    dydt = -y*(1-x)*((1/2)-x)
    vettore_derivate = np.array([dxdt,dydt])
    return vettore_derivate 
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
# PUNTI FISSI
# Per trovare automaticamente gli zeri del sistema di equazioni differenziali

# define the system in this way (asuming a predator-prey-system with no negative values)
# to avoid interference x = u (predator) and y = v (prey)
# Qua devi riscrivere le eq del sistema con i nuovi nomi delle variabili per evitare di usare le stesse lettere già usate per non creare confusioni. Visto che queste saranno simboli
# Se ci fosse bisogno di una terza variabile per una terza eq. differenziale usa "w" e "W"
u, v = sm.symbols('u, v', negative=False)
U = -u*(1-v)
#V = -v*(1-u)*((1/2)-u) # non funziona bene, dà un sacco di zeri
#V = -v*(1-u)*((1-u-u)/2) # questo funziona bene, ma non so perchè
V = -v*(1-u)*( sm.Rational(1, 2) - u) # anche questo funziona bene
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
x0_1 = 1.1
y0_1 = 1.1
#  dall'alto a sinistra va all'origine
x0_2 = 0.12
y0_2 = 2
# da in basso a destra verso l'origine
x0_3 = 2     
y0_3 = 0.2
# cerchio grande
x0_4 = 0.7   
y0_4 = 1
# sella braccio stabile va all'origine. In realtà questo è inutile perchè il braccio che scende da in alto a sinistra poi fa questo giro lungo il pesce torna al nodo sella passando sotto e poi continua verso l'origine
x0_5 = 0.499
y0_5 = .999
#  dall'alto a sinistra va all'origine
x0_6 = 0.1
y0_6 = 2   
# parte in alto a sinistra fa il giro lungo i cerchi
x0_7 = 0.15  
y0_7 = 2
# parte in alto a sinistra fa il giro lungo cerchi
x0_8 = 0.2  
y0_8 = 2
# sella braccio instabile, fa il giro e si richiude sotto. In realtà questo è inutile perchè il braccio che scende da in alto a sinistra poi fa questo giro lungo il pesce torna al nodo sella passando sotto e poi continua verso l'origine
x0_9 = 0.51
y0_9 = 1.001
# braccio stabile della sella che parte in alto a sinistra, è stato difficilissimo da trovare! N.B.: A 0.11 e 2 fa ancora il giro lungo dal centro invece che andare al bracio
x0_10 = 0.10514  # x0_10=0.105 potrebbe quasi andare bene, se invece metti 0.1051 non fa tutto il gro del pesce ma dall'alto va verso il basso tipo una cuspide, io penso che il vero punto sella dovrebbe fare il giro del pesce: dall'alto arriva al punto va a destra fa il pesce, arriva sotto al pesce e poi va all'origine
y0_10 = 2   
# alto sinistra verso il basso che va all'origine
x0_11 = 0.08 
y0_11 = 2
# esterno sotto
x0_12 = 2  
y0_12 = 0.5

#dt = 0.01
#numero_di_punti_tempo_voluti = 5000  # deve essere un numero tale che tempo_finale venga un numero intero
dt=0.01
numero_di_punti_tempo_voluti = 5000
tempo_finale = numero_di_punti_tempo_voluti * dt
t = np.linspace(0,tempo_finale,numero_di_punti_tempo_voluti) # per fare il grafico dopo, l'ultimo è proprio tempo_finale

# ovviamente tutto ciò si poteva saltare mettendo direttamente i numeri al posto di x0_i e y0_i definendoli qui e non sopra che è ridondante
x_1_flusso, y_1_flusso = flusso_traiettoria_singola(x0_1,y0_1,dt,numero_di_punti_tempo_voluti,sistema_II_semplificato)
x_2_flusso, y_2_flusso = flusso_traiettoria_singola(x0_2,y0_2,dt,numero_di_punti_tempo_voluti,sistema_II_semplificato)
x_3_flusso, y_3_flusso = flusso_traiettoria_singola(x0_3,y0_3,dt,numero_di_punti_tempo_voluti,sistema_II_semplificato)
x_4_flusso, y_4_flusso = flusso_traiettoria_singola(x0_4,y0_4,dt,numero_di_punti_tempo_voluti,sistema_II_semplificato)
x_5_flusso, y_5_flusso = flusso_traiettoria_singola(x0_5,y0_5,dt,numero_di_punti_tempo_voluti,sistema_II_semplificato)
x_6_flusso, y_6_flusso = flusso_traiettoria_singola(x0_6,y0_6,dt,numero_di_punti_tempo_voluti,sistema_II_semplificato)
x_7_flusso, y_7_flusso = flusso_traiettoria_singola(x0_7,y0_7,dt,numero_di_punti_tempo_voluti,sistema_II_semplificato)
x_8_flusso, y_8_flusso = flusso_traiettoria_singola(x0_8,y0_8,dt,numero_di_punti_tempo_voluti,sistema_II_semplificato)
x_9_flusso, y_9_flusso = flusso_traiettoria_singola(x0_9,y0_9,dt,numero_di_punti_tempo_voluti,sistema_II_semplificato)
x_10_flusso, y_10_flusso = flusso_traiettoria_singola(x0_10,y0_10,dt,numero_di_punti_tempo_voluti,sistema_II_semplificato)
x_11_flusso, y_11_flusso = flusso_traiettoria_singola(x0_11,y0_11,dt,numero_di_punti_tempo_voluti,sistema_II_semplificato)
x_12_flusso, y_12_flusso = flusso_traiettoria_singola(x0_12,y0_12,dt,numero_di_punti_tempo_voluti,sistema_II_semplificato)
#%%
# GRAFICO
font = {'family':'serif','color':'darkred','weight':'normal','size': 10.7,}
# plt.title('Damped exponential decay', fontdict=font)
# plt.text(2, 0.65, r'$\cos(2 \pi t) \exp(-t)$', fontdict=font)
# plt.xlabel('time (s)', fontdict=font)
# plt.ylabel('voltage (mV)', fontdict=font)


# grafico temporale con una bella dinamica ma forse non è quello che vuole Turchetti, oppure ne vuole tante a confronto 
plt.plot(t, x_1_flusso,'r', label='predatore $x(t)$')
plt.plot(t, y_1_flusso,'black', label='preda $y(t)$')
plt.legend(loc='lower right', prop={'size': 6})
plt.title('andamento temporale',fontdict = font)
plt.xlabel('tempo t')
plt.ylabel('$x, y$')
plt.grid()
#plt.savefig('Figure 2. Wage share and employment, basic Goodwin model.jpg', dpi=850, transparent=False)
plt.show()


# FRECCE DEL CAMPO VETTORIALE

# creo la griglia, in realtà saranno 2 matrici una con coordinata x e l'altra y così prese elemento per elemento danno il punto (x_i,y_i)
x_griglia = np.linspace(0, 2, 14)
y_griglia = np.linspace(0, 2, 14)

# il comando sotto è interessante perchè assegna automaticamente i termini in ordine: es. gaia=np.array([3,4]) poi con hippie1,hippie2=gaia si ha hippie1=gaia[0]=3 e hippie2=gaia[1]=4
X1 , Y1  = np.meshgrid(x_griglia, y_griglia)    # create a grid
# ritorna una prima matrice che ha come righe il vettore x ripetuto tante volte quante la lunghezza del vettore y, cioè len(y) volte, e una seconda matrice che ha il vettore y come colonne ripetuto tante volte quanta la lunghezza del vettore x, len(x)
DX1, DY1 = sistema_II_semplificato([X1, Y1])    # compute growth rate on the grid
# la moltiplicazione * usata in "sistema_II_semplificato" è quella normale e se X come in questo caso è composto da 2 matrici, quando spacchetta X assegna la prima matrice (cioè le coordinate x su tutta la griglia) a x e la seconda (cioè la griglia per le coordinate y) a y, quando poi si trova a fare x*y moltiplica elemento per elemento entrambe le matrici senza fare il prodotto matriciale! Perciò le due matrici devono avere entrambe la stessa dimensione, tutte e due devono essere nxn
M = (np.hypot(DX1, DY1))                        # norm growth rate 
M[ M == 0] = 1.                                 # avoid zero division errors 
DX1 /= M                                        # normalize each arrows
DY1 /= M
# se fai M.shape dà (20,20) è infatti una matrice 20x20 siccome hai usato due volte x = np.linspace(0, 2, 20) o x e y, ma sono uguali x=y=np.linspace(0, 2, 20)

#plt.quiver(X1, Y1, DX1, DY1, color='deepskyblue',pivot='mid')
#plt.quiver(X1, Y1, DX1, DY1, M, pivot='mid') # M gli dà una scala di colori che sono i numeri


plt.plot(x_1_flusso, y_1_flusso,'b')
plt.plot(x_2_flusso, y_2_flusso,'b')
plt.plot(x_3_flusso, y_3_flusso,'b')
plt.plot(x_4_flusso, y_4_flusso,'b')
#plt.plot(x_5_flusso, y_5_flusso,'b') # inutile
plt.plot(x_6_flusso, y_6_flusso,'b')
plt.plot(x_7_flusso, y_7_flusso,'b')
plt.plot(x_8_flusso, y_8_flusso,'b')
#plt.plot(x_9_flusso, y_9_flusso,'b') # inutile
plt.plot(x_10_flusso, y_10_flusso,'b')
plt.plot(x_11_flusso, y_11_flusso,'b')
plt.plot(x_12_flusso, y_12_flusso,'b')

# Nel grafico che sto creando metto anche tutti punti fissi. Ci sono 2 modi possibili per farlo
# versione "automatica" però non puoi variare i colori:
#for point in equilibria:
#    plt.plot(point[0],point[1],"royalblue", marker = "o", markersize = 7.0)

# versione "manuale":
plt.plot(0,0,'ob')
plt.plot(0.5,1,'o',color='blue')
plt.plot(1,1,'o',color='royalblue')   # condizione iniziale o meglio punto del centro dei cicli. Basta allontanarsi poco da questo punto con la condizione iniziale che si formano i cycle!
#plt.plot(omega0eq1,lanbda0eq1,'og', label='Basic Goodwin Equilibrium') # condizione iniziale o meglio punto del centro dei cicli. Basta allontanarsi poco da questo punto con la condizione iniziale che si formano i cycle!

ip=500
ip1=100
ip2=300
ip3=1820
ip4=2900
ip5=4250
plt.arrow(x_1_flusso[ip], y_1_flusso[ip],x_1_flusso[ip+1]-x_1_flusso[ip], y_1_flusso[ip+1]-y_1_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_2_flusso[ip], y_2_flusso[ip],x_2_flusso[ip+1]-x_2_flusso[ip], y_2_flusso[ip+1]-y_2_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_3_flusso[ip1], y_3_flusso[ip1],x_3_flusso[ip1+1]-x_3_flusso[ip1], y_3_flusso[ip1+1]-y_3_flusso[ip1],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_4_flusso[ip], y_4_flusso[ip],x_4_flusso[ip+1]-x_4_flusso[ip], y_4_flusso[ip+1]-y_4_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
#plt.arrow(x_5_flusso[ip], y_5_flusso[ip],x_5_flusso[ip+1]-x_5_flusso[ip], y_5_flusso[ip+1]-y_5_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_6_flusso[ip], y_6_flusso[ip],x_6_flusso[ip+1]-x_6_flusso[ip], y_6_flusso[ip+1]-y_6_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_7_flusso[ip], y_7_flusso[ip],x_7_flusso[ip+1]-x_7_flusso[ip], y_7_flusso[ip+1]-y_7_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_8_flusso[ip], y_8_flusso[ip],x_8_flusso[ip+1]-x_8_flusso[ip], y_8_flusso[ip+1]-y_8_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
#plt.arrow(x_9_flusso[ip], y_9_flusso[ip],x_9_flusso[ip+1]-x_9_flusso[ip], y_9_flusso[ip+1]-y_9_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_10_flusso[ip], y_10_flusso[ip],x_10_flusso[ip+1]-x_10_flusso[ip], y_10_flusso[ip+1]-y_10_flusso[ip],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_10_flusso[ip3], y_10_flusso[ip3],x_10_flusso[ip3+1]-x_10_flusso[ip3], y_10_flusso[ip3+1]-y_10_flusso[ip3],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_10_flusso[ip4], y_10_flusso[ip4],x_10_flusso[ip4+1]-x_10_flusso[ip4], y_10_flusso[ip4+1]-y_10_flusso[ip4],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_10_flusso[ip5], y_10_flusso[ip5],x_10_flusso[ip5+1]-x_10_flusso[ip5], y_10_flusso[ip5+1]-y_10_flusso[ip5],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_11_flusso[ip2], y_11_flusso[ip2],x_11_flusso[ip2+1]-x_11_flusso[ip2], y_11_flusso[ip2+1]-y_11_flusso[ip2],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(x_12_flusso[ip1], y_12_flusso[ip1],x_12_flusso[ip1+1]-x_12_flusso[ip1], y_12_flusso[ip1+1]-y_12_flusso[ip1],shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(0.3, 0, -0.01, 0,shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(0.8, 0, -0.01, 0,shape='full', lw=0, length_includes_head=True, head_width=.04)
plt.arrow(0, 0.8, 0, -0.01,shape='full', lw=0, length_includes_head=True, head_width=.05)
plt.arrow(0, 1.4, 0, -0.01,shape='full', lw=0, length_includes_head=True, head_width=.05)

plt.legend(loc='best', prop={'size': 7})
plt.title('spazio delle fasi II sistema',fontdict = font)
plt.xlabel('predatore x')
plt.ylabel('preda y')
plt.axis([0, 2, 0, 2])
plt.grid()
#plt.savefig('Figure 3. Cyclical and equilibrium time paths with Basic Goodwin Equilibrium.jpg', dpi=550, transparent=False)
plt.show()
#%%
# ANDAMENTO TEMPORALE PARAGONE DIVERSE DINAMICHE COME PARE VOLERE TURCHETTI
# Forse aggiungere legenda col punto iniziale o informazioni ad esempio "dinamica sella"

plt.plot(t, x_4_flusso,'r',t, x_2_flusso,'--r',label='predatore $x(t)$')
plt.plot(t, y_4_flusso,'k',t, y_2_flusso,'--k',label='preda $y(t)$')
plt.legend(loc='best',prop={'size': 7})
#plt.plot(t, x_1_flusso,'r', label='predatore $x(t)$')
#plt.plot(t, y_1_flusso,'black', label='preda $y(t)$')
#plt.legend(loc='lower right', prop={'size': 6})
plt.title('andamento temporale',fontdict = font)
plt.xlabel('tempo t')
plt.ylabel('$x, y$')
plt.grid()
#plt.savefig('Figure 2. Wage share and employment, basic Goodwin model.jpg', dpi=850, transparent=False)
plt.show()


plt.plot(t, x_9_flusso,'r',t, x_10_flusso,'--r',label='predatore $x(t)$')
plt.plot(t, y_9_flusso,'k',t, y_10_flusso,'--k',label='preda $y(t)$')
plt.legend(loc='best',prop={'size': 7})
plt.title('andamento temporale',fontdict = font)
plt.xlabel('tempo t')
plt.ylabel('$x, y$')
plt.grid()
#plt.savefig('Figure 2. Wage share and employment, basic Goodwin model.jpg', dpi=850, transparent=False)
plt.show()
#%%
#%%
# DISEGNO NULLCLINE e GRAFICO GLI ZERI DELLE FUNZIONI 

# nullcline e zeri della funzione f(x) che compare in dy/dt come vuole Turchetti
assex=np.linspace(-2,2,100)
plt.plot(assex,-(1-assex)*((1/2)-assex),'blue',[0.5,.5],[-1,2],'r',[-.5,2],[1,1],'r',[1,1],[-1,2],'r',[-1,2],[0,0],'k',[0.5,1],[0,0],'or')
plt.axis([-0.5, 2, -1, 1.3])
#plt.xlim(-1,2)
#plt.ylim(-1,1)
plt.grid()
plt.show()

# Faccio vedere a scopo dimostrativo che il campo vettoriale è perpendicolare alle nullcline vicino ad esse
# Draw Nullclines and Quiver plot
# plot nullclines
plt.plot([0,2],[1,1], 'r-', lw=2, label='x-nullcline')
plt.plot([1,1],[0,2], 'b-', lw=2, label='y-nullcline')
plt.plot([0.5,0.5],[0,2], 'b-', lw=2,)

# plot fixed points
for point in equilibria:
    plt.plot(point[0],point[1],"red", marker = "o", markersize = 10.0)
plt.title("Quiverplot with nullclines")
plt.legend(loc='best')
plt.quiver(X1, Y1, DX1, DY1, color='deepskyblue',pivot='mid')
plt.grid()
plt.show()

#%%
#%%
# CALCOLO AUTOVALORI DEI PUNTI FISSI PER STUDIARE STABILITA'. N.B.: NON FIDARSI DI QUELLO CHE DICE CHE SONO. Inoltre devi inserire manualmente i punti fissi trovati sopra e la matrice Jacobiana generica
# non funziona molto bene e l'altro modo invece non funziona proprio! Siamo messi bene...

from scipy import sqrt

def eigenvalues(x,y):
    # -x*(1 - y)
    a11 = -(1 - y)                        # differentiated with respect to x
    a12 = x                               # differentiated with respect to y
    # -y*(1 - x)*(1/2 - x)
    a21 = -y*(-x**2 + (3/2)*x - (1/2))    # differentiated with respect to x
    a22 = -(1 - x)*((1/2) - x)            # differentiated with respect to y

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
        print('The fixed point in %s, %s is unstable' % (x,y))
    print('----------------------------')
    return lambda1 , lambda2

equilibria_manuale=np.array([[0,0],[0.5,1],[1,1]]) # perchè con equilibria non funziona siccome ci sono tutti i decimali, dice che non può fare sqrt() di un float...non so perchè! 
# iterate through list of fixed points
for x,y in equilibria_manuale:
    eigenvalues(x,y)    

#for item in equilibria:
#    eigenvalues(item[0],item[1])
#%%
#%%
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