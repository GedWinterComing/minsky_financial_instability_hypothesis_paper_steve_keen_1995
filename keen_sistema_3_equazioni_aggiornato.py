# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 18:46:32 2020

@author: Gabriele
"""

#%%
# VERSIONE INGEGNERIZZATA PER FARE QUANTE PIU' COSE AUTOMATICAMENTE ED
# EVITARE IL PROBLEMA DI NON TROVARE DEI PUNTI FISSI COME ACCADUTO NELLA RELAZIONE CON LA SADDLE

#%%
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sympy as sm
#from scipy.integrate import odeint
#%%
# Per disegnare frecce anche nei grafici 3D, poichè plt.arrow() funziona solo per caso 2D, vuole cioè 4 argomenti invece di 6; in più la sintassi è anche differente...
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)
#%%
from numpy import linalg as LA
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
def flusso_traiettoria_singola(x_iniz,y_iniz,z_iniz,dt,punti_di_tempo,FUNZA):
    X = np.array([x_iniz,y_iniz,z_iniz])
    risultati_integrazione = X
    for r in range(1,punti_di_tempo):
        X = runge_kutta(FUNZA,X,dt)
        risultati_integrazione = np.vstack((risultati_integrazione,X))
    x_flusso = risultati_integrazione[:,0]
    y_flusso = risultati_integrazione[:,1]
    z_flusso = risultati_integrazione[:,2]
    return x_flusso, y_flusso, z_flusso

    #np.vstack((A,b)) questo comando aggiunge il vettore b, inteso come vettore riga, in fondo alla matrice A aggiungendone una riga. Quindi se A è una matrice (m,n), np.vstack((A,b)) restituisce una matrice (m+1,n). Esiste anche il comando analogo per le colonne: np.hstack((A,b)) restituisce una matrice (m,n+1) con l'aggiunta di b inteso come vettore colonna
    
    #risultati_integrazione=np.array([risultati_integrazione,X])
    # QUESTO COMANDO NON FUNZIONAVA poichè non hai capito come funzionano gli array in Python!
    # Di fatto creavi una matrice con 2 sole righe invece di aggiungere il nuovo X come ultima riga della matrice precedente. La prima maxi riga era risultati_integrazione e la seconda era X. Cioè stavi creando un array composto da 2 array risultati_integrazione e X, invece di aggiungere l'array X agli ALTRI precedenti array!!!! risultati_integrazione viene considerato come un array nel suo complesso invece che come una collezione di tanti array
#%%
# SISTEMA PRINCIPALE DOVE CI SONO TUTTE LE ISTRUZIONI
# caso modello Minsky di Steve Keen (3D cioè 3 equazioni)

def sistema_finanziario_3_eq_minsky(X):
    #spacchettamento del vettore    
    omega = X[0]
    lanbda = X[1]
    d = X[2]
    #b = X[3]
    
    r = zeta + phi*d
    pi = 1 - omega - r*d
    
    domegadt = omega * (functionalform(A,B,C,D,lanbda) - alfa)
    dlanbdadt = lanbda * (((functionalform(E,F,G,H,(pi/ni)))/ni) - alfa - beta - gamma)
    dddt = r*d - pi + (ni - d) * (((functionalform(E,F,G,H,(pi/ni)))/ni) - gamma)
    #dddt = r*d - pi - (ni - d) * (((functionalform(E,F,G,H,(pi/ni)))/ni) - gamma) # col -(ni -d) come scritto nel paper, ma secondo me è sbagliato, facendo i calcoli ci vuole un +
    #dbdt = (phi*d + r) * (b - pi - (ni - d) * (((functionalform(E,F,G,H,(pi/ni)))/ni) - gamma))
    
    vettore_derivate = np.array([domegadt, dlanbdadt, dddt])
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
#%%
# PARAMETRI VARIABILI ED EQUAZIONI DEI PUNTI FISSI

#zeta = 4.808/100
#zeta=4.1/100
#zeta=3.1/100
zeta=2.5/100
#zeta = 0.046 # base rate cioè il tasso di interesse fissato dalla banca centrale
# a 0.026 sembra essere ancora stabile, punto di equilibrio con i 3 finiti, con phi a 0.005
#phi = 0.23/100
phi = 0
#phi = 0.005 # sensibilità al debito a pag.627 nel grafico scrive: .5% credo sia quindi lo 0.005


#lanbda0 = 0.35
#omega0 = 0.493
#omega0 = 0.96  # condizioni inizili preferite: le uniche mostrate da Keen
#lanbda0 = 0.9  # condizioni inizili preferite: le uniche mostrate da Keen

#d0 = 0  # condizioni inizili preferite: le uniche mostrate da Keen
#d0 = 0.1

######
# RETTA DI PUNTI INSTABILI CHE SI VERIFICA SOLO PER 2 r
# Indago punti fissi (0,lambda,d_{1 lambdanullo}) che non ha lo stesso d_1 chiamato "deq1" del punto fisso (omega_1, lambda_1, d_1) però ha gli stessi 2 pi greCo sopra chiamato pi1 che ha il segno - e quello con il segno +
# Questi punti esistono solo se vale la condizione: (1-ni*k^{-1} [ni*(alfa+beta+gamma)] )/r = (ni*k^{-1} [ni*(alfa+beta+gamma)]-ni*(alfa+beta))/(r-(alfa+beta)
# che risolta per r fornisce i 2 r_1 associati ciascuno al rispettivo pi:
# (nella relazione si potrebbe far vedere quella per il segno -)
pi1meno=((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ) # sotto è chiamato pi1
pi1più=((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) )
r_1lambdanullomeno = ((alfa+beta)*(pi1meno-1))/(2*pi1meno - 1 - ni*(alfa+beta))
r_1lambdanullopiù = ((alfa+beta)*(pi1più-1))/(2*pi1più - 1 - ni*(alfa+beta))
d_1lambdanullomeno = (pi1meno - ni*(alfa+beta))/(r_1lambdanullomeno - (alfa+beta))  # usato dt=0.001
d_1lambdanullopiù = (pi1più - ni*(alfa+beta))/(r_1lambdanullopiù - (alfa+beta)) # usato dt=0.001
# questo d_1 NON è deq1=1.2368049752642498
# poi dopo come punto iniziale si potrà fare d_1 + 0.01 per avere la dinamica del punto fisso instabile che agisce come "source"
#zeta = r_1lambdanullomeno
#zeta = r_1lambdanullopiù
######
# PUNTO DI EQUILIBRIO COME CONDIZIONE INIZIALE (quindi NO dinamica) NEL CASO phi=0 cioè r=zeta, NON funziona se r=zeta+phi*d perchè devi ritrovare le condizioni per omega1 e d1 in quanto, sostituendo r, comparirà un d1 da entrambi i lati 
r0=zeta # vedi più sopra dove viene scelto il zeta
pi1=((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ) # esiste anche un altro pi con segno + alla radice, però produce un deq1<0 tipo -378
omegaeq1=1-pi1-r0*((pi1-ni*(alfa+beta))/(r0-(alfa+beta)))
lanbdaeq1=(B/C) - ( np.sqrt(A/(alfa + D)) )/C
deq1=(pi1-ni*(alfa+beta))/(r0-(alfa+beta))
#omega0=omegaeq1
#lanbda0=lanbdaeq1
#d0=deq1
#omega0=round(omegaeq1,9)  # Non funziona perchè sono sufficienti quelle poche cifre tagliate via per creare spirali infinitesime quindi meglio direttamente prendere tutti i valori
#lanbda0=round(lanbdaeq1,9)
#d0=round(deq1,9)
# ovvero questi sono, però attento che dipendono dai valori di zeta scelto: (per zeta=4.1/100)
#omega0=0.81042224079154401
#lanbda0=0.96588112190688502
#d0=1.2368049752642498

# IN REALTA' E' pi1=((ni*F)/G) -o+ (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) )
# E QUINDI CI SONO DUE COPPIE DI (omegaeq1,deq1) PERO' CON pi1 CON SEGNO + VIENE deq1<0 ad esempio in questo caso deq1=-26.79 PERO' FORSE DEVE ESSERE RIPORTATO NEL GRAFICO SPAZIO DELLE FASI!

#(sul foglio inizialmente lo avevi trovato così omegaeq1=((alfa+beta)*(1-r0*ni))/(alfa+beta-r0) - r0/(alfa+beta-r0) + ((2*r0-alfa-beta)/(alfa+beta-r0))*( ((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ) )    )
#(sul foglio inizialmente lo avevi trovato così deq1=(1-ni*(alfa+beta))/(2*r0-alfa-beta) - (1/(alfa+beta-r0))*( ((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ) ) - (1/(2*r0-alfa-beta))*( ( ((alfa+beta)*(1-r0*ni))/(alfa+beta-r0) )-(r0/(alfa+beta-r0)) )    )
######
#
######
# PUNTO DI EQUILIBRIO COME CONDIZIONE INIZIALE (quindi NO dinamica) NEL CASO r=zeta+phi*d
#pi1=((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) )
#lanbdaeq1=(B/C) - ( np.sqrt(A/(alfa + D)) )/C
#deq1_con_phi_segno_positivo=(-(zeta-alfa-beta) + np.sqrt( (zeta-alfa-beta)**2 - 4*phi*(ni*(alfa+beta)-pi1) ))/(2*phi)
#deq1_con_phi_segno_negativo=(-(zeta-alfa-beta) - np.sqrt( (zeta-alfa-beta)**2 - 4*phi*(ni*(alfa+beta)-pi1) ))/(2*phi)
#omegaeq1_con_phi_segno_positivo=1-pi1-zeta*deq1_con_phi_segno_positivo-phi*deq1_con_phi_segno_positivo**2
#omegaeq1_con_phi_segno_negativo=1-pi1-zeta*deq1_con_phi_segno_negativo-phi*deq1_con_phi_segno_negativo**2

#omega0=omegaeq1_con_phi_segno_positivo
#omega0=omegaeq1_con_phi_segno_negativo
#lanbda0=lanbdaeq1
#d0=deq1_con_phi_segno_positivo
#d0=deq1_con_phi_segno_negativo
# ovvero questi sono, però attento che dipendono dai valori di zeta e phi scelti: (per zeta=3.1/100 e phi = 0.23/100)
#omega0=0.49094930433037387  # segno positivo
#omega0=0.84053219651826938  # segno negativo
#lanbda0=0.96588112190688502
#d0=7.6262637044876511  # segno positivo
#d0=0.63460586072974201  # segno negativo
#pi1=0.13886875522262176
######
#%%
# PUNTI FISSI

# Per trovare automaticamente gli zeri del sistema di equazioni differenziali
# define the system in this way (asuming a system with no negative values)
# to avoid interference use u, v, w and U, V, W instead of omega, lanbda, d
# r = zeta + phi*d

u, v, w = sm.symbols('u, v, w', negative=False)
U = u * (functionalform(A,B,C,D,v) - alfa)
V = v * (((functionalform(E,F,G,H,((1-u-zeta*w)/ni)))/ni) - alfa - beta - gamma)
W = 2*zeta*w - 1 + u + (ni - w) * (((functionalform(E,F,G,H,((1-u-zeta*w)/ni)))/ni) - gamma)

#    dddt = r*d - pi + (ni - d) * (((functionalform(E,F,G,H,(pi/ni)))/ni) - gamma)
#    dddt = 2*r*d - 1 + omega + (ni - d) * (((functionalform(E,F,G,H,(pi/ni)))/ni) - gamma)

# use sympy's way of setting equations to zero
UEqual = sm.Eq(U, 0)
VEqual = sm.Eq(V, 0)
WEqual = sm.Eq(W, 0)

# compute fixed points
equilibria = sm.solve( (UEqual, VEqual, WEqual), u, v, w )
print(equilibria)
# NON FUNZIONA PROBABILMENTE POICHE' UNA DELLE SOLUZIONI E' INFINITO
#%%
# PUNTI FISSI DEL SISTEMA CAMBIO VARIABILI u=1/d

# Per trovare automaticamente gli zeri del sistema di equazioni differenziali
# define the system in this way (asuming a system with no negative values)
# to avoid interference use u, v, w and U, V, W instead of omega, lanbda, d
# r = zeta + phi*d

u, v, w = sm.symbols('u, v, w', negative=False)
U = u * (functionalform(A,B,C,D,v) - alfa)
V = v * (((functionalform(E,F,G,H,((1-u-(zeta/w))/ni)))/ni) - alfa - beta - gamma)
W = w**2 - u*w**2 - 2*zeta*w - (ni*w**2 - w) * (((functionalform(E,F,G,H,((1-u-(zeta/w))/ni)))/ni) - gamma)

#    dddt = r*d - pi + (ni - d) * (((functionalform(E,F,G,H,(pi/ni)))/ni) - gamma)
#    dddt = 2*r*d - 1 + omega + (ni - d) * (((functionalform(E,F,G,H,(pi/ni)))/ni) - gamma)

# use sympy's way of setting equations to zero
UEqual = sm.Eq(U, 0)
VEqual = sm.Eq(V, 0)
WEqual = sm.Eq(W, 0)

# compute fixed points
equilibria = sm.solve( (UEqual, VEqual, WEqual), u, v, w )
print(equilibria)
# NON FUNZIONA PROBABILMENTE POICHE' UNA DELLE SOLUZIONI E' INFINITO, PERO' QUA DOVREBBE ESSERE 0 quindi andare bene
#%%
# PUNTI FISSI SOLO DELLA TERZA EQUAZIONE cerco cioè solo d_0 del punto fisso (0,0,d_0)

# Per trovare automaticamente gli zeri del sistema di equazioni differenziali
# define the system in this way (asuming a system with no negative values)
# to avoid interference use u, v, w and U, V, W instead of omega, lanbda, d
# r = zeta + phi*d

w = sm.symbols('w', negative=False)
W = 2*zeta*w - 1 + (ni - w) * (((functionalform(E,F,G,H,((1-zeta*w)/ni)))/ni) - gamma)

#    dddt = r*d - pi + (ni - d) * (((functionalform(E,F,G,H,(pi/ni)))/ni) - gamma)
#    dddt = 2*r*d - 1 + omega + (ni - d) * (((functionalform(E,F,G,H,(pi/ni)))/ni) - gamma)

# use sympy's way of setting equations to zero
WEqual = sm.Eq(W, 0)

# compute fixed points
equilibria_d_0 = sm.solve( (WEqual), w )
print(equilibria_d_0)
#%%
# TESTO CHE i d_0 TROVATI (per il valore di r, o zeta, scelto sopra) DEL PUNTO FISSO (0,0,d_0) SIANO EFFETTIVAMENTE SOLUZIONI della TERZA EQUAZIONE dd/dt=0

# inserisco manualmente i punti fissi equilibria_d_0 calcolati per r = zeta = r_1lambdanullomeno
debito0annullaterzaeq1 = 8.62402480510675    # se fossero diversi poichè hai cambiato zeta e non stai più usando r_1lambdanullomeno,
debito0annullaterzaeq2 = 11.7241196652734    # dovrai ovviamente inserire manualmente quelli giusti
debito0annullaterzaeq3 = 17.4452497910951    # Infatti sotto ho messo il controllo che l'r usato sia quello della retta (0,lambda,d_1)
debito0annullaterzaeq_vettore = np.array([debito0annullaterzaeq1,debito0annullaterzaeq2,debito0annullaterzaeq3])
terzaeqsiannullaperd0_punto_di_domanda_vettore = 2*zeta*debito0annullaterzaeq_vettore - 1 + (ni - debito0annullaterzaeq_vettore) * (((functionalform(E,F,G,H,((1 - zeta*debito0annullaterzaeq_vettore)/ni)))/ni) - gamma)
# un risultato del tipo -1.35724765e-14 significa -1.347*10^{-14} cioè in pratica è zero. Nella terminologia di Python "e" sighifica "moltiplicato per 10^{}"
r_1lambdanullomeno==zeta
#%%
#%%
# FACCIO I GRAFICI PER IL CASO r=4.9..% CIOE' QUELLO PER IL QUALE ESISTE LA RETTA (0,lambda,d_1) DI PUNTI FISSI
#%%
#%%       
# INTEGRAZIONE

# definisco le condizioni iniziali di ogni flusso ma è ridondante perchè potrei direttamente farlo nella chiamata di flusso_traiettoria_singola(x_0,y_0,..,..) mettendo direttamente i numeri

# quello che parte dalla retta di punti fissi instabili discostandosi di un epsilon piccolo, sennò ci rimarrebbe
x0_1 = 0
y0_1 = 0.8
z0_1 = d_1lambdanullomeno + 0.01
# un altro sempre che parte dalla retta
x0_2 = 0 
y0_2 = 0.3 
z0_2 = d_1lambdanullomeno + 0.01
# parte a destra e va al secondo dei punti (0,0,d_0) che è quello FISSO
x0_3 = 0.14
y0_3 = 0
z0_3 = 11.7241196652734
# Va a -oo quasi direttamente. NON LO GRAFICO perchè abbastanza inutile da mostrare forse convertire questo usandolo per mostrare il braccio di qualche sella
#x0_3 = 0.3
#y0_3 = 0.3
#z0_3 = 0.3
x0_4 = 0.5
y0_4 = 0.5
z0_4 = 0.5

# QUELLO SCELTO DA STEVE KEEN
x0_5 = 0.96
y0_5 = 0.9
z0_5 = 0
# parte vicino a quello di Keen ma con un debito prossimo al punto fisso stabile che però siccome r>4.6 potrebbe non essere raggiunto e andare a +oo, crollo economia
x0_6 = 0.96
y0_6 = 0.9
z0_6 = deq1 - 0.1
x0_7 = omegaeq1
y0_7 = lanbdaeq1 
z0_7 = deq1 - 0.2
# fa casino nel vortice, credo
x0_8 = 1.4
y0_8 = 0.7
z0_8 = 0
# va a +oo
x0_9 = 0.6
y0_9 = 0.6
z0_9 = 0

# fa un mulinello e va a +oo
x0_10 = 1
y0_10 = 0.99
z0_10 = 1
# parte a destra e vicino al terzo dei punti (0,0,d_0) che è instabile
x0_11 = 0.14
y0_11 = 0
z0_11 = 17.4452497910951
# provo a popolare l'area a destra. NON LO GRAFICO PERCHE' MI SBALLA LA SCALA di omega COMUNQUE VA A +OO
#x0_11 = 0.8
#y0_11 = 0.2
#z0_11 = 0
# parte a destra vicino al primo dei punti (0,0,d_0) che è instabile e va a quello stabile nel mezzo
x0_12 = 0.14
y0_12 = 0
z0_12 = 8.62402480510675
# provo a popolare l'area a destra. NON LO GRAFICO
#x0_12 = 1.2
#y0_12 = 0.8
#z0_12 = 0
# popoliamo area in basso a destra
x0_13 = 1.2
y0_13 = 0.4
z0_13 = 3
# partire vicino al più piccolo dei punti (0,0,d0) per capire cosa succede
x0_14 = 0
y0_14 = 0.9
z0_14 = 8.6224
# va a -oo partendo da sotto il punto di equilibrio (omega1,lambda1,d1)
x0_15 = omegaeq1
y0_15 = lanbdaeq1
z0_15 = 1.5
# il punto fisso che era equilibrio desiderabile credo sia instabile perchè partendoci vicino va +oo
x0_16 = omegaeq1
y0_16 = lanbdaeq1
z0_16 = deq1 + 0.01
# provi a popolare la zona del flusso 6 facendolo partire poco più avanti del 6
x0_17 = 1.4
y0_17 = 0.96
z0_17 = deq1 - 0.1
# lo metto sulla retta rossa a lambda=0, va a +oo
x0_18 = 0
y0_18 = 0
z0_18 = d_1lambdanullomeno + 0.01
# ottimi tenere
x0_19 = 1.5
y0_19 = 0
z0_19 = 0
x0_20 = 1.5
y0_20 = 0.5
z0_20 = 0.5


#dt = 0.01
#numero_di_punti_tempo_voluti = 5000  # deve essere un numero tale che tempo_finale venga un numero intero
#dt = 0.01
dt=0.03 # QUELLO USATO! GIUSTO!
#dt=0.005
#dt=0.004 # o forse meglio 0.002
#numero_di_punti_tempo_voluti = 50000
#numero_di_punti_tempo_voluti = 25000
numero_di_punti_tempo_voluti = 15000
tempo_finale = numero_di_punti_tempo_voluti * dt
t = np.linspace(0,tempo_finale,numero_di_punti_tempo_voluti) # per fare il grafico dopo, l'ultimo è proprio tempo_finale

# ovviamente tutto ciò si poteva saltare mettendo direttamente i numeri al posto di x0_i e y0_i definendoli qui e non sopra che è ridondante
x_1_flusso, y_1_flusso, z_1_flusso = flusso_traiettoria_singola(x0_1,y0_1,z0_1,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_2_flusso, y_2_flusso, z_2_flusso = flusso_traiettoria_singola(x0_2,y0_2,z0_2,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_3_flusso, y_3_flusso, z_3_flusso = flusso_traiettoria_singola(x0_3,y0_3,z0_3,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_4_flusso, y_4_flusso, z_4_flusso = flusso_traiettoria_singola(x0_4,y0_4,z0_4,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_5_flusso, y_5_flusso, z_5_flusso = flusso_traiettoria_singola(x0_5,y0_5,z0_5,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_6_flusso, y_6_flusso, z_6_flusso = flusso_traiettoria_singola(x0_6,y0_6,z0_6,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_7_flusso, y_7_flusso, z_7_flusso = flusso_traiettoria_singola(x0_7,y0_7,z0_7,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_8_flusso, y_8_flusso, z_8_flusso = flusso_traiettoria_singola(x0_8,y0_8,z0_8,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_9_flusso, y_9_flusso, z_9_flusso = flusso_traiettoria_singola(x0_9,y0_9,z0_9,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_10_flusso, y_10_flusso, z_10_flusso = flusso_traiettoria_singola(x0_10,y0_10,z0_10,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_11_flusso, y_11_flusso, z_11_flusso = flusso_traiettoria_singola(x0_11,y0_11,z0_11,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_12_flusso, y_12_flusso, z_12_flusso = flusso_traiettoria_singola(x0_12,y0_12,z0_12,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_13_flusso, y_13_flusso, z_13_flusso = flusso_traiettoria_singola(x0_13,y0_13,z0_13,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_14_flusso, y_14_flusso, z_14_flusso = flusso_traiettoria_singola(x0_14,y0_14,z0_14,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_15_flusso, y_15_flusso, z_15_flusso = flusso_traiettoria_singola(x0_15,y0_15,z0_15,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_16_flusso, y_16_flusso, z_16_flusso = flusso_traiettoria_singola(x0_16,y0_16,z0_16,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_17_flusso, y_17_flusso, z_17_flusso = flusso_traiettoria_singola(x0_17,y0_17,z0_17,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_18_flusso, y_18_flusso, z_18_flusso = flusso_traiettoria_singola(x0_18,y0_18,z0_18,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_19_flusso, y_19_flusso, z_19_flusso = flusso_traiettoria_singola(x0_19,y0_19,z0_19,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_20_flusso, y_20_flusso, z_20_flusso = flusso_traiettoria_singola(x0_20,y0_20,z0_20,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
#%%
# GRAFICO SPAZIO DELLE FASI
font = {'family':'serif','color':'darkred','weight':'normal','size': 10.7,}
# plt.text(2, 0.65, r'$\cos(2 \pi t) \exp(-t)$', fontdict=font)


# FRECCE DEL CAMPO VETTORIALE

# creo la griglia, in realtà saranno 2 matrici una con coordinata x e l'altra y così prese elemento per elemento danno il punto (x_i,y_i)
x_griglia = np.linspace(0, 2, 6) # buono anche solo 14 però poi le frecce non sono tutte ortogonali come dovrebbe invece essere
y_griglia = np.linspace(0, 2, 6)
z_griglia = np.linspace(0, 2, 6)

# il comando sotto è interessante perchè assegna automaticamente i termini in ordine: es. gaia=np.array([3,4]) poi con hippie1,hippie2=gaia si ha hippie1=gaia[0]=3 e hippie2=gaia[1]=4
X1 , Y1, Z1  = np.meshgrid(x_griglia, y_griglia, z_griglia)    # create a grid
# ritorna una prima matrice che ha come righe il vettore x ripetuto tante volte quante la lunghezza del vettore y, cioè len(y) volte, e una seconda matrice che ha il vettore y come colonne ripetuto tante volte quanta la lunghezza del vettore x, len(x)
DX1, DY1, DZ1 = sistema_finanziario_3_eq_minsky([X1, Y1, Z1])    # compute growth rate on the grid
# la moltiplicazione * usata in "sistema_II_semplificato" è quella normale e se X come in questo caso è composto da 2 matrici, quando spacchetta X assegna la prima matrice (cioè le coordinate x su tutta la griglia) a x e la seconda (cioè la griglia per le coordinate y) a y, quando poi si trova a fare x*y moltiplica elemento per elemento entrambe le matrici senza fare il prodotto matriciale! Perciò le due matrici devono avere entrambe la stessa dimensione, tutte e due devono essere nxn
M = (np.hypot(DX1, DY1, DZ1))                        # norm growth rate 
M[ M == 0] = 1.                                 # avoid zero division errors 
DX1 /= M                                        # normalize each arrows
DY1 /= M
# se fai M.shape dà (20,20) è infatti una matrice 20x20 siccome hai usato due volte x = np.linspace(0, 2, 20) o x e y, ma sono uguali x=y=np.linspace(0, 2, 20)


fig = plt.figure(num=None, figsize=(8, 6), dpi=800, facecolor='w', edgecolor='k')
ax = Axes3D(fig)


# ATTIVARE SE VUOI VEDERE LE FRECCE DEL CAMPO VETTORIALE ma non penso che Turchetti le voglia
#ax.quiver(X1, Y1, Z1, DX1, DY1, DZ1, color='deepskyblue', length=0.1, pivot='middle') # anche senza length=0.1 magari
#plt.quiver(X1, Y1, Z1, DX1, DY1, DZ1, M, pivot='mid') # M gli dà una scala di colori che sono i numeri

# per non graficare tutti i punti sennò vanno a -oo e +oo 
#mostrafinoa = 13540
mostrafinoa = 7540
mostrafinoa1 = 400 # per flussi che divergono subito BUONO PARE
#mostrafinoa1 = 3400 # POI CANCELLARE L'HAI MESSO SOLO PER TROVARE QUEL SEGMENTO CHE PARTE DALLA FINE DELLA RETTA DOVE C'E' PUNTO (omega_1,lambda_1,d_1) E CAPIRE DOVE VA A FINIRE. Sei arrivato al 6 flusso da commentare e vedere se quello sparisce che vuoi te E' il 7
mostrafinoa2 = 700 # per flussi che hanno bisogno di farsi vedere un poco di più
mostrafinoa3 = 3080 # per il flusso 14
mostrafinoa4 = 150 # per il flusso 11 che va a +oo pare
mostrafinoa5 = 230 # per il flusso 17 che va a +oo pare
mostrafinoa6 = 860 # per il flusso 13 che va a +oo pare

#ax.plot_wireframe(x_1_flusso, y_1_flusso, z_1_flusso, color='green')
ax.plot_wireframe(x_1_flusso[0:mostrafinoa1], y_1_flusso[0:mostrafinoa1], z_1_flusso[0:mostrafinoa1], color='green')
#ax.plot_wireframe(x_2_flusso, y_2_flusso, z_2_flusso, color='green')
ax.plot_wireframe(x_2_flusso[0:mostrafinoa1], y_2_flusso[0:mostrafinoa1], z_2_flusso[0:mostrafinoa1], color='green')
ax.plot_wireframe(x_3_flusso, y_3_flusso, z_3_flusso, color='green')
#ax.plot_wireframe(x_3_flusso[0:mostrafinoa1], y_3_flusso[0:mostrafinoa1], z_3_flusso[0:mostrafinoa1], color='green')
#ax.plot_wireframe(x_4_flusso, y_4_flusso, z_4_flusso, color='green')
ax.plot_wireframe(x_4_flusso[0:mostrafinoa1], y_4_flusso[0:mostrafinoa1], z_4_flusso[0:mostrafinoa1], color='green')
#ax.plot_wireframe(x_5_flusso, y_5_flusso, z_5_flusso, color='green')
ax.plot_wireframe(x_5_flusso[0:mostrafinoa], y_5_flusso[0:mostrafinoa], z_5_flusso[0:mostrafinoa], color='green')
#ax.plot_wireframe(x_6_flusso, y_6_flusso, z_6_flusso, color='green')
ax.plot_wireframe(x_6_flusso[0:mostrafinoa1-140], y_6_flusso[0:mostrafinoa1-140], z_6_flusso[0:mostrafinoa1-140], color='green')
ax.plot_wireframe(x_7_flusso, y_7_flusso, z_7_flusso, color='green')
#ax.plot_wireframe(x_7_flusso[0:mostrafinoa], y_7_flusso[0:mostrafinoa], z_7_flusso[0:mostrafinoa], color='green')
#ax.plot_wireframe(x_8_flusso, y_8_flusso, z_8_flusso, color='green')
#ax.plot_wireframe(x_8_flusso[0:mostrafinoa2], y_8_flusso[0:mostrafinoa2], z_8_flusso[0:mostrafinoa2], color='green')
#ax.plot_wireframe(x_9_flusso, y_9_flusso, z_9_flusso, color='green')
ax.plot_wireframe(x_9_flusso[0:mostrafinoa1], y_9_flusso[0:mostrafinoa1], z_9_flusso[0:mostrafinoa1], color='green')
#ax.plot_wireframe(x_10_flusso, y_10_flusso, z_10_flusso, color='green')
ax.plot_wireframe(x_10_flusso[0:mostrafinoa6+1500], y_10_flusso[0:mostrafinoa6+1500], z_10_flusso[0:mostrafinoa6+1500], color='green')
#ax.plot_wireframe(x_11_flusso, y_11_flusso, z_11_flusso, color='green')
#ax.plot_wireframe(x_11_flusso[0:mostrafinoa5], y_11_flusso[0:mostrafinoa5], z_11_flusso[0:mostrafinoa5], color='green')
#ax.plot_wireframe(x_12_flusso, y_12_flusso, z_12_flusso, color='green')
#ax.plot_wireframe(x_12_flusso[0:mostrafinoa1], y_12_flusso[0:mostrafinoa1], z_12_flusso[0:mostrafinoa1], color='green')
#ax.plot_wireframe(x_13_flusso, y_13_flusso, z_13_flusso, color='green')
ax.plot_wireframe(x_13_flusso[0:mostrafinoa6], y_13_flusso[0:mostrafinoa6], z_13_flusso[0:mostrafinoa6], color='green')
#ax.plot_wireframe(x_14_flusso, y_14_flusso, z_14_flusso, color='green')
ax.plot_wireframe(x_14_flusso[0:mostrafinoa3], y_14_flusso[0:mostrafinoa3], z_14_flusso[0:mostrafinoa3], color='green')
#ax.plot_wireframe(x_15_flusso, y_15_flusso, z_15_flusso, color='green')
ax.plot_wireframe(x_15_flusso[0:mostrafinoa1], y_15_flusso[0:mostrafinoa1], z_15_flusso[0:mostrafinoa1], color='green')
#ax.plot_wireframe(x_16_flusso, y_16_flusso, z_16_flusso, color='green')
ax.plot_wireframe(x_16_flusso[0:mostrafinoa1], y_16_flusso[0:mostrafinoa1], z_16_flusso[0:mostrafinoa1], color='green')
#ax.plot_wireframe(x_17_flusso, y_17_flusso, z_17_flusso, color='green')
ax.plot_wireframe(x_17_flusso[0:mostrafinoa5], y_17_flusso[0:mostrafinoa5], z_17_flusso[0:mostrafinoa5], color='green')
#ax.plot_wireframe(x_18_flusso, y_18_flusso, z_18_flusso, color='green')
ax.plot_wireframe(x_18_flusso[0:mostrafinoa1], y_18_flusso[0:mostrafinoa1], z_18_flusso[0:mostrafinoa1], color='green')
#ax.plot_wireframe(x_19_flusso, y_19_flusso, z_19_flusso, color='green')
ax.plot_wireframe(x_19_flusso[0:mostrafinoa1], y_19_flusso[0:mostrafinoa1], z_19_flusso[0:mostrafinoa1], color='green')
#ax.plot_wireframe(x_20_flusso, y_20_flusso, z_20_flusso, color='green')
ax.plot_wireframe(x_20_flusso[0:mostrafinoa6+500], y_20_flusso[0:mostrafinoa6+500], z_20_flusso[0:mostrafinoa6+500], color='green')
# ALLA FINE TOLTI FLUSSI 8,11 e 12 E RISPETTIVE FRECCE. Invece 3 cambiato e lasciato, converge al fisso d_0 stabile, se lo sposti di poco va a -oo


# Nel grafico che sto creando metto anche tutti punti fissi. Ci sono 2 modi possibili per farlo
# versione "automatica" però non puoi variare i colori ma qua non funziona perchè non riesce a trovarli:
#for point in equilibria:
#    plt.plot(point[0],point[1],"royalblue", marker = "o", markersize = 7.0)

# versione "manuale":
ax.scatter(omegaeq1, lanbdaeq1, deq1, 'o',s=30, c='royalblue',depthshade=True)
ax.plot_wireframe([0, 0], [0, 1], [d_1lambdanullomeno, d_1lambdanullomeno], lw=2, color='red') # i punti fissi instabili della soluzione (0,lambda,d_1)
ax.scatter([0, 0, 0], [0, 0, 0], [8.62402480510675, 11.7241196652734, 17.4452497910951], 'o',s=30, c='red',depthshade=True) # i punti fissi instabili della soluzione (0,0,d_0)

#plt.plot(0,0,'ob')
#plt.plot(1 - ((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),(B/C) - ( np.sqrt(A/(alfa + D)) )/C,'o',color='royalblue')   # condizione iniziale o meglio punto del centro dei cicli. Basta allontanarsi poco da questo punto con la condizione iniziale che si formano i cycle!
#plt.plot(1 - ((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),(B/C) - ( np.sqrt(A/(alfa + D)) )/C,'o',color='royalblue')
#plt.plot(1 - ((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),(B/C) + ( np.sqrt(A/(alfa + D)) )/C,'o',color='royalblue') 
#plt.plot(1 - ((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),(B/C) + ( np.sqrt(A/(alfa + D)) )/C,'o',color='royalblue')


#ip=1300 # vedere se va bene
ip=280
ip1=100
ip2=300
ip3=1820
ip4=3840
ip5=200
ip6=2300 # per il flusso 10


# E' diversa dalla funzione arrow 2D standard, quella usava x,y e dx, dy questa usa il punto iniziale della freccia (x[i],y[i],z[i]) e quello finale (x[i+1],y[i+1],z[i+1]) inoltre li vuole a coppie come vettori cioè ciascuno così [..,..] e non di seguito arrow(x[i],y[i],z[i],x[i+1],y[i+1],z[i+1]) come lo stile dell'altra
a1 = Arrow3D([x_1_flusso[ip],x_1_flusso[ip+1]],[y_1_flusso[ip],y_1_flusso[ip+1]],[z_1_flusso[ip],z_1_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green") # arrowstyle="-|>" per fare la freccia piena
ax.add_artist(a1) # serve per disegnare la freccia
a2 = Arrow3D([x_2_flusso[ip],x_2_flusso[ip+1]],[y_2_flusso[ip],y_2_flusso[ip+1]],[z_2_flusso[ip],z_2_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green") # arrowstyle="-|>" per fare la freccia piena
ax.add_artist(a2) # serve per disegnare la freccia
#a3 = Arrow3D([x_3_flusso[ip2],x_3_flusso[ip2+1]],[y_3_flusso[ip2],y_3_flusso[ip2+1]],[z_3_flusso[ip2],z_3_flusso[ip2+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green") # arrowstyle="-|>" per fare la freccia piena
#ax.add_artist(a3) # serve per disegnare la freccia
a4 = Arrow3D([x_4_flusso[ip],x_4_flusso[ip+1]],[y_4_flusso[ip],y_4_flusso[ip+1]],[z_4_flusso[ip],z_4_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green") # arrowstyle="-|>" per fare la freccia piena
ax.add_artist(a4) # serve per disegnare la freccia
a5 = Arrow3D([x_5_flusso[ip],x_5_flusso[ip+1]],[y_5_flusso[ip],y_5_flusso[ip+1]],[z_5_flusso[ip],z_5_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green") # arrowstyle="-|>" per fare la freccia piena
ax.add_artist(a5) # serve per disegnare la freccia
a6 = Arrow3D([x_6_flusso[ip1],x_6_flusso[ip1+1]],[y_6_flusso[ip1],y_6_flusso[ip1+1]],[z_6_flusso[ip1],z_6_flusso[ip1+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a6)
a7 = Arrow3D([x_7_flusso[ip],x_7_flusso[ip+1]],[y_7_flusso[ip],y_7_flusso[ip+1]],[z_7_flusso[ip],z_7_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a7)
#a8 = Arrow3D([x_8_flusso[ip],x_8_flusso[ip+1]],[y_8_flusso[ip],y_8_flusso[ip+1]],[z_8_flusso[ip],z_8_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
#ax.add_artist(a8)
a9 = Arrow3D([x_9_flusso[ip],x_9_flusso[ip+1]],[y_9_flusso[ip],y_9_flusso[ip+1]],[z_9_flusso[ip],z_9_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a9)
a10 = Arrow3D([x_10_flusso[ip6],x_10_flusso[ip6+1]],[y_10_flusso[ip6],y_10_flusso[ip6+1]],[z_10_flusso[ip6],z_10_flusso[ip6+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a10)
#a11 = Arrow3D([x_11_flusso[ip],x_11_flusso[ip+1]],[y_11_flusso[ip],y_11_flusso[ip+1]],[z_11_flusso[ip],z_11_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
#ax.add_artist(a11)
#a12 = Arrow3D([x_12_flusso[ip],x_12_flusso[ip+1]],[y_12_flusso[ip],y_12_flusso[ip+1]],[z_12_flusso[ip],z_12_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
#ax.add_artist(a12)
a13 = Arrow3D([x_13_flusso[ip],x_13_flusso[ip+1]],[y_13_flusso[ip],y_13_flusso[ip+1]],[z_13_flusso[ip],z_13_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a13)
a14 = Arrow3D([x_14_flusso[ip],x_14_flusso[ip+1]],[y_14_flusso[ip],y_14_flusso[ip+1]],[z_14_flusso[ip],z_14_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a14)
a15 = Arrow3D([x_15_flusso[ip],x_15_flusso[ip+1]],[y_15_flusso[ip],y_15_flusso[ip+1]],[z_15_flusso[ip],z_15_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a15)
a16 = Arrow3D([x_16_flusso[ip],x_16_flusso[ip+1]],[y_16_flusso[ip],y_16_flusso[ip+1]],[z_16_flusso[ip],z_16_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a16)
a17 = Arrow3D([x_17_flusso[ip1],x_17_flusso[ip1+1]],[y_17_flusso[ip1],y_17_flusso[ip1+1]],[z_17_flusso[ip1],z_17_flusso[ip1+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a17)
a18 = Arrow3D([x_18_flusso[ip],x_18_flusso[ip+1]],[y_18_flusso[ip],y_18_flusso[ip+1]],[z_18_flusso[ip],z_18_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a18)
a19 = Arrow3D([x_19_flusso[ip],x_19_flusso[ip+1]],[y_19_flusso[ip],y_19_flusso[ip+1]],[z_19_flusso[ip],z_19_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a19)
a20 = Arrow3D([x_20_flusso[ip],x_20_flusso[ip+1]],[y_20_flusso[ip],y_20_flusso[ip+1]],[z_20_flusso[ip],z_20_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a20)


# LE FRECCE DEGLI ASSI QUA NON SERVONO
#plt.arrow(0.5, 0, -0.01, 0,shape='full', lw=0, length_includes_head=True, head_width=.04)
#plt.arrow(1.1, 0, -0.01, 0,shape='full', lw=0, length_includes_head=True, head_width=.04)
#plt.arrow(0, 0.9, 0, -0.01,shape='full', lw=0, length_includes_head=True, head_width=.05)
#plt.arrow(0, 0.4, 0, -0.01,shape='full', lw=0, length_includes_head=True, head_width=.05)


ax.set_title('$r=$''%g' % zeta, fontdict = font)
#ax.set_title('Modello Minsky basso r spazio delle fasi \n $\zeta=$' + '%g' % zeta + ', $\phi=$ ' + '%g' % phi, fontdict = font)
ax.set_xlabel('quota salari $\omega$')
ax.set_ylabel('occupazione $\lambda$')
ax.set_zlabel('quota debito privato $d$')
#ax.set_xlim(-0.5,2)
#ax.set_ylim(-0.5,1.5)
ax.set_zlim(-20,30)
#ax.set_aspect("equal") # serve ad usare un'unica scala vedere se usare o meno, far la prova se piace o no: se attiva schiacci la figura

#plt.axis([0, 1.4, 0, 2.1]) # X=1.4 buono # PER FARE ANCHE PARTE SOPRA CON OCCUPATI>POPOLAZIONE cioì lambda>1
#plt.axis([0, 1.2, .9, 1.1]) # zoom sulla parte incasinata
#plt.axis([0, 1.4, 0, 1.0]) # PER FARE GRAFICO SOLO DELLA PARTE SENSATA cioè con lambda NON superiore 1
#plt.savefig('relazione corretta Minsky alto r con retta lambda spazio delle fasi.jpg', dpi=850, transparent=False)
#%%
#%% # PROVA GRAFICO e PER CAPIRE DOVE VANNO LE ORBITE POI SI POTREBBE ANCHE ELIMINARE  ???

#x0_5 = 0.96
#y0_5 = 0.9
#z0_5 = 1.7

# indago punti fissi (0,lambda,d_{1 speciale}) che non ha lo stesso d_1 sopra chiamato deq1 del punto fisso (omega_1, lambda_1, d_1) però ha gli stessi 2 pi grego sopra chiamato pi1 che ha il segno - e quello con il segno +
# la condizione per r è
pi1meno=((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) )
pi1più=((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) )
r_1lambdanullomeno = ((alfa+beta)*(pi1meno-1))/(2*pi1meno - 1 - ni*(alfa+beta))
r_1lambdanullopiù = ((alfa+beta)*(pi1più-1))/(2*pi1più - 1 - ni*(alfa+beta))
#x0_5 = 0
#y0_5 = 0.95
#z0_5 = (pi1meno - ni*(alfa+beta))/(r_1lambdanullomeno - (alfa+beta)) + 0.01 # usato dt=0.001
#z0_5 = (pi1più - ni*(alfa+beta))/(r_1lambdanullopiù - (alfa+beta)) + 0.01 # usato dt=0.001
#zeta=r_1lambdanullomeno
#zeta=r_1lambdanullopiù
#questo z NON è deq1=1.2368049752642498

# 8.62402480510675; 11.7241196652734; 17.4452497910951 sono i tre punti (0,0,d_0) per r=z=r_1lambdanullomeno
# 13.2796458147542; 22.4204384888053; 35.3726429691678 sono i tre punti (0,0,d_0) per r=z=2.5%
x0_45 = 0.96
y0_45 = 0.9
z0_45 = 0

x0_altro = 0.96
y0_altro = 0.9
z0_altro = 0


#x0_45 = 0.96
#y0_45 = 0.9
#z0_45 = 0 # con dt = 0.020196 e r_1lambdanullomeno pare diventare negativo forse andare a meno infinito e poi dopo a + infinito
#z0_5=9.80
#zeta = 4.6/100
#phi = 0
#dt=0.001
#dt=0.005
dt = 0.03
#dt = 0.020196
#dt=0.00028
#numero_di_punti_tempo_voluti = 15000
numero_di_punti_tempo_voluti = 25000
tempo_finale = numero_di_punti_tempo_voluti * dt
t = np.linspace(0,tempo_finale,numero_di_punti_tempo_voluti) # per fare il grafico dopo, l'ultimo è proprio tempo_finale
x_griglia = np.linspace(0, 1, 6) # buono anche solo 14 però poi le frecce non sono tutte ortogonali come dovrebbe invece essere
y_griglia = np.linspace(0, 1, 6)
z_griglia = np.linspace(0, 1, 6)
X1 , Y1, Z1  = np.meshgrid(x_griglia, y_griglia, z_griglia)    # create a grid
# ritorna una prima matrice che ha come righe il vettore x ripetuto tante volte quante la lunghezza del vettore y, cioè len(y) volte, e una seconda matrice che ha il vettore y come colonne ripetuto tante volte quanta la lunghezza del vettore x, len(x)
DX1, DY1, DZ1 = sistema_finanziario_3_eq_minsky([X1, Y1, Z1])    # compute growth rate on the grid

x_45_flusso, y_45_flusso, z_45_flusso = flusso_traiettoria_singola(x0_45,y0_45,z0_45,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)

font = {'family':'serif','color':'darkred','weight':'normal','size': 10.7}
fig = plt.figure(num=None, figsize=(8, 6), dpi=800, facecolor='w', edgecolor='k')
ax = Axes3D(fig)
ax.set_xlabel('quota salari $\omega$')
ax.set_ylabel('occupazione $\lambda$')
ax.set_zlabel('quota debito privato $d$')

x_altro_flusso, y_altro_flusso, z_altro_flusso = flusso_traiettoria_singola(x0_altro,y0_altro,z0_altro,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
#ax.plot_wireframe(x_altro_flusso, y_altro_flusso, z_altro_flusso, color='purple')

#ax.plot_wireframe(x_45_flusso, y_45_flusso, z_45_flusso, color='green')
# PER MOSTRARE SOLO PARTE INIZIALE DEL FLUSSO
arrivare = 13000
ax.plot_wireframe(x_45_flusso[0:arrivare], y_45_flusso[0:arrivare], z_45_flusso[0:arrivare], color='green')

# 8.62402480510675; 11.7241196652734; 17.4452497910951 sono i tre punti (0,0,d_0) per r=r_1lambdanullomeno
#ax.scatter(0, 0, 8.62402480510675, 'o', s=30, color='blue', depthshade=True)
#ax.scatter(0, 0, 11.7241196652734, 'o', s=30, color='royalblue', depthshade=True)
#ax.scatter(0, 0, 17.4452497910951, 'o', s=30, color='blue', depthshade=True)
# 13.2796458147542; 22.4204384888053; 35.3726429691678 sono i tre punti (0,0,d_0) per r=z=2.5%
#ax.scatter(0, 0, 13.2796458147542, 'o', s=30, color='blue', depthshade=True)
#ax.scatter(0, 0, 22.4204384888053, 'o', s=30, color='royalblue', depthshade=True)
#ax.scatter(0, 0, 35.3726429691678, 'o', s=30, color='blue', depthshade=True)
#ax.scatter(x0_45, y0_45, z0_45, 'o', s=20, color='green', depthshade=True)
#ax.scatter(omegaeq1, lanbdaeq1, deq1, 'o',s=30, c='green',depthshade=True) # primq era red
#ax.plot_wireframe([0, 0], [0, 1], [(pi1meno - ni*(alfa+beta))/(r_1lambdanullomeno - (alfa+beta)), (pi1meno - ni*(alfa+beta))/(r_1lambdanullomeno - (alfa+beta))], lw=2, color='red') # i punti fissi instabili della soluzione (0,lambda,d_1)
#ax.quiver(X1, Y1, Z1, DX1, DY1, DZ1, color='deepskyblue', length=0.1, pivot='middle')

#ax.set_xlim(-0.5,2)
#ax.set_ylim(-0.5,1.5)
#ax.set_zlim(-20,40)


ip=350
#ip=numero_di_punti_tempo_voluti-12000


# E' diversa dalla funzione arrow 2D standard, quella usava x,y e dx, dy questa usa il punto iniziale della freccia (x[i],y[i],z[i]) e quello finale (x[i+1],y[i+1],z[i+1]) inoltre li vuole a coppie come vettori cioè ciascuno così [..,..] e non di seguito arrow(x[i],y[i],z[i],x[i+1],y[i+1],z[i+1]) come lo stile dell'altra
#a = Arrow3D([x_45_flusso[ip],x_45_flusso[ip+1]-x_45_flusso[ip]],[y_45_flusso[ip],y_45_flusso[ip+1]-y_45_flusso[ip]],[z_45_flusso[ip],z_45_flusso[ip+1]-z_45_flusso[ip]], mutation_scale=2, lw=1, arrowstyle="->",color="r")
a = Arrow3D([x_45_flusso[ip],x_45_flusso[ip+1]],[y_45_flusso[ip],y_45_flusso[ip+1]],[z_45_flusso[ip],z_45_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green") # arrowstyle="-|>" per fare la freccia piena
ax.add_artist(a) # serve per disegnare la freccia
#ax.set_aspect("equal") # serve ad usare un'unica scala vedere se usare o meno, far la prova se piace o no: se attiva schiacci la figura
ax.set_title('$r=$''%g' % zeta, fontdict = font)
#plt.savefig('relazione corretta Minsky alto r Figure?.jpg', dpi=850, transparent=False)
#ax.legend()
plt.show()

fig, ax1 = plt.subplots()
ax1.plot(t, x_45_flusso,'r',label='quota salari $\omega(t)$')
ax1.plot(t, y_45_flusso,'k',label='occupazione $\lambda(t)$')
#plt.legend(loc='lower right', prop={'size': 6})
plt.title('$r=$''%g' % zeta, fontdict = font)
#plt.title('grafico doppia scala con debito divergente a destra',fontdict = font)
ax1.set_xlabel('tempo t')
ax1.set_ylabel('$\omega, \lambda$')
plt.grid()
#ax1.set_ylim(0,1.2)  # POI RIMETTERE
# scala del debito privato divergente creando una nuova figura
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
colore = 'blue'
ax2.plot(t, z_45_flusso, color=colore,label='quota debito privato $d(t)$')
#ax2.legend(loc='center',prop={'size': 6})
ax1.plot(np.nan, colore,label='quota debito $d(t)$')  # Creo un falso grafico con np.nan solo per potere collocare la legenda del debito privato nella prima, infatti è ax1
ax1.legend(loc='best',prop={'size': 8})
ax1.set_xlim(0,700)
ax2.set_ylabel('$d$', color=colore)  # we already handled the x-label with ax1
ax2.set_ylim(-100,800)  # POI RIMETTERE
ax2.tick_params(axis='y', labelcolor=colore)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.savefig('relazione corretta Minsky alto r andamento temporale Figure?.jpg', dpi=650, transparent=False)
plt.show()
#%%
#%%
# ANDAMENTO TEMPORALE PARAGONE DIVERSE DINAMICHE COME PARE VOLERE TURCHETTI
# Forse aggiungere legenda col punto iniziale o informazioni ad esempio "dinamica sella"

colore = 'blue'

fig, ax1 = plt.subplots()
ax1.plot(t, x_7_flusso,'r',label='quota salari $\omega(t)$')
ax1.plot(t, y_7_flusso,'k',label='occupazione $\lambda(t)$')
#plt.legend(loc='lower right', prop={'size': 6})
#plt.title('converge al punto fisso stabile che è uno dei 3 (0,0,d0)',fontdict = font)
ax1.set_xlabel('tempo t')
ax1.set_ylabel('$\omega, \lambda$')
plt.grid()
ax1.set_ylim(-0.5,2)
# scala del debito privato divergente creando una nuova figura
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(t, z_7_flusso, color=colore,label='quota debito privato $d(t)$')
#ax2.legend(loc='center',prop={'size': 6})
ax1.plot(np.nan, colore,label='quota debito $d(t)$')  # Creo un falso grafico con np.nan solo per potere collocare la legenda del debito privato nella prima, infatti è ax1
ax1.legend(loc='best',prop={'size': 8})
ax1.set_xlim(0,200)
ax2.set_ylabel('$d$', color=colore)  # we already handled the x-label with ax1
ax2.set_ylim(-10,20)
ax2.tick_params(axis='y', labelcolor=colore)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.savefig('relazione corretta Minsky alto r con retta lambda andamento temporale 1.jpg', dpi=650, transparent=False)
plt.show()

fig, ax1 = plt.subplots()
ax1.plot(t, x_10_flusso,'r',label='quota salari $\omega(t)$')
ax1.plot(t, y_10_flusso,'k',label='occupazione $\lambda(t)$')
#plt.legend(loc='lower right', prop={'size': 6})
#plt.title('vortice a +oo, quello scelto da STEVE KEEN',fontdict = font)
ax1.set_xlabel('tempo t')
ax1.set_ylabel('$\omega, \lambda$')
plt.grid()
ax1.set_ylim(0,1.5)
# scala del debito privato divergente creando una nuova figura
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(t, z_10_flusso, color=colore,label='quota debito privato $d(t)$')
#ax2.legend(loc='center',prop={'size': 6})
ax1.plot(np.nan, colore,label='quota debito $d(t)$')  # Creo un falso grafico con np.nan solo per potere collocare la legenda del debito privato nella prima, infatti è ax1
ax1.legend(loc='best',prop={'size': 8})
ax1.set_xlim(0,200)
ax2.set_ylabel('$d$', color=colore)  # we already handled the x-label with ax1
ax2.set_ylim(0,250)
ax2.tick_params(axis='y', labelcolor=colore)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.savefig('relazione corretta Minsky alto r con retta lambda andamento temporale 2.jpg', dpi=650, transparent=False)
plt.show()

fig, ax1 = plt.subplots()
ax1.plot(t, x_20_flusso,'r',label='quota salari $\omega(t)$')
ax1.plot(t, y_20_flusso,'k',label='occupazione $\lambda(t)$')
#plt.legend(loc='lower right', prop={'size': 6})
#plt.title('ghirigoro poi va a -oo',fontdict = font)
ax1.set_xlabel('tempo t')
ax1.set_ylabel('$\omega, \lambda$')
plt.grid()
ax1.set_ylim(0,2)
# scala del debito privato divergente creando una nuova figura
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(t, z_20_flusso, color=colore,label='quota debito privato $d(t)$')
#ax2.legend(loc='center',prop={'size': 6})
ax1.plot(np.nan, colore,label='quota debito $d(t)$')  # Creo un falso grafico con np.nan solo per potere collocare la legenda del debito privato nella prima, infatti è ax1
ax1.legend(loc='best',prop={'size': 8})
ax1.set_xlim(0,150)
ax2.set_ylabel('$d$', color=colore)  # we already handled the x-label with ax1
ax2.set_ylim(-250,10)
ax2.tick_params(axis='y', labelcolor=colore)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.savefig('relazione corretta Minsky alto r con retta lambda andamento temporale 3.jpg', dpi=650, transparent=False)
plt.show()

fig, ax1 = plt.subplots()
ax1.plot(t, x_16_flusso,'r',label='quota salari $\omega(t)$')
ax1.plot(t, y_16_flusso,'k',label='occupazione $\lambda(t)$')
#plt.legend(loc='lower right', prop={'size': 6})
#plt.title('parte da quello che era il punto fisso 1 stabile e va a +oo',fontdict = font)
ax1.set_xlabel('tempo t')
ax1.set_ylabel('$\omega, \lambda$')
plt.grid()
ax1.set_ylim(0,1)
# scala del debito privato divergente creando una nuova figura
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(t, z_16_flusso, color=colore,label='quota debito privato $d(t)$')
#ax2.legend(loc='center',prop={'size': 6})
ax1.plot(np.nan, colore,label='quota debito $d(t)$')  # Creo un falso grafico con np.nan solo per potere collocare la legenda del debito privato nella prima, infatti è ax1
ax1.legend(loc='best',prop={'size': 8})
ax1.set_xlim(0,100)
ax2.set_ylabel('$d$', color=colore)  # we already handled the x-label with ax1
ax2.set_ylim(0,400)
ax2.tick_params(axis='y', labelcolor=colore)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.savefig('relazione corretta Minsky alto r con retta lambda andamento temporale 4.jpg', dpi=650, transparent=False)
plt.show()
# IN TOTALE QUINDI NELLA RELAZIONE NE METTI 4 DI GRAFICI TEMPORALI
#%%
#%%
# DISEGNO NULLCLINE e GRAFICO GLI ZERI DELLE FUNZIONI 

# nullcline e zeri della funzione f(omega) che compare in dlambda/dt come vuole Turchetti
assex=np.linspace(-1,2,1000)

#    domegadt = omega * (functionalform(A,B,C,D,lanbda) - alfa)
#    dlanbdadt = lanbda * (((functionalform(E,F,G,H,((1-omega)/ni)))/ni) - alfa - beta - gamma)

#omegaeq = 1 - ((ni*F)/G) +o- (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) )   # è un centro
#lanbdaeq = (B/C) -o+ ( np.sqrt(A/(alfa + D)) )/C   # è un centro, per cui functionalform(A,B,C,D,lanbda0)=alfa anche se dopo tanti zeri differisce per cui potrebbe non andare bene

# Grafico pi - f(pi)
plt.plot(assex, ((functionalform(E,F,G,H,((assex)/ni)))/ni) - alfa - beta - gamma,'blue')
plt.plot([0,0],[-10,20],'k')
plt.plot([-10,20],[0,0],'k')
plt.plot([(ni*F)/G,(ni*F)/G],[-10,20],'k')
plt.plot(((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),0,'ob')
plt.plot(((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),0,'ob')
plt.legend(loc='best',prop={'size': 6})
#plt.title('grafico pi f(pi) NON è la nullcline di $\dot{\lambda}$',fontdict = font)
plt.xlabel('$\pi$')
plt.ylabel('$f(\pi)$')
plt.grid()
plt.xlim(-0.2,1.1)
plt.ylim(-0.4,2.2)
#plt.savefig('relazione corretta Minsky grafico pi f(pi).jpg', dpi=650, transparent=False)
plt.show()

# Grafico lambda - g(lambda)
plt.plot(assex, functionalform(A,B,C,D,assex) - alfa,'red')
plt.plot([0,0],[-1,3],'k')
plt.plot([-1,3],[0,0],'k')
plt.plot([1,1],[-1,3],'k')
plt.plot(lanbdaeq1,0,'or')
plt.legend(loc='best',prop={'size': 6})
#plt.title('grafico lambda g(lambda) NON è la nullcline di $\dot{\omega}$',fontdict = font)
plt.xlabel('$\lambda$')
plt.ylabel('$g(\lambda)$')
plt.grid()
plt.xlim(-0.1,1.2)
plt.ylim(-0.2,1.2)
#plt.savefig('relazione corretta Minsky grafico lambda g(lambda).jpg', dpi=650, transparent=False)
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
#plt.savefig('relazione corretta Minsky grafico lambda w[lambda].jpg', dpi=650, transparent=False)
plt.show()

# Grafico pi - k[pi/ni]
plt.plot(assex, functionalform(E,F,G,H,(assex/ni)),'red')
plt.plot([0,0],[-1,3],'k')
plt.plot([-1,3],[0,0],'k')
plt.plot([(ni*F)/G,(ni*F)/G],[-10,20],'k')
plt.plot(((ni*F)/G) + (ni/G)*( np.sqrt(E/H) ),0,'or')
plt.plot(((ni*F)/G) - (ni/G)*( np.sqrt(E/H) ),0,'or')
plt.legend(loc='best',prop={'size': 6})
#plt.title('funzione risposta capitalisti grafico pi k[pi/ni]',fontdict = font)
plt.xlabel('$\pi$')
plt.ylabel('$k[\pi / \\nu]$')
plt.grid()
plt.xlim(-0.1,1.2)
plt.ylim(-0.2,1.2)
#plt.savefig('relazione corretta Minsky grafico omega k[pi/nu].jpg', dpi=650, transparent=False)
plt.show()

# Faccio vedere a scopo dimostrativo che il campo vettoriale è perpendicolare alle nullcline vicino ad esse
# Draw Nullclines and Quiver plot
# plot nullclines
plt.plot([-5,5],[(B/C) - ( np.sqrt(A/(alfa + D)) )/C ,(B/C) - ( np.sqrt(A/(alfa + D)) )/C ], 'b-', lw=2, label='$\dot{\omega}$ nullcline')
plt.plot([-5,5],[(B/C) + ( np.sqrt(A/(alfa + D)) )/C ,(B/C) + ( np.sqrt(A/(alfa + D)) )/C ], 'b-', lw=2)
plt.plot([((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ), ((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) )],[-5,5], 'r-', lw=2, label='$\dot{\lambda}$ nullcline')
plt.plot([((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ), ((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) )],[-5,5], 'r-', lw=2)
plt.plot(functionalform(A,B,C,D,assex) - alfa,assex,'green') # questa in realtà andrebbe tolta perchè N.B.: g(lambda) NON è la nullcline, essa sono i due valori di lambda* infatti la nullcline di lambda non dipende da omega, è un retta orizzontale quindi
plt.plot(assex, ((functionalform(E,F,G,H,(assex/ni)))/ni) - alfa - beta - gamma,'green') #  questa in realtà andrebbe tolta perchè N.B.: f(omega) NON è la nullcline, essa è il valore di omega* infatti la nullcline di omega non dipende da lambda, è un retta verticale quindi
plt.plot([0,0],[-5,5],'k',lw=2)
plt.plot([-5,5],[0,0],'k',lw=2)
plt.plot([-5,5],[1,1],'--k')
plt.plot([(ni*F)/G,(ni*F)/G],[-5,5],'--k')
plt.legend(loc='best',prop={'size': 9})
plt.xlabel('$\pi$')
plt.ylabel('$\lambda$')
plt.xlim(-0.1,1.2)
plt.ylim(-0.1,1.2)
# plot fixed points
# qua non funziona più il metodo di trovare gli equilibri automaticamente per cui devi fare manualmente
#for point in equilibria:
#    plt.plot(point[0],point[1],"blue", marker = "o", markersize = 6.0)
plt.plot(((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),0,'or')
plt.plot(((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),0,'or')
plt.plot(0,(B/C) - (1/C)*( np.sqrt(A/(D + (alfa ))) ),'ob')
plt.plot(0,(B/C) + (1/C)*( np.sqrt(A/(D + (alfa ))) ),'ob')
plt.title("Quiverplot with nullclines",fontdict = font)
#plt.quiver(X1, Y1, DX1, DY1, color='deepskyblue',pivot='mid')
plt.grid()
# carina da tenere che si vede meglio anche se NON LA DEVI METTERE NELLA RELAZIONE
#plt.savefig('relazione corretta Minsky nullcline spazio delle fasi.jpg', dpi=650, transparent=False)
plt.show()

# Grafico omega,lambda - f(omega),g(lambda)
plt.plot(assex,((functionalform(E,F,G,H,(assex/ni)))/ni) - alfa - beta - gamma,'blue',label='$f(\pi)$')
plt.plot(assex,functionalform(A,B,C,D,assex)-alfa,'red', label='$g(\lambda)$')
plt.plot([-1,2],[0,0],'k', [0,0],[-1,2],'k',[((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) ),lanbdaeq1],[0,0,0],'or')
# oppure tutto insieme ma non funziona il comando label non so perchè!
#plt.plot(assex,functionalform(A,B,C,D,assex)-alfa,'red', assex,(((1-assex)/ni) - alfa - beta - gamma),'blue', [-1,2],[0,0],'k', [0,0],[-1,2],'k',[omega0eq1,lanbda0eq1],[0,0],'or')
plt.legend(loc='upper center',prop={'size': 8})
plt.plot([(ni*F)/G,(ni*F)/G],[-5,5],'--k')
plt.plot([1,1],[-5,5],'--k')
plt.axis([-0.2, 1.2, -0.2, 1.3])
plt.xlabel('$\pi, \lambda$')
plt.ylabel('$f(\pi), g(\lambda)$')
#plt.xlim(-1,2)
#plt.ylim(-1,1)
plt.grid()
plt.show()
#%%
#%%
# AUTOVALORI JACOBIANA SISTEMA 3 EQUAZIONI MINSKY con r=zeta dove N.B.: zeta è quello che stato dichiarato globalmente all'inizio del programma, quindi stai attento a cosa compili!

def jacobiana_sistema_finanziario_3_eq_minsky(omegaeq,lanbdaeq,deq):
    r = zeta
    pi = 1 - omegaeq - r*deq
    A11 = functionalform(A,B,C,D,lanbdaeq) - alfa
    A12 = omegaeq * ( (2*A*C)/(B-C*lanbdaeq)**3 )
    A13 = 0
    A21 = (-(lanbdaeq/ni)) * ( ((2*E*G)/ni)/(F-G*((pi)/ni))**3 )
    A22 = (functionalform(E,F,G,H,((pi)/ni))/ni) - alfa - beta - gamma
    A23 = (-(lanbdaeq/ni)) * ( r*((2*E*G)/ni)/(F-G*((pi)/ni))**3 )
    A31 = 1 - ((ni-deq)/ni) * ( ((2*E*G)/ni)/(F-G*((pi)/ni))**3 )
    A32 = 0
    A33 = 2*r - ( (functionalform(E,F,G,H,((pi)/ni))/ni) - gamma ) - ((ni-deq)/ni) * ( r*((2*E*G)/ni)/(F-G*((pi)/ni))**3 )
    #A11=round(A11,9) # sennò non viene zero perchè tiene 17 cifre decimali!
    #A12=round(A12,9)
    #A21=round(A21,9)
    #A22=round(A22,9)
    #A22=round(A22) # perchè sennò viene -0.0 quando omegaeq è quello che annulla "quarto". N.B.: TTENZIONE CHE IN TUTTI GLI ALTRI CASI TI FA VENIRE RISULTATI SBAGLIATI!!! Meglio lasciare sempre commentato, non usare
    return np.array([[A11, A12, A13],[A21, A22, A23],[A31, A32, A33]])

def studia_i_vari_punti_fissi(x,y,z):
    J = jacobiana_sistema_finanziario_3_eq_minsky(x,y,z)
    traccia_di_J = np.trace(J)
    determinante_di_J = np.linalg.det(J)
    autovalori_di_J, autovettori_di_J = LA.eig(J)
    print('Sto analizzando il punto fisso  %s, %s, %s' % (x,y,z)) 
    print('\nLa matrice Jacobiana del modello finanziario Minsky a 3 equazioni di Steve Keen è:\n',J,'\n')
    print('Il determinante det(J) è:\n', determinante_di_J)
    print('La traccia tr(J) è:\n', traccia_di_J)
    print('Gli autovalori della matrice J sono:\n', autovalori_di_J)
    print('Una possibile scelta per i 3 autovettori di J è:\n',autovettori_di_J)
    print('----------------------------')
    print('\n')
    return autovalori_di_J, determinante_di_J, traccia_di_J
#%%
# Indago stabilità punti fissi per il caso r GRANDE con retta punti (0,lambda,d_1) cioè per r_1lambdanullomeno

debito0annullaterzaeq1 = 8.62402480510675    # se fossero diversi poichè hai cambiato zeta e non stai più usando r_1lambdanullomeno,
debito0annullaterzaeq2 = 11.7241196652734    # dovrai ovviamente inserire manualmente quelli giusti
debito0annullaterzaeq3 = 17.4452497910951    # Infatti sotto ho messo il controllo che l'r usato sia quello della retta (0,lambda,d_1)

# inserisco manualmente i punti fissi in una matrice dove ogni riga è un punto fisso
matrice_punti_fissi = np.array([ [0,0,debito0annullaterzaeq1],[0,0,debito0annullaterzaeq2],[0,0,debito0annullaterzaeq3],[omegaeq1,lanbdaeq1,deq1],[0,0.4,d_1lambdanullomeno],[0,0.8,d_1lambdanullomeno] ])
# metto dei separatori per far capire dove inizia la stampa dei risultati
print('\n')
print('----------------------------')
print('----------------------------')
print('----------------------------')
print('\n')
# itero su tutti i punti fissi inseriti manualmente
for x,y,z in matrice_punti_fissi:
    studia_i_vari_punti_fissi(x,y,z)    
#%%    
#%%
# Indago stabilità punti fissi per il caso r PICCOLO cioè r=2.5% e NON c'è la retta punti (0,lambda,d_1). r NON è r_1lambdanullomeno

debito0annullaterzaeq1 = 13.2796458147542  
debito0annullaterzaeq2 = 22.4204384888053
debito0annullaterzaeq3 = 35.3726429691678

# inserisco manualmente i punti fissi in una matrice dove ogni riga è un punto fisso
matrice_punti_fissi = np.array([ [0,0,debito0annullaterzaeq1],[0,0,debito0annullaterzaeq2],[0,0,debito0annullaterzaeq3],[omegaeq1,lanbdaeq1,deq1] ])
# metto dei separatori per far capire dove inizia la stampa dei risultati
print('\n')
print('----------------------------')
print('----------------------------')
print('----------------------------')
print('\n')
# itero su tutti i punti fissi inseriti manualmente
for x,y,z in matrice_punti_fissi:
    studia_i_vari_punti_fissi(x,y,z)
#%%
#%%
#%%
#%%
#%%
#%%
#%%
#%%
#                                     ALTRO CASO. GRAFICI r PICCOLO (r=zeta=0.025)
#%%
# TESTO CHE i d_0 TROVATI (per il valore di r, o zeta, scelto sopra) DEL PUNTO FISSO (0,0,d_0) SIANO EFFETTIVAMENTE SOLUZIONI della TERZA EQUAZIONE dd/dt=0

# inserisco manualmente i punti fissi equilibria_d_0 calcolati per r = zeta = 0.025
debito0annullaterzaeq1 = 13.2796458147542    # se fossero diversi poichè hai cambiato zeta e non stai più usando zeta=0.025,
debito0annullaterzaeq2 = 22.4204384888053    # dovrai ovviamente inserire manualmente quelli giusti
debito0annullaterzaeq3 = 35.3726429691678    # Infatti sotto ho messo il controllo che l'r usato NON sia quello della retta (0,lambda,d_1)
debito0annullaterzaeq_vettore = np.array([debito0annullaterzaeq1,debito0annullaterzaeq2,debito0annullaterzaeq3])
terzaeqsiannullaperd0_punto_di_domanda_vettore = 2*zeta*debito0annullaterzaeq_vettore - 1 + (ni - debito0annullaterzaeq_vettore) * (((functionalform(E,F,G,H,((1 - zeta*debito0annullaterzaeq_vettore)/ni)))/ni) - gamma)
# un risultato del tipo -1.35724765e-14 significa -1.347*10^{-14} cioè in pratica è zero. Nella terminologia di Python "e" sighifica "moltiplicato per 10^{}"
2.5/100==zeta
#%%
#%%
# FACCIO I GRAFICI PER IL CASO r=2.5%, NON ESISTE LA RETTA (0,lambda,d_1) DI PUNTI FISSI
#%%
#%%       
# INTEGRAZIONE
# definisco le condizioni iniziali di ogni flusso ma è ridondante perchè potrei direttamente farlo nella chiamata di flusso_traiettoria_singola(x_0,y_0,..,..) mettendo direttamente i numeri

# converge al secondo (0,0,d_0) che è stabile partendoci da sopra. Ci arriva anche da (0.2,0.2,22.4204384888053+3)
x0_1 = 0 # 0.2 e 0.2 ci arriva
y0_1 = 0.6
z0_1 = 22.4204384888053+6.5
# non è quello di Keen, fa vortice da sotto il punto fisso stabile e ci arriva. Disegnare tutto flusso
#x0_1 = 0.9
#y0_1 = 0.96
#z0_1 = 0
# fa vortice e converge a punto fisso stabile, fa un vortice da sotto che si bacia con quello di Keen. Disegnare tutto flusso
x0_2 = 0.8
y0_2 = 0.8
z0_2 = 0
# fa vortice e converge a punto fisso stabile. Disegnare tutto flusso
# per omega iniziale appena superiore di 1.3 va a debito +oo, quindi questo è l'ultimo stabile
x0_3 = 1.3
y0_3 = 0.9
z0_3 = 0

# INFO UTILI sulla dinamica di altri flussi che divergono però sballano la scala di omega e di lambda perchè li fanno crescere enormemente. Nessuna interpretazione economica possibile! 
# (0.7,0.7,0) il debito va a -oo ma NON LO GRAFICO perchè lambda va +oo e sballa la scala del grafico
# (0.75,0.75,0) il debito va a +oo ma NON LO GRAFICO perchè omega va +oo e sballa la scala del grafico
# (0.76,0.76,0) CURIOSAMENTE anche se è solo di 0.01 più grande del precedente, si comporta come quello (0.7,0.7,0):il debito va a -oo ma NON LO GRAFICO perchè lambda va +oo e sballa la scala del grafico

# va a -oo
x0_4 = 0.5
y0_4 = 0.5
z0_4 = 0.5
# QUELLO SCELTO DA STEVE KEEN. Fa vortice da sopra il punto fisso stabile e ci converge 
x0_5 = 0.96
y0_5 = 0.9
z0_5 = 0
# va a -oo come tutti i punti (0.#,0.#,0.#) da qua fino a 0
x0_6 = 0.6
y0_6 = 0.6
z0_6 = 0.6
# va a -oo
x0_7 = 0.4
y0_7 = 0.4 
z0_7 = 0.4
# va a -oo
x0_8 = 0.3
y0_8 = 0.3
z0_8 = 0.3
# va a -oo
x0_9 = 0.2
y0_9 = 0.2
z0_9 = 0.2

#[13.2796458147542, 22.4204384888053, 35.3726429691678]

# va a +oo. Mostrare flusso fino a 650
x0_10 = 1
y0_10 = 0.6
z0_10 = 13.2796458147542
# parte sopra al terzo dei punti (0,0,d_0) che è instabile e va a +oo comportandosi bene
x0_11 = 0
y0_11 = 0
z0_11 = 35.3726429691678 + 0.1
# parte a destra del terzo dei punti (0,0,d_0) che è instabile e va a +oo
x0_12 = 1
y0_12 = 0
z0_12 = 35.3726429691678
# va a +oo popola un po' l'area
x0_13 = 0.7
y0_13 = 0.7
z0_13 = 35.3726429691678 - 15
# parte a destra all'altezza del secondo dei punti (0,0,d_0) che siccome è fisso pensavo convergesse ad esso invece va a +oo bene
x0_14 = 1
y0_14 = 0
z0_14 = 22.4204384888053
# va a -oo partendo da vicino al secondo punto dei (0,0,d_0) che è stabile ma non lo attrae. Mostrare flusso fino a 2100, è molto lento
x0_15 = 0.3
y0_15 = 0
z0_15 = 22.4204384888053
# va al secondo dei punti (0,0,d_0) che è stabile.
x0_16 = 0.28
y0_16 = 0
z0_16 = 22.4204384888053
# questo invece va al secondo punto dei (0,0,d_0) che è stabile.
#CURIOSAMENTE basta partire da omega=0.29 o da lambda=0.4 che non converge al punto fisso stabile!!!
x0_17 = 0.28
y0_17 = 0.3
z0_17 = 22.4204384888053
# va a -oo comportandosi bene, curiosamente partendo da lambda appena più piccolo cioè 0.3 sarebbe andato al punto fisso stabile invece
x0_18 = 28
y0_18 = 0.4
z0_18 = 22.4204384888053
# va a -oo. Mostrare flusso fino a 1700
x0_19 = 0.3
y0_19 = 0.3
z0_19 = 13.2796458147542
# va a +oo. Mostrare flusso fino a 650
x0_20 = 1
y0_20 = 0.4
z0_20 = 13.2796458147542


#dt = 0.01
#numero_di_punti_tempo_voluti = 5000  # deve essere un numero tale che tempo_finale venga un numero intero
#dt = 0.01
dt=0.03 # QUELLO USATO! GIUSTO!
#dt=0.005
#dt=0.004 # o forse meglio 0.002
#numero_di_punti_tempo_voluti = 50000
#numero_di_punti_tempo_voluti = 25000
numero_di_punti_tempo_voluti = 15000
tempo_finale = numero_di_punti_tempo_voluti * dt
t = np.linspace(0,tempo_finale,numero_di_punti_tempo_voluti) # per fare il grafico dopo, l'ultimo è proprio tempo_finale

# ovviamente tutto ciò si poteva saltare mettendo direttamente i numeri al posto di x0_i e y0_i definendoli qui e non sopra che è ridondante
x_1_flusso, y_1_flusso, z_1_flusso = flusso_traiettoria_singola(x0_1,y0_1,z0_1,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_2_flusso, y_2_flusso, z_2_flusso = flusso_traiettoria_singola(x0_2,y0_2,z0_2,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_3_flusso, y_3_flusso, z_3_flusso = flusso_traiettoria_singola(x0_3,y0_3,z0_3,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_4_flusso, y_4_flusso, z_4_flusso = flusso_traiettoria_singola(x0_4,y0_4,z0_4,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_5_flusso, y_5_flusso, z_5_flusso = flusso_traiettoria_singola(x0_5,y0_5,z0_5,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_6_flusso, y_6_flusso, z_6_flusso = flusso_traiettoria_singola(x0_6,y0_6,z0_6,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_7_flusso, y_7_flusso, z_7_flusso = flusso_traiettoria_singola(x0_7,y0_7,z0_7,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_8_flusso, y_8_flusso, z_8_flusso = flusso_traiettoria_singola(x0_8,y0_8,z0_8,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_9_flusso, y_9_flusso, z_9_flusso = flusso_traiettoria_singola(x0_9,y0_9,z0_9,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_10_flusso, y_10_flusso, z_10_flusso = flusso_traiettoria_singola(x0_10,y0_10,z0_10,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_11_flusso, y_11_flusso, z_11_flusso = flusso_traiettoria_singola(x0_11,y0_11,z0_11,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_12_flusso, y_12_flusso, z_12_flusso = flusso_traiettoria_singola(x0_12,y0_12,z0_12,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_13_flusso, y_13_flusso, z_13_flusso = flusso_traiettoria_singola(x0_13,y0_13,z0_13,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_14_flusso, y_14_flusso, z_14_flusso = flusso_traiettoria_singola(x0_14,y0_14,z0_14,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_15_flusso, y_15_flusso, z_15_flusso = flusso_traiettoria_singola(x0_15,y0_15,z0_15,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_16_flusso, y_16_flusso, z_16_flusso = flusso_traiettoria_singola(x0_16,y0_16,z0_16,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_17_flusso, y_17_flusso, z_17_flusso = flusso_traiettoria_singola(x0_17,y0_17,z0_17,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_18_flusso, y_18_flusso, z_18_flusso = flusso_traiettoria_singola(x0_18,y0_18,z0_18,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_19_flusso, y_19_flusso, z_19_flusso = flusso_traiettoria_singola(x0_19,y0_19,z0_19,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
x_20_flusso, y_20_flusso, z_20_flusso = flusso_traiettoria_singola(x0_20,y0_20,z0_20,dt,numero_di_punti_tempo_voluti,sistema_finanziario_3_eq_minsky)
#%%
# GRAFICO SPAZIO DELLE FASI
font = {'family':'serif','color':'darkred','weight':'normal','size': 10.7,}
# plt.text(2, 0.65, r'$\cos(2 \pi t) \exp(-t)$', fontdict=font)


# FRECCE DEL CAMPO VETTORIALE

# creo la griglia, in realtà saranno 2 matrici una con coordinata x e l'altra y così prese elemento per elemento danno il punto (x_i,y_i)
x_griglia = np.linspace(0, 2, 6) # buono anche solo 14 però poi le frecce non sono tutte ortogonali come dovrebbe invece essere
y_griglia = np.linspace(0, 2, 6)
z_griglia = np.linspace(0, 2, 6)

# il comando sotto è interessante perchè assegna automaticamente i termini in ordine: es. gaia=np.array([3,4]) poi con hippie1,hippie2=gaia si ha hippie1=gaia[0]=3 e hippie2=gaia[1]=4
X1 , Y1, Z1  = np.meshgrid(x_griglia, y_griglia, z_griglia)    # create a grid
# ritorna una prima matrice che ha come righe il vettore x ripetuto tante volte quante la lunghezza del vettore y, cioè len(y) volte, e una seconda matrice che ha il vettore y come colonne ripetuto tante volte quanta la lunghezza del vettore x, len(x)
DX1, DY1, DZ1 = sistema_finanziario_3_eq_minsky([X1, Y1, Z1])    # compute growth rate on the grid
# la moltiplicazione * usata in "sistema_II_semplificato" è quella normale e se X come in questo caso è composto da 2 matrici, quando spacchetta X assegna la prima matrice (cioè le coordinate x su tutta la griglia) a x e la seconda (cioè la griglia per le coordinate y) a y, quando poi si trova a fare x*y moltiplica elemento per elemento entrambe le matrici senza fare il prodotto matriciale! Perciò le due matrici devono avere entrambe la stessa dimensione, tutte e due devono essere nxn
M = (np.hypot(DX1, DY1, DZ1))                        # norm growth rate 
M[ M == 0] = 1.                                 # avoid zero division errors 
DX1 /= M                                        # normalize each arrows
DY1 /= M
# se fai M.shape dà (20,20) è infatti una matrice 20x20 siccome hai usato due volte x = np.linspace(0, 2, 20) o x e y, ma sono uguali x=y=np.linspace(0, 2, 20)


fig = plt.figure(num=None, figsize=(8, 6), dpi=800, facecolor='w', edgecolor='k')
ax = Axes3D(fig)


# ATTIVARE SE VUOI VEDERE LE FRECCE DEL CAMPO VETTORIALE ma non penso che Turchetti le voglia
#ax.quiver(X1, Y1, Z1, DX1, DY1, DZ1, color='deepskyblue', length=0.1, pivot='middle') # anche senza length=0.1 magari
#plt.quiver(X1, Y1, Z1, DX1, DY1, DZ1, M, pivot='mid') # M gli dà una scala di colori che sono i numeri

# per non graficare tutti i punti sennò vanno a -oo e +oo 
mostrafinoa = 2600
mostrafinoa1 = 400 # per flussi che divergono subito BUONO PARE
mostrafinoa2 = 700 # per flussi che hanno bisogno di farsi vedere un poco di più
mostrafinoa3 = 3080 # per il flusso 14
mostrafinoa4 = 150 # per il flusso 11 che va a +oo pare
mostrafinoa5 = 230+200 # per il flusso 17 che va a +oo pare
mostrafinoa6 = 860 # per il flusso 13 che va a +oo pare

ax.plot_wireframe(x_1_flusso, y_1_flusso, z_1_flusso, color='green')
#ax.plot_wireframe(x_1_flusso[0:mostrafinoa1], y_1_flusso[0:mostrafinoa1], z_1_flusso[0:mostrafinoa1], color='green')
ax.plot_wireframe(x_2_flusso, y_2_flusso, z_2_flusso, color='green')
#ax.plot_wireframe(x_2_flusso[0:mostrafinoa1], y_2_flusso[0:mostrafinoa1], z_2_flusso[0:mostrafinoa1], color='green')
ax.plot_wireframe(x_3_flusso, y_3_flusso, z_3_flusso, color='green')
#ax.plot_wireframe(x_3_flusso[0:mostrafinoa1], y_3_flusso[0:mostrafinoa1], z_3_flusso[0:mostrafinoa1], color='green')
#ax.plot_wireframe(x_4_flusso, y_4_flusso, z_4_flusso, color='green')
ax.plot_wireframe(x_4_flusso[0:mostrafinoa1+100], y_4_flusso[0:mostrafinoa1+100], z_4_flusso[0:mostrafinoa1+100], color='green')
ax.plot_wireframe(x_5_flusso, y_5_flusso, z_5_flusso, color='green')
#ax.plot_wireframe(x_5_flusso[0:mostrafinoa], y_5_flusso[0:mostrafinoa], z_5_flusso[0:mostrafinoa], color='green')
#ax.plot_wireframe(x_6_flusso, y_6_flusso, z_6_flusso, color='green')
ax.plot_wireframe(x_6_flusso[0:mostrafinoa5], y_6_flusso[0:mostrafinoa5], z_6_flusso[0:mostrafinoa5], color='green')
#ax.plot_wireframe(x_7_flusso, y_7_flusso, z_7_flusso, color='green')
ax.plot_wireframe(x_7_flusso[0:mostrafinoa5], y_7_flusso[0:mostrafinoa5], z_7_flusso[0:mostrafinoa5], color='green')
#ax.plot_wireframe(x_8_flusso, y_8_flusso, z_8_flusso, color='green')
ax.plot_wireframe(x_8_flusso[0:mostrafinoa5], y_8_flusso[0:mostrafinoa5], z_8_flusso[0:mostrafinoa5], color='green')
#ax.plot_wireframe(x_9_flusso, y_9_flusso, z_9_flusso, color='green')
ax.plot_wireframe(x_9_flusso[0:mostrafinoa5], y_9_flusso[0:mostrafinoa5], z_9_flusso[0:mostrafinoa5], color='green')
#ax.plot_wireframe(x_10_flusso, y_10_flusso, z_10_flusso, color='green')
ax.plot_wireframe(x_10_flusso[0:mostrafinoa5], y_10_flusso[0:mostrafinoa5], z_10_flusso[0:mostrafinoa5], color='green')
#ax.plot_wireframe(x_11_flusso, y_11_flusso, z_11_flusso, color='green')
ax.plot_wireframe(x_11_flusso[0:mostrafinoa5-100], y_11_flusso[0:mostrafinoa5-100], z_11_flusso[0:mostrafinoa5-100], color='green')
#ax.plot_wireframe(x_12_flusso, y_12_flusso, z_12_flusso, color='green')
ax.plot_wireframe(x_12_flusso[0:mostrafinoa1-150], y_12_flusso[0:mostrafinoa1-150], z_12_flusso[0:mostrafinoa1-150], color='green')
#ax.plot_wireframe(x_13_flusso, y_13_flusso, z_13_flusso, color='green')
ax.plot_wireframe(x_13_flusso[0:mostrafinoa5], y_13_flusso[0:mostrafinoa5], z_13_flusso[0:mostrafinoa5], color='green')
#ax.plot_wireframe(x_14_flusso, y_14_flusso, z_14_flusso, color='green')
ax.plot_wireframe(x_14_flusso[0:mostrafinoa5], y_14_flusso[0:mostrafinoa5], z_14_flusso[0:mostrafinoa5], color='green')
#ax.plot_wireframe(x_15_flusso, y_15_flusso, z_15_flusso, color='green')
ax.plot_wireframe(x_15_flusso[0:mostrafinoa], y_15_flusso[0:mostrafinoa], z_15_flusso[0:mostrafinoa], color='green')
ax.plot_wireframe(x_16_flusso, y_16_flusso, z_16_flusso, color='green')
#ax.plot_wireframe(x_16_flusso[0:mostrafinoa1], y_16_flusso[0:mostrafinoa1], z_16_flusso[0:mostrafinoa1], color='green')
ax.plot_wireframe(x_17_flusso, y_17_flusso, z_17_flusso, color='green')
#ax.plot_wireframe(x_17_flusso[0:mostrafinoa5], y_17_flusso[0:mostrafinoa5], z_17_flusso[0:mostrafinoa5], color='green')
#ax.plot_wireframe(x_18_flusso, y_18_flusso, z_18_flusso, color='green')
ax.plot_wireframe(x_18_flusso[0:mostrafinoa5], y_18_flusso[0:mostrafinoa5], z_18_flusso[0:mostrafinoa5], color='green')
#ax.plot_wireframe(x_19_flusso, y_19_flusso, z_19_flusso, color='green')
ax.plot_wireframe(x_19_flusso[0:1700], y_19_flusso[0:1700], z_19_flusso[0:1700], color='green')
#ax.plot_wireframe(x_20_flusso, y_20_flusso, z_20_flusso, color='green')
ax.plot_wireframe(x_20_flusso[0:650], y_20_flusso[0:650], z_20_flusso[0:650], color='green')

# Nel grafico che sto creando metto anche tutti punti fissi. Ci sono 2 modi possibili per farlo
# versione "automatica" però non puoi variare i colori ma qua non funziona perchè non riesce a trovarli:
#for point in equilibria:
#    plt.plot(point[0],point[1],"royalblue", marker = "o", markersize = 7.0)

# versione "manuale":
ax.scatter(omegaeq1, lanbdaeq1, deq1, 'o',s=30, c='royalblue',depthshade=True)
ax.scatter([0, 0, 0], [0, 0, 0], [13.2796458147542, 22.4204384888053, 35.3726429691678], 'o',s=30, c='red',depthshade=True) # i punti fissi instabili della soluzione (0,0,d_0)

#ip=1300 # vedere se va bene
ip=280+20
ip1=100
ip2=300
ip3=1820
ip4=1000 # per il flusso 16
ip5=200
ip6=2300


# E' diversa dalla funzione arrow 2D standard, quella usava x,y e dx, dy questa usa il punto iniziale della freccia (x[i],y[i],z[i]) e quello finale (x[i+1],y[i+1],z[i+1]) inoltre li vuole a coppie come vettori cioè ciascuno così [..,..] e non di seguito arrow(x[i],y[i],z[i],x[i+1],y[i+1],z[i+1]) come lo stile dell'altra
a1 = Arrow3D([x_1_flusso[ip],x_1_flusso[ip+1]],[y_1_flusso[ip],y_1_flusso[ip+1]],[z_1_flusso[ip],z_1_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green") # arrowstyle="-|>" per fare la freccia piena
ax.add_artist(a1) # serve per disegnare la freccia
a2 = Arrow3D([x_2_flusso[ip],x_2_flusso[ip+1]],[y_2_flusso[ip],y_2_flusso[ip+1]],[z_2_flusso[ip],z_2_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green") # arrowstyle="-|>" per fare la freccia piena
ax.add_artist(a2) # serve per disegnare la freccia
a3 = Arrow3D([x_3_flusso[ip5],x_3_flusso[ip5+1]],[y_3_flusso[ip5],y_3_flusso[ip5+1]],[z_3_flusso[ip5],z_3_flusso[ip5+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green") # arrowstyle="-|>" per fare la freccia piena
ax.add_artist(a3) # serve per disegnare la freccia
a4 = Arrow3D([x_4_flusso[ip],x_4_flusso[ip+1]],[y_4_flusso[ip],y_4_flusso[ip+1]],[z_4_flusso[ip],z_4_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green") # arrowstyle="-|>" per fare la freccia piena
ax.add_artist(a4) # serve per disegnare la freccia
a5 = Arrow3D([x_5_flusso[ip],x_5_flusso[ip+1]],[y_5_flusso[ip],y_5_flusso[ip+1]],[z_5_flusso[ip],z_5_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green") # arrowstyle="-|>" per fare la freccia piena
ax.add_artist(a5) # serve per disegnare la freccia
a6 = Arrow3D([x_6_flusso[ip1],x_6_flusso[ip1+1]],[y_6_flusso[ip1],y_6_flusso[ip1+1]],[z_6_flusso[ip1],z_6_flusso[ip1+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a6)
a7 = Arrow3D([x_7_flusso[ip],x_7_flusso[ip+1]],[y_7_flusso[ip],y_7_flusso[ip+1]],[z_7_flusso[ip],z_7_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a7)
a8 = Arrow3D([x_8_flusso[ip],x_8_flusso[ip+1]],[y_8_flusso[ip],y_8_flusso[ip+1]],[z_8_flusso[ip],z_8_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a8)
a9 = Arrow3D([x_9_flusso[ip],x_9_flusso[ip+1]],[y_9_flusso[ip],y_9_flusso[ip+1]],[z_9_flusso[ip],z_9_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a9)
a10 = Arrow3D([x_10_flusso[ip],x_10_flusso[ip+1]],[y_10_flusso[ip],y_10_flusso[ip+1]],[z_10_flusso[ip],z_10_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a10)
a11 = Arrow3D([x_11_flusso[ip],x_11_flusso[ip+1]],[y_11_flusso[ip],y_11_flusso[ip+1]],[z_11_flusso[ip],z_11_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a11)
a12 = Arrow3D([x_12_flusso[ip1],x_12_flusso[ip1+1]],[y_12_flusso[ip1],y_12_flusso[ip1+1]],[z_12_flusso[ip1],z_12_flusso[ip1+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a12)
a13 = Arrow3D([x_13_flusso[ip],x_13_flusso[ip+1]],[y_13_flusso[ip],y_13_flusso[ip+1]],[z_13_flusso[ip],z_13_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a13)
a14 = Arrow3D([x_14_flusso[ip],x_14_flusso[ip+1]],[y_14_flusso[ip],y_14_flusso[ip+1]],[z_14_flusso[ip],z_14_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a14)
a15 = Arrow3D([x_15_flusso[ip6],x_15_flusso[ip6+1]],[y_15_flusso[ip6],y_15_flusso[ip6+1]],[z_15_flusso[ip6],z_15_flusso[ip6+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a15)
a16 = Arrow3D([x_16_flusso[ip4],x_16_flusso[ip4+1]],[y_16_flusso[ip4],y_16_flusso[ip4+1]],[z_16_flusso[ip4],z_16_flusso[ip4+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a16)
a17 = Arrow3D([x_17_flusso[ip5],x_17_flusso[ip5+1]],[y_17_flusso[ip5],y_17_flusso[ip5+1]],[z_17_flusso[ip5],z_17_flusso[ip5+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a17)
a18 = Arrow3D([x_18_flusso[ip],x_18_flusso[ip+1]],[y_18_flusso[ip],y_18_flusso[ip+1]],[z_18_flusso[ip],z_18_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a18)
a19 = Arrow3D([x_19_flusso[ip],x_19_flusso[ip+1]],[y_19_flusso[ip],y_19_flusso[ip+1]],[z_19_flusso[ip],z_19_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a19)
a20 = Arrow3D([x_20_flusso[ip],x_20_flusso[ip+1]],[y_20_flusso[ip],y_20_flusso[ip+1]],[z_20_flusso[ip],z_20_flusso[ip+1]], mutation_scale=18, lw=2, arrowstyle="->",color="green")
ax.add_artist(a20)


ax.set_title('$r=$''%g' % zeta, fontdict = font)
#ax.set_title('Modello Minsky basso r spazio delle fasi \n $\zeta=$' + '%g' % zeta + ', $\phi=$ ' + '%g' % phi, fontdict = font)
ax.set_xlabel('quota salari $\omega$')
ax.set_ylabel('occupazione $\lambda$')
ax.set_zlabel('quota debito privato $d$')
ax.set_xlim(-0.2,1.5)
ax.set_ylim(-0.2,1)
ax.set_zlim(-20,40)

# ZOOM sulla parte incasinata
#ax.set_xlim(-0.5,2)
#ax.set_ylim(-0.5,1.5)
#ax.set_zlim(-10,10)

ax.set_aspect("equal") # serve ad usare un'unica scala vedere se usare o meno, far la prova se piace o no: se attiva schiacci la figura
#plt.savefig('relazione corretta Minsky r basso spazio delle fasi equalizzato.jpg', dpi=850, transparent=False)
#%%
# ANDAMENTO TEMPORALE PARAGONE DIVERSE DINAMICHE COME PARE VOLERE TURCHETTI
# Forse aggiungere legenda col punto iniziale o informazioni ad esempio "dinamica sella"

colore = 'blue'

fig, ax1 = plt.subplots()
ax1.plot(t, x_2_flusso,'r',label='quota salari $\omega(t)$')
ax1.plot(t, y_2_flusso,'k',label='occupazione $\lambda(t)$')
#plt.legend(loc='lower right', prop={'size': 6})
#plt.title('converge al punto fisso stabile (omega_1,lambda_1,d_1)',fontdict = font)
ax1.set_xlabel('tempo t')
ax1.set_ylabel('$\omega, \lambda$')
plt.grid()
ax1.set_ylim(0.,1.2)
# scala del debito privato divergente creando una nuova figura
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(t, z_2_flusso, color=colore,label='quota debito privato $d(t)$')
#ax2.legend(loc='center',prop={'size': 6})
ax1.plot(np.nan, colore,label='quota debito $d(t)$')  # Creo un falso grafico con np.nan solo per potere collocare la legenda del debito privato nella prima, infatti è ax1
ax1.legend(loc='best',prop={'size': 8})
ax1.set_xlim(0,300)
ax2.set_ylabel('$d$', color=colore)  # we already handled the x-label with ax1
ax2.set_ylim(-5,5)
ax2.tick_params(axis='y', labelcolor=colore)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.savefig('relazione corretta Minsky r basso andamento temporale 1.jpg', dpi=650, transparent=False)
plt.show()

fig, ax1 = plt.subplots()
ax1.plot(t, x_3_flusso,'r',label='quota salari $\omega(t)$')
ax1.plot(t, y_3_flusso,'k',label='occupazione $\lambda(t)$')
#plt.legend(loc='lower right', prop={'size': 6})
#plt.title('converge al punto fisso stabile (omega_1,lambda_1,d_1) ed l'omega iniziale più estremo tra i possibili',fontdict = font)
ax1.set_xlabel('tempo t')
ax1.set_ylabel('$\omega, \lambda$')
plt.grid()
ax1.set_ylim(0,1.5)
# scala del debito privato divergente creando una nuova figura
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(t, z_3_flusso, color=colore,label='quota debito privato $d(t)$')
#ax2.legend(loc='center',prop={'size': 6})
ax1.plot(np.nan, colore,label='quota debito $d(t)$')  # Creo un falso grafico con np.nan solo per potere collocare la legenda del debito privato nella prima, infatti è ax1
ax1.legend(loc='best',prop={'size': 8})
ax1.set_xlim(0,300)
ax2.set_ylabel('$d$', color=colore)  # we already handled the x-label with ax1
ax2.set_ylim(0,5)
ax2.tick_params(axis='y', labelcolor=colore)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.savefig('relazione corretta Minsky r basso andamento temporale 2.jpg', dpi=650, transparent=False)
plt.show()

fig, ax1 = plt.subplots()
ax1.plot(t, x_17_flusso,'r',label='quota salari $\omega(t)$')
ax1.plot(t, y_17_flusso,'k',label='occupazione $\lambda(t)$')
#plt.legend(loc='lower right', prop={'size': 6})
#plt.title('converge al secondo punto fisso dei (0,0,d_0) che è stabile',fontdict = font)
ax1.set_xlabel('tempo t')
ax1.set_ylabel('$\omega, \lambda$')
plt.grid()
ax1.set_ylim(0,1)
# scala del debito privato divergente creando una nuova figura
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(t, z_17_flusso, color=colore,label='quota debito privato $d(t)$')
#ax2.legend(loc='center',prop={'size': 6})
ax1.plot(np.nan, colore,label='quota debito $d(t)$')  # Creo un falso grafico con np.nan solo per potere collocare la legenda del debito privato nella prima, infatti è ax1
ax1.legend(loc='best',prop={'size': 8})
ax1.set_xlim(0,100)
ax2.set_ylabel('$d$', color=colore)  # we already handled the x-label with ax1
ax2.set_ylim(0,30)
ax2.tick_params(axis='y', labelcolor=colore)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.savefig('relazione corretta Minsky r basso andamento temporale 3.jpg', dpi=650, transparent=False)
plt.show()

fig, ax1 = plt.subplots()
ax1.plot(t, x_4_flusso,'r',label='quota salari $\omega(t)$')
ax1.plot(t, y_4_flusso,'k',label='occupazione $\lambda(t)$')
#plt.legend(loc='lower right', prop={'size': 6})
#plt.title('va a -oo',fontdict = font)
ax1.set_xlabel('tempo t')
ax1.set_ylabel('$\omega, \lambda$')
plt.grid()
ax1.set_ylim(0,1)
# scala del debito privato divergente creando una nuova figura
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(t, z_4_flusso, color=colore,label='quota debito privato $d(t)$')
#ax2.legend(loc='center',prop={'size': 6})
ax1.plot(np.nan, colore,label='quota debito $d(t)$')  # Creo un falso grafico con np.nan solo per potere collocare la legenda del debito privato nella prima, infatti è ax1
ax1.legend(loc='best',prop={'size': 8})
ax1.set_xlim(0,100)
ax2.set_ylabel('$d$', color=colore)  # we already handled the x-label with ax1
ax2.set_ylim(-400,50)
ax2.tick_params(axis='y', labelcolor=colore)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
#plt.savefig('relazione corretta Minsky r basso andamento temporale 4.jpg', dpi=650, transparent=False)
plt.show()
# IN TOTALE QUINDI NELLA RELAZIONE NE METTI 4 DI GRAFICI TEMPORALI
#%%
#%%
#%%
#%%
#%%
#########################
#%%
#########################
#%%
# CONTROLLO CONDIZIONI DI EQUILIBRIO PUNTO FISSO VALORI FINITI ECONOMICAMENTE AUSPICABILE
# (per alcune cose è una ripetizione delle formul del punto di equilibrio messe sopra, effettivamente si potrebbe anche cancellare questa parte)
zeta=2.5/100
#zeta=15.2/100
r0=zeta
#r0 = r_1lambdanullomeno
#r0 = r_1lambdanullopiù   # Attento che devi usare il pi1 col segno + se usi questo! Anche questo non rispetto le conzioni per cui il punto fisso è instabile.
pi1=((ni*F)/G) - (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) )
#pi1=((ni*F)/G) + (ni/G)*( np.sqrt(E/(H + ni*(alfa + beta + gamma))) )
omegaeq1=1-pi1-r0*((pi1-ni*(alfa+beta))/(r0-(alfa+beta)))
lanbdaeq1=(B/C) - ( np.sqrt(A/(alfa + D)) )/C
deq1=(pi1-ni*(alfa+beta))/(r0-(alfa+beta))
controllare2=r0*((ni-deq1)*( ((2*E*G)/ni)/(F-G*(pi1/ni))**3 -ni)) # sì, è maggiore di zero poichè viene 0.0868958 per zeta=0.025= 2.5/100
controllare1=r0*(2-((ni-deq1)/ni)*( ((2*E*G)/ni)/(F-G*(pi1/ni))**3)) # controllare se (alfa+beta)=0.05=5/100 è MAGGIORE di questo? Sì, lo è perchè viene addirittura negativo: -0.042834 
print((alfa+beta)>controllare1) # è True
print(controllare2) # Dovrebbe essere positiva affinchè il punto fisso (omega_1,lambda_1,d_1) sia stabile. Se viene NEGATIVA è instabile per quel tasso di interesse r scelto!
#%%
#########################
#%%
#########################
#%%
#%%
#%%
#%%
#%%
#%%
#%%
#%%
# RIFACCIO DEI GRAFICI DELLA VECCHIA RELAZIONE MA TOGLIENDO TITOLO

# N.B.: si utilizzano i flussi creati sopra nella parte di prova, fissare z appropriato all'inizio

#quota_salari_su_output = risultati_integrazione[:,0]
#occupazione_su_forza_lavoro = risultati_integrazione[:,1]
#debito_privato_su_output = risultati_integrazione[:,2]
#quota_banchieri_su_output = (zeta + phi*debito_privato_su_output)*debito_privato_su_output
#quota_profitti_su_output = 1 - quota_salari_su_output - quota_banchieri_su_output
quota_banchieri_su_output = (zeta)*z_45_flusso
quota_profitti_su_output = 1 - x_45_flusso - quota_banchieri_su_output

#plt.plot(x_45_flusso[0:13000], y_45_flusso[0:13000])
plt.plot(x_45_flusso, y_45_flusso)
plt.legend(loc='best', prop={'size': 6})
plt.title('$r=$''%g' % zeta, fontdict = font)
plt.xlabel('quota salari $\omega$')
plt.ylabel('occupazione $\lambda$')
#plt.ylim(0.75,1)
#plt.xlim(0.78,0.97)
plt.grid()
#plt.savefig('relazione corretta Minsky basso r Cyclical and equilibrium time paths. Stability at low interest.jpg', dpi=550, transparent=False)
plt.show()
#%%
# GRAFICO SPECIALE con anche la funzione di investimento dei capitalisti (che NON sono gli Investimenti Netti I)

plt.plot(t, functionalform(E,F,G,H,(quota_profitti_su_output/ni)),'darkgreen', label='investimenti $k[\\frac{\pi}{\\nu}](t)$')
plt.plot(t, quota_profitti_su_output,'blue', label='quota profitti $\pi(t)$')
plt.plot(t, quota_banchieri_su_output,'r', label='quota banchieri $b(t)$')
#plt.legend(loc='lower left', prop={'size': 6})
plt.legend()
plt.title('$r=$''%g' % zeta, fontdict = font)
plt.xlabel('tempo t')
plt.ylabel('$k[\\frac{\pi}{\\nu}], \pi, b$')
plt.xlim(0,200)
plt.grid()
#plt.savefig('relazione corretta Minsky basso r Figure 5. Profit, investment, and bank share at low interest ???.jpg', dpi=450, transparent=False)
plt.show()