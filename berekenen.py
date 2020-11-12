
#%% hulp functie
import time

class timer():
    def __init__(self):
        self.begin_tijd = 0
        self.eind_tijd = 0
           
    def __enter__(self):
        self.begin_tijd = time.time()
        
    def __exit__(self, *args):
        self.eind_tijd = time.time()
        print('Time elapsed [s]: {}'.format(self.eind_tijd-self.begin_tijd))


#%%



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df_knp = pd.read_csv('knp.csv', delimiter=';')
df_lei = pd.read_csv('leidingen.csv', delimiter=';')


#%% Inlezen van Knooppunten
knp = df_knp['knp'].values




knp_1 = df_lei['knp 1'].values
knp_2 = df_lei['knp 2'].values
bob_1 = df_lei['bob 1'].values
bob_2 = df_lei['bob2'].values

assert (~np.isnan(bob_1) +  ~np.isnan(bob_2)).all(), 'Een leiding mist beide bob waardes.'

if (np.isnan(bob_1) ^ np.isnan(bob_2)).any():
    nanbobs = 0
    if np.isnan(bob_1).any():
        nanbob_1 = np.isnan(bob_1)
        bob_1[nanbob_1] = bob_2[nanbob_1]
        nanbobs += sum(nanbob_1)
    if np.isnan(bob_2).any():
        nanbob_2 = np.isnan(bob_2)
        bob_2[nanbob_2] = bob_1[nanbob_2]
        nanbobs += sum(nanbob_2)
        
#%% Opzetten BOB matrix
    


n = len(knp) # aantal knopen

Inf = np.Infinity #Begin waardes bij de knooppunten

putbodem = np.full(n, Inf, dtype=np.float64)
bob = np.full((n,n), Inf, dtype=np.float64)

for startput, eindput, bob_1_, bob_2_ in zip(knp_1, knp_2, bob_1, bob_2):
    
    for j, put in enumerate(knp):
        if put == startput:
            j1 = j
            if putbodem[j] > bob_1_:
                putbodem[j] = bob_1_
        if put == eindput:
            j2 = j
            if putbodem[j] > bob_2_:
                putbodem[j] = bob_2_
    
    if bob_1_ > bob_2_:
        bob[j1, j2] = bob_1_
    else:
        bob[j1, j2] = bob_2_
    
    # Bob matrix is symmetrisch
    bob[j2, j1] = bob[j1,j2]
    
    

#%% Dijkstra algorithme voor waterstand


wachtrij = np.full(n, True, dtype=bool)
vorige_ind = np.full(n, -1, dtype=np.int64) #Numpy does not allow NaN values in integer matrices/vectors
vorige_knp = np.full(n, -1, dtype=np.int64)
waterstand = np.full(n, Inf, dtype=np.float64)

beginput = 201829

for i, putnummer in enumerate(knp):
    if putnummer == beginput:
        waterstand[i] = putbodem[i] - 0.01
# Correctie met 0.01 voor drempelbepaling



# Dijkstra algorithme
while True:
    dist = Inf
    
    for i, knp_i in enumerate(knp):
        if wachtrij[i]:
            if waterstand[i] < dist:
                dist = waterstand[i]
                u = i
                knp_u = knp_i
    
    
    if dist == Inf:
        break
    
    
    wachtrij[u] = False
    
    for j, knp_ in enumerate(knp):
        if wachtrij[j] and bob[u, j] != Inf:
          
            if waterstand[u] > bob[u, j]:
                alt = waterstand[u]
            else:
                alt = bob[u,j]
            
            if alt < waterstand[j]:
                waterstand[j] = alt
                vorige_ind[j] = u
                vorige_knp[j] = knp_u
    


#%% Verlorenberging Waterstand bij begin en eindpunt leiding toevoegen
def leiding_waterstand(knp, waterstand, knp_1, knp_2):
    """
    Voegt de waterstand uit de dijkstra berekening toe aan een lijst met
    waterstanden voor de leidingen. Knooppunt 1 is het beginpunt van een
    leiding. Knooppunt 2 is het eindpunt van een leiding.

    Parameters
    ----------
    knp : np.array
        Lijst met alle knooppunten.
    waterstand : np.array
        Lijst met alle waterstanden uit de dijkstra voor elk individueel
        knooppunt.
    knp_1 : np.array
        Eerste knooppunten lijst voor alle leidingen.
    knp_2 : np.array
        Tweede knooppunten lijst voor alle leidingen.

    Returns
    -------
    np.array
        Twee lijsten met waterstanden voor het eerste en tweede knooppunt van 
        een leiding.
    """
    knp_1_ws_list = []
    knp_2_ws_list = []
    for knp_1_, knp_2_ in zip(knp_1, knp_2):
        knp_1_ws_list.append(waterstand[knp == knp_1_][0])
        knp_2_ws_list.append(waterstand[knp == knp_2_][0])
    return np.array(knp_1_ws_list), np.array(knp_2_ws_list)

knp_1_ws, knp_2_ws = leiding_waterstand(knp, waterstand, knp_1, knp_2)

#%% Inlezen Leidingen tabel

lengte = df_lei['lengte'].values
breedte = df_lei['breedte'].values
vorm = df_lei['vorm'].values
hoogte = df_lei['hoogte'].values

# Bij ronde buizen is de hoogte==breedte
hoogte[vorm==0] = breedte[vorm==0]


#%% Volume berekeningen voor de Bergingskromme
def volume_berekenen(ws_1, ws_2, bob_1, bob_2, h, b, l, vorm):
    """
    Berekend het volume water in een leiding waarvan de waterstanden aan beide 
    kanten bekent is. Als de leiding niet volledig gevuld is, wordt de
    waterstand linear geinterpoleerd over de buis en in segmenten uitgerekend.

    Parameters
    ----------
    ws_1 : float
        Waterstand aan het eerste uiteinde van de leiding [mNAP].
    ws_2 : float
        Waterstand aan het tweede uiteinde van de leiding [mNAP].
    bob_1 : float
        Bovenkant Onderkant Buis 1 [mNAP].
    bob_2 : TYPE
        Bovenkant Onderkant Buis 2 [mNAP].
    h : float
        Hoogte leiding [m].
    b : TYPE
        Breedte leiding [m].
    l : TYPE
        lengte leiding [m].
    vorm : int
        Vorm van de leiding.
        0: cilinder vormig (rond)
        1: eivormig
        2: 
            ...

    Returns
    -------
    volume : float
        Volume water in de leiding.

    """
    volume = 0
    bbb_1 = bob_1 + h
    bbb_2 = bob_2 + h
    
    if (ws_1 >= bbb_1) & (ws_2 >= bbb_2):
        if vorm==0:
            A_lei = np.pi*(.5*b)**2
            volume = A_lei * l
        elif vorm==1:
            A_lei = 0.510*h**2
            volume = A_lei * l
        #elif vorm==2:
            
    elif ws_1 >= bob_1 or ws_2 >= bob_2:
        N_seg = int(round(l*10) - 1)
        dL = l/N_seg
        
        bob_lin = np.linspace(min([bob_1, bob_2]), max([bob_1, bob_2]), N_seg + 1) #BoB linear geinterpoleerd over lengte buis
        bob_seg = (bob_lin[1::] + bob_lin[0:-1])/2
        if ws_1 == ws_2:
            ws_seg = np.ones(shape=bob_seg.shape)*ws_1
        else:
            ws_lin = np.linspace(min([ws_1, ws_2]), max([ws_1, ws_2]), N_seg + 1) # Waterstand linear geinterpoleerd over lengte buis
            ws_seg = (ws_lin[1::] + ws_lin[0:-1])/2
        segment_whs = ws_seg - bob_seg #Waterhoogtes per segment van de buis
        segment_whs[segment_whs < 0] = 0 #Negatieve waterhoogtes kennen we niet
        
        if vorm==0:
            for segment_wh in segment_whs:
                R = .5*b #Radius buis
                if segment_wh >= 2*R:
                    A_water = np.pi*R**2
                    volume += dL*A_water
                elif segment_wh == 0:
                    pass
                elif segment_wh >= R:
                    theta = 2*np.arccos((segment_wh-R)/R) # Radialen van de leiding die helemaal gevuld zijn met water
                    A_water = np.pi*R**2 - R**2/2*(theta - np.sin(theta)) # A_totaal - A_zonder_water
                    volume += dL*A_water # Inhoud segment
                elif segment_wh < R:
                    theta = 2*np.arccos((R-segment_wh)/R)
                    A_water = R**2/2*(theta - np.sin(theta))
                    volume += dL*A_water
        elif vorm==1:
            for segment_wh in segment_whs:
                A_totaal = 0.510*h**2
                if segment_wh >= h:
                    volume += dL*A_totaal
                else:
                    # Eerste zes polynomale coefficienten die de waterhoogte vs. wateroppervlak benaderen
                    x = [0.2425, 2.1142, -2.1657, 0.8277, 0.8396, -0.8558]
                    wh = segment_wh/h #genormaliseerde hoogte [0,1] m
                    # genormaliseerd oppervlak [0,1] m2
                    A_rel_vol = wh*x[0] + wh**2*x[1] + wh**3*x[2] + wh**4*x[3] + wh**5*x[4] + wh**6*x[5] 
                    A_vol = A_totaal * A_rel_vol # Gevuld oppervlakte
                    volume += dL*A_vol
    return volume


riool_waterstanden = np.arange(np.min(putbodem), np.max(waterstand), 0.01) # min: riool is leeg, max: riool is compleet gevuld
array_wat_vol_lei = [] # Volume gevuld bij een rioolwaterstand.
array_ver_berg = [] # Verloren berging bij een rioolwaterstand

with timer():
    vol_volume = np.full(len(bob_1), 0, dtype=float)
    bool_vol = np.full(len(bob_1), False, dtype=bool)
    vol_verloren_volume = np.full(len(bob_1), 0, dtype=float)
    bool_vol_verloren = np.full(len(bob_1), False, dtype=bool)
    
    for r_ws in riool_waterstanden:
        
        volume_water_lei = 0 # Totale volume water in leiding
        verloren_berging = 0 # Volume water dat verloren berging is
        
        for i, (knp_1_ws_, knp_2_ws_, bob_1_, bob_2_, l, b, vorm_, h) in enumerate(zip(knp_1_ws, knp_2_ws, bob_1, bob_2, lengte, breedte, vorm, hoogte)):
            
            # Als de riool waterstand hoger is dan verloren berging waterstand dan h
            if r_ws >= knp_1_ws_:
                ws_ver_berg_1 = knp_1_ws_
            else:
                ws_ver_berg_1 = r_ws
            if r_ws >= knp_2_ws_:
                ws_ver_berg_2 = knp_2_ws_
            else:
                ws_ver_berg_2 = r_ws
            
            if (r_ws >= bob_1_ + h) & (r_ws >= bob_2_ + h):
                if bool_vol[i]:
                    volume_water_lei += vol_volume[i]
                else:
                    vol_volume[i] = volume_berekenen(r_ws, r_ws, bob_1_, bob_2_, h, b, l, vorm_)
                    bool_vol[i] = True
                    volume_water_lei += vol_volume[i]
            else:
                volume_water_lei += volume_berekenen(r_ws, r_ws, bob_1_, bob_2_, h, b, l, vorm_)
            
            if (r_ws >= knp_1_ws_) & (r_ws >= knp_2_ws_):
                if bool_vol_verloren[i]:
                    verloren_berging += vol_verloren_volume[i]
                else:
                    vol_verloren_volume[i] = volume_berekenen(ws_ver_berg_1, ws_ver_berg_2, bob_1_, bob_2_, h, b, l, vorm_)
                    bool_vol_verloren[i] = True
                    verloren_berging += vol_verloren_volume[i]
            else:
                verloren_berging += volume_berekenen(ws_ver_berg_1, ws_ver_berg_2, bob_1_, bob_2_, h, b, l, vorm_)
            
            
        array_wat_vol_lei.append(volume_water_lei)
        array_ver_berg.append(verloren_berging)
            
    array_wat_vol_lei = np.array(array_wat_vol_lei)
    array_ver_berg = np.array(array_ver_berg)
    
array_berg = array_wat_vol_lei - array_ver_berg # totaal volume - verloren berging
#%% plot bergingskromme

plt.figure(figsize=(20,10))
plt.plot(array_ver_berg, riool_waterstanden, array_wat_vol_lei, riool_waterstanden, array_berg, riool_waterstanden)
plt.legend(['Verloren berging [m3]', 'Totale inhoud [m3]', 'Berging [m3]'])

#%% Afvoerend oppervlak
df_afv = pd.read_csv('afv.csv', delimiter=';')
nlei = len(knp_1)

geslhell = np.nansum(df_afv['gesl hell'])
geslvlak = np.nansum(df_afv['gesl vlak'])
gesluitg = np.nansum(df_afv['gesl uitg'])
openhell = np.nansum(df_afv['open hell'])
openvlak = np.nansum(df_afv['open vlak'])
openuitg = np.nansum(df_afv['open uitg'])
dakhell = np.nansum(df_afv['dak hell'])
dakvlak = np.nansum(df_afv['dak vlak'])
dakuitg = np.nansum(df_afv['dak uitg'])
onvhell = np.nansum(df_afv['onv hell'])
onvvlak = np.nansum(df_afv['onv vlak'])
onvuitg = np.nansum(df_afv['onv uitg'])

totaalvo = np.nansum([geslhell, geslvlak, gesluitg, openhell, openvlak, openuitg, dakhell, dakvlak, dakuitg, onvhell, onvvlak, onvuitg])

