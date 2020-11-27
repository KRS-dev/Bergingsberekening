# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 10:57:22 2020

@author: Desktop
"""
def bergingsberekening(
        put_1,
        put_2,
        bob_1,
        bob_2,
        lengte,
        vorm,
        hoogte,
        breedte,
        begin_put,
        overstorten
        ):
    '''
    Parameters
    ----------
    put_1 : np.array of int
        beginput.
    put_2 : np.array of int
        eindput.
    bob_1 : np.array of float
        bovenkant onderkant buis begin.
    bob_2 : np.array of float
        bob eind.
    lengte : np.array of float
        lengte buis.
    vorm : np.array of varchar
        Vorm van de leiding.
        'RO': cilinder vormig (rond)
        'EI': eivormig
        'RE': rechthoekig
        'MU': muilvormig
    hoogte : np.array of float
        hoogte/diameter buis.
    breedte : np.array of float, nullable
        breedte buis.
    begin_put : int
        Afvoer put waar het gemaal geplaatst is.
    overstorten : list or np.array
        Lijst aan overstortdrempel niveaus in het rioleringstelsel .

    Returns
    -------
    None.

    '''
    putten = np.unique(np.concatenate([put_1, put_2])) #Unieke putten
    
    n = len(putten)
    
    Inf = np.Infinity
    
    putbodem = np.full(n, Inf, dtype=np.float64) # Laagste bob bij een put
    
    bob_graaf = np.full((n,n), Inf, dtype=np.float64)  
    # Vierkante matrix van putten x putten, put[i] heeft een connectie met put[j]
    # dus er is non-infinite bobwaarde in de bob_graaf[i,j] = max(bob_1, bob_2)
    # het waterniveau in het systeem moet hoger zijn dan max(bob_1, bob_2) om
    # van put[i] naar put[j] te stromen
    
    for startput, eindput, bob_1_, bob_2_ in zip(put_1, put_2, bob_1, bob_2):
    
        for j, put in enumerate(putten):
            if put == startput:
                j1 = j
                if putbodem[j] > bob_1_:
                    putbodem[j] = bob_1_
            if put == eindput:
                j2 = j
                if putbodem[j] > bob_2_:
                    putbodem[j] = bob_2_
        
        
        bob_graaf[j1, j2] = max([bob_1_, bob_2_])
        # Bob graafmatrix is symmetrisch
        bob_graaf[j2, j1] = bob_graaf[j1,j2]
        
        
    wachtrij = np.full(n, True, dtype=bool)
    vorige_i = np.full(n, -1, dtype=np.int64) 
    vorige_put = np.full(n, -1, dtype=np.int64)
    waterstand = np.full(n, Inf, dtype=np.float64) #Verloren waterstand
    
    # Zorg dat de bemalingsput/beginput leeg is, waterstand == bob
    # De -0.01 extra werd door Albert kemeling aanbevolen
    waterstand[putten == beginput] = putbodem[putten == beginput] -.01
    ###Dijkstra
    while True:
        dist = Inf
        
        for i, put_i in enumerate(putten):
            if wachtrij[i]:
                if waterstand[i] < dist:
                    dist = waterstand[i]
                    u = i
                    put_u = put_i
        
        
        if dist == Inf:
            break
        
        
        wachtrij[u] = False
        
        for j, put_ in enumerate(putten):
            if wachtrij[j] and bob_graaf[u, j] != Inf:
              
                if waterstand[u] > bob_graaf[u, j]:
                    alt = waterstand[u]
                else:
                    alt = bob_graaf[u,j]
                
                if alt < waterstand[j]:
                    waterstand[j] = alt
                    vorige_i[j] = u
                    vorige_put[j] = put_u
    
    

    
    
    ### bergingskromme
    
    #Verloren waterstand wordt bij beginput en eindput van elk streng gezocht
    put_1_ws, put_2_ws = leiding_waterstand(putten, waterstand, put_1, put_2)
    
    N_streng = len(bob_1)
    
    # De volgende arrays zijn ervoor om te verkomen dat het volle (verloren) volume
    # van strengen niet meerdere keren berekend wordt.
    # elke streng krijgt een ruimte om het volle volume op te slaan en de 
    # boolean: bool_vol en bool_vol_verloren zorgen ervoor dat ze bij de volgende
    # iteratie overgeslagen worden.
    vol_volume = np.full(N_streng, 0, dtype=float) 
    bool_vol = np.full(N_streng, False, dtype=bool)
    vol_verloren_volume = np.full(N_streng, 0, dtype=float)
    bool_vol_verloren = np.full(N_streng, False, dtype=bool)
    
    list_volume_water = []
    list_volume_verlorenberging = []
    
    
    riool_waterstanden = np.arange(np.nanmin(putbodem), np.nanmax(waterstand), 0.01) # min: riool is leeg, max: riool is compleet gevuld
    
    for r_ws in riool_waterstanden:
        
        volume_water_lei = 0 # Totale volume water in de strengen bij water niveau r_ws
        verloren_berging = 0 # Volume verloren berging
        
        for i, (put_1_ws_, put_2_ws_, bob_1_, bob_2_, l, b, vorm_, h) in enumerate(zip(put_1_ws, put_2_ws, bob_1, bob_2, lengte, breedte, vorm, hoogte)):
            
            # Als de riool waterstand hoger is dan verloren berging waterstand dan h
            if r_ws >= put_1_ws_:
                ws_ver_berg_1 = put_1_ws_
            else:
                ws_ver_berg_1 = r_ws
                
            if r_ws >= put_2_ws_:
                ws_ver_berg_2 = put_2_ws_
            else:
                ws_ver_berg_2 = r_ws
            
            # Komt r_ws boven BBB uit? dan is de streng vol. 
            # Reken anders het gedeeltelijk volume uit
            if (r_ws >= bob_1_ + h) & (r_ws >= bob_2_ + h):
                if bool_vol[i]:
                    volume_water_lei += vol_volume[i]
                else:
                    vol_volume[i] = volume_berekenen(r_ws, r_ws, bob_1_, bob_2_, h, b, l, vorm_)
                    bool_vol[i] = True
                    volume_water_lei += vol_volume[i]
            else:
                volume_water_lei += volume_berekenen(r_ws, r_ws, bob_1_, bob_2_, h, b, l, vorm_)
            
            # Komt r_ws boven de verloren waterstand uit? Dan is de verloren berging compleet gevuld
            # Reken anders de gedeeltelijke verloren berging uit.
            if (r_ws >= put_1_ws_) & (r_ws >= put_2_ws_):
                if bool_vol_verloren[i]:
                    verloren_berging += vol_verloren_volume[i]
                else:
                    vol_verloren_volume[i] = volume_berekenen(ws_ver_berg_1, ws_ver_berg_2, bob_1_, bob_2_, h, b, l, vorm_)
                    bool_vol_verloren[i] = True
                    verloren_berging += vol_verloren_volume[i]
            else:
                verloren_berging += volume_berekenen(ws_ver_berg_1, ws_ver_berg_2, bob_1_, bob_2_, h, b, l, vorm_)
            
        
        list_volume_water.append(volume_water_lei)
        list_volume_verlorenberging.append(verloren_berging)
    
    # Maak arrays van de lijst data
    array_volume_water = np.array(list_volume_water)
    array_verlorenberging = np.array(list_volume_verlorenberging)
    
    array_berging = array_volume_water - array_verlorenberging # totaal volume - verloren berging
    
    
    ### Plotting
    fig = plt.figure(figsize=(20,12))
    
    plt.plot(array_volume_water, riool_waterstanden, array_verlorenberging, riool_waterstanden, array_berging, riool_waterstanden)
    plt.legend(['Totale inhoud [$m^3$]', 'Verloren berging [$m^3$]',  'Berging [$m^3$]'])
    
    for overstort in overstorten:
        plt.plot([np.min(array_volume_water), np.max(array_volume_water)], [overstort, overstort])
    plt.xlabel('$m^3$')
    plt.ylabel('$mNAP$')
    
    
    b = riool_waterstanden < np.min(overstorten)
    
    volume_ovs = np.max(array_volume_water[b])
    verloren_berg_ovs = np.max(array_verlorenberging[b])
    berg_ovs = np.max(array_berging[b])
    
    text = ('Overstortdrempel= {:.1f} $mNAP$ \nOnder overstort:\nTotaal volume= {:.2f} $m^3$'.format(np.min(overstorten), volume_ovs) +
    '\nVerloren berging = {:.2f} $m^3$\nBerging = {:.2f} $m^3$'.format(verloren_berg_ovs, berg_ovs))
    
    ax = plt.gca()
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.70, 0.20, text, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    
    return fig    