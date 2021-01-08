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
        'RO': cilinder vormig (rond)
        'EI': eivormig
        'RE': rechthoekig
        'MU': muilvormig
    Returns
    -------
    volume : float
        Volume water in de leiding.
    """
    volume = 0
    bbb_1 = bob_1 + h
    bbb_2 = bob_2 + h
    
    if (ws_1 >= bbb_1) & (ws_2 >= bbb_2):
        if vorm=='RO':
            A_lei = np.pi*(.5*h)**2
            volume = A_lei * l
        elif vorm=='EI':
            A_lei = 0.510*h**2
            volume = A_lei * l
        elif vorm=='RE':
            A_lei = h*b
            volume = A_lei * l
        elif vorm=='MU':
            A_lei = 0.7645 * h * b
            volume = A_lei * l
            
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
        segment_whs = ws_seg - bob_seg #Waterdieptes per segment van de buis
        segment_whs[segment_whs < 0] = 0 #Negatieve waterdieptes kennen we niet
        
        if vorm=='RO':
            for segment_wh in segment_whs:
                R = .5*h #Radius buis
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
        elif vorm=='EI':
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
        elif vorm=='RE':
            for segment_wh in segment_whs:
                A_totaal = 0.510*h**2
                if segment_wh >= h:
                    volume += dL*h
                else:
                    volume += dL*segment_wh
        elif vorm=='MU':
            for segment_wh in segment_whs:
                A_totaal = 0.7645 * h * b
                if segment_wh >= h:
                    volume += dL*A_totaal
                else:
                    # Eerste zes polynomale coefficienten die de waterhoogte vs. wateroppervlak benaderen
                    x = [1.4114, -6.9511, 12.33, -11.014, 5.0266, 0.1994]
                    wh = segment_wh/h #genormaliseerde diepte [0,1] m/m
                    # genormaliseerd oppervlak [0,1] 
                    A_rel_vol = wh*x[0] + wh**2*x[1] + wh**3*x[2] + wh**4*x[3] + wh**5*x[4] + wh**6*x[5] 
                    A_vol = A_totaal * A_rel_vol # Gevuld oppervlakte
                    volume += dL*A_vol

    return volume

def bergingsberekening(
        bemalingsgebied_naam,
        bemalingsgebied_id,
        put_1,
        put_2,
        bob_1,
        bob_2,
        lengte,
        vorm,
        hoogte,
        breedte,
        afvoerputten,
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
    afvoerputten : np.array of int
        Afvoerput ids waar de gemalen staan.
    overstorten : list or np.array
        Lijst aan overstortdrempel niveaus in het rioleringstelsel .
    Returns
    -------
    fig 
        Figure of the bergingskromme.
        
    '''
    putten = np.unique(np.concatenate([put_1, put_2])) #Unieke putten
    
    n = len(putten)
    N_streng = len(bob_1)
    
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
        
        if np.isnan(bob_1_) & np.isnan(bob_2_):
            # Bob graafmatrix is symmetrisch
            bob_graaf[j1, j2] = bob_graaf[j2, j1] = np.nan
        else:
            bob_graaf[j1, j2] = bob_graaf[j2, j1] = np.nanmax([bob_1_, bob_2_])

        
        
    wachtrij = np.full(n, True, dtype=bool)
    vorige_i = np.full(n, -1, dtype=np.int64) 
    vorige_put = np.full(n, -1, dtype=np.int64)
    waterstand = np.full(n, Inf, dtype=np.float64) #Verloren waterstand
    
    check = [(putten == afvoerput).any() for afvoerput in afvoerputten]
    assert any(check), 'Geen van de afvoerputten {} is onderdeel van stelsel {}.'.format(afvoerputten, bemalingsgebied_naam)
    
    # Zorg dat de bemalingsput/beginput leeg is, waterstand == bob
    # De -0.01 extra werd door Albert kemeling aanbevolen
    for afvoerput in afvoerputten:
        waterstand[putten == afvoerput] = putbodem[putten == afvoerput] -.01
    
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
        
        for j, put_j in enumerate(putten):
            if wachtrij[j] and bob_graaf[u, j] != Inf:
                
                if np.isnan(bob_graaf[u, j]):
                    bob_graaf[u, j] = waterstand[u]


                if waterstand[u] > bob_graaf[u, j]:
                    alt = waterstand[u]
                else:
                    alt = bob_graaf[u,j]
                
                if alt < waterstand[j]:
                    waterstand[j] = alt
                    vorige_i[j] = u
                    vorige_put[j] = put_u
    
    #Verloren waterstand wordt bij beginput en eindput van elk streng gezocht
    put_1_ws, put_2_ws = leiding_waterstand(putten, waterstand, put_1, put_2)
    
    ### vervalputten
    
    verval = np.full(n, -1, dtype=np.float)
    verval_strengen = np.empty(n, dtype=object)
    # Dijkstra geeft een stroomvolgorde aan (vorige_put naar put) voor
    # het kritieke pad.
    # Hiermee kunnen we gemakkelijk het verval per put berekenen op het
    # kritieke pad.
    kritieke_pad_strengen = np.full(N_streng, False, dtype=bool)
    
    for i, (vp, hp) in enumerate(zip(vorige_put, putten)):
        ind_1 = (put_1 == hp) & (put_2 == vp)
        ind_2 = (put_2 == hp) & (put_1 == vp)
        if ind_1.any():
            bob_temp = bob_2[ind_1]
            ws = put_2_ws[ind_1]
            if bob_temp > ws:
                verval[i] = bob_temp - ws
                verval_strengen[i] = str(put_1[ind_1][0]) + '-' + str(put_2[ind_1][0])
            kritieke_pad_strengen[ind_2] = True
        elif ind_1.any():
            bob_temp = bob_1[ind_2]
            ws = put_1_ws[ind_2]
            if bob_temp > ws:
                verval[i] = bob_temp - ws
                verval_strengen[i] = str(put_1[ind_2][0]) + '-' + str(put_2[ind_2][0])
            kritieke_pad_strengen[ind_1] = True
    
    # Voor de strengen die niet op het kritieke pad liggen is er de volgende
    # berekening.
    a = ~kritieke_pad_strengen
    for bob_1_, bob_2_, put_1_, put_2_, ws_1, ws_2 in zip(bob_1[a], bob_2[a], put_1[a], put_2[a], put_1_ws[a], put_2_ws[a]):
        # De stroomrichting wordt door de bob's bepaald, als de bob's gelijk 
        # zijn dan wordt het verval naar beide putten uitgerekend.
        if bob_2_ > bob_1_:
            hp = put_1_
            bob_hp = bob_1_
            ws_hp = ws_1
            if bob_hp > ws_hp:
                ind = putten == hp
                verval_hoogte = bob_hp - ws_hp
                if verval[ind] < verval_hoogte:
                    verval[ind] = verval_hoogte
                    verval_strengen[ind] = str(put_1_) + '-' + str(put_2_)
        elif bob_1_ > bob_2_:
            hp = put_2_
            bob_hp = bob_2_
            ws_hp = ws_2
            if bob_hp > ws_hp:
                ind = putten == hp
                verval_hoogte = bob_hp - ws_hp
                if verval[ind] < verval_hoogte:
                    verval[ind] = verval_hoogte
                    verval_strengen[ind] = str(put_1_) + '-' + str(put_2_)
        elif bob_1_ == bob_2_:
            if bob_1_ > ws_1:
                ind = putten == put_1_
                verval_hoogte = bob_1_ - ws_1
                if verval[ind] < verval_hoogte:
                    verval[ind] = verval_hoogte
                    verval_strengen[ind] = str(put_1_) + '-' + str(put_2_)
            if bob_2_ > ws_2:
                ind = putten == put_2_
                verval_hoogte = bob_2_ - ws_2
                if verval[ind] < verval_hoogte:
                    verval[ind] = verval_hoogte
                    verval_strengen[ind] = str(put_1_) + '-' + str(put_2_)
                
    
    ### bergingskromme
    
    # De volgende arrays zijn ervoor om te verkomen dat het volle (verloren) volume
    # van strengen niet meerdere keren berekend wordt.
    # elke streng krijgt een ruimte om het volle volume op te slaan en de 
    # boolean: bool_vol en bool_vol_verloren zorgen ervoor dat ze bij de volgende
    # iteratie overgeslagen worden.
    vol_volume = np.full(N_streng, np.nan, dtype=float) 
    bool_vol = np.full(N_streng, False, dtype=bool)
    vol_verloren_volume = np.full(N_streng, np.nan, dtype=float)
    bool_vol_verloren = np.full(N_streng, False, dtype=bool)
    
    list_volume_water = []
    list_volume_verlorenberging = []
  
    riool_waterstanden = np.arange(np.nanmin(putbodem), np.nanmax(waterstand) + .01, 0.01) # min: riool is leeg, max: riool is compleet gevuld
    
    
    for r_ws in riool_waterstanden:
        
        volume_water_lei = 0 # Totale volume water in de strengen bij water niveau r_ws
        verloren_berging = 0 # Volume verloren berging
        
        for i, (put_1_ws_, put_2_ws_, bob_1_, bob_2_, l, b, vorm_, h) in enumerate(zip(put_1_ws, put_2_ws, bob_1, bob_2, lengte, breedte, vorm, hoogte)):
            
            # Als beide bob's missen dan berekenen we het volume niet voor deze buis
            # Als 1 bob mist, bob_1 == bob_2, de tweede lijn is misschien overbodig
            # vanwege de check die eerder gedaan wordt maar geeft wel extra robuustheid
            # voor veranderingen.
            if np.isnan(bob_1_) & np.isnan(bob_2_):
                continue # Skip this loop
            elif np.isnan(bob_1_) | np.isnan(bob_2_):
                bob_1_ = bob_2_ = np.nanmax([bob_1_, bob_2_])

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
            if (r_ws >= (bob_1_ + h)) & (r_ws >= (bob_2_ + h)):
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
    
    info_berging = np.column_stack([np.array(riool_waterstanden), array_volume_water, array_verlorenberging, array_berging])
    info_putten = np.column_stack([putten, vorige_put, putbodem, waterstand, verval, verval_strengen])
    info_strengen = np.column_stack([put_1, put_2, bob_1, bob_2, put_1_ws, put_2_ws,lengte, hoogte, vol_volume, vol_verloren_volume, vol_verloren_volume/vol_volume])
    
    ### Plotting
    fig = plt.figure(figsize=(20,12))
    
    plt.plot(array_volume_water, riool_waterstanden, array_verlorenberging, riool_waterstanden, array_berging, riool_waterstanden)
    plt.legend(['Totale inhoud [$m^3$]', 'Verloren berging [$m^3$]',  'Berging [$m^3$]'])
    plt.xlabel('$m^3$')
    plt.ylabel('$mNAP$')
    plt.title(str(bemalingsgebied_id) + ', ' + bemalingsgebied_naam)
    
    if len(overstorten) > 0:
        b = riool_waterstanden < np.min(overstorten)
        volume_ovs = np.nanmax(array_volume_water[b]) #hoogste watervolume onder de laagste overstort
        verloren_berg_ovs = np.nanmax(array_verlorenberging[b])
        berg_ovs = np.nanmax(array_berging[b])
        for overstort in overstorten:
            plt.plot([np.min(array_volume_water), np.max(array_volume_water)], [overstort, overstort])

        text = ('Overstortdrempel= {:.1f} $mNAP$ \nOnder overstort:\nTotaal volume= {:.0f} $m^3$'.format(np.min(overstorten), volume_ovs) +
        '\nVerloren berging = {:.0f} $m^3$\nBerging = {:.0f} $m^3$'.format(verloren_berg_ovs, berg_ovs))
    else:
        volume_max = np.nanmax(array_volume_water) #hoogste watervolume onder de laagste overstort
        verloren_berg_max = np.nanmax(array_verlorenberging)
        berg_max = np.nanmax(array_berging)
        
        text = ('Max volume= {:.0f} $m^3$'.format(volume_max) +
        '\nMax Verloren berging = {:.0f} $m^3$\nMax Berging = {:.0f} $m^3$'.format(verloren_berg_max, berg_max))
        
    ax = plt.gca()
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.70, 0.20, text, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    
    return fig, info_berging, info_putten, info_strengen

#%%
import cx_Oracle as ora
import os, traceback, time
import numpy as np
import matplotlib.pyplot as plt
import logging


#%%

folder = r'\\rotterdam.local\vdfs001\HOMEDIRECTORY\UserHome01\932829\Mijn Documenten\Bergingsberekening\Test 7 albert'

## TO DO: sql folder
sql_folder = r'\\rotterdam.local\vdfs001\HOMEDIRECTORY\UserHome01\932829\Mijn Documenten\Bergingsberekening\Test 04'

sql_bemaling_link = 'bemalingsgebieden.sql'
sql_strengen_link = 'strengen.sql'
sql_overstorten_link = 'overstorten.sql'

sql_bemaling = open(os.path.join(sql_folder, sql_bemaling_link)).read()
sql_strengen = open(os.path.join(sql_folder, sql_strengen_link)).read()
sql_overstorten = open(os.path.join(sql_folder, sql_overstorten_link)).read()


#%% Alle verschillende bemalingsgebieden + afvoerput ids ophalen
with ora.connect('wmorplg', 'raadplegen_wmo', 'WMO.PROD') as con:
    cur = con.cursor()
    
    
    cur.execute(sql_bemaling)
    
    list = cur.fetchall()
    columns = np.array(cur.description)[:, 0]
    bemaling_temp = np.array(list)
    bemaling_ids = bemaling_temp[:,0].astype(np.int)
    bemalingsgebieden = []
    # Bemalingsgebieden staan dubbel voor meer dan 1 afvoerput
    # Hier worden de afvoerputten bij  hetzelfde bemalingsgebied gemerged
    for gemaal_id in np.unique(bemaling_ids):
        i = bemaling_ids == gemaal_id
        
        bemalingsgebied_naam = bemaling_temp[i][0, 1]
        afvoerputten = bemaling_temp[i][:,2].astype(np.int)
        bemalingsgebieden.append([gemaal_id, bemalingsgebied_naam, afvoerputten])

#%% Loop over alle belangrijke bemalingsgebieden

for i in range(0, 100):
    print('--------------\nBergings berekening voor {}, {} \n'.format(bemalingsgebieden[i][0], bemalingsgebieden[i][1]))
    
    bemalingsgebied = bemalingsgebieden[i]
    gemaal_id = str(bemalingsgebied[0])
    gemaal_naam = bemalingsgebied[1].replace('.','').replace('/', '-')
    
    subfolder_name = gemaal_id.zfill(4) + '_' + gemaal_naam
    subfolder = os.path.join(folder, subfolder_name)
    if os.path.exists(subfolder) is False:
        os.mkdir(subfolder)
    
    log = logging.getLogger(subfolder)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fileHandler = logging.FileHandler(os.path.join(subfolder, subfolder_name + '.log'), mode='w')
    fileHandler.setFormatter(formatter)
    log.setLevel(logging.INFO)
    log.addHandler(fileHandler)    
    
    log.info('Start bergingsberekening voor %s', gemaal_naam)
    
    try:
        
        
        afvoerputten = bemalingsgebied[2]
        assert len(afvoerputten)> 0, 'Er zijn geen afvoerputten van dit stelsel gevonden'
        
        # Query van de overstorten
        with ora.connect('wmorplg', 'raadplegen_wmo', 'WMO.PROD') as con:
            cur = con.cursor()
            
            
            cur.execute(sql_overstorten, [gemaal_id])
            
            list = cur.fetchall()
            ovs_col = np.array(cur.description)[:, 0]
            ovs = np.array(list)
        
        #TO DO: warning not exception let op geen overstort
        if ovs.size != 0:
            ovs = ovs[:,2]
            overstorten = [float(ov) for ov in ovs if ov != None]
            log.info('Er zijn %i overstorten gevonden.', len(overstorten))
            #ovs[:,2][ovs[:,2] != None].astype(np.float)
        else:
            overstorten = []
            log.warning('Er zijn geen overstorten gevonden.')
        
        def queryLeidingen(gemaal_id):
            with ora.connect('wmorplg', 'raadplegen_wmo', 'WMO.PROD') as con:
                cur = con.cursor()
                
                cur.execute(sql_strengen, [gemaal_id])
                fetched = cur.fetchall()
                col_lei = np.array(cur.description)[:,0]
                array_lei = np.array(fetched)
                return array_lei, col_lei

        # Query van de leidingen
        array_lei, col_lei = queryLeidingen(gemaal_id)
        
        assert array_lei.size !=0, 'Er zijn geen strengen van dit stelsel gevonden.'
        
        log.info('Er zijn %i strengen gevonden', len(array_lei[:,1]))
        
        # Initialisatie van de belangrijke arrays
        put_1 = array_lei[:,2].astype(np.int)
        put_2 = array_lei[:,3].astype(np.int)
        bob_1 = array_lei[:,4].astype(np.float)
        bob_2 = array_lei[:,6].astype(np.float)
        lengte = array_lei[:,8].astype(np.float)
        vorm = array_lei[:,9]
        hoogte = array_lei[:,10].astype(np.float)/1000
        breedte = array_lei[:,11].astype(np.float)/1000
        diameter = array_lei[:,12].astype(np.float)/1000
        hoogte[~np.isnan(diameter)] = diameter[~np.isnan(diameter)]
        
        
        # Als een van de twee bob waardes mist wordt die vervangen door de andere.
        if (np.isnan(bob_1) | np.isnan(bob_2)).any():
            nanbobs = 0
            nanbob_1 = np.isnan(bob_1)
            nanbob_2 = np.isnan(bob_2)
            if nanbob_1.any():
                bob_1[nanbob_1] = bob_2[nanbob_1]
                nanbobs += sum(nanbob_1)
            if nanbob_2.any():
                bob_2[nanbob_2] = bob_1[nanbob_2]
                nanbobs += sum(nanbob_2)
            
            streng_ids = array_lei[:,0][(nanbob_1 | nanbob_2)]
            log.warning('Er zijn %i streng(en) die een of twee bob\'s missen', len(streng_ids))
            log.info('Onvolledige strengen:\n' + np.array2string(streng_ids, max_line_width=40, separator=';'))
                   
            # Beide bob waardes van een leiding missen, deze kunnen essentieel zijn
            # voor een bergingsberekening. Het algorithme geeft dan verkeerde waardes
            # uit.
            both_nanbobs = np.sum(nanbob_1 + nanbob_2)
            #assert (nanbob_2 & nanbob_2).all(), 'Er zijn {} strengen die beide bob waardes missen.'.format(both_nanbobs)


        dict_h = {
            'bemalingsgebied_id':gemaal_id,
            'bemalingsgebied_naam':gemaal_naam,
            'put_1': put_1,
            'put_2': put_2,
            'bob_1': bob_1,
            'bob_2': bob_2,
            'lengte': lengte,
            'vorm': vorm,
            'hoogte': hoogte,
            'breedte': breedte,
            'afvoerputten': afvoerputten,
            'overstorten': overstorten
        }
        log.info('Berekening Berging')
        fig, info_berging, info_putten, info_strengen = bergingsberekening(**dict_h)
        log.info('Bergingsberekening wordt opgeslagen in %s', subfolder)
        
        np.savetxt(os.path.join(subfolder, gemaal_id.zfill(4) + '_berging.txt'), info_berging, delimiter=';', fmt='%4.2f;%.2f;%.2f;%.2f', header='Riool Waterstand[mNAP];Volume Water[m3];Verloren Berging[m3];Berging[m3]')
        np.savetxt(os.path.join(subfolder, gemaal_id.zfill(4) + '_putten.txt'), info_putten, delimiter=';', fmt='%i;%i;%.2f;%.2f;%.2f;%s', header='Put;Vorige Put;Putbodem[mNAP];Waterstand[mNAP];Verval Hoogte[m];Verval Streng')
        np.savetxt(os.path.join(subfolder, gemaal_id.zfill(4) + '_strengen.txt'), info_strengen, delimiter=';', fmt='%i;%i;%.2f;%.2f;%.2f;%.2f;%.2f;%.2f;%.3f;%.3f;%.2f', header='put 1;put 2;bob 1[mNAP];bob 2[mNAP];ws 1[mNAP];ws 2[mNAP];Lengte[m];Hoogte[m];Totaal Volume[m3];Verloren Volume[m3];Ratio')
        fig.savefig(os.path.join(subfolder, gemaal_id.zfill(4) + '_' + gemaal_naam + '.png'))
        
        log.info('Succes!')
        
        
    except AssertionError as e:
        print(e)
        log.exception('Assertion Error')
        log.info('Bergingsberekening afgebroken')
    except Exception as e:
        print(traceback.format_exc())
        log.exception('Exception Occurred')
        log.info('Bergingsberekening afgebroken')

    print('--------------\n')