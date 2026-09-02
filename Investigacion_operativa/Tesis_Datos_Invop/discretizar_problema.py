import numpy as np
import scipy.sparse as sp
import data_generator as gen
import gurobipy as gp
from gurobipy import GRB, quicksum
import time

def construir_grilla(coord_x,coord_y):
    """
    La grilla se compone de nodos organizados en filas y columnas
    y formará rectángulos. El modelo no requiere que los nodos de cada fila
    y de cada columna sean equidistantes.
    """
    grilla = []
    for y0 in coord_y:
        grilla.append([])
        for x0 in coord_x:
            grilla[-1].append(np.array([x0,y0]))
    return grilla

#funciona

def vecinos_de_nodo_intragrilla(ind_fila,ind_col,fil,col,grilla):
    """
    Esta función toma un nodo de la grilla y lista sus vecinos de la grilla,
    así como la distancia a cada uno.
    Así se obtiene luego las únicas aristas que podrá recorrer el vehículo base.
    PROCEDEMOS A ELIMINAR LAS DIAGONALES
    """
    vecinos = []
    nodo = grilla[ind_fila][ind_col]
    if ind_fila > 0:
        vecino = grilla[ind_fila-1][ind_col]
        vecinos.append((vecino[0],vecino[1],np.linalg.norm(nodo-vecino)))
    if ind_fila < fil-1:
        vecino = grilla[ind_fila+1][ind_col]
        vecinos.append((vecino[0],vecino[1],np.linalg.norm(nodo-vecino)))
    if ind_col > 0:
        vecino = grilla[ind_fila][ind_col-1]
        vecinos.append((vecino[0],vecino[1],np.linalg.norm(nodo-vecino)))
    if ind_col < col-1:
        vecino = grilla[ind_fila][ind_col+1]
        vecinos.append((vecino[0],vecino[1],np.linalg.norm(nodo-vecino)))
    return vecinos

#funciona

def construir_vecinos_intragrafo(grafo):
    """
    Construimos un diccionario de vecinos de los nodos de un grafo ya discretizado.
    Aún no tienen los índices, sino que están dados por sus coordenadas en el plano.
    También se agrega la longitud de la arista, que será luego el costo de atravesarla.
    """
    dic_vecinos = dict()
    dic_dos_sentidos = dict()
    for lista_nodos in grafo: #Ahora cada lista tiene sólo dos nodos.
        dic_vecinos[(lista_nodos[0][0],lista_nodos[0][1])] = []
        dic_vecinos[(lista_nodos[1][0],lista_nodos[1][1])] = []
        dic_dos_sentidos[(lista_nodos[0][0],lista_nodos[0][1])] = []
        dic_dos_sentidos[(lista_nodos[1][0],lista_nodos[1][1])] = []
    for lista_nodos in grafo:
        #Ponemos los vecinos con longitud de la arista
        longitud = np.linalg.norm(lista_nodos[0]-lista_nodos[1])
        #Armamos las vecindades en un solo sentido para luego ubicar con mayor facilidad los dos sentidos de cada arista
        nodo = (lista_nodos[0][0],lista_nodos[0][1])
        vecino = (lista_nodos[1][0],lista_nodos[1][1],longitud)
        dic_vecinos[nodo].append(vecino)
        dic_dos_sentidos[nodo].append(vecino)
        dic_dos_sentidos[(vecino[0],vecino[1])].append((nodo[0],nodo[1],longitud))
    return dic_vecinos, dic_dos_sentidos

#funciona

def construir_lista_vecinos_intragrafos(lista_grafos):
    """
    Construimos las vecindades que constituyen aristas de recorrido obligatorio en cada grafo.
    """
    lista_vecindades_intragrafo = []
    lista_vecindades_dos_sentidos = []
    for grafo in lista_grafos:
        dic_vecinos, dic_dos_sentidos = construir_vecinos_intragrafo(grafo)
        lista_vecindades_intragrafo.append(dic_vecinos)
        lista_vecindades_dos_sentidos.append(dic_dos_sentidos)
    return lista_vecindades_intragrafo, lista_vecindades_dos_sentidos

#funciona

def construir_completacion_grafo(dic_dos_sentidos):
    dic_completacion = dict()
    for nodo in dic_dos_sentidos:
        nodos_vecinos = [(vecino[0],vecino[1]) for vecino in dic_dos_sentidos[nodo]]
        resto_nodos = []
        for otro_nodo in dic_dos_sentidos:
            if otro_nodo != nodo and not otro_nodo in nodos_vecinos:
                distancia = np.linalg.norm(np.array([nodo[0],nodo[1]])-np.array([otro_nodo[0],otro_nodo[1]]))
                resto_nodos.append((otro_nodo[0],otro_nodo[1],distancia))
        dic_completacion[nodo] = resto_nodos
    return dic_completacion

#funciona

def construir_lista_completaciones(lista_dic_dos_sentidos):
    completaciones = []
    for dic in lista_dic_dos_sentidos:
        completaciones.append(construir_completacion_grafo(dic))
    return completaciones

#funciona

def construir_dic_vecinos_intragrilla(grilla):
    fil = len(grilla)
    col = len(grilla[0])
    salidas_por_nodo_intragrilla = dict()
    for ind_fil in range(fil):
        for ind_col in range(col):
            nodo_tupla = (grilla[ind_fil][ind_col][0],grilla[ind_fil][ind_col][1])
            salidas_por_nodo_intragrilla[nodo_tupla] = vecinos_de_nodo_intragrilla(ind_fil,ind_col,fil,col,grilla)
    return salidas_por_nodo_intragrilla

#funciona

#Procedemos a construir las vecindades de grilla completa
#dando la posibilidad de elegir si tomamos distancia 1 o distancia 2
def construir_dic_grilla_completa(grilla,metrica):
    """
    metrica es 1 para Manhattan y 2 para euclídea
    """
    lista_nodos_tupla = []
    filas = len(grilla)
    cols = len(grilla[0])
    dic_intragrilla = dict()
    for fil in range(filas):
        for col in range(cols):
            nodo_tupla = (grilla[fil][col][0],grilla[fil][col][1])
            dic_intragrilla[nodo_tupla] = []
            lista_nodos_tupla.append(nodo_tupla)
    cant_nodos = filas*cols
    for j in range(cant_nodos):
        for k in range(j+1,cant_nodos):
            nodo1 = np.array([lista_nodos_tupla[j][0],lista_nodos_tupla[j][1]])
            nodo2 = np.array([lista_nodos_tupla[k][0],lista_nodos_tupla[k][1]])
            distancia = np.linalg.norm(nodo1-nodo2,ord=metrica)
            dic_intragrilla[lista_nodos_tupla[j]].append((lista_nodos_tupla[k][0],lista_nodos_tupla[k][1],distancia))
            dic_intragrilla[lista_nodos_tupla[k]].append((lista_nodos_tupla[j][0],lista_nodos_tupla[j][1],distancia))
    return dic_intragrilla

#funciona

def diccionario_grilla(grilla,tipo_grilla,metrica):
    """
    tipo_grilla es 'damero' o 'completa'; metrica es 1 o 2
    """
    if tipo_grilla == 'damero':
        return construir_dic_vecinos_intragrilla(grilla)
    elif tipo_grilla == 'completa':
        return construir_dic_grilla_completa(grilla,metrica)
    else:
        print('mal escrito el tipo_grilla, es damero o completa')
        return dict()

def construir_vecindades_grilla_a_grafo(dic_intragrilla,dic_intragrafo):
    conexiones = dict()
    #Sólo los nodos de grilla son claves, así se construye vecindades en un solo sentido
    for nodo_grilla in dic_intragrilla:
        conexiones[nodo_grilla] = []
    for nodo_grilla in dic_intragrilla:
        for nodo_grafo in dic_intragrafo:
            distancia = np.linalg.norm(np.array([nodo_grilla[0],nodo_grilla[1]])-np.array([nodo_grafo[0],nodo_grafo[1]]))
            conexiones[nodo_grilla].append((nodo_grafo[0],nodo_grafo[1],distancia))
    return conexiones

#funciona

def lista_conexiones_grilla_grafos(dic_intragrilla,lista_dic_intragrafo):
    lista_conexiones = []
    for dic_intragrafo in lista_dic_intragrafo:
        lista_conexiones.append(construir_vecindades_grilla_a_grafo(dic_intragrilla,dic_intragrafo))
    return lista_conexiones

#funciona

def construir_diccionario_indices(dic_intragrilla,lista_dic_intragrafo):
    indices = dict()
    nodos_por_indice = []
    ind_actual = 0
    for nodo in dic_intragrilla:
        indices[nodo] = ind_actual
        nodos_por_indice.append(nodo)
        ind_actual += 1
    for dic_intragrafo in lista_dic_intragrafo:
        for nodo in dic_intragrafo:
            indices[nodo] = ind_actual
            nodos_por_indice.append(nodo)
            ind_actual += 1
    return indices,nodos_por_indice

#funciona

def construir_listas_arcos(dic_intragrilla,lista_dic_intragrafo,lista_completaciones,vecindades_grilla_grafo):
    dic_indices,nodos_por_indice = construir_diccionario_indices(dic_intragrilla,lista_dic_intragrafo)
    lista_arcos_dron = []
    lista_arcos_grilla = []
    arcos_por_grafo = []
    #lista de G listas, donde la g-ésima es la lista de arcos del grafo g completo
    listas_pares_por_grafo = []
    #lista de G listas, cada una lista para el grafo g obligatorio los pares de arcos
    #en ambos sentidos de un mismo par de nodos.
    salidas_a_grafos_por_nodo_grilla = []
    llegadas_desde_grafos_por_nodo_grilla = []
    #listas de G diccionarios
    salidas_a_grilla_por_nodo_grilla = dict()
    llegadas_desde_grilla_por_nodo_grilla = dict()
    lista_salidas_por_nodo_de_grafo = []
    lista_llegadas_por_nodo_de_grafo = []
    #Listas de G diccionarios, donde el g-ésimo tiene como claves a los nodos del g-ésimo grafo
    #Reúnen tanto las vecindades de grafo completo como las de grilla
    G = len(lista_dic_intragrafo)
    for nodo in dic_intragrilla:
        salidas_a_grilla_por_nodo_grilla[nodo] = []
        llegadas_desde_grilla_por_nodo_grilla[nodo] = []
    ind_actual_grilla = 0
    for nodo in dic_intragrilla:
        for vecino in dic_intragrilla[nodo]:
            vec = (vecino[0],vecino[1])
            lista_arcos_grilla.append((dic_indices[nodo],dic_indices[vec],vecino[2]))
            salidas_a_grilla_por_nodo_grilla[nodo].append(ind_actual_grilla)
            llegadas_desde_grilla_por_nodo_grilla[vec].append(ind_actual_grilla)
            ind_actual_grilla += 1
    #Completamos las aristas de grilla, pasamos a las del dron
    ind_actual_dron = 0
    #Comenzamos agregando las aristas de recorrido obligatorio
    for g in range(G):
        dic_intragrafo = lista_dic_intragrafo[g]
        salidas_por_nodo = dict() #agregar luego a la lista de diccionarios correspondiente
        llegadas_por_nodo = dict() #agregar luego a la lista de diccionarios correspondiente
        pares_en_grafo = [] #luego agregar a lista_pares_por_grafo
        arcos_en_grafo = [] #luego agregar a arcos_por_grafo
        for nodo in dic_intragrafo:
            salidas_por_nodo[nodo] = []
            llegadas_por_nodo[nodo] = [] 
        for nodo in dic_intragrafo:
            for vecino in dic_intragrafo[nodo]:
                vec = (vecino[0],vecino[1])
                indices_par = []
                lista_arcos_dron.append((dic_indices[nodo],dic_indices[vec],vecino[2]))
                arcos_en_grafo.append(ind_actual_dron)
                indices_par.append(ind_actual_dron)
                salidas_por_nodo[nodo].append(ind_actual_dron)
                llegadas_por_nodo[vec].append(ind_actual_dron)
                ind_actual_dron += 1
                lista_arcos_dron.append((dic_indices[vec],dic_indices[nodo],vecino[2]))
                arcos_en_grafo.append(ind_actual_dron)
                indices_par.append(ind_actual_dron)
                salidas_por_nodo[vec].append(ind_actual_dron) #acá se hace en el otro sentido también
                #porque dic_intragrafo tenía las vecindades en un solo sentido
                llegadas_por_nodo[nodo].append(ind_actual_dron)
                pares_en_grafo.append(indices_par)
                ind_actual_dron += 1
        listas_pares_por_grafo.append(pares_en_grafo)
        dic_completacion = lista_completaciones[g]
        for nodo in dic_completacion:
            for vecino in dic_completacion[nodo]:
                vec = (vecino[0],vecino[1])
                lista_arcos_dron.append((dic_indices[nodo],dic_indices[vec],vecino[2]))
                arcos_en_grafo.append(ind_actual_dron)
                salidas_por_nodo[nodo].append(ind_actual_dron)
                llegadas_por_nodo[vec].append(ind_actual_dron)
                ind_actual_dron += 1
        arcos_por_grafo.append(arcos_en_grafo)
        dic_grilla_grafo = vecindades_grilla_grafo[g]
        salidas_de_grilla = dict()
        llegadas_a_grilla = dict()
        for nodo_grilla in dic_grilla_grafo:
            salidas_de_grilla[nodo_grilla] = []
            llegadas_a_grilla[nodo_grilla] = []
            for vecino_grafo in dic_grilla_grafo[nodo_grilla]:
                vec = (vecino_grafo[0],vecino_grafo[1])
                lista_arcos_dron.append((dic_indices[nodo_grilla],dic_indices[vec],vecino_grafo[2]))
                salidas_de_grilla[nodo_grilla].append(ind_actual_dron)
                llegadas_por_nodo[vec].append(ind_actual_dron)
                ind_actual_dron += 1
                lista_arcos_dron.append((dic_indices[vec],dic_indices[nodo_grilla],vecino_grafo[2]))
                llegadas_a_grilla[nodo_grilla].append(ind_actual_dron)
                salidas_por_nodo[vec].append(ind_actual_dron)
                ind_actual_dron += 1
        salidas_a_grafos_por_nodo_grilla.append(salidas_de_grilla)
        llegadas_desde_grafos_por_nodo_grilla.append(llegadas_a_grilla)
        lista_salidas_por_nodo_de_grafo.append(salidas_por_nodo)
        lista_llegadas_por_nodo_de_grafo.append(llegadas_por_nodo)
    return lista_arcos_dron,lista_arcos_grilla,arcos_por_grafo,listas_pares_por_grafo,salidas_a_grafos_por_nodo_grilla,llegadas_desde_grafos_por_nodo_grilla,salidas_a_grilla_por_nodo_grilla,llegadas_desde_grilla_por_nodo_grilla,lista_salidas_por_nodo_de_grafo,lista_llegadas_por_nodo_de_grafo,dic_indices,nodos_por_indice

#funciona

def problema_discretizado(dic_intragrilla,lista_grafos):
    #lo definimos en función de las vecindades de grilla para que el código
    #no dependa de que el mapa de arcos para la base forme una grilla
    lista_dic_intragrafo, lista_dic_dos_sentidos = construir_lista_vecinos_intragrafos(lista_grafos)
    lista_dic_completaciones = construir_lista_completaciones(lista_dic_dos_sentidos)
    lista_conexiones_grilla_grafo = lista_conexiones_grilla_grafos(dic_intragrilla,lista_dic_intragrafo)
    lista_arcos_dron,lista_arcos_grilla,arcos_por_grafo,listas_pares_por_grafo,salidas_a_grafos_por_nodo_grilla,llegadas_desde_grafos_por_nodo_grilla,salidas_a_grilla_por_nodo_grilla,llegadas_desde_grilla_por_nodo_grilla,lista_salidas_por_nodo_de_grafo,lista_llegadas_por_nodo_de_grafo,dic_indices,nodos_por_indice = construir_listas_arcos(dic_intragrilla,lista_dic_intragrafo,lista_dic_completaciones,lista_conexiones_grilla_grafo)
    return lista_arcos_dron,lista_arcos_grilla,arcos_por_grafo,listas_pares_por_grafo,salidas_a_grafos_por_nodo_grilla,llegadas_desde_grafos_por_nodo_grilla,salidas_a_grilla_por_nodo_grilla,llegadas_desde_grilla_por_nodo_grilla,lista_salidas_por_nodo_de_grafo,lista_llegadas_por_nodo_de_grafo,dic_indices,nodos_por_indice, lista_dic_completaciones, lista_dic_intragrafo

def longitud_total_grafo(grafo_adaptado):
    longitud = 0
    for arco in grafo_adaptado:
        longitud += np.linalg.norm(arco[0]-arco[1])
    return longitud

def modelar_instancia(lista_grafos_importados,tipo_grilla,metrica,filas=11,columnas=11,lim_inf=0,lim_sup=100,lim_izq=0,lim_der=100,orig=(50,50),dest=(50,50),G=5,D=2,v_D=3,v_M=1.5,N_D=40):
    """
    tipo_grilla es 'completa' o 'damero'
    metrica es 1 o 2
    """
    coord_x = np.linspace(lim_izq,lim_der,columnas)
    coord_y = np.linspace(lim_inf,lim_sup,filas)
    grilla = construir_grilla(coord_x,coord_y)
    dic_intragrilla = diccionario_grilla(grilla,tipo_grilla,metrica)
    lista_grafos = gen.crear_lista_grafos_adaptados(lista_grafos_importados)
    lista_arcos_dron,lista_arcos_grilla,arcos_por_grafo,listas_pares_por_grafo,salidas_a_grafos_por_nodo_grilla,llegadas_desde_grafos_por_nodo_grilla,salidas_a_grilla_por_nodo_grilla,llegadas_desde_grilla_por_nodo_grilla,lista_salidas_por_nodo_de_grafo,lista_llegadas_por_nodo_de_grafo,dic_indices,nodos_por_indice, lista_dic_completaciones, lista_dic_intragrafo = problema_discretizado(dic_intragrilla,lista_grafos)
    cant_nodos_grilla = len(dic_intragrilla)
    cant_nodos_por_grafo = [len(dic_intragrafo) for dic_intragrafo in lista_dic_intragrafo]
    cant_aristas_por_grafo_completo = [k*(k-1) for k in cant_nodos_por_grafo]
    aristas_totales_dron = len(lista_arcos_dron)
    aristas_totales_grilla = len(lista_arcos_grilla)
    indice_orig = dic_indices[orig]
    indice_dest = dic_indices[dest]
    T = 2*G
    #empezamos a armar el modelo
    t_const_0 = time.time()
    m = gp.Model("mapa_discretizado")
    #definimos las variables
    Z = m.addMVar(shape=(G,D,aristas_totales_dron),vtype=GRB.BINARY, name="Z_{t,d,arc}") #El primer D es T//2, la cantidad de etapas
    M = m.addMVar(shape=(T,aristas_totales_grilla),vtype=GRB.BINARY, name="M_{t,d,arc}")
    L = m.addMVar(shape=(T,cant_nodos_grilla),vtype=GRB.BINARY, name="L_{t,v}")
    R = m.addMVar(shape=(T,cant_nodos_grilla),vtype=GRB.BINARY, name="R_{t,v}")
    MU = m.addMVar(shape=(T,D,G),vtype=GRB.BINARY,name="Mu_{t,d,g}")    
    B = m.addMVar(shape=T,vtype=GRB.BINARY,name="Beta_t")
    K = m.addMVar(shape=T,vtype=GRB.INTEGER,name="k_t")
    DIST = m.addMVar(shape=T,vtype=GRB.CONTINUOUS,name="d_t")
    """
    Restricciones elementales
    """
    m.addConstrs(L[t,range(cant_nodos_grilla)].sum() == 1 for t in range(T))
    m.addConstrs(R[t,range(cant_nodos_grilla)].sum() == 1 for t in range(T))
    m.addConstr(L[0,indice_orig] == 1) #(3)
    m.addConstr(R[T-1,indice_dest] == 1) #(4)
    m.addConstrs((R[t,:] - L[t+1,:] == 0 for t in range(T-1))) #(5)
    #(6)
    m.addConstrs(Z[k,d,np.array(lista_salidas_por_nodo_de_grafo[g][nodo],dtype=int)].sum()-Z[k,d,np.array(lista_llegadas_por_nodo_de_grafo[g][nodo],dtype=int)].sum() == 0 for g in range(G) for k in range(G) for d in range(D) for nodo in lista_salidas_por_nodo_de_grafo[g])
    for g in range(G): #(7)
        pares = listas_pares_por_grafo[g]
        for par in pares:
            m.addConstr(Z[:,:,par[0]].sum()+Z[:,:,par[1]].sum() >= 1)
    #(8),(9)
    m.addConstrs(cant_aristas_por_grafo_completo[g]*MU[2*k+1,d,g]-Z[k,d,arcos_por_grafo[g]].sum() >= 0 for g in range(G) for k in range(G) for d in range(D))
    m.addConstrs(MU[2*k+1,d,g]-Z[k,d,np.array(arcos_por_grafo[g],dtype=int)].sum() <= 0 for g in range(G) for k in range(G) for d in range(D))
    m.addConstr(sum(MU[2*k,:,:].sum() for k in range(G)) == 0)
    m.addConstrs(MU[:,:,g].sum() == 1 for g in range(G))
    #(12)
    m.addConstrs(MU[2*k+1,d,:].sum() <= 1 for k in range(G) for d in range(D)) #(12)
    distancias_grilla = np.array([lista_arcos_grilla[a][2] for a in range(aristas_totales_grilla)])
    m.addConstrs(sum(M[t,:]*distancias_grilla[:]) - v_M*DIST[t] <= 0 for t in range(T))
    distancias_dron = np.array([lista_arcos_dron[a][2] for a in range(aristas_totales_dron)])
    #(14),(15)
    m.addConstrs(sum(Z[k,d,:]*distancias_dron[:]) <= v_D*N_D for d in range(D) for k in range(G))
    m.addConstrs(sum(Z[k,d,:]*distancias_dron[:]) - v_D*DIST[2*k+1] <= 0 for d in range(D) for k in range(G))
    """
    Restricciones de conectividad para los drones   
    """
    #(16),(17)
    m.addConstrs(sum(Z[k,d,np.array(salidas_a_grafos_por_nodo_grilla[g][nodo],dtype=int)].sum() for nodo in dic_intragrilla) - MU[2*k+1,d,g] == 0 for g in range(G) for d in range(D) for k in range(G))
    m.addConstrs(sum(Z[k,d,np.array(llegadas_desde_grafos_por_nodo_grilla[g][nodo],dtype=int)].sum() for nodo in dic_intragrilla) - MU[2*k+1,d,g] == 0 for g in range(G) for d in range(D) for k in range(G))
    m.addConstrs(Z[k,d,np.array(salidas_a_grafos_por_nodo_grilla[g][nodo],dtype=int)].sum() - L[2*k+1,dic_indices[nodo]] <= 0 for nodo in dic_intragrilla for g in range(G) for k in range(G) for d in range(D))
    m.addConstrs(Z[k,d,np.array(llegadas_desde_grafos_por_nodo_grilla[g][nodo],dtype=int)].sum() - R[2*k+1,dic_indices[nodo]] <= 0 for nodo in dic_intragrilla for g in range(G) for k in range(G) for d in range(D))
    if tipo_grilla == 'damero':
        """
        Vehículo base: grilla damero    
        """
        P = m.addMVar(shape=(T,cant_nodos_grilla),vtype=GRB.INTEGER,name="p_v^t")
        #(20),(21),(22)
        m.addConstrs(M[t,np.array(salidas_a_grilla_por_nodo_grilla[nodos_por_indice[v]],dtype=int)].sum() - M[t,np.array(llegadas_desde_grilla_por_nodo_grilla[nodos_por_indice[v]],dtype=int)].sum() - L[t,v] + R[t,v] == 0 for v in range(cant_nodos_grilla) for t in range(T)) #(20)
        m.addConstrs(P[t,v] >= 0 for v in range(cant_nodos_grilla) for t in range(T))
        m.addConstrs(P[t,v] - M[t,range(aristas_totales_grilla)].sum() <= 0 for v in range(cant_nodos_grilla) for t in range(T))
        #(23)
        m.addConstrs(P[t,lista_arcos_grilla[a][0]]-P[t,lista_arcos_grilla[a][1]]+cant_nodos_grilla*M[t,a] <= cant_nodos_grilla-1 for a in range(aristas_totales_grilla) for t in range(T)) #(23)
    if tipo_grilla == 'completa':
        """
        Vehículo base: grilla completa
        """
        #(24),(25)
        #m.addConstrs(M[t,a]-L[t,lista_arcos_grilla[a][0]]-R[t,lista_arcos_grilla[a][1]] >= -3/2 for t in range(T) for a in range(aristas_totales_grilla)) #(24)
        #m.addConstrs(M[t,a]-0.5*L[t,lista_arcos_grilla[a][0]]-0.5*R[t,lista_arcos_grilla[a][1]] <= 0 for t in range(T) for a in range(aristas_totales_grilla)) #(25)
        m.addConstrs(M[t,np.array(salidas_a_grilla_por_nodo_grilla[nodos_por_indice[v]],dtype=int)].sum() - M[t,np.array(llegadas_desde_grilla_por_nodo_grilla[nodos_por_indice[v]],dtype=int)].sum() - L[t,v] + R[t,v] == 0 for v in range(cant_nodos_grilla) for t in range(T))
        m.addConstrs(M[t,:].sum() <= 1 for t in range(T))
    """
    Reducción de soluciones redundantes
    """
    #(26)
    m.addConstrs(MU[2*k+1,d,:].sum() - MU[2*k+1,d+1,:].sum() >= 0 for d in range(D-1) for k in range(G)) #(26)
    m.addConstrs(K[t] - MU[t,:,:].sum() == 0 for t in range(T)) #(27)
    m.addConstr(B[0]==0) #(28)
    m.addConstrs(G*B[t]-K[:t-1].sum() <= 0 for t in range(1,T)) #(29)
    m.addConstrs(B[t]-K[:t-1].sum() >= 1-G for t in range(1,T)) #(30)
    m.addConstrs(MU[2*k+1,:,:].sum() + B[2*k+1] >= 1 for k in range(G)) #(31)
    """
    Función objetivo
    """
    m.setObjective(quicksum(DIST[t] for t in range(T)),sense=GRB.MINIMIZE)
    t_const_f = time.time()
    tiempo_construccion = t_const_f-t_const_0
    if tiempo_construccion < 1800:
        m.setParam('timeLimit',1800-tiempo_construccion)
    #le damos media hora para construir y resolver
    return m, tiempo_construccion, Z, M, L
#devolviendo las variables que me interesa usar después para que Python no las borre
#después de optimizar














































