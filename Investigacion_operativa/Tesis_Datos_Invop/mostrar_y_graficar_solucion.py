import numpy as np
import matplotlib.pyplot as plt
import discretizar_problema as disc

def graficar_linea(fig,ax,nodo_1,nodo_2,color):
    #nodo_1,nodo_2 están dados como tuplas de dos coordenadas
    ax.plot([nodo_1[0],nodo_2[0]],[nodo_1[1],nodo_2[1]],linewidth=0.5,color=color)

def graficar_vecindades_grilla(fig,ax,dic_intragrilla,dic_indices):
    for nodo in dic_intragrilla:
        for vecino in dic_intragrilla[nodo]:
            vec = (vecino[0],vecino[1])
            if dic_indices[nodo] < dic_indices[vec]:
                graficar_linea(fig,ax,nodo,vec,'gray')

def graficar_arista(fig,ax,extremos_arista,color_arista,lim_izq,lim_der,lim_inf,lim_sup):
    long_hori_grilla = lim_der-lim_izq
    long_verti_grilla = lim_sup-lim_inf
    plt.xlim(lim_izq-long_hori_grilla/20, lim_der+long_hori_grilla/20)
    plt.ylim(lim_inf-long_verti_grilla/20, lim_sup+long_verti_grilla/20)
    ax.plot([extremos_arista[0][0],extremos_arista[1][0]],[extremos_arista[0][1],extremos_arista[1][1]],color=color_arista,linewidth=1)

def graficar_lista_puntos(fig,ax,lista_puntos,color_puntos,lim_izq,lim_der,lim_inf,lim_sup):
    long_hori_grilla = lim_der-lim_izq
    long_verti_grilla = lim_sup-lim_inf
    plt.xlim(lim_izq-long_hori_grilla/20, lim_der+long_hori_grilla/20)
    plt.ylim(lim_inf-long_verti_grilla/20, lim_sup+long_verti_grilla/20)
    lista_x = []
    lista_y = []
    for punto in lista_puntos:
        lista_x.append(punto[0])
        lista_y.append(punto[1])
    ax.scatter(lista_x,lista_y,s=10,color=color_puntos)

def graficar_grafos(fig,ax,lista_grafos,lista_colores_grafos,lim_izq,lim_der,lim_inf,lim_sup):
    cant_grafos = len(lista_grafos)
    long_hori_grilla = lim_der-lim_izq
    long_verti_grilla = lim_sup-lim_inf
    plt.xlim(lim_izq-long_hori_grilla/20, lim_der+long_hori_grilla/20)
    plt.ylim(lim_inf-long_verti_grilla/20, lim_sup+long_verti_grilla/20)
    for k in range(cant_grafos):
        for lista_puntos in lista_grafos[k]:
            graficar_lista_puntos(fig,ax,lista_puntos,lista_colores_grafos[k],lim_izq,lim_der,lim_inf,lim_sup)
            graficar_arista(fig,ax,[lista_puntos[0],lista_puntos[-1]],lista_colores_grafos[k],lim_izq,lim_der,lim_inf,lim_sup)

def graficar_completaciones_grafos(fig,ax,lista_dic_completaciones,dic_indices,color_completaciones,lim_izq,lim_der,lim_inf,lim_sup):
    long_hori_grilla = lim_der-lim_izq
    long_verti_grilla = lim_sup-lim_inf
    plt.xlim(lim_izq-long_hori_grilla/20, lim_der+long_hori_grilla/20)
    plt.ylim(lim_inf-long_verti_grilla/20, lim_sup+long_verti_grilla/20)
    for dic_completacion in lista_dic_completaciones:
        for nodo in dic_completacion:
            for vecino in dic_completacion[nodo]:
                vec = (vecino[0],vecino[1])
                if dic_indices[nodo] < dic_indices[vec]:
                    graficar_linea(fig,ax,nodo,vec,color_completaciones)

def listas_variables_activadas(Z,M,L):
    M_activados = []
    Z_activados = []
    L_activados = []
    G,D,aristas_totales_dron = Z.shape
    T,aristas_totales_grilla = M.shape
    etapas,nodos_grilla = L.shape
    for v in range(nodos_grilla):
        for t in range(etapas):
            if L.X[t,v] == 1:
                L_activados.append((t,v))
    for k in range(G):
        for d in range(D):
            for a in range(aristas_totales_dron):
                if Z.X[k,d,a] == 1:
                    Z_activados.append((k,d,a))
    for t in range(T):
        for a in range(aristas_totales_grilla):
            if M.X[t,a] == 1:
                M_activados.append((t,a))
    return Z_activados,M_activados,L_activados

def imprimir_caminos(Z_activados,M_activados,lista_arcos_grilla,lista_arcos_dron,D,G):
    for k in range(G):
        print('base en etapa ',2*k)
        aristas_imprimir = []
        for M_act in M_activados:
            if M_act[0] == 2*k:
                arista_actual = lista_arcos_grilla[M_act[1]]
                aristas_imprimir.append((arista_actual[0],arista_actual[1]))
        print(aristas_imprimir)
        print('base en etapa ',2*k+1)
        aristas_imprimir = []
        for M_act in M_activados:
            if M_act[0] == 2*k+1:
                arista_actual = lista_arcos_grilla[M_act[1]]
                aristas_imprimir.append((arista_actual[0],arista_actual[1]))
        print(aristas_imprimir)
        for d in range(D):
            print('dron ',d,' en etapa ',2*k+1)
            aristas_imprimir = []
            for Z_act in Z_activados:
                if Z_act[0] == k and Z_act[1] == d:
                    arista_actual = lista_arcos_dron[Z_act[2]]
                    aristas_imprimir.append((arista_actual[0],arista_actual[1]))
            print(aristas_imprimir)

def graficar_aristas_recorridas(fig,ax,L_activados,Z_activados,M_activados,lim_izq,lim_der,lim_inf,lim_sup,lista_colores_por_dron,lista_arcos_grilla,lista_arcos_dron,lista_nodos_por_indice,etapa=-1):
    """
    Si etapa = -1, se grafica todo. Si etapa >= 0, se grafica sólo esa etapa.
    """
    lista_partidas = []
    for lanzamiento in L_activados:
        lista_partidas.append(lista_nodos_por_indice[lanzamiento[1]])
    graficar_lista_puntos(fig,ax,lista_partidas,'black',lim_izq,lim_der,lim_inf,lim_sup)
    for M_act in M_activados:
        ind_arista = M_act[1]
        nodo1 = lista_nodos_por_indice[lista_arcos_grilla[ind_arista][0]]
        nodo2 = lista_nodos_por_indice[lista_arcos_grilla[ind_arista][1]]
        t = M_act[0]
        if etapa == t or etapa == -1:
            ax.plot([nodo1[0],nodo2[0]],[nodo1[1],nodo2[1]],linewidth=1,color='black')
    if etapa % 2 == 1: #Acá ya se incluye la opción etapa = -1
        k = (etapa-1) // 2
        for Z_act in Z_activados:
            if etapa == -1:
                k = Z_act[0]
            if Z_act[0] == k:    
                ind_arista = Z_act[2]
                nodo1 = lista_nodos_por_indice[lista_arcos_dron[ind_arista][0]]
                nodo2 = lista_nodos_por_indice[lista_arcos_dron[ind_arista][1]]
                ax.plot([nodo1[0],nodo2[0]],[nodo1[1],nodo2[1]],linestyle='--',linewidth=1,color=lista_colores_por_dron[Z_act[1]])