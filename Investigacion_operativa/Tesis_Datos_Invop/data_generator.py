"""
Created on Thu Nov  5 11:31:51 2020

@author: rafaelblanquero
"""

# import ujson
import pickle

from neighbourhood import *
from data import Data

r = 1
alpha_list = [True, False]
endurances = [20, 30, 40, 50, 60]
drones = [1, 2, 3]
initialization = True    

def generar_instancias():
    seed = np.random.seed(2)
    list_1 = [4, 6, 8, 10, 12]
    list_2 = [4, 4, 6, 6, 8, 8, 10, 10, 12, 12]
    lists = [list_1, list_2]
    instances = {}
    grid = True
    time_limit = 3600
    for lista, j in zip(lists, range(len(lists))):
        for i in range(5):
            m = len(lista)
            grid = True
            data = Data([], graphs_number=m, grid_mode=grid, time_limit=time_limit, alpha=True, fleet_size = 1, time_endurance=20, seed=seed, drone_speed = 3, truck_speed = 1)
            data.generate_grid()
            data.generate_graphs(lista)
            instances[(i, len(lista))] = data
    return instances

def adaptar_grafo_importado(grafo_APV):
    lista_nodos = grafo_APV.V
    cant_nodos = len(lista_nodos)
    matriz_adyacencias = grafo_APV.A
    lista_arcos = []
    for i in range(cant_nodos):
        for j in range(i,cant_nodos):
               if matriz_adyacencias[i,j] > 0: #esto indica que la arista existe
                    lista_arcos.append([lista_nodos[i],lista_nodos[j]])
    return lista_arcos

def crear_lista_grafos_adaptados(lista_grafos_importados):
    lista_grafos = []
    for grafo_APV in lista_grafos_importados:
        lista_grafos.append(adaptar_grafo_importado(grafo_APV))
    return lista_grafos

def crear_lista_instancias_adaptadas(dic_instancias):
    nro_instancia = [0,1,2,3,4]
    cant_grafos = [5,10]
    instancias_adaptadas = dict()
    for G in cant_grafos:
        for inst in nro_instancia:
            instancias_adaptadas[(inst,G)] = crear_lista_grafos_adaptados(dic_instancias[(inst,G)].instances)
    return instancias_adaptadas
