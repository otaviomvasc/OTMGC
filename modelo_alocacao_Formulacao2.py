import mip
from itertools import product
import matplotlib.pyplot as plt
from math import sqrt
import numpy as np

"""
Código para modelo de alocação de facilites de acordo com a segunda formulação.
Variável Z:

Z (cliente, posição i no vetor de distâncias do cliente i até as plantas P ordenado do menor para o maior.
Criar matriz D por clientes com lista ordenada.

"""

#Pré-OTM
c_instalacao = [430, 452, 382, 242, 366, 259, 449, 500, 371, 297]
pos_clientes = np.array([(25, 36), (27, 7), (3, 23), (17, 28), (33, 21), (32, 40), (26, 41), (20, 14), (31, 49), (23,
31), (38, 29), (14, 34), (45, 17), (9, 4), (19, 1), (36, 6), (7, 26), (46, 35), (30, 20), (18, 16), (40,
27), (43, 46), (42, 11), (5, 8), (10, 24), (4, 39), (47, 47), (22, 3), (11, 48), (16, 19)])

pos_plantas = np.array([(19, 24), (20, 40), (14, 23), (39, 12), (33, 43), (10, 26), (41, 18), (25, 42), (36, 33), (8, 7)])
demanda = [10, 7, 6, 10, 4, 5, 3, 4, 3, 1, 10, 5, 8, 2, 2, 3, 3, 1, 2, 9, 7, 9, 5, 9, 4, 4, 10, 7, 10, 5]

n_plantas = len(c_instalacao)
n_clientes = len(demanda)

plantas = range(n_plantas)
clientes = range(n_clientes)



#Ao invés de usar matriz vou usar dicionário. É mais lento, porém acredito que código fica mais fácil de padronizar!
m_custo = {(j, i): sqrt((pos_clientes[i][0] - pos_plantas[j][0]) ** 2 + (pos_clientes[i][1] - pos_plantas[j][1]) ** 2)
            for j in plantas for i in clientes}

distancias_d = dict()
distancias_dict_valores = dict()
for cliente in clientes:
    dict_aux = {cb: m_custo[cb] for cb in m_custo if cb[1] == cliente}
    distancias_d[cliente] = [i[0] for i in sorted(dict_aux, key=dict_aux.get, reverse=False)]
    distancias_dict_valores[cliente] = {i:dict_aux[i] for i in sorted(dict_aux, key=dict_aux.get, reverse=False) }


#chamada do modelo
model = mip.Model('Problema_Localizacao', mip.MINIMIZE)

#Criação das variáveis
#Variavel binária de escolha de planta

#teste para ver semelhança com o pulp!
var_atv_planta = {planta: model.add_var(var_type=mip.BINARY, name=f'atv_planta_{planta}') for planta in plantas}

#Variável da porcentagem da demanda do cliente C atendida pela planta P
var_atend_demanda = {(planta_cliente): model.add_var(var_type=mip.CONTINUOUS, ub=1,
                                                     name=f"% da demanda da planta {planta_cliente[0]} - atendendo cliente {planta_cliente[1]}")
                for planta_cliente in product(plantas, clientes)}

#Variavel z - planta localizada até a distancia máxima distancias_d(cliente, facilities)7
dict_teste = dict()
z = dict()
for cliente in distancias_d:
    z[cliente] = [model.add_var(var_type=mip.BINARY, name=f'planta_ate {(cliente, planta)}') for planta
                  in distancias_d[cliente]]


model.objective = mip.minimize(
    #Custo de ativação de planta
    mip.xsum(var_atv_planta[pl] * c_instalacao[pl] for pl in plantas)
    #custo de transporte!
    + mip.xsum(
        (m_custo[(distancias_d[cliente][0], cliente)] +
        #somatório de diferenças consecutivas!

        sum(m_custo[(distancias_d[cliente][k+1], cliente)] - m_custo[(distancias_d[cliente][k], cliente)]#TODO: Porque se usar o mip.xsum aqui da erro?
         for k in range(len(distancias_d[cliente])-1) if k < n_plantas) #somátorio de todos os k.Confirmar se é isso mesmo!
        )#Fim do somatorio de dentro
        * z[cliente][k1] #multiplicação pela variável Z
        for cliente in clientes for k1 in range(len(distancias_d[cliente])-1) if k1 < n_plantas
    )
    )


#A planta escolhida vai zerar a variável z posição 0 daquele cliente!
for cliente in clientes:
    model += (
            z[cliente][0] +
            mip.xsum(var_atv_planta[pl] for pl in distancias_d[cliente]) >= 1
    )



for cliente in clientes:
    for k in plantas:
        if k == 0:
            continue
        model += (
                z[cliente][k] +
                mip.xsum(var_atv_planta[pl] for pl in distancias_d[cliente]) >= z[cliente][k-1]
        )


#Chamadado Solver
status = model.optimize()

list_z = list()
for c in clientes:
    for vari in z[c]:
        if vari.x > 0:
            list_z.append(vari)
#Pós-Otimização!

if status == mip.OptimizationStatus.OPTIMAL:
    print(f'Custo de Instalação = {sum(var_atv_planta[p].x * c_instalacao[p] for p in plantas)}')
    print(f'Custo de Transporte = {round(sum(var_atend_demanda[planta_cliente].x * demanda[planta_cliente[1]] * m_custo[planta_cliente] for planta_cliente in product(plantas, clientes)),2)}')
    print(f'Custo Total = {model.objective_value}')

    fig, ax = plt.subplots()

    #marcando os clientes
    plt.scatter(pos_clientes[:, 0], pos_clientes[:, 1], marker="o", color='black', s=10, label="clientes")

    for i in clientes:
        plt.text(pos_clientes[i,0], pos_clientes[i,1], "{:d}".format(i + 1))

    for (i, j) in [(i, j) for (i, j) in product(plantas, clientes) if var_atend_demanda[(i, j)].x >= 1e-6]:
        plt.plot((pos_clientes[j][0], pos_plantas[i][0]), (pos_clientes[j][1], pos_plantas[i][1]), linestyle="--", color="black")

    plt.scatter(pos_plantas[:, 0], pos_plantas[:, 1], marker="^", color='black', s=100, label="plantas")

    for j in plantas:
        plt.text(pos_plantas[j][0] + .5, pos_plantas[j][1], "{:d}".format(j+1))

    plt.legend()
    plt.plot()
    plt.show()

