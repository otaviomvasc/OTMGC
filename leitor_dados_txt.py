import copy
import os
from collections import defaultdict
import tsplib95
from itertools import product
from math import sqrt
import numpy as np

"""
Classes que vão abrir arquivos txt de cada problema e retornar dados para o modelo
em forma de dicionário a principio
"""

class leitor_dados:
    def __init__(self, problema):
        dir_name = os.path.dirname(os.path.realpath(__file__))
        if problema == 'atribuicao_generalizada':
            self.path = dir_name + "\\instancias\\atribuicao_generalizada"

        if problema == 'TSP':
            self.path = dir_name + "\\instancias\\TSP"


        if problema == 'Alocacao':
            self.path = dir_name + "\\instancias\\Alocacao"



    def lista_txt(self):
        lista_arquivos = os.listdir(self.path)
        return lista_arquivos

    def le_dados_atribuicao_gen(self):
        arquivos = self.lista_txt()
        #arquivos = ['gap_6.txt','gap_7.txt']
        dict_end = defaultdict()
        dados_exportacao = defaultdict()
        for arquivo in arquivos:
            count_l = 0
            with open(self.path + "\\" + arquivo, "r") as arq:
                linha = arq.readline()
                if count_l == 0:
                    #primeira linha. Traz a quantidade de problemas no txt!
                    qntd_problemas = int(linha.split()[0])
                    #dict_end['qntd_problemas'] = qntd_problemas
                    count_l += 1
                    n = 1
                    le_equipes_tarefas = True
                    le_custo_alocacao_tarefa_equipe = False
                    le_consumo_capacidade = False
                    capacidade_equipes = False

                while True:
                    linha = arq.readline().split()
                    if len(linha) == 2 and le_equipes_tarefas:
                        #Significa que a posiçao 0 da lista é o numero de equipese 1 é o número de tarefas
                        dict_problemas = defaultdict()
                        dict_problemas["equipes"] = int(linha[0])
                        dict_problemas["tarefas"] = int(linha[1])
                        le_equipes_tarefas = False
                        count_l += 1
                        le_custo_alocacao_tarefa_equipe = True
                        equip = 0
                        list_aux = list()
                        dict_alocacao_equipe = dict()
                        continue

                    if le_custo_alocacao_tarefa_equipe:
                        list_aux.extend(linha)
                        if len(list_aux) == dict_problemas["tarefas"]:
                            dict_alocacao_equipe[equip] = list_aux[:]
                            equip += 1
                            list_aux.clear()
                            if equip == dict_problemas['equipes']:
                                #Significa que já leu o custo de alocacao para todas as equipes!
                                dict_problemas['custo_alocacao_equipes'] = dict_alocacao_equipe
                                le_custo_alocacao_tarefa_equipe = False
                                le_consumo_capacidade = True
                                equip = 0
                                dict_consumo_tarefa = dict()
                                #list_aux.clear()
                                continue

                    if le_consumo_capacidade:
                        list_aux.extend(linha)
                        if len(list_aux) == dict_problemas["tarefas"]:
                            dict_consumo_tarefa[equip] = list_aux[:]
                            equip += 1
                            list_aux.clear()
                            if equip == dict_problemas['equipes']:
                                dict_problemas['consumo_capacidade'] = dict_consumo_tarefa
                                le_consumo_capacidade = False
                                equip = 0
                                lista_capacidade_equipes = list()
                                capacidade_equipes = True
                                continue

                    if capacidade_equipes:
                        lista_capacidade_equipes.extend(linha)
                        capacidade_equipes = False
                        dict_problemas['capacidade_equipe'] = lista_capacidade_equipes[:]
                        lista_capacidade_equipes.clear()
                        dict_end[n] = copy.deepcopy(dict_problemas)
                        n = n+1
                        le_equipes_tarefas = True
                        if n > qntd_problemas:
                            arq.close()
                            break
            dados_exportacao[arquivo] = copy.deepcopy(dict_end)
            dict_end.clear()

        return dados_exportacao

    def le_dados_TSP(self):
        def formata_matriz_diagonal(pr):
            full_list_dist = [dist for _ in pr.edge_weights for dist in _]
            list_aux = list()
            list_fim = list()
            #Formatando na matriz
            for vl in full_list_dist:
                list_aux.append(vl)
                if vl == 0:
                    list_fim.append(list_aux[:])
                    list_aux.clear()

            #formatando os pares do dicionário - Método que vou usar no código!
            #Optei por criar o dicionário assim porque fica mais fácil de trabalhar com problemas assimétricos

            dict_return = dict()
            for pl in product(range(pr.dimension), range(pr.dimension)):
                try:
                    dict_return[pl] = list_fim[pl[0]][pl[1]]
                except IndexError:
                    dict_return[pl] = list_fim[pl[1]][pl[0]]

            return dict_return

        def formata_coords_geograficas(pr):
            nodes = range(pr.dimension)
            dict_c = {par: (pr.node_coords[par[0] + 1][0] - pr.node_coords[par[1] + 1][0]) ** 2 +
                       (pr.node_coords[par[0] + 1][1] - pr.node_coords[par[1] + 1][1]) ** 2
                      for par in product(nodes, nodes)}

            return dict_c

        def formata_full_matrix(pr):
            nodes = range(pr.dimension - 1)
            dict_c = {par: pr.edge_weights[par[0]][par[1]] for par in product(range(pr.dimension),range(pr.dimension))}
            return dict_c

        lista_arquivos = os.listdir(self.path)
        dados = dict()
        for arq in lista_arquivos:
            problem = tsplib95.load(self.path + "//" + arq)
            if problem.edge_weight_type == "EXPLICIT" and problem.edge_weight_format == 'LOWER_DIAG_ROW':
                problem.edge_weights = formata_matriz_diagonal(pr=problem) #Rodar meu problema só a partir dessa matriz não é problemático?
            elif problem.edge_weight_type == "EXPLICIT" and problem.edge_weight_format == 'FULL_MATRIX':
                problem.edge_weights = formata_full_matrix(pr=problem)
            elif problem.edge_weight_type == "EUC_2D":
                problem.edge_weights = formata_coords_geograficas(pr=problem)

            dados[arq] = copy.deepcopy(problem)
        return dados

    def le_dados_alocacao_facilities(self):
        lista_arquivos = os.listdir(self.path)
        dados = dict()
        count_l = 0
        dict_aux = dict()
        for arquivo in lista_arquivos:
            with open(self.path + "\\" + arquivo, "r") as arq:
                linha = arq.readline().split()
                n_plantas = int(linha[0])
                n_clientes = int(linha[1])
                le_capacidade_custo_abr = True
                le_demanda_cliente = False
                le_distancia_planta_cliente = False

                lista_cap = list()
                lista_custo_abertura = list()
                demanda_cliente = list()
                matriz_cliente_planta_aux = list()
                matriz_cliente_planta = list()

                while True:
                    linha = arq.readline().split()

                    if le_capacidade_custo_abr:
                        lista_cap.append(int(linha[0]))
                        lista_custo_abertura.append(int(linha[1][0:len(linha[1]) - 1]))
                        if len(lista_cap) == int(n_plantas):
                            #Significa que já pegou todas as capacidades e custos de abertura
                            le_capacidade_custo_abr = False
                            le_demanda_cliente = True
                            continue

                    if le_demanda_cliente:
                        demanda_cliente.append(linha[0])
                        le_demanda_cliente = False
                        le_distancia_planta_cliente = True
                        continue

                    if le_distancia_planta_cliente:
                        matriz_cliente_planta_aux.extend(linha)
                        if len(matriz_cliente_planta_aux) == int(n_plantas):
                            matriz_cliente_planta.append(matriz_cliente_planta_aux[:])
                            matriz_cliente_planta_aux.clear()
                            if len(matriz_cliente_planta) == int(n_clientes):
                                arq.close()
                                break
                            le_distancia_planta_cliente = False
                            le_demanda_cliente = True

                mp_transp = np.transpose(matriz_cliente_planta)
                dados[arquivo] = {
                    "plantas" : n_plantas,
                    "clientes" : n_clientes,
                    "capacidade_plantas" : {p: lista_cap[p] for p in range(n_plantas)},
                    "custo_abertura" : {p: lista_custo_abertura[p] for p in range(n_plantas)},
                    "demanda_clientes" : {c: int(demanda_cliente[c]) for c in range(n_clientes)},
                    "matriz_custo" : {(p_c): float(mp_transp[p_c[0],p_c[1]] )for p_c in product(range(n_plantas), range(n_clientes))}
                }

        return dados




if __name__ == '__main__':
    prob = leitor_dados(problema='Alocacao').le_dados_alocacao_facilities() #Criar método para ler ddados e esse método chama a leitura de dados internos!
