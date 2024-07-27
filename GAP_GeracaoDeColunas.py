import time

import mip
from itertools import product
import matplotlib.pyplot as plt
from math import sqrt
import numpy as np
import copy
import numpy.ma as ma
import math
import random

random.seed(0)

class GeracaoColunasGAP():
    def __init__(self, dados):
        self.n_equipes = dados['equipes']
        self.n_tarefas = dados["tarefas"]
        self.pr =  [pr for pr in product(range(self.n_equipes), range(self.n_tarefas))]
        # self.custo_alocacao = {eq: [int(custo_eq) for custo_eq in custos] for eq, custos in
        #                        dados['custo_alocacao_equipes'].items()}

        # self.consumo_cap = {eq: [int(cap_eq) for cap_eq in consumo] for eq, consumo in
        #                     dados['consumo_capacidade'].items()}
        self.custo_alocacao = [int(j) for i in dados['custo_alocacao_equipes'].values() for j in i]
        self.consumo_cap = [int(j) for i in dados['consumo_capacidade'].values() for j in i]
        self.capacidades = [int(c) for c in dados["capacidade_equipe"]]
        self.BigM = 15000
        

    
    def otimiza(self):
        def resolve_PMDW_GAP(colunas, self):
            """
            Método para criar e resolver o problema mestre de Dantizg-Wolfe para o GAP.
            retorna valores dos coeficientes duais da restrições 
            """
            self.model = mip.Model("GAP_GC", sense=mip.MINIMIZE)
            self.model.verbose = 0
            self.var_lambda = {i: self.model.add_var(f'col_{i}', var_type=mip.CONTINUOUS) for i in range(len(colunas))}
            var_atr = self.model.add_var("BigM", var_type=mip.CONTINUOUS)
            self.model.objective = mip.xsum([sum(self.custo_alocacao[pr] * colunas[col][pr] for pr in range(len(colunas[col]))) * self.var_lambda[col] for col in range(1, len(colunas))]) \
                             + (var_atr * self.BigM) + self.var_lambda[0] * 1000
            #restrições
            #
            init = 0
            end = self.n_tarefas - 1
            for eq in range(self.n_equipes):
                self.model += mip.xsum([sum(self.consumo_cap[pr] * colunas[col][pr] for pr in range(len(colunas[col])) if pr >= init and pr <= end) 
                                  * self.var_lambda[col] for col in range(len(colunas))]) <= self.capacidades[eq] + (var_atr * sum(self.capacidades))
                init = end + 1
                end = init +  self.n_tarefas - 1

            self.model += mip.xsum(self.var_lambda[i] for i in range(len(colunas))) + var_atr == 1
            status = self.model.optimize()

            #preciso encontrar valor único de U. Como tenho várias restrições de capacidade (5 por coluna), preciso achar um valor único de U. 
            #Provavelmente vou tentar usar relaxação lagrangeana. Mas para inicio vou considerar que esse U vai ser a soma de todos os Us das colunas.

            
            U0 = self.model.constrs[-1].pi
            #Teste do U sendo a soma de todos enquanto não coloco a relaxação lagrangeana!!!
            U = [self.model.constrs[cr].pi for cr in range(len(self.model.constrs) -1)] 
            print([v.x for v in self.var_lambda.values()])
            return U0, U

        def resolve_SubProblema_GAP_(U, U0, self):
            self.model_sub = mip.Model("Subproblema_GAP_GC", sense=mip.MINIMIZE)
            self.model_sub.verbose = 0
            var_x = {pr: self.model_sub.add_var(f'atr_{pr}', var_type=mip.CONTINUOUS) for pr in self.pr}
            p1 = mip.xsum([self.custo_alocacao[pr] * var_x[self.pr[pr]] for pr in range(len(self.pr))])
            p2 = list()
            for eq in range(self.n_equipes):
                 p2.append(mip.xsum([var_x[(eq, t)] * self.custo_alocacao[self.pr.index((eq,t))] for t in range(self.n_tarefas)]) * int(U[eq]))

            self.model_sub.objective = p1 - mip.xsum(p2) - U0
            
            for tar in range(self.n_tarefas):
                self.model_sub += mip.xsum([var_x[(eq,tar)] for eq in range(self.n_equipes)]) == 1
            
            self.model_sub.optimize()
            col = [v.x for v in var_x.values()]
            return col, self.model_sub.objective_value
        
        colunas = [[0 for pr in self.pr]]
        for tar in range(self.n_tarefas): 
            eq = random.randint(0, self.n_equipes-1)
            pos = self.pr.index((eq,tar))
            colunas[0][pos] = 1
        while True:
            #Criar Problema mestre de dantizg-wolfe que retorne os coeficientes duais das restrições
            U0, U = resolve_PMDW_GAP(colunas=colunas, self=self)
            #Criar subproblema que retorne a coluna (X) e o valor da F.o
            col, FO = resolve_SubProblema_GAP_(U, U0, self)
            #Verificar critérios de parada (custo relativo negativo)
            colunas.append(col)
            print(FO)
            if FO > 0:
                c_final = [0 for _ in self.pr]
                #Achar a Solução ótima através da combinação convexa dos pontos?
                for col in range(len(colunas)):
                    try:
                        c_final += np.array(colunas[col]) *  self.var_lambda[col].x
                    except:
                        c_final += np.array(colunas[col]) *  0
                break
