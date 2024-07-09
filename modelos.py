import time

import mip
from itertools import product
import matplotlib.pyplot as plt
from math import sqrt
import numpy as np
import copy
import numpy.ma as ma
import math

class ModelAtribuicaoGeneralizada:
    def __init__(self, dados):
        self.n_equipes = dados['equipes']
        self.n_tarefas = dados["tarefas"]
        self.custo_alocacao = {eq: [int(custo_eq) for custo_eq in custos] for eq, custos in  dados['custo_alocacao_equipes'].items()}
        self.consumo_cap = {eq: [int(cap_eq) for cap_eq in consumo] for eq,consumo in dados['consumo_capacidade'].items()}
        self.capacidades = [int(c) for c in dados["capacidade_equipe"]]

    def otimiza(self):
        modelo = mip.Model('Problema_Atribuicao', mip.MINIMIZE)
        var_atr_equipes = {(eq, tarefa): modelo.add_var(f'alocacao equipe {eq} tarefa {tarefa}', var_type=mip.BINARY)
                           for eq, tarefa in product(range(self.n_equipes), range(self.n_tarefas))}

        #Função objetivo de custo!
        modelo.objective = mip.xsum(
            var_atr_equipes[(eq, tarefa)] * self.custo_alocacao[eq][tarefa]
            for eq, tarefa in product(range(self.n_equipes), range(self.n_tarefas))
        )

        #Restrições:
        #Todas as tarefas precisam ser atendidas
        for tarefa in range(self.n_tarefas):
            modelo += mip.xsum(var_atr_equipes[(eq, tarefa)] for eq in range(self.n_equipes)) == 1


        #Capacidade das equipes deve ser respeitada!
        for equips in range(self.n_equipes):
            modelo += mip.xsum(var_atr_equipes[(equips, tar)] * self.consumo_cap[equips][tar] for tar in range(self.n_tarefas)) <= self.capacidades[equips]

        init = time.time()
        #chama solve!
        status = modelo.optimize()
        end = time.time()

        if status == mip.OptimizationStatus.OPTIMAL:
            for eq, tarefa in product(range(self.n_equipes), range(self.n_tarefas)):
                if var_atr_equipes[(eq, tarefa)].x > 0:
                    print(var_atr_equipes[(eq, tarefa)].name)

        #Monta exportação!
        valor_fo = modelo.objective_value
        tamanho_instacia = f'{self.n_equipes} e {self.n_tarefas}'
        tempo_gasto =end-init

        return {"valor_fo": valor_fo, "tamanho_instacia": tamanho_instacia, "tempo": tempo_gasto }

class ModelTSP:
    def __init__(self, dados):
        self.n = dados.dimension
        self.matriz_c = dados.edge_weights

    def otimiza(self, subrota_TMZ=False, subrota_fluxo_ficticio=False, sub_rota_iterativa=True, max_seconds = 3600):
        def rest_subrota_TMZ(model, max_seconds):
            var_seq = {node: model.add_var(name=f'seq {node}', var_type=mip.INTEGER)
                       for node in range(self.n)}

            for n,j in product(range(1, self.n), range(1, self.n)):
                if n == j:
                    continue
                model += var_seq[n] - var_seq[j] + self.n * var_arco[(n,j)] <= self.n - 1

            status = model.optimize(max_seconds=max_seconds)

        def rest_sub_rota_iterativa(model, max_seconds):
            """
                Pseudocódigo para retirada subrota por restrição de fluxo
                Roda solve
                while true:
                    Pega todos os pontos visitados
                    verifica se há subrota
                    se sim:
                        pega arcos de cada subrota
                        cria restrição onde a soma dos arcos da subrota precisa ser menor que o tamanho desse conjunto -1.
                        ex:
                        Subrota com os pontos (1,3,5):
                        sum (x13 + x35 + x51 <= len(1,3,5) - 1

                    senão:
                        retorna resultado porque não houve subrota!
                        sai do while

            """

            while True:
                model.optimize(max_seconds=max_seconds)
                #checa se há subrota:
                pt = 0
                rotas_totais = list()
                rotas_aux = list()
                while True:
                    n_arc = next(v for v in var_arco if v[0] == pt and var_arco[v].x > 0)
                    if n_arc[1] in [arc[0] for arc in rotas_aux]:
                        rotas_aux.append(n_arc)
                        rotas_totais.append(rotas_aux[:])
                        if sum(len(sub_rotas) for sub_rotas in rotas_totais) == self.n:
                            #todas as subrotas já foram mapeadas!
                            rotas_aux.clear()
                            break

                        pt = next(n for n in range(self.n) if n not in [b[1] for a in rotas_totais for b in a])
                        rotas_aux.clear()
                        #rotas_aux.append(pt)
                        continue
                        # Fecha subrota!!!
                    rotas_aux.append(n_arc)
                    pt = copy.deepcopy(n_arc[1])

                if len(rotas_totais[0]) == self.n: #Não gerou subrotas
                    break

                #Cria restrição de fluxo para cada subrota!
                for sub_r in rotas_totais:
                   model += mip.xsum(var_arco[p] for p in sub_r) <= len(sub_r) - 1

        def rest_subrota_fluxo_ficticio(model, max_seconds):
            #Criar variável de fluxo (i,j)
            var_fluxo = {pr: model.add_var(name=f'fluxo {pr}', var_type=mip.INTEGER)
                         for pr in product(range(self.n), range(self.n)) if pr[0] != pr[1]}

            #restrição que garanta que o veículo saia da planta com n-1 unidades de fluxo.
            model += (mip.xsum(var_fluxo[(0, j)] for j in range(1, self.n)) == self.n-1)


            #Restrição de balanço de fluxo:
            for node in range(1, self.n):
                model += (
                    mip.xsum(var_fluxo[(i, node)] for i in range(self.n) if i != node) - mip.xsum(var_fluxo[(node, i)] for i in range(self.n) if i != node) == 1
                )

            #Ativação do fluxo x(i,j)
            for pr in product(range(self.n), range(self.n)):
                if pr[0] == pr[1]:
                    continue
                model += (
                    var_fluxo[pr] <= (self.n -1) * var_arco[pr]
                )

            status = model.optimize(max_seconds=max_seconds)


        model = mip.Model('TSP', mip.MINIMIZE)
        var_arco = {pr:model.add_var(name=f'arco {pr[0]} - {pr[1]}', var_type=mip.BINARY)
                    for pr in product(range(self.n), range(self.n)) if pr[0] != pr[1]}


        model.objective = mip.minimize(
            mip.xsum(var_arco[pr] * self.matriz_c[pr] for pr in product(range(self.n), range(self.n)) if pr[0] != pr[1])
        )

        #Restrições:
        #Cada nó é origem uma vez!
        init = time.time()
        for n in range(self.n):
            model += mip.xsum(var_arco[(n,j)] for j in range(self.n) if n != j) == 1

        #Cada nó é destino uma vez!
        for n in range(self.n):
            model += mip.xsum(var_arco[(j,n)] for j in range(self.n) if n != j) == 1

        if subrota_TMZ:
            #Restrição MTZ para retirada de sub-rota. Considerei nó 0 como origem!
            rest_subrota_TMZ(model=model, max_seconds=max_seconds)
            #return model.objective_value

        elif subrota_fluxo_ficticio:
            rest_subrota_fluxo_ficticio(model=model, max_seconds=max_seconds)
            #return model.objective_value

        elif sub_rota_iterativa():
            rest_sub_rota_iterativa(model=model, max_seconds=max_seconds)
            #return model.objective_value

        end = time.time()
        valor_fo = model.objective_value
        tamanho_instacia = f'{self.n}'
        tempo_gasto = end - init

        return {"valor_fo": valor_fo, "tamanho_instacia": tamanho_instacia, "tempo": tempo_gasto}

class AlocacaoFacilities():
    def __init__(self, dados):
        self.c_instalacao = dados['custo_abertura']
        self.demanda = dados['demanda_clientes']
        self.n_plantas = dados['plantas']
        self.n_clientes = dados['clientes']
        self.m_custo = dados['matriz_custo']
        self.plantas = range(self.n_plantas)
        self.clientes = range(self.n_clientes)


    def otimiza(self):
        distancias_d = dict()
        distancias_dict_valores = dict()
        for cliente in self.clientes:
            dict_aux = {cb: self.m_custo[cb] for cb in self.m_custo if cb[1] == cliente}
            distancias_d[cliente] = [i[0] for i in sorted(dict_aux, key=dict_aux.get, reverse=False)]
            distancias_dict_valores[cliente] = {i: dict_aux[i] for i in
                                                sorted(dict_aux, key=dict_aux.get, reverse=False)}

        # chamada do modelo
        model = mip.Model('Problema_Localizacao', mip.MINIMIZE)

        # Criação das variáveis
        # Variavel binária de escolha de planta

        var_atv_planta = {planta: model.add_var(var_type=mip.BINARY, name=f'atv_planta_{planta}') for planta in self.plantas}


        # Variavel z - planta localizada até a distancia máxima distancias_d(cliente, facilities)7
        dict_teste = dict()
        z = dict()
        for cliente in distancias_d:
            z[cliente] = [model.add_var(var_type=mip.CONTINUOUS, lb=0, name=f'planta_ate {(cliente, planta)}') for
                          planta
                          in distancias_d[cliente]]

        model.objective = mip.minimize(
            # Custo de ativação de planta
            mip.xsum(var_atv_planta[pl] * self.c_instalacao[pl] for pl in self.plantas) +
            # custo de transporte!
            mip.xsum(
                (self.m_custo[(distancias_d[cliente][0], cliente)] +
                 (
                     mip.xsum((self.m_custo[(distancias_d[cliente][k + 1], cliente)] - self.m_custo[
                         (distancias_d[cliente][k], cliente)]) * z[cliente][k]
                              for k in range(len(distancias_d[cliente]) - 1) if k < self.n_plantas)
                 )
                 ) #* self.demanda[cliente]

                for cliente in self.clientes
            )
        )


        # A planta escolhida vai zerar a variável z posição 0 daquele cliente!
        for cliente in self.clientes:
            model += (z[cliente][0] + var_atv_planta[distancias_d[cliente][0]] >= 1)

        for cliente in self.clientes:
            for k in self.plantas:
                if k == 0:
                    continue

                model += (z[cliente][k] + var_atv_planta[distancias_d[cliente][k]] >= z[cliente][k - 1])

        # Chamadado Solver]
        init = time.time()
        status = model.optimize()

        end = time.time()
        valor_fo = model.objective_value
        tamanho_instacia = f'{self.n_plantas} x {self.n_clientes}'
        tempo_gasto = end - init

        return {"valor_fo": valor_fo, "tamanho_instacia": tamanho_instacia, "tempo": tempo_gasto}

class AlocacaoFacilitiesPrimal():
    """
    CLasse utilizada para calcular o valor do subproblema primal de alocação de facilities com Y fixado
    para garantir que o desenvolvimento do meu dual estará correto, já que ambos precisam ter a mesma F.O
    dado o mesmo Y
    """

    def __init__(self, dados):
        self.c_instalacao = dados['custo_abertura']
        self.demanda = dados['demanda_clientes']
        self.n_plantas = dados['plantas']
        self.n_clientes = dados['clientes']
        self.m_custo = dados['matriz_custo']
        self.plantas = range(self.n_plantas)
        self.clientes = range(self.n_clientes)


    def otimiza(self):
        distancias_d = dict()
        distancias_dict_valores = dict()
        for cliente in self.clientes:
            dict_aux = {cb: self.m_custo[cb] for cb in self.m_custo if cb[1] == cliente}
            distancias_d[cliente] = [i[0] for i in sorted(dict_aux, key=dict_aux.get, reverse=False)]
            distancias_dict_valores[cliente] = {i: dict_aux[i] for i in
                                                sorted(dict_aux, key=dict_aux.get, reverse=False)}

        # chamada do modelo
        model = mip.Model('Problema_Localizacao', mip.MINIMIZE)

        # Criação das variáveis
        # Variavel binária de escolha de planta
        lista_y_inteiros = [1,4,8]
        var_atv_planta = {planta: (1 if planta in lista_y_inteiros else 0) for planta in self.plantas}


        # Variavel z - planta localizada até a distancia máxima distancias_d(cliente, facilities)7
        dict_teste = dict()
        z = dict()
        for cliente in distancias_d:
            z[cliente] = [model.add_var(var_type=mip.CONTINUOUS, lb=0, name=f'planta_ate {(cliente, planta)}') for
                          planta
                          in distancias_d[cliente]]

        model.objective = mip.minimize(
            # Custo de ativação de planta - Retirada da atv da planta!
            #mip.xsum(var_atv_planta[pl] * self.c_instalacao[pl] for pl in self.plantas) +
            # custo de transporte!
            mip.xsum(
                (self.m_custo[(distancias_d[cliente][0], cliente)] +
                 (
                     mip.xsum((self.m_custo[(distancias_d[cliente][k + 1], cliente)] - self.m_custo[
                         (distancias_d[cliente][k], cliente)]) * z[cliente][k]
                              for k in range(len(distancias_d[cliente]) - 1) if k < self.n_plantas)
                 )
                 ) #* self.demanda[cliente]

                for cliente in self.clientes
            )
        )


        # A planta escolhida vai zerar a variável z posição 0 daquele cliente!
        for cliente in self.clientes:
            model += (z[cliente][0] + var_atv_planta[distancias_d[cliente][0]] >= 1)

        for cliente in self.clientes:
            for k in self.plantas:
                if k == 0:
                    continue

                model += (z[cliente][k] + var_atv_planta[distancias_d[cliente][k]] >= z[cliente][k - 1])

        # Chamadado Solver]
        init = time.time()
        status = model.optimize()

        end = time.time()
        valor_fo = model.objective_value
        tamanho_instacia = f'{self.n_plantas} x {self.n_clientes}'
        tempo_gasto = end - init

        return {"valor_fo": valor_fo, "tamanho_instacia": tamanho_instacia, "tempo": tempo_gasto}

class AlocacaoFacilitiesDual():
    """
    CLasse utilizada para calcular o valor do subproblema DUAL de alocação de facilities com Y fixado nos mesmos valores do primal
    para garantir que o desenvolvimento do meu dual estará correto, já que ambos precisam ter a mesma F.O
    dado o mesmo Y
    """

    def __init__(self, dados):
        self.c_instalacao = dados['custo_abertura']
        self.demanda = dados['demanda_clientes']
        self.n_plantas = dados['plantas']
        self.n_clientes = dados['clientes']
        self.m_custo = dados['matriz_custo']
        self.plantas = range(self.n_plantas)
        self.clientes = range(self.n_clientes)

    def otimiza(self):
        distancias_d = dict()
        distancias_dict_valores = dict()
        for cliente in self.clientes:
            dict_aux = {cb: self.m_custo[cb] for cb in self.m_custo if cb[1] == cliente}
            distancias_d[cliente] = [i[0] for i in sorted(dict_aux, key=dict_aux.get, reverse=False)]
            distancias_dict_valores[cliente] = {i: dict_aux[i] for i in
                                                sorted(dict_aux, key=dict_aux.get, reverse=False)}

        # chamada do modelo
        model = mip.Model('Problema_Localizacao', mip.MINIMIZE)

        # Criação das variáveis
        # Variavel binária de escolha de planta
        lista_y_inteiros = [1, 4, 8]
        var_atv_planta = {planta: (1 if planta in lista_y_inteiros else 0) for planta in self.plantas} #Variável Y nos artigos
        v = dict()
        # for cliente in self.clientes:
        #     v[cliente] = [model.add_var()]

class EgapLangrange():
    def __init__(self, dados):
        self.n_equipes = dados['equipes']
        self.n_tarefas = dados["tarefas"]
        self.custo_alocacao = {eq: [int(custo_eq) for custo_eq in custos] for eq, custos in
                               dados['custo_alocacao_equipes'].items()}
        self.consumo_cap = {eq: [int(cap_eq) for cap_eq in consumo] for eq, consumo in
                            dados['consumo_capacidade'].items()}
        self.capacidades = [int(c) for c in dados["capacidade_equipe"]]

    def otimiza(self):
        def calcula_atr_menor_custo(u_, custo ):
            var_y_ij = {p: 0 for p in product(range(self.n_equipes), range(self.n_tarefas))}
            for i in range(self.n_equipes):
                v_custo_tarefa = [self.custo_alocacao[i][j] * self.beta + u_[i, j] for j in range(self.n_tarefas)]
                j = np.argmin(v_custo_tarefa)
                var_y_ij[(i, j)] = 1
                # Calculo do Custo - Poderia fazer de uma vez só, mas vou devagar para entender melhor o passo a passo - Jogar para método depois!
            custo_prob_y = sum([(self.custo_alocacao[i][j] * self.beta + u_[i, j]) * var_y_ij[i, j] for i in range(self.n_equipes) for j in range(self.n_tarefas)])
            return var_y_ij, custo_prob_y

        def calcula_atr_menor_custo_por_tarefa(u_, custo ):
            var_y_ij_x = {p: 0 for p in product(range(self.n_equipes), range(self.n_tarefas))}
            for tar in range(self.n_tarefas):
                v_custo_tarefa = [custo[tar][i] * self.beta + u_[i, tar] for i in range(self.n_equipes)]
                eq = np.argmin(v_custo_tarefa)
                var_y_ij_x[(eq, tar)] = 1
                # Calculo do Custo - Poderia fazer de uma vez só, mas vou devagar para entender melhor o passo a passo - Jogar para método depois!
            custo_prob_y = sum([(self.custo_alocacao[i][j] * self.beta + u_[i, j]) * var_y_ij_x[i, j] for i in range(self.n_equipes) for j in range(self.n_tarefas)])
            return var_y_ij_x, custo_prob_y

        def calcula_knapsack_por_recurso(_u, eq):
            model = mip.Model(name="Knapsack", sense=mip.MINIMIZE)
            model.verbose = 0
            var_x = {tar: model.add_var(name=f'tar {tar}', var_type=mip.BINARY) for tar in range(self.n_tarefas)}
            model.objective = mip.xsum((self.alfa * self.custo_alocacao[eq][tar] - _u[(eq, tar)]) * var_x[tar] for tar in range(self.n_tarefas))  # TODO: não seja burro e use matrizes para custo!!!
            model += mip.xsum(self.consumo_cap[eq][tar] * var_x[tar] for tar in range(self.n_tarefas)) <= self.capacidades[eq]
            model.optimize()
            obj = model.objective_value
            return obj, var_x

        def calcula_knapsack_por_recurso_usando_obj_add_var(_u, eq):
            model = mip.Model(name="Knapsack", sense=mip.MINIMIZE)
            model.verbose = 0
            var_x = {tar: model.add_var(name=f'tar {tar}', var_type=mip.BINARY, obj=self.alfa * self.custo_alocacao[eq][tar] - _u[(eq, tar)])
                     for tar in range(self.n_tarefas)}
            # Tem diferença em usar o obj do add_var ou a função objetivo? - Testar
            #model.objective = mip.xsum((self.alfa * self.custo_alocacao[eq][tar] - _u[(eq, tar)]) * var_x[tar] for tar in
                                       #range(self.n_tarefas))
            model += mip.xsum(self.custo_alocacao[eq][tar] * var_x[tar] for tar in range(self.n_tarefas)) <= \
                     self.capacidades[eq]
            st = model.optimize()
            obj = model.objective_value
            return obj, var_x

        def print_stats(h, ub, infimum, lb, gap, mu, norm):
            print("{:4d} ".format(h), end='')
            print("{:12,.2f} ".format(ub), end='')
            print("{:12,.2f} ".format(infimum), end='')
            print("{:12,.2f} ".format(lb), end='')
            print("{:12,.2f} %".format(gap), end='')
            print("{:12,.8f} ".format(mu), end='')
            print("{:12,.2f} ".format(norm, end=''))

        def resolve_s_inicial():
            model = mip.Model(name='U inicial', sense=mip.MINIMIZE)
            model.verbose = 0
            var_x = {p: model.add_var(name=f'x{p}', var_type=mip.BINARY) for p in product(range(self.n_equipes),range(self.n_tarefas))}
            model.objective = mip.xsum(self.custo_alocacao[p[0]][p[1]] * var_x[p] for p in product(range(self.n_equipes),range(self.n_tarefas)))
            for j in range(self.n_equipes):
                model += mip.xsum(var_x[i,j] for i in range(self.n_equipes)) == 1

            model.optimize()
            return {p: var_x[p].x for p in product(range(self.n_equipes),range(self.n_tarefas))}

        #Iniciar variáveis de controle e bounds!
        lb = 0
        custo_por_tarefa = {tarefa: [self.custo_alocacao[i][tarefa] for i in range(self.n_equipes)] for tarefa in range(self.n_tarefas) } #Criar heurística para criar um UB melhor!
        ub = sum(max(custo_por_tarefa[tarefa]) for tarefa in range(self.n_tarefas))
        same = 0
        same_limit = 15
        h_limite = 200
        mu=2
        #u_new = {p: 0 for p in product(range(self.n_equipes), range(self.n_tarefas))} # Criar U
        u_new = resolve_s_inicial()

        h = 0
        gap = 100.0
        self.alfa = .5
        self.beta = .5

        while True:
            #Resolve problema Y - Argmin
            u = copy.deepcopy(u_new)
            #Calculo do Y
            Y_ij, FO_y = calcula_atr_menor_custo_por_tarefa(u_=u, custo=custo_por_tarefa)

            #Resolve Problema X - Quantas atividades a equipe I pode receber por tarefa J?
            FO_x = 0
            var_x_ij = {p: 0 for p in product(range(self.n_equipes), range(self.n_tarefas))}
            for eq in range(self.n_equipes):
                #Resolver problema por método exato!!

                FO_x_r, x_ = calcula_knapsack_por_recurso(_u=u, eq=eq)
                FO_x += FO_x_r
                for j in range(self.n_tarefas):
                    var_x_ij[eq,j] = x_[j].x

            #Vetor que junta as duas soluções é y-x segundo a restrição e o artigo!!
            vetor_u_final = {p: Y_ij[p] - var_x_ij[p] for p in product(range(self.n_equipes), range(self.n_tarefas))}
            FO_t = FO_x + FO_y



            if lb < FO_t:
                lb = FO_t
                same = 0 #Iterações sem melhora!!
            else:
                same += 1

            if (same > same_limit):
                mu = mu / 2.0 #máximo de iterações sem melhora eu divido o parametro de convergencia
                same = 0


            #Calculo da norma do U-final
            norm = sum(vetor_u_final[p] ** 2 for p in product(range(self.n_equipes), range(self.n_tarefas))) #Tiro a raiz ou não??

            h += 1
            if (norm < 0.001) or (mu < 0.0005) or (h > h_limite) or (gap < 0.0001):
            #if ([float(i) for i in Y_ij.values()] == var_x_ij.values()) or (h > h_limite) or (mu < 0.0005):
                break

            step = mu * (ub- lb) / norm
            gap = 100.0 * (ub - lb) / ub
            u_new = {p: u[p] + (step * vetor_u_final[p]) for p in product(range(self.n_equipes), range(self.n_tarefas))}
            ax = sum(abs(u_new[p] - u[p]) for p in product(range(self.n_equipes), range(self.n_tarefas)))
            if ax < 0.00001:
                break


            print_stats(h, ub, FO_t, lb, gap, mu, norm)
        print(f"Y ótimo: {[y for y in Y_ij if Y_ij[y] > 0]}")
        print(f"X ótimo: {[x for x in var_x_ij if var_x_ij[x] > 0]}")

        h_i = 400
        c = 0
        while h_i > c:
            eq_violada = [eq for eq in range(self.n_equipes) if sum(Y_ij[(eq, j)] * self.consumo_cap[eq][j] for j in range(self.n_tarefas)) > self.capacidades[eq]]
            eq_apta = [eq for eq in range(self.n_equipes) if sum(Y_ij[(eq, j)] * self.consumo_cap[eq][j] for j in range(self.n_tarefas)) <= self.capacidades[eq]]

            # ordenar equipes aptas pela disponibilidade de capacidade!
            eq_apta_s = sorted(eq_apta, key=lambda x: self.capacidades[x] - sum(Y_ij[(x, j)] * self.consumo_cap[x][j] for j in range(self.n_tarefas)), reverse=True)
            eq_violada_s = sorted(eq_violada, key=lambda x: self.capacidades[x] - sum(Y_ij[(x, j)] * self.consumo_cap[x][j] for j in range(self.n_tarefas)))
            for eq_v in eq_violada_s:
                #atividade que consome mais tempo?
                dt = np.array([Y_ij[(eq_v, j)] * self.consumo_cap[eq][j] for j in range(self.n_tarefas)])
                atv = np.argmin(ma.masked_where(dt==0, dt))
                #passo minha atividade para a equipe mais disponível!
                Y_ij[(eq_apta_s[0], atv)] = 1
                #atribuo a atv na equipe com mais capacidade!
                Y_ij[(eq_v, atv)] = 0

            c += 1
            if len(eq_violada) == 0:
                break
        FO_fim = sum([Y_ij[p] * self.custo_alocacao[p[0]][p[1]] for p in product(range(self.n_equipes), range(self.n_tarefas))])
        print_stats(h, ub, FO_fim, lb, gap, mu, norm)
        for tar in range(self.n_tarefas):
            eq = next(eq for eq in range(self.n_equipes) if Y_ij[eq, tar] > 0)
            print(f'Tarefa {tar} - Equipe {eq}')

        print(f'Capacidades disponiveis:: {self.capacidades}')
        a = {eq: sum(Y_ij[(eq, j)] * self.consumo_cap[eq][j] for j in range(self.n_tarefas)) for eq in range(self.n_equipes)}
        print(f'Capacidades utilizadas = {a}')
        print(f'Função Objetivo Calculada por Custo Xij = {FO_fim}')
        return {"valor_fo": FO_fim, "tamanho_instacia": 0, "tempo": 0 }













