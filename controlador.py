from leitor_dados_txt import leitor_dados
from modelos import *
import time
import pandas as pd
from Decomposicao_Benders import *

class Controlador:
    def __init__(self, problema, arquivo=None, num=None):
        self.problema = problema

        if self.problema == "atribuicao_generalizada":
            self.otimiza_atribuicao_generalizada(arquivo=arquivo, num=num)

        elif self.problema == 'TSP':
            self.otimiza_TSP(arquivo=arquivo)

        elif self.problema == 'Alocacao':
            self.otimiza_alocacao(arquivo=arquivo)

        elif self.problema == 'AlocacaoBENDERS':
            self.otimiza_Benders(arquivo=arquivo)

        elif self.problema == "AlocacaoBENDERSBendersBranchAndCut":
            self.otimiza_Benders_BranchAndCut(arquivo=arquivo)

        elif self.problema == "GapLagrange":
            self.otimiza_Lagrange(arquivo=arquivo)

        else:
            print(f'O Problema {problema} está incorreto. '
                  f'Tente: atribuicao_generalizada '
                  f'ou TSP'
                  f'ou Alocacao')

    def otimiza_Lagrange(self, arquivo=None, num=None):
        dados_result = list()
        leitor = leitor_dados(problema="atribuicao_generalizada")
        dados = leitor.le_dados_atribuicao_gen()
        # TODO: testar rodada do modelo!!
        if arquivo:
            dados = dados[arquivo]
            if num:
                dados = dados[int(num)]
        dados_consolidados = list()
        for arq in dados:
            for prs in dados[arq]:
                model = EgapLangrange(dados=dados[arq][prs])
                #inicio = time.time()
                resposta = model.otimiza()
                #fim = time.time() - inicio
                resposta['problema'] = 'atr_generalizada'
                resposta['instancia'] = f'{arq} - {prs}'
                # dict_results[f"{arq} - {prs}"] = resposta
                resposta['problema'] = 'atr_generalizada'
                resposta['instancia'] = f'{arq} - {prs}'
                resposta['tempo'] = 0
                dados_result.append(copy.deepcopy(resposta))
        return dados_result


    def otimiza_atribuicao_generalizada(self, arquivo=None, num=None):
        dados_result = list()
        leitor = leitor_dados(problema=self.problema)
        dados = leitor.le_dados_atribuicao_gen()
        #TODO: testar rodada do modelo!!
        if arquivo:
            dados = dados[arquivo]
            if num:
                dados = dados[int(num)]
        dados_consolidados = list()
        for arq in dados:
            for prs in dados[arq]:
                model = ModelAtribuicaoGeneralizada(dados=dados[arq][prs])
                #inicio = time.time()
                resposta = model.otimiza()
                fim = 0 - 0
                resposta['problema'] = 'atr_generalizada'
                resposta['instancia'] = f'{arq} - {prs}'
                #dict_results[f"{arq} - {prs}"] = resposta
                resposta['problema'] = 'atr_generalizada'
                resposta['instancia'] = f'{arq} - {prs}'
                resposta['tempo'] = fim
                dados_result.append(copy.deepcopy(resposta))

        return dados_result

    def otimiza_TSP(self, arquivo=None):
        dados_result = list()
        leitor = leitor_dados(problema=self.problema)
        dados = leitor.le_dados_TSP()
        if arquivo:
            dados = {arquivo: dados[arquivo]}
        for instancia in dados:
            print('-' * 90)
            print(f'Rodando o Arquivo {instancia}')
            print('-' * 90)
            inicio = time.time()
            model = ModelTSP(dados=dados[instancia])
            resposta = model.otimiza(subrota_fluxo_ficticio=True, subrota_TMZ=False, sub_rota_iterativa=False)
            resposta['problema'] = 'TSP'
            resposta['instancia'] = instancia
            fim=time.time()
            print('-' * 90)
            print(f'{dados[instancia].name} = {resposta}')
            print(f'Tempo de resposta: {fim - inicio}')
            print('-' * 90)
            resposta['tempo'] = fim - inicio
            dados_result.append(copy.deepcopy(resposta))

        return dados_result

    def otimiza_Benders(self, arquivo=None):
        def cria_matriz_custo_ordenada(d):
            #TODO: SERÁ QUE COMPENSA TRABALHAR COM 2 DICTS AGORA QUE ELES VÃO SER CONSULTADOS E ITERADOS MAIS VEZES PELO BENDERS?
            #TODO: preciso retirar as instâncias com distâncias repetidas!!
            distancias_d = dict()
            distancias_dict_valores = dict()
            for cliente in range(d['clientes']):
                dict_aux = {cb: d['matriz_custo'][cb] for cb in d['matriz_custo'] if cb[1] == cliente}
                distancias_d[cliente] = [i[0] for i in sorted(dict_aux, key=dict_aux.get, reverse=False)]
                distancias_dict_valores[cliente] = {i: dict_aux[i] for i in
                                                    sorted(dict_aux, key=dict_aux.get, reverse=False)}

            return distancias_d, distancias_dict_valores
        #TODO:os quatro metodos dessa classe são quase a mesma coisa. Juntar para 1 só.
        dados_result = list()
        self.problema = "Alocacao"
        leitor = leitor_dados(problema=self.problema)
        dados = leitor.le_dados_alocacao_facilities()
        if arquivo:
            dados = {arquivo: dados[arquivo]}
        for instancia in dados:
            print('-' * 90)
            print(f'Rodando o Arquivo {instancia}')
            print('-' * 90)
            inicio = time()
            D_Ordenado, _ = cria_matriz_custo_ordenada(d=dados[instancia])
            dados[instancia]["D_ord"] = D_Ordenado
            model = Benders(dat=dados[instancia])
            resposta, fo = model.run()
            tempo_total = time() - inicio
            dados_result.append({"Cenario": instancia, "tempo": tempo_total, "FO": fo})

        return dados_result

    def otimiza_Benders_BranchAndCut(self, arquivo=None):
        def cria_matriz_custo_ordenada(d):
            # TODO: SERÁ QUE COMPENSA TRABALHAR COM 2 DICTS AGORA QUE ELES VÃO SER CONSULTADOS E ITERADOS MAIS VEZES PELO BENDERS?
            distancias_d = dict()
            distancias_dict_valores = dict()
            for cliente in range(d['clientes']):
                dict_aux = {cb: d['matriz_custo'][cb] for cb in d['matriz_custo'] if cb[1] == cliente}
                distancias_d[cliente] = [i[0] for i in sorted(dict_aux, key=dict_aux.get, reverse=False)]
                distancias_dict_valores[cliente] = {i: dict_aux[i] for i in
                                                    sorted(dict_aux, key=dict_aux.get, reverse=False)}

            return distancias_d, distancias_dict_valores

        # TODO:os quatro metodos dessa classe são quase a mesma coisa. Juntar para 1 só.
        dados_result = list()
        self.problema = "Alocacao"
        leitor = leitor_dados(problema=self.problema)
        dados = leitor.le_dados_alocacao_facilities()
        if arquivo:
            dados = {arquivo: dados[arquivo]}
        for instancia in dados:
            print('-' * 90)
            print(f'Rodando o Arquivo {instancia}')
            print('-' * 90)
            inicio = time()
            D_Ordenado, _ = cria_matriz_custo_ordenada(d=dados[instancia])
            dados[instancia]["D_ord"] = D_Ordenado
            model = BendersBranchAndCut(dat=dados[instancia])
            resposta, fo = model.run()
            tempo_total = time() - inicio
            dados_result.append({"Cenario": instancia, "tempo": tempo_total, "FO": fo})

        return dados_result

"""
a = [11, 11, 0, 24, 11, 0, 16, 5, 11, 11, 3, 10, 5, 0, 6, 11, 3, 10, 3, 6, 3, 6, 10, 0, 11, 10, 12, 10, 10, 0, 0, 10, 0, 16, 11, 11, 24, 5, 23, 24, 10, 3, 11, 6, 22, 23, 23, 6, 24, 11]
array([ 0,  3,  5,  6, 10, 11, 12, 16, 22, 23, 24])
"""
if __name__ == '__main__':
    #NOME DO PROBLEMA PRECISA SER O MESMO DA PASTA COM AS INSTANCIAS!!
    #

    #Passar nome do arquivo para rodar exclusivamente algum. Precisa ser do mesmo problema!
    #arquivo = None
    exporta = True
    #num: Cada arquivo de cada problema tem um número de problemas. usar para passar algum especifico!

    # problema = "AlocacaoBENDERS"
    # controle_aloc = Controlador(problema=problema).otimiza_Benders()
    arquivo = "gap_1.txt"
    problema = "atribuicao_generalizada"
    controle_atr = Controlador(problema=problema).otimiza_atribuicao_generalizada()

    problema  = "GapLagrange"
    lagrange = Controlador(problema=problema)
    for i in range(len(controle_atr)):
        print(f'Exato - Heuristica: {controle_atr[i]["valor_fo"] - lagrange[i]["valor_fo"]}')
    b=0
    # problema = "AlocacaoBENDERSBendersBranchAndCut"
    # controle_aloc_2 = Controlador(problema=problema).otimiza_Benders_BranchAndCut()

    # df_benders = pd.DataFrame(controle_aloc)
    # df_benders['Método'] = 'Benders_MultiCut'
    #
    # df_benders_Branch_and_cut = pd.DataFrame(controle_aloc_2)
    # df_benders_Branch_and_cut['Método'] = 'Benders_BranchAnCut'
    #
    # #df_fim = pd.concat([df_benders, df_benders_Branch_and_cut])
    # df_fim = df_benders.merge(df_benders_Branch_and_cut, on='Cenario')
    # df_fim['Diff_tempo'] = df_fim.tempo_x - df_fim.tempo_y
    # df_fim['Diff_FO'] = df_fim.FO_x - df_fim.FO_y
    # df_fim.to_excel("Respostas_tp2.xlsx")

    # problema = "atribuicao_generalizada"
    # controle_atr = Controlador(problema=problema).otimiza_atribuicao_generalizada()
    #
    # problema = "Alocacao"
    # controle_aloc = Controlador(problema=problema).otimiza_alocacao()
    #
    #
    # problema = "TSP"
    # controle_tsp = Controlador(problema=problema).otimiza_TSP()

    # if exporta:
    #     list_fim = list()
    #     list_fim.extend(controle_aloc)
    #     list_fim.extend(controle_atr)
    #     list_fim.extend(controle_tsp)
    #     df = pd.DataFrame(list_fim)
    #     df.to_excel("output.xlsx")

