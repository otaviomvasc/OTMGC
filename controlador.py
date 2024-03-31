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
        else:
            print(f'O Problema {problema} está incorreto. '
                  f'Tente: atribuicao_generalizada '
                  f'ou TSP'
                  f'ou Alocacao')


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
                inicio = time.time()
                resposta = model.otimiza()
                fim = time.time() - inicio
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
            distancias_d = dict()
            distancias_dict_valores = dict()
            for cliente in range(d['clientes']):
                dict_aux = {cb: d['matriz_custo'][cb] for cb in d['matriz_custo'] if cb[1] == cliente}
                distancias_d[cliente] = [i[0] for i in sorted(dict_aux, key=dict_aux.get, reverse=False)]
                distancias_dict_valores[cliente] = {i: dict_aux[i] for i in
                                                    sorted(dict_aux, key=dict_aux.get, reverse=False)}

            return distancias_d, distancias_dict_valores
        #TODO:os três metodos dessa classe são quase a mesma coisa. Juntar para 1 só.
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
            #inicio = time.time()
            D_Ordenado, _ = cria_matriz_custo_ordenada(d=dados[instancia])
            dados[instancia]["D_ord"] = D_Ordenado
            model = Benders(dat=dados[instancia])
            resposta = model.run()
            resposta['problema'] = 'alocacao'
            resposta['instancia'] = instancia
            #dict_results[instancia] = resposta
            #fim = time.time()
            # print('-' * 90)
            # print(f'{instancia} = {resposta}')
            # print(f'Tempo de resposta: {fim - inicio}')
            # print('-' * 90)
            # resposta['tempo'] = fim - inicio
            dados_result.append(copy.deepcopy(resposta))

        return dados_result


if __name__ == '__main__':
    #NOME DO PROBLEMA PRECISA SER O MESMO DA PASTA COM AS INSTANCIAS!!
    #

    #Passar nome do arquivo para rodar exclusivamente algum. Precisa ser do mesmo problema!
    #arquivo = None
    exporta = True
    #num: Cada arquivo de cada problema tem um número de problemas. usar para passar algum especifico!

    problema = "AlocacaoBENDERS"
    controle_aloc = Controlador(problema=problema).otimiza_Benders()

    # problema = "AlocacaoPrimal"
    # controle_aloc = Controlador(problema=problema).otimiza_alocacao()
    #
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

