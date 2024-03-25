from leitor_dados_txt import leitor_dados
from modelos import *
import time

class Controlador:
    def __init__(self, problema, arquivo=None, num=None):
        self.problema = problema

        if self.problema == "atribuicao_generalizada":
            self.otimiza_atribuicao_generalizada(arquivo=arquivo, num=num)

        elif self.problema == 'TSP':
            self.otimiza_TSP(arquivo=arquivo)

        elif self.problema == 'Alocacao':
            self.otimiza_alocacao(arquivo=arquivo)

        else:
            print(f'O Problema {problema} está incorreto. '
                  f'Tente: atribuicao_generalizada '
                  f'ou TSP'
                  f'ou Alocacao')


    def otimiza_atribuicao_generalizada(self, arquivo=None, num=None):
        dict_results= dict()
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
                resposta = model.otimiza()
                dict_results[f"{arq} - {prs}"] = resposta
                resposta[arquivo] = arq
                resposta[problema] = prs
                dados_consolidados.append(resposta)
        return dict_results

    def otimiza_TSP(self, arquivo=None):
        dict_results = dict()
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
            dict_results[instancia] = resposta
            fim=time.time()
            print('-' * 90)
            print(f'{dados[instancia].name} = {resposta}')
            print(f'Tempo de resposta: {fim - inicio}')
            print('-' * 90)

        return dict_results
    def otimiza_alocacao(self,arquivo=None):
        #TODO:os três metodos dessa classe são quase a mesma coisa. Juntar para 1 só.
        dict_results = dict()
        leitor = leitor_dados(problema=self.problema)
        dados = leitor.le_dados_alocacao_facilities()
        if arquivo:
            dados = {arquivo: dados[arquivo]}
        for instancia in dados:
            print('-' * 90)
            print(f'Rodando o Arquivo {instancia}')
            print('-' * 90)
            inicio = time.time()
            model = AlocacaoFacilities(dados=dados[instancia])
            resposta = model.otimiza()
            dict_results[instancia] = resposta
            fim = time.time()
            print('-' * 90)
            print(f'{instancia} = {resposta}')
            print(f'Tempo de resposta: {fim - inicio}')
            print('-' * 90)

        return dict_results


if __name__ == '__main__':
    #NOME DO PROBLEMA PRECISA SER O MESMO DA PASTA COM AS INSTANCIAS!!
    #

    #Passar nome do arquivo para rodar exclusivamente algum. Precisa ser do mesmo problema!
    #arquivo = None

    #num: Cada arquivo de cada problema tem um número de problemas. usar para passar algum especifico!
    problema = "Alocacao"
    controle_aloc = Controlador(problema=problema).otimiza_alocacao()

    problema = "atribuicao_generalizada"
    controle_atr = Controlador(problema=problema).otimiza_atribuicao_generalizada()

    problema = "TSP"
    controle_tsp = Controlador(problema=problema).otimiza_TSP()

    b=0


