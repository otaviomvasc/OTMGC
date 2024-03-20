from leitor_dados_txt import leitor_dados
from modelos import *

class Controlador:
    def __init__(self, problema, arquivo=None, num=None):
        self.problema = problema

        if self.problema == "atribuicao_generalizada":
            self.otimiza_atribuicao_generalizada(arquivo=arquivo, num=num)

        if self.problema == 'TSP':
            self.otimiza_TSP(arquivo=arquivo)

    def otimiza_atribuicao_generalizada(self, arquivo=None, num=None):
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
                resposta[arquivo] = arq
                resposta[problema] = prs
                dados_consolidados.append(resposta)
        return dados_consolidados

    def otimiza_TSP(self, arquivo=None):
        leitor = leitor_dados(problema=self.problema)
        dados = leitor.le_dados_TSP()
        if arquivo:
            dados = dados[arquivo]
        for instancia in dados:
            model = ModelTSP(dados=dados[instancia])
            resposta = model.otimiza()
            print(f'{dados[instancia].name} = {resposta}')


if __name__ == '__main__':
    #NOME DO PROBLEMA PRECISA SER O MESMO DA PASTA COM AS INSTANCIAS!!
    #problema = "atribuicao_generalizada"

    #Passar nome do arquivo para rodar exclusivamente algum. Precisa ser do mesmo problema!
    #arquivo = None

    #num: Cada arquivo de cada problema tem um número de problemas. usar para passar algum especifico!
    #controle = Controlador(problema=problema)

    problema = "TSP"
    controle = Controlador(problema=problema)
