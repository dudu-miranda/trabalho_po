
from copy import *

class Matriz(object):
	"""Classe que implementará as operações básicas de matrizes."""
	def __init__(self, listas):

		self.m = listas
		self.lin = len(listas)
		self.col = len(listas[0])


	#Função que retorna a forma escrita de uma matriz
	def __str__(self):
		string = ''
		for i in self.m:
			string += '| '
			for j in i:
				string += str(j) + ' '
			string += '|\n'
		return string
		

	#Cria matriz cheia de zeros
	def criaMatriz0(self,linha,coluna):
		matriz = []
		for i in range(linha):
			matriz.append([])
			for j in range(coluna):
				matriz[i].append(0)

		return Matriz(matriz)


	#Função que cria um vetor a partir de uma coluna em especifico da matriz
	def getColuna(self,ind):
		col = []
		for i in range(self.lin):
			col.append(self.m[i][ind])

		return Matriz(col)

	#Função que cria uma matriz identidade de tamanho N
	def criaIdentidade(self, n):

		m = []
		for i in range(0,n):
			m.append([])
			for j in range(0,n):
				if(i==j):
					m[i].append(1)
				else:
					m[i].append(0)

		return m


	#Função que multiplica uma matriz por um escalar OK
	def multEscalar(self, escalar):

		temp = []
		for i in range(0, len(self.m)):
			temp.append([])

			for j in range(0, len(self.m[i])):
				temp[i].append(self.m[i][j] * escalar)

		return Matriz(temp)


	#Função que transpôe uma matriz OK
	def transpostaMatriz(self):

		temp = []

		for i in range(0, len(self.m[0])):
			temp.append([])

			for j in range(0, len(self.m)):
				temp[i].append(self.m[j][i])

		return Matriz(temp)


	#Função que multiplica duas matrizes OK
	def multiplicaMatriz(self,matrizB):

		if(self.col!=matrizB.lin):
			print('Impossível multiplicar essas matrizes!')
			exit()

		temp = []
		for i in range(0, self.lin):
			temp.append([])

			for j in range(0, matrizB.col):
				temp[i].append(0)

				for k in range(0, self.col):
					temp[i][j] += self.m[i][k] * matrizB.m[k][j]

		return Matriz(temp)


	#Função que retorna o produto escalar de dois vetores OK
	def produtoEscalar(self,vetor):

		return self.multiplicaMatriz(vetor)[0][0]


	#Função que será utilizada para calculo da matriz inversa parte da decomposição LU (Acho que ta ok)
	def decomposicaoLU(self,n,A):

		pivot = []

		#Inicialização ordenada de pivot[]
		for i in range(0,n):
			pivot.append(i)

		for j in range(0,n-2):
		
			#Escolha do pivô
			p, maior = j, abs(A[0][j])
			for k in range(j+1,n):
				if(abs(A[k][j]) > maior):
					maior = abs(A[k][j])
					p = k

			if(p != j):
				A[j], A[p] = A[p], A[j]
				#for k in range(0, n):
					#Troca das linhas p e j (PS: Ver se da pra trocar uma linha de uma vez dps)
				#	A[j][k], A[p][k] = A[p][k], A[j][k]

				#Armazena as permutas de B
				pivot[j], pivot[p] = pivot[p], pivot[j]

			if(abs(A[j][j]) != 0):
				for i in range(j+1, n):
					#Pivoteamento por eliminação de Gauss
					multiplicador = A[i][j] / A[j][j]
					#Multiplicadores m[i][j] de L
					A[i][j] = A[i][j] * multiplicador

					for k in range(j+1, n):
						#Coeficientes U[i][j] de U
						A[i][k] = A[i][k] - multiplicador * A[j][k]

		#print(pivot)
		#print(str(A[0]))
		#print(str(A[1]))
		#print(str(A[2]))

		#Retorna matriz A = L/U e vetor Pivot[]
		return A, pivot


	#1 parte do algoritimo que fará a inversa da matriz
	def substuicaoSucessivaPivotal(self, n, L, b, pivot):
		#Solução do sistema Triangular inferior
		k = pivot[0]
		y = [0]*(n)
		y[0] = b[k]

		for i in range(1,n):
			soma = 0

			for j in range(0, i-1):
				soma += L[i][j] * y[j]

			#Obtém indice de b para o i-ésimo pivô
			k = pivot[i]
			y[i] = b[k] - soma

		return y


	#Segunda parte do algoritimo que fará a inversa da matriz
	def substituicaoRetroativa(self, n, U, y):

		x = [0]*(n)
		x[n-1] = y[n-1] / U[n-1][n-1]

		for i in range(n-2, 0, -1):
			soma = 0

			for j in range(i+1, n):
				soma += U[i][j] * x[j]

			x[i] = (y[i] - soma) / U[i][i]
		
		return x


	#Função que usa a decomposição LU para retornar a matriz inversa
	def matrizInversa(self):

		#LU, pivo = self.decomposicaoLU(self.lin,self.m)
		identidade = self.criaIdentidade(self.lin)
		inversa = []

		for i in range(0,self.lin):
			#y = self.substuicaoSucessivaPivotal(self.lin, LU, identidade[i], pivo)
			#inversa.append( self.substituicaoRetroativa(self.lin, LU, y) )
			X = deepcopy(self.m)
			
			for j in range(0, len(X)):
				X[j].append(identidade[i][j])

			pedaco = self.gaussJordan(X)
			col = []
			for el in X:
				col.append(el[-1])

			inversa.append(col)

		return Matriz(inversa).transpostaMatriz()


	#Codigo do gauss jordan de pivoteamento para matriz inversa
	def gaussJordan(self, m):
		"""
		Retorna True caso sucesso e False caso 'm' seja singular.
		Por parâmetros tem m que é a lista de listas que contém a matriz
		"""
		
		singular = 1.0/(10**10)

		(lin, col) = (len(m), len(m[0]))
		for i in range(0,lin):
			limiteLinha = i
			for j in range(i+1, lin):    # Acha o maior pivô
				if abs(m[j][i]) > abs(m[limiteLinha][i]):
					limiteLinha = j
			(m[i], m[limiteLinha]) = (m[limiteLinha], m[i])
			if abs(m[i][i]) <= singular:     # Singular?
				return False
			
			for j in range(i+1, lin):    # Elimina coluna i
				c = m[j][i] / m[i][i]
				for x in range(i, col):
					m[j][x] -= m[i][x] * c
		
		for i in range(lin-1, 0-1, -1): # Substituição retroativa
			c  = m[i][i]
			for j in range(0,i):
				for x in range(col-1, i-1, -1):
					m[j][x] -=  m[i][x] * m[j][i] / c
			m[i][i] /= c
			for x in range(lin, col):       # Normaliza linha i
				m[i][x] /= c

		return True

		
	def gaussJordan2(self, A, m, n):

		for i in range(0, m-1):
			for j in range(i+1, m):
				multiplicador = A[j][i] / A[i][i]

				for k in range(0, n):
					A[j][k] = A[j][k] - multiplicador * A[i][k]

		for i in range(m-1,1,-1):
			for j in range(i-1, 0, -1):
				
				multiplicador = A[j][i] / A[i][i]
				for k  in range(0,n):
					A[j][k] = A[j][k] - multiplicador * A[i][k]

		for i in range(0,m):
			for j in range(0,n):
				A[i][j] = A[i][j] / A[i][i]

		return A




























class Simplex(object):
	"""docstring for Simplex"""
	#Model criaModelo(Matrix A, Matrix b, Matrix c, Matrix base, Matrix naoBase, Matrix artificiais)
	def __init__(self, arquivo):	

		self.A, self.b, self.custo, self.base, self.naoBase, self.artificiais = self.carregaModelo(arquivo)

		#print(self.A)
		#print(self.b)
		#print(self.custo)
		#print(self.base)
		#print(self.naoBase)
		#print(self.artificiais)

		
		self.A          	 # Matriz A.
		self.b = self.b.transpostaMatriz()          	 # Vetor de igualdade das restricoes.
		self.custo = self.custo.transpostaMatriz()    	 # Vetor de custo (min).
		self.base = self.base.transpostaMatriz()   	     # Vetor com as variaveis na base.
		self.naoBase 	     # Vetor com as variaveis nao-base.
		self.artificiais	 # Vetor com as variaveis artificiais.
		self.B = []   	     # Matriz das colunas de cada variavel base.
		self.solucao = None	 # Vetor solucao das variaveis presentes na base.
		
		self.B = self.A.criaMatriz0(self.A.lin,self.A.lin)

		self.simplex()


	# Sub-rotinas:
	def carregaModelo(self, arquivo):
		# Abre o arquivo para leitura:
		parq = open(arquivo, "r");

		# Caso o arquivo nao exista, retorna NULL:
		if(parq == None):
			print("Arquivo nao encontrado.\n")
			return None

		# Determinar quantidade de linhas e colunas;
		temp = parq.readline()
		i,j,a = temp.split()
		i, j, a = int(i), int(j), int(a)

		# Cria vetor b:
		b = Matriz([[0]*i]*1)
		temp = parq.readline().split()
		for k in range(0,len(temp)):
			b.m[0][k] = float(temp[k])

		# Cria vetor de custo:
		c = Matriz([[0]*j]*1)
		temp = parq.readline().split()
		for k in range(len(temp)):
			c.m[0][k] = float(temp[k])

		# Cria vetor base:
		base = Matriz([[0]*i]*1)
		temp = parq.readline().split()
		for k in range(i):
			base.m[0][k] = int(temp[k])- 1

		# Cria vetor nao-base:
		naoBase = Matriz([[0]*(j - i)]*1)
		temp = parq.readline().split()
		for k in range(j - i):
			naoBase.m[0][k] = int(temp[k]) - 1	

		# Cria vetor de variaveis artificiais:
		artificiais = Matriz([[0]*a]*1)
		if(a > 0):
			temp = parq.readline().split()
			for k in range(a):
				artificiais.m[0][k] = valor - 1	

		# Cria matriz A:
		A = []
		# Le e salva cada elemento na matriz:
		for i in range(0,i):
			A.append(parq.readline().split())

		#Converte esses elementos para inteiro
		for linha in range(len(A)):
			for el in range(len(A[0])):
				A[linha][el] = int(A[linha][el])


		parq.close()
		A = Matriz(A)

		return A, b, c, base, naoBase, artificiais




	def simplex(self):
		''' 
		 * Executa o simplex e atualiza o campo 'solucao' do modelo.
		 * Retorna -1 se o modelo nao possuir solucao;
		 * Retorna -2 se a solucao for ilimitada;
		 * Retorna -3 se houverem multiplas solucoes;
		 * Retorna >= 0 se houver uma unica solucao otima. Retorna numero de iteracoes.
		 * [DEBUG] matImprime(Matrix m) mostra qualquer matriz/vetor no console.
		 * [DEBUG] imprimeModelo(Model m) mostra todos os dados do simplex na iteracao atual.
		'''

		m = self.A.lin
		n = self.A.col
		iteracoes = 0

		while True:
			# Calcula a B^-1:
			self.geraMatrizBase();
			invB = self.B.matrizInversa()

			# (REQUISITO 03) Calcular a SBF inicial (invB x b):
			self.solucao = invB.multiplicaMatriz(self.b)
			# Calcula o valor da funcao objetivo:
			objetivo = self.calcObjetivo();

			# Monta o vetor de custo basico:
			custoBase = self.A.criaMatriz0(self.base.lin,1)
			
			for i in range(m):
				indice = self.base.m[i][0]
				valor = self.custo.m[indice][0]
				custoBase.m[i][0] = valor

			# Determina qual variavel entra na base:
			custoBaseT = custoBase.transpostaMatriz()
			custoReduzido = self.A.criaMatriz0(1, n)
			menorCusto = 99999999
			indiceMenor = -1
			#Breakpoint

			for i in range(self.naoBase.lin):

				indice = self.naoBase.m[i][0]

				# Custo do indice na funcao objetivo:
				custo = self.custo.m[indice][0]

				# Coluna do indice nao basico da matriz A:
				Aj = self.A.getColuna(indice)

				# (REQUISITO 05) Direcao:
				direcao = invB.multiplicaMatriz(Aj)

				# resOp1 = custoBase x invB:
				resOp1 = custoBaseT.multiplicaMatriz(invB)

				# resOp2 = custoBase x invB x Aj:
				resOp2 = resOp1,multiplicaMatriz(Aj)
				
				# (REQUISITO 04) Custo reduzido:
				custo = custo - resOp2.m[0][0]
				custoReduzido.m[0][indice] = custo
				
				# Trata valores de custo muito pequenos e negativos (i.e: -1e-16).
				# (REQUISITO 09) Degeneracao - Garantir que custos proximos de 0 sejam +0.
				# Arredonda o valor para 5 casas decimais.
				fac = pow(10, 5)
				custo = round(custo * fac) / fac

				# (REQUISITO 09) Degeneracao - Custos reduzidos iguais a zero nao entram na base.
				# Se o custo for negativo:
				if(custo < 0):
					if(custo < menorCusto):
						# Atualiza variavel candidata para entrar na base:
						indiceMenor = indice
						menorCusto = custo
					
				
				
				# [DEBUG] (REQUISITO 05) Mostra direcao encontrada:
				# printf("Direcao factivel %d, custo reduzido %lf\n", indice, custo);
				# matImprime(direcao);

				# Libera matrizes utilizadas nas operacoes:
				Aj = None
				resOp1 = None
				resOp2 = None
				direcao = None
			

			# [DEBUG] Mostra o custo reduzido encontrado:
			# printf("Custo reduzido:\n");
			# matImprime(custoReduzido);

			# Nao houve nenhum indice com custo reduzido negativo.
			# Solucao otima encontrada:
			if indiceMenor == -1:
				objetivo = self.calcObjetivo()

				# (REQUISITO 09) - Solucao inexistente (Variaveis artificiais presentes na solucao):
				if self.artificialNaBase():
					printf("[Simplex] Solucao inexistente.\n")
					return -1
				
				# (REQUISITO 09) - Multiplos Otimos:
				elif self.existeNaoBase0(custoReduzido):
					printf("[Simplex] Multiplos Otimos.\nIteracoes = %d\nCusto otimo = %lf\n", iteracoes, objetivo);
					return -3;
				
				# (REQUISITO 09) - Solucao Otima Unica:
				else:
					printf("[Simplex] Solucao Otima Unica.\nIteracoes = %d\nCusto otimo = %lf\n", iteracoes, objetivo);
					return iteracoes;
				
			

			custoBaseT = None
			custoReduzido = None

			# [DEBUG] Indice da variavel a entrar na base:
			# printf("Entra na base %d, custo = %lf\n", indiceMenor, menorCusto);

			# (REQUISITO 09) - Solucao Ilimitada:
			Aj = self.A.getColuna(indiceMenor)
			u = invB.multiplicaMatriz(Aj)
			if not self.existePositivo(u):
				printf("[Simplex] Solucao Ilimitada.\nCusto otimo = -Infinito\n");
				return -2;
			
			Aj = None
			u = None

			# (REQUISITO 06) Determina o valor de theta:
			theta = 99999999
			indice = -1;

			for i in range(m):
			
				if u.m[i][0] > 0:

					razao = self.solucao.m[i][0] / u.m[i][0]

					if(razao < theta):
						theta = razao
						indice = self.base.m[i][0]

			# [DEBUG] Indice da variavel a sair da base:
			# printf("Variavel sai da base: %d, theta = %lf\n", indice, theta);

			# (REQUISITO 07) Atualiza a base:
			for i in range(m):

				if self.base.m[i][0] == indice:
					self.solucao.m[i][0] = theta
					self.base.m[i][0] = indiceMenor
				

			# (REQUISITO 07) Atualiza o conjunto de variaveis nao-basicas:
			for i in range(n - m):
			
				if self.naoBase.m[i][0] == indiceMenor:
					self.naoBase.m[i][0] = indice
				
			# Libera a matriz inversa utilizada para a iteracao atual:
			invB = None

			iteracoes+= 1

		return -1;


	#Calcula a função objetivo do modelo
	def calcObjetivo(self):

		objt = 0
		for i in range(self.base.lin):
			base = self.base.m[i][0]
			objt += self.custo.m[base][0] * self.solucao.m[i][0]

		return objt


	#Função para imprimir o modelo no seu estado atual
	def imprimeModelo(self):
		print("\nMatriz A:\n");
		print(self.A);

		print("\nVetor b:\n");
		print(self.b);

		print("\nVetor custo:\n");
		print(self.custo);

		print("\nVetor base (indices contam a partir de 0):\n");
		print(self.base);

		print("\nVetor nao-base (indices contam a partir de 0):\n");
		print(self.naoBase);

		print("\nMatriz basica:\n");
		print(self.B);


	# Checa se existe um valor positivo na linha 'x':
	def existePositivo(self,x):

		for i in range(x.lin):
			if x[i][0] > 0:
				return True

		return False


	#Checa se existe uma variavel artificial na base:
	def artificialNaBase(self):

		for i in range(self.artificiais.lin):
			var = self.artificiais.m[i][0]
			for j in range(self.base.lin):
				if(var == self.base[j][0]):
					return 1

		return 0


	# Copia as colunas da base 'modelo -> base' para a matriz B:
	def geraMatrizBase(self):

		valor = 0
		for i in range(self.A.lin):
			base = self.base.m[i][0]

			for j in range(self.A.lin):
				valor = self.A.m[j][base]
				self.B.m[j][i] = valor


	#Checa se existe uma variavel nao base com custo reduzido zero:
	def existeNaoBase(self, custoReduzido):

		if custoReduzido == None:
			return 0

		for i in range(custoReduzido.col):
			custo = custoReduzido.m[0][i]
		
			# Trata valores muito pequenos de custo (i.e: 1e-16).
			# Arredonda o valor para 5 casas decimais:
			fac = pow(10, 5);
			custo = round(custo * fac) / fac;

			# Encontrou um custo reduzido igual a zero:
			if(custo == 0):
				# Olha se o indice 'i' existe na base:
				existe = 0
				for j in range(self.base.lin):
					if self.base.m[j][0] == i:
						existe = 1
						break 
				
				if(not existe):
					return 1;

		return 0


	#Saida do programa
	def outputModelo(self, iteracoes, entrada, saida):
		if((self == None) or (self.solucao == None)):
			return

		arq = open(saida, "w")

		if iteracoes == -1:
			arq.write("Modelo: "+str(entrada)+"\nSolucao inexistente\n")
			
		elif iteracoes == -2:
			arq.write("Modelo: "+str(entrada)+"\nSolucao ilimitada\n")
		
		elif iteracoes == -3:
			arq.write("Modelo: "+str(entrada)+"\nMultiplos otimos\n")

		else:
			arq.write("Modelo: "+str(entrada)+"\nSolucao unica otima encontrada\n\nIteracoes: "+str(iteracoes))

		arq.write("\nOtimo: "+str(self.calcObjetivo)+"\n\n")

		# Exibe o valor das variaveis:
		for i in range(self.custo.lin):
			valor = 0;

			# Procura se a variavel 'i' esta na base:
			for j in range(self.base.lin):
				base = self.base.m[j][0]
				if(i == base):
					# Pega o valor da base na solucao:
					valor = self.solucao.m[j][0]
				
			arq.write("x[" + str(i+1) + "] = " + str(valor) + "\n")		

		arq.close();
	

















"""
1 2 3
4 5 6
7 8 9

1 4 7
2 5 8 
3 6 9
"""

def menuzinho():

	s = Simplex()

	op=1
	while op!=0:
		print('Digite 0 para sair do programa')
		print('Digite 1 para carregar o arquivo do modelo.\n')
		print('Digite 2 para modificar o modelo.\n')
		op = int(input())
		if(op==1):
			matriz = carregarArquivo()
			s.calcula()
			s.mostra()
		elif(op==2):
			matriz = modificaModelo()
			s.calcula()
			s.mostra()

#menuzinho()


s = Simplex("goldbard-pg104.mod")