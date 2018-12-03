
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
		for i in m:
			for j in i:
				string += str(j)
			string += '\n'
		return string
		

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

		return temp


	#Função que transpôe uma matriz OK
	def transpostaMatriz(self):

		temp = []

		for i in range(0, len(self.m[0])):
			temp.append([])

			for j in range(0, len(self.m)):
				temp[i].append(self.m[j][i])

		return temp


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

		return temp


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
	def __init__(self, A, b, c, base, naoBase, artificiais):	

		self.A = A          			 # Matriz A.
		self.b = b            			 # Vetor de igualdade das restricoes.
		self.custo = custo        		 # Vetor de custo (min).
		self.base = base   			     # Vetor com as variaveis na base.
		self.naoBase = naoBase  	     # Vetor com as variaveis nao-base.
		self.artificiais = artificiais	 # Vetor com as variaveis artificiais.
		self.B = []   			         # Matriz das colunas de cada variavel base.
		self.solucao = None    			 # Vetor solucao das variaveis presentes na base.
		
		self.B = [[0]*self.A.lin]*self.A.lin


	#Calcula a função objetivo do modelo
	def calcObjetivo(self):

		objt = 0

		for i in range(self.base.lin):
			base = self.base[i][0]
			objt += self.custo[base][i] * self.solucao[i][0]

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
				return 1

		return 0

	#Checa se existe uma variavel artificial na base:
	def artificialNaBase(self):

		for i in range(self.artificiais.lin):
			var = self.artificiais.m[i][0]
			for j in range(self.base.lin):
				if(var == self.base[j][0]):
					return 1

		return 0

	def geraMatrizBase(self) {
	# Copia as colunas da base 'modelo -> base' para a matriz B:

	valor = 0

	for i in range(self.A.lin):
		base = self.base.m[i][0]

		for j in range(self.A.lin):
			valor = self.A.m[j][base]
			self.B.m[j][i] = valor































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