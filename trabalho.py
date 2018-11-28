

class Matriz(object):
	"""Classe que implementará as operações básicas de matrizes."""
	def __init__(self, listas):

		self.m = listas
		self.lin = len(listas)
		self.col = len(listas[0])
		

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


	def produtoEscalar(self,vetor):

		return self.multiplicaMatriz(vetor)[0][0]
		
"""
1 2 3
4 5 6
7 8 9

1 4 7
2 5 8 
3 6 9
"""

m1 = [ [1], [1], [2], [3] ]
ident = [ [1,2,3,4] ]
ident1 = Matriz(ident)

mat1 = Matriz(m1)
print(str(ident1.produtoEscalar(mat1)))