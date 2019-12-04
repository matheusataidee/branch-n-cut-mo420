#coding: utf-8

import sys
import glob

if '__main__' == __name__:

	for url_file in sorted(glob.glob('./solutions/*')):

		temp = url_file.split('.')

		if(temp[3] == 'sol'):

			rfile = open(url_file, 'r')

			lines = rfile.readlines()
			for i in list(range(len(lines))): lines[i] = lines[i].split('\n')
			for i in list(range(len(lines))): lines[i] = lines[i][0]

			lines[0] = lines[0].split('_')

			n = int(lines[0][2])

			vertices_grafo = []
			arestas_grafo = []

			for i in list(range(1, n)):
				temp = lines[i].split(' ')
				arestas_grafo.append((int(temp[0]), int(temp[1])))

			if n - 1 == len(arestas_grafo):

				vertices_grafo.append(arestas_grafo[0][0])
				vertices_grafo.append(arestas_grafo[0][1])
				arestas_grafo.remove((arestas_grafo[0][0], arestas_grafo[0][1]))

				find = True
				while find and len(arestas_grafo) > 0:

					find = False
					arestas_remove = []
					vertices_add = []

					for aresta in arestas_grafo:

						origem = aresta[0]
						destino = aresta[1]
						entry = 0

						if origem in vertices_grafo:
							vertices_add.append(destino)
							arestas_remove.append((origem, destino))
							find = True
							entry = entry + 1

						if destino in vertices_grafo:
							vertices_add.append(origem)
							arestas_remove.append((origem, destino))
							find = True
							entry = entry + 1

						if(entry == 2): arestas_remove.remove((origem, destino))

					for vertice in vertices_add: vertices_grafo.append(vertice)
					for aresta in arestas_remove: arestas_grafo.remove(aresta)

				vertices_grafo = sorted(set(vertices_grafo))

				if len(vertices_grafo) == n and len(arestas_grafo) == 0: print('{}: o grafo esta conectado.'.format(url_file))
				else: print('{}: o grafo nao esta conectado.'.format(url_file))

			else: print('{}: numero incorreto de arestas.'.format(url_file))
