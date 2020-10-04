import networkx as nx
import matplotlib.pyplot as plt
from math import *
import csv
import sys
from pathlib import Path
import copy


WEG = 10**100
RAD = pi / 180
R = 6371


def Floyd_Warshall(matrix, n):                        #Floyd_Warshall algorithm with saving matrix for restoring path
    global lenght, previous
    lenght = [[matrix[i][j] for j in range(n)] for i in range(n)] 
    previous = [[i for j in range(n)] for i in range(n)] 
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if lenght[i][k] + lenght[k][j] < lenght[i][j]:
                    lenght[i][j] = lenght[i][k] + lenght[k][j]
                    previous[i][j] = previous[k][j]


def path(i, j, prev):                               #restoring path from vertex i to j
    ways = [] 
    while j != i: 
        ways.append(j) 
        j = prev[i][j] 
    ways.append(i)
    return ways[::-1]


def res_path(i, j, matrix):                       #reserve path
    n = len(matrix)
    way = path(i, j, previous_main)
    res_mat = copy.deepcopy(matrix)
    ind = 0
    const_way = copy.deepcopy(way)
    if len(way) <= 2:                           #building reserve_matrix in case we cannot repeat vertexes and edges
        res_mat[i][j] = WEG                     
        res_mat[j][i] = WEG
        ind = 1
    else:
        way = way[1:-1]
        for elem in way:
            for k in range(n):
                res_mat[elem][k] = WEG
                res_mat[k][elem] = WEG

    Floyd_Warshall(res_mat, n)                #building the shortest path in case we cannot repeat vertexes 
    wei = lenght[i][j]
    if wei != WEG:
        res_path = path(i, j, previous)
        return [res_path, wei]
    elif ind == 1:
        return ['no', '-']
    '''
    if it is impossible to build way with unic vertexes and edges,
	go through all the variants with one vertex used
	and will build list of answers for all this versions
	'''
    ans = []                               

    for l in range(len(way)):
        res_mat = copy.deepcopy(matrix)
        res_mat[const_way[l]][const_way[l+1]] = WEG                #excluding the used edges from the remaining vertex
        res_mat[const_way[l+1]][const_way[l]] = WEG
        res_mat[const_way[l+2]][const_way[l+1]] = WEG
        res_mat[const_way[l+1]][const_way[l+2]] = WEG
        new_way = way[:l] + way[l+1:]
        for elem in new_way:
            for k in range(n):
                res_mat[elem][k] = WEG         #building reserve_matrix in case we can repeat vertex number way[k]
                res_mat[k][elem] = WEG         #where 'way' is shortest path between i and j vertexes
        Floyd_Warshall(res_mat, n)                 #building a path
        wei = lenght[i][j]
        res_path = path(i, j, previous)
        add = []
        add.append(wei)
        add.append(res_path)
        ans.append(add)
    ans.sort()
    if ans[0][0] == WEG:
        return ['no', '-']
    return ans[0][::-1]



line = sys.argv                                     #input with checking
if (len(line) <= 2):
	print("Incorrect input. Please, write options you need")
	sys.exit(1)
if line[1] != '-t':
	print("input name of file only after -t")
	sys.exit(1)
if (Path(line[2]).exists()):                        #create list of edges and nodes
	name = line[2].split('.')
	if name[1] == 'gml':                     
		G = nx.read_gml(line[2])
		edge = list(G.edges())
		nod = list(G.nodes())
	else:
		G = nx.read_graphml(line[2])
		prom_nod = list(G.nodes(data = 'label'))    
		nod = []
		count = len(prom_nod)
		for i in range(count):
			nod.append(prom_nod[i][1])
		prom_edge = list(G.edges())
		count = len(prom_edge)
		edge = []
		for i in range(count):
			edge.append((nod[int(prom_edge[i][0])], nod[int(prom_edge[i][1])]))
else:
	print("Incorrect name of file")
	sys.exit(1)


#topology parsing
Long = list(G.nodes(data = 'Longitude'))  
Lat = list(G.nodes(data = 'Latitude'))
sity = []          #will contain elements in format ['City', Longitude, Latitude]
l = len(nod)
for i in range(l):
	add = []
	add.append(nod[i])
	add.append(Long[i][1])
	add.append(Lat[i][1])
	sity.append(add)

matrix = [[WEG for i in range(l)] for j in range(l)]    #matrix of weight
for i in range(l):
	matrix[i][i] = 0

#start of first CSV file generation 
data = [['Node 1 (id)','Node 1 (label)', 'Node 1 (longitude)', 'Node 1 (lalitude)', 'Node 2 (id)','Node 2 (label)', 'Node 2 (longitude)', 'Node 2 (lalitude)', 'Distance (km)', 'Delay (mks)']]
for i in range(len(edge)):                    #calculating distances
	data_add = []                             #will contain new string for first CSV file
	j = 0
	while edge[i][0] != sity[j][0] and j < len(edge):
		j += 1
	k = 0
	while edge[i][1] != sity[k][0] and k < len(edge):
		k += 1

	if j == len(edge) or k == len(edge) or sity[j][2] == None or sity[j][1] == None or sity[k][2] == None or sity[k][1] == None:
		continue

	data_add.append(str(j))
	data_add += sity[j]
	data_add.append(str(k))
	data_add += sity[k]                                   

	'''
	https://www.kobzarev.com/programming/calculation-of-distances-between-cities-on-their-coordinates/  - formula for distance between 2 points
	difference is because of need to convert the angle to radians
	''' 
 
	leng = atan2(sqrt((cos(sity[k][2] * RAD)*sin((sity[j][1] - sity[k][1])* RAD))**2 + (cos(sity[j][2] * RAD) * sin(sity[k][2] * RAD) - 
    	sin(sity[j][2] * RAD) * cos(sity[k][2] * RAD) * cos((sity[j][1] - sity[k][1])* RAD))**2),
    	sin(sity[j][2] * RAD) * sin(sity[k][2] * RAD) + cos(sity[k][2] * RAD) * cos(sity[j][2] * RAD) * cos((sity[j][1] - sity[k][1])* RAD))
	leng *= R
	delay = 4.8 * leng                 #delay in mks - weght of node
	matrix[k][j] = delay               #filling the weight matrix
	matrix[j][k] = delay

	data_add.append(leng)
	data_add.append(delay)
	data.append(data_add)              #add a new string in CSV file


#first CSV file 
FILENAME = name[0] + '_topo.csv'
with open(FILENAME, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerows(data)

l = len(nod)
Floyd_Warshall(matrix, l)                    
global lenght_main, previous_main
lenght_main = copy.deepcopy(lenght)      #saving matrix of shortest distances
previous_main = copy.deepcopy(previous)          #saving matrix for recovering paths

#start of second CSV file generation 
data = [['Node 1 (id)', 'Node 2 (id)', 'Path type', 'Path', 'Delay (mks)']]

for i in range(l):                       #add a new strings in CSV file
	for j in range(l):
		if i != j:
			way = path(i, j, previous_main)
			if lenght_main[i][j] == WEG:              
				add = [i, j, 'main', 'no', '-']
			else:
				add = [i, j, 'main', way, lenght_main[i][j]]
			data.append(add)
			if len(line) > 3:                               #add reserve ways in file if need
				reserve_way = res_path(i, j, matrix)
				add = [i, j, 'reserv', reserve_way[0], reserve_way[1]]
				data.append(add)


gr = nx.Graph()                             #create graph that we will use for drawing
gr.add_nodes_from(nod)
for elem in edge:
	gr.add_edge(elem[0],elem[1],color='black',weight=1)

#second CSV file 
FILENAME = name[0] + '_routes.csv'
with open(FILENAME, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerows(data)

if len(line) == 4:                                      #chacking for correctness of input
	if line[3] != '-r':
		print("Incorrect input. \nPlease write '-r' in case you need reserve path, \nor write start and end of path")
		sys.exit(1)
elif ((len(line) == 7) or (len(line) == 8)):
	if line[3] != '-s' or line[5] != '-d':
		print('Incorrect format of input')
		sys.exit(1)
	elif (int(line[4]) > len(sity)) or (int(line[4]) < 0) or (int(line[6]) > len(sity)) or (int(line[6]) < 0):
		print('Incorrect number (ID) of sity')
		sys.exit(1)
	else:                                                                         #output the path between nodes if need
		main_way = path(int(line[4]), int(line[6]), previous_main)
		if lenght_main[int(line[4])][int(line[6])] == WEG:
			print('There is no path')
		else:
			print('Main path:')
			print(*main_way)
			for i in range(len(main_way) - 1):
				gr[nod[main_way[i]]][nod[main_way[i + 1]]]['color'] = 'red'        #coloring and highlight edges in graph 
				gr[nod[main_way[i]]][nod[main_way[i + 1]]]['weight'] = 2
			if len(line) == 8:                                                     #output the reserve path between nodes if need
				print('Reserve path:')
				reserve_way = res_path(int(line[4]), int(line[6]), matrix)
				reserve_way = reserve_way[0]
				if reserve_way == 'no':
					print('there is no reserve way')
				else:
					for i in range(len(reserve_way) - 1):
						gr[nod[reserve_way[i]]][nod[reserve_way[i + 1]]]['color'] = 'green'      #coloring and highlight reserve edges in graph
						gr[nod[reserve_way[i]]][nod[reserve_way[i + 1]]]['weight'] = 2
					print(*reserve_way)
elif len(line) > 3:
	print('Incorrect input')
	sys.exit(1)

pos = nx.circular_layout(gr)                                                    #drawing graph
e = gr.edges()
colors = [gr[u][v]['color'] for u,v in e]
w = [gr[u][v]['weight'] for u,v in e]

nx.draw(gr, pos, edges = e, edge_color = colors, width = w, with_labels = True, horizontalalignment = 'left', node_size = 10)
plt.show() 

sys.exit(0)