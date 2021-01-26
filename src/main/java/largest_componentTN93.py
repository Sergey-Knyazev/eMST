#USE AS FOLLOWS
# python largest_componentTN93.py tn93path.csv output_largestcomponentpath.csv


import os
import sys
import networkx as nx
import pandas as pd
import csv

TN93path = sys.argv[1]
outpath = sys.argv[2]

#Load the TN93 graph into networkx
GraphData = open(TN93path, 'r')
next(GraphData, None)#for ignoring the  first line
Graphtype = nx.Graph()
G= nx.parse_edgelist(GraphData, delimiter=',', create_using = Graphtype, nodetype=str, data=(('weight',float),))

#split the graph into individal connected components
connected_components = sorted(nx.connected_components(G), key=len, reverse=True)
Connected_SubGrpahs = [G.subgraph(c).copy() for c in connected_components] #Connected components in list format in decreasing order

# #write all components of the graphs to csv files
# for i in range(len(Connected_SubGrpahs)):
#   nx.readwrite.write_edgelist(Connected_SubGrpahs[i], TN93path+'_largest_comp', delimiter=',', data=True, encoding='utf-8')

# write largest component of the graph to csv file
fields = ['Source', 'target', 'weight']  
with open(outpath, 'w') as csvfile:   
    csvwriter = csv.writer(csvfile)   
    csvwriter.writerow(fields)  
    for u,v,a in Connected_SubGrpahs[0].edges(data=True):
        csvwriter.writerow([u, v, a['weight']])     
# nx.readwrite.write_edgelist(Connected_SubGrpahs[0], TN93path[:-4]+'_largest_comp.csv', delimiter=',', data=True, encoding='utf-8')



# # create dictionary of metadata of the form where each key - value pair is name:{'continent':continent, 'country':country, 'division':division ......other attributes and their value ........} 
# MetaData_Index_num = pd.read_csv('../AlignedSequencesFull/metadata.csv')
# MetaData = MetaData_Index_num.set_index('strain')
# MetaData_Dict = MetaData.to_dict('Index')
# print('metadata loaded')
# epsilon = [('0', 0), ('0.00001', 0.00001), ('0.00002', 0.00002), ('0.00005', 0.00005), ('0.0001', 0.0001), ('0.0002', 0.0002), ('0.0005', 0.0005), ('0.001', 0.001), ('0.002', 0.002), ('0.005', 0.005), ('0.01', 0.01), ('0.02', 0.02), ('0.05', 0.05), ('0.1', 0.1)]
  
# with open('./nextstrain_168/maskedAssortativityTable.csv', 'w') as csvfile:
#   writer = csv.writer(csvfile)
#   writer.writerow(['Epsilon','Degree', 'Continent', 'Country', 'Division'])
#   for e in epsilon:
#     print(str(e))
#     # create graph 
#     GraphData = open('./nextstrain_168/168_nextstrain_e'+ e[0] +'.csv', 'r')
#     next(GraphData, None)#for ignoring the  first line
#     Graphtype = nx.Graph()
#     G= nx.parse_edgelist(GraphData, delimiter=',', create_using = Graphtype, nodetype=str, data=(('weight',float),))

#     # for all nodes in graph, get attributes using Metadata_Dict
#     for node in G.nodes:
#       nx.set_node_attributes(G, {node : MetaData_Dict[node]})

#     deg_assortativity = nx.degree_assortativity_coefficient(G)
#     continent_assortavity = nx.attribute_assortativity_coefficient(G, 'region')
#     country_assortativity = nx.attribute_assortativity_coefficient(G, 'country')
#     division_assortavity = nx.attribute_assortativity_coefficient(G, 'division')

#     # write results for this graph in csv
#     writer.writerow([e[1], deg_assortativity, continent_assortavity, country_assortativity ,division_assortavity])

# df = pd.read_csv('./nextstrain_168/maskedAssortativityTable.csv')
# print(df.head)

# # Make 4 graphs and save them
# deg = df.plot(x = 'Epsilon', y = 'Degree', logx = True)
# cont = df.plot(x = 'Epsilon', y = 'Continent', logx = True)
# count = df.plot(x = 'Epsilon', y = 'Country', logx = True)
# div = df.plot(x = 'Epsilon', y = 'Division', logx = True)

# figdeg = deg.get_figure()
# figdeg.savefig("./nextstrain_168/degree_graph.jpg")

# figcont = cont.get_figure()
# figcont.savefig("./nextstrain_168/continent_graph.jpg")

# figcount = count.get_figure()
# figcount.savefig("./nextstrain_168/country_graph.jpg")

# figdiv = div.get_figure()
# figdiv.savefig("./nextstrain_168/div_graph.jpg")