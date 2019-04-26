# nn
Nearest Neighbor Network is eMST.

Having genetic distances between viral samples, we want to build a graph where samples are nodes, and genetically related samples are connected by edges. A Minimum Spanning Tree is a good starting point for these purposes. Unfortunately, if distances between viral samples are not unique, there may be many solutions for MST. Having all them determined, we can obtain a joint of Minimal Spanning Trees (jMST). It is a graph that contains all edges that can belong to any MST built on given viral samples. Luckily, jMST contains all edges which connect the closest nodes. But what if we may want to have edges for the 2nd, 3rd, and so on neighbor for every node? In that case, we may specify an epsilon – a value that sets the limit for MST weight increase. Having MST and an epsilon e, we are allowed to remove any edge. Suppose that we chose an edge with weight w. After removal of this edge, we have two disjointed trees which we are going to join back by using another appropriate edge. Any edge is appropriate if it connects these two trees back and have a weight less then w*(1+e). Joining of all possible such trees gives us a joint of epsilon Minimal Spanning Trees (eMST).

Installation:
mvn clean install
