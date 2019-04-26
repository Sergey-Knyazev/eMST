import java.util.*;

class NearestNeighbourGraph {
    // The variable stores the number of nodes in the graph
    private int V = 0;
    List<HashSet<Integer>> nearest_neighbour_graph(double[][] graph, int[] mst_parents, double epsilon) {
        V = graph.length;
        // build mst
        List<HashSet<Integer>> mst = new ArrayList<HashSet<Integer>>();
        for(int i=0; i<V; ++i)
            mst.add(new HashSet<Integer>());
        for(int i=1; i<V; ++i) {
            mst.get(i).add(mst_parents[i]);
            mst.get(mst_parents[i]).add(i);
        }
        List<HashSet<Integer>> nng = new ArrayList<HashSet<Integer>>();
        for(int i=0; i<V; ++i)
            nng.add(new HashSet<Integer>());
        /*
        // Add nearest neighbours to every node
        for(int i=0; i<V; ++i) {
            double closest_neighbour_dist = Double.MAX_VALUE;
            for(int j=0; j<V; ++j)
                if(i!=j && graph[i][j]>0 && graph[i][j] < closest_neighbour_dist)
                    closest_neighbour_dist = graph[i][j];
            for(int j=0; j<V; ++j)
                if(i!=j && graph[i][j]>=closest_neighbour_dist && graph[i][j] <= closest_neighbour_dist*(1.0+epsilon)) {
                    nng.get(i).add(j);
                    nng.get(j).add(i);
                }
        }
        // Add MST to NN
        for(int i=0; i<V; ++i) {
            for(int a: mst.get(i)) {
                nng.get(i).add(a);
                nng.get(a).add(i);
            }
        }*/
        // Add bridges between NN components
        double[][] longest_edge = new double[V][V];
        for(int i=0; i<V; ++i)
            for(int j=0; j<V; ++j)
                longest_edge[i][j] = 0;
        for(int i=0; i<V; ++i)
            bfs_update_matrix(mst, graph, i, longest_edge);
        for(int i=0; i<V; ++i) {
            for(int j=0; j<V; ++j) {
                if(graph[i][j] > 0 && graph[i][j] <= longest_edge[i][j]*(1.0+epsilon)) {
                    nng.get(i).add(j);
                    nng.get(j).add(i);
                }
            }
        }
        return nng;
    }
    private void bfs_update_matrix(List<HashSet<Integer>> mst, double[][] weights, int root, double[][] longest_edge) {
        boolean[] visited = new boolean[V];
        Queue<Integer> queue = new LinkedList<Integer>();
        queue.add(root);
        while(!queue.isEmpty()) {
            int v = queue.remove();
            visited[v] = true;
            for(Integer u: mst.get(v)) {
                if(visited[u]) continue;
                queue.add(u);
                longest_edge[root][u] = longest_edge[u][root] =
                        Math.max(weights[v][u], Math.max(longest_edge[root][u], longest_edge[root][v]));
            }
        }
    }
}
