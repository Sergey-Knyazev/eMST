import java.util.*;

class NearestNeighbourGraph {
    // The variable stores the number of nodes in the graph
    private int V = 0;
    List<HashSet<Integer>> nearest_neighbour_graph(ArrayList<ArrayList<Double>> graph, int[] mst_parents, double epsilon) {
        V = graph.size();
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
        ArrayList<ArrayList<Double>> longest_edge = new ArrayList<ArrayList<Double>>();
        // ArrayList<Double> alist = new ArrayList<Double>();
        for(int i=0; i<V; ++i){
            // alist = new ArrayList<Double>() ;
            longest_edge.add(new ArrayList<Double>());
            for(int j=0; j<V; ++j)
                // alist.add(0.0);
                longest_edge.get(i).add(0.0);
        }
        for(int i=0; i<V; ++i){
            bfs_update_matrix(mst, graph, i, longest_edge);
        }
        for(int i=0; i<V; ++i) {
            for(int j=0; j<V; ++j) {
                if(graph.get(i).get(j) > 0 && graph.get(i).get(j) <= longest_edge.get(i).get(j)*(1.0+epsilon)) {
                    nng.get(i).add(j);
                    nng.get(j).add(i);
                }
            }
        }
        return nng;
    }
    private void bfs_update_matrix(List<HashSet<Integer>> mst, ArrayList<ArrayList<Double>> weights, int root, ArrayList<ArrayList<Double>> longest_edge) {
        boolean[] visited = new boolean[V]; //False by default
        Queue<Integer> queue = new LinkedList<Integer>();
        queue.add(root);
        while(!queue.isEmpty()) {
            int v = queue.remove();
            visited[v] = true;
            for(Integer u: mst.get(v)) {
                if(visited[u]) continue;
                queue.add(u);
                longest_edge.get(root).set(u, Math.max(weights.get(v).get(u), Math.max(longest_edge.get(root).get(u), longest_edge.get(root).get(v))));
                longest_edge.get(u).set(root, Math.max(weights.get(v).get(u), Math.max(longest_edge.get(root).get(u), longest_edge.get(root).get(v))));
            }
        }
    }
}
