import java.io.*;
import java.util.*;

interface nng_strategy{
    public void make_eMST(List<String[]>edges, HashMap<String, Integer>node_indices, HashMap<Integer, String>node_names, double epsilon, File outputFile);
}


class NearestNeighbourGraph_matrix implements nng_strategy{
    // The variable stores the number of nodes in the graph
    private int V = 0;

    public void make_eMST(List<String[]>edges, HashMap<String, Integer>node_indices, HashMap<Integer, String>node_names, double epsilon, File outputFile){
        double[][] graph = get_distance_matrix(edges, node_indices);
        List<HashSet<Integer>> nng = build_nearest_neighbour_graph(graph, epsilon);
        try{
            export_graph(nng, graph, node_names, outputFile);
            }
        catch (Exception e) {
                e.printStackTrace();
            }
    }


    public List<HashSet<Integer>> nearest_neighbour_graph(double[][] graph, int[] mst_parents, double epsilon) {
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

    private List<HashSet<Integer>> build_nearest_neighbour_graph(double[][] graph, double epsilon) {
        MST_matrix t = new MST_matrix();
        int[] mst_parents = t.primMST(graph);

        return nearest_neighbour_graph(graph, mst_parents, epsilon);
    }
    private double[][] get_distance_matrix(List<String[]> edges, HashMap<String, Integer> node_indices) {
        double[][] distance_matrix = new double[node_indices.size()][node_indices.size()];
        for(int i=0; i<distance_matrix.length; ++i)
            for(int j=0; j<distance_matrix.length; ++j)
                distance_matrix[i][j] = -1.0;
        for(String[] a: edges) {
            int u = node_indices.get(a[0]);
            int v = node_indices.get(a[1]);
            double dist = Double.parseDouble(a[2]);
            distance_matrix[u][v] = distance_matrix[v][u] = dist;
        }
        return distance_matrix;
    }
    private void export_graph(List<HashSet<Integer>> adj_list, double[][] graph,
                                     HashMap<Integer, String> node_names,
                                     File file_name) throws FileNotFoundException {
        PrintWriter f = new PrintWriter(file_name);
        f.println("Source,Target,Dist");
        for(int i=0; i<adj_list.size(); ++i)
            for (int j : adj_list.get(i))
                if (i < j) f.println(String.format("%s,%s,%f", node_names.get(i), node_names.get(j), graph[i][j]));
        f.close();
    }
}

class NearestNeighbourGraph_list implements nng_strategy{
    // The variable stores the number of nodes in the graph
    private int V = 0;

    public void make_eMST(List<String[]>edges, HashMap<String, Integer>node_indices, HashMap<Integer, String>node_names, double epsilon, File outputFile){
            ArrayList<ArrayList<Double>> graph = get_distance_graph(edges, node_indices);
            List<HashSet<Integer>> nng = build_nearest_neighbour_graph(graph, epsilon);
            try{
                export_graph(nng, graph, node_names, outputFile);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            
    }

    public List<HashSet<Integer>> nearest_neighbour_graph(ArrayList<ArrayList<Double>> graph, int[] mst_parents, double epsilon) {
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

    private List<HashSet<Integer>> build_nearest_neighbour_graph(ArrayList<ArrayList<Double>> graph, double epsilon) {
        MST_list t = new MST_list();
        int[] mst_parents = t.primMST(graph);

        return nearest_neighbour_graph(graph, mst_parents, epsilon);
    }
    private ArrayList<ArrayList<Double>> get_distance_graph(List<String[]> edges, HashMap<String, Integer> node_indices) {
        ArrayList<ArrayList<Double>> distance_graph = new ArrayList<ArrayList<Double>>(node_indices.size());
        for(int i=0; i<node_indices.size(); ++i){
            ArrayList<Double> alist = new ArrayList<>(node_indices.size());
            for(int j=0; j<node_indices.size(); ++j)
                alist.add(-1.0);
            distance_graph.add(alist);
        }
        //now we have the distance matrix in the form of arraylists 
        for(String[] a: edges) {
            int u = node_indices.get(a[0]);
            int v = node_indices.get(a[1]);
            double dist = Double.parseDouble(a[2]);
            distance_graph.get(u).set(v, dist);
            distance_graph.get(v).set(u, dist);
        }
        return distance_graph;
    }

    private void export_graph(List<HashSet<Integer>> adj_list, ArrayList<ArrayList<Double>> graph,
                                     HashMap<Integer, String> node_names,
                                     File file_name) throws FileNotFoundException {
        PrintWriter f = new PrintWriter(file_name);
        f.println("Source,Target,Dist");
        for(int i=0; i<adj_list.size(); ++i)
            for (int j : adj_list.get(i))
                if (i < j) f.println(String.format("%s,%s,%f", node_names.get(i), node_names.get(j), graph.get(i).get(j)));
        f.close();
    }
}
