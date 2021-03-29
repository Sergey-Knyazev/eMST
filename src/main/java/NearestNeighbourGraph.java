import java.io.*;
import java.util.*;

import io.vavr.Tuple2;

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
            // ArrayList<ArrayList<Integer>> graph_node_connectivity = get_node_con_graph(edges, node_indices);

            // Make a graph object haviong the graph in adjacency list format, node_indices and node_names.
            Graph multi_comp_graph = new Graph(graph, node_indices, node_names);
            // get the connected components of the Graph object
            ConnectedComponents conn_comp_obj = new ConnectedComponents();
            List<Graph> graph_connected_components = conn_comp_obj.getConnectedComponents(multi_comp_graph);

            // NOTE: What can be done is to shift some functionss from this file(NNG class) to graph class
            for(Graph graph_subcomponent: graph_connected_components){
                List<HashSet<Integer>> nng = build_nearest_neighbour_graph(graph_subcomponent.get_graph(), epsilon);
                try{
                    export_graph(nng, graph_subcomponent.get_graph(), graph_subcomponent.get_node_names(), outputFile);
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }

            // List<HashSet<Integer>> nng = build_nearest_neighbour_graph(graph, epsilon);
            // try{
            //     export_graph(nng, graph, node_names, outputFile);
            // }
            // catch (Exception e) {
            //     e.printStackTrace();
            // }
            
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
        PrintWriter f = new PrintWriter(new FileOutputStream(file_name,true));
        if(file_name.length()==0){
            f.println("Source,Target,Dist");
        }
        for(int i=0; i<adj_list.size(); ++i)
            for (int j : adj_list.get(i))
                if (i < j) f.println(String.format("%s,%s,%f", node_names.get(i), node_names.get(j), graph.get(i).get(j)));
        f.close();
    }
}

class NearestNeighbourGraph_fasta{

    private int V = 0;
    private int edgeThreshold = 10000;

    public void make_eMST(double epsilon, File outputFile, File alnfile){

        try {
            HashMap<Integer, String> node_seqnames = get_node_seqnames(alnfile);
            HashMap<Integer, Seq> node_sequences = get_node_sequences(alnfile);
            
            // TODO: get consensus using node_sequences hashmap, find positions of 
            // change from all other sequeces and store it in anotehr hashmap node_diff_consensus<Integer, ArrayList<Integer>>, 
            // then into the build and export function, pass node_seqnames as before, node_diff_consensus, epsilon and outputfile as before
            
            Seq consensus = get_consensus_strain(node_sequences);
            HashMap<Integer, ArrayList<Integer>> diff_consensus = get_diff_consensus(node_sequences, consensus);
            
            buildandexport_nearest_neighbour_graph(diff_consensus, node_sequences, node_seqnames, epsilon, outputFile);
        
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private HashMap<Integer, ArrayList<Integer>> get_diff_consensus(HashMap<Integer, Seq> sequences, Seq consensus){
        
        //the hashmap storing the poitions of differences of each sequence from the consensus
        HashMap<Integer, ArrayList<Integer>> diff_consensus = new HashMap<Integer, ArrayList<Integer>>();

        int[] consensus_seqenc = consensus.getSeq_enc();
        int t = 0;
        for(Map.Entry<Integer,Seq> entry : sequences.entrySet()){
            int[] seq = entry.getValue().getSeq_enc();
            ArrayList<Integer> store_diffs = new ArrayList<Integer>();

            for(int i = 0; i<seq.length; i++){
                if(seq[i]!=consensus_seqenc[i]){
                    store_diffs.add(i);
                }
            }
            diff_consensus.put(t, store_diffs);
            t+=1;
        }

        return diff_consensus;
    }

    private Seq get_consensus_strain(HashMap<Integer, Seq> sequences){
        int seqlen = sequences.get(0).get_seq_len();
        ArrayList<ArrayList<Integer>> counts = new ArrayList<ArrayList<Integer>>(seqlen);
        for(int i = 0; i<seqlen; i++){
            // if(i%1000 == 0){
            //     System.out.println(i);
            // }
            ArrayList<Integer> base_count = new ArrayList<Integer>(5);
            for(int j=0; j<5; j++){
                base_count.add(0);
            }
            counts.add(base_count);//5 because we have to keep track of counts of all 4 base pairs + 1(N) and then make a consensus at the end
        }
        for(Map.Entry<Integer,Seq> entry : sequences.entrySet()){
            int[] seq = entry.getValue().getSeq_enc();
            for(int i = 0; i<seq.length; i++){
                //seq[i] is the i'th base pair -1,2,3 or 4 and if its for eg 2 then we have to update the counts apropriately
                int newval = counts.get(i).get(seq[i])+1;
                counts.get(i).set(seq[i], newval);
            }
        }
        //we have the counts now select the max number base pair for all positions of the genome
        int[] consensus = new int[seqlen];
        for(int i = 0; i<seqlen; i++){
            consensus[i] = counts.get(i).indexOf(Collections.max(counts.get(i)));
        }
        Seq consensus_seq = new Seq("consensus", "");
        consensus_seq.set_seq_enc(consensus);
        consensus_seq.set_sequence();

        return consensus_seq;
    }

    private HashMap<Integer, String> get_node_seqnames(File alnFile) throws FileNotFoundException{
        
        Scanner sc = new Scanner(alnFile);

        HashMap<Integer, String> seqs = new HashMap<Integer, String>();
        String name="", seq="";
        while(sc.hasNextLine()) {
            String line = sc.nextLine().trim();
            if(line.length() == 0) continue;
            if(line.charAt(0)=='>') {
                if (name.length()!=0) seqs.put(seqs.size(), name);
                name = line.substring(1);
                seq="";
            }
            else seq=seq.concat(line);
        }
        if(name.length()!=0) seqs.put(seqs.size(), name);
        
        sc.close();

        return seqs;

    }
    private HashMap<Integer, Seq> get_node_sequences(File alnFile) throws FileNotFoundException{
        
        Scanner sc = new Scanner(alnFile);

        HashMap<Integer, Seq> seqs = new HashMap<Integer, Seq>();
        String name="", seq="";
        while(sc.hasNextLine()) {
            String line = sc.nextLine().trim();
            if(line.length() == 0) continue;
            if(line.charAt(0)=='>') {
                if (name.length()!=0) seqs.put(seqs.size(), new Seq(name, seq));
                name = line.substring(1);
                seq="";
            }
            else seq=seq.concat(line);
        }
        if(name.length()!=0) seqs.put(seqs.size(), new Seq(name, seq));
        
        sc.close();

        return seqs;

    }
    public void nearest_neighbour_graph(HashMap<Integer, ArrayList<Integer>> diff_consensus, HashMap<Integer, Seq> node_sequences, HashMap<Integer, String> node_seqnames, int[] mst_parents, double epsilon, File outputFile) throws FileNotFoundException{
        V = node_sequences.size();
        // build mst
        List<HashSet<Integer>> mst = new ArrayList<HashSet<Integer>>();
        for(int i=0; i<V; ++i)
            mst.add(new HashSet<Integer>());

        for(int i=1; i<V; ++i) {
            mst.get(i).add(mst_parents[i]);
            mst.get(mst_parents[i]).add(i);
        }

        PrintWriter f = new PrintWriter(outputFile);
        f.print("Source,Target,Dist\n");
        f.close();

        for(int i=0; i<V; ++i){
            System.out.println(i);
            ArrayList<Double> longest_edges = new ArrayList<Double>();
            for(int j=0; j<V; j++){
                longest_edges.add(0.0);
            }
            bfs_update_longedges(i, mst, longest_edges, diff_consensus, node_sequences);

            for(int j=0; j<V;j++){

                double dist_i_j = distance_hamming_using_consensus(i, j, diff_consensus, node_sequences);

                // if(node_seqnames.get(i).equals("Switzerland/VD0503/2020") && node_seqnames.get(j).equals("Switzerland/BE2536/2020")){
                //     System.out.println("yes I got it - i = " + i + "j = " + j + "dist = " + dist_i_j + "longestedge = " + longest_edges.get(j));
                // }

                // System.out.println(dist_i_j + "," + longest_edges.get(j));
                if(dist_i_j > 0 && dist_i_j <= (1.0 + epsilon)*longest_edges.get(j)){
                    write_edge(i, j, node_seqnames, dist_i_j, outputFile);
                }
            }

        }
    }

    private void write_edge(int i, int j, HashMap<Integer, String> node_seqnames, double distance, File file_name){
        try {
            FileWriter f = new FileWriter(file_name, true);
            // if(node_seqnames.get(i).equals("Switzerland/VD0503/2020") && node_seqnames.get(j).equals("Switzerland/BE2536/2020")){
            //     System.out.println("yes I got it - i = " + i + "j = " + j + "dist = " + distance);
            // }
            if (i < j && distance<edgeThreshold){
                f.write(String.format("%s,%s,%f\n", node_seqnames.get(i), node_seqnames.get(j), distance));
            }
            f.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
    private void bfs_update_longedges(int root, List<HashSet<Integer>> mst, ArrayList<Double> longest_edges, HashMap<Integer,ArrayList<Integer>> diff_consensus, HashMap<Integer, Seq> node_sequences) {
        
        boolean[] visited = new boolean[V]; //False by default
        Queue<Integer> queue = new LinkedList<Integer>();
        queue.add(root);
        while(!queue.isEmpty()) {
            int v = queue.remove();
            visited[v] = true;
            for(Integer u: mst.get(v)) {
                if(visited[u]) continue;
                queue.add(u);
                longest_edges.set(u, Math.max(distance_hamming_using_consensus(u, v, diff_consensus, node_sequences), Math.max(longest_edges.get(u), longest_edges.get(v))));
                // longest_edge.get(u).set(root, Math.max(weights.get(v).get(u), Math.max(longest_edge.get(root).get(u), longest_edge.get(root).get(v))));
            }
        }
    }

    private double distance_hamming_using_consensus(int u, int v, HashMap<Integer, ArrayList<Integer>> diff_consensus, HashMap<Integer, Seq> node_sequences){
        ArrayList<Integer> diff_1 = diff_consensus.get(u);
        ArrayList<Integer> diff_2 = diff_consensus.get(v);

        Seq seq_1 = node_sequences.get(u);
        Seq seq_2 = node_sequences.get(v);

        double hamdist = 0;

        int i=0;
        int j=0;
        // get distance using diff_1, diff_2 and node_sequences(which will be used when 2 numbers are same in the diff arrays, to check is the bas pairs are same or not at that place)
        while(i<diff_1.size() && j<diff_2.size()){
            if(diff_1.get(i)<diff_2.get(j)){
                i++;
                hamdist++;
            }
            else if(diff_1.get(i)>diff_2.get(j)){
                j++;
                hamdist++;
            }
            else{
                //positions of dissimilarity with consensus is same so we have to check if the base pairs are different or same at this position
                if(seq_1.getSeq_enc()[diff_1.get(i)] != seq_2.getSeq_enc()[diff_2.get(j)]){
                    hamdist++;
                }
                i++;
                j++;
            }
        }

        if(i==diff_1.size()){
            hamdist+=(diff_2.size()-j);
        }
        else{
            hamdist+=(diff_1.size()-i);
        }

        return hamdist;
    }

    // private double tn93_distance(int u, int v, HashMap<Integer, Seq> node_sequences){

    //     Seq s1 = node_sequences.get(u);
    //     Seq s2 = node_sequences.get(v);
    //     return TN93.tn93(s1.getSeq_enc(), s2.getSeq_enc());
    // }
    private void buildandexport_nearest_neighbour_graph(HashMap<Integer, ArrayList<Integer>> diff_consensus, HashMap<Integer, Seq> node_sequences, HashMap<Integer, String> node_seqnames, double epsilon,
            File outputFile) throws FileNotFoundException{
        MST_fasta t = new MST_fasta();
        int[] mst_parents = t.primMST(diff_consensus, node_sequences);
        nearest_neighbour_graph(diff_consensus, node_sequences, node_seqnames, mst_parents, epsilon, outputFile);
    }

}
