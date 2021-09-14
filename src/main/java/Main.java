import com.univocity.parsers.csv.CsvParserSettings;
import com.univocity.parsers.csv.CsvParser;
import io.vavr.Tuple2;
import picocli.CommandLine;
import TN93.TN93;

// import org.jgrapht.*;
// import org.jgrapht.graph.*;
// import org.jgrapht.nio.*;
// import org.jgrapht.nio.dot.*;
// import org.jgrapht.traverse.*;
// import org.jgrapht.alg.connectivity.*;


import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.ArrayList;


@CommandLine.Command(name = "nn", mixinStandardHelpOptions = true, version = "0.0")
public class Main implements Runnable{
    @CommandLine.Option(names={"-i", "--inFile"}, description="input file with edges in CSV format or fasta file(will work only if a=2)",
            paramLabel = "FILE", required=true)
    private File inputFile = new File("infile.defaultfile");
    @CommandLine.Option(names={"-o", "--outFile"},
            description="output file with edges of nearest neighbour graph in CSV format",
            paramLabel = "FILE")
    private File outputFile = new File("outfile.defaultfile");
    @CommandLine.Option(names={"-e", "--epsilon"},
            description="epsilon - an edge will be added in nn if its weight is only 1+e times greater than needed",
            paramLabel = "e")
    private double epsilon = 0.0;
    @CommandLine.Option(names={"-t", "--threshold"},
    description="threshold - remove all edges with length > threshold.",
    paramLabel = "t")
    private double edge_threshold = -1.0;   
    @CommandLine.Option(names={"-a", "--algorithm"},
            description="algorithm - whether to use matrices(a=0), Lists(a=1), on-the-fly-algo(a=2), print the number of components in the graph(a=3), or calculate the pairwise hamming distances file(a=4)",
            paramLabel = "a")
    private int algorithm = 1;
    @CommandLine.Option(names={"-d", "--distance"},
            description="distance - Hamming distance(d=1) [default], TN93(d=1)",
            paramLabel = "d")
    private int distance_metric = 1;

    public void run() {

        try {

            // double[][] graph = get_distance_matrix(edges, node_indices);
            if(algorithm == 0){
                System.out.println("You're using the matrix algorithm(a=0)");

                List<String[]> edges = load_edges(inputFile);
                Tuple2<HashMap<String, Integer>, HashMap<Integer, String>> nodes = get_node_indices(edges);
                HashMap<String, Integer> node_indices = nodes._1;
                HashMap<Integer, String> node_names = nodes._2;

                NearestNeighbourGraph_matrix g = new NearestNeighbourGraph_matrix();
                g.make_eMST(edges, node_indices, node_names, epsilon, outputFile);
            }
            else if(algorithm == 1){
                System.out.println("You're using the list algorithm(a=1)");

                // List<String[]> all_edges = load_edges(inputFile);
                // Tuple2<HashMap<String, Integer>, HashMap<Integer, String>> all_nodes = get_node_indices(all_edges);
                // HashMap<String, Integer> all_node_indices = all_nodes._1;
                // HashMap<Integer, String> all_node_names = all_nodes._2;

                // Make a list of list of edges, ie each component of the graph will be a lius tof edge and all components will be a part of a single list
                // List<List<String[]>> edges_multi_components = get_multiple_components(all_edges, all_node_indices, all_node_names);
                
                // for(List<String[]> comp_edges: edges_multi_components){

                List<String[]> edges = load_edges(inputFile);

                // Tuple2<HashMap<String, Integer>, HashMap<Integer, String>> nodes = get_node_indices(comp_edges);
                Tuple2<HashMap<String, Integer>, HashMap<Integer, String>> nodes = get_node_indices(edges);
                HashMap<String, Integer> node_indices = nodes._1;
                HashMap<Integer, String> node_names = nodes._2;

                NearestNeighbourGraph_list g = new NearestNeighbourGraph_list();
                // g.make_eMST(comp_edges, node_indices, node_names, epsilon, outputFile);
                g.make_eMST(edges, node_indices, node_names, epsilon, outputFile);

                // }
            }
            else if(algorithm == 2){
                NearestNeighbourGraph_fasta g = new NearestNeighbourGraph_fasta();
                System.out.println("You're inputting a fasta file and using the on-the-fly algorithm(a=2)");
                g.make_eMST(epsilon, outputFile, inputFile, distance_metric, edge_threshold);
            }
            else if(algorithm == 3){

                // returns the number of connected components in the graph

                List<String[]> edges = load_edges(inputFile);
                Tuple2<HashMap<String, Integer>, HashMap<Integer, String>> nodes = get_node_indices(edges);
                HashMap<String, Integer> node_indices = nodes._1;
                HashMap<Integer, String> node_names = nodes._2;

                ArrayList<ArrayList<Double>> graph = get_distance_graph(edges, node_indices);
                Graph multi_comp_graph = new Graph(graph, node_indices, node_names);
                // get the number of connected components of the Graph object
                ConnectedComponents conn_comp_obj = new ConnectedComponents();
                int num_components = conn_comp_obj.numberOfConnecedComponents(multi_comp_graph);

                System.out.println(num_components);
            }
            else if(algorithm == 4){
                Ham_dist testobj = new Ham_dist();
                testobj.setInputFile(inputFile);
                // File file = new File("./outputTN93test.csv");
                testobj.setOutputFile(outputFile);
                testobj.hammingFasta();
            }
            else if(algorithm == 5){
                TN93 testobj = new TN93();
                testobj.setInputFile(inputFile);
                // File file = new File("./outputTN93test.csv");
                testobj.setOutputFile(outputFile);
                testobj.tn93Fasta();
            }

        }
        catch(FileNotFoundException e) {
            e.printStackTrace();
        }
    }
    
    public static void main(String[] args) {
        CommandLine.run(new Main(), System.out, args);
    }

    // private static List<List<String[]>> get_multiple_components(List<String[]> edges, HashMap<String, Integer> node_indices, HashMap<Integer, String> node_names){

    //     ArrayList<ArrayList<Double>> full_graph = get_full_distance_graph(edges, node_indices);
    //     List<List<String[]>> multi_component_graph = connectedComponents(full_graph, node_names);

    //     return multi_component_graph;

    // }
    // private static List<List<String[]>> DFSUtil(int v, boolean[] visited, ArrayList<ArrayList<Double>> full_graph, HashMap<Integer, String> node_names, List<String[]> component)
    // {
    //     // Mark the current node as visited and print it
    //     visited[v] = true;
    //     System.out.print(v + " ");
    //     // Recur for all the vertices
    //     // adjacent to this vertex
    //     for (int x : full_graph.get(v)) {
    //         if (!visited[x])

    //             component = DFSUtil(x, visited, full_graph, node_names, component);
    //     }

    //     return component;
    // }
    // private static List<List<String[]>> connectedComponents(ArrayList<ArrayList<Double>> full_graph, HashMap<Integer, String> node_names)
    // {
    //     List<List<String[]>> multi_component_graph = new ArrayList<List<String[]>>();
    //     // Mark all the vertices as not visited
    //     boolean[] visited = new boolean[node_names.size()];
    //     for (int v = 0; v < V; ++v) {
    //         if (!visited[v]) {
    //             List<String[]> component_v = new ArrayList<String[]>();
    //             // print all reachable vertices
    //             // from v
    //             List<String[]> component = DFSUtil(v, visited, full_graph, node_names, component_v);
    //             // System.out.println();
    //             multi_component_graph.add(component_v);

    //         }
    //     }

    //     return multi_component_graph;
    // }



    // ***********************************************



    // private void get_graph_string(List<HashSet<Integer>> adj_list, ArrayList<ArrayList<Double>> graph, HashMap<Integer, String> node_names){
    //     for(int i=0; i<adj_list.size(); ++i)
    //         for (int j : adj_list.get(i))
    //             if (i < j) f.println(String.format("%s,%s,%f", node_names.get(i), node_names.get(j), graph.get(i).get(j)));
    // }

    private static ArrayList<ArrayList<Double>> get_full_distance_graph(List<String[]> edges, HashMap<String, Integer> node_indices) {
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

    private static List<String[]> load_edges(File file_name) throws FileNotFoundException{
        CsvParserSettings settings= new CsvParserSettings();
        settings.getFormat().setLineSeparator("\n");
        settings.setHeaderExtractionEnabled(true);
        CsvParser parser = new CsvParser(settings);
        return parser.parseAll(new FileInputStream(file_name));
    }
    private static Tuple2<HashMap<String, Integer>, HashMap<Integer,String>> get_node_indices(List<String[]> edges) {
        HashMap<String, Integer> node_indices = new HashMap<String, Integer>();
        HashMap<Integer, String> node_names = new HashMap<Integer, String>();
        for(String[] a: edges) {
            for (int i: new int[]{0, 1}) {
               if(node_indices.containsKey(a[i])) continue;
               int idx = node_indices.size();
               node_indices.put(a[i], idx);
               node_names.put(idx, a[i]);
            }
        }
        return new Tuple2<HashMap<String, Integer>, HashMap<Integer, String>>(node_indices, node_names);
    }

    private static ArrayList<ArrayList<Double>> get_distance_graph(List<String[]> edges, HashMap<String, Integer> node_indices) {
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
}
