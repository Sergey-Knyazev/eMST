import com.univocity.parsers.csv.CsvParserSettings;
import com.univocity.parsers.csv.CsvParser;
import io.vavr.Tuple2;
import picocli.CommandLine;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

@CommandLine.Command(name = "nn", mixinStandardHelpOptions = true, version = "0.0")
public class Main implements Runnable{
    @CommandLine.Option(names={"-i", "--inFile"}, description="input file with edges in CSV format",
            paramLabel = "FILE", required=true)
    private File inputFile;
    @CommandLine.Option(names={"-o", "--outFile"},
            description="output file with edges of nearest neighbour graph in CSV format",
            paramLabel = "FILE", required=true)
    private File outputFile;
    @CommandLine.Option(names={"-e", "--epsilon"},
            description="epsilon - an edge will be added in nn if its weight is only 1+e times greater than needed",
            paramLabel = "e")
    private double epsilon = 0.0;

    public void run() {
        try {
            List<String[]> edges = load_edges(inputFile);
            Tuple2<HashMap<String, Integer>, HashMap<Integer, String>> nodes = get_node_indices(edges);
            HashMap<String, Integer> node_indices = nodes._1;
            HashMap<Integer, String> node_names = nodes._2;
            double[][] graph = get_distance_matrix(edges, node_indices);
            List<HashSet<Integer>> nng = build_nearest_neighbour_graph(graph, epsilon);
            export_graph(nng, graph, node_names, outputFile);
        }
        catch(FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        CommandLine.run(new Main(), System.out, args);
    }
    private static List<HashSet<Integer>> build_nearest_neighbour_graph(double[][] graph, double epsilon) {
        MST t = new MST();
        int[] mst_parents = t.primMST(graph);
        NearestNeighbourGraph g = new NearestNeighbourGraph();
        return g.nearest_neighbour_graph(graph, mst_parents, epsilon);
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
    private static double[][] get_distance_matrix(List<String[]> edges, HashMap<String, Integer> node_indices) {
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
    private static void export_graph(List<HashSet<Integer>> adj_list, double[][] graph,
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
