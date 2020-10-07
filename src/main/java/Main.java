import com.univocity.parsers.csv.CsvParserSettings;
import com.univocity.parsers.csv.CsvParser;
import io.vavr.Tuple2;
import picocli.CommandLine;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.ArrayList;

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
    @CommandLine.Option(names={"-a", "--algorithm"},
            description="algorithm - whether to use matrices(a=0) or Lists(a=1)",
            paramLabel = "a")
    private int algorithm = 1;

    public void run() {
        try {
            List<String[]> edges = load_edges(inputFile);
            Tuple2<HashMap<String, Integer>, HashMap<Integer, String>> nodes = get_node_indices(edges);
            HashMap<String, Integer> node_indices = nodes._1;
            HashMap<Integer, String> node_names = nodes._2;
            // double[][] graph = get_distance_matrix(edges, node_indices);
            if(algorithm == 0){
                NearestNeighbourGraph_matrix g = new NearestNeighbourGraph_matrix();
                g.make_eMST(edges, node_indices, node_names, epsilon, outputFile);
                System.out.println("You're using the matrix algorithm(a=0)");
            }
            else{
                NearestNeighbourGraph_list g = new NearestNeighbourGraph_list();
                g.make_eMST(edges, node_indices, node_names, epsilon, outputFile);
                System.out.println("You're using the list algorithm(a=1)");
            }

        }
        catch(FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        CommandLine.run(new Main(), System.out, args);
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

}
