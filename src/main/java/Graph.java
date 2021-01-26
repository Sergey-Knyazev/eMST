import java.util.ArrayList;
import java.io.*;
import java.util.*;

class Graph {
    int V;
    ArrayList<ArrayList<Double>> full_graph;
    HashMap<String, Integer> node_indices;
    HashMap<Integer, String> node_names;

    Graph(ArrayList<ArrayList<Double>> full_graph, HashMap<String, Integer> node_indices, HashMap<Integer, String> node_names){
        this.node_indices = node_indices;
        this.node_names = node_names;
        this.full_graph = full_graph;
        this.V = node_indices.size();
    }

    public ArrayList<ArrayList<Double>> get_graph(){
        return this.full_graph;
    }

    public HashMap<Integer, String> get_node_names(){
        return this.node_names;
    }

    public HashMap<String, Integer> get_node_indices(){
        return this.node_indices;
    }
}

