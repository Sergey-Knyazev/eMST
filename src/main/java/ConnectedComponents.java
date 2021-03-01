import java.util.*;

class ConnectedComponents {


    public List<Graph> getConnectedComponents(Graph graph){
        ArrayList<ArrayList<Double>> full_graph = graph.get_graph();
        HashMap<String, Integer> full_node_indices = graph.get_node_indices();
        HashMap<Integer, String> full_node_names = graph.get_node_names();


        List<Graph> multi_component_graph = new ArrayList<Graph>();
        // Mark all the vertices as not visited
        boolean[] visited = new boolean[full_node_names.size()];
        for (int v = 0; v < full_node_names.size(); ++v) {
            if (!visited[v]) {
                Graph component = BFS(v, visited, full_graph, full_node_names, full_node_indices);
                multi_component_graph.add(component);
            }
        }

        return multi_component_graph;

    }

    public int numberOfConnecedComponents(Graph graph){

        List<Graph> multi_component_graph = getConnectedComponents(graph);
        int num_components = multi_component_graph.size();
        return num_components;
    }

    private static Graph BFS(int v, boolean[] visited, ArrayList<ArrayList<Double>> full_graph, HashMap<Integer, String> node_names, HashMap<String, Integer> node_indices)
    {

        HashMap<String, Integer> sub_nodeindices = new HashMap<String, Integer>();
        HashMap<Integer, String> sub_nodenames = new HashMap<Integer, String>();
        ArrayList<Integer> nodelist = new ArrayList<Integer>();

        // Adding v to the list of nodes
        // System.out.println("-------------");
        sub_nodenames.put(0, node_names.get(v));
        sub_nodeindices.put(node_names.get(v), 0);
        nodelist.add(v);

        int nodenum = 1;
        Queue<Integer> queue = new LinkedList<Integer>();
        queue.add(v);
        visited[v]=true;
        while(!queue.isEmpty()) {
            int x = queue.remove();
            // visited[x] = true;
            for(int u = 0; u < node_names.size(); u++) {
                if(full_graph.get(x).get(u) < 0.0) continue; // no need to add u to the queue if its not connected to x
                if(visited[u]) continue;
                queue.add(u);
                visited[u] = true;
                sub_nodenames.put(nodenum, node_names.get(u));
                sub_nodeindices.put(node_names.get(u), nodenum);
                nodelist.add(u);
                nodenum++;
                // System.out.println(node_names.get(u));
            }
        }
        // System.out.println(sub_nodenames);
        // System.out.println(nodelist);


        ArrayList<ArrayList<Double>> subgraph = new ArrayList<ArrayList<Double>>(nodelist.size());
        for(int i = 0; i<nodelist.size(); ++i){
            ArrayList<Double> alist = new ArrayList<>(nodelist.size());
            for(int j=0; j<nodelist.size(); ++j)
                alist.add(-1.0);
            subgraph.add(alist);
        }

        // Finally use the collected nodes (of the connected subgraph) to get the adjacency list form of the subgraph 
        for(int i = 0; i<nodelist.size(); i++){
            for(int j=0; j<nodelist.size();j++){
                subgraph.get(i).set(j, full_graph.get(nodelist.get(i)).get(nodelist.get(j)));
                subgraph.get(j).set(i, full_graph.get(nodelist.get(i)).get(nodelist.get(j)));
            }
        }
        
        Graph component = new Graph(subgraph, sub_nodeindices, sub_nodenames);
            
        return component;
    }
}
