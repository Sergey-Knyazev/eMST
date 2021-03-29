// A Java program for Prim's Minimum Spanning Tree (MST) algorithm.
// The program is for adjacency matrix representation of the graph

import java.lang.*;
import java.util.ArrayList;

import java.util.HashMap;


class MST_matrix{
    // Number of vertices in the graph
    private int V=0;

    // A utility function to find the vertex with minimum key
    // value, from the set of vertices not yet included in MST
    private int minKey(double key[], Boolean mstSet[])
    {
        // Initialize min value
        double min = Integer.MAX_VALUE;
        int min_index=-1;

        for (int v = 0; v < V; v++)
            if (!mstSet[v] && key[v] < min)
            {
                min = key[v];
                min_index = v;
            }

        return min_index;
    }
    // Function to construct and print MST for a graph represented
    //  using adjacency matrix representation
    int[] primMST(double graph[][])
    {
        V = graph.length;
        // Array to store constructed MST
        int parent[] = new int[V];

        // Key values used to pick minimum weight edge in cut
        double key[] = new double [V];

        // To represent set of vertices not yet included in MST
        Boolean mstSet[] = new Boolean[V];

        // Initialize all keys as INFINITE
        for (int i = 0; i < V; i++)
        {
            key[i] = Double.MAX_VALUE;
            mstSet[i] = false;
        }

        // Always include first 1st vertex in MST.
        key[0] = 0.0;     // Make key 0 so that this vertex is
        // picked as first vertex
        parent[0] = -1; // First node is always root of MST

        // The MST will have V vertices
        for (int count = 0; count < V-1; count++)
        {
            // Pick thd minimum key vertex from the set of vertices
            // not yet included in MST
            int u = minKey(key, mstSet);

            // Add the picked vertex to the MST Set
            mstSet[u] = true;

            // Update key value and parent index of the adjacent
            // vertices of the picked vertex. Consider only those
            // vertices which are not yet included in MST
            for (int v = 0; v < V; v++)

                // graph[u][v] is non zero only for adjacent vertices of m
                // mstSet[v] is false for vertices not yet included in MST
                // Update the key only if graph[u][v] is smaller than key[v]
                if (graph[u][v]>=0 && !mstSet[v] &&
                        graph[u][v] <  key[v])
                {
                    parent[v]  = u;
                    key[v] = graph[u][v];
                }
        }
        return parent;
    }
}

class MST_list{
    // Number of vertices in the graph
    private int V=0;

    // A utility function to find the vertex with minimum key
    // value, from the set of vertices not yet included in MST
    // can be improved using priority queue
    private int minKey(double key[], Boolean mstSet[])
    {
        // Initialize min value
        double min = Double.MAX_VALUE;
        int min_index=-1;

        for (int v = 0; v < V; v++)
            if (!mstSet[v] && key[v] < min)
            {
                min = key[v];
                min_index = v;
            }

        return min_index;
    }

    // Function to construct and print MST for a graph represented
    //  using adjacency matrix representation is converted to array list representation
    int[] primMST(ArrayList<ArrayList<Double>> graph)
    {
        //size of graph is the total number of nodes in the graph
        V = graph.size();
        // Array to store constructed MST
        int parent[] = new int[V];

        // Key values used to pick minimum weight edge in cut
        double key[] = new double [V];

        // To represent set of vertices not yet included in MST
        Boolean mstSet[] = new Boolean[V];

        // Initialize all keys as INFINITE
        for (int i = 0; i < V; i++)
        {
            key[i] = Double.MAX_VALUE;//sort of infinite value.This means the max value of integer in java
            mstSet[i] = false;
        }

        // Always include first 1st vertex in MST.
        key[0] = 0.0;     // Make key 0 so that this vertex is picked as first vertex
        parent[0] = -1; // First node is always root of MST

        // The MST will have V vertices and 1 vertex is already present in parent so we need to run only V-1 times
        for (int count = 0; count < V-1; count++)
        {
            // Pick thd minimum key vertex from the set of vertices
            // not yet included in MST
            int u = minKey(key, mstSet);

            // Add the picked vertex to the MST Set
            mstSet[u] = true;

            // Update key value and parent index of the adjacent
            // vertices of the picked vertex. Consider only those
            // vertices which are not yet included in MST
            for (int v = 0; v < V; v++)

                // graph[u][v] is non zero only for adjacent vertices of m
                // mstSet[v] is false for vertices not yet included in MST
                // Update the key only if graph[u][v] is smaller than key[v]
                // the following statement effectively gets the vertex in mstSet which is closest to u
                if (graph.get(u).get(v)>=0 && !mstSet[v] && graph.get(u).get(v) <  key[v]){
                    parent[v]  = u;
                    key[v] = graph.get(u).get(v);
                }
        }
        return parent;
    }
}

class MST_fasta{
    // Number of vertices in the graph
    private int V=0;

    // A utility function to find the vertex with minimum key
    // value, from the set of vertices not yet included in MST
    // can be improved using priority queue
    private int minKey(double key[], Boolean mstSet[])
    {
        // Initialize min value
        double min = Double.MAX_VALUE;
        int min_index=-1;

        for (int v = 0; v < V; v++)
            if (!mstSet[v] && key[v] < min)
            {
                min = key[v];
                min_index = v;
            }

        return min_index;
    }

    // Function to construct and print MST for a graph represented
    //  using adjacency matrix representation is converted to array list representation
    int[] primMST(HashMap<Integer, ArrayList<Integer>> diff_consensus, HashMap<Integer, Seq> node_sequences)
    {
        //size of graph is the total number of nodes in the graph
        V = node_sequences.size();
        // Array to store constructed MST
        int parent[] = new int[V];

        // Key values used to pick minimum weight edge in cut
        double key[] = new double [V];

        // To represent set of vertices not yet included in MST
        Boolean mstSet[] = new Boolean[V];

        // Initialize all keys as INFINITE
        for (int i = 0; i < V; i++)
        {
            key[i] = Double.MAX_VALUE;//sort of infinite value.This means the max value of integer in java
            mstSet[i] = false;
        }

        // Always include first 1st vertex in MST.
        key[0] = 0.0;     // Make key 0 so that this vertex is picked as first vertex
        parent[0] = -1; // First node is always root of MST

        // The MST will have V vertices and 1 vertex is already present in parent so we need to run only V-1 times
        for (int count = 0; count < V-1; count++)
        {
            // Pick thd minimum key vertex from the set of vertices
            // not yet included in MST
            int u = minKey(key, mstSet);

            // Add the picked vertex to the MST Set
            mstSet[u] = true;

            // Update key value and parent index of the adjacent
            // vertices of the picked vertex. Consider only those
            // vertices which are not yet included in MST
            double dist = 0;
            for (int v = 0; v < V; v++){

                // graph[u][v] is non zero only for adjacent vertices of m
                // mstSet[v] is false for vertices not yet included in MST
                // Update the key only if graph[u][v] is smaller than key[v]
                // the following statement effectively gets the vertex in mstSet which is closest to u
                
                // calculate the hamming distance using arrays having difference from consensus
                dist = distance_hamming_using_consensus(u, v, diff_consensus, node_sequences);
                if (dist >=0 && !mstSet[v] && dist <  key[v]){
                    parent[v]  = u;
                    key[v] = dist;
                }
            }
        }
        return parent;
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
}
// This code is contributed by Aakash Hasija and modified by Harman Singh