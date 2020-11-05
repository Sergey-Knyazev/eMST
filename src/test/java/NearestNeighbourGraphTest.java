import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class NearestNeighbourGraphTest {
    @Test
    void nearest_neighbour_graph() {
        ArrayList<ArrayList<Double>> a = new ArrayList<ArrayList<Double>>(5);
        a.add(new ArrayList<Double>(Arrays.asList(-1.0, 2.0, -1.0, 6.0, -1.0)));
        a.add(new ArrayList<Double>(Arrays.asList(2.0, -1.0, 3.0, 8.0, 5.0)));
        a.add(new ArrayList<Double>(Arrays.asList(-1.0, 3.0, -1.0, -1.0, 7.0)));
        a.add(new ArrayList<Double>(Arrays.asList(6.0, 8.0, -1.0, -1.0, 9.0)));
        a.add(new ArrayList<Double>(Arrays.asList(-1.0, 5.0, 7.0, 9.0, -1.0)));
                // {{-1, 2, -1, 6, -1},            //   2      3
                //         {2, -1, 3, 8, 5},       //a------b-----c
                //         {-1, 3, -1, -1, 7},     //  \    | \   |
                //         {6, 8, -1, -1, 9},      //    \6 |8  \5|7
                //         {-1, 5, 7, 9, -1},      //       d-----e
                // };                              //          9
        int[] b = {-1, 0, 1, 0, 1};     //MST
        NearestNeighbourGraph_list c = new NearestNeighbourGraph_list();
        List<HashSet<Integer>> d = new ArrayList<HashSet<Integer>>();
        d.add(new HashSet<Integer>());
        for (int i = 1; i < b.length; ++i) {
            d.add(new HashSet<Integer>());
            d.get(i).add(b[i]);
            d.get(b[i]).add(i);
        }
        assertEquals(d, c.nearest_neighbour_graph(a, b, 0.3));
        // double[][] e =
        //         {{-1, 2, -1, 6, -1},
        //                 {2, -1, 3, 8, 5},
        //                 {-1, 3, -1, -1, 5},
        //                 {6, 8, -1, -1, 9},
        //                 {-1, 5, 5, 9, -1},
        //         };
        ArrayList<ArrayList<Double>> e = new ArrayList<ArrayList<Double>>(5);
        e.add(new ArrayList<Double>(Arrays.asList(-1.0, 2.0, -1.0, 6.0, -1.0)));
        e.add(new ArrayList<Double>(Arrays.asList(2.0, -1.0, 3.0, 8.0, 5.0)));
        e.add(new ArrayList<Double>(Arrays.asList(-1.0, 3.0, -1.0, -1.0, 7.0)));
        e.add(new ArrayList<Double>(Arrays.asList(6.0, 8.0, -1.0, -1.0, 9.0)));
        e.add(new ArrayList<Double>(Arrays.asList(-1.0, 5.0, 7.0, 9.0, -1.0)));

        d.get(2).add(4);
        d.get(4).add(2);
        assertEquals(d, c.nearest_neighbour_graph(e, b, 0.0));
    }

    @Test
    void nearest_neighbour_graph2() {
        NearestNeighbourGraph_list c = new NearestNeighbourGraph_list();
        // double[][] f =
        //         {{-1, 1, 2,-1},
        //          { 1,-1,-1, 2},
        //          { 2,-1,-1, 1},
        //          {-1, 2, 1,-1}};
        ArrayList<ArrayList<Double>> f = new ArrayList<ArrayList<Double>>(4);
        f.add(new ArrayList<Double>(Arrays.asList(-1.0, 1.0, 2.0, -1.0)));
        f.add(new ArrayList<Double>(Arrays.asList(1.0, -1.0, -1.0, 2.0)));
        f.add(new ArrayList<Double>(Arrays.asList(2.0, -1.0, -1.0, 1.0)));
        f.add(new ArrayList<Double>(Arrays.asList(-1.0, 2.0, 1.0, -1.0)));
        
        int[] g = {-1,0,3,1};
        List<HashSet<Integer>> h = new ArrayList<HashSet<Integer>>();
        for(int i: g) {
            h.add(new HashSet<Integer>());
        }
        for(int i=1; i<g.length; ++i) {
            h.get(i).add(g[i]);
            h.get(g[i]).add(i);
        }
        h.get(0).add(2);
        h.get(2).add(0);
        assertEquals(h, c.nearest_neighbour_graph(f, g, 0.0));
    }
}