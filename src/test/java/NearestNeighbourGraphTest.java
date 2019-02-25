import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class NearestNeighbourGraphTest {
    @Test
    void nearest_neighbour_graph() {
        double[][] a =
                {{-1, 2, -1, 6, -1},            //   2      3
                        {2, -1, 3, 8, 5},       //a------b-----c
                        {-1, 3, -1, -1, 7},     //  \    | \   |
                        {6, 8, -1, -1, 9},      //    \6 |8  \5|7
                        {-1, 5, 7, 9, -1},      //       d-----e
                };                              //          9
        int[] b = {-1, 0, 1, 0, 1};     //MST
        NearestNeighbourGraph c = new NearestNeighbourGraph();
        List<HashSet<Integer>> d = new ArrayList<HashSet<Integer>>();
        d.add(new HashSet<Integer>());
        for (int i = 1; i < b.length; ++i) {
            d.add(new HashSet<Integer>());
            d.get(i).add(b[i]);
            d.get(b[i]).add(i);
        }
        assertEquals(d, c.nearest_neighbour_graph(a, b, 0.3));
        double[][] e =
                {{-1, 2, -1, 6, -1},
                        {2, -1, 3, 8, 5},
                        {-1, 3, -1, -1, 5},
                        {6, 8, -1, -1, 9},
                        {-1, 5, 5, 9, -1},
                };
        d.get(2).add(4);
        d.get(4).add(2);
        assertEquals(d, c.nearest_neighbour_graph(e, b, 0.0));
    }

    @Test
    void nearest_neighbour_graph2() {
        NearestNeighbourGraph c = new NearestNeighbourGraph();
        double[][] f =
                {{-1, 1, 2,-1},
                 { 1,-1,-1, 2},
                 { 2,-1,-1, 1},
                 {-1, 2, 1,-1}};
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