import org.junit.jupiter.api.Test;
import java.util.ArrayList;
import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

class MSTTest {
    
    @Test
    void primMST() {
        MST_list m = new MST_list();
        // double[][] a =
        //         {{-1, 2,-1, 6,-1},
        //          { 2,-1, 3, 8, 5},
        //          {-1, 3,-1,-1, 7},
        //          { 6, 8,-1,-1, 9},
        //          {-1, 5, 7, 9,-1},
        // };
        ArrayList<ArrayList<Double>> a = new ArrayList<ArrayList<Double>>(5);
        a.add(new ArrayList<Double>(Arrays.asList(-1.0, 2.0, -1.0, 6.0, -1.0)));
        a.add(new ArrayList<Double>(Arrays.asList(2.0, -1.0, 3.0, 8.0, 5.0)));
        a.add(new ArrayList<Double>(Arrays.asList(-1.0, 3.0, -1.0, -1.0, 7.0)));
        a.add(new ArrayList<Double>(Arrays.asList(6.0, 8.0, -1.0, -1.0, 9.0)));
        a.add(new ArrayList<Double>(Arrays.asList(-1.0, 5.0, 7.0, 9.0, -1.0)));
        int[] b = {-1, 0, 1, 0, 1};
        assertArrayEquals(b, m.primMST(a));
    }
}