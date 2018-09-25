import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class MSTTest {
    
    @Test
    void primMST() {
        MST m = new MST();
        double[][] a =
                {{-1, 2,-1, 6,-1},
                 { 2,-1, 3, 8, 5},
                 {-1, 3,-1,-1, 7},
                 { 6, 8,-1,-1, 9},
                 {-1, 5, 7, 9,-1},
        };
        int[] b = {-1, 0, 1, 0, 1};
        assertArrayEquals(b, m.primMST(a));
    }
}