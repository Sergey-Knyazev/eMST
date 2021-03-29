import java.util.HashMap;
import java.util.Map;

public class Seq {
    final static int A=0,C=1,G=2,T=3,N=4;
    static Map<Character, Integer> nucl;
    static {
        nucl = new HashMap<Character, Integer>();
        nucl.put('A', A);
        nucl.put('C', C);
        nucl.put('G', G);
        nucl.put('T', T);
        nucl.put('N', N);
    }
    static Map<Integer, Character> inv_nucl;
    static {
        inv_nucl = new HashMap<Integer, Character>();
        inv_nucl.put(A, 'A');
        inv_nucl.put(C, 'C');
        inv_nucl.put(G, 'G');
        inv_nucl.put(T, 'T');
        inv_nucl.put(N, 'N');
    }

    private String name;
    private String seq;
    private int[] seq_enc;
    Seq(String name, String seq) {
        this.name = name;
        this.seq = seq.toUpperCase();
        seq_enc = new int[seq.length()];
        for(int i=0;i<seq.length();++i) {
            char c = seq.charAt(i);
            if(!nucl.containsKey(c)) seq_enc[i] = nucl.get('N');
            else seq_enc[i] = nucl.get(c);
        }
    }
    public String getName() {
        return name;
    }

    public String getSeq() {
        return seq;
    }

    public int[] getSeq_enc() {
        return seq_enc;
    }

    public int get_seq_len(){
        return seq_enc.length;
    }

    public void set_seq_enc(int[] seq_enc){
        this.seq_enc = seq_enc;
    }

    public void set_sequence(){
        this.seq = "";
        for(int i=0;i<get_seq_len();++i) {
            int base_num = this.seq_enc[i];
            this.seq = this.seq + inv_nucl.get(base_num);
        }
    }
}