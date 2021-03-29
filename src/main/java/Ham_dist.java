import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.Observable;

import static java.lang.Math.log;

public class Ham_dist extends Observable {
    private float edgeThreshold = 10000;
    private File inputFile;
    private File outputFile;
    
    public void setEdgeThreshold(float edgeThreshold) {
        this.edgeThreshold = edgeThreshold;
    }


    public void setInputFile(File inputFile) {
        this.inputFile = inputFile;
    }

    public void setOutputFile(File outputFile) {
        this.outputFile = outputFile;
    }

    public void hammingFasta() {
        PrintWriter f = null;
        try {
            LinkedList<Seq> seqs = read_fasta(inputFile);
            // System.out.print(String.format("%d% numseq ", seqs.size()));
            double[][] dist = hamming(seqs);
            f = new PrintWriter(outputFile);
            f.println("Source,Target,Dist");
            for (int i = 0; i < dist.length; ++i) {
                for (int j = 0; j < dist.length ; ++j) {
                    if (dist[i][j] <= edgeThreshold && j>i) f.println(
                            String.format("%s,%s,%.11f", seqs.get(i).getName(), seqs.get(j).getName(), dist[i][j]));
                }
            }
        }
        catch(FileNotFoundException e) {
            e.printStackTrace();
        }
        finally {
            if(f != null) f.close();
        }
    }

    public double[][] hamming(LinkedList<Seq> seqs) {
        double[][] dist = new double[seqs.size()][seqs.size()];
        long startTime = System.nanoTime(), estimatedTime;
        int pairs_count = (dist.length * dist.length - dist.length)/2;
        int current_pair = 0;
        for (int i = 0; i < dist.length; ++i) {
            for (int j = 0; j < i; ++j) {
                current_pair++;
                dist[i][j] = dist[j][i] = hamming(seqs.get(i).getSeq_enc(), seqs.get(j).getSeq_enc());
            }
        }
        setChanged();
        notifyObservers(100);
        return dist;
    }

    public static double hamming(int[] s1, int[] s2) {
        if(s1.length != s2.length){
            System.out.println("Seq lengths not equal");
        }
        int length = s1.length;
        int dist = 0;
        for(int i=0; i<length; ++i) {
            if(s1[i] != s2[i]){
                dist = dist + 1;  
            }
        }

        return dist;

    }

    public static LinkedList<Seq> read_seqs(Scanner sc) {
        LinkedList<Seq> seqs = new LinkedList<Seq>();
        String name="", seq="";
        while(sc.hasNextLine()) {
            String line = sc.nextLine().trim();
            if(line.length() == 0) continue;
            if(line.charAt(0)=='>') {
                if (name.length()!=0) seqs.add(new Seq(name, seq));
                name = line.substring(1);
                seq="";
            }
            else seq=seq.concat(line);
        }
        if(name.length()!=0) seqs.add(new Seq(name, seq));
        return seqs;
    }

    private static LinkedList<Seq> read_fasta(File inputFile) throws FileNotFoundException {
        Scanner sc = new Scanner(inputFile);
        LinkedList<Seq> a = read_seqs(sc);
        sc.close();
        return a;
    }
}
