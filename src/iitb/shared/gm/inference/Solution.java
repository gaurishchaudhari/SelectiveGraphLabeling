package iitb.shared.gm.inference;

public class Solution extends Labeling implements Comparable<Solution> {
    public Solution next;
    public int[] indexInfo;
    
    public Solution() {
    }
    public Solution(int n) {
        super(n);
    }
    public Solution(Solution s) {
        copy(s);
    }
    public Solution(int[] labelings, double score) {
        super(labelings,score);
    }
    public void copy(Solution s) {
        if (labeling==null) labeling = new int[s.labeling.length];
        for(int i=0;i < s.labeling.length;i++)
            labeling[i] = s.labeling[i];
        score = s.score;
    }
    public String indexHash() {
        StringBuffer s = new StringBuffer("");
        for(int i=0;i < indexInfo.length;i++) {
            s.append(indexInfo[i] + " ");
        }
        return s.toString();
    }
    public int compareTo(Solution s1) {
        double diff = score - s1.score;
        if(diff < -0.000001d)
            return 1;
        if(diff > 0.000001d)
            return -1;
        assert(labeling.length == s1.labeling.length);
        for (int i = 0; i < labeling.length; i++)
            if (labeling[i] > s1.labeling[i])
                return -1;
            else if(labeling[i] < s1.labeling[i])
                return 1;
        return 0;
    }
    public boolean equals(Object o) {
        return (compareTo((Solution)o) == 0);
    }
    
    public int getLabel(int var) 
    {
    	return labeling[var];
    }
    
    public void copyLabels(Solution sol) {
        for (int i = 0; i < labeling.length; i++)
            if (labeling[i] != -1) sol.labeling[i] = labeling[i];
    }
    public String toString(){
        String str = super.toString();
        if (next != null) str+= "\n"+next.toString();
        return str;
    }
   
}
