package etmo.metaheuristics.maoeac;


import etmo.core.Solution;
import etmo.util.comparators.RankComparator;


import java.util.Comparator;

public class SumValueComparator implements Comparator {


    /**
     * stores a comparator for check the rank of solutions
     */
    private static final Comparator comparator = new RankComparator();

    /**
     * Compares two solutions.
     *
     * @param o1
     *            Object representing the first <code>Solution</code>.
     * @param o2
     *            Object representing the second <code>Solution</code>.
     * @return -1, or 0, or 1 if o1 is less than, equal, or greater than o2,
     *         respectively.
     */
    public int compare(Object o1, Object o2) {
        if (o1 == null)
            return 1;
        else if (o2 == null)
            return -1;

        int flagComparatorRank = comparator.compare(o1, o2);
        if (flagComparatorRank != 0)
            return flagComparatorRank;

        double fitness1 = ((Solution) o1).getSumValue();
        double fitness2 = ((Solution) o2).getSumValue();
        if (fitness1 < fitness2) {
            return -1;
        }

        if (fitness1 > fitness2) {
            return 1;
        }

        return 0;
    } // compare

}
