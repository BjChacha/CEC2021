package etmo.util.comparators;

import etmo.core.Solution;
import java.util.Comparator;

public class PdfComparator implements Comparator {

    @Override
    public int compare(Object o1, Object o2) {
        if (o1 == null)
            return 1;
        else if (o2 == null)
            return -1;

        double distance1 = ((Solution) o1).getPdf();
        double distance2 = ((Solution) o2).getPdf();
        if (distance1 > distance2)
            return -1;

        if (distance1 < distance2)
            return 1;

        return 0;
    }
}
