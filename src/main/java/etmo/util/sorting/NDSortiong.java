package etmo.util.sorting;

import etmo.core.ProblemSet;
import etmo.core.SolutionSet;
import etmo.util.Distance;
import etmo.util.PORanking;
import etmo.util.comparators.CrowdingComparator;
import etmo.util.comparators.LocationComparator;

public class NDSortiong {
    public static void sort(SolutionSet pop, ProblemSet problemSet, int taskID) {
        for (int i = 0; i < pop.size(); i++)
            pop.get(i).setLocation(Integer.MAX_VALUE);

        boolean selec[] = new boolean[problemSet.getTotalNumberOfObjs()];
        int objStart = problemSet.get(taskID).getStartObjPos();
        int objEnd = problemSet.get(taskID).getEndObjPos();
        for (int i = objStart; i <= objEnd; i++) {
            selec[i] = true;
        }
        Distance distance = new Distance();
        PORanking pr = new PORanking(pop, selec);
        int loc = 0;
        for (int i = 0; i < pr.getNumberOfSubfronts(); i++) {
            SolutionSet front = pr.getSubfront(i);
            distance.crowdingDistanceAssignment(front, problemSet.getTotalNumberOfObjs(), selec);
            front.sort(new CrowdingComparator());
            for (int j = 0; j < front.size(); j++) {
                if (loc < front.get(j).getLocation())
                    front.get(j).setLocation(loc);
                loc++;
            }
        }
        pop.sort(new LocationComparator());
    }
}
