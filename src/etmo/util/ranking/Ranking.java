package etmo.util.ranking;

import etmo.core.SolutionSet;

public interface Ranking {
	public SolutionSet getSubfront(int layer);
	public int getNumberOfSubfronts();
}
