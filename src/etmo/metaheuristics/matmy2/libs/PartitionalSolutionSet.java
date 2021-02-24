package etmo.metaheuristics.matmy2.libs;

import etmo.core.SolutionSet;
import etmo.metaheuristics.maoeac.Utils;

public class PartitionalSolutionSet {
	SolutionSet solutionSet_;
	int numberOfClusters;
	int[][] neighbor;

	public PartitionalSolutionSet(SolutionSet solutionSet_,int numberOfObjectives){
		this.solutionSet_ = solutionSet_;
		this.numberOfClusters = numberOfObjectives;
		neighbor = new int[numberOfObjectives][solutionSet_.size()/numberOfObjectives];
	}

	public SolutionSet[] partitional(){
		SolutionSet[] sols = new SolutionSet[numberOfClusters];
		int sbSize = solutionSet_.size() / numberOfClusters;
		boolean[] isAsocciated = new boolean[solutionSet_.size()];
		for(int i=0;i<numberOfClusters;i++){
			sols[i] = new SolutionSet(sbSize);
		}
		int[] permutation = new int[numberOfClusters];
		etmo.metaheuristics.maoeac.Utils.randomPermutation(permutation, numberOfClusters);
		for(int i=0;i<numberOfClusters;i++){
			int n = permutation[i];
			int currentSize = solutionSet_.size();
			for(int p=0;p<solutionSet_.size();p++){
				if(isAsocciated[p]){
					currentSize--;
				}
			}
			double[] x = new double[currentSize];
			int[] idx = new int[currentSize];
			int t = 0;
			for(int j=0;j<solutionSet_.size();j++){
				if(!isAsocciated[j]){
					x[t] = Math.acos(Math.abs(solutionSet_.get(j).getNormalizedObjective(n)/solutionSet_.get(j).getDistanceToIdealPoint()));
					idx[t] = j;
					t++;
				}
			}
			// find 'niche' nearest neighboring subproblems
			Utils.minFastSort(x, idx, currentSize, sbSize);

			System.arraycopy(idx, 0, neighbor[n], 0, sbSize);
			for(int k=0;k<neighbor[n].length;k++){
				isAsocciated[neighbor[n][k]] = true;
			}
		}
		for(int i=0;i<numberOfClusters;i++){
			for(int j=0;j<sbSize;j++){
				sols[i].add(solutionSet_.get(neighbor[i][j]));
			}
		}
		return sols;
	}

}
