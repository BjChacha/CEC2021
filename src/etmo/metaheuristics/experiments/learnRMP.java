package etmo.metaheuristics.experiments;

import etmo.core.SolutionSet;
import etmo.core.Variable;
import etmo.util.JMException;
import etmo.util.PseudoRandom;

import java.util.Random;

/*
 * 
*/
public class learnRMP {
	SolutionSet[] subPops_;
	int[] numVars_;
	//int maxVar;
	
	public learnRMP(SolutionSet[] subPops, int[] numVars) {
		this.subPops_ = subPops;
		this.numVars_ = numVars;
		/*
		 * maxVar = numVars_[0]; for(int i=1;i<numVars_.length;i++) { if(maxVar <
		 * numVars_[i]) { maxVar = numVars_[i]; } }
		 */
	}
	
	public double[][] learning() throws JMException{
		int numTasks = subPops_.length;
		double[][] rmpMatrix = new double[numTasks][numTasks];
		//Initialize the rmp matrix
		for(int i=0;i<numTasks;i++) {
			for(int j=0;j<numTasks;j++) {
				if(i == j) {
					rmpMatrix[i][j] = 1.0; 
				}else {
					rmpMatrix[i][j] = 0.0; 
				}
			}	
		}
		
		//add noise and build probabilistic models
		int[] numSamples = new int[numTasks];
		double[][] means = new double[numTasks][];
		double[][] stdevs = new double[numTasks][];
		double[][][] data = new double[numTasks][][];
		for(int i=0;i<numTasks;i++) {
			numSamples[i] = subPops_[i].size();
			int numRandSamples = (int)Math.floor(0.1*numSamples[i]);
			data[i] = new double[numSamples[i]][numVars_[i]];
			double[][] data_noise = new double[numSamples[i]+numRandSamples][numVars_[i]];
			for(int p=0;p<numSamples[i];p++) {
				Variable[] decisionVariables = subPops_[i].get(p).getDecisionVariables();
				for(int d=0;d<numVars_[i];d++) {
					data_noise[p][d] = decisionVariables[d].getValue();
					data[i][p][d] = decisionVariables[d].getValue();
				}
			}
			//Add noise into the data
			for(int p=0;p<numRandSamples;p++) {
				for(int d=0;d<numVars_[i];d++) {
					data_noise[p+numSamples[i]][d] = PseudoRandom.randDouble();
				}
			}
			
			//Get the mean and standard deviation for univariate distribution
			means[i] = new double[numVars_[i]];
			stdevs[i] = new double[numVars_[i]];
			int numbs = numSamples[i] + numRandSamples; 
			for(int d=0;d<numVars_[i];d++) {
				means[i][d] = 0.0;
				for(int p=0;p<numbs;p++) {
					means[i][d] += data_noise[p][d]; 
				}
				means[i][d] = means[i][d]/numbs;
			}
			
			for(int d=0;d<numVars_[i];d++) {
				stdevs[i][d] = 0.0;
				for(int p=0;p<numbs;p++) {
					stdevs[i][d] += (data_noise[p][d]-means[i][d])*(data_noise[p][d]-means[i][d]); 
				}
				stdevs[i][d] = Math.sqrt(stdevs[i][d]/(numbs-1));
			}
		}
		
		//learning the rmp matrix
		for(int i=0;i<numTasks;i++) {
			for(int j=i+1;j<numTasks;j++) {
				int dim = Math.min(numVars_[i], numVars_[j]);
				double[][][] probMatrixs = new double[2][][];
				probMatrixs[0] = new double [numSamples[i]][2];
				probMatrixs[1] = new double [numSamples[j]][2];
				for(int p=0;p<numSamples[i];p++) {
					probMatrixs[0][p][0] = 1;
					probMatrixs[0][p][1] = 1;
					for(int d=0;d<dim;d++) {
						probMatrixs[0][p][0] = probMatrixs[0][p][0]*normpdf(data[i][p][d], means[i][d], stdevs[i][d]);
						probMatrixs[0][p][1] = probMatrixs[0][p][1]*normpdf(data[i][p][d], means[j][d], stdevs[j][d]);
					}
				}
				for(int p=0;p<numSamples[j];p++) {
					probMatrixs[1][p][0] = 1;
					probMatrixs[1][p][1] = 1;
					for(int d=0;d<dim;d++) {
						probMatrixs[1][p][0] = probMatrixs[1][p][0]*normpdf(data[j][p][d], means[i][d], stdevs[i][d]);
						probMatrixs[1][p][1] = probMatrixs[1][p][1]*normpdf(data[j][p][d], means[j][d], stdevs[j][d]);
					}
				}
			    double[] testNumbs = new double[51];
			    for(int t=0;t<=50;t++) {
			    	testNumbs[t] = t*(1.0/50);
			    }
			    double minTest = loglik(testNumbs[0], probMatrixs, numTasks);
			    for(int t=1;t<=50;t++) {
			    	double value = loglik(testNumbs[t], probMatrixs, numTasks);
			    	if(minTest > value) {
			    		minTest = value;
			    	}
			    } 
				rmpMatrix[i][j] = Math.max(0, minTest + normRand(0,0.01));
				rmpMatrix[j][i] = rmpMatrix[i][j];
			}
		}
		
		return rmpMatrix;
	}
	
	/*
	 * Normal probability density function, returns the pdf of the normal
	 * distribution with mean mu and standard deviation sigma, evaluated at the value of x
	 */
	
	public double normpdf(double x, double mu, double sigma) {
		double y = 0.0;
		y = Math.exp(-0.5*Math.pow((x-mu)/sigma, 2))/(sigma*Math.sqrt(2*Math.PI));
		return y;
	}
	
	public double normRand(double mu, double sigma) {
		Random rand = new Random();
		double rnd = rand.nextGaussian();
		rnd = sigma*rnd + mu;
		return rnd;
	}
	
	public double loglik(double rmp, double[][][] probMatrix, int numTasks) {
		double y = 0.0;
		double[][][] matrixs = new double [2][][];
		for(int i=0;i<2;i++) {
			int len = probMatrix[i].length;
			matrixs[i] = new double[len][2];
			for(int d=0;d<len;d++) {
				for(int j=0;j<2;j++) {
					matrixs[i][d][j] = probMatrix[i][d][j];
				}
			}
		}
		for(int i=0;i<2;i++) {
			int len = probMatrix[i].length;
			for(int j=0;j<2;j++) {
				if(i == j) {
					for(int d=0;d<len;d++) {
						matrixs[i][d][j] = matrixs[i][d][j]*(1-(0.5*(numTasks-1)*rmp/numTasks));
					}
				}else {
					for(int d=0;d<len;d++) {
						matrixs[i][d][j] = matrixs[i][d][j]*0.5*(numTasks-1)*rmp/numTasks;
					}
				}
				
			}
			double[] sumx = new double[len];
			for(int d=0;d<len;d++) {
				sumx[d] = 0.0;
				for(int j=0;j<2;j++) {
					sumx[d] += matrixs[i][d][j];
				}
				sumx[d] = -Math.log(sumx[d]);
				y += sumx[d];
			}
		}
		return y;
	}

}
