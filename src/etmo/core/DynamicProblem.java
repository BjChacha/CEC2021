package etmo.core;

import etmo.util.JMException;

public class DynamicProblem extends Problem{
	/**
	 * Defines the change times
	 */
	protected double t_;
	/**
	 * Defines the change cycle
	 */
	protected double fc_;
	/**
	 * Defines the change ratio
	 */
	protected double sc_;
	
	/*
	 * gc indicates the generation counter
	 * fc indicates the frequency of change
	 * sc indicates the severity of change
	*/
	
	public void updateTime(int gc)
	{
		this.t_ = (1.0/sc_)*Math.floor(gc/fc_);
	}
	
	public void resetTime()
	{
		this.t_ = 0.0;
	}
	@Override
	public void evaluate(Solution solution) throws JMException {
		// TODO Auto-generated method stub

	}

	@Override
	public void dynamicEvaluate(Solution solution, int currentGeneration) throws JMException {
		// TODO Auto-generated method stub
		
	}

}
