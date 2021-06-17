//  LZ09_F3.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package etmo.problems.LZ09;

import etmo.core.Problem;
import etmo.core.Solution;
import etmo.core.Variable;
import etmo.util.JMException;

import java.util.Vector;

/** 
 * Class representing problem LZ09_F3 
 */
public class LZ09_F3 extends Problem {
    LZ09 LZ09_ ;

    public LZ09_F3(int numberOfObjectives, int numberOfVariables) {

			Integer ptype=21;
			Integer dtype=1;
			Integer ltype=23;

     numberOfVariables_  = numberOfVariables;
     numberOfObjectives_ = numberOfObjectives;
     numberOfConstraints_= 0;
     problemName_        = "LZ09_F3";

   	 LZ09_  = new LZ09(numberOfVariables_,
   			               numberOfObjectives_,
   			               ptype,
   			               dtype,
   			               ltype) ;

     lowerLimit_ = new double[numberOfVariables_];
     upperLimit_ = new double[numberOfVariables_];
     for (int var = 0; var < numberOfVariables_; var++){
       lowerLimit_[var] = 0.0;
       upperLimit_[var] = 1.0;
     } //for

    } // LZ09_F2
   
   /** 
    * Evaluates a solution 
    * @param solution The solution to evaluate
     * @throws JMException 
    */    
    public void evaluate(Solution solution) throws JMException {
      Variable[] gen  = solution.getDecisionVariables();
      
      Vector<Double> x = new Vector<Double>(numberOfVariables_) ;
      Vector<Double> y = new Vector<Double>(numberOfObjectives_);
          
      for (int i = 0; i < numberOfVariables_; i++) {
      	x.addElement(gen[i].getValue());
      	y.addElement(0.0) ;
      } // for
        
      LZ09_.objective(x, y) ;
      
      for (int i = 0; i < numberOfObjectives_; i++)
        solution.setObjective(i, y.get(i)); 
    } // evaluate
    @Override
    public void dynamicEvaluate(Solution solution, int currentGeneration) throws JMException {

    }
} // LZ09_F3


