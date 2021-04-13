package etmo.metaheuristics.matmy2.libs;

import etmo.core.Solution;
import etmo.core.SolutionSet;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class HierarchicalClustering1 {
    List<SolutionSet> list = new <SolutionSet>ArrayList();
    public HierarchicalClustering1(List list){
        this.list = list;
    }
    public List<SolutionSet> clusteringAnalysis(int clusteringSize){
        int size = list.size();
        double minAngle = Double.MAX_VALUE;
        double angle = 0.0;
        double[] minAngles = new double[size];//每个类与之最近类之间的角度
        int[] minIndexs = new int[size];//每个类与之最近的类
        int[] index = new int[4];//当前最相似的两个类
        index[0] = index[1] = index[2] = index[3] = -1;
        for(int i=0; i<size; i++){
            double min = Double.MAX_VALUE;
            for(int j=0;j<size;j++){
                if(i != j){
                    angle = computeAngle(list.get(i).getCentroidVector(),list.get(j).getCentroidVector());
                    //angle = computeAngle(list.get(i).getCentroid(),list.get(j).getCentroid());
                    if(min > angle){
                        min = angle;
                        index[0] = i;
                        index[1] = j;
                    }
                }
            }

            minAngles[i] =  min;
            minIndexs[i] = index[1];
            if(minAngle > min){
                minAngle = min;
                index[2] = index[0];
                index[3] = index[1];
            }

//            if (index[0] == -1 && index[1] == -1 && index[2] == -1 && index[3] == -1 && i == size - 1)
//                System.out.println("DEBUG");

        }
        while(size > clusteringSize){
            SolutionSet sols = (list.get(index[2]).union(list.get(index[3])));
            list.get(index[2]).setRemove(true);
            list.remove(index[3]);
            list.add(index[3], sols);
            /*
             * 更新把index[2]个体当作角度最近个体i的最近个体indexs[i];
             */
            for(int i=0;i<list.size();i++){
                if(minIndexs[i]==index[2] && !list.get(i).isRemove()){
                    double min = Double.MAX_VALUE;
                    int sb = -1;
                    double ss = 0.0;
                    for(int j=0;j<list.size();j++){
                        if(!list.get(j).isRemove() && i!=j){
                            //ss = computeAngle(list.get(i).getCentroid(),list.get(j).getCentroid());
                            ss = computeAngle(list.get(i).getCentroidVector(),list.get(j).getCentroidVector());
                            if(min > ss){
                                min = ss;
                                sb = j;
                            }//if
                        }//if
                    }//for
                    minAngles[i] = min;
                    minIndexs[i] = sb;
                }//if
            }//for
            /*
             * 更新把index[3]个体当作角度最近个体i的最近个体indexs[i];
             */
            for(int i=0;i<list.size();i++){
                if(minIndexs[i]==index[3] && !list.get(i).isRemove()){
                    double min = Double.MAX_VALUE;
                    int sb = -1;
                    double ss = 0.0;
                    for(int j=0;j<list.size();j++){
                        if(!list.get(j).isRemove() && i!=j){
                            //ss = computeAngle(list.get(i).getCentroid(),list.get(j).getCentroid());
                            ss = computeAngle(list.get(i).getCentroidVector(),list.get(j).getCentroidVector());
                            if(min > ss){
                                min = ss;
                                sb = j;
                            }//if
                        }//if
                    }//for
                    minAngles[i] = min;
                    minIndexs[i] = sb;
                }//if
            }//for
            /*
             * 更新当前最近角度的两个个体的index;
             */
            double sAngle = Double.MAX_VALUE;
            for(int p=0;p<list.size();p++){
                if(!list.get(p).isRemove()){
                    if(sAngle > minAngles[p]){
                        sAngle = minAngles[p];
                        index[2] = p;
                        index[3] = minIndexs[p];
                    }
                }
            }

            size--;
        }//while
		/*int sss = 0;
		for(int i=0;i<list.size();i++){
			if(!list.get(i).isRemove()){
				sss++;
			}
		}*/
		/*if(sss != 16){
			System.out.println("sss = " + sss);
			System.exit(0);
		}*/
        Iterator<SolutionSet> iterator = list.iterator();
        while(iterator.hasNext()){
            if(iterator.next().isRemove()){
                iterator.remove();
            }
        }
        return this.list;
    }

    /*
     * 求两个个体之间的角度值
     */
    public double computeAngle(Solution s1, Solution s2){
        double angle = 0.0;//所求两个向量的角度
        double distanceToidealPoint1 = s1.getDistanceToIdealPoint();//S1到理想点的距离
        double distanceToidealPoint2 = s2.getDistanceToIdealPoint();//S2到理想点的距离

        double innerProduc = 0.0; //两个向量的内积
        for(int i=0; i<s1.getNumberOfObjectives(); i++){
            innerProduc += s1.getNormalizedObjective(i) * s2.getNormalizedObjective(i);
        }
        double value = Math.abs(innerProduc/(distanceToidealPoint1*distanceToidealPoint2));
        angle = Math.acos(value);
        return angle;
    }//computeAngle


}
