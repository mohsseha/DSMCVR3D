package edu.mit.husain.DSMCVR3D;

import java.io.Serializable;
import static edu.mit.husain.DSMCVR3D.myUtils.*;

public class Result implements Serializable{
	/**
	 * Here are a bunch of important fields:
	 */
	private static final long serialVersionUID = 7L;
	public final Parameters params;
	public final Cell[][][][] trans;
	public final Cell[][][] SS;
	
	/**get the transient results at a particular time
	 * @param t	time 1st time is 1 not zero (to conform to MMA's system)
	 * @return
	 */ 
	public Cell[][][] getTrans(int t){return trans[t-1];}
	/**get a transeint cross section at a particular postion:
	 * 
	 */
	public Cell[] getTrans(int i,int j, int k){
		Cell[] res=new Cell[trans.length];
		for(int t=0;t<trans.length;t++){
			res[t]=trans[t][i-1][j-1][k-1];
		}
		return res;
	}
	/**get the average of transient result range
	 * 
	 */
	public Cell[][][] getTransAvg(int start,int end){
		Cell[][][] res=new Cell[params.cellx][params.celly][params.cellz];
		for(int i=0;i<params.cellx;i++)for(int j=0;j<params.celly;j++)for(int k=0;k<params.cellz;k++){
			res[i][j][k]=new Cell(null,i,j,k,params);
		}
		for(int t=start;t<end;t++){//we are using MMA method of counting postion in time shifting is done in getTrans(i)
			addSecondCellArrayToFirst(res,getTrans(t));
		}
		return res;
	}
	
	
	/**loads a bunch of files and averages them, then outputs them
	 * @param filenames of restart files
	 */	
	public static void main(String[] args) {
		Result output=loadAndSumResults(args);
		System.out.print(output);
		saveObject("output.result.avg",output);//the diff name is to allow wildcards to work without averaging previous results
	}
	/** load all the results in the array args and average them and then return the average to caller
	 * @param args	filenames
	 * @return Result object that has all of the results averaged. 
	 */
	public static Result loadAndSumResults(String[] args){
		Result output=(Result)loadObjectFromFile(args[0]);//first file
		for(int i=1;i<args.length;i++) output.add((Result)loadObjectFromFile(args[i]));
		return output;
	}
	/**adds an other result to this result. remember this is not averaging but just adding 
	 * 
	 * @param other result
	 */

	private void add(Result other) {
		//sanity check to make sure we do not average result that are from different runs:
		if(!other.params.equals(params)){
			System.err.println("Result:add:Error can't add two runs with different parameters");
			System.err.println("here is my params:\n***********"+params+"\n***********");
			System.err.println("here are the other  params:\n***********"+other.params+"\n***********");			
			System.exit(1);
		}
		for(int t=0;t<params.endTransSampling;t++){
			addSecondCellArrayToFirst(trans[t],other.trans[t]);
		}
		//do SS:
		addSecondCellArrayToFirst(SS,other.SS);
	}




	/** 
	 * Default constructor 
	 * @param calcParams reference to problem parameters
	 * @param gas 
	 */
	Result(Gas gas, Parameters calcParams){
		params=calcParams;//local Reference to parameters to make things conveinent here
		//Create data structures for stroring calculation result. 
		//first we calculate the transient cells:
		trans=new Cell[params.endTransSampling][params.cellx][params.celly][params.cellz];
		for(int t=0;t<params.endTransSampling;t++)
			for(int i=0;i<params.cellx;i++)
				for(int j=0;j<params.celly;j++)
					for(int k=0;k<params.cellz;k++){
						trans[t][i][j][k]=new Cell(null,i,j,k,params);//I dont' want cells with 100s ofk of reserved memory unsed
					}
		//now we do the same for the SS cells:
		SS=new Cell[params.cellx][params.celly][params.cellz];
		for(int i=0;i<params.cellx;i++)for(int j=0;j<params.celly;j++)for(int k=0;k<params.cellz;k++){
			SS[i][j][k]=new Cell(null,i,j,k,params);
		}
	}

	/** 
	 * @param gas in case we need to pass along details relating to the walls and wall collisions not just the cell properties. 
	 * @param time the time that the measurement was taken
	 */
	public void add(Cell[][][] mesh,int time){
		if(time<params.endTransSampling)addSecondCellArrayToFirst(trans[time],mesh);
		else if(time>params.startSSSampling)addSecondCellArrayToFirst(SS,mesh);
	}

	private void addSecondCellArrayToFirst(Cell[][][] cells, Cell[][][] mesh) {
		for(int i=0;i<params.cellx;i++)   // cells.length
			for(int j=0;j<params.celly;j++) // cells[i].length
				for(int k=0;k<params.cellz;k++){ // cells[i][j].length
					cells[i][j][k].add(mesh[i][j][k]);
				}
	}


	/** 
	 * human readable version of the results
	 */
	public String toString(){
		String result="\nTransient Results:\n";
		for(int t=0;t<trans.length;t++)result+=printMesh(trans[t]);
		result+="\nSSresult:\n";
		result+=printMesh(SS);
		return result;
	}



	/** 
	 * human readable version of the string we get
	 * @param Cell[][][] the cells we want to print
	 * 
	 */

	private  String printMesh(Cell[][][] cells){
		String result="";
		for(int i=0;i<params.cellx;i++)   // cells.length
			for(int j=0;j<params.celly;j++) // cells[i].length
				for(int k=0;k<params.cellz;k++){ // cells[i][j].length
					result+=cells[i][j][k]+"\t";
				}
		return result;
	}



	/** MMA function to return the transient ux variance of NE and VR 
	 * 
	 * @return double[2] results
	 */
	public double[] MMATransVariance(){
		double[] res=new double[2];
		int len=trans.length*params.cellx*params.celly*params.cellz;//length of data
		double[] uxSamples=new double[len],uxVRSamples=new double[len];
		System.err.print("length is="+len);
		int s=0;
		for(int time=0;time<trans.length;time++){
			for(int i=0;i<params.cellx;i++)   
				for(int j=0;j<params.celly;j++)
					for(int k=0;k<params.cellz;k++){ 
						Cell acell=trans[time][i][j][k];
						uxSamples[s]=acell.uxne();
						uxVRSamples[s]=acell.ux();
						s++;
					}
		}
		res[0]=variance(uxSamples);
		res[1]=variance(uxVRSamples);
		return res;
	}
}



