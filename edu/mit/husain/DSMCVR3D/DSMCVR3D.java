package edu.mit.husain.DSMCVR3D;



import java.io.*;
import java.net.InetAddress;
import java.net.UnknownHostException;
//import static edu.mit.husain.DSMCVR3D.myUtils.*;//my personal shortcuts

/** (C) 2008 All Rights Reserved 
 * 
 * 
 * @author Husain A. Al-Mohssen
 *
 */
public class DSMCVR3D extends myUtils implements Serializable{

	/**
	 * 
	 */
	private static final long serialVersionUID = 20L,uSecBetweenWrites=3*60*1000;//2mins
	private long timeLastSave=System.currentTimeMillis(),calcStartTime=timeLastSave,timeLastRestart=timeLastSave;
	final private String restartFile,resultFile;
	private final static String inputFilename="input.txt";//never changes
	public int currentStep,ensemble;
	Parameters params;
	public Gas gas;
	public Result result;//this contains all of our calculation results


	/**
	 * @param args is either empty => start a calculation or a filename which will restart the calculation from the file name
	 */
	public static void main(String[] args) {
		DSMCVR3D calc;//this calculation 
		if(args.length>0){//need to load the reload file:
			calc=loadFromFile(args[0]);
		}else{//load and parse the input file
			calc=new DSMCVR3D();
		}
		System.err.println("DSMCVR3D -Slow KDTree- Ver "+serialVersionUID+"\n"+calc.params);//write the state to output 
		calc.run();
		System.out.print(calc.result);//Human readable results 
		saveResultsAndRestart(calc);
	}


	private static void saveResultsAndRestart(DSMCVR3D calc) {
		System.err.print("Writing Output to"+calc.resultFile);
		saveObject(calc.resultFile,calc.result);		
		saveObject(calc.restartFile,calc);
		if(calc.params.deleteRestart)deleteFileIfExist(calc.restartFile);	
	}


	public static DSMCVR3D loadFromFile(String file) {
		Object obj=loadObjectFromFile(file);//load the object
		try {
			if (obj instanceof DSMCVR3D)return (DSMCVR3D) obj;
			else {
				System.err.println(file+" is not correct time it is "+obj.getClass().toString());
				throw new Exception();
			}
		} catch (Exception e) {
			System.err.println("loadFromFile: not correct type of ojbect  in "+file);
			e.printStackTrace();
		}
		return null;
	}

	DSMCVR3D(){
		String basename="null";
		currentStep=0;ensemble=0;
		params=new Parameters(inputFilename);
		gas=new Gas(params);
		result=new Result(gas,params);
		try {
			basename=InetAddress.getLocalHost().getHostName()+"."+calcStartTime;
		} catch (UnknownHostException e) {
			e.printStackTrace();
		}
		restartFile=basename+".restart";
		resultFile=basename+".result";
		//wow that was easy :)
	}

	/**
	 * Start running our calculation from whatever point we are at based on our global variable state
	 */
	void run() {
		for(;ensemble<params.totalEnsembles;ensemble++){
			for(;currentStep<params.totalSteps;currentStep++){
				saveRestartAndResultToDisk();//save the current calculation state to disk so we can restart out calculation at any point. 
				//work is done here:
				gas.step();
				result.add(gas.cells, currentStep);//this will still capture the first (time=0) step because sort is always before collide
				printExecutionProgress();
				//save the state periodically
			}
			currentStep=0;//needed so we can do the next ensemble calc.....
			gas.resetInitialConditions();
		}
	}

	/* this variable is in a funny place. it is a global variable but only used by the method beneith it*/
	static long timeLastProgress=0;
	private void printExecutionProgress() {
		double usedRamInMB=(Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/1e6;
		long now=System.currentTimeMillis();
		if(now-timeLastProgress>60*1000) {//every minute 
			double varWSample=variance(this.gas.w);
			double[] MX0=(this.gas.cells[0][0][0]).MX();
			double progressratio=(double)(ensemble+currentStep/(double)params.totalSteps)/(((double)(params.totalEnsembles)));
			System.err.print("Progress: step="+currentStep+" ensemble="+ensemble+"	("+
					twoSig(100.0*progressratio)+"%)            started "+twoSig((now-calcStartTime)/60.e3)+"mins ago"+
					"             RAM Used="+twoSig(usedRamInMB)+"MB \t\t VarW="+varWSample+"\t"+"MX[0][][]="+twoSig(MX0[0])+"\tWtMX"+twoSig(MX0[1])+"\r");//tapes at end to make things pretty 
			timeLastProgress=now;
		}
	}



	private void saveRestartAndResultToDisk() {
		if((System.currentTimeMillis()-timeLastSave)<uSecBetweenWrites)return;
		timeLastSave=System.currentTimeMillis();
		saveObject(resultFile,result);
		//Since we are in the mood for cleaning things our we might as 
		//well call garbage collection (not guaranteed to do anything)
		System.gc();
		if((System.currentTimeMillis()-timeLastRestart)<48*48*uSecBetweenWrites)return;//every 4 hours  
		timeLastRestart=System.currentTimeMillis();		
		saveObject(restartFile,this);
	}

}
