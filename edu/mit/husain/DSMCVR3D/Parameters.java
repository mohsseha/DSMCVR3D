package edu.mit.husain.DSMCVR3D;

import java.io.*;
import java.util.Properties;
import static java.lang.Math.*;//so u can say thing like sin(x)//it seems that this class is really slow ... 
import static edu.mit.husain.DSMCVR3D.myUtils.*;//my personal shortcuts


public class Parameters implements Serializable {
	private static final long serialVersionUID = 28L;
	final int totalSteps,totalEnsembles,endTransSampling,startSSSampling,defPkVCell;
	final public int cellx,celly,cellz,pkPerCell,PKS_EPS;	
	final public double Lx, Ly,Lz,n0Ref,dt,wallSpeed,mass,dia,k,TWallTop,TWallBottom,TRef,prefactorRef,CbarRef,neff,mfpRef,MX;
	final public double mpvTop,mpvRef,mpvBottom;
	final public boolean homolleIC,regularDSMC,deleteRestart;
	final public long debugSeed;
	
	public Parameters(String inputFileName) {
		//		Load the configuration from a config file
		Properties properties = new Properties();
		try {
			properties.load(new FileInputStream(inputFileName));
		} catch (IOException e) {e.printStackTrace();}
		k=1.3806503e-23;

		//these should fail if the variable is not defined in the input file:
		totalSteps =parseInt(properties.getProperty("totalSteps"));
		totalEnsembles =parseInt(properties.getProperty("totalEnsembles"));
		endTransSampling =parseInt(properties.getProperty("endTransSampling"));
		int tmp =parseInt(properties.getProperty("startSSSampling"));
		if(tmp==-1)startSSSampling=endTransSampling+1;
		else startSSSampling=tmp;
		cellx =parseInt(properties.getProperty("cellx"));
		celly =parseInt(properties.getProperty("celly"));
		cellz =parseInt(properties.getProperty("cellz"));
		defPkVCell=parseInt(properties.getProperty("defPkVCell"));

		pkPerCell =parseInt(properties.getProperty("pkPerCell")); 
		dia =parseDouble(properties.getProperty("dia"));
		n0Ref =parseDouble(properties.getProperty("n0Ref"));

		homolleIC=parseBoolean(properties.getProperty("homolleIC"));
		regularDSMC=parseBoolean(properties.getProperty("regularDSMC"));
		//dynamicallyAdj=parseBoolean(properties.getProperty("dynamicallyAdj"))&&(!regularDSMC);//false if we are doing DSMC
		deleteRestart=parseBoolean(properties.getProperty("deleteRestart"));

		mfpRef=1.0 / (sqrt(2.0) * PI * dia * dia * n0Ref);
		Lx =mfpRef*parseDouble(properties.getProperty("NormalizedLx"));
		Ly =mfpRef*parseDouble(properties.getProperty("NormalizedLy"));
		Lz =mfpRef*parseDouble(properties.getProperty("NormalizedLz"));
		double relativedt=parseDouble(properties.getProperty("relativeDt"));
		mass =parseDouble(properties.getProperty("mass"));
		TRef =parseDouble(properties.getProperty("TRef")); 
		double dT=parseDouble(properties.getProperty("dT"));
		TWallTop=TRef+dT;
		TWallBottom=TRef-dT;

		mpvRef=sqrt(2*k*TRef/mass);
		mpvTop=sqrt(2*k*TWallTop/mass);
		mpvBottom=sqrt(2*k*TWallBottom/mass);



		CbarRef = 1.12838 * mpvRef;
		wallSpeed=mpvRef*parseDouble(properties.getProperty("NormalizedWallSpeed"));

		if(pkPerCell<20){System.err.print("dude we don't have enough particles per cell for the run to be valid!\n");System.exit(1);}

		//calculated variables:
		prefactorRef=pow(PI,-3./2.)*pow(mpvRef,-3.);
		MX=7*mpvRef;
		dt=relativedt*mfpRef/CbarRef;
		double minCellLength=min(min(Lx/cellx,Ly/celly),Lz/cellz);//smallest cell dimention
		if(dt*CbarRef>2*minCellLength)System.err.println("WARNING:Parameters:\t your steps are larger than 2*cell length... will proceed...");


		int gasSize=pkPerCell*cellx*celly*cellz;//note this is exact number of pk to be created
		neff=Lx*Ly*Lz*n0Ref/gasSize;//total # of gas molecules/#num of computational particles

		PKS_EPS=Integer.parseInt(properties.getProperty("PKS_EPS"));

		debugSeed=System.currentTimeMillis();
		//initialize random seed:
		tmp=Integer.parseInt(properties.getProperty("randomSeed"));
		if(tmp!=-1)System.err.println("debug:setting randomSeed to "+tmp);
		if(tmp==-1)r.setSeed(debugSeed);
		else r.setSeed(tmp);
	}

	//these methods exist because the properties class is retarted 
	private int parseInt(String property) {return Integer.parseInt(property.trim());}
	private boolean parseBoolean(String property) {return Boolean.parseBoolean(property.trim());}
	private double parseDouble(String property) {return Double.parseDouble(property.trim());}

	public String toString(){ 
		String result="Mode:";
		if(regularDSMC)result+="Regular DSMC\n"; 
		else {result+="VR DSMC (Stabalized)\n";
			result+="\t pk/eps="+PKS_EPS+",eps/mpv="+twoSig(eps()/mpvRef)+",eps="+twoSig(eps())+"m/s";
			result+="\tf is Dynamically adjusted reference state\n";//always the case 
		}
		result+="Steps:\t totSteps="+totalSteps+",totalEnsembles="+ totalEnsembles+"\t\t,endTransSampling="+endTransSampling+",startSSSampling="+startSSSampling+"\t pk/cell="+pkPerCell+"\n";
		if(homolleIC)result+="Running Mode: \t\t Homolle IC problem. This problem is 0D and has a Kn=0.0\n";
		else{
			result+="Kn:\t in x,y,z directions\t"+twoSig(mfpRef/Lx)+","+twoSig(mfpRef/Ly)+","+twoSig(mfpRef/Lz)+"\n";
			result+="Cells:\t x,y,z="+cellx+","+celly+","+cellz+"\tCell Knx,y,z="+twoSig((Lx/cellx)/mfpRef)+","+twoSig((Ly/celly)/mfpRef)+","+twoSig((Lz/cellz)/mfpRef)+"\n";
		}
		result+="dt:\t "+dt+"\tsec.\tdt/MTBC="+twoSig(dt/(mfpRef/CbarRef))+"\t total run time="+twoSig(totalSteps*dt/(mfpRef/CbarRef))+"MTBC\n";
		result+="BCs:\twallSpeed="+twoSig(wallSpeed)+"\tTRef="+twoSig(TRef)+",TWallTop="+twoSig(TWallTop)+",TWallBottom="+twoSig(TWallBottom)+"\t\t\tCbar="+twoSig(CbarRef);
		return result;
	}


	/**  It's tells you when to parameters are the same. 
	 * Java is a bit restarted this way ....
	 */
	public boolean equals(Object obj) {
		if(obj==null)return false;
		return this.toString().equals(obj.toString());
	}

	public int getPkPerCell() {
		return pkPerCell;
	}

	/** epsilon used in NN calculation at iteration i
	 * 
	 * @param i current iteration
	 * @return epsilon
	 */
	final private double eps(final int i) {return eps((double)i);}
	final private double eps(final double i){
		return mpvRef*pow(i/(0.26*pkPerCell),1.0/3.0);// it's easy to argue the form of this formula the 0.26 was "experimentally" determined
	}

	/** so everyone can know what the epsilon of this ell is:
	 */
	public double eps(){
		return  eps(PKS_EPS);
	}
}


