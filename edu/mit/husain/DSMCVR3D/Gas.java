package edu.mit.husain.DSMCVR3D;


import java.io.Serializable;
import static edu.mit.husain.DSMCVR3D.myUtils.*;
import static java.lang.Math.*;

public class Gas implements Serializable {

	private static final long serialVersionUID = 16L;
	// a bunch of static variables defining the main particle properties
	final public double[] x, y, z, Cx, Cy, Cz, w;
	final Parameters params;
	final Cell[][][] cells;// mesh of current cells and parameters
	//final double SigmaW0; //find this once and forget about it. 
	final double SigmaW0; //find this once and forget about it. 
	double SigmaWLastStep,SigmaWinTop=0,SigmaWinBottom=0;//1st variableneeded for COM and last 2 are needed for no-flux BC. 
	final ArrayListInt topPks,bottomPks;//particles that hit the top and bottom BC in this Time step
	
	
	public Gas(Parameters pars) {
		params = pars;
		final int cellx = pars.cellx, celly = pars.celly, cellz = pars.cellz;
		int numPk = pars.cellx * pars.celly * pars.cellz * pars.pkPerCell;// total pks' needed 
		SigmaWLastStep=SigmaW0=(double)numPk;//set point for number of equlibrium particles in the simulation//last step is the same because we start everything @ equlibrum
		x = new double[numPk];
		y = new double[numPk];
		z = new double[numPk];
		Cx = new double[numPk];
		Cy = new double[numPk];
		Cz = new double[numPk];
		w = new double[numPk]; 
		cells = new Cell[cellx][celly][cellz];
		for (int i = 0; i < params.cellx; i++)
			for (int j = 0; j < params.celly; j++)
				for (int k = 0; k < params.cellz; k++) {
					cells[i][j][k] = new Cell(this,i,j,k,params);
				}
		resetInitialConditions();
		topPks=new ArrayListInt(pars.pkPerCell);
		bottomPks=new ArrayListInt(pars.pkPerCell);
		
		//make initial weight 
		System.err.println("debug:Gas():Reminder we only fix initial conditions for domain wide not cell wide");//because sorting now is a pain due to my small brain;//TODO: delme jsut for debugging
		double wSum=0;
		for(int i=0;i<numPk;i++)wSum+=w[i];
		final double a=((double)numPk)/wSum;//being paranoid with casting and ('s
		for(int i=0;i<numPk;i++)w[i]*=a;
	}
	/** this method resets the initial conditions based on our configuration
	 */
	public void resetInitialConditions() {
		// uniformly distribute the location of all particles:
		for (int i = 0; i < x.length; i++) {
			x[i] = params.Lx * nextDouble();
			y[i] = params.Ly * nextDouble();
			z[i] = params.Lz * nextDouble();
		}
		
		//TODO: DELETE ALL OF THIS STUFF
		//			setupHomolleIC();
		if (params.homolleIC)die(13323,"pelase kill me please .... ahhhwaaahhh");
		
		equilIC();
		//has to be done so we can measure initial step (ie before 1st collision)
		sort();
	}

	private void equilIC() {
		final double pre = params.mpvRef / sqrt(2.);
		for (int i = 0; i < x.length; i++) {
			w[i] = 1.0;
			Cx[i] = pre * nextGaussian();
			Cy[i] = pre * nextGaussian();
			Cz[i] = pre * nextGaussian();
		}
	}

	/**
	 * Method used to setup the homole relaxation problem. This is the problem
	 * we discussed in our RGD26 Paper <i>Al-Mohssen, H. A., Hadjiconstantinou,
	 * N.G.; Yet another variance reduction method for direct Monte Carlo
	 * simulations of low-signal flows, 26th International Symposium on Rarefied
	 * Gas Dynamics, July 21-25, 2008.</i>
	 */
//	private void setupHomolleIC() {
//		double vw;
//		// this is equibilent to the code in my original DSMC test code without
//		// all the cosines and sines.
//		final double pre = params.mpvRef / sqrt(2.);
//		for (int i = 0; i < x.length; i++) {
//			Cx[i] = pre * nextGaussian();
//			Cy[i] = pre * nextGaussian();
//			Cz[i] = pre * nextGaussian();
//			// remember we are dowing two Gaussians which are offset by the wall
//			// velocity
//			if (nextBoolean())
//				vw = params.wallSpeed;
//			else
//				vw = -params.wallSpeed;
//			Cx[i] += vw;
//			double feqval=nfmb(Cx[i], Cy[i], Cz[i],0,0,0,params.n0Ref,params.mass,params.TRef);
//			w[i] = feqval/ fNEHomolle(Cx[i], Cy[i], Cz[i]);
//			
//		}
//
//	}





	public void step() {		
		if(!params.homolleIC)move();
		collide();//sort is called insider collide
	}
	static State cellState, refState;
	private void collide() {
		//never collide before sorting:
		sort();
		if(cellState==null)cellState=new State();if(refState==null)refState=new State();
		refState.makeRef(params);//global reference state
		for (int i = 0; i < params.cellx; i++)
			for (int j = 0; j < params.celly; j++)
				for (int k = 0; k < params.cellz; k++) {final Cell cell=cells[i][j][k];
					if(params.regularDSMC){//what to do when we are pure DSMC:
						cell.collideUsingKDE(); 
					}else{//VRDSMC:
						//Find initial weight:
						final double WCellInit=summUpWeights(cell.pkPtrs);
						
						//Actual cell collision
						cellState.VRmeasureAndAssign(cell);cell.adjRefPDF(refState,cellState,this);//we need to use dynamically adjusted reference state when colliding 
						cell.collideUsingKDE(); 
						cell.adjRefPDF(cellState,refState,this);//Move things back 
						
						//Conserve Wt's for this cell:
						makeWtSumInArrayList(cell.pkPtrs,WCellInit);
						
					}
				}
	}
	/** multiply Wt's by a scalar 
	 * 
	 * @param pkList points to which weights. 
	 * @param Wbar value to scale by
	 */
	final private void scaleWtsInArrayList(ArrayListInt pkList,
			final double Wbar) {
			for(int i=0;i<pkList.length;i++){
				final int pkj=pkList.get(i);//pkj
				w[pkj]*=Wbar;
			}
	}
	final private double summUpWeights(ArrayListInt ptrs){
		double res=0.000000;
		for(int i=0;i<ptrs.length;i++){
			final int pkj=ptrs.get(i);//pkj
			res+=w[pkj];
		}
		return res;
	}
	
	final private void makeWtSumInArrayList(ArrayListInt ptrs,final double targetWtSum){
		final double currWtSum=summUpWeights(ptrs);
		scaleWtsInArrayList(ptrs,targetWtSum/currWtSum);
	}
	
	
	
	/** sort particles into cells. This function is responsible for collecting statistics that are also stored in cells
	 * 
	 */
	private void sort() {
		final double dx=cells[0][0][0].dx(),dy=cells[0][0][0].dy(),dz=cells[0][0][0].dz();
		final double a=SigmaW0/SigmaWLastStep;//THis is more to prevent a very small random walk that is caused by machine precision 
		if(abs(a-1.0)>1e-12&&!params.regularDSMC)System.err.println("Debug:Gas.sort():excessive global weight adjustment:\ta="+a);//TODO: DEL THIS LINE THIS IS ONLY FOR DEBUGGING
		SigmaWLastStep=0;//to do 
		//first Let's zero out all the pointer arrays:
		for (int i = 0; i < params.cellx; i++)
			for (int j = 0; j < params.celly; j++)
				for (int k = 0; k < params.cellz; k++) {
					cells[i][j][k].reset();
					cells[i][j][k].incSteps();
				}
		for(int ic=0;ic<x.length;ic++){
			SigmaWLastStep+=(w[ic]*=a);//2 things: 1. scale values of w to push total weight back to set point and 2. sum value for next step.  
			int i=(int)((x[ic])/dx);
			int j=(int)((y[ic])/dy);
			int k=(int)((z[ic])/dz);
			if(Double.isNaN(w[ic])){//TODO: debug del this block
				System.err.print("dude we are NAn!!!");
			}
			cells[i][j][k].addpk(Cx[ic],Cy[ic],Cz[ic],w[ic],ic);//adds particle to the cell and updates all the cell's property 
			if(params.regularDSMC)w[ic]=1.0;
		}
		
	}


	final private void move() {
		if(!params.regularDSMC)resetFluxCounters();
		for (int pk = 0; pk < x.length; pk++) {
			step(pk, params.dt);
		}
		//Adjust the flux for the boundary condition: this ensure conservation of mass at the edge. 
		if(!params.regularDSMC){
			if(topPks.length>0)makeWtSumInArrayList(topPks,SigmaWinTop);
			if(bottomPks.length>0)makeWtSumInArrayList(bottomPks,SigmaWinBottom);
		}
		
	}

	
	
	
	final private void resetFluxCounters() {
		SigmaWinTop=SigmaWinBottom=0.0;
		topPks.reset();
		bottomPks.reset();
	}
	/**
	 * have particle pk move dt. This implementation is for a 1D couette flow.
	 * 
	 * @param pk
	 *            particle pointer
	 * @param dt
	 *            time to integrate it. needed to do things recursively.
	 */
	void step(int pk, double dt) { // take a particle and integrate it until u
		// are done with dt!
		final double preTop = params.mpvTop/sqrt(2.),preBottom = params.mpvBottom/sqrt(2.);
		double tw;//time until we hit wall
		final double Ly = params.Ly, mass=params.mass;
		final double mpvTop = params.mpvTop,mpvBottom = params.mpvBottom;
		double cx = Cx[pk], cy = Cy[pk], cz = Cz[pk], y0 = y[pk];
		// I'm ignoring the x0 and z0 axes since we are only doing a 1D problem
		final double vw = params.wallSpeed;
		// which wall will we hit?
		if (cy > 0) {// we hit top wall potentially
			tw = (Ly - y0) / cy;
			if (dt > tw) {// we hit the top wall
				y0 = Ly - 1.e-12; // we set ourselfs sligly below top wall
				// reset velocities:
				/* **************************************************************
				 * *********TOP WALL INTERACTION HERE!!!***********
				 **************************************************************/
				cx = preTop * nextGaussian() + vw;// wall is moving
				cz = preTop * nextGaussian();
				cy = -mpvTop * nextBiased(); // pointing down
	
				SigmaWinTop+=w[pk];
				//w=\alpha (g_eq,i/g_ne,i)* Sqrt(T_NE/T_eq); alpha will be added when we are finished with the step. 
				w[pk]=(nfmb(cx,cy,cz,0,0,0,1.0,mass,params.TRef)/nfmb(cx,cy,cz,vw,0,0,1.0,mass,params.TWallTop))*sqrt(params.TWallTop/params.TRef);
				//w[pk]*=(nfmb(cx,cy,cz,0,0,0,1.0,mass,params.TRef)/nfmb(cx,cy,cz,vw,0,0,1.0,mass,params.TWallTop))*sqrt(params.TWallTop/params.TRef);
				if(!params.regularDSMC) topPks.add(pk);

				y[pk] = y0;
				Cx[pk] = cx;
				Cy[pk] = cy;
				Cz[pk] = cz;
				step(pk, dt - tw);
				return;
			}
		} else {// we potentially hit bottom wall
			tw = -y0 / cy;// positive time
			if (dt > tw) {// we hit the bottom wall
				y0 = 1.e-12;
				// reset velocities:
				// *****************************************************
				// ********** LOWER WALL INTERACTION HERE!!!***********
				// ***************************************************
				// **************************** bottom wall not moving!
				cx = preBottom * nextGaussian() - vw;// wall is moving
				cz = preBottom * nextGaussian();
				cy = mpvBottom * nextBiased(); // pointing down
				
				
				SigmaWinBottom+=w[pk];
				//w=\alpha (g_eq,i/g_ne,i)* Sqrt(T_NE/T_eq); alpha will be added when we are finished with the step. 
//				w[pk]*=(nfmb(cx,cy,cz,0,0,0,1.0,mass,params.TRef)/nfmb(cx,cy,cz,-vw,0,0,1.0,mass,params.TWallBottom))*sqrt(params.TWallBottom/params.TRef);
				w[pk]=(nfmb(cx,cy,cz,0,0,0,1.0,mass,params.TRef)/nfmb(cx,cy,cz,-vw,0,0,1.0,mass,params.TWallBottom))*sqrt(params.TWallBottom/params.TRef);
				if(!params.regularDSMC)bottomPks.add(pk);
				
				y[pk] = y0;
				Cx[pk] = cx;
				Cy[pk] = cy;
				Cz[pk] = cz;
				step(pk, dt - tw);
				return;
			}
		}

		y0 += cy * dt;
		y[pk] = y0;
		Cx[pk] = cx;
		Cy[pk] = cy;
		Cz[pk] = cz;

		return;
	}

}
