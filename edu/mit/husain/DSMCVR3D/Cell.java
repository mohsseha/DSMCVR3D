package edu.mit.husain.DSMCVR3D;

import static edu.mit.husain.DSMCVR3D.myUtils.*;
import static java.lang.Math.*;

import java.io.Serializable;



final public class Cell implements Serializable{

	/* and the test branch delme DELME
	 * this is the 	MASTER  branch
	 */
	private static final long serialVersionUID = 29L;
	final ArrayListInt pkPtrs;
	final Parameters pars;
	private final double x0,y0,z0,x1,y1,z1;//edges of the cell
	private double carryOver,MX,WtMX;//Number of collision carry overs from last step
	static final int MAX_HIST=30;
	final private int[] KDEHistogram=new int[MAX_HIST];//should be much smaller but JIC 
	private int stepCollisions;//number of collisions that actually happned in last time collision() was called
	//	transient static KDTree NNEngine;
	private transient static GaussianNNEngine NNEngine;
	private Gas gas;

	final static SlowPoint tmpA=new SlowPoint(0,0,0),tmpB=new SlowPoint(0,0,0);//I hate using the SlowPoint Class but if it's static and final we shoudl be fairly efficiant 
	public final Sums NE,Eq;

	/** this method should allow us to add the statistic from two cell. meant mainly to be used when aggregating results
	 * @param other the other cell*/
	final public void add(Cell other){//we don't need to aggrage all variables
		MX=max(MX,other.MX);//Propagate the worst case sernario
		WtMX=max(WtMX,other.WtMX);//
		stepCollisions+=other.stepCollisions;	
		arrayAddTo(KDEHistogram,other.KDEHistogram);
		assert pars==other.pars :"Cell:add: very basic sanity check!";
		NE.add(other.NE);
		Eq.add(other.Eq);
	}


	/** Create cell between these set of points and with this initial capacity
	 * 
	 * @param i, j ,k where the cell is going to be positioned
	 * @param pars parameters of the problem 
	 */
	public Cell(Gas g,int i, int j , int k, Parameters params, int initSize) {
		pars=params;
		gas=g;
		if(!pars.regularDSMC&&NNEngine==null)NNEngine=new GaussianNNEngine(gas,pars);//create on object that is re-used across all cells;

		NE=new Sums(pars);Eq=new Sums(pars);
		//initialize the other variables:
		carryOver=0;MX=6.0*pars.mpvRef;//completely arbitrary
		stepCollisions=0;//number of collisions that actually happned in last time collision() was called

		pkPtrs=new ArrayListInt(initSize);

		x0=i*(pars.Lx/pars.cellx);x1=(i+1)*(pars.Lx/pars.cellx);
		y0=j*(pars.Ly/pars.celly);y1=(j+1)*(pars.Ly/pars.celly);
		z0=k*(pars.Lz/pars.cellz);z1=(k+1)*(pars.Lz/pars.cellz);
		reset();
	}
	/** this constructor creates a cell with no memory used in the pkPtrs arraylist. 
	 */
	public Cell(Gas g,int i,int j,int k, Parameters params){
		this(g,i,j,k,params,0);
	}


	/** cell width
	 * @return dx cell of a cell
	 */
	final public double dx() {return x1-x0;}
	/** cell hight
	 * @return dy cell of a cell
	 */
	final public double dy() {return y1-y0;}
	/** cell depth
	 * @return dz cell of a cell
	 */
	final public double dz() {return z1-z0;}




	// this is the code that returns VR cell properties:

		
		/** returns the VR number density of the cell in units of (real) pk/m^3	 */
		public double n0(){      
			return (1.0/volume())*(pars.neff*(NE.pkSum()-Eq.pkSum() + pars.pkPerCell));}
		/** returns the regular(NE) number density of the cell in units of (real) pk/m^3	 */
		public double rhoNE(){return (1.0/volume())*(pars.neff*( NE.pkSum()));}

		/** returns VR ux */
		public double ux(){			
			double CapitalDeltau = (NE.Cx()-Eq.Cx());
			return (pars.neff/(volume()*n0()))*
			((Eq.Nsamples/Eq.Nsteps)*CapitalDeltau) ;//+ 
					//( pars.n0Ref*volume())/n0()*0.0);//u0=0 because reference state has 0 mean
		}
		/** returns Non-Equilibrium ux */
		public double uxne(){
			return (pars.neff/(volume()*rhoNE()))* 
			((Eq.Nsamples/Eq.Nsteps)*NE.Cx());
		}

		/** returns VR uy */
		public double uy(){
			double CapitalDeltau = (NE.Cy()-Eq.Cy());
			return (pars.neff/(volume()*n0()))*
			((Eq.Nsamples/Eq.Nsteps)*CapitalDeltau);//u0=0 because reference state has 0 mean
		}
		/** returns Non-Equilibrium uy */
		public double uyne(){
			return (pars.neff/(volume()*rhoNE()))*
			((Eq.Nsamples/Eq.Nsteps)*NE.Cy());
		}

		/** returns VR uz */
		public double uz(){
			double CapitalDeltau = (NE.Cz()-Eq.Cz());
			return (pars.neff/(volume()*n0()))*
			((Eq.Nsamples/Eq.Nsteps)*CapitalDeltau);//u0=0 because reference state has 0 mean
		}
		/** returns Non-Equilibrium uz */
		public double uzne(){
			return (pars.neff/(volume()*rhoNE()))*
			((Eq.Nsamples/Eq.Nsteps)*NE.Cz());
		}


		/** VR Temperature*/	
		public double T(){
			final double k=pars.k,mass=pars.mass, Nsamples=Eq.Nsamples,Nsteps = Eq.Nsteps,Neq = pars.n0Ref,NEff = pars.neff,Teq=pars.TRef;
			double CapitalDeltac2 = ((NE.Cx2() + NE.Cy2() + NE.Cz2()) 
					- (Eq.Cx2() + Eq.Cy2() + Eq.Cz2()));
			CapitalDeltac2 = (Nsamples*CapitalDeltac2)/Nsteps;
			final double n0VR=n0();
			double CapitalSigmaU2VR = (ux()*ux() + uy()*uy() + uz()*uz());
			double res=-(CapitalSigmaU2VR*mass)/(3.*k) + (Neq*Teq)/n0VR + 
			(CapitalDeltac2*mass*NEff)/(3.*k*volume()*n0VR);
			if(res<0){ res=Tne();
			}
			return res;
		}
		
		
		private void printWts() {//TODO: DELETE THIS METHOD
			for(int i=0;i<gas.w.length;i++)System.err.println(gas.w[i]);
		}


		/** NE Temperature*/
		public double Tne(){
			final double k=pars.k,mass=pars.mass, Nsamples=Eq.Nsamples,Nsteps = Eq.Nsteps,NEff = pars.neff;
			double CapitalSigmac2NE = (Nsamples/Nsteps)*(NE.Cx2() + NE.Cy2() + NE.Cz2());
			double CapitalSigmaU2NE = pow(uxne(),2.0)+pow(uyne(),2.0)+pow(uzne(),2.0);
			return (mass*(-CapitalSigmaU2NE + (CapitalSigmac2NE*NEff)/(volume()*rhoNE())))/(3.*k);
		}

		/** return VR Shear in the xy plane*/
		public double pixy(){
			final double k=pars.k,mass=pars.mass, Nsamples=Eq.Nsamples,Nsteps = Eq.Nsteps,Neq = pars.n0Ref,NEff = pars.neff,Teq=pars.TRef,vol=volume();
			double CapitalDeltaxy = NE.CxCy()-Eq.CxCy();
			return 218232.0;//TODO:this is almost certainly wrong! why is there a NE property here? ~!!!(mass*((Nsamples*CapitalDeltaxy*
//					NEff)/Nsteps
//					- vol*ux()*uy()*rhoNE()))/vol;
		}


		/** return NE Shear in the xy plane*/
		public double pixyne(){
			final double mass=pars.mass, Nsamples=Eq.Nsamples,Nsteps = Eq.Nsteps,NEff = pars.neff,vol=volume();
			double CapitalDeltaxyne = NE.CxCy();
			return 218232.5;//TODO: this is also wrong .....possibley.....(mass*((Nsamples*CapitalDeltaxyne*NEff)/Nsteps- 
//					vol*uxne()*uyne()*rhoNE()))/vol;
		}


		/** return VR heat-flux*/
		public double qy(){
			final double k=pars.k,mass=pars.mass, Nsamples=Eq.Nsamples,Nsteps = Eq.Nsteps,Neq = pars.n0Ref,NEff = pars.neff,Teq=pars.TRef,vol=volume();
			double Myx2NE = (NEff*Nsamples/Nsteps)*NE.CyCx2();
			double Myx2Eq = (NEff*Nsamples/Nsteps)*Eq.CyCx2();
			double Myy2NE = (NEff*Nsamples/Nsteps)*NE.CyCy2();
			double Myy2Eq = (NEff*Nsamples/Nsteps)*Eq.CyCy2();
			double Myz2NE = (NEff*Nsamples/Nsteps)*NE.CyCz2();
			double Myz2Eq = (NEff*Nsamples/Nsteps)*Eq.CyCz2();

			double MxyNE = (NEff*Nsamples/Nsteps)*NE.CxCy();
			double MxyEq = (NEff*Nsamples/Nsteps)*Eq.CxCy();

			double Mx2NE = (NEff*Nsamples/Nsteps)*NE.Cx2();
			double Mx2Eq = (NEff*Nsamples/Nsteps)*Eq.Cx2();
			double My2NE = (NEff*Nsamples/Nsteps)*NE.Cy2();
			double My2Eq = (NEff*Nsamples/Nsteps)*Eq.Cy2();
			double Mz2NE = (NEff*Nsamples/Nsteps)*NE.Cz2();
			double Mz2Eq = (NEff*Nsamples/Nsteps)*Eq.Cz2();

			double MCxNE = (NEff*Nsamples/Nsteps)*NE.Cx(); 
			double MCxEq = (NEff*Nsamples/Nsteps)*Eq.Cx(); 
			double MCyNE = (NEff*Nsamples/Nsteps)*NE.Cy(); 
			double MCyEq = (NEff*Nsamples/Nsteps)*Eq.Cy(); 
			double WsumEq= NEff*Eq.pkSum();

			double Mux = ux(), Muy = uy();
			double R = 2.*MCxNE*Mux*Muy - Muy*(pow(Mux,2.0) + pow(Muy,2.0)) + 
			MCyNE*(pow(Mux,2.0) + 3.*pow(Muy,2.)) - Muy*Mx2NE - 2*Mux*MxyNE - 
			3.*Muy*My2NE + Myx2NE + Myy2NE + Myz2NE - Muy*Mz2NE;
			double WR = 2*MCxEq*Mux*Muy + MCyEq*(pow(Mux,2) + 3*pow(Muy,2)) - 
			Muy*Mx2Eq - 2*Mux*MxyEq - 3*Muy*My2Eq + Myx2Eq + 
			Myy2Eq + Myz2Eq -Muy*Mz2Eq -Muy*(pow(Mux,2) + 
					pow(Muy,2))*WsumEq;
			double phieq=(mass*((-5*k*Neq*Teq*vol*uy())/mass - Neq*uy()*(pow(ux(),2) + pow(uy(),2))))/(2.*vol);
			return (mass*(R - WR))/(2.*vol) + phieq;
		}

		/** return NE heat-flux in y direction*/
		public double qyne(){
			final double k=pars.k,mass=pars.mass, Nsamples=Eq.Nsamples,Nsteps = Eq.Nsteps,Neq = pars.n0Ref,NEff = pars.neff,Teq=pars.TRef,vol=volume();
			double Myx2NE = (NEff*Nsamples/Nsteps)*NE.CyCx2();
			double Myy2NE = (NEff*Nsamples/Nsteps)*NE.CyCy2();
			double Myz2NE = (NEff*Nsamples/Nsteps)*NE.CyCz2();

			double MxyNE = (NEff*Nsamples/Nsteps)*NE.CxCy();

			double Mx2NE = (NEff*Nsamples/Nsteps)*NE.Cx2();
			double My2NE = (NEff*Nsamples/Nsteps)*NE.Cy2();
			double Mz2NE = (NEff*Nsamples/Nsteps)*NE.Cz2();

			double MCxNE = (NEff*Nsamples/Nsteps)*NE.Cx(); 
			double MCyNE = (NEff*Nsamples/Nsteps)*NE.Cy(); 

			double Mux = uxne(), Muy = uyne();
			double R = 2.*MCxNE*Mux*Muy - Muy*(pow(Mux,2.0) + pow(Muy,2.0)) + 
			MCyNE*(pow(Mux,2.0) + 3.*pow(Muy,2.)) - Muy*Mx2NE - 2*Mux*MxyNE - 
			3.*Muy*My2NE + Myx2NE + Myy2NE + Myz2NE - Muy*Mz2NE;
			return (mass*(R))/(2.*vol);
		}


	/**this should be consistent with the calculate candidates method below in addition to how Neff is calcualated . 
	 * 
	 * @return volume of the cell
	 */
	final public double volume() {
		return dz()*dy()*dx();
	}
	/** number of collisions that happened*/
	public double collisions(){return stepCollisions;}
	/** MX*/
	final public double[] MX(){double[] a={MX,WtMX};return a;}
	/** Calculate the number of collision candidates we need for a collision in this cell
	 * 
	 * @return an integer number of candidates. carryOver will move things from one call to the function to the next to tke care of the non integer arithemtic
	 */
	public int calculateCandidatesIncludingCarryOver() {
		double exactCanidate,dia=pars.dia,neff=pars.neff,volume=volume(),dt=pars.dt;
		int cellNum=pkPtrs.length;
		//God I hope this is exactly right. note that you have to be careful with the casting and such. 
		exactCanidate = ((double)(cellNum) * (double)(cellNum)  * PI *                   //careful about the FREAKING 1.0!!!!
				dia * dia * MX * neff *										// eventually SHould be:cellNum * (cellNum - 1.0)
				dt) / (2.000000 * volume);
		exactCanidate += carryOver; //carryOver takes care of the cases when the number of collisions is not an integer 
		int canidates = (int) exactCanidate;
		carryOver = exactCanidate - canidates;//what's left is used for next time we call this function

		return canidates;
	}


	final public void collideUsingKDE() {
		final double[] Cx=gas.Cx,Cy=gas.Cy,Cz=gas.Cz,w=gas.w;
		final double eps=pars.eps();
		stepCollisions = 0;

		if(!pars.regularDSMC)NNEngine.rebuild(pkPtrs);

		int candidates=calculateCandidatesIncludingCarryOver(); 
		for (int ic = 0; ic < candidates; ic++) {
			//these are the offsets inside the pkPtrs Array. If I know to need the exact offset in the global x array I do pkPtrs.get(offA)
			final int offA = nextInt(pkPtrs.length),gA=pkPtrs.get(offA);//gA and gB are the global offsets of the particles 
			final int offB = nextInt(pkPtrs.length),gB=pkPtrs.get(offB);
			final double vij = vij(gA,gB);

			//MX updates:
			final double mxwt=max(1.0,max(w[gA],w[gB]));//the 1.0 is in case we are running a regular DSMC calculation
			if(mxwt*vij>MX){
				if(mxwt>WtMX)WtMX=mxwt;
				MX=1.1*mxwt*vij;candidates=calculateCandidatesIncludingCarryOver();}//updating candidates is required otherwise this particular step would not be executed properly
			final double PNE=vij/MX;

			if (nextDouble() < PNE) {//We use NONequiblrum Collision rates				
				final double whatA=wKDE(offA,eps),whatB=wKDE(offB,eps);
				if(!pars.regularDSMC){NNEngine.cache(offA);NNEngine.cache(offB);}
				stepCollisions++;	
				//Dear Lord: I humbly ask you to make the following lines correct, Amen++++
				calculatePostC(Cx,Cy,Cz,gA,gB,vij);
				tmpA.assignVsAtPos(Cx,Cy,Cz,gA);
				tmpB.assignVsAtPos(Cx,Cy,Cz,gB);
				if(!pars.regularDSMC){NNEngine.update(offA);NNEngine.update(offB);}
				w[gA]=(w[gB]=(whatA*whatB));
			}else{//we have to adjust the weights when the step is rejected:
				//I'm being aggressive with final because I'm not sure (but hope) it might help generate slightly faster code. 
				final double PeqA=w[gB]*vij/MX;
				final double PeqB=w[gA]*vij/MX;

				//Dear Lord: and this please++:
				w[gA]*=((1.0-PeqA)/(1.0-PNE));
				w[gB]*=((1.0-PeqB)/(1.0-PNE));
			}
		}//went through all collisions in this cell
	}



	private double vij(final int aA,final int aB) {
		final double[] Cx=gas.Cx,Cy=gas.Cy,Cz=gas.Cz;
		return abs(Cx[aA] - Cx[aB],
				Cy[aA] - Cy[aB],
				Cz[aA] - Cz[aB]);
	}




	final private double wKDE(final int offj, final double eps) {
		final double[] w=gas.w;
		if(pars.regularDSMC)return 1.0;//when doing regular DSMC weights are not used., 
		final ArrayListInt nnList=NNEngine.NN(offj,eps);
		KDEHistogram[min(nnList.length,MAX_HIST-1)]++;
		return avg(w,nnList);
	}







	final private static void calculatePostC(final double[] cx,final double[] cy,final double[] cz,
			final int a, final int b,final double vr) {
		final double angx, angy, angz, vcmx, vcmy, vcmz, sinTheta, cosTheta, phi;

		vcmx = 0.5 * (cx[a] + cx[b]);
		vcmy = 0.5 * (cy[a] + cy[b]);
		vcmz = 0.5 * (cz[a] + cz[b]);

		//calculate a uniformly distributed scattering angle:
		cosTheta = 1 - 2 * nextDouble();
		sinTheta = sqrt(1 - cosTheta * cosTheta);
		phi = 2 * PI * nextDouble();
		angx =  cosTheta;
		angy =  sinTheta * cos(phi);
		angz =  sinTheta * sin(phi);


		final double aX,aY,aZ,bX,bY,bZ;
		aX = vcmx + 0.5 * vr*angx;
		aY = vcmy + 0.5 * vr*angy;
		aZ = vcmz + 0.5 * vr*angz;

		bX = vcmx - 0.5 * vr*angx;
		bY = vcmy - 0.5 * vr*angy;
		bZ = vcmz - 0.5 * vr*angz;

		tmpA.assignVs(aX,aY,aZ);
		tmpB.assignVs(bX,bY,bZ);

	}


	/**reset all the variables in this cell. The goal here is to have a cell that contains all the information calculated withint a timestep 
	 * for individual integration steps or contains the state of many averages if it is part of results..... 
	 */
	final public void reset() {
		stepCollisions=0;
		arrayAssignTo(KDEHistogram,0);
		pkPtrs.reset();
		NE.reset();Eq.reset();
	}
	/** add particle to this cell and update the cell's properties to reflect his:
	 * @param cx,cy,cz the position of pk in velocity space
	 * @param w wt of pk
	 * @param ic location of particle 
	 */
	final public void addpk(final double cx,final double cy,final double cz,final double w,final int loc) {
		if(Double.isNaN(w))die(1337,"\n\n\n\t\tDied because w is NaN; probably caused by a small number of particles in a cell;pk is"+loc);
		pkPtrs.add(loc);
		NE.addpk(cx,cy,cz,1.0);
		Eq.addpk(cx,cy,cz,w);
	}
	/** this is the function to modify to change the format of the cell result output. 
	 * currently it outputs the variance of the weights only.
	 */
	public String toString(){
		return (Eq.varW())+" ";
	}


/** 
Adjust the refernce state to the one given. States are assuming VR porerpties
 * 	this uses the funciton in utilities...
 * */
	void adjRefPDF(State from, State to, Gas gas) {
		final double mass0=pars.mass;
		final double[] cx=gas.Cx,cy=gas.Cy,cz=gas.Cz,w=gas.w;
		for(int i=0;i<pkPtrs.length;i++){//go through all pk's in this cell ...
			final int pk=pkPtrs.get(i);
			w[pk]*=nfmb(cx[pk],cy[pk],cz[pk],to.ux,to.uy,to.uz,to.n0,mass0,to.T0)/nfmb(cx[pk],cy[pk],cz[pk],from.ux,from.uy,from.uz,from.n0,mass0,from.T0);
		}
	}






	/**return the average cell epsilon used in the last collision
	 * 
	 */
	public double avgKDENN(){
		double isum=0.,sum=0.;
		for(int i=0;i<KDEHistogram.length;i++){
			sum+=KDEHistogram[i];
			isum+=i*KDEHistogram[i];
		}
		double ibar=isum/sum;
		return ibar;
	}
	
	public int[] KDENNHistogram(){return KDEHistogram;}

	public void incSteps() {
		NE.incSteps();Eq.incSteps();		
	}



}




final class Sums implements Serializable{//data structure for the properties we want to keep track of everything
	private static final long serialVersionUID = 6L;
	Parameters params;
	public long Nsamples,Nsteps;//Arguably this should be part of the cell....
	private double w,x,y,xy,z,x2,y2,z2,x3,y3,z3,x4,y4,z4,w2;
	private double yx2,yy2,yz2;//the 3 variables to track the heat flux 
	@SuppressWarnings("unused")
	private Sums(){}//you can't create a Sums that does not have a reference to a params
	public void incSteps() {Nsteps++;	}
	public Sums(Parameters pars){params=pars;reset();}
	final public void reset(){Nsteps=Nsamples=0;w2=w=x=y=xy=z=x2=y2=z2=x3=y3=z3=x4=y4=z4=0;
		yx2=yy2=yz2=0;
	}//reset counters; 
	final public void add(Sums other) {
		Nsamples+=other.Nsamples;Nsteps+=other.Nsteps;
		w+=other.w;x+=other.x;y+=other.y;z+=other.z;
		xy+=other.xy;
		x2+=other.x2;y2+=other.y2;z2+=other.z2;
		x3+=other.x3;y3+=other.y3;z3+=other.z3;
		x4+=other.x4;y4+=other.y4;z4+=other.z4;
		yx2+=other.yx2;yy2+=other.yy2;yz2+=other.yz2;
		w2+=other.w2;
	}
	/** add pk with these properties to the appropriate sum
	 * 
	 * @param cx 
	 * @param cy
	 * @param cz
	 * @param w weight
	 */
	final public void addpk(final double cx, final double cy, final double cz, final double wi) {
		Nsamples++;
		//Variance in Weight:
		w2+=sq(wi);

		//density:
		w+=wi;
		//measure mean
		x+=wi*cx;y+=wi*cy;z+=wi*cz;
		//mixed moment:
		xy+=wi*cx*cy;
		//measure 2nd moment
		final double xsq=sq(cx),ysq=sq(cy),zsq=sq(cz);
		x2+=xsq*wi;y2+=ysq*wi;z2+=zsq*wi;
		//c3
		x3+=cx*xsq*wi;y3+=cy*ysq*wi;z3+=cz*zsq*wi;
		//measure c4
		final double xsqsq=sq(xsq),ysqsq=sq(ysq),zsqsq=sq(zsq);
		x4+=xsqsq*wi;y4+=ysqsq*wi;z4+=zsqsq*wi;

		// c3 to measure heat flow eventually
		yx2+=wi*cy*xsq;yy2+=wi*cy*ysq;yz2+=wi*cy*zsq;
	}
	/** Average pk sum over all Nsteps samples	 */
	public double pkSum(){return w/Nsteps;}
	/** Sum of Cx averaged over all pk's (I know it's confusing). Should re-write so that we just return the sum and have the caller figure out what to do with Nsamples*/
	public double Cx(){return x/Nsamples;}
	public double Cy(){return y/Nsamples;}
	public double Cz(){return z/Nsamples;}
	public double CxCy(){return xy/Nsamples;}
	/** approximation of T don't USE!
	 * @return rough estimate of temperature without density effects. DO NOT USE 
	 */
//	public double SloppyT(){
//		final double ux=ux(),uy=uy(),uz=uz(),m=params.mass,k=params.k;//nasty uy bug
//		final double tmp=m/(3.*k),Ti;
//		final double cx2=Cx2(),cy2=Cy2(),cz2=Cz2(),c2=cx2+cy2+cz2;
//		Ti=tmp*(c2-sq(ux)-sq(uy)-sq(uz));
//		return Ti/(1.0-1.0/Nsamples);}
	public double Cx2(){return x2/Nsamples;}//needed for the Homolle problem 
	public double Cy2(){return y2/Nsamples;}//needed for the Homolle problem 
	public double Cz2(){return z2/Nsamples;}//needed for the homolle problem 
	public double Cx3(){return x3/Nsamples;}//needed for the homolle problem 
	public double Cy3(){return y3/Nsamples;}//needed for the homolle problem 
	public double Cz3(){return z3/Nsamples;}//needed for the homolle problem 
	public double Cx4(){return x4/Nsamples;}//needed for the homolle problem 
	public double Cy4(){return y4/Nsamples;} 
	public double Cz4(){return z4/Nsamples;} 
	public double CyCx2(){return yx2/Nsamples;}
	public double CyCy2(){return yy2/Nsamples;}
	public double CyCz2(){return yz2/Nsamples;}
	public double varW(){return (w2-Nsamples*sq(w/Nsamples))/(Nsamples-1);}//unbias estimator of variance 
	public double w() {		return w;}

}
