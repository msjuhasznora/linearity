package StochasticVariability;

// NO FILE OUTPUT, hopefully faster, for determining EXPERIMENTAL POE

// a refactored version of Sadegh's code

// We introduce a new cell type: a cell of "type 3" represents an infected cell that does not produce virus.
// The category is introduced in order to track the number of infected cells that got infected because of
// one specific cell.

// init description todo

import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GridWindow;
import HAL.Tools.FileIO;

import java.io.File;
import java.util.Random;

import static HAL.Util.*;

public class StochasticVariability {

	public static String experimentalOrTheoretical = "theoretical";
	public static int numberOfExperiments = 1000;

	public static void main(String[] args) {

		int numberOfTicks = 500; // per experiment. will be increased during runtime if necessary

		// this data set will store the number of new infections at each time step for each experiment
		// double[][] branchingProcessData = new double[numberOfExperiments][numberOfTicks];

		int nrExtinct = 0;

		java.util.Date now = new java.util.Date();
		java.text.SimpleDateFormat dateFormat = new java.text.SimpleDateFormat("yyyy-MM-dd_HH-mm-ss");
		String date_time = dateFormat.format(now);

		String projPath = PWD() + "/stochasticVariability";
		final String output_dir = projPath + "/output/" + experimentalOrTheoretical + "/" + date_time + "/";

		new File(output_dir).mkdirs();
		//FileIO outfile = new FileIO(output_dir.concat("/").concat("Out").concat(".csv"), "w");
		//FileIO dataLastTick = new FileIO(output_dir.concat("/").concat("dataLastTick").concat(".csv"), "w");

		int x = 200;
		int y = 200;

		int visScale = 2;
		GridWindow win = new GridWindow("Cellular state space, virus concentration.", x*2, y, visScale, true);

		System.out.println("Calculations for PoE type: " + experimentalOrTheoretical);

		for (int experimentIterator = 0; experimentIterator < numberOfExperiments; experimentIterator++){

			// System.out.println("Experiment ID: " + experimentIterator);

			StochasticExperiment experiment = new StochasticExperiment(x, y, experimentalOrTheoretical);
			experiment.Init(x, y);

			double[] cellCountResult = experiment.RunExperiment(numberOfTicks, win, experimentalOrTheoretical);

			if (experimentalOrTheoretical.equals("theoretical")){

				double numberOfType3Cells = cellCountResult[3];
				System.out.println("Number of offsprings: " + numberOfType3Cells);

			} else if (experimentalOrTheoretical.equals("experimental")){

				double remainingHealthyCells = cellCountResult[0];
				boolean isExtinct = (remainingHealthyCells >= (x * y) * 0.99);
				if (isExtinct) {
					nrExtinct += 1;
				}

			} else {
				System.out.println("PoE should be either experimental or theoretical.");
			}
		}

		if(experimentalOrTheoretical.equals("experimental")){
			System.out.println(StochasticExperiment.virusRemovalRate + "VRR. Extinct: " +  nrExtinct + " out of " + numberOfExperiments + " experiments.");
		}

		//for (int experimentIterator = 0; experimentIterator < branchingProcessData.length; experimentIterator++) {

		//for (int tick = 0; tick < branchingProcessData[0].length ; tick++){
		//	outfile.Write(branchingProcessData[experimentIterator][tick] +",");
		//}

		//outfile.Write("\n");
		//dataLastTick.Write(branchingProcessData[experimentIterator][branchingProcessData[0].length-1] + "\n");

		//}

		//outfile.Close();
		//dataLastTick.Close();

		win.Close();
	}
}


class StochasticExperiment extends AgentGrid2D<Cells>{

	public PDEGrid2D virusCon;
	public Random random;
	public String experimentalOrTheoretical;

	public double[] cellularVirusCon = new double[length];

	// SARS-CoV-2 parameters
	public static double virusRemovalRate = 0.1 * Math.pow(10,-3); // 1.67
	public double virusMax = 3.72 * Math.pow(10,-3);  // f_{i,j}
	public double virusDiffCoeff = 0.2; // D_V [sigma^2 / min]
	public double deathProb = 7.2 * Math.pow(10,-4);
	public double infectionRate = 1.01 * Math.pow(10,-7);

	public StochasticExperiment(int x, int y, String experimentalOrTheoretical){
		super(x, y, Cells.class);
		this.experimentalOrTheoretical = experimentalOrTheoretical;
		this.random = new Random();
		virusCon = new PDEGrid2D(xDim, yDim);
		virusCon.Update();
	}

	void Init(int x, int y){

		int initialPlace = length / 2 - length % 2 + x / 2 - x % 2; // right in the middle of the ABM

		for (int i = 0; i < length; i++){

			Cells c = NewAgentSQ(i);

			if(i == initialPlace) {
				c.CellInit(false,true,false,false);
			} else {
				c.CellInit(true,false,false,false);
			}
		}
	}

	double[] RunExperiment(int numberOfTicks, GridWindow win, String experimentalOrTheoretical){

		int numberOfEffectiveTicks = numberOfTicks;

		double[] cellCount = CountCells();

		for (int tick = 0; tick < numberOfEffectiveTicks; tick ++){

			TimeStep();
			//DrawModel(win);

			double totalVirusCon = TotalVirusCon();
			cellCount = CountCells();

			if (experimentalOrTheoretical.equals("experimental") && (cellCount[0] < xDim * yDim * 0.99)){
				break;
			}

			boolean doesItNeedMoreTime = ((tick == numberOfEffectiveTicks - 1) && ((totalVirusCon > 0.0001) || (cellCount[1] > 0)));
			if (doesItNeedMoreTime) {
				//System.out.println("Number of effective ticks: " + numberOfEffectiveTicks);
				numberOfEffectiveTicks += 100;
				//System.out.println("Number of effective ticks: " + numberOfEffectiveTicks);
				//System.out.println(totalVirusCon);
				//double[][] branchingProcessDataExt = new double[StochasticVariability.numberOfExperiments][numberOfEffectiveTicks];
				//branchingProcessData = branchingProcessDataExt;
				// fix copying the data from the previous one if needed, we don't need it now though
			}

			//branchingProcessData[experimentIterator][tick] = cellCount[3];
		}

		//System.out.println("Number of effective ticks: " + numberOfEffectiveTicks);
		return cellCount;
	}

	double[] CountCells(){

		double healthyCells = 0, infectedCells = 0, deadCells = 0, type3Cells = 0;
		double[] cellCount = new double[4];

		for (Cells cell: this){
			if (cell.CellType == 0){
				healthyCells += 1;
			} else if (cell.CellType == 1){
				infectedCells += 1;
			} else if (cell.CellType == 2){
				deadCells += 1;
			} else if (cell.CellType == 3){
				type3Cells += 1;
			}

		}

		cellCount[0] = healthyCells;
		cellCount[1] = infectedCells;
		cellCount[2] = deadCells;
		cellCount[3] = type3Cells;

		return cellCount;

	}

	double randomGenerator() {
		return random.nextDouble();
	}

	public void TimeStep(){

		TimeStepVirus();
		TimeStepCells();

	}

	void TimeStepCells(){

		for (Cells cell : this){
			cell.CellStepInfection(this.experimentalOrTheoretical);
		}

		for (Cells cell : this) {
			cell.CellStepDeath();
		}


	}

	void TimeStepVirus(){

		// decay of the virus
		for (Cells cell : this){
			virusCon.Add(cell.Isq(), -virusRemovalRate * virusCon.Get(cell.Isq()));
		}
		virusCon.Update();

		// virus production
		for (Cells cell : this){
			if (cell.CellType == 1){ // infected cell
				double addedVirusCon = VirusSource();
				double currentVirusCon = virusCon.Get(cell.Isq());
				double newVirusCon = addedVirusCon + currentVirusCon;
				virusCon.Set(cell.Isq(), newVirusCon);
			}
		}

		virusCon.DiffusionADI(virusDiffCoeff);
		virusCon.Update();

	}

	double TotalVirusCon() {

		double totalVirusCon = 0;
		for (int i = 0; i < length; i++){
			cellularVirusCon[i] = virusCon.Get(i);
		}

		for (double virusConInCell : cellularVirusCon ){
			totalVirusCon = totalVirusCon + virusConInCell;
		}

		return totalVirusCon;

	}

	double VirusSource(){

		return virusMax;

	}

	public void DrawModel(GridWindow vis){

		for (int i = 0; i < length; i++) {

			Cells drawMe = GetAgent(i);

			if (drawMe == null){
				vis.SetPix(i, RGB256(255, 255, 255));
			} else{

				if (drawMe.CellType == 0){ // healthy cells
					vis.SetPix(i, RGB256(0, 255, 51));
				}
				else if (drawMe.CellType == 1){ // infected cells
					vis.SetPix(i, RGB256(255, 0, 0));
				}
				else if (drawMe.CellType == 2){
					vis.SetPix(i, RGB256(0, 0, 0));
				}
				else if (drawMe.CellType == 3){
					vis.SetPix(i, RGB256(180, 3, 255));
				}

				vis.SetPix(ItoX(i) + xDim, ItoY(i), HeatMapRGB(100 * virusCon.Get(i)));

			}
		}
	}

}

class Cells extends AgentSQ2Dunstackable<StochasticExperiment>{

	int CellType;

	// cell types: 0 - healthy, 1 - infected, 2 - dead, 3 - an infected cell that does not produce virus

	public void CellInit(boolean isHealthy, boolean isInfected, boolean isDead, boolean isType3){

		if(isHealthy == true){
			this.CellType = 0;
		} else if(isInfected == true){
			this.CellType = 1;
		} else if(isDead == true) {
			this.CellType = 2;
		} else if(isType3 == true) {
			this.CellType = 3;
		}
	}

	public void CellStepInfection(String experimentalOrTheoretical){

		double virusConAtCell = G.virusCon.Get(Isq());

		double infectionProb = G.infectionRate * G.xDim * G.yDim * virusConAtCell;
		if (G.randomGenerator() < infectionProb) {

			if (this.CellType == 0){ // healthy cell
				if (experimentalOrTheoretical.equals("theoretical")) {
					this.CellType = 3; // a "type 3" cell is an infected cell that does not produce virus
				} else {
					this.CellType = 1;
				}
			}
		}
	}

	public void CellStepDeath(){

		if (this.CellType == 1) {

			if(G.randomGenerator() < G.deathProb){
				this.CellType = 2;
			}
		}
	}

}