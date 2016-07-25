import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import gurobi.GRBLinExpr;
import gurobi.GRBModel;
import gurobi.GRBVar;

public class ClassicalHLP_Main {
	private static final double[][] coordinates = MyArray.read("coordinates.txt");
	private static double[][] distances = Distance.get(coordinates);
	private static final double[][] tmpFlows = MyArray.read("w.txt");
	private static int nVar = tmpFlows.length;
	private static double[][] flows = new double[nVar][nVar];
	private static final double[][] fixedCosts = MyArray.read("fixedcharge.txt");
	private static final int P = 3;
	private static final double alpha = 0.2;

	public static void main(String[] args) throws GRBException {

		// Filling in the flows matrix assymetrically
		for (int i = 0; i < nVar; i++) {
			for (int j = 0; j < nVar; j++) {
				flows[i][j] = tmpFlows[i][j] + tmpFlows[j][i];
			}
		}

		GRBEnv env = new GRBEnv("RpHND.log");
		GRBModel model = new GRBModel(env);

		// Create variables
		GRBVar[][][][] x = new GRBVar[nVar][nVar][nVar][nVar];
		GRBVar[] y = new GRBVar[nVar];

		for (int i = 0; i < nVar; i++) {
			for (int j = i + 1; j < nVar; j++) {
				for (int k = 0; k < nVar; k++) {
					for (int m = 0; m < nVar; m++) {
						String varName = "x" + i + "_" + k + "_" + m + "_" + j;
						x[i][k][m][j] = model.addVar(0.0, 1.0, flows[i][j] * Cikmj(i, k, m, j), GRB.CONTINUOUS,
								varName);
					}
				}
			}
		}

		for (int i = 0; i < nVar; i++) {
			y[i] = model.addVar(0, 1.0, fixedCosts[i][0], GRB.CONTINUOUS, "y" + i);
		}

		// Integrate new variables
		model.update();

		// Adding constraints
		// Constraint 2
		GRBLinExpr con2 = new GRBLinExpr();
		for (int i = 0; i < nVar; i++) {
			con2.addTerm(1, y[i]);
		}
		model.addConstr(con2, GRB.EQUAL, P, "u2");

		// Constraint 3
		for (int i = 0; i < nVar; i++) {
			for (int j = i + 1; j < nVar; j++) {

				GRBLinExpr con3 = new GRBLinExpr();
				for (int k = 0; k < nVar; k++) {
					for (int m = 0; m < nVar; m++) {
						con3.addTerm(1, x[i][k][m][j]);
					}
				}
				model.addConstr(con3, GRB.EQUAL, 1, "u3_" + i + "_" + j);
			}
		}

		// Constraint 4
		for (int i = 0; i < nVar; i++) {
			for (int j = i + 1; j < nVar; j++) {
				for (int k = 0; k < nVar; k++) {

					GRBLinExpr con4 = new GRBLinExpr();
					for (int m = 0; m < nVar; m++) {
						con4.addTerm(1, x[i][k][m][j]);
					}
					con4.addTerm(-1, y[k]);
					model.addConstr(con4, GRB.LESS_EQUAL, 0, "u4_" + i + "_" + j + "_" + k);
				}
			}
		}

		// Constraint 5
		for (int i = 0; i < nVar; i++) {
			for (int j = i + 1; j < nVar; j++) {
				for (int m = 0; m < nVar; m++) {

					GRBLinExpr con5 = new GRBLinExpr();
					for (int k = 0; k < nVar; k++) {
						con5.addTerm(1, x[i][k][m][j]);
					}
					con5.addTerm(-1, y[m]);
					model.addConstr(con5, GRB.LESS_EQUAL, 0, "u5_" + i + "_" + j + "_" + m);
				}
			}
		}
		
		// optimize model
		model.optimize();
		model.write("d:/classicalModel.lp");
		printSol(model);
		/*for ( GRBVar var : model.getVars()){
			System.out.println(var.get(GRB.StringAttr.VarName) + ": " + var.get(GRB.DoubleAttr.Obj));
		}*/
	}
	
	/**
	 * prints the solution
	 * @param model
	 * @throws GRBException
	 */
	private static void printSol(GRBModel model) throws GRBException{
		for (GRBVar var : model.getVars()) 
			if (var.get(GRB.DoubleAttr.X)>0 /*&& var.get(GRB.StringAttr.VarName).contains("y")*/) 
				System.out.println(var.get(GRB.StringAttr.VarName) + " : " 
						+ var.get(GRB.DoubleAttr.X) + " - "
						+ var.get(GRB.DoubleAttr.Obj));
	}

	/**
	 * Cikmj
	 * 
	 */
	private static double Cikmj(int i, int k, int m, int j) {
		return distances[i][k] + (1 - alpha) * distances[k][m] + distances[m][j];
	}
}
