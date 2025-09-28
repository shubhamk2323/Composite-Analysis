import java.util.*;

public class Macromechanics{

    public static class LaminateStiffnessProperties {
        public double Ex;
        public double Ey;
        public double Gxy;
        public double Exb;
        public double Eyb;
    
        public LaminateStiffnessProperties(double Ex, double Ey, double Gxy, double Exb, double Eyb) {
            this.Ex = Ex;
            this.Ey = Ey;
            this.Gxy = Gxy;
            this.Exb = Exb;
            this.Eyb = Eyb;
        }
    }

    public static class MaterialProperties {
        int id;
        double E1, E2, G12, v12, v21;
        double Lt, Lc, Tt, Tc, s;

        public MaterialProperties(int id, double E1, double E2, double G12, double v12,
        double Lt, double Lc, double Tt, double Tc, double s) {
            this.id= id;
            this.E1 = E1;
            this.E2 = E2;
            this.G12 = G12;
            this.v12 = v12;
            this.v21 = (E2 / E1) * v12;
            this.Lt= Lt;
            this.Lc= Lc;
            this.Tt= Tt;
            this.Tc= Tc;
            this.s= s;
        }
    }

    public static class layup{
        int id;
        int angle;
        double t;
        public layup(int id, int angle, double t){
            this.id= id;
            this.angle= angle;
            this.t= t;
        }
    }

    public static class load{
        double nx, ny, nxy, mx, my, mxy;
        public load(double nx, double ny, double nxy, double mx,
        double my, double mxy){
            this.nx= nx;
            this.ny= ny;
            this.nxy= nxy;
            this.mx= mx;
            this.my= my;
            this.mxy= mxy;
        }
    }

    public static double[][] Smat(MaterialProperties m){
        double s11= 1/(m.E1);
        double s12= -(m.v12)/(m.E1);
        double s21= s12;
        double s22= 1/(m.E2);
        double s66= 1/(m.G12);
        double[][] smat= new double[3][3];
        smat[0]= new double[]{s11, s12, 0};
        smat[1]= new double[]{s21, s22, 0};
        smat[2]= new double[]{0, 0, s66};
        return smat;
    }

    public static double[][] Qmat(MaterialProperties m){
        double[][] A= Smat(m);

        double[][] inv= new double[3][3];

        double det = A[0][0]*A[1][1]*A[2][2]
                   + A[0][1]*A[1][2]*A[2][0]
                   + A[0][2]*A[1][0]*A[2][1]
                   - A[0][2]*A[1][1]*A[2][0]
                   - A[0][0]*A[1][2]*A[2][1]
                   - A[0][1]*A[1][0]*A[2][2];

        inv[0][0]= (A[1][1]*A[2][2]- A[1][2]*A[2][1]) / det;
        inv[0][1]= (A[0][2]*A[2][1]- A[0][1]*A[2][2]) / det;
        inv[0][2]= (A[0][1]*A[1][2]- A[0][2]*A[1][1]) / det;

        inv[1][0]= (A[1][2]*A[2][0]- A[1][0]*A[2][2]) / det;
        inv[1][1]= (A[0][0]*A[2][2]- A[0][2]*A[2][0]) / det;
        inv[1][2]= (A[0][2]*A[1][0]- A[0][0]*A[1][2]) / det;

        inv[2][0]= (A[1][0]*A[2][1]- A[1][1]*A[2][0]) / det;
        inv[2][1]= (A[0][1]*A[2][0]- A[0][0]*A[2][1]) / det;
        inv[2][2]= (A[0][0]*A[1][1]- A[0][1]*A[1][0]) / det;

        return inv;

    }

    public static double[][] getSbarFromS(double[][] S, layup l) {
        double theta= Math.toRadians(l.angle);
        double c= Math.cos(theta);
        double s= Math.sin(theta);
        double c2= c*c;
        double s2= s*s;
        double sc= s*c;
        double[][] R = {{1, 0, 0}, {0, 1, 0}, {0, 0, 2}};
        double[][] Rinv = {{1, 0, 0}, {0, 1, 0}, {0, 0, 0.5}};

        double[][] T = {
            {c2, s2, 2*sc},
            {s2, c2, -2*sc},
            {-sc, sc, c2-s2}
        };

        double[][] Tinv = {
            {c2, s2, -2*sc},
            {s2, c2, 2*sc},
            {sc, -sc, c2 - s2}
        };

        double[][] RTinv= multiplyMatrices(R, Tinv);
        double[][] RTRinv= multiplyMatrices(RTinv, Rinv);
        double[][] A= multiplyMatrices(RTRinv, S);
        double[][] Sbar= multiplyMatrices(A, T);
        return Sbar;
    }

    public static double[][] multiplyMatrices(double[][] A, double[][] B) {
        int rowsA= A.length;
        int colsA= A[0].length;
        int colsB= B[0].length;
        double[][] result= new double[rowsA][colsB];

        for(int i=0; i<rowsA; i++) {
            for (int j=0; j<colsB; j++) {
                for (int k=0; k<colsA; k++) {
                    result[i][j]+= A[i][k]*B[k][j];
                }
            }
        }
        return result;
    }

    public static double[][] getQbarFromQ(double[][] Q, layup l) {
        double theta= Math.toRadians(l.angle);
        double c= Math.cos(theta);
        double s= Math.sin(theta);
        double c2= c * c;
        double s2= s * s;
        double sc= s * c;

        double[][] R= {{1, 0, 0},{0, 1, 0},{0, 0, 2}};
        double[][] Rinv= {{1, 0, 0},{0, 1, 0},{0, 0, 0.5}};
        double[][] T= {
            {c2, s2, 2 * sc},
            {s2, c2, -2 * sc},
            {-sc, sc, c2 - s2}
        };
    
        double[][] Tinv= {
            {c2, s2, -2 * sc},
            {s2, c2, 2 * sc},
            {sc, -sc, c2 - s2}
        };
        double[][] TinvQ= multiplyMatrices(Tinv, Q);
        double[][] A= multiplyMatrices(TinvQ, R);
        double[][] B= multiplyMatrices(A, T);
        double[][] Qbar= multiplyMatrices(B, Rinv);
        return Qbar;
    }

    public static double[][] computeAMatrix(List<layup> layups, Map<Integer, MaterialProperties> materialMap) {
        double[][] A = new double[3][3];

        double totalThickness = 0;
        for (layup l : layups) {
            totalThickness += l.t;
        }

        double zBottom = -totalThickness / 2;

        for (layup l : layups) {
            double zTop = zBottom + l.t;

            MaterialProperties m = materialMap.get(l.id);
            double[][] Q = Qmat(m);              
            double[][] Qbar = getQbarFromQ(Q, l);

            double dz = zTop - zBottom;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    A[i][j] += Qbar[i][j] * dz;
                }
            }

            zBottom = zTop;
        }

        return A;
    }

    public static double[][] computeBMatrix(List<layup> layups, Map<Integer, MaterialProperties> materialMap) {
        double[][] B = new double[3][3];

        double totalThickness = 0;
        for (layup l : layups) {
            totalThickness += l.t;
        }

        double zBottom = -totalThickness / 2;

        for (layup l : layups) {
            double zTop = zBottom + l.t;

            MaterialProperties m = materialMap.get(l.id);
            double[][] Q = Qmat(m);
            double[][] Qbar = getQbarFromQ(Q, l);

            double zTop2 = zTop * zTop;
            double zBottom2 = zBottom * zBottom;
            double dz2 = zTop2 - zBottom2;

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    B[i][j] += 0.5 * Qbar[i][j] * dz2;
                }
            }

            zBottom = zTop;
        }

        return B;
    }

    public static double[][] computeDMatrix(List<layup> layups, Map<Integer, MaterialProperties> materialMap) {
        double[][] D = new double[3][3];

        double totalThickness = 0;
        for (layup l : layups) {
            totalThickness += l.t;
        }

        double zBottom = -totalThickness / 2;

        for (layup l : layups) {
            double zTop = zBottom + l.t;

            MaterialProperties m = materialMap.get(l.id);
            double[][] Q = Qmat(m);
            double[][] Qbar = getQbarFromQ(Q, l);

            double zTop3 = zTop * zTop * zTop;
            double zBottom3 = zBottom * zBottom * zBottom;
            double dz3 = zTop3 - zBottom3;

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    D[i][j] += (1.0 / 3.0) * Qbar[i][j] * dz3;
                }
            }

            zBottom = zTop;
        }

        return D;
    }


    public static double[][] makeABD(double[][] A, double[][] B, double[][] D) {
        double[][] abd = new double[6][6];

        // Fill A
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                abd[i][j] = A[i][j];

        // Fill B
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                abd[i][j + 3] = B[i][j];
                abd[i + 3][j] = B[i][j];
            }

        // Fill D
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                abd[i + 3][j + 3] = D[i][j];

        return abd;
    }



    public static double[][] invert6x6(double[][] matrix) {
        int n = 6;
        double[][] a = new double[n][n];
        double[][] inverse = new double[n][n];

        for (int i = 0; i < n; i++) {
            System.arraycopy(matrix[i], 0, a[i], 0, n);
            inverse[i][i] = 1;
        }

        for (int i = 0; i < n; i++) {
            double pivot = a[i][i];
            if (pivot == 0) {
                throw new ArithmeticException("Matrix is singular and cannot be inverted.");
            }

            for (int j = 0; j < n; j++) {
                a[i][j] /= pivot;
                inverse[i][j] /= pivot;
            }

            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double factor = a[k][i];
                    for (int j = 0; j < n; j++) {
                        a[k][j] -= factor * a[i][j];
                        inverse[k][j] -= factor * inverse[i][j];
                    }
                }
            }
        }

        return inverse;
    }

    // public static LaminateStiffnessProperties computeApparentStiffness(double[][] A, double[][] D, List<layup> layups) {
    //     double h = layups.stream().mapToDouble(l -> l.t).sum();

    //     double A11inv = 1.0 / A[0][0];
    //     double A22inv = 1.0 / A[1][1];
    //     double A66inv = 1.0 / A[2][2];

    //     double D11inv = 1.0 / D[0][0];
    //     double D22inv = 1.0 / D[1][1];

    //     double Ex = (1.0 / h) * (1.0 / A11inv);
    //     double Ey = (1.0 / h) * (1.0 / A22inv);
    //     double Gxy = (1.0 / h) * (1.0 / A66inv);

    //     double Exb = (12.0 / Math.pow(h, 3)) * (1.0 / D11inv);
    //     double Eyb = (12.0 / Math.pow(h, 3)) * (1.0 / D22inv);

    //     return new LaminateStiffnessProperties(Ex, Ey, Gxy, Exb, Eyb);
    // }

    public static LaminateStiffnessProperties computeApparentStiffness(double[][] A, double[][] D, List<layup> layups) {
        // Step 1: Calculate total thickness
        double totalThickness = 0;
        for (layup l : layups) {
            totalThickness += l.t;
        }

        // Step 2: Normalize A matrix by total thickness
        double[][] A_bar = new double[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                A_bar[i][j] = A[i][j] / totalThickness;
            }
        }

        // Step 3: Inverse of A̅ and D matrices
        double[][] AbarInv = invert3x3Matrix(A_bar);
        double[][] DInv = invert3x3Matrix(D);

        // Step 4: Calculate properties
        double Ex = 1.0 / AbarInv[0][0];      // Longitudinal stiffness
        double Ey = 1.0 / AbarInv[1][1];      // Transverse stiffness
        double Gxy = 1.0 / AbarInv[2][2];     // Shear stiffness

        double Exb = 12.0 / DInv[0][0];       // Bending stiffness in x
        double Eyb = 12.0 / DInv[1][1];       // Bending stiffness in y

        // Step 5: Return result
        return new LaminateStiffnessProperties(Ex, Ey, Gxy, Exb, Eyb);
    }


    public static double[][] invert3x3Matrix(double[][] m) {
        double[][] inv = new double[3][3];

        double det = m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1])
                   - m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0])
                   + m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);

        if (Math.abs(det) < 1e-12) throw new RuntimeException("Matrix not invertible");

        inv[0][0] =  (m[1][1]*m[2][2] - m[1][2]*m[2][1]) / det;
        inv[0][1] = -(m[0][1]*m[2][2] - m[0][2]*m[2][1]) / det;
        inv[0][2] =  (m[0][1]*m[1][2] - m[0][2]*m[1][1]) / det;

        inv[1][0] = -(m[1][0]*m[2][2] - m[1][2]*m[2][0]) / det;
        inv[1][1] =  (m[0][0]*m[2][2] - m[0][2]*m[2][0]) / det;
        inv[1][2] = -(m[0][0]*m[1][2] - m[0][2]*m[1][0]) / det;

        inv[2][0] =  (m[1][0]*m[2][1] - m[1][1]*m[2][0]) / det;
        inv[2][1] = -(m[0][0]*m[2][1] - m[0][1]*m[2][0]) / det;
        inv[2][2] =  (m[0][0]*m[1][1] - m[0][1]*m[1][0]) / det;

        return inv;
    }

    public static void printMatrix(double[][] matrix) {
        for (double[] row : matrix) {
            for (double val : row) {
                System.out.printf("%12.3e ", val);
            }
            System.out.println();
        }
    }

    public static double[] matrixVectorMultiply(double[][] matrix, double[] vector) {
        int rows = matrix.length;
        int cols = matrix[0].length;
        double[] result = new double[rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i] += matrix[i][j] * vector[j];
            }
        }
        return result;
    }


    public static class PlyStressStrain {
        public double[][] strainLower;
        public double[][] strainUpper;
        public double[][] strainMid;
        public double[][] stressLower;
        public double[][] stressUpper;
        public double[][] stressMid;

        public PlyStressStrain(int n) {
            strainLower = new double[n][3];
            strainUpper = new double[n][3];
            strainMid = new double[n][3];
            stressLower = new double[n][3];
            stressUpper = new double[n][3];
            stressMid = new double[n][3];
        }
    }

    public static PlyStressStrain calculateStressStrain(
            int n,
            double[][] ABD_inv,
            double[][][] QbarAll,
            double[] lowerBounds,
            double[] upperBounds,
            load inputLoad
    ) {
        // Step 1: Create force vector
        double[][] forceVector = {
                {inputLoad.nx},
                {inputLoad.ny},
                {inputLoad.nxy},
                {inputLoad.mx},
                {inputLoad.my},
                {inputLoad.mxy}
        };

        // Step 2: midplane strain and curvature
        double[][] strainCurvature = multiplyMatrices(ABD_inv, forceVector); // 6x1

        PlyStressStrain result = new PlyStressStrain(n);

        for (int i = 0; i < n; i++) {
            double zLower = lowerBounds[i];
            double zUpper = upperBounds[i];
            double zMid = (zLower + zUpper) / 2.0;

            // Compute strain at each surface
            for (int j = 0; j < 3; j++) {
                result.strainLower[i][j] = strainCurvature[j][0] + zLower * strainCurvature[j + 3][0];
                result.strainUpper[i][j] = strainCurvature[j][0] + zUpper * strainCurvature[j + 3][0];
                result.strainMid[i][j]   = strainCurvature[j][0] + zMid   * strainCurvature[j + 3][0];
            }

            // Multiply with Q̄ to get stress
            result.stressLower[i] = matrixVectorMultiply(QbarAll[i], result.strainLower[i]);
            result.stressUpper[i] = matrixVectorMultiply(QbarAll[i], result.strainUpper[i]);
            result.stressMid[i]   = matrixVectorMultiply(QbarAll[i], result.strainMid[i]);
        }

        return result;
    }

    public static double[] transformStress(double[] stress, double angleDegrees) {
        double theta = Math.toRadians(angleDegrees);
        double c = Math.cos(theta);
        double s = Math.sin(theta);
        double c2 = c * c;
        double s2 = s * s;
        double sc = s * c;

        // Transformation matrix
        double[][] T = {
                {c2, s2, 2 * sc},
                {s2, c2, -2 * sc},
                {-sc, sc, c2 - s2}
        };

        // Perform matrix-vector multiplication
        return matrixVectorMultiply(T, stress);
    }

    public static double[] transformStrain(double[] globalStrain, int angleDegrees) {
        double theta = Math.toRadians(angleDegrees);
        double c = Math.cos(theta);
        double s = Math.sin(theta);
        double c2 = c * c;
        double s2 = s * s;
        double sc = s * c;

        double[][] R = {
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 2}
        };
        double[][] Rinv = {
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 0.5}
        };

        double[][] T = {
            {c2, s2, 2 * sc},
            {s2, c2, -2 * sc},
            {-sc, sc, c2 - s2}
        };

        // Multiply R * T
        double[][] RT = multiplyMatrices(R, T);
        // Multiply (R * T) * Rinv
        double[][] RTRinv = multiplyMatrices(RT, Rinv);

        // Now do matrix-vector multiplication with global strain
        double[] localStrain = new double[3];
        for (int i = 0; i < 3; i++) {
            localStrain[i] = 0;
            for (int j = 0; j < 3; j++) {
                localStrain[i] += RTRinv[i][j] * globalStrain[j];
            }
        }

        return localStrain;
    }



    public static void ans() {
        Scanner sc = new Scanner(System.in);

        // Input Material Properties
        System.out.println("Enter Material Properties:");
        System.out.print("ID: ");
        int id = sc.nextInt();
        System.out.print("E1 :");
        double E1 = sc.nextDouble();
        System.out.print("E2 :");
        double E2 = sc.nextDouble();
        System.out.print("G12 :");
        double G12 = sc.nextDouble();
        System.out.print("v12: ");
        double v12 = sc.nextDouble();
        System.out.print("Lt: ");
        double Lt = sc.nextDouble();
        System.out.print("Lc: ");
        double Lc = sc.nextDouble();
        System.out.print("Tt: ");
        double Tt = sc.nextDouble();
        System.out.print("Tc: ");
        double Tc = sc.nextDouble();
        System.out.print("s: ");
        double ss = sc.nextDouble();

        MaterialProperties m = new MaterialProperties(id, E1, E2, G12, v12, Lt, Lc, Tt, Tc, ss);

        // Input Layup
        System.out.println("\nEnter Layup:");
        System.out.print("ID: ");
        int layupId = sc.nextInt();
        System.out.print("Angle (degrees): ");
        int angle = sc.nextInt();
        System.out.print("Thickness (mm): ");
        double t = sc.nextDouble();
        layup l = new layup(layupId, angle, t);

        // Calculate S
        double[][] S = Smat(m);

        // Calculate Q
        double[][] Q = Qmat(m);

        // Calculate S̄ and Q̄
        double[][] Sbar = getSbarFromS(S, l);
        double[][] Qbar = getQbarFromQ(Q, l);

        // Print S
        System.out.println("\nMatrix S:");
        printMatrix(S);

        // Print S̄
        System.out.println("\nMatrix S̄:");
        printMatrix(Sbar);

        // Print Q
        System.out.println("\nMatrix Q:");
        printMatrix(Q);

        // Print Q̄
        System.out.println("\nMatrix Q̄:");
        printMatrix(Qbar);

        System.out.print("\nEnter angle θ (in degrees): ");
        double theta = Math.toRadians(sc.nextDouble());
        double c = Math.cos(theta);
        double s = Math.sin(theta);
        double c2 = c * c;
        double s2 = s * s;
        double c4 = c2 * c2;
        double s4 = s2 * s2;
        double sin2theta = Math.sin(2 * theta);
        double sin2theta2 = sin2theta * sin2theta;
        double cos2theta2 = Math.cos(2 * theta) * Math.cos(2 * theta);

        // 1 / Exx
        double inv_Exx = c4 / E1 + s4 / E2 + 0.25 * (1.0 / G12 - 2 * v12 / E1) * sin2theta2;
        double Exx = 1.0 / inv_Exx;

        // 1 / Eyy
        double inv_Eyy = s4 / E1 + c4 / E2 + 0.25 * (1.0 / G12 - 2 * v12 / E1) * sin2theta2;
        double Eyy = 1.0 / inv_Eyy;

        // 1 / Gxy
        double midz= 1.0 / E1 + 2 * v12 / E1 + 1.0 / E2;
        double inv_Gxy = midz - ((midz - (1.0 / G12)) * cos2theta2);
        double Gxy = 1.0 / inv_Gxy;

        // νxyθ
        double term = (v12 / E1 - 0.25 * sin2theta2 * (1.0 / E1 + 2 * v12 / E1 + 1.0 / E2 - 1.0 / G12));
        double vxy_theta = Exx * term;

        // νyxθ
        double vyx_theta = (Eyy / Exx) * vxy_theta;

        // Output
        System.out.printf("Exx (θ): %.3f GPa\n", Exx / 1e9);
        System.out.printf("Eyy (θ): %.3f GPa\n", Eyy / 1e9);
        System.out.printf("Gxy (θ): %.3f GPa\n", Gxy / 1e9);
        System.out.printf("νxy (θ): %.5f\n", vxy_theta);
        System.out.printf("νyx (θ): %.5f\n", vyx_theta);

    }

    public static void answ() {
        Scanner sc = new Scanner(System.in);

        System.out.print("Enter E1");
        double E1 = sc.nextDouble();

        System.out.print("Enter E2");
        double E2 = sc.nextDouble();

        System.out.print("Enter G12");
        double G12 = sc.nextDouble();

        System.out.print("Enter v12: ");
        double v12 = sc.nextDouble();

        double v21 = (E2 / E1) * v12;

        MaterialProperties mat = new MaterialProperties(1, E1, E2, G12, v12, 0, 0, 0, 0, 0);

        System.out.print("Enter σx ");
        double sx = sc.nextDouble() ;

        System.out.print("Enter σy ");
        double sy = sc.nextDouble() ;

        System.out.print("Enter τxy ");
        double txy = sc.nextDouble() ;

        System.out.print("Enter fiber angle θ (degrees): ");
        double angle = sc.nextDouble();

        double[] stressXY = {sx, sy, txy};
        double[] stress12 = transformStress(stressXY, angle);
        

        double[][] S = {
            {1 / E1, -v12 / E1, 0},
            {-v12 / E1, 1 / E2, 0},
            {0, 0, 1 / G12}
        };

        double[] strain12 = matrixVectorMultiply(S, stress12);

        double[] strainXY = transformStrain(strain12, (int)angle);



        System.out.println("\nStress in material direction :");
        System.out.printf("sig1   = %.3e%n", stress12[0]);
        System.out.printf("sig2   = %.3e%n", stress12[1]);
        System.out.printf("gamma12  = %.3e%n", stress12[2]);


        System.out.println("\nStrain in given direction (ε1, ε2, γ12):");
        System.out.printf("ε1   = %.3e%n", strainXY[0]);
        System.out.printf("ε2   = %.3e%n", strainXY[1]);
        System.out.printf("γ12  = %.3e%n", strainXY[2]);


        System.out.println("\nStrain in material direction (ε1, ε2, γ12):");
        System.out.printf("ε1   = %.3e%n", strain12[0]);
        System.out.printf("ε2   = %.3e%n", strain12[1]);
        System.out.printf("γ12  = %.3e%n", strain12[2]);

        sc.close();
    }

    public static void printPlyStrengthFromUser() {
        Scanner sc = new Scanner(System.in);

        // --- 1. Material Properties Input ---
        System.out.print("Enter E1 (Pa): ");
        double E1 = sc.nextDouble();

        System.out.print("Enter E2 (Pa): ");
        double E2 = sc.nextDouble();

        System.out.print("Enter G12 (Pa): ");
        double G12 = sc.nextDouble();

        System.out.print("Enter v12: ");
        double v12 = sc.nextDouble();

        System.out.print("Enter Longitudinal Tensile Strength Lt (Pa): ");
        double Lt = sc.nextDouble();

        System.out.print("Enter Longitudinal Compressive Strength Lc (Pa): ");
        double Lc = sc.nextDouble();

        System.out.print("Enter Transverse Tensile Strength Tt (Pa): ");
        double Tt = sc.nextDouble();

        System.out.print("Enter Transverse Compressive Strength Tc (Pa): ");
        double Tc = sc.nextDouble();

        System.out.print("Enter Shear Strength s (Pa): ");
        double s = sc.nextDouble();

        System.out.print("Enter fiber angle θ (degrees): ");
        double angle = sc.nextDouble();

        // --- 2. Create MaterialProperties object ---
        MaterialProperties mat = new MaterialProperties(1, E1, E2, G12, v12, Lt, Lc, Tt, Tc, s);

        // --- 3. Local stress vectors ---
        double[] tensionLocal = {Lt, 0, 0};
        double[] compressionLocal = {-Lc, 0, 0};
        double[] transTensionLocal = {0, Tt, 0};
        double[] transCompressionLocal = {0, -Tc, 0};
        double[] shearLocal = {0, 0, s};

        // --- 4. Transform to global (XY) direction ---
        double[] tensionXY = transformStress(tensionLocal, angle);
        double[] compressionXY = transformStress(compressionLocal, angle);
        double[] transTensionXY = transformStress(transTensionLocal, angle);
        double[] transCompressionXY = transformStress(transCompressionLocal, angle);
        double[] shearXY = transformStress(shearLocal, angle);

        // --- 5. Print output ---
        System.out.println("\nAllowable strengths (transformed to global XY coordinates):");
        System.out.printf("Tensile in X-dir     = %.2f MPa\n", tensionXY[0] / 1e6);
        System.out.printf("Compressive in X-dir = %.2f MPa\n", compressionXY[0] / 1e6);
        System.out.printf("Tensile in Y-dir     = %.2f MPa\n", transTensionXY[1] / 1e6);
        System.out.printf("Compressive in Y-dir = %.2f MPa\n", transCompressionXY[1] / 1e6);
        System.out.printf("Shear in XY-dir      = %.2f MPa\n", shearXY[2] / 1e6);

        sc.close();
    }

    public static void computeABDwithMappedMaterials() {
        Scanner sc = new Scanner(System.in);

        // Step 1: Input number of materials
        System.out.print("Enter number of unique materials: ");
        int numMaterials = sc.nextInt();

        Map<Integer, MaterialProperties> materialMap = new HashMap<>();

        // Input each material
        for (int i = 0; i < numMaterials; i++) {
            System.out.println("\n--- Material " + (i + 1) + " ---");
            System.out.print("Material ID: ");
            int id = sc.nextInt();

            System.out.print("E1 (Pa): ");
            double E1 = sc.nextDouble();

            System.out.print("E2 (Pa): ");
            double E2 = sc.nextDouble();

            System.out.print("G12 (Pa): ");
            double G12 = sc.nextDouble();

            System.out.print("v12: ");
            double v12 = sc.nextDouble();

            

            materialMap.put(id, new MaterialProperties(id, E1, E2, G12, v12, 0, 0, 0, 0, 0));
        }

        // Step 2: Input number of layups
        System.out.print("\nEnter number of layup plies: ");
        int numLayers = sc.nextInt();

        List<layup> l = new ArrayList<>();
        double totalThickness = 0;

        for (int i = 0; i < numLayers; i++) {
            System.out.println("\n--- Ply " + (i + 1) + " ---");

            System.out.print("Enter material ID: ");
            int matId = sc.nextInt();

            if (!materialMap.containsKey(matId)) {
                System.out.println("Invalid material ID. Exiting.");
                return;
            }

            System.out.print("Enter fiber angle (deg): ");
            int angle = sc.nextInt();

            System.out.print("Enter thickness (m): ");
            double t = sc.nextDouble();

            l.add(new layup(matId, angle, t));
            totalThickness += t;
        }

        // Step 3: Compute z-coordinates for plies
        double[] z = new double[numLayers + 1];
        z[0] = -totalThickness / 2;

        for (int i = 1; i <= numLayers; i++) {
            z[i] = z[i - 1] + l.get(i - 1).t;
        }


        // Step 4: Calculate A, B, D matrices
        double[][] A = computeAMatrix(l, materialMap);
        double[][] B = computeBMatrix(l, materialMap);
        double[][] D = computeDMatrix(l, materialMap);
        double[][] ABD = makeABD(A, B, D);

        // Step 6: Invert ABD
        double[][] ABDinv = invert6x6(ABD);
        System.out.println("\n--- ABD Matrix ---"); printMatrix(ABD);
        System.out.println("\n--- Inverse ABD Matrix ---"); printMatrix(ABDinv);

        sc.close();
    }

    public static void strengthTransformUserInput() {
        Scanner sc = new Scanner(System.in);
        sc.useLocale(Locale.US); // To ensure dot-based decimal

        // Input
        System.out.print("Enter Xt (Pa): ");
        double Xt = sc.nextDouble();

        System.out.print("Enter Yt (Pa): ");
        double Yt = sc.nextDouble();

        System.out.print("Enter Xc (Pa): ");
        double Xc = sc.nextDouble();

        System.out.print("Enter Yc (Pa): ");
        double Yc = sc.nextDouble();

        System.out.print("Enter shear strength S (Pa): ");
        double S = sc.nextDouble();

        System.out.print("Enter fiber angle θ (in degrees): ");
        double theta = sc.nextDouble();

        // Transform
        double thetaRad = Math.toRadians(theta);
        double sin = Math.sin(thetaRad);
        double cos = Math.cos(thetaRad);

        double[][] T = {
            {cos * cos, sin * sin, 2 * sin * cos},
            {sin * sin, cos * cos, -2 * sin * cos},
            {-sin * cos, sin * cos, cos * cos - sin * sin}
        };

        double[][] Tinv = invert3x3Matrix(T);

        double[] sigma12t = {Xt, Yt, S};
        double[] sigmaXYt = matrixVectorMultiply(Tinv, sigma12t);

        double[] sigma12c = {Xc, Yc, S};
        double[] sigmaXYc = matrixVectorMultiply(Tinv, sigma12c);

        // Output in scientific notation
        System.out.println("\n--- Transformed Strengths in Global XY Direction ---");

        System.out.printf("Tensile Limit (σx_t, σy_t, τxy_t): %n");
        System.out.printf("σx_t  = %.6e Pa%n", sigmaXYt[0]);
        System.out.printf("σy_t  = %.6e Pa%n", sigmaXYt[1]);
        System.out.printf("τxy_t = %.6e Pa%n", sigmaXYt[2]);

        System.out.printf("\nCompressive Limit (σx_c, σy_c, τxy_c): %n");
        System.out.printf("σx_c  = %.6e Pa%n", sigmaXYc[0]);
        System.out.printf("σy_c  = %.6e Pa%n", sigmaXYc[1]);
        System.out.printf("τxy_c = %.6e Pa%n", sigmaXYc[2]);

        sc.close();
    }



    public static void fullLaminateAnalysis() {
        Scanner sc = new Scanner(System.in);

        // Step 1: Material input
        System.out.print("Enter number of materials: ");
        int numMaterials = sc.nextInt();
        Map<Integer, MaterialProperties> materialMap = new HashMap<>();

        for (int i = 0; i < numMaterials; i++) {
            System.out.println("\n--- Material " + (i + 1) + " ---");
            System.out.print("Enter ID: ");
            int id = sc.nextInt();
            System.out.print("E1 (Pa): "); double E1 = sc.nextDouble();
            System.out.print("E2 (Pa): "); double E2 = sc.nextDouble();
            System.out.print("G12 (Pa): "); double G12 = sc.nextDouble();
            System.out.print("v12: "); double v12 = sc.nextDouble();
            System.out.print("Xt (Pa): "); double Xt = sc.nextDouble();
            System.out.print("Xc (Pa): "); double Xc = sc.nextDouble();
            System.out.print("Yt (Pa): "); double Yt = sc.nextDouble();
            System.out.print("Yc (Pa): "); double Yc = sc.nextDouble();
            System.out.print("S (Pa): "); double S = sc.nextDouble();
            materialMap.put(id, new MaterialProperties(id, E1, E2, G12, v12, Xt, Xc, Yt, Yc, S));
        }

        // Step 2: Layup input
        System.out.print("\nEnter number of layers: ");
        int numLayers = sc.nextInt();
        List<layup> layups = new ArrayList<>();
        double totalThickness = 0;

        for (int i = 0; i < numLayers; i++) {
            System.out.println("\n--- Ply " + (i + 1) + " ---");
            System.out.print("Material ID: ");
            int id = sc.nextInt();
            if (!materialMap.containsKey(id)) {
                System.out.println("Invalid ID. Exiting.");
                return;
            }
            System.out.print("Angle (deg): ");
            int angle = sc.nextInt();
            System.out.print("Thickness (m): ");
            double t = sc.nextDouble();
            totalThickness += t;
            layups.add(new layup(id, angle, t));
        }

        // Step 3: Compute z-coordinates
        double[] lowerBounds = new double[numLayers];
        double[] upperBounds = new double[numLayers];
        double zBottom = -totalThickness / 2;
        for (int i = 0; i < numLayers; i++) {
            lowerBounds[i] = zBottom;
            upperBounds[i] = zBottom + layups.get(i).t;
            zBottom = upperBounds[i];
        }

        // Step 4: Compute Q̄ for all plies
        double[][][] QbarAll = new double[numLayers][3][3];
        for (int i = 0; i < numLayers; i++) {
            MaterialProperties mat = materialMap.get(layups.get(i).id);
            double[][] Q = Qmat(mat);
            QbarAll[i] = getQbarFromQ(Q, layups.get(i));
        }

        // Step 5: Compute A, B, D and ABD⁻¹
        double[][] A = computeAMatrix(layups, materialMap);
        double[][] B = computeBMatrix(layups, materialMap);
        double[][] D = computeDMatrix(layups, materialMap);
        double[][] ABD = makeABD(A, B, D);
        double[][] ABD_inv = invert6x6(ABD);

        // Step 6: Load input
        System.out.println("\nEnter loads:");
        System.out.print("Nx (N/m): "); double nx = sc.nextDouble();
        System.out.print("Ny (N/m): "); double ny = sc.nextDouble();
        System.out.print("Nxy (N/m): "); double nxy = sc.nextDouble();
        System.out.print("Mx (N-m/m): "); double mx = sc.nextDouble();
        System.out.print("My (N-m/m): "); double my = sc.nextDouble();
        System.out.print("Mxy (N-m/m): "); double mxy = sc.nextDouble();

        load inputLoad = new load(nx, ny, nxy, mx, my, mxy);

        // Step 7: Run stress-strain calculation
        PlyStressStrain result = calculateStressStrain(
                numLayers, ABD_inv, QbarAll, lowerBounds, upperBounds, inputLoad);

        // Step 8: Output results
        System.out.println("\n--- Midplane Strains & Curvatures ---");
        double[][] forceVector = {
                {inputLoad.nx}, {inputLoad.ny}, {inputLoad.nxy},
                {inputLoad.mx}, {inputLoad.my}, {inputLoad.mxy}
        };
        double[][] strainCurvature = multiplyMatrices(ABD_inv, forceVector);
        for (int i = 0; i < 3; i++) {
            System.out.printf("ε%d = %.6e%n", i + 1, strainCurvature[i][0]);
        }
        for (int i = 3; i < 6; i++) {
            System.out.printf("κ%d = %.6e%n", i - 2, strainCurvature[i][0]);
        }

        for (int i = 0; i < numLayers; i++) {
            System.out.printf("\n--- Ply %d ---\n", i + 1);
            System.out.printf("Strain Lower: ε1=%.3e, ε2=%.3e, γ12=%.3e%n",
                    result.strainLower[i][0], result.strainLower[i][1], result.strainLower[i][2]);
            System.out.printf("Strain Mid:   ε1=%.3e, ε2=%.3e, γ12=%.3e%n",
                    result.strainMid[i][0], result.strainMid[i][1], result.strainMid[i][2]);
            System.out.printf("Strain Upper: ε1=%.3e, ε2=%.3e, γ12=%.3e%n",
                    result.strainUpper[i][0], result.strainUpper[i][1], result.strainUpper[i][2]);

            System.out.printf("Stress Lower: σ1=%.3e, σ2=%.3e, τ12=%.3e%n",
                    result.stressLower[i][0], result.stressLower[i][1], result.stressLower[i][2]);
            System.out.printf("Stress Mid:   σ1=%.3e, σ2=%.3e, τ12=%.3e%n",
                    result.stressMid[i][0], result.stressMid[i][1], result.stressMid[i][2]);
            System.out.printf("Stress Upper: σ1=%.3e, σ2=%.3e, τ12=%.3e%n",
                    result.stressUpper[i][0], result.stressUpper[i][1], result.stressUpper[i][2]);
        }



        


        sc.close();
    }




    public static void main(String[] args){
        System.out.println("1. Micromechanics");
        System.out.println("2. S, Sbar, Q, Qbar Matrices and Modulus in give direction calculation");
        System.out.println("3. Transform Stress, Strain");
        System.out.println("4. Strength in any given direction");
        System.out.println("5. ABD, ABD inv matrix calculation");
        System.out.println("6. Mid plane strains, curvatures and stress strain values for top, middle, bottom surfaces of ply");
        Scanner sc= new Scanner(System.in);
        while(true){
            int b= sc.nextInt();
            if(b==1){
                Micromechanics.answer();
            }
            else if(b==2){
                ans();
            }
            else if(b==3){
                answ();
            }
            else if(b==4){
                printPlyStrengthFromUser();
            }
            else if(b==5){
                computeABDwithMappedMaterials();
            }
            else if(b==6){
                fullLaminateAnalysis();
            }
            else{
                break;
            }
        }
    }


    // public static void main(String[] args) {
    // // Step 1: Material Definition
    //     Macromechanics.MaterialProperties carbonEpoxy = new Macromechanics.MaterialProperties(
    //         1, 181e9, 10.3e9, 7.17e9, 0.28,
    //         1e9, 1e9, 1e9, 1e9, 1e9
    //     );

    //     Macromechanics.MaterialProperties boronEpoxy = new Macromechanics.MaterialProperties(
    //         2, 204e9, 18.5e9, 5.59e9, 0.23,
    //         1e9, 1e9, 1e9, 1e9, 1e9
    //     );

    //     Macromechanics.MaterialProperties glassEpoxy = new Macromechanics.MaterialProperties(
    //         3, 38.6e9, 8.27e9, 4.14e9, 0.26,
    //         1e9, 1e9, 1e9, 1e9, 1e9
    //     );

        

    //     // Step 2: Layup Definition
    //     List<Macromechanics.layup> layups = new ArrayList<>();
    //     layups.add(new Macromechanics.layup(1, 30, 0.006));
    //     layups.add(new Macromechanics.layup(1, 45, 0.006));
    //     layups.add(new Macromechanics.layup(2, 0, 0.005));
    //     layups.add(new Macromechanics.layup(2, 60, 0.005));
    //     layups.add(new Macromechanics.layup(3, -45, 0.008));
    //     layups.add(new Macromechanics.layup(3, 45, 0.008));

    //     // Step 3: Material Map
    //     Map<Integer, Macromechanics.MaterialProperties> matMap = new HashMap<>();
    //     matMap.put(1, carbonEpoxy);
    //     matMap.put(2, boronEpoxy);
    //     matMap.put(3, glassEpoxy);

    //     // Step 4: Compute A, B, and D matrices
    //     double[][] A = Macromechanics.computeAMatrix(layups, matMap);
    //     double[][] B = Macromechanics.computeBMatrix(layups, matMap);
    //     double[][] D = Macromechanics.computeDMatrix(layups, matMap);

    //     // Step 6: Compute apparent stiffness properties
    //     LaminateStiffnessProperties props = Macromechanics.computeApparentStiffness(A, D, layups);

    //     System.out.println("Apparent Laminate Stiffness Properties:");
    //     System.out.println("----------------------------------------");
    //     System.out.printf("Ex  (in-plane modulus in x)      : %12.3e Pa\n", props.Ex);
    //     System.out.printf("Ey  (in-plane modulus in y)      : %12.3e Pa\n", props.Ey);
    //     System.out.printf("Gxy (in-plane shear modulus)     : %12.3e Pa\n", props.Gxy);
    //     System.out.printf("Exb (bending modulus in x)       : %12.3e Pa\n", props.Exb);
    //     System.out.printf("Eyb (bending modulus in y)       : %12.3e Pa\n", props.Eyb);

    // }


}

