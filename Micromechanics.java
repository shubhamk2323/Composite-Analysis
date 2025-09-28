import java.util.*;

public class Micromechanics{

    public static class FiberInfo{
        double E1, E2, G12, v12, v21, Vf, Xf;
    
        public FiberInfo(double E1, double E2, double G12, double v12, double Vf, double Xf) {
            this.E1 = E1;
            this.E2 = E2;
            this.G12 = G12;
            this.v12 = v12;
            this.Vf = Vf;
            this.v21 = (E2 / E1) * v12;
            this.Xf= Xf;
        }
    
    }
    
    public static class MatrixInfo{
    
        double Em, Gm, Vm, Xm;
    
        public MatrixInfo(double Em, double Gm, double Vm, double Xm){
            this.Em= Em;
            this.Gm= Gm;
            this.Vm= Vm;
            this.Xm= Xm;
        }
    
    }

    public static class CompositeMicroMech{

        double Ec1, Ec2, Gc12, vc12, vc21, Xt;

        public CompositeMicroMech(double Ec1, double Ec2, double Gc12, double vc12, double vc21, double Xt){
            this.Ec1= Ec1;
            this.Ec2= Ec2;
            this.Gc12= Gc12;
            this.vc12= vc12;
            this.vc21 = vc21;
            this.Xt= Xt;
        }

    }

    public static CompositeMicroMech CompositeInfoCalc(FiberInfo x, MatrixInfo y){
        double Ec1 = (x.E1 * x.Vf) + (y.Em * (1 - x.Vf));
        double Ec2 = 1 / ((x.Vf / x.E2) + ((1 - x.Vf) / y.Em));
        double Gc12 = 1 / ((x.Vf / x.G12) + ((1 - x.Vf) / y.Gm));
        double vc12 = (x.v12 * x.Vf) + (y.Vm * (1 - x.Vf));
        double vc21 = (Ec2 / Ec1) * vc12;
        double Xt;
        double ef= x.Xf/x.E1;
        double em= y.Xm/y.Em;
        if(ef<em){
            Xt= ((1-x.Vf)*(y.Xm)) + ((x.Vf)*(x.E1)*ef);
            // Xt= (1-x.Vf)*(y.Xm);
        }
        else{
            Xt= ((1-x.Vf)*(y.Xm)) + ((x.Vf)*(x.E1)*ef);
            //Xt= (1-x.Vf)*(y.Xm);
        }

        return new CompositeMicroMech(Ec1, Ec2, Gc12, vc12, vc21, Xt);

    }


    public static void answer() {
        Scanner sc = new Scanner(System.in);

        System.out.println("Enter Fiber Properties:");
        System.out.print("E1 (GPa): ");
        double E1 = sc.nextDouble();
        System.out.print("E2 (GPa): ");
        double E2 = sc.nextDouble();
        System.out.print("G12 (GPa): ");
        double G12 = sc.nextDouble();
        System.out.print("v12: ");
        double v12 = sc.nextDouble();
        System.out.print("Vf (volume fraction): ");
        double Vf = sc.nextDouble();
        System.out.print("Xf (tensile strength MPa): ");
        double Xf = sc.nextDouble();

        System.out.println("\nEnter Matrix Properties:");
        System.out.print("Em (GPa): ");
        double Em = sc.nextDouble();
        System.out.print("Gm (GPa): ");
        double Gm = sc.nextDouble();
        System.out.print("Vm (Poisson ratio): ");
        double Vm = sc.nextDouble();
        System.out.print("Xm (tensile strength MPa): ");
        double Xm = sc.nextDouble();

        FiberInfo fiber = new FiberInfo(E1, E2, G12, v12, Vf, Xf);
        MatrixInfo matrix = new MatrixInfo(Em, Gm, Vm, Xm);

        CompositeMicroMech composite = CompositeInfoCalc(fiber, matrix);

        System.out.println("\n----- Composite Properties (Micromechanics) -----");
        System.out.printf("E_c1 (Longitudinal Modulus) = %.3e GPa\n", composite.Ec1);
        System.out.printf("E_c2 (Transverse Modulus) = %.3e GPa\n", composite.Ec2);
        System.out.printf("G_c12 (Shear Modulus) = %.3e GPa\n", composite.Gc12);
        System.out.printf("v_c12 (Poisson Ratio) = %.3e", composite.vc12);
        System.out.printf("v_c21 (Poisson Ratio) = %.3e", composite.vc21);
        System.out.printf("Xt (Tensile Strength) = %.3e MPa\n", composite.Xt);
    }



    public static void main(String[] args) {
        answer();
    }

}