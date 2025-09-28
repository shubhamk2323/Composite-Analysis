
import java.util.*;

public class Dataclasses{

    public class FiberInfo{
        double E1;
        double E2;
        double G12;
        double v12;
        double v21;
        double Vf;
        double Xf;
    
        public FiberInfo(double E1, double E2, 
        double G12, double v12, double Vf, double Xf) {
            this.E1 = E1;
            this.E2 = E2;
            this.G12 = G12;
            this.v12 = v12;
            this.Vf = Vf;
            this.v21 = (E2 / E1) * v12;
            this.Xf= Xf;
        }
    
    }
    
    public class MatrixInfo{
    
        double Em;
        double Gm;
        double Vm;
        double Xm
    
        public MatrixInfo(double Em, double Gm, double Vm, double Xm){
            this.Em= Em;
            this.Gm= Gm;
            this.Vm= Vm;
            this.Xm= Xm
        }
    
    }

    public class CompositeMicroMech{

        double Ec1;
        double Ec2;
        double Gc12;
        double vc12;
        double vc21;
        double Xt;

        public CompositeMicroMech(double Ec1, double Ec2, double Gc12,
        double vc12, double Xt){
            this.Ec1= Ec1;
            this.Ec2= Ec2;
            this.Gc12= Gc12;
            this.vc12= vc12;
            this.v21 = (E2 / E1) * v12;
            this.Xt= Xt;
        }

    }
    
}   