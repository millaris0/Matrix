//#include "matrix_class.cpp"
//#include <cmath>
//
//double* jacobi(double** a, double* y, int n){
//    double* res = new double[n];
//    int i, j;
//    for(i = 0; i<n; i++){
//        res[i] = y[i] / a[i][i];
//    }
//    double sum = 0, eps = 0.0001;
//    double* Xn = new double[n];
//
//    do{
//        for(i=0; i<n; i++){
//            Xn[i] = y[i] / a[i][i];
//            for(j=0; j<n; j++){
//                if(i == j)
//                    continue;
//                else{
//                    Xn[i] -= a[i][j] / a[i][i] * res[j];
//                }
//            }
//        }
//        bool flag = true;
//        for(i=0; i<n; i++){
//            if(abs(Xn[i] - res[i]) > eps){
//                flag = false;
//                break;
//            }
//        }
//
//        for(i=0; i<n; i++){
//            res[i] = Xn[i];
//        }
//        if (flag)
//            break;
//
//    } while (1);
//    return res;
//}