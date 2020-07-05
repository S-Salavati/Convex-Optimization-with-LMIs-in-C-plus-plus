#include <iostream>
#include <iomanip>
#include <cmath>
#include <tuple>
#include "fusion.h"

using namespace mosek::fusion;
using namespace monty;
using namespace std;

void System_Dynamics (double rho_1, double rho_2, double delta_1, double W_e, double W_u, Matrix::t &A_p, Matrix::t &B_p, Matrix::t &D_p, Matrix::t &C_y, Matrix::t &C_z, Matrix::t &D_zu, Matrix::t &D_zd)
{
    A_p = Matrix::dense(new_array_ptr<double,2>({ {-1/rho_1,0.0}, {-1.0,-delta_1} }));
    B_p = Matrix::dense(new_array_ptr<double,2>({ {rho_2/rho_1}, {0.0} }));
    D_p = Matrix::dense(new_array_ptr<double,2>({ {0.0,0.01}, {1.0,0.0} }));
    C_y = Matrix::dense(new_array_ptr<double,2>({ {1.0,0.0}, {0.0,1.0} }));
    C_z = Matrix::dense(new_array_ptr<double,2>({ {0.0,W_e}, {0.0,0.0} }));
    D_zu = Matrix::dense(new_array_ptr<double,2>({ {0.0}, {W_u} }));
    D_zd = Matrix::dense(new_array_ptr<double,2>({ {0.0,0.0}, {0.0,0.0} }));
 //   return {A_p,B_p,D_p,C_y,C_z,D_zu,D_zd};
}

void Input_Dist_Dynamics (Matrix::t &W, Matrix::t &V, Matrix::t &H, Matrix::t &J)
{
    W = Matrix::dense(new_array_ptr<double,2>({ {0.,0.0035*M_PI}, {-0.0035*M_PI,0} }));
    V = Matrix::dense(new_array_ptr<double,2>({ {1000.0,100.0} }));
    H = Matrix::dense(new_array_ptr<double,2>({ {0.01}, {0.01} }));
    J = Matrix::dense(new_array_ptr<double,2>({ {0.0} }));
   // J=0.0;
 //   return {W,H,V,J};
}

void Constant_Mats (int n_p, int n_d, double c_w_1, double c_w_2, Matrix::t &C_ew, Matrix::t &E, Matrix::t &G_A, Matrix::t &G_B, Matrix::t &G_D)
{
    C_ew = Matrix::dense(new_array_ptr<double,2>({ {c_w_1,0.0},{0.0, c_w_2} }));
    E = Matrix::dense(new_array_ptr<double,2>({ {0.02,0.0},{0.0, 0.0} }));
    G_A = Matrix::dense(new_array_ptr<double,2>({ {0.005,0.0},{0.0, 0.0} }));
    G_B = Matrix::dense(new_array_ptr<double,2>({ {0.005},{0.0} }));        
    G_D = Matrix::dense(new_array_ptr<double,2>({ {0.0, 0.0},{0.0, 0.0} }));    
//    return {C_ew, E, G_A, G_B, G_D};
}


int main(int argc, char ** argv)
{
    Model::t M = new Model("LPV_sat"); auto _M = finally([&]() { M->dispose(); });
    M->setLogHandler([](const string & msg) { cout << msg << flush; });


    double rho_1_vec[] = {100.0, 150.0, 200.0};
    double rho_2_vec[] = {0.45, 0.55, 0.66};
    double kappa=50.0, epsilon=5.0, u_bar=120.0;
    double delta_1=1.0e-8, W_e=0.1, W_u=0.3, delta=4.0e-5, c_w_1=0.01, c_w_2=0.01, inequality = 1.0e-8;
    double v_1_val = 0.06 , v_2_val = 2.0e-4, v_3_val = 0.5, theta_max=8.0, theta_prime = 1.0;
//     cout<<v_2<<'\n';

    Matrix::t A_p,B_p,D_p,C_y,C_z,D_zu,D_zd, W, V, H, J, C_ew, E, G_A, G_B, G_D;
  //  double J;
    System_Dynamics (rho_1_vec[0], rho_2_vec[0], delta_1, W_e, W_u, A_p,B_p,D_p,C_y,C_z,D_zu,D_zd);

    Input_Dist_Dynamics(W, V, H, J);

    int n_p = A_p->numColumns();
    int n_u = B_p->numColumns();
    int n_y = C_y->numRows();
    int n_d = D_p->numColumns();
    int n_omega = W->numRows();
    int n_d_i = V->numRows();
    int n_nu = H->numColumns();                                    

    auto V_bar= Expr::hstack(Expr::constTerm(Matrix::sparse(n_u,n_p)), Expr::constTerm(Matrix::sparse(n_u,n_p)), Expr::constTerm(V));
    auto J_bar=Expr::hstack(Expr::constTerm(Matrix::sparse(n_u,n_p)), Expr::constTerm(J)); 
  
    Constant_Mats (n_p, n_d, c_w_1, c_w_2,C_ew, E, G_A, G_B, G_D);
    int n_C_ew = C_ew->numRows();                                                         

    auto P_0 = M->variable( new_array_ptr<int,1>({2*n_p+n_omega,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto P_1 = M->variable( new_array_ptr<int,1>({2*n_p+n_omega,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto P_2 = M->variable( new_array_ptr<int,1>({2*n_p+n_omega,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto P_3 = M->variable( new_array_ptr<int,1>({2*n_p+n_omega,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto P_4 = M->variable( new_array_ptr<int,1>({2*n_p+n_omega,2*n_p+n_omega}), Domain::unbounded() ) ;

    auto Q_0 = M->variable( new_array_ptr<int,1>({2*n_p+n_omega,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto Q_1 = M->variable( new_array_ptr<int,1>({2*n_p+n_omega,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto Q_2 = M->variable( new_array_ptr<int,1>({2*n_p+n_omega,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto Q_3 = M->variable( new_array_ptr<int,1>({2*n_p+n_omega,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto Q_4 = M->variable( new_array_ptr<int,1>({2*n_p+n_omega,2*n_p+n_omega}), Domain::unbounded() ) ;

    auto S_0 = M->variable( new_array_ptr<int,1>({2*n_p+n_omega,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto S_1 = M->variable( new_array_ptr<int,1>({2*n_p+n_omega,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto S_2 = M->variable( new_array_ptr<int,1>({2*n_p+n_omega,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto S_3 = M->variable( new_array_ptr<int,1>({2*n_p+n_omega,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto S_4 = M->variable( new_array_ptr<int,1>({2*n_p+n_omega,2*n_p+n_omega}), Domain::unbounded() ) ;
    
    auto R = M->variable(Domain::inPSDCone(2*n_p+n_omega)) ;        
    auto S1 = M->variable( new_array_ptr<int,1>({2*n_p+n_omega,2*n_p+n_omega}), Domain::unbounded() ) ;

    auto gamma_sq = M->variable("gamma_sq", Domain::greaterThan(0.));
    auto beta = M->variable("beta", Domain::greaterThan(0.));
        
    auto X_0 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto X_1 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto X_2 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto X_3 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto X_4 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
           
    auto Y_0 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto Y_1 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto Y_2 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto Y_3 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto Y_4 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
        
    auto Z_0 = M->variable( new_array_ptr<int,1>({n_omega,n_omega}), Domain::unbounded() ) ;
    auto Z_1 = M->variable( new_array_ptr<int,1>({n_omega,n_omega}), Domain::unbounded() ) ;
    auto Z_2 = M->variable( new_array_ptr<int,1>({n_omega,n_omega}), Domain::unbounded() ) ;
    auto Z_3 = M->variable( new_array_ptr<int,1>({n_omega,n_omega}), Domain::unbounded() ) ;
    auto Z_4 = M->variable( new_array_ptr<int,1>({n_omega,n_omega}), Domain::unbounded() ) ;   
        
    auto A_K_tilde_0 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto A_K_tilde_1 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto A_K_tilde_2 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto A_K_tilde_3 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto A_K_tilde_4 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
        
    auto B_K_tilde_0 = M->variable( new_array_ptr<int,1>({n_p,n_y}), Domain::unbounded() ) ;
    auto B_K_tilde_1 = M->variable( new_array_ptr<int,1>({n_p,n_y}), Domain::unbounded() ) ;
    auto B_K_tilde_2 = M->variable( new_array_ptr<int,1>({n_p,n_y}), Domain::unbounded() ) ;
    auto B_K_tilde_3 = M->variable( new_array_ptr<int,1>({n_p,n_y}), Domain::unbounded() ) ;
    auto B_K_tilde_4 = M->variable( new_array_ptr<int,1>({n_p,n_y}), Domain::unbounded() ) ;

    auto L_K_tilde_0 = M->variable( new_array_ptr<int,1>({n_omega,n_p}), Domain::unbounded() ) ;
    auto L_K_tilde_1 = M->variable( new_array_ptr<int,1>({n_omega,n_p}), Domain::unbounded() ) ;
    auto L_K_tilde_2 = M->variable( new_array_ptr<int,1>({n_omega,n_p}), Domain::unbounded() ) ;
    auto L_K_tilde_3 = M->variable( new_array_ptr<int,1>({n_omega,n_p}), Domain::unbounded() ) ;
    auto L_K_tilde_4 = M->variable( new_array_ptr<int,1>({n_omega,n_p}), Domain::unbounded() ) ;

    auto L_d_tilde_0 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;
    auto L_d_tilde_1 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;
    auto L_d_tilde_2 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;
    auto L_d_tilde_3 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;
    auto L_d_tilde_4 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;
             
    auto L_y_tilde_0 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;
    auto L_y_tilde_1 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;
    auto L_y_tilde_2 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;
    auto L_y_tilde_3 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;
    auto L_y_tilde_4 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;

    auto C_theta_tilde_theta_0 = M->variable( new_array_ptr<int,1>({n_u,n_p}), Domain::unbounded() ) ;
    auto C_theta_tilde_theta_1 = M->variable( new_array_ptr<int,1>({n_u,n_p}), Domain::unbounded() ) ;
    auto C_theta_tilde_theta_2 = M->variable( new_array_ptr<int,1>({n_u,n_p}), Domain::unbounded() ) ;
    auto C_theta_tilde_theta_3 = M->variable( new_array_ptr<int,1>({n_u,n_p}), Domain::unbounded() ) ;
    auto C_theta_tilde_theta_4 = M->variable( new_array_ptr<int,1>({n_u,n_p}), Domain::unbounded() ) ;
    auto C_theta_tilde_theta_5 = M->variable( new_array_ptr<int,1>({n_u,n_p}), Domain::unbounded() ) ;
    auto C_theta_tilde_theta_6 = M->variable( new_array_ptr<int,1>({n_u,n_p}), Domain::unbounded() ) ;

    auto A_theta_tilde_theta_0 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto A_theta_tilde_theta_1 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto A_theta_tilde_theta_2 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto A_theta_tilde_theta_3 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto A_theta_tilde_theta_4 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto A_theta_tilde_theta_5 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;
    auto A_theta_tilde_theta_6 = M->variable( new_array_ptr<int,1>({n_p,n_p}), Domain::unbounded() ) ;

    auto B_theta_tilde_theta_0 = M->variable( new_array_ptr<int,1>({n_p,n_y}), Domain::unbounded() ) ;
    auto B_theta_tilde_theta_1 = M->variable( new_array_ptr<int,1>({n_p,n_y}), Domain::unbounded() ) ;
    auto B_theta_tilde_theta_2 = M->variable( new_array_ptr<int,1>({n_p,n_y}), Domain::unbounded() ) ;
    auto B_theta_tilde_theta_3 = M->variable( new_array_ptr<int,1>({n_p,n_y}), Domain::unbounded() ) ;
    auto B_theta_tilde_theta_4 = M->variable( new_array_ptr<int,1>({n_p,n_y}), Domain::unbounded() ) ;
    auto B_theta_tilde_theta_5 = M->variable( new_array_ptr<int,1>({n_p,n_y}), Domain::unbounded() ) ;
    auto B_theta_tilde_theta_6 = M->variable( new_array_ptr<int,1>({n_p,n_y}), Domain::unbounded() ) ;                

    auto D_K_theta_0 = M->variable( new_array_ptr<int,1>({n_u,n_y}), Domain::unbounded() ) ;
    auto D_K_theta_1 = M->variable( new_array_ptr<int,1>({n_u,n_y}), Domain::unbounded() ) ;
    auto D_K_theta_2 = M->variable( new_array_ptr<int,1>({n_u,n_y}), Domain::unbounded() ) ;
    auto D_K_theta_3 = M->variable( new_array_ptr<int,1>({n_u,n_y}), Domain::unbounded() ) ;
    auto D_K_theta_4 = M->variable( new_array_ptr<int,1>({n_u,n_y}), Domain::unbounded() ) ;
    auto D_K_theta_5 = M->variable( new_array_ptr<int,1>({n_u,n_y}), Domain::unbounded() ) ;
    auto D_K_theta_6 = M->variable( new_array_ptr<int,1>({n_u,n_y}), Domain::unbounded() ) ;   

    auto L_K_theta_tilde_theta_0 = M->variable( new_array_ptr<int,1>({n_omega,n_p}), Domain::unbounded() ) ;
    auto L_K_theta_tilde_theta_1 = M->variable( new_array_ptr<int,1>({n_omega,n_p}), Domain::unbounded() ) ;
    auto L_K_theta_tilde_theta_2 = M->variable( new_array_ptr<int,1>({n_omega,n_p}), Domain::unbounded() ) ;
    auto L_K_theta_tilde_theta_3 = M->variable( new_array_ptr<int,1>({n_omega,n_p}), Domain::unbounded() ) ;
    auto L_K_theta_tilde_theta_4 = M->variable( new_array_ptr<int,1>({n_omega,n_p}), Domain::unbounded() ) ;
    auto L_K_theta_tilde_theta_5 = M->variable( new_array_ptr<int,1>({n_omega,n_p}), Domain::unbounded() ) ;
    auto L_K_theta_tilde_theta_6 = M->variable( new_array_ptr<int,1>({n_omega,n_p}), Domain::unbounded() ) ;

    auto L_y_theta_tilde_theta_0 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;
    auto L_y_theta_tilde_theta_1 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;
    auto L_y_theta_tilde_theta_2 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;
    auto L_y_theta_tilde_theta_3 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;
    auto L_y_theta_tilde_theta_4 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;
    auto L_y_theta_tilde_theta_5 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;
    auto L_y_theta_tilde_theta_6 = M->variable( new_array_ptr<int,1>({n_omega,n_y}), Domain::unbounded() ) ;        
        
    auto E_K_tilde_theta_0 = M->variable( new_array_ptr<int,1>({n_p,n_u}), Domain::unbounded() ) ;
    auto E_K_tilde_theta_1 = M->variable( new_array_ptr<int,1>({n_p,n_u}), Domain::unbounded() ) ;
    auto E_K_tilde_theta_2 = M->variable( new_array_ptr<int,1>({n_p,n_u}), Domain::unbounded() ) ;
    auto E_K_tilde_theta_3 = M->variable( new_array_ptr<int,1>({n_p,n_u}), Domain::unbounded() ) ;
    auto E_K_tilde_theta_4 = M->variable( new_array_ptr<int,1>({n_p,n_u}), Domain::unbounded() ) ;
    auto E_K_tilde_theta_5 = M->variable( new_array_ptr<int,1>({n_p,n_u}), Domain::unbounded() ) ;
    auto E_K_tilde_theta_6 = M->variable( new_array_ptr<int,1>({n_p,n_u}), Domain::unbounded() ) ;  

    auto F_K_tilde_theta_0 = M->variable( new_array_ptr<int,1>({n_omega,n_u}), Domain::unbounded() ) ;
    auto F_K_tilde_theta_1 = M->variable( new_array_ptr<int,1>({n_omega,n_u}), Domain::unbounded() ) ;
    auto F_K_tilde_theta_2 = M->variable( new_array_ptr<int,1>({n_omega,n_u}), Domain::unbounded() ) ;
    auto F_K_tilde_theta_3 = M->variable( new_array_ptr<int,1>({n_omega,n_u}), Domain::unbounded() ) ;
    auto F_K_tilde_theta_4 = M->variable( new_array_ptr<int,1>({n_omega,n_u}), Domain::unbounded() ) ;
    auto F_K_tilde_theta_5 = M->variable( new_array_ptr<int,1>({n_omega,n_u}), Domain::unbounded() ) ;
    auto F_K_tilde_theta_6 = M->variable( new_array_ptr<int,1>({n_omega,n_u}), Domain::unbounded() ) ;  

    auto t_bar_theta_0 = M->variable("t_bar_theta_0", n_u, Domain::unbounded());
    auto t_bar_theta_1 = M->variable("t_bar_theta_1", n_u, Domain::unbounded()); 
    auto t_bar_theta_2 = M->variable("t_bar_theta_2", n_u, Domain::unbounded());
    auto t_bar_theta_3 = M->variable("t_bar_theta_3", n_u, Domain::unbounded());
    auto t_bar_theta_4 = M->variable("t_bar_theta_4", n_u, Domain::unbounded());
    auto t_bar_theta_5 = M->variable("t_bar_theta_5", n_u, Domain::unbounded());
    auto t_bar_theta_6 = M->variable("t_bar_theta_6", n_u, Domain::unbounded());

//     auto x = M->variable(n);

    auto T_bar_theta_0 = Expr::add(
                            make_shared<ndarray<Expression::t,1>>(
                            shape(n_u), 
                            [t_bar_theta_0,n_u](int i) { return Expr::mul( t_bar_theta_0->index(i), 
                                 Matrix::sparse(n_u, n_u, 
                                                new_array_ptr<int, 1>({i}), 
                                                new_array_ptr<int, 1>({i}), 
                                                new_array_ptr<double, 1>({1.0})) 
                             );}));

    auto T_bar_theta_1 = Expr::add(
                            make_shared<ndarray<Expression::t,1>>(
                            shape(n_u), 
                            [t_bar_theta_1,n_u](int i) { return Expr::mul( t_bar_theta_1->index(i), 
                                 Matrix::sparse(n_u, n_u, 
                                                new_array_ptr<int, 1>({i}), 
                                                new_array_ptr<int, 1>({i}), 
                                                new_array_ptr<double, 1>({1.0})) 
                             );}));
                             
    auto T_bar_theta_2 = Expr::add(
                            make_shared<ndarray<Expression::t,1>>(
                            shape(n_u), 
                            [t_bar_theta_2,n_u](int i) { return Expr::mul( t_bar_theta_2->index(i), 
                                 Matrix::sparse(n_u, n_u, 
                                                new_array_ptr<int, 1>({i}), 
                                                new_array_ptr<int, 1>({i}), 
                                                new_array_ptr<double, 1>({1.0})) 
                             );}));
                             
    auto T_bar_theta_3 = Expr::add(
                            make_shared<ndarray<Expression::t,1>>(
                            shape(n_u), 
                            [t_bar_theta_3,n_u](int i) { return Expr::mul( t_bar_theta_3->index(i), 
                                 Matrix::sparse(n_u, n_u, 
                                                new_array_ptr<int, 1>({i}), 
                                                new_array_ptr<int, 1>({i}), 
                                                new_array_ptr<double, 1>({1.0})) 
                             );}));
                             
    auto T_bar_theta_4 = Expr::add(
                            make_shared<ndarray<Expression::t,1>>(
                            shape(n_u), 
                            [t_bar_theta_4,n_u](int i) { return Expr::mul( t_bar_theta_4->index(i), 
                                 Matrix::sparse(n_u, n_u, 
                                                new_array_ptr<int, 1>({i}), 
                                                new_array_ptr<int, 1>({i}), 
                                                new_array_ptr<double, 1>({1.0})) 
                             );}));
                             
                             
    auto T_bar_theta_5 = Expr::add(
                            make_shared<ndarray<Expression::t,1>>(
                            shape(n_u), 
                            [t_bar_theta_5,n_u](int i) { return Expr::mul( t_bar_theta_5->index(i), 
                                 Matrix::sparse(n_u, n_u, 
                                                new_array_ptr<int, 1>({i}), 
                                                new_array_ptr<int, 1>({i}), 
                                                new_array_ptr<double, 1>({1.0})) 
                             );}));
                             
    auto T_bar_theta_6 = Expr::add(
                            make_shared<ndarray<Expression::t,1>>(
                            shape(n_u), 
                            [t_bar_theta_6,n_u](int i) { return Expr::mul( t_bar_theta_6->index(i), 
                                 Matrix::sparse(n_u, n_u, 
                                                new_array_ptr<int, 1>({i}), 
                                                new_array_ptr<int, 1>({i}), 
                                                new_array_ptr<double, 1>({1.0})) 
                             );}));

//            auto T_bar_theta_0 = M->variable( new_array_ptr<int,1>({n_u,n_u}), Domain::unbounded() ) ;
//            Matrix::t T_bar_theta_0 = Matrix::diag(new_array_ptr<int,1>({n_u,n_u}) t_bar_theta_0);
//            auto T_bar_theta_1 = M->variable( new_array_ptr<int,1>({n_u,n_u}), Domain::unbounded() ) ;
//           auto T_bar_theta_2 = M->variable( new_array_ptr<int,1>({n_u,n_u}), Domain::unbounded() ) ;
//          auto T_bar_theta_3 = M->variable( new_array_ptr<int,1>({n_u,n_u}), Domain::unbounded() ) ;
//         auto T_bar_theta_4 = M->variable( new_array_ptr<int,1>({n_u,n_u}), Domain::unbounded() ) ;
//         auto T_bar_theta_5 = M->variable( new_array_ptr<int,1>({n_u,n_u}), Domain::unbounded() ) ;
//         auto T_bar_theta_6 = M->variable( new_array_ptr<int,1>({n_u,n_u}), Domain::unbounded() ) ;
        
    auto G_hat_theta_0 = M->variable( new_array_ptr<int,1>({n_u,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto G_hat_theta_1 = M->variable( new_array_ptr<int,1>({n_u,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto G_hat_theta_2 = M->variable( new_array_ptr<int,1>({n_u,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto G_hat_theta_3 = M->variable( new_array_ptr<int,1>({n_u,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto G_hat_theta_4 = M->variable( new_array_ptr<int,1>({n_u,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto G_hat_theta_5 = M->variable( new_array_ptr<int,1>({n_u,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto G_hat_theta_6 = M->variable( new_array_ptr<int,1>({n_u,2*n_p+n_omega}), Domain::unbounded() ) ;  

    auto G_1_hat_theta_0 = M->variable( new_array_ptr<int,1>({n_u,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto G_1_hat_theta_1 = M->variable( new_array_ptr<int,1>({n_u,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto G_1_hat_theta_2 = M->variable( new_array_ptr<int,1>({n_u,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto G_1_hat_theta_3 = M->variable( new_array_ptr<int,1>({n_u,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto G_1_hat_theta_4 = M->variable( new_array_ptr<int,1>({n_u,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto G_1_hat_theta_5 = M->variable( new_array_ptr<int,1>({n_u,2*n_p+n_omega}), Domain::unbounded() ) ;
    auto G_1_hat_theta_6 = M->variable( new_array_ptr<int,1>({n_u,2*n_p+n_omega}), Domain::unbounded() ) ;  
    
    auto LMI_1 = Var::vstack(Var::hstack(R, S1), Var::hstack(S1->transpose(), R));
    M->constraint(LMI_1, Domain::inPSDCone((*LMI_1->getShape())[0]));
    M->constraint( Expr::sub(delta,beta), Domain::greaterThan(0.));
    for (auto rho_1 : rho_1_vec)
        {
        for (auto rho_2 : rho_2_vec)
            {
            System_Dynamics (rho_1, rho_2, delta_1, W_e, W_u, A_p,B_p,D_p,C_y,C_z,D_zu,D_zd);
           
            auto X= Expr::add(X_0, Expr::add(Expr::mul(rho_1,X_1), Expr::add(Expr::mul(0.5*pow(rho_1,2),X_2), Expr::add(Expr::mul(rho_2,X_3), Expr::mul(0.5*pow(rho_2,2),X_4)))));
            auto Y= Expr::add(Y_0, Expr::add(Expr::mul(rho_1,Y_1), Expr::add(Expr::mul(0.5*pow(rho_1,2),Y_2), Expr::add(Expr::mul(rho_2,Y_3), Expr::mul(0.5*pow(rho_2,2),Y_4)))));
            auto Z= Expr::add(Z_0, Expr::add(Expr::mul(rho_1,Z_1), Expr::add(Expr::mul(0.5*pow(rho_1,2),Z_2), Expr::add(Expr::mul(rho_2,Z_3), Expr::mul(0.5*pow(rho_2,2),Z_4)))));
            M->constraint(Expr::sub(Expr::transpose(X),X), Domain::equalsTo(0.0));
            M->constraint(Expr::sub(Expr::transpose(Y),Y), Domain::equalsTo(0.0));
            M->constraint(Expr::sub(Expr::transpose(Z),Z), Domain::equalsTo(0.0));
                            
            auto L_d_tilde= Expr::add(L_d_tilde_0, Expr::add(Expr::mul(rho_1,L_d_tilde_1), Expr::add(Expr::mul(0.5*pow(rho_1,2),L_d_tilde_2), Expr::add(Expr::mul(rho_2,L_d_tilde_3), Expr::mul(0.5*pow(rho_2,2),L_d_tilde_4)))));      
            auto A_K_tilde= Expr::add(A_K_tilde_0, Expr::add(Expr::mul(rho_1,A_K_tilde_1), Expr::add(Expr::mul(0.5*pow(rho_1,2),A_K_tilde_2), Expr::add(Expr::mul(rho_2,A_K_tilde_3), Expr::mul(0.5*pow(rho_2,2),A_K_tilde_4)))));    
            auto B_K_tilde= Expr::add(B_K_tilde_0, Expr::add(Expr::mul(rho_1,B_K_tilde_1), Expr::add(Expr::mul(0.5*pow(rho_1,2),B_K_tilde_2), Expr::add(Expr::mul(rho_2,B_K_tilde_3), 
                 Expr::mul(0.5*pow(rho_2,2),B_K_tilde_4)))));            
            auto L_K_tilde= Expr::add(L_K_tilde_0, Expr::add(Expr::mul(rho_1,L_K_tilde_1), Expr::add(Expr::mul(0.5*pow(rho_1,2),L_K_tilde_2), Expr::add(Expr::mul(rho_2,L_K_tilde_3), Expr::mul(0.5*pow(rho_2,2),L_K_tilde_4)))));    
            auto L_y_tilde= Expr::add(L_y_tilde_0, Expr::add(Expr::mul(rho_1,L_y_tilde_1), Expr::add(Expr::mul(0.5*pow(rho_1,2),L_y_tilde_2), Expr::add(Expr::mul(rho_2,L_y_tilde_3),
                 Expr::mul(0.5*pow(rho_2,2),L_y_tilde_4)))));            
                   
            auto P_bar= Expr::add(P_0, Expr::add(Expr::mul(rho_1,P_1), Expr::add(Expr::mul(0.5*pow(rho_1,2),P_2), Expr::add(Expr::mul(rho_2,P_3), Expr::mul(0.5*pow(rho_2,2),P_4)))));
            auto Q_bar= Expr::add(Q_0, Expr::add(Expr::mul(rho_1,Q_1), Expr::add(Expr::mul(0.5*pow(rho_1,2),Q_2), Expr::add(Expr::mul(rho_2,Q_3), Expr::mul(0.5*pow(rho_2,2),Q_4)))));
            auto S_bar= Expr::add(S_0, Expr::add(Expr::mul(rho_1,S_1), Expr::add(Expr::mul(0.5*pow(rho_1,2),S_2), Expr::add(Expr::mul(rho_2,S_3), Expr::mul(0.5*pow(rho_2,2),S_4)))));
                   
            auto B_w_hat = Expr::vstack(Expr::hstack(Expr::constTerm(D_p), Expr::mul(B_p,Expr::constTerm(J))), Expr::hstack(Expr::mul(X,D_p), Expr::mul(Expr::mul(X,B_p),J)), Expr::hstack(Expr::mul(Expr::mul(L_d_tilde, C_y),D_p),
                 Expr::add(Expr::mul(Z,H), Expr::mul(Expr::mul(Expr::mul(L_d_tilde, C_y),B_p),J))));
            auto C_hat = Expr::vstack(Expr::hstack(Expr::mul(C_z,Y), Expr::constTerm(C_z), Expr::mul(D_zu,Expr::constTerm(V))), Expr::hstack(Expr::constTerm(Matrix::sparse(n_C_ew,n_p)), Expr::constTerm(Matrix::sparse(n_C_ew,n_p)),
                 Expr::constTerm(C_ew)));
            auto D_psi = Expr::vstack(Expr::mul(Expr::constTerm(-1),D_zu), Expr::constTerm(Matrix::sparse(n_C_ew,n_u)));  
            auto D_w = Expr::vstack(Expr::hstack(Expr::constTerm(D_zd), Expr::mul(D_zu, Expr::constTerm(J))), Expr::constTerm(Matrix::sparse(n_C_ew,n_d+n_nu)));
            auto E_hat = Expr::vstack(Expr::hstack(Expr::constTerm(E), Expr::constTerm(Matrix::sparse(n_p,n_p)), Expr::constTerm(Matrix::sparse(n_p,n_p))), Expr::hstack(Expr::mul(X,E), Expr::constTerm(Matrix::sparse(n_p,n_p)),
                 Expr::constTerm(Matrix::sparse(n_p,n_p))), Expr::hstack(Expr::mul(Expr::mul(L_d_tilde, C_y),E), Expr::constTerm(Matrix::sparse(n_omega,n_p)), Expr::constTerm(Matrix::sparse(n_omega,n_p))));
            M->constraint(Expr::sub(P_bar, Expr::mul(inequality, Expr::constTerm(Matrix::eye((*P_bar->getShape())[0])))), Domain::inPSDCone((*P_bar->getShape())[0] )); 
            M->constraint(Expr::sub(Q_bar, Expr::mul(inequality, Expr::constTerm(Matrix::eye((*Q_bar->getShape())[0])))), Domain::inPSDCone((*Q_bar->getShape())[0] )); 
            M->constraint(Expr::sub(S_bar, Expr::mul(inequality, Expr::constTerm(Matrix::eye((*S_bar->getShape())[0])))), Domain::inPSDCone((*S_bar->getShape())[0] )); 
            
            for (auto rho_1_d : rho_1_vec)
                {
                for (auto rho_2_d : rho_2_vec)
                    {
                    auto L_K_theta_tilde_theta= Expr::add(L_K_theta_tilde_theta_0, Expr::add(Expr::mul(rho_1, L_K_theta_tilde_theta_1), Expr::add(Expr::mul(rho_2, L_K_theta_tilde_theta_2), Expr::add(Expr::mul(rho_1, Expr::mul(rho_2,
                         L_K_theta_tilde_theta_3)), Expr::add(Expr::mul(rho_1_d, L_K_theta_tilde_theta_4), Expr::add(Expr::mul(rho_2_d, L_K_theta_tilde_theta_5), Expr::mul(rho_1_d, Expr::mul(rho_2_d, L_K_theta_tilde_theta_6))))))));
                    
                    auto A_theta_tilde_theta= Expr::add(A_theta_tilde_theta_0, Expr::add(Expr::mul(rho_1, A_theta_tilde_theta_1), Expr::add(Expr::mul(rho_2, A_theta_tilde_theta_2), Expr::add(Expr::mul(rho_1, Expr::mul(rho_2,
                         A_theta_tilde_theta_3)), Expr::add(Expr::mul(rho_1_d, A_theta_tilde_theta_4), Expr::add(Expr::mul(rho_2_d, A_theta_tilde_theta_5), Expr::mul(rho_1_d, Expr::mul(rho_2_d, A_theta_tilde_theta_6))))))));

                    auto B_theta_tilde_theta= Expr::add(B_theta_tilde_theta_0, Expr::add(Expr::mul(rho_1, B_theta_tilde_theta_1), Expr::add(Expr::mul(rho_2, B_theta_tilde_theta_2), Expr::add(Expr::mul(rho_1, Expr::mul(rho_2,
                         B_theta_tilde_theta_3)), Expr::add(Expr::mul(rho_1_d, B_theta_tilde_theta_4), Expr::add(Expr::mul(rho_2_d, B_theta_tilde_theta_5), Expr::mul(rho_1_d, Expr::mul(rho_2_d, B_theta_tilde_theta_6))))))));
                         
                    auto C_theta_tilde_theta= Expr::add(C_theta_tilde_theta_0, Expr::add(Expr::mul(rho_1, C_theta_tilde_theta_1), Expr::add(Expr::mul(rho_2, C_theta_tilde_theta_2), Expr::add(Expr::mul(rho_1, Expr::mul(rho_2,
                         C_theta_tilde_theta_3)), Expr::add(Expr::mul(rho_1_d, C_theta_tilde_theta_4), Expr::add(Expr::mul(rho_2_d, C_theta_tilde_theta_5), Expr::mul(rho_1_d, Expr::mul(rho_2_d, C_theta_tilde_theta_6))))))));
                         
                    auto T_bar_theta= Expr::add(T_bar_theta_0, Expr::add(Expr::mul(rho_1, T_bar_theta_1), Expr::add(Expr::mul(rho_2, T_bar_theta_2), Expr::add(Expr::mul(rho_1, Expr::mul(rho_2, T_bar_theta_3)),
                         Expr::add(Expr::mul(rho_1_d, T_bar_theta_4), Expr::add(Expr::mul(rho_2_d, T_bar_theta_5), Expr::mul(rho_1_d, Expr::mul(rho_2_d, T_bar_theta_6))))))));
                    
                    auto E_K_tilde_theta= Expr::add(E_K_tilde_theta_0, Expr::add(Expr::mul(rho_1, E_K_tilde_theta_1), Expr::add(Expr::mul(rho_2, E_K_tilde_theta_2), Expr::add(Expr::mul(rho_1, Expr::mul(rho_2, E_K_tilde_theta_3)),
                         Expr::add(Expr::mul(rho_1_d, E_K_tilde_theta_4), Expr::add(Expr::mul(rho_2_d, E_K_tilde_theta_5), Expr::mul(rho_1_d, Expr::mul(rho_2_d, E_K_tilde_theta_6))))))));
                        
                    auto F_K_tilde_theta= Expr::add(F_K_tilde_theta_0, Expr::add(Expr::mul(rho_1, F_K_tilde_theta_1), Expr::add(Expr::mul(rho_2, F_K_tilde_theta_2), Expr::add(Expr::mul(rho_1, Expr::mul(rho_2, F_K_tilde_theta_3)),
                         Expr::add(Expr::mul(rho_1_d, F_K_tilde_theta_4), Expr::add(Expr::mul(rho_2_d, F_K_tilde_theta_5), Expr::mul(rho_1_d, Expr::mul(rho_2_d, F_K_tilde_theta_6))))))));
                         
                    auto G_hat_theta= Expr::add(G_hat_theta_0, Expr::add(Expr::mul(rho_1, G_hat_theta_1), Expr::add(Expr::mul(rho_2, G_hat_theta_2), Expr::add(Expr::mul(rho_1, Expr::mul(rho_2, G_hat_theta_3)),
                         Expr::add(Expr::mul(rho_1_d, G_hat_theta_4), Expr::add(Expr::mul(rho_2_d, G_hat_theta_5), Expr::mul(rho_1_d, Expr::mul(rho_2_d, G_hat_theta_6))))))));
                        
                    auto G_1_hat_theta= Expr::add(G_1_hat_theta_0, Expr::add(Expr::mul(rho_1, G_1_hat_theta_1), Expr::add(Expr::mul(rho_2, G_1_hat_theta_2), Expr::add(Expr::mul(rho_1, Expr::mul(rho_2, G_1_hat_theta_3)),
                         Expr::add(Expr::mul(rho_1_d, G_1_hat_theta_4), Expr::add(Expr::mul(rho_2_d, G_1_hat_theta_5), Expr::mul(rho_1_d, Expr::mul(rho_2_d, G_1_hat_theta_6))))))));
                         
                    auto D_K_theta= Expr::add(D_K_theta_0, Expr::add(Expr::mul(rho_1, D_K_theta_1), Expr::add(Expr::mul(rho_2, D_K_theta_2), Expr::add(Expr::mul(rho_1, Expr::mul(rho_2, D_K_theta_3)), Expr::add(Expr::mul(rho_1_d, 
                         D_K_theta_4), Expr::add(Expr::mul(rho_2_d, D_K_theta_5), Expr::mul(rho_1_d, Expr::mul(rho_2_d, D_K_theta_6))))))));
                        
                    auto L_y_theta_tilde_theta= Expr::add(L_y_theta_tilde_theta_0, Expr::add(Expr::mul(rho_1, L_y_theta_tilde_theta_1), Expr::add(Expr::mul(rho_2, L_y_theta_tilde_theta_2), Expr::add(Expr::mul(rho_1, Expr::mul(rho_2,
                         L_y_theta_tilde_theta_3)), Expr::add(Expr::mul(rho_1_d, L_y_theta_tilde_theta_4), Expr::add(Expr::mul(rho_2_d, L_y_theta_tilde_theta_5), Expr::mul(rho_1_d, Expr::mul(rho_2_d, L_y_theta_tilde_theta_6))))))));
                    
                    auto P_bar_delayed = Expr::add(P_0, Expr::add(Expr::mul(rho_1_d, P_1), Expr::add(Expr::mul(0.5*pow(rho_1_d, 2), P_2), Expr::add(Expr::mul(rho_2_d, P_3), Expr::mul(0.5*pow(rho_2_d,2), P_4)))));
                    
                    auto Q_bar_delayed = Expr::add(Q_0, Expr::add(Expr::mul(rho_1_d, Q_1), Expr::add(Expr::mul(0.5*pow(rho_1_d, 2), Q_2), Expr::add(Expr::mul(rho_2_d, Q_3), Expr::mul(0.5*pow(rho_2_d,2), Q_4)))));
                    
                    auto S_h           = Expr::add(S_0, Expr::add(Expr::mul(rho_1_d, S_1), Expr::add(Expr::mul(0.5*pow(rho_1_d, 2), S_2), Expr::add(Expr::mul(rho_2_d, S_3), Expr::mul(0.5*pow(rho_2_d,2), S_4)))));
                 
                    auto A_d_hat = Expr::vstack(Expr::hstack(Expr::mul(B_p,C_theta_tilde_theta), Expr::mul(B_p,Expr::mul(D_K_theta,C_y)), Expr::constTerm(Matrix::sparse(n_p,n_omega))), Expr::hstack(A_theta_tilde_theta,
                         Expr::mul(B_theta_tilde_theta, C_y), Expr::constTerm(Matrix::sparse(n_p,n_omega))), Expr::hstack(L_K_theta_tilde_theta, Expr::mul(L_y_theta_tilde_theta, C_y), Expr::constTerm(Matrix::sparse(n_p,n_omega))));
            
                    auto K_bb_theta = Expr::hstack(C_theta_tilde_theta, Expr::mul(D_K_theta,C_y), Expr::constTerm(Matrix::sparse(n_u,n_omega)));
               
                    auto B_psi_hat_theta = Expr::vstack(Expr::mul(-1,Expr::mul(B_p,T_bar_theta)), E_K_tilde_theta, F_K_tilde_theta);
            
                    auto C_d_hat = Expr::vstack(Expr::hstack(Expr::mul(D_zu,C_theta_tilde_theta), Expr::mul(D_zu,Expr::mul(D_K_theta,C_y)), Expr::constTerm(Matrix::sparse(n_p,n_omega))), Expr::constTerm(Matrix::sparse(n_C_ew,
                         2*n_p+n_omega)));
                         
                    M->constraint(Expr::sub(P_bar_delayed, Expr::mul(inequality, Expr::constTerm(Matrix::eye((*P_bar_delayed->getShape())[0])))), Domain::inPSDCone((*P_bar_delayed->getShape())[0] )); 
                    M->constraint(Expr::sub(Q_bar_delayed, Expr::mul(inequality, Expr::constTerm(Matrix::eye((*Q_bar_delayed->getShape())[0])))), Domain::inPSDCone((*Q_bar_delayed->getShape())[0] )); 
                    M->constraint(Expr::sub(S_h, Expr::mul(inequality, Expr::constTerm(Matrix::eye((*S_h->getShape())[0])))), Domain::inPSDCone((*S_h->getShape())[0] ));
                                  
                    auto L_13 = Expr::add(R,Expr::sub(A_d_hat,S1));
                    auto L_14 = S1;
                    auto L_15 = Expr::add(B_psi_hat_theta, Expr::add(Expr::transpose(G_hat_theta), Expr::sub(Expr::mul(2, Expr::transpose(V_bar)), Expr::transpose(K_bb_theta))));
                    auto L_16 = B_w_hat;
                    auto L_17 = Expr::transpose(C_hat);
                    auto L_18 = E_hat;
                    auto L_19 = Expr::mul(epsilon, Expr::vstack(Expr::hstack(Expr::mul(G_A,Y), Expr::constTerm(G_A), Expr::mul(G_B,Expr::constTerm(V))), Expr::constTerm(Matrix::sparse(2*n_p,2*n_p+n_omega))));
                    
                    auto L_22 = Expr::add(Expr::mul(pow(theta_max,2), R), Expr::mul(-2*kappa, Expr::vstack(Expr::hstack(Y, Expr::constTerm(Matrix::eye(n_p)), Expr::constTerm(Matrix::sparse(n_p, n_omega))), 
                         Expr::hstack(Expr::constTerm(Matrix::eye(n_p)), X, Expr::constTerm(Matrix::sparse(n_p, n_omega))), Expr::hstack(Expr::constTerm(Matrix::sparse(n_omega,2*n_p)), Z))));
                    auto L_23 = Expr::mul(kappa, A_d_hat);
                    auto L_24 = Expr::constTerm(Matrix::sparse(2*n_p+n_omega,2*n_p+n_omega));
                    auto L_25 = Expr::mul(kappa, B_psi_hat_theta);
                    auto L_26 = Expr::mul(kappa, B_w_hat);
                    auto L_27 = Expr::constTerm(Matrix::sparse(2*n_p+n_omega,(*C_hat->getShape())[0]));
                    auto L_28 = Expr::mul(kappa,E_hat);
                    auto L_29 = Expr::constTerm(Matrix::sparse(2*n_p+n_omega,3*n_p));

                    auto L_34 = Expr::sub(R, S1->transpose());
                    auto L_35 = Expr::add(Expr::transpose(G_1_hat_theta), Expr::transpose(K_bb_theta));
                    auto L_36 = Expr::constTerm(Matrix::sparse(2*n_p+n_omega, (*B_w_hat->getShape())[1]));
                    auto L_37 = Expr::transpose(C_d_hat);
                    auto L_38 = Expr::constTerm(Matrix::sparse(2*n_p+n_omega,3*n_p));
                    auto L_39 = Expr::mul(epsilon, Expr::transpose(Expr::vstack(Expr::hstack(Expr::mul(G_B,C_theta_tilde_theta), Expr::mul(G_B, Expr::mul(D_K_theta, C_y)), Expr::constTerm(Matrix::sparse(n_p, n_omega))), 
                         Expr::constTerm(Matrix::sparse(2*n_p,2*n_p+n_omega)))));
                         
                    auto L_44 = Expr::mul(-1, Expr::add(S_h, R));
                    auto L_45 = Expr::constTerm(Matrix::sparse(2*n_p+n_omega, (*T_bar_theta->getShape())[1]));
                    auto L_46 = Expr::constTerm(Matrix::sparse(2*n_p+n_omega, (*B_w_hat->getShape())[1]));
                    auto L_47 = Expr::constTerm(Matrix::sparse(2*n_p+n_omega, (*C_hat->getShape())[0]));
                    auto L_48 = Expr::constTerm(Matrix::sparse(2*n_p+n_omega, 3*n_p));
                    auto L_49 = Expr::constTerm(Matrix::sparse(2*n_p+n_omega, 3*n_p));
                        
                    auto L_55 = Expr::mul(-4, T_bar_theta);
                    auto L_56 = Expr::mul(2, J_bar);
                    auto L_57 = Expr::hstack(Expr::mul(-1, Expr::mul(T_bar_theta, D_zu->transpose())), Expr::constTerm(Matrix::sparse(n_u, n_C_ew)));
                    auto L_58 = Expr::constTerm(Matrix::sparse((*T_bar_theta->getShape())[0], 3*n_p));
                    auto L_59 = Expr::mul(-epsilon, Expr::hstack(Expr::mul(T_bar_theta, G_B ->transpose()), Expr::constTerm(Matrix::sparse(n_u, 2*n_p))));
                  

                  
                    auto L_66 = Expr::mul(-1, Expr::constTerm(Matrix::eye((*B_w_hat->getShape())[1])));
                    auto L_67 = Expr::transpose(D_w);
                    auto L_68 = Expr::constTerm(Matrix::sparse((*D_w->getShape())[1], 3*n_p));
                    auto L_69 = Expr::mul(epsilon, Expr::transpose(Expr::vstack(Expr::hstack(Expr::constTerm(G_D), Expr::mul(G_B, Expr::constTerm(J))), Expr::constTerm(Matrix::sparse(2*n_p, n_d+n_nu)))) );
                    
                    // cout<<(*L_16->getShape())[0]<<'\n'<<(*L_26->getShape())[0]<<'\n'<<(*L_36->getShape())[0]<<'\n'<<(*L_46->getShape())[0]<<'\n'<<(*L_56->getShape())[0]<<'\n'<<(*L_66->getShape())[1]<<'\n'<<(*L_67->getShape())[1]<<'\n'<<(*L_68->getShape())[1]<<'\n'<<(*L_69->getShape())[1]<<'\n';

                   // cout<<(*L_15->getShape())[0]<<'\n'<<(*L_25->getShape())[0]<<'\n'<<(*L_35->getShape())[0]<<'\n'<<(*L_45->getShape())[0]<<'\n'<<(*L_55->getShape())[1]<<'\n'<<(*L_56->getShape())[1]<<'\n'<<(*L_57->getShape())[1]<<'\n'<<(*L_58->getShape())[1]<<'\n'<<(*L_59->getShape())[1]<<'\n';
              //      cout<<(*L_15->getShape())[0]<<'\t'<<(*L_25->getShape())[0]<<'\t'<<(*L_35->getShape())[0]<<'\t'<<(*L_45->getShape())[0]<<'\t'<<(*L_55->getShape())[1]<<'\t'<<(*L_56->getShape())[1]<<'\t'<<(*L_57->getShape())[1]<<'\t'<<(*L_58->getShape())[1]<<'\t'<<(*L_59->getShape())[1]<<'\n';
                    
                //    cout<<(*L_16->getShape())[0]<<'\t'<<(*L_26->getShape())[0]<<'\t'<<(*L_36->getShape())[0]<<'\t'<<(*L_46->getShape())[0]<<'\t'<<(*L_56->getShape())[0]<<'\t'<<(*L_66->getShape())[1]<<'\t'<<(*L_67->getShape())[1]<<'\t'<<(*L_68->getShape())[1]<<'\t'<<(*L_69->getShape())[1]<<'\n';

                    
                    auto L_77 = Expr::mul(-1, Expr::mul(gamma_sq, Matrix::eye((*D_w->getShape())[0])));
                    auto L_78 = Expr::constTerm(Matrix::sparse((*D_w->getShape())[0],3*n_p));
                    auto L_79 = Expr::constTerm(Matrix::sparse((*D_w->getShape())[0],3*n_p));
                    
                    auto L_88= Expr::mul(-epsilon, Expr::constTerm(Matrix::eye(3*n_p)));
                    auto L_89= Expr::constTerm(Matrix::sparse(3*n_p,3*n_p));
                        
                    auto L_99= Expr::mul(-epsilon, Expr::constTerm(Matrix::eye(3*n_p)));
                    
                    for (auto v_1 : {-v_1_val, v_1_val})
                        {
                        for (auto v_2 : {-v_2_val, v_2_val})
                            {
                            auto P_dot = Expr::add(Expr::mul(v_1, Expr::add(P_1, Expr::mul(rho_1, P_2))), Expr::mul(v_2, Expr::add(P_3, Expr::mul(rho_2, P_4))));
                            
                            auto L_d_tilde_dot = Expr::add(Expr::mul(v_1, Expr::add(L_d_tilde_1, Expr::mul(rho_1, L_d_tilde_2))), Expr::mul(v_2, Expr::add(L_d_tilde_3, Expr::mul(rho_2, L_d_tilde_4))));

                            auto A_hat = Expr::vstack(Expr::hstack(Expr::mul(A_p,Y), Expr::constTerm(A_p), Expr::mul(B_p, Expr::constTerm(V))), 
                                                      Expr::hstack(A_K_tilde, Expr::add(Expr::mul(X, A_p), Expr::mul(B_K_tilde, C_y)), Expr::mul(Expr::mul(X, B_p), V)),
                                                      Expr::hstack(L_K_tilde, Expr::add(Expr::mul(Expr::mul(L_d_tilde, C_y),A_p), Expr::mul(Expr::sub(L_d_tilde_dot, L_y_tilde), C_y)), Expr::add(Expr::mul(Z, W),
                                                                  Expr::mul(Expr::mul(Expr::mul(L_d_tilde, C_y), B_p), V))) );
                         
                            auto L_11 = Expr::add(A_hat, Expr::add(Expr::transpose(A_hat), Expr::add(Q_bar, Expr::add(S_bar, Expr::sub(P_dot,R)))));
                            auto L_12 = Expr::add(Expr::sub(P_bar, Expr::vstack(Expr::hstack(Y, Expr::constTerm(Matrix::eye(n_p)), Expr::constTerm(Matrix::sparse(n_p, n_omega))), 
                                                                                Expr::hstack(Expr::constTerm(Matrix::eye(n_p)), X, Expr::mul(Expr::mul(X, B_p), V)),
                                                                                Expr::hstack(Expr::constTerm(Matrix::sparse(n_omega, 2*n_p)), Z))), Expr::mul(kappa, Expr::transpose(A_hat)) );
                            for (auto v_3 : {-v_3_val, v_3_val})
                                {
                                auto L_33 = Expr::add(Expr::mul(-2, R), Expr::add(S1, Expr::sub(Expr::transpose(S1), Expr::mul(1-v_3*theta_prime, Q_bar_delayed) )));
                                
                                auto L = Expr::vstack(
                                                      Expr::hstack(L_11, L_12, Expr::hstack(L_13, L_14, Expr::hstack(L_15, L_16, Expr::hstack(L_17, L_18, L_19)))), 
                                                      Expr::hstack(Expr::transpose(L_12), L_22, Expr::hstack(L_23, L_24, Expr::hstack(L_25, L_26, Expr::hstack(L_27, L_28, L_29)))),
                                                      Expr::vstack(
                                                                   Expr::hstack(Expr::transpose(L_13), Expr::transpose(L_23), Expr::hstack(L_33, L_34, Expr::hstack(L_35, L_36, Expr::hstack(L_37, L_38, L_39)))),
                                                                   Expr::hstack(Expr::transpose(L_14), Expr::transpose(L_24), Expr::hstack(Expr::transpose(L_34), L_44, Expr::hstack(L_45, L_46, Expr::hstack(L_47, L_48, L_49)))),
                                                                   Expr::vstack(
                                                                                Expr::hstack(Expr::transpose(L_15), Expr::transpose(L_25), Expr::hstack(Expr::transpose(L_35), Expr::transpose(L_45), Expr::hstack(L_55, L_56,
                                                                                                                                                                                                     Expr::hstack(L_57, L_58, L_59)))),
                                                                                Expr::hstack(Expr::transpose(L_16), Expr::transpose(L_26), Expr::hstack(Expr::transpose(L_36), Expr::transpose(L_46), Expr::hstack(Expr::transpose(L_56),
                                                                                                                                                                                                 L_66, Expr::hstack(L_67, L_68, L_69)))),
                                                                                Expr::vstack(
                                                                                            Expr::hstack(Expr::transpose(L_17), Expr::transpose(L_27), Expr::hstack(Expr::transpose(L_37), Expr::transpose(L_47), 
                                                                                                                                      Expr::hstack(Expr::transpose(L_57), Expr::transpose(L_67), Expr::hstack(L_77, L_78, L_79)))),
                                                                                            Expr::hstack(Expr::transpose(L_18), Expr::transpose(L_28), Expr::hstack(Expr::transpose(L_38), Expr::transpose(L_48),
                                                                                                                     Expr::hstack(Expr::transpose(L_58), Expr::transpose(L_68), Expr::hstack(Expr::transpose(L_78), L_88,  L_89)))), 
                                                                                            Expr::hstack(Expr::transpose(L_19), Expr::transpose(L_29), Expr::hstack(Expr::transpose(L_39), Expr::transpose(L_49),
                                                                                                       Expr::hstack(Expr::transpose(L_59), Expr::transpose(L_69), Expr::hstack(Expr::transpose(L_79), Expr::transpose(L_89), L_99))))
                                                        ))));
                                                                  
                               M->constraint(Expr::sub(Expr::mul(-inequality, Expr::constTerm(Matrix::eye((*L->getShape())[0]))), L), Domain::inPSDCone((*L->getShape())[0] ));     
                              
                                }
                            }                        
                        }                     
                   
                    auto LMI_sat_d = Expr::vstack(Expr::hstack(beta, Expr::sub(K_bb_theta, G_1_hat_theta)), 
                                                  Expr::hstack(Expr::transpose(Expr::sub(K_bb_theta, G_1_hat_theta)), Expr::mul(pow(u_bar, 2), P_bar_delayed)));
                                            //  (Expr::sub(K_bb_theta, G_1_hat_theta))->slice( new_array_ptr<int,1>({0,0}), new_array_ptr<int,1>({0,2}) )   
                    M->constraint(LMI_sat_d, Domain::inPSDCone((*LMI_sat_d->getShape())[0] ));  
                    
                    auto LMI_sat = Expr::vstack(Expr::hstack(beta, Expr::sub(K_bb_theta, G_hat_theta)), 
                                                  Expr::hstack(Expr::transpose(Expr::sub(K_bb_theta, G_hat_theta)), Expr::mul(pow(u_bar, 2), P_bar)));

                    M->constraint(LMI_sat, Domain::inPSDCone((*LMI_sat->getShape())[0] )); 
                     
                    }
                }
                 
             }
        }
        M->objective(ObjectiveSense::Minimize, gamma_sq);
        //M->setSolverParam("intpntMaxIterations", 1);
        M->setSolverParam("intpntCoTolDfeas", 1.0e-4);
        M->setSolverParam("intpntCoTolRelGap", 1.0e-4);
        M->setSolverParam("intpntCoTolPfeas", 1.0e-4);


        M->solve();

    // This will print the value of gamma (note the use of level() method)
        cout << pow((*gamma_sq->level())[0],0.5);
      //  cout<<"Log1"<<'\n';
    return 0;
}
