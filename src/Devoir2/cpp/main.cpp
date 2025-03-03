#include <iostream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <fstream>
#include <string>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

void read_input_file(const std::string& filename, double& C_e, double& D_eff, double& R, int& N, double& k, double& T, int& K, double& C_i) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Failed to open input file." << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        std::string key;
        char equalSign;
        if (iss >> key >> equalSign) {
            if (key == "C_e") {
                iss >> C_e;
            } else if (key == "D_eff") {
                iss >> D_eff;
            }
            else if (key == "R") {
                iss >> R;
            }
            else if (key == "N") {
                iss >> N;
            }
            else if (key == "k") {
                iss >> k;
            }
            else if (key == "T") {
                iss >> T;
            }
            else if (key == "K") {
                iss >> K;
            }
            else if (key == "C_i") {
                iss >> C_i;
            }
        }
    }


    std::cout << "Read parameters from input file:\n";
    std::cout << "C_e = " << C_e << std::endl;
    std::cout << "D_eff = " << D_eff << std::endl;
    std::cout << "R = " << R << std::endl;
    std::cout << "N = " << N << std::endl;
    std::cout << "k = " << k << std::endl;
    std::cout << "T = " << T << std::endl;
    std::cout << "K = " << K << std::endl;
    std::cout << "C_i = " << C_i << std::endl;


    inputFile.close();
}

class Parameters {
public:
    double C_e, D_eff, R, k, T, C_i;
    int N, K;

    Parameters(double C_e, double D_eff, double R, int N, double k, double T, int K, double C_i)
        : C_e(C_e), D_eff(D_eff), R(R), k(k), T(T), C_i(C_i), N(N), K(K) {}
};
    
VectorXd solveSparseSystem(const MatrixXd& A, const VectorXd& B) {
    // Convert the dense matrix A to a sparse matrix
    SparseMatrix<double> sparseA = A.sparseView();

    // Solve the system using SparseLU
    SparseLU<SparseMatrix<double>> solver;
    solver.compute(sparseA);

    if (solver.info() != Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        return VectorXd();
    }

    VectorXd x = solver.solve(B);

    if (solver.info() != Success) {
        std::cerr << "Solving failed!" << std::endl;
        return VectorXd();
    }

    return x;
}

VectorXd solveSparseSystemBiCGSTAB(const MatrixXd& A, const VectorXd& B) {
    // Convert the dense matrix A to a sparse matrix
    SparseMatrix<double> sparseA = A.sparseView();

    // Create a BiCGSTAB solver
    BiCGSTAB<SparseMatrix<double>, IncompleteLUT<double>> solver;

    // Set tolerance and maximum iterations (optional)
    solver.setTolerance(1e-11);
    solver.setMaxIterations(1000);

    // Compute the solution
    solver.compute(sparseA);

    if (solver.info() != Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        return VectorXd();
    }

    VectorXd x = solver.solve(B);

    if (solver.info() != Success) {
        std::cerr << "Solving failed!" << std::endl;
        return VectorXd();
    }

    // std::cout << "BiCGSTAB solver iterations: " << solver.iterations() << std::endl;
    // std::cout << "Estimated error: " << solver.error() << std::endl;

    return x;
}


Eigen::ArrayXXd concentration_anal_MMS(Parameters& prm) {
    /*
    Function to compute the analytical solution for salt concentration in a pillar.

    Parameters:
    ----------
    prm : Object of class Parametres
        - D_eff : Effective diffusion coefficient of salt in concrete [1/s]
        - R : Radius of the pillar [m]
        - N : Number of spatial nodes [-]
        - k : Reaction constant [1/s]
        - C_e : Salt concentration at the surface of the pillar [mol/m^3]
        - K : Number of time steps [-]
        - T : Total simulation time [s]

    Returns:
    -------
    C_anal : Eigen::MatrixXd
        Matrix of salt concentration values according to position and time [mol/L]
    */

    int N = prm.N;
    int K = prm.K;

    // Spatial and temporal discretization
    Eigen::ArrayXd r = Eigen::ArrayXd::LinSpaced(N, 0, prm.R);
    Eigen::ArrayXd t = Eigen::ArrayXd::LinSpaced(K, 0, prm.T);

    // Creating a meshgrid equivalent in Eigen
    auto R = r.replicate(1, K).transpose();  // Replicating r as columns
    auto T = t.replicate(1, N);  // Replicating t as rows

    // Computing the concentration matrix
    auto C_anal = prm.C_e*(exp(1.0e-10*T) - 1)*cos(0.5*M_PI*R.pow(2)/pow(prm.R, 2)) + prm.C_e;

    return C_anal;
}

Eigen::VectorXd source_MMS(Parameters& prm, Eigen::VectorXd& r, double& t) {
    /*
    Fonction qui calcule le terme source selon la solution MMS

    Parameters
    ----------
    prm : Objet class parametres()
        - D_eff : Coefficient de diffusion effectif du sel dans le béton [1/s]
        - R : Rayon du pilier [m]
        - N : Nombre de noeuds [-]
        - k : Constante de réaction [1/s]
        - Ce : Concentration en sel à la surface du poteau [mol/m^3]
        - K : Nombre de pas de temps [-]
        - T : Temps total de simulation [s]
    r : float
        Position radiale à laquelle le terme source est calculé [m]
    t : float
        Temps à laquelle le terme source est calculé [s]

    Returns
    -------
    f : float
        Valeur du terme source à la position r et au temps t [mol/m^3/s]

    */
    double C_e = prm.C_e;
    double R = prm.R;
    double k = prm.k;
    double Deff = prm.D_eff;
    double pi = M_PI;
    Eigen::ArrayXd f = C_e*(-M_PI*Deff*(1 - exp(1.0e-10*t))*(2*pow(R, 2)*sin((1.0/2.0)*M_PI*r.array().pow(2)/pow(R, 2)) 
    + M_PI*r.array().pow(2)*cos((1.0/2.0)*M_PI*r.array().pow(2)/pow(R, 2))) 
    + pow(R, 4)*(k*(-(1 - exp(1.0e-10*t))*cos((1.0/2.0)*M_PI*r.array().pow(2)/pow(R, 2)) + 1) 
    + 1.0e-10*exp(1.0e-10*t)*cos((1.0/2.0)*M_PI*r.array().pow(2)/pow(R, 2))))/pow(R, 4);

    return f.matrix();
}

double source_zero(Parameters& prm, double& r, double& t) {
    /*
    Fonction qui calcule le terme source selon la solution MMS

    Parameters
    ----------
    prm : Objet class parametres()
        - D_eff : Coefficient de diffusion effectif du sel dans le béton [1/s]
        - R : Rayon du pilier [m]
        - N : Nombre de noeuds [-]
        - k : Constante de réaction [1/s]
        - Ce : Concentration en sel à la surface du poteau [mol/m^3]
        - K : Nombre de pas de temps [-]
        - T : Temps total de simulation [s]
    r : float
        Position radiale à laquelle le terme source est calculé [m]
    t : float
        Temps à laquelle le terme source est calculé [s]

    Returns
    -------
    f : float
        Valeur du terme source à la position r et au temps t [mol/m^3/s]

    */

    return 0;
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::MatrixXd> concentration_mdf_o2(Parameters& prm, 
    Eigen::VectorXd (*source)(Parameters&, Eigen::VectorXd&, double&)) {
    /*
    Fonction qui simule, selon la méthode des différence finies avec un 
    schéma d'approximation d'ordre 2, la concentration de sel dans un pilier dans le temps
    

    Parameters
    ----------
    prm : Objet class parametres()
        - D_eff : Coefficient de diffusion effectif du sel dans le béton [1/s]
        - R : Rayon du pilier [m]
        - N : Nombre de noeuds [-]
        - k : Constante de réaction [1/s]
        - C_e : Concentration en sel à la surface du poteau [mol/m^3]
        - K : Nombre de pas de temps [-]
        - T : Temps total de simulation [s]
        - C_i : Concentration initiale [mol/m^3]
    
    source : Fonction
        Fonction qui calcule le terme source selon la solution

    Returns
    -------
    r : Vecteur (array) de dimension N composé de la position radiale à 
        laquelle les concentrations sont calculées, où N le nombre de noeuds.
    dr : float représentant de la pas de discrétisation spatial selon le rayon 
        du pilier
    C : Vecteur (array) composée de la concentration en sel selon la 
        position [mol/L]
    */
    int N = prm.N;
    int K = prm.K;

    Eigen::VectorXd r = Eigen::VectorXd::LinSpaced(N, 0, prm.R);
    Eigen::VectorXd t = Eigen::VectorXd::LinSpaced(K, 0, prm.T);

    double dr = prm.R / (N - 1);
    double dt = prm.T / (K - 1);
    std::cout << "dr :" << dr << std::endl;
    std::cout << "dt :" << dt << std::endl;

    Eigen::MatrixXd C_all = Eigen::MatrixXd::Zero(K, N);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N, N);
    Eigen::VectorXd B = Eigen::VectorXd::Zero(N);

    // Boundary conditions
    A(0, 0) = -3;
    A(0, 1) = 4;
    A(0, 2) = -1;
    A(N-1, N-1) = 1;
    B(N-1) = prm.C_e;

    // Initial condition
    Eigen::VectorXd C_t = Eigen::VectorXd::Constant(N, prm.C_i);
    std::cout << "C_i :" << prm.C_i << std::endl;
    C_all.row(0) = C_t;
    
    for (int i=1; i<N-1; i++) {
        A(i, i-1) = prm.D_eff * dt * (1 / (dr * dr) - 1 / (r(i) * 2 * dr));
        A(i, i) = -2 * prm.D_eff * dt / (dr * dr) - 1 - prm.k * dt;
        A(i, i+1) = prm.D_eff * dt * (1 / (dr * dr) + 1 / (r(i) * 2 * dr));
    }

    for (int k=1; k<K; k++) {
        // std::cout << "Time step: " << k << std::endl;
        B(Eigen::seq(1, Eigen::last-1)) = -C_t(Eigen::seq(1, Eigen::last-1)) - dt * source(prm, r, t(k))(Eigen::seq(1, Eigen::last-1));
        // C_t = A.partialPivLu().solve(B);
        // C_t = solveSparseSystem(A, B);
        C_t = solveSparseSystemBiCGSTAB(A, B);
        C_all.row(k) = C_t;
    }

    return {r, t, C_all};
}

// Function to compute the discrete L1 norm of the error
double erreur_L1(const MatrixXd& u, const MatrixXd& u_ref) {
    double L1 = (u - u_ref).array().abs().sum() / u.size();
    return L1;
}

// Function to compute the discrete L2 norm of the error
double erreur_L2(const MatrixXd& u, const MatrixXd& u_ref) {
    double L2 = std::sqrt((u - u_ref).array().square().sum() / u.size());
    return L2;
}

// Function to compute the discrete Linf norm of the error
double erreur_inf(const MatrixXd& u, const MatrixXd& u_ref) {
    double L_inf = (u - u_ref).array().abs().maxCoeff();
    return L_inf;
}

void save_to_csv(const std::string& filename, const Eigen::MatrixXd& r, const Eigen::MatrixXd& t, const Eigen::MatrixXd& C_all) {
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    // Set precision to maximum for double values
    file << std::setprecision(15);

    // Write the header
    file << "Time/Position";
    for (int j = 0; j < r.size(); ++j) {
        file << "," << r(j);
    }
    file << "\n";

    // Write the data (each row corresponds to a time step)
    for (int i = 0; i < t.size(); ++i) {
        file << t(i);  // Write the time value
        for (int j = 0; j < r.size(); ++j) {
            file << "," << C_all(i, j);
        }
        file << "\n";
    }

    file.close();
    std::cout << "Data saved to " << filename << std::endl;
}

int main(int argc, char* argv[]) {

    auto time_start = std::chrono::high_resolution_clock::now();

    // Ensure the input file is passed as a command line argument
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <parameters_filename>" << std::endl;
        return 1;
    }

    std::string input_filename = argv[1]; // The input filename is passed as the first argument

    // Read parameters from the input file
    double C_e, D_eff, R, k, T, C_i;
    int N, K;

    read_input_file(input_filename, C_e, D_eff, R, N, k, T, K, C_i);

    Parameters prm(C_e, D_eff, R, N, k, T, K, C_i);

    // Résolution du problème
    auto C_anal_MMS = concentration_anal_MMS(prm);
    auto [r, t, C_mdf_o2] = concentration_mdf_o2(prm, source_MMS);

    save_to_csv("output.csv", r, t, C_mdf_o2);

    auto time_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = time_end - time_start;
    std::cout << "Elapsed time: " << elapsed_time.count() << " s" << std::endl;

    // std::cout << "C_anal_MMS \n" << C_anal_MMS << std::endl;
    // std::cout << "C_mdf_o2 \n" << C_mdf_o2 << std::endl;

    double L1 = erreur_L1(C_mdf_o2, C_anal_MMS.matrix());
    double L2 = erreur_L2(C_mdf_o2, C_anal_MMS.matrix());
    double L_inf = erreur_inf(C_mdf_o2, C_anal_MMS.matrix());

    std::cout << std::fixed << std::setprecision(15);
    std::cout << "L1 error = " << L1 << std::endl;
    std::cout << "L2 error = " << L2 << std::endl;
    std::cout << "L_inf error = " << L_inf << std::endl;
}