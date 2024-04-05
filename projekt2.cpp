#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>

using namespace std;

const int N = 977;
const double A1 = 11.0;
const double Da1 = 3.0;
const double a2 = -1.0;
const double a3 = -1.0;
const int maxIter = 1000;
const double tol = 1e-9;

void fillMatrix(double** A, double* b, double a1, int N)
{
    // wypełnianie wektora b
    for (int i = 0; i < N; i++)
    {
        b[i] = sin((i) * (2 + 1));

        // wypełnianie macierzy A
        for (int j = 0; j < N; j++)
        {
            if (i == j) // na diagonali macierzy A jest wartość a1
            {
                A[i][j] = a1;
            }
            else if (i == j + 1 || i == j - 1) // nad i pod diagonali są wartości a2
            {
                A[i][j] = a2;
            }
            else if (i == j + 2 || i == j - 2) // dwie kolejne wartości nad/pod diagonali są a3
            {
                A[i][j] = a3;
            }
            else // pozostałe wartości są zerami
            {
                A[i][j] = 0.0;
            }
        }
    }
}

double jacobi(double** A, double* b, double* x, int N, bool save) {
    double* prev_x = new double[N]; // wektor poprzedniego wyniku
    double* res = new double[N]; // wektor residuum
    double* tmp_x; // tymczasowy wektor 
    double* normRes = new double[maxIter]; // wektor do przechowywania norm residuum
    int iterations = 0; // liczba wykonanych iteracji

    for (int i = 0; i < N; i++) {
        prev_x[i] = 0; // początkowe wartości wektora prev_x
        x[i] = 1; // początkowe wartości wektora x
    }

    auto start = chrono::high_resolution_clock::now(); // pomiar czasu rozpoczyna się tutaj

    // wykonywanie iteracji
    while (iterations < maxIter) {
        // obliczanie nowych wartości wektora x
        for (int i = 0; i < N; i++) {
            double s = 0;

            // obliczenie sumy A_ij * x_j^(k) dla j != i
            for (int j = 0; j < N; j++) {
                if (i != j) {
                    s += A[i][j] * prev_x[j];
                }
            }
            // obliczenie nowej wartości x_i^(k+1) na podstawie wzoru
            x[i] = (b[i] - s) / A[i][i];
        }

        // sprawdzenie warunku stopu
        for (int i = 0; i < N; i++) {
            // obliczenie wektora residuum r^(k)
            res[i] = b[i];
            for (int j = 0; j < N; j++) {
                res[i] -= A[i][j] * x[j];
            }
        }

        // obliczenie normy residuum ||r^(k)||_2
        double norm_res = 0;
        for (int i = 0; i < N; i++) {
            norm_res += pow(res[i], 2);
        }
        norm_res = sqrt(norm_res);
        normRes[iterations] = norm_res;
        if (norm_res < tol) {
            break;
        }

        // zapisanie wyniku z tego przebiegu iteracji jako poprzedniego wyniku
        tmp_x = prev_x;
        prev_x = x;
        x = tmp_x;

        iterations++;
    }

    auto end = chrono::high_resolution_clock::now(); // pomiar czasu kończy się tutaj

    // wyświetlenie wyników
    cout << "Jacobi method:" << endl;
    cout << "Number of iterations: " << iterations << endl;
    //chrono::duration<double> diff = end - start;
    auto diff = chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "Elapsed time: " << diff.count() << " milliseconds" << endl;

    if (save) {
        ofstream outfile("residuumJacobi.txt");

        // Zapisanie wektora z normami residuum
        for (int i = 0; i <= iterations; ++i) {
            outfile << normRes[i] << ",";
        }
        outfile << endl;

        outfile.close();
    }


    return diff.count();
}

double gaussSeidel(double** A, double* b, double* x, int N, bool save) {
    double* prev_x = new double[N]; // wektor poprzedniego wyniku
    double* res = new double[N]; // wektor residuum
    double* tmp_x; // tymczasowy wektor
    double* normRes = new double[maxIter]; // wektor do przechowywania norm residuum
    int iterations = 0; // liczba wykonanych iteracji

    for (int i = 0; i < N; i++) {
        prev_x[i] = 0; // początkowe wartości wektora prev_x
        x[i] = 1; // początkowe wartości wektora x
    }

    auto start = chrono::high_resolution_clock::now();

    // wykonywanie iteracji
    while (iterations < maxIter) {
        // obliczanie nowych wartości wektora x
        for (int i = 0; i < N; i++) {
            double s = 0;

            // obliczanie sumy elementów powyżej i-tej diagonali
            for (int j = 0; j < i; j++) {
                s += A[i][j] * x[j];
            }

            // obliczanie sumy elementów poniżej i-tej diagonali z poprzedniej iteracji
            for (int j = i + 1; j < N; j++) {
                s += A[i][j] * prev_x[j];
            }

            // obliczenie i-tej współrzędnej wektora x na podstawie wzoru
            x[i] = (b[i] - s) / A[i][i];
        }

        // sprawdzenie warunku stopu
        for (int i = 0; i < N; i++) {
            res[i] = b[i];
            for (int j = 0; j < N; j++) {
                res[i] -= A[i][j] * x[j];
            }
        }

        // obliczenie normy wektora residuum
        double norm_res = 0;
        for (int i = 0; i < N; i++) {
            norm_res += pow(res[i], 2);
        }

        norm_res = sqrt(norm_res);
        normRes[iterations] = norm_res;
        if (norm_res < tol) {
            break;
        }

        // zapisanie wyniku z tego przebiegu iteracji jako poprzedniego wyniku
        tmp_x = prev_x;
        prev_x = x;
        x = tmp_x;

        iterations++;
    }

    auto end = chrono::high_resolution_clock::now(); // pomiar czasu kończy się tutaj

    // wyświetlenie wyników
    cout << "Gauss-Seidel method:" << endl;
    cout << "Number of iterations: " << iterations << endl;
    //chrono::duration<double> diff = end - start;
    auto diff = chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "Elapsed time: " << diff.count() << " milliseconds" << endl;
    
    if (save) {
        ofstream outfile("residuumGauss.txt");

        // Zapisanie wektora norm residuum
        for (int i = 0; i < iterations; ++i) {
            outfile << normRes[i] << ",";
        }
        outfile << endl;

        outfile.close();
    }


    return diff.count();
}

double LU(double** A, double* b, double* x, int N) {
    double** L = new double* [N];
    double** U = new double* [N];

    auto start = chrono::high_resolution_clock::now();

    for (int i = 0; i < N; i++) {
        L[i] = new double[N];
        U[i] = new double[N];
        for (int j = 0; j < N; j++) {
            L[i][j] = 0;
            U[i][j] = 0;
        }
    }

    // kopiowanie macierzy A do macierzy U
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // przypisanie wartości macierzy A do macierzy U
            U[i][j] = A[i][j];
        }
    }

    // tworzenie macierzy L i U
    for (int j = 0; j < N; j++) {
        // wartość na przekątnej macierzy L jest równa 1
        L[j][j] = 1;
        for (int i = j + 1; i < N; i++) {
            // obliczenie wartości elementu macierzy L na podstawie elementów macierzy U
            L[i][j] = U[i][j] / U[j][j];
            for (int k = j; k < N; k++) {
                // aktualizacja wartości elementu macierzy U zgodnie ze wzorem
                U[i][k] -= L[i][j] * U[j][k];
            }
        }
    }

    // rozwiązanie układu równań
    double* y = new double[N];
    for (int i = 0; i < N; i++) {
        double s = 0;
        // obliczenie wartości wektora y na podstawie macierzy L i wektora b
        for (int j = 0; j < i; j++) {
            s += L[i][j] * y[j];
        }
        y[i] = b[i] - s;
    }
    for (int i = N - 1; i >= 0; i--) {
        double s = 0;
        // obliczenie wartości wektora x na podstawie macierzy U i wektora y
        for (int j = i + 1; j < N; j++) {
            s += U[i][j] * x[j];
        }
        x[i] = (y[i] - s) / U[i][i];
    }

    // obliczenie residuum
    double* r = new double[N];
    for (int i = 0; i < N; i++) {
        r[i] = b[i];
        for (int j = 0; j < N; j++) {
            r[i] -= A[i][j] * x[j];
        }
    }
    double res_norm = 0;
    for (int i = 0; i < N; i++) {
        res_norm += pow(r[i], 2);
    }
    res_norm = sqrt(res_norm);

    auto end = chrono::high_resolution_clock::now(); // pomiar czasu kończy się tutaj

    // wyświetlenie wyników
    cout << "LU method:" << endl;
    cout << "Residuum norm: " << res_norm << endl;
    //chrono::duration<double> diff = end - start;
    auto diff = chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "Elapsed time: " << diff.count() << " milliseconds" << endl;

    return diff.count();
}

void displayMatrix(double** A, double* b) {
    std::cout << "Macierz A:\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << A[i][j] << " ";
        }
        std::cout << "\n";
        std::cout << " | " << b[i] << std::endl;
    }
}


int main()
{
    // Tworzenie macierzy A[N][N] i wektorów b[N], x[N]
    double** A = new double* [N];
    for (int i = 0; i < N; i++)
    {
        A[i] = new double[N];
    }
    double* b = new double[N];
    double* x = new double[N];

    // Wypełnienie macierzy A[N][N] i wektora b[N] wartościami
    fillMatrix(A, b, A1, N);

    // Zadanie B - Rozwiązanie układu równań metodami Jacobiego i Gaussa-Seidla
    std::cout << "Zadanie B" << std::endl;
    jacobi(A, b, x, N, true);
    gaussSeidel(A, b, x, N, true);

    // Zadanie C - zamienic 5 parametr w metodach aby zapisc wyniki z zadania b
    std::cout << "Zadanie C" << std::endl;
    fillMatrix(A, b, Da1, N);
    jacobi(A, b, x, N, false);
    gaussSeidel(A, b, x, N, false);

    //Zadnie D
    std::cout << "Zadanie D" << std::endl;
    fillMatrix(A, b, Da1, N);
    LU(A, b, x, N);
    
    //Zadanie E
    std::cout << "\nZadanie E" << std::endl;
    int nValues[10] = { 100, 200, 300, 500, 800, 1200, 1500, 2000, 2500, 3000 };
    double timeJacobi[10];
    double timeGauss[10];
    double timeLU[10];

    for (int i = 0; i < 10; i++) {
        std::cout << "N = "<< nValues[i] << std::endl;

        double** A = new double* [nValues[i]];
        for (int j = 0; j < nValues[i]; j++)
        {
            A[j] = new double[nValues[i]];
        }
        double* b = new double[nValues[i]];
        double* x = new double[nValues[i]];

        fillMatrix(A, b, A1, nValues[i]);
        timeJacobi[i] = jacobi(A, b, x, nValues[i], false);
        timeGauss[i] = gaussSeidel(A, b, x, nValues[i], false);
        timeLU[i] = LU(A, b, x, nValues[i]);
    }

    // Otwarcie pliku do zapisu
    ofstream outfile("tables.txt");

    // Zapisanie timeJacobi
    for (int i = 0; i < 10; ++i) {
        outfile << timeJacobi[i] << ",";
    }
    outfile << endl;

    // Zapisanie tabeli timeGauss
    for (int i = 0; i < 10; ++i) {
        outfile << timeGauss[i] << ",";
    }
    outfile << endl;

    // Zapisanie tabeli timeLU
    for (int i = 0; i < 10; ++i) {
        outfile << timeLU[i] << ",";
    }
    outfile << endl;

    // Zamknięcie pliku
    outfile.close();

    for (int i = 0; i < N; i++)
    {
        delete[] A[i];
    }
    delete[] A;
    delete[] b;
    delete[] x;

    return 0;
}


