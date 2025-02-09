#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

class Matrix {
public:
    Matrix(int i, int j) : i(i), j(j), Matrixbuffer(std::make_unique<double[]>(i * j)) {
    }

    void setI(int i) {
        Matrix::i = i;
    }

    void setJ(int j) {
        Matrix::j = j;
    }

    int getJ() const {
        return j;
    }

    int getI() const {
        return i;
    }

    void insert(int element, int pos_x, int pos_y) {
        bool valid = ValidPos(pos_x, pos_y);
        if (!valid) {
            std::cout << "Unvalid Position." << std::endl;
            return;
        }
        Matrixbuffer[pos_y * i + pos_x] = element;
    }

    int getElement(int posx, int posy) const {
        bool valid = ValidPos(posx, posy);
        if (valid)
            return Matrixbuffer[posy * i + posx];
        else
            return 0;
    }


    Matrix operator+(const Matrix &othermatrix) {
        if (i == othermatrix.getI() && j == othermatrix.getJ()) {
            Matrix C(i, j);
            int totalElements = i * j;
            // Accesso diretto ai buffer interni
            double *aData = Matrixbuffer.get();
            double *bData = othermatrix.Matrixbuffer.get();
            double *cData = C.Matrixbuffer.get();
            for (int idx = 0; idx < totalElements; ++idx) {
                cData[idx] = aData[idx] + bData[idx];
            }
            C.PrintMatrix();
            return C;
        }
        std::cout << "nothing is happening..." << std::endl;
        return Matrix(0, 0);
    }

    Matrix operator-(const Matrix &othermatrix) {
        if (i == othermatrix.getI() && j == othermatrix.getJ()) {
            Matrix C(i, j);
            int totalElements = i * j;
            // Accesso diretto ai buffer interni
            double *aData = Matrixbuffer.get();
            double *bData = othermatrix.Matrixbuffer.get();
            double *cData = C.Matrixbuffer.get();
            for (int idx = 0; idx < totalElements; ++idx) {
                cData[idx] = aData[idx] - bData[idx];
            }
            C.PrintMatrix();
            return C;
        }
        std::cout << "nothing is happening..." << std::endl;
        return Matrix(0, 0);
    }


    int Rank() const {
        int m = i; // numero di righe
        int n = j; // numero di colonne

        // Copia della matrice in un vettore di double per una maggiore precisione
        std::vector<double> A(m * n);
        for (int r = 0; r < m; r++) {
            for (int c = 0; c < n; c++) {
                A[r * n + c] = static_cast<double>( getElement(r, c) );
            }
        }

        int rank = 0;
        int currentRow = 0;
        const double eps = 1e-9; // tolleranza per considerare un valore come zero

        // Loop sulle colonne
        for (int col = 0; col < n && currentRow < m; col++) {
            // Trova la riga pivot per la colonna col, cioè la prima riga (a partire da currentRow)
            // in cui l'elemento è "non nullo"
            int pivotRow = currentRow;
            while (pivotRow < m && std::fabs(A[pivotRow * n + col]) < eps) {
                pivotRow++;
            }
            if (pivotRow == m) {
                // Nessun pivot in questa colonna: passa alla successiva
                continue;
            }
            // Se il pivot non è già in currentRow, scambia le righe
            if (pivotRow != currentRow) {
                for (int k = col; k < n; k++) {
                    std::swap(A[currentRow * n + k], A[pivotRow * n + k]);
                }
            }
            // Abbiamo trovato un pivot: incrementa il conteggio
            rank++;

            // Elimina gli elementi sotto il pivot nella colonna col
            for (int r = currentRow + 1; r < m; r++) {
                // Se l'elemento è già zero, non serve eliminare
                if (std::fabs(A[r * n + col]) < eps)
                    continue;
                double factor = A[r * n + col] / A[currentRow * n + col];
                for (int c = col; c < n; c++) {
                    A[r * n + c] -= factor * A[currentRow * n + c];
                }
            }
            // Passa alla riga successiva per il prossimo pivot
            currentRow++;
        }
        std::cout << "Rank:"<< rank << std::endl;
        return rank;
    }


    void PrintMatrix(){
        for(int l = 0; l<getI() ; l++){
            for(int k = 0;k<getJ();k++){
                std::cout<< getElement(l,k) << " ";
            }
            std::cout << std::endl;
        }
    }

private:
    int i;
    int j;
    std::unique_ptr<double[]> Matrixbuffer;
private:
    bool ValidPos(int pos_x, int pos_y) const {
        if (pos_x < 0 || pos_y < 0 || pos_x >= i || pos_y >= j)
            return false;
        return true;
    }
};

int main() {
    Matrix matrix(3,3);
    Matrix matrix2(3,3);

    //first Matrix
    matrix.insert(1,0,0);
    matrix.insert(2,0,1);
    matrix.insert(3,0,2);
    matrix.insert(4,1,0);
    matrix.insert(5,1,1);
    matrix.insert(6,1,2);
    matrix.insert(7,2,0);
    matrix.insert(8,2,1);
    matrix.insert(9,2,2);

    //the second
    matrix2.insert(9, 0, 0);
    matrix2.insert(8, 0, 1);
    matrix2.insert(7, 0, 2);
    matrix2.insert(6, 1, 0);
    matrix2.insert(5, 1, 1);
    matrix2.insert(4, 1, 2);
    matrix2.insert(3, 2, 0);
    matrix2.insert(2, 2, 1);
    matrix2.insert(1, 2, 2);


    //first Matrix
    std::cout << "First Matrix:" << " " << std::endl;
    matrix.PrintMatrix(); std::cout << " " << std::endl;
    std::cout << "Second Matrix:" << " " << std::endl;
    matrix2.PrintMatrix();

    std::cout << "Sum of both:" << std::endl;
    // Uso della funzione statica per ottenere la somma
    Matrix result = matrix + matrix2;
    matrix.Rank();




}
