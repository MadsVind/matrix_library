#include "main.hpp"

typedef std::vector<std::vector<std::vector<float>>> Weights;

std::vector<float> createRandomVector(int size) {
    std::vector<float> randomVector;
    for (int i = 0; i < size; i++) {
        randomVector.push_back(rand() % 100 / 100.0f);
    }
    return randomVector;
}

Weights createWeights(std::vector<int> layers) {
    Weights weights;
    for (int i = 0; i < layers.size() - 1; i++) {
        std::vector<std::vector<float>> layerWeights;
        for (int j = 0; j < layers[i]; j++) {
            layerWeights.push_back(createRandomVector(layers[i + 1]));
        }
        weights.push_back(layerWeights);
    }
    return weights;
}

std::vector<float> createBias(std::vector<int> layers) {
    std::vector<float> bias;
    for (int i = 1; i < layers.size(); i++) {
        bias.push_back(rand() % 100 / 100.0f);
    }
    return bias;
}

void printWeights(Weights weights) {
    for (int i = 0; i < weights.size(); i++) {
        for (int j = 0; j < weights[i].size(); j++) {
            for (int k = 0; k < weights[i][j].size(); k++) {
                std::cout << weights[i][j][k] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}


std::vector<float> sigmoid(std::vector<float> input) {
    std::vector<float> output;
    for (int i = 0; i < input.size(); i++) {
        output.push_back(1 / (1 + exp(-input[i])));
    }
    return output;
}

//std::vector<float> feedForward(std::vector<float> input, Weights weights, std::vector<float> bias) {
//    if (input.size() != weights[0].size()) {
//        std::cout << "Input size does not match input layer size" << std::endl;
//        return {};
//    }
//    std::vector<float> output;
//    for (int i = 0; i < weights.size(); i++) {
//        
//    }
//}


int main() {
    Matrix<double> matrix1;
    Matrix<double> matrix2;
    Matrix<double> matrix3;
    Matrix<double> matrix4;
    Matrix<double> matrix5;
    Matrix<double> matrix6;

    matrix1.addRow({2, 0, 2, 0.6})
           .addRow({3, 3, 4, -2})
           .addRow({5, 5, 4, 2})
           .addRow({-1, -2, 3.4, -1});

    matrix2.addRow({1, 0, 0})
           .addRow({0, 1, 0})
           .addRow({0, 0, 1});

    matrix3.addRow({0, 2, 2})
           .addRow({1, 2, 3})
           .addRow({4, 4, 4});

    matrix4.addRow({1, 2})
           .addRow({5, 4});

    matrix5.addRow({1, 1, 0})
           .addRow({1, 0, 1})
           .addRow({0, 1, 1})
           .addRow({1, 0, 1});

    matrix6.addRow({1, 1})
           .addRow({1, -1});

    Matrix<double>::Qr qu = matrix5.qr();

    //std::cout << "\n\n Determinant: " << matrix5.determinant() << "\n";

    //qr.orthogonal.print();
    //std::cout << "\n";
    //qr.upper.print();


    Matrix<double> q = qu.orthogonal;
    Matrix<double> r = qu.upper;
    
//    q.addRow({0.577, 0.516, -0.632})
//     .addRow({0.577, -0.258, 0.316})
//     .addRow({0, 0.775, 0.632})
//     .addRow({0.577, -0.258, 0.316});
//
//    r.addRow({1.732, 0.577, 1.155})
//     .addRow({0,     1.291, 0.258})
//     .addRow({0,     0,     1.265});

    std::cout << "reconstructed A\n\n";

    (q * r).print(0.001);

    std::cout << "\nr:\n";
    r.print();
    
    std::cout << "\nq.T:\n";
    q.transpose().print();
    
    std::cout << "\nq:\n";
    q.print();

    
    std::cout << "\n";

    (q.transpose() * q).print(0.001);

    std::cout << "3\n\n";

    return 0;
}
