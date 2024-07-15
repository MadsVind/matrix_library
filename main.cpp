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

    matrix1.addRow({1, 2, 3})
           .addRow({4, 5, 6})
           .addRow({7, 8, 9});

    matrix2.addRow({1, 2, 3, 4})
           .addRow({5, 6, 7, 8})
           .addRow({9, 10, 11, 12});

    Matrix<double>::Qr qr = matrix2.decompQR();

    qr.orthogonal.print();
    std::cout << "\n";
    qr.upper.print();

}
