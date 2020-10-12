#include <iostream>
#include <map>
#include <vector>

using namespace std;
class Network;
class FFLayer;
class Neuron;
class ConvolutionalLayer;
class Kernel;

class Network{
    vector<ConvolutionalLayer> convolutionalLayers;
    vector<FFLayer*> ffLayers;
public:
    Network(int convLayerCount, int ffLayersCount){
        for(int i = 0; i < convLayerCount; i ++){
            
        }
    }
};
class ConvolutionalLayer{
    vector<Kernel*> kernels;
    Network* net;
    int idInNet;
public:
    ConvolutionalLayer(Network *net, int id){
        this->net = net;
        this->idInNet = id;
    }
};
class Kernel{
    vector<vector<float>> weights;
    int idInLayer;
    int kernelSize;
    ConvolutionalLayer *layer;
public:
    Kernel(ConvolutionalLayer* layer, int id, int kernelSize){
        this->layer = layer;
        this->idInLayer = id;
        this->kernelSize = kernelSize;
    }
    void setWeightsZ(vector<float> weightSet, int z){
        this->weights.at(z) = weightSet;
    }
    float getValue(int x, int z){
        return weights.at(z).at(x);
    }
    vector<string> getConvolutionResult(vector<string> input){
        vector<float> result;
        for(int x = 0; x < input.at(0).size(); x++){
            float multidimensionalSum = 0;
            for(int z = 0; z < input.size(); z ++){
                vector<float> weightSet = weights.at(z);
                for (int i = 0; i < weightSet.size(); i ++){
                    multidimensionalSum += ((int) input.at(z).at(x+i))*weightSet.at(i);
                }
            }
            result.push_back(multidimensionalSum);
        }
    }
};
class FFLayer{

};
class Neuron{
    float weight;
    float bias;
    int idInLayer;
    Network * net;
    FFLayer * layer;
};

int main(){
    Network net();
    return 0;
}