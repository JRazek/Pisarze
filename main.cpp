#include <iostream>
#include <map>
#include <vector>

using namespace std;
class Network;
class FFLayer;
class Neuron;
class ConvolutionalLayer;
class Connection;

class Kernel;

class Neuron{
    float bias;
    int idInLayer;
    FFLayer * layer;
    vector<Connection*> connections;
public:
    Neuron(int idInLayer, FFLayer * l){
        this->idInLayer = idInLayer;
        this->layer = l;
    }
    void setWeight(Connection *c){
        this->connections.push_back(c);
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
    void initRandom(){

    }
};
class FFLayer{
    Network* net;
    int idInNet;
    vector<Neuron *> neurons;
    bool inputLayer = false;
public:
    FFLayer(Network *net, int id){
        this->net = net;
        this->idInNet = id;
    }
    void setAsInputLayer(){
        this->inputLayer = true;
    }
    void initRandom(){
        for(int i = 0; i < neurons.size(); i++){
            if(!inputLayer){
                //Connection * conn = new Connection();
                //neurons.at(i)->setWeight();
                int prevLayer = idInNet -1;

            }
        }
    }
};
class Network{
    vector<ConvolutionalLayer*> convolutionalLayers;
    vector<FFLayer*> ffLayers;
public:
    Network(int convLayersCount, int ffLayersCount){
        for(int i = 0; i < convLayersCount; i ++){
            convolutionalLayers.push_back(new ConvolutionalLayer(this,i));
        }
        for(int i = convLayersCount; i < ffLayersCount; i ++){
            ffLayers.push_back(new FFLayer(this, i));
        }
    }
    vector<FFLayer *> getFFLayers(){
    	return ffLayers;
    }
    vector<ConvolutionalLayer *> getConvolutionLayers(){
    	return convolutionalLayers;
    }
    void initRandom(){
        for(int i = 0; i < convolutionalLayers.size(); i++){
            convolutionalLayers.at(i)->initRandom();
        }
        for(int i = 0; i < ffLayers.size(); i++){
            ffLayers.at(i)->initRandom();
        }
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
class Connection{
    float weight;
    Neuron * inputNeuron;
    Neuron * outputNeuron;
public:
    Connection(Neuron * inputNeuron, Neuron * outputNeuron, float w){
        this->inputNeuron = inputNeuron;
        this->outputNeuron = outputNeuron;
        this->weight = w;
    }
};
int main(){
    Network net();
    return 0;
}