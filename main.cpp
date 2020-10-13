#include <iostream>
#include <map>
#include <vector>
#include <math.h>

#include <stdlib.h>

using namespace std;
class Network;
class FFLayer;
class Neuron;
class ConvolutionalLayer;
class Connection;

class Kernel;
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
    void setWeight(float w){
        this->weight = w;
    }
    float getWeight(){
        return weight;
    }
    Neuron * getInputNeuron(){
        return inputNeuron;
    }
};
class Neuron{
    const float e = 2.71828;
    float bias;
    int idInLayer;
    FFLayer * layer;
    vector<Connection*> connections;
    float output;
    void countActivation(){
        float sum = 0;
        for(int i = 0; i < connections.size(); i ++){
            Connection *c = connections.at(i);
            sum += c->getInputNeuron()->getOutput()*c->getWeight();
        }
        this->output = (1.f/1.f-pow(e, -sum));
    }
    float getOutput(){
        return output;
    }
public:
    Neuron(int idInLayer, FFLayer * l){
        this->idInLayer = idInLayer;
        this->layer = l;
    }
    void addConnection(Connection * c){
        this->connections.push_back(c);
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
    void initRandom();
};
class FFLayer{
    Network* net;
    int idInNet;
    vector<Neuron *> neurons;
public:
    FFLayer(Network *net, int id){
        this->net = net;
        this->idInNet = id;
    }
    vector<Neuron *> getNeurons(){
        return neurons;
    }
    void initRandom();
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

    static double randZeroToOne(){
        return rand() / (RAND_MAX + 1.);
    }
    vector<FFLayer *> getFFLayers(){
    	return ffLayers;
    }
    vector<ConvolutionalLayer *> getConvolutionLayers(){
    	return convolutionalLayers;
    }
    void initRandom();
};
void ConvolutionalLayer::initRandom() {
    if(idInNet != 0) {
        net->getConvolutionLayers().at(idInNet - 1);
    }
}
void FFLayer::initRandom() {
    if(idInNet != 0) {
       vector<Neuron *> prevV = net->getFFLayers().at(idInNet-1)->getNeurons();
       for(int i = 0; i < neurons.size(); i ++){
           Neuron * n = neurons.at(i);
           for(int j = 0; j < prevV.size(); j++){
               Neuron * pn = prevV.at(j);
               Connection * c = new Connection(n, pn, net->randZeroToOne());
                n->addConnection(c);
           }
       }
    }
}
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

int main(){
    Network net();
    return 0;
}