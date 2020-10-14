#include <iostream>
#include <map>
#include <vector>
#include <random>
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

public:
    Neuron(int idInLayer, FFLayer * l){
        this->idInLayer = idInLayer;
        this->layer = l;
    }
    void addConnection(Connection * c){
        this->connections.push_back(c);
    }
    void setLayer(FFLayer * l){
        this->layer = l;
    }
    vector<Connection*> getConnections(){
        return connections;
    }
    void setWeight(Connection *c){
        this->connections.push_back(c);
    }
    void countActivation(){
        float sum = 0;
        for(int i = 0; i < connections.size(); i ++){
            Connection *c = connections.at(i);
            sum += (c->getInputNeuron()->getOutput() * c->getWeight());
         //   cout<<c->getWeight();
        }
        float tmp = (1.f + pow(e, -sum));
        this->output = (1.f/tmp);
        cout<<"";
    }
    void setOutput(float o){
        this->output = o;
    }
    float getOutput(){
        return output;
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
    void initRandom(int neuronsCount);
    void setNet(Network * n){
        this->net = n;
    }
    void run();
};
class Network{
    vector<ConvolutionalLayer*> convolutionalLayers;
    vector<FFLayer*> ffLayers;
    vector<float> input;

public:
    Network(int convLayersCount, int ffLayersCount){
        for(int i = 0; i < convLayersCount; i ++){
            convolutionalLayers.push_back(new ConvolutionalLayer(this,i));
        }
        for(int i = convLayersCount; i < ffLayersCount; i ++){
            ffLayers.push_back(new FFLayer(this, i));
        }
    }
    Network(Network * p1, Network * p2){
        if(p1->ffLayers.size() == p2->ffLayers.size())
            for(int i = 0; i < p1->ffLayers.size(); i ++){
                bool r = std::rand() % 2;
                FFLayer * l = r ? p1->ffLayers.at(i) : p2->ffLayers.at(i);
                for(int j = 0; j < l->getNeurons().size(); j ++){
                    l->getNeurons().at(j)->setLayer(l);
                }
                this->ffLayers.push_back(l);
            }
    }
    void run(vector<float> input){
        this->input = input;
        //tmp until convolution;
        //todo
        if(input.size() == ffLayers.at(0)->getNeurons().size()){
            for(int i = 0; i < ffLayers.size(); i ++){
                ffLayers.at(i)->run();
            }
        }else{
            cout<<"ERROR in = " << input.size() <<" exp = " << ffLayers.at(0)->getNeurons().size();
        }
    }
    vector<float> getInput(){
        return input;
    }
    void setFFLayer(int x, FFLayer * l){
        this->ffLayers[x] = l;
    }
    static float random(float min, float max){
        float r = (max - min) * ((float)rand() / RAND_MAX) + min;
      //  cout<<"\nrand = " << r ;
        return r;
    }
    vector<FFLayer *> getFFLayers(){
    	return ffLayers;
    }
    void mutate(float mutationRate){
        for(int i = 0; i < ffLayers.size(); i ++){
            for(int j = 0; j < ffLayers.at(i)->getNeurons().size();j++){
                vector<Connection *> conns = ffLayers.at(i)->getNeurons().at(j)->getConnections();
                for(int k = 0; k < conns.size(); k ++){
                    if(random(0, 1000) <= mutationRate*10){
                        conns.at(k)->setWeight(conns.at(k)->getWeight() + Network::random(-1,1));
                    }
                }
            }
        }
    }
    vector<ConvolutionalLayer *> getConvolutionLayers(){
    	return convolutionalLayers;
    }
    void initRandom(){
        for(int i = 0; i < ffLayers.size(); i ++){
            FFLayer * l = ffLayers.at(i);
            l->initRandom(16);
        }
    }

    vector<float> getResult(){
        vector<float> result;
        vector<Neuron *> neurons = ffLayers.at(ffLayers.size()-1)->getNeurons();
        for(int i =0; i < neurons.size(); i ++){
            result.push_back(neurons.at(i)->getOutput());
        }
        return result;
    }
};
void ConvolutionalLayer::initRandom() {
    if(idInNet != 0) {
        net->getConvolutionLayers().at(idInNet - 1);
    }
}
void FFLayer::initRandom(int neuronsCount) {
    for(int i = 0; i < neuronsCount; i ++){
        Neuron * n = new Neuron(i, this);
        neurons.push_back(n);
    }
    if(idInNet != 0) {
       vector<Neuron *> prevV = net->getFFLayers().at(idInNet-1)->getNeurons();
       for(int i = 0; i < neuronsCount; i ++){
           Neuron * n = neurons.at(i);
           for(int j = 0; j < prevV.size(); j++){
               Neuron * pn = prevV.at(j);
               Connection * c = new Connection(pn, n, net->random(0,1));
                n->addConnection(c);
           }
       }
    }
}
void FFLayer::run(){
    if(idInNet == 0){
        vector<float> output = net->getInput();
        for(int i = 0; i < neurons.size(); i ++){
            neurons.at(i)->setOutput(output.at(i));
        }
    }else {
        for (int i = 0; i < neurons.size(); i++) {
            Neuron *n = neurons.at(i);
            n->countActivation();
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
                    multidimensionalSum += ((float ) input.at(z).at(x+i))*weightSet.at(i);
                }
            }
            result.push_back(multidimensionalSum);
        }
    }
};
class Population{
    vector<Network *> species;
    int speciesPerGeneration;
    float mutationRate;

    vector<float> curInput;
    vector<float> currTarget;
    float currPopulationLoss;
    int ffLayersCount;
    map<Network *, float> adaptationLevel;
    float bestSpecieScore;

public:
    Population(int speciesPerGeneration, float mutationRate, int ffLayersCount){
        this->speciesPerGeneration = speciesPerGeneration;
        this->mutationRate = mutationRate;
        this->ffLayersCount = ffLayersCount;
    }
    float getMutationRate(){
        return mutationRate;
    }
    float getBestSpecieScore(){
        return bestSpecieScore;
    }
    void init(){
        for(int i = 0; i < speciesPerGeneration; i ++){
            Network * specie = new Network(0, ffLayersCount);
            specie->initRandom();
            this->species.push_back(specie);
        }
    }

    void runTraining(vector<float> input, vector<float> target){
        float totalLoss = 0;

        float bestScore = 0;
        for(int i = 0; i < species.size(); i ++){
            Network * specie = species.at(i);
            specie->run(input);
            vector<float> result = specie->getResult();
            if(result.size() == target.size()){
                float loss = 0;
                for(int j = 0; j < result.size(); j ++){
                    loss += pow((result.at(j) - target.at(j)), 2);
                }
                adaptationLevel[specie] = (1.f/(loss))*100.f;
                if(loss < bestScore){
                    bestScore = loss;
                }
                bestSpecieScore = bestScore;
                totalLoss += loss;
            }else{
                cout<<"error";
            }
        }
        this->currPopulationLoss = totalLoss;
    }
    struct RouletteElement{
        float start;
        float end;
        Network * specie;
    };
    void cross(){
        float sum = 0;
        vector<RouletteElement*> roulette;
        for (std::map<Network*, float>::iterator it=adaptationLevel.begin(); it!=adaptationLevel.end(); ++it){
            RouletteElement * r = new RouletteElement();
            r->specie = it->first;
            r->start = sum;
            sum += it->second;
            r->end = sum;
            roulette.push_back(r);
        }
        vector<Network *> children;
        for(int i = 0; i < speciesPerGeneration; i ++){
            float randomSpecie1 = Network::random(0, sum);
            float randomSpecie2;
            Network *n1 = nullptr;
            Network *n2 = nullptr;
            while(randomSpecie2 = Network::random(0, sum) && randomSpecie2 == randomSpecie1);
            for(int j = 0; j < roulette.size(); j++){
                float s = roulette.at(j)->start;
                float e = roulette.at(j)->end;
                if(s <= randomSpecie1 && e >= randomSpecie1){
                    n1 = roulette.at(j)->specie;
                }
                if(s <= randomSpecie2 && e >= randomSpecie2){
                    n2 = roulette.at(j)->specie;
                }
                if(n1 != nullptr && n2 != nullptr) break;
            }
            Network *child = new Network(n1,n2);
            child->mutate(mutationRate);
            children.push_back(child);
        }
        free(species);
        this->species = children;
        this->adaptationLevel.empty();
    }
    float getCurrPopulationLoss(){
        return currPopulationLoss;
    }
};

int main(){
    srand((unsigned int)time(NULL));
    Population * population = new Population(20,2,4);
    population->init();

    int maxIterations = 100000000;
    int age = 1;
    vector<float> input = {0.44, 0.55, 0.1, 0.6, 0.44, 0.55, 0.1, 0.6, 0.44, 0.55, 0.1, 0.6, 0.44, 0.55, 0.1, 0.6};
    vector<float> output = {0.0f, 0.0f, 0.0f, 0.0f, 1.f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,0.0f, 0.0f,};
    while(age <= maxIterations){
        population->runTraining(input, output);
        age++;
        if(age % 25 == 0){
            population->cross();
        }
        if(age % 100 == 0){
            cout << "Specie Loss = " << population->getBestSpecieScore()<<"\n";
            cout << "Loss = " << population->getCurrPopulationLoss()<<"\n";
        }
    }
    return 0;
}