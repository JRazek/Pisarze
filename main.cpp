#include <iostream>
#include <map>
#include <vector>
#include <random>
#include <math.h>
#include <fstream>

#include <stdlib.h>

using namespace std;
class Network;
class FFLayer;
class Neuron;
class Connection;

class Connection{
    float weight;
    Neuron * inputNeuron;
    Neuron * outputNeuron;

    bool chainSet = false;
    float chain;
public:
    Connection(Neuron * inputNeuron, Neuron * outputNeuron, float w){
        this->inputNeuron = inputNeuron;
        this->outputNeuron = outputNeuron;
        this->weight = w;
    }
    void setWeight(float w){
        this->weight = w;
    }
    void setChain(float chain){
        this->chain = chain;
        chainSet = true;
    }
    void resetChain(){
        this->chainSet = false;
    }
    bool isChainSet(){
        return chainSet;
    }
    float getChain(){
        return chain;
    }
    float getWeight(){
        return weight;
    }
    Neuron * getInputNeuron(){
        return inputNeuron;
    }
    Neuron * getOutputNeuron(){
        return outputNeuron;
    }
};

static float sigmoid(float x){
    static const float e = 2.71828;
    return (1.f/(1.f + pow(e, -x)));
}
class Neuron{
    float bias;
    int idInLayer;
    FFLayer * layer;
    vector<Connection*> inputConnections;
    vector<Connection*> outputConnections;
    float output;
    float netVal;

public:
    Neuron(int idInLayer, FFLayer * l){
        this->idInLayer = idInLayer;
        this->layer = l;
    }
    void addInputConnection(Connection * c){
        this->inputConnections.push_back(c);
    }
    void addOutputConnection(Connection * c){
        this->outputConnections.push_back(c);
    }
    void setLayer(FFLayer * l){
        this->layer = l;
    }
    FFLayer * getLayer(){
        return layer;
    }
    vector<Connection*> getInputConnections(){
        return inputConnections;
    }
    vector<Connection*> getOutputConnections(){
        return outputConnections;
    }
    void setWeight(Connection *c){
        this->inputConnections.push_back(c);
    }
    void countActivation(){
        float sum = 0;
        for(int i = 0; i < inputConnections.size(); i ++){
            Connection *c = inputConnections.at(i);
            sum += (c->getInputNeuron()->getOutput() * c->getWeight());
        }
        this->netVal = sum;
        this->output = sigmoid(sum);
        cout<<"";
    }
    float getNetVal(){
        return netVal;
    }
    void setOutput(float o){
        this->output = o;
    }
    int getIndex(){
        return idInLayer;
    }
    float getOutput(){
        return output;
    }
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
    int getIndexInNet(){
        return idInNet;
    }
    void run();
};
class Network{
    vector<FFLayer*> ffLayers;
    vector<float> input;
    vector<Connection *> connList;

public:
    Network(int ffLayersCount){
        for(int i = 0; i < ffLayersCount; i ++){
            ffLayers.push_back(new FFLayer(this, i));
        }
    }
    vector<Connection * > getConnList(){
        return connList;
    }
    void addToConnList(Connection *c){
        this->connList.push_back(c);
    }
    vector<float> getOutputFromLastLayer(){
        vector<float> out;
        FFLayer *last = ffLayers.at(ffLayers.size()-1);
        for(int i = 0; i < last->getNeurons().size(); i ++){
            out.push_back(last->getNeurons().at(i)->getOutput());
        }
        return out;
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
        return r;
    }
    vector<FFLayer *> getFFLayers(){
    	return ffLayers;
    }
    void initRandom(){
        for(int i = 0; i < ffLayers.size(); i ++){
            FFLayer * l = ffLayers.at(i);
            int neurons = 16;
            if(i == ffLayers.size()-1)
                neurons = 3;
            else if(i == 0)
                neurons = 40;
            l->initRandom(neurons);
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
               n->addInputConnection(c);
               pn->addOutputConnection(c);
               net->addToConnList(c);
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
class BackProp{
    const float learningRate = 1;
    Network * net;
    float getChain(Connection *c, vector<float> expected) {
        float chain = 1;
        float act = sigmoid(c->getOutputNeuron()->getNetVal());
        chain *= act * (1.f - act);
        float tmp = 0;
        if(c->isChainSet())
            return c->getChain();

        if (c->getOutputNeuron()->getLayer()->getIndexInNet() != net->getFFLayers().size() - 1) {
            vector<Connection *> connList = c->getOutputNeuron()->getOutputConnections();
            for (int i = 0; i < connList.size(); i++) {
                Connection *nc = connList.at(i);
                tmp += nc->getWeight() * getChain(nc, expected);
            }
            chain *= tmp;
        } else{
            chain *= 2*(c->getOutputNeuron()->getOutput() - expected.at(c->getOutputNeuron()->getIndex()));
        }
        c->setChain(chain);
        //todo reset chains in neuron!!
        return chain;
    }

public:
    BackProp(Network * net){
        this->net = net;
    }
    float getLoss(vector<float> expected){
        vector<Neuron *> neurons = net->getFFLayers().at(net->getFFLayers().size()-1)->getNeurons();
        float loss = 0;
        for(int i = 0; i < neurons.size(); i++){
            loss += pow(neurons.at(i)->getOutput() - expected.at(i), 2);
        }
        return loss;
    }
    void learn(vector<float> expected){
        if(net->getFFLayers().at(net->getFFLayers().size() - 1)->getNeurons().size() != expected.size()){
            cout<<"ERROR!!";
            return;
        }
        for(int i = 0; i < net->getConnList().size(); i ++){
            net->getConnList().at(i)->resetChain();
        }
        for(int i = 0; i < net->getConnList().size(); i ++){
            Connection * c = net->getConnList().at(i);
            float delta = learningRate * c->getInputNeuron()->getOutput() * getChain(c, expected);
            c->setWeight(c->getWeight() - delta);
            cout<<"";
        }
    }
};
class TextManager{
    const string path = "/home/user/CLionProjects/Pisarze/texts/";
    fstream mic = fstream (path + "Mickiewicz.txt");
    fstream sien = fstream (path + "Sienkiewicz.txt");
    fstream prus = fstream (path + "Prus.txt");
    const int batchSize = 40;
public:
    struct RandomText{
        vector<float> text;
        int author;
    };
    RandomText getRandomText(){
        mic.seekg(0);
        sien.seekg(0);
        prus.seekg(0);
        string text;
        int author = Network::random(0, 3);
        RandomText rt;
        if(!(mic.good() && sien.good() && prus.good())){
            cout<<"ERROR";
        }
        if(author == 0){
            //micus
            getline(mic, text);
        } else if(author == 1){
            //prusak
            getline(sien, text);
        } else if(author == 2){
            //twoj stary pijany
            getline(prus, text);
        }
        int randStart = Network::random(0, text.size() - batchSize);
        vector<float> cutText;
        if(randStart < 0){
            cout<<"";
        }
        for(int i = 0; i < batchSize; i ++){
            cutText.push_back(text[randStart + i]/255.f);
        }
        cout<<"";
        rt.text = cutText;
        rt.author = author;
        return rt;
    }
};
class Accuracy{
    int timesMeasuredFromLastReset = 0;
    int timesOKFromLastReset = 0;
public:
    void measure(TextManager::RandomText rt, vector<float> prediction){
        float max = 0;
        int author = -1;
        for(int i = 0; i < prediction.size(); i ++){
            if(prediction[i] > max) {
                max = prediction[i];
                author = i;
            }
        }
        if(author == rt.author){
            timesOKFromLastReset++;
        }
        timesMeasuredFromLastReset++;
    }
    float getAccuracyPercentage(){
        return ((float)timesOKFromLastReset)/((float)timesMeasuredFromLastReset)*100.f;
    }
    void reset(){
        timesMeasuredFromLastReset = 0;
        timesOKFromLastReset = 0;
    }
};

int main(){

    srand((unsigned int)time(NULL));
    Network * network = new Network(10);
    network->initRandom();

    BackProp *backProp = new BackProp(network);
    TextManager textManager;
    Accuracy acc;
   // vector<float> input = {0.44, 0.55, 0.1, 0.6, 0.44, 0.55, 0.1, 0.6, 0.44, 0.55, 0.1, 0.6, 0.44, 0.55, 0.1, 0.6};
    int resetRate = 1000;
    for(int i = 0; i < 100000; i ++){
        TextManager::RandomText randomText = textManager.getRandomText();
        vector<float> expected = {0.0f, 0.0f, 0.0f};
        expected[randomText.author] = 0.99f;
        network->run(randomText.text);
        backProp->learn(expected);
        acc.measure(randomText, network->getOutputFromLastLayer());
        if((i + 1) % resetRate == 0){
            cout << "Loss = " << backProp->getLoss(expected)<<"\n";
            cout << "Accuracy = " << acc.getAccuracyPercentage()<<"%\n";
            acc.reset();
        }
    }
    return 0;
}