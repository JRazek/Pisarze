#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
using namespace std;

vector<string> split(string str, char divider){
    vector<string> result;

    string currWord = "";
    for(int i = 0; i < str.size(); i ++){
        currWord+=str[i];
        if(str[i] == divider || str.size()-1 == i){
            result.push_back(currWord);
            currWord = "";
        }
    }
    return result;
}
string mickiewicz[] = {"tadeusz","wojski","gerwazy","telimena","klucznik","tadeusza","podkomorzy","wtenczas","litwie","asesor","zosia","horeszkow","jacek","soplica","maciej","soplicow","moskali","protazy","dabrowski","strzelcy","bernardyn","dobrzynscy","biezy","chrzciciel","soplicowie","stolnik","wojskiego","maciek","charty","jegry","krzykneli","zosie","dobrzynski","mopanku","rzecze","kropiciel","kwestarz","rapier","bernardyna","buchman","konewka","kropic","stolnika","telimenie","telimeny","hrabie","jankiel","karabele","klucznika","macka","majorze","moskal","moskala","polowanie","powiecie","scyzoryk","strzelby","chartow","dobrzynskich","dobyl","hejze","jegrow","litwa","mysliwi","soplice","stola","litewskich","prusak","ramiony","rebajlo","rykowa","soplicy","szabli","szczwacze","zascianku","ascka","asesora","chrzciciela","dobrzynie","drazki","dziedzic","konopie","kusego","litwy","lowach","moskalow","mysliwcow","obadwa","ochmistrzyni","pleban","pluta","podkomorzego","ptastwo","skoluba","sokol","soplicowa","wloscian","wszerz","zareczyny","brzytewka",};
string sienkiewicz[] = {"winicjusz","petroniusz","ligia","cezara","cezar","winicjusza","ligii","albowiem","ligie","chilo","ursus","petroniusza","piotr","chrzescijanie","apostol","eunice","ancjum","nerona","pomponia","atrium","pomponii","niewolnicy","chilona","ostrianum","tygellinus","aulusow","bogowie","winicjuszowi","poppea","aulus","ligio","apostola","niemal","chilon","glauka","tygellina","cyrku","naokol","uczynic","ursusa","pretorianow","tygellin","uczynil","tarsu","chrzescijaninem","glaukus","arenie","aulusa","cezarem","cezarowi","krotona","chryzotemis","marku","winicjuszu","augustianie","poppei","triclinium","arene","augustianow","niewolnik","kryspus","lektyki","achai","cezarze","czynil","greka","odrzekla","petroniuszu","kroto","linus","linusa","nazariusz","winicjuszem","amfiteatrze","ligow","miriam","palatynie","seneka","wiezieniu","petroniuszowi","westynus","cubiculum","zatybrze","amfiteatr","areny","odnalezc","poczucie","willi","dziewice","grecji","pretorianie","chrzescijanami","grecyna","neapolis","plaucjusz","centurion","igrzyska","kryspa","legii","lektyke",};
string prus[] = {"wokulski","izabela","wokulskiego","rubli","rzecki","baron","ignacy","izabeli","stach","ochocki","baronowa","wasowska","starski","stawska","szlangbaum","wokulskiemu","stawskiej","prezesowa","szuman","naturalnie","wokulskim","adwokat","hrabina","maruszewicz","fotelu","geist","lecki","baronowej","izabele","mraczewski","florentyna","klejn","krzeszowska","odparla","stacha","kupiec","lisiecki","procent","kamienice","rzeckiego","prezesowej","stanislawie","kobiete","staruszka","szlangbauma","mincel","misiewiczowa","klacz","starskiego","subiekt","lecka","paryzu","spolki","stosunki","wirski","hrabiny","krzeszowskiej","suzin","leckich","przedpokoju","geista","radca","stawka","meliton","felicja","salonie","wasowskiej","diabla","mezczyzna","stachu","wtracila","szumana","gabinetu","wariat","wyobraz","adwokata","leckiego","maruszewicza","ochockiego","racji","raczek","spolke","krzeszowski","rezultacie","spolka","helunia","ignacego","leckiej","mikolaj","szescdziesiat","diabli","hopfera","kupca","lalka","pociag","kupcem","miesiecy","mraczewskiego","pietnascie","starskim",};


string toLower(string data){
    std::for_each(data.begin(), data.end(), [](char & c) {
        c = ::tolower(c);
    });
    return data;
}
bool isInList(char c){
    string restrictedChars = " ,.\"!';?";
    for(int i = 0; i < restrictedChars.size(); i++){
        if(restrictedChars[i] == c)
            return true;
    }
    return false;
}
bool isInMap(map<string, int> m, string key){
    if (m.find(key) == m.end())
        return false;
    return true;
}
int main() {
    string line;
    getline(cin, line);
    vector<string> args = split(line, ' ');
    int n = stoi(args[0]);

    for(int i = 0; i < n; i ++){
        map<string, int> words;//int for time of occurrences
        getline(cin, line);

        string currWord = "";
        for(int i = 0; i < line.size(); i ++){
            char c = line[i];
            if(isInList(c)){
                if(currWord.size() >= 5) {
                    if (isInMap(words, currWord)) {
                        words[currWord]++;
                    } else {
                        words[currWord] = 1;
                    }
                }
                currWord = "";
            } else{
                currWord += c;
            }
        }
        for (const auto& [key, value] : words) {
            std::cout << key << " x" << value << std::endl;
        }
        cout<<"\n\n\n\n";
    }
    return 0;
}
