#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <bitset>
#include <cmath>

#define FASTA "/home/zuzu/CLionProjects/akwb3/sample.fasta"
#define QUAL "/home/zuzu/CLionProjects/akwb3/sample.qual"

using namespace std;
struct sequence {
    string name;
    string nucleotydes;
    vector<long int> qual;
};

struct node {
    string name;
    string nucleotydes;
    vector<long int> qual;
    unsigned int number_of_sequence = 0;
    unsigned int position = 0;
    vector<node *> connected_nodes;

    bool operator==(node &other) {
        return name == other.name && number_of_sequence == other.number_of_sequence && position == other.position;
    }
};

bool find_node(vector<node *> &connected, node &node1) {
    for (auto &it : connected) {
        if (*it == node1) return true;
    }
    return false;
}

vector<sequence> read_from_file() {
    vector<sequence> all_sequences;
    sequence one_sequence;
    ifstream file(FASTA);
    while (file.good()) {
        getline(file, one_sequence.name);
        getline(file, one_sequence.nucleotydes);
        all_sequences.push_back(one_sequence);
    }
    file.close();

    file.open(QUAL);
    string name_of_seq;
    string line;
    unsigned int temp_qual;
    while (getline(file, name_of_seq))
        for (auto &sequence_r : all_sequences)
            if (name_of_seq == sequence_r.name) {
                getline(file, line);
                istringstream l(line);
                while (l >> temp_qual) sequence_r.qual.push_back(temp_qual);
                if (sequence_r.nucleotydes.size() != sequence_r.qual.size()) {
                    cout << "Wielkość sekwencji nukleotydów nie zgadza się z liczbą oceny wiarygodności nukleotydów"
                         << endl;
                    cout << "Sekwencja: " << sequence_r.name << endl << "Nukleotydów: " << sequence_r.nucleotydes.size()
                         << endl
                         << "Wiarygodności: "
                         << sequence_r.qual.size() << endl;
                    exit(1);
                }
                break;
            }
    file.close();

    cout << endl << endl;
    return all_sequences;
}

vector<vector<node>> make_graph(unsigned int &n, int q, vector<sequence> &all_sequences) {
    vector<node> nodes_in_seq;
    vector<vector<node>> graph;
    unsigned int number = 0;
    unsigned int pos = 0;
    for (auto &sequence_r : all_sequences) {
        pos = 0;
        nodes_in_seq.clear();
        auto it = sequence_r.nucleotydes.begin() + n;
        auto sec_it = sequence_r.nucleotydes.begin();
        for (unsigned int i = 0; i < (sequence_r.nucleotydes.length() - n + 1); i++) {
            node temp_node;
            temp_node.number_of_sequence = number;
            temp_node.name = sequence_r.name;
            temp_node.nucleotydes = string(sec_it, it);
            sec_it++;
            it++;
            temp_node.position = pos;
            temp_node.qual = vector<long>(sequence_r.qual.begin() + pos, sequence_r.qual.begin() + pos + n);
            pos++;
            nodes_in_seq.push_back(temp_node);
        }
        number++;
        graph.push_back(nodes_in_seq);
    }

    for (auto &seq1 : graph)
        for (auto &seq : graph)
            for (auto &node_of_seq : seq1)
                for (auto &node_of_other : seq)
                    if (node_of_seq.number_of_sequence != node_of_other.number_of_sequence &&
                        (node_of_seq.nucleotydes == node_of_other.nucleotydes)) {
                        if (!find_node(node_of_seq.connected_nodes, node_of_other))
                            node_of_seq.connected_nodes.push_back(&node_of_other);
                        if (!find_node(node_of_other.connected_nodes, node_of_seq))
                            node_of_other.connected_nodes.push_back(&node_of_seq);
                    }


    for (auto &sequences : graph)
        for (auto &note1 : sequences)
            for (unsigned int quality_it = 0; quality_it < note1.qual.size(); quality_it++)
                if (note1.qual[quality_it] < q) {
                    note1.nucleotydes.erase(note1.nucleotydes.begin() + quality_it);
                    note1.qual.erase(note1.qual.begin() + quality_it);
                    quality_it--;
                }

    for (auto &seq1 : graph)
        for (auto &seq : graph)
            for (auto &node_of_seq : seq1)
                for (auto &node_of_other : seq)
                    if (node_of_seq.number_of_sequence != node_of_other.number_of_sequence &&
                        (node_of_seq.nucleotydes == node_of_other.nucleotydes ||
                         node_of_seq.nucleotydes.find(node_of_other.nucleotydes) != string::npos)) {
                        if (!find_node(node_of_seq.connected_nodes, node_of_other))
                            node_of_seq.connected_nodes.push_back(&node_of_other);
                        if (!find_node(node_of_other.connected_nodes, node_of_seq))
                            node_of_other.connected_nodes.push_back(&node_of_seq);
                    }

    return graph;
}

bool check_candidate(vector<node> &temp_clique, node &candidate) {
    for (auto &it : temp_clique)
        if (!find_node(it.connected_nodes, candidate))
            return false;
    return true;
}

bool check_positions_in_last_clique(vector<node> &last_clique, node &candidate) {
    for (auto &it : last_clique)
        if (it.number_of_sequence == candidate.number_of_sequence && it.position + 1 != candidate.position)
            return false;

    return true;
}

vector<vector<node>> find_series_of_cliques(vector<vector<node> *> &graph) {
    vector<vector<node>> temp_series;
    vector<vector<node>> best_series;
    vector<node> temp_clique;
    for (auto &current : *(graph[0])) {
        temp_clique.push_back(current);
        for (unsigned int i = 1; i < graph.size(); i++)
            for (auto &candidate : *(graph[i]))
                if (temp_series.empty()) {
                    if (check_candidate(temp_clique, candidate)) {
                        temp_clique.push_back(candidate);
                        break;
                    }
                } else {
                    if (check_candidate(temp_clique, candidate) &&
                        check_positions_in_last_clique(temp_series.back(), candidate)) {
                        temp_clique.push_back(candidate);
                        break;
                    }
                }

        if (temp_clique.size() == graph.size()) {
            temp_series.push_back(temp_clique);
            temp_clique.clear();
            if (temp_series.size() > best_series.size())
                best_series = temp_series;
        } else {
            temp_clique.clear();
            temp_series.clear();
        }
    }
    return best_series;
}


vector<vector<unsigned int>> comb(unsigned int &K) {
    string temp;
    vector<vector<unsigned int>> results;
    for (unsigned int i = 0; i < pow(2, 7); i++) {
        bitset<7> t(i);
        if (t.count() == K) {
            vector<unsigned int> v1;
            for (unsigned int j = 0; j < 7; j++)
                if (t[j]) v1.push_back(j);
            results.push_back(v1);
        }
    }
    return results;
}

void print_best_clique(vector<sequence> seq, vector<vector<node>> best, unsigned int &n) {
    string temp;
    cout << "Najlepsza seria składa się z klik o " << best[0].size() << " wierzchołkach" << endl;
    cout << "Długość serii: " << best.size() << endl;
    for (auto &it : best[0]) {
        cout << "W sekwencji nr: " << it.number_of_sequence << " od pozycji nr: " << it.position << " motyw: ";
        temp = seq[it.number_of_sequence].nucleotydes.substr(it.position, (best.size() + n - 1));
        cout << temp << endl;
    }

}


int main() {
    unsigned int n;
    unsigned int q;
    while (1){
        cout << "Podaj rozmiar okna (od 4 do 7): " << endl;
        cin >> n;
        if(n < 4 || n > 7) continue;
        else break;
    }
    cout << endl << "Podaj próg wiarygodności " << endl;
    cin >> q;
    cout << endl;
    vector<sequence> seq = read_from_file();
//    cout << "read_from_file completed" << endl;
    vector<vector<node>> graph = make_graph(n, q, seq); // n, q, seq
//    cout << "make_graph completed" << endl;
    vector<vector<node> *> temp_graph;
    vector<vector<node>> temp_series;
    vector<vector<node>> best_of_comb;
    vector<vector<node>> best;

    for (unsigned int k = 7; k >= 4; k--) {
//        cout << "k=" << k << " started" << endl;
        int completed = 0;
        vector<vector<unsigned int>> combinations = comb(k);
        for (auto &i : combinations) {
            for (auto &j: i) {
                temp_graph.push_back(&(graph[j]));
            }
            temp_series = find_series_of_cliques(temp_graph);
            temp_graph.clear();
            completed++;
//            cout << "find_series_of_cliques completed " << completed << endl;
            if (temp_series.size() > best_of_comb.size()) best_of_comb = temp_series;
        }
        if (best_of_comb.size() > best.size()) best = best_of_comb;
    }
    if (best.empty())
        cout << endl << "W podanych sekwencjach nie odnaleziono motywu (serii klik)" << endl;
    else
        print_best_clique(seq, best, n);

    return 0;
}