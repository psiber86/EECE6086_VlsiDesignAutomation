#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

int main(int argc, char* argv[])
{
    int line;
    int linenum = 0;
    int cellcount = 0;
    int netcount = 0;
    std::vector<int> a, b;
    int cell;

    if (argc < 2) {
        std::cout << "Usage: klalgo <filename>" << std::endl;
        exit(0);
    }

    //get netlist
    std::ifstream infile(argv[1]);
    if (infile.is_open()) {
        while (infile.good()) {
            if (linenum == 0) {
                infile >> cellcount;
            } else if (linenum == 1) {
                infile >> netcount;
            } else {
                infile >> cell;
                a.push_back(cell);
                infile >> cell;
                b.push_back(cell);
            }
            linenum++;
        }
        std::cout << cellcount << std::endl;
        std::cout << netcount << std::endl;
        for (int i = 0; i < a.size(); i++) {
            std::cout << a[i] << " " << b[i] << std::endl;
        }
        infile.close();
    } else {
        std::cout << "Unable to open file" << std::endl;
    }

    //TODO: call KL algorithm

    return 0;
}
