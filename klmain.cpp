#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cstring>

//global variables
bool debug = false;
int cellcount = 0;
int netcount = 0;
int cutset;

enum {
    SETA = 0,  
    SETB = 1
}; 

enum {
    UNLOCKED = 0,
    LOCKED
};

void debug_print_matrices(std::map<int, std::map<int, int> > mat, std::string label)
{
    std::cout << "Printing matrix " <<  label << std::endl;
    for (std::map<int, std::map<int, int> >::iterator i1 = mat.begin(); 
            i1 != mat.end(); i1++)
    {
        for (std::map<int, int>::iterator i2 = i1->second.begin(); 
             i2 != i1->second.end(); i2++) 
        {
            printf("rowA=%i <-> colB=%i\n", i1->first, i2->first);
        }
    }
}

std::vector<int> compDVals(std::map<int, std::map<int, int> > a, 
                           std::map<int, std::map<int, int> > b,
                           std::map<int, std::map<int, int> > &c,
                           std::map<int, int> setmap,
                           std::map<int, int> lockmap) 
{
    std::vector<int> d;
    int intsum;
    int extsum;

    //reset cutset
    cutset = 0;

    for (std::map<int, std::map<int, int> >::iterator i1 = a.begin(); 
            i1 != a.end(); i1++)
    {
        if (lockmap[i1->first] == UNLOCKED) {
            intsum = 0;
            extsum = 0;
            for (std::map<int, int>::iterator i2 = i1->second.begin(); 
                 i2 != i1->second.end(); i2++) 
            {
                if (setmap[i2->first] == SETA) {
                    intsum++;
                } else if (setmap[i2->first] == SETB) {
                    extsum++;
                    cutset++;

                    //construct c matrix of A<->B links
                    c[i1->first][i2->first] = 1;
                }
            }
            if(debug) printf("cell %i: %i int links; %i ext links; cost = %i\n", i1->first, intsum, extsum, extsum-intsum);
            d.push_back(extsum - intsum);
        } else {
            if (debug) printf("cell %i in set a is locked. skipping\n", i1->first);
        }
    }
    
    for (std::map<int, std::map<int, int> >::iterator i1 = b.begin(); 
            i1 != b.end(); i1++)
    {
        if (lockmap[i1->first] == UNLOCKED) {
            intsum = 0;
            extsum = 0;
            for (std::map<int, int>::iterator i2 = i1->second.begin(); 
                 i2 != i1->second.end(); i2++) 
            {
                if (setmap[i2->first] == SETB) {
                    intsum++;
                } else if (setmap[i2->first] == SETA) {
                    extsum++;
                    cutset++;

                    //construct c matrix of A<->B links
                    c[i2->first][i1->first] = 1;
                }
            }
            if(debug) printf("cell %i: %i int links; %i ext links; cost = %i\n", i1->first, intsum, extsum, extsum-intsum);
            d.push_back(extsum - intsum);
        } else {
            if (debug) printf("cell %i in set b is locked. skipping\n", i1->first);
        }
    }

    return d;
}

int main(int argc, char* argv[])
{
    int line;
    int linenum = 0;
    bool afull = false;
    std::map<int, std::map<int, int> > a, b, c;
    std::map<int, int> set_hashmap;
    std::map<int, int> lock_hashmap;
    int gmax = -1000;
    int count = 0;

    if (argc < 2) {
        std::cout << "Usage: klalgo <filename> [-d]" << std::endl;
        exit(0);
    } else if (argc == 3 && std::string(argv[2]) == "-d") {
        debug = true;
    }

    //get netlist, populate data structs
    std::ifstream infile(argv[1]);
    if (infile.is_open()) {
        while (1) {
            int tmp1, tmp2;
            infile >> tmp1;
            if (infile.eof()) {
                break;
            } else if (linenum == 0) {
                cellcount = tmp1;
            } else if (linenum == 1) {
                netcount = tmp1;
            } else if (!afull) {
                if (set_hashmap.count(tmp1) == 0) {
                    set_hashmap[tmp1] = SETA;
                    if (set_hashmap.size() == cellcount/2) { afull = true; }
                }
                infile >> tmp2;
                a[tmp1][tmp2] = 1;
            } else {
                if (set_hashmap.count(tmp1) == 0) {
                    set_hashmap[tmp1] = SETB;
                }                    
                infile >> tmp2;
                b[tmp1][tmp2] = 1;
            }
            linenum++;
        }
//        if (debug) {
//            printf("cell count: %i\n", cellcount); 
//            printf("net count: %i\n", netcount);
//            for (std::map<int, std::map<int, int> >::iterator i1 = a.begin(); 
//                    i1 != a.end(); i1++)
//            {
//                for (std::map<int, int>::iterator i2 = i1->second.begin(); 
//                     i2 != i1->second.end(); i2++) 
//                {
//                    printf("set a, cell %i connects to %i\n", i1->first, i2->first );
//                }
//            }
//            for (std::map<int, std::map<int, int> >::iterator i1 = b.begin(); 
//                    i1 != b.end(); i1++)
//            {
//                for (std::map<int, int>::iterator i2 = i1->second.begin(); 
//                     i2 != i1->second.end(); i2++) 
//                {
//                    printf("set b, cell %i connects to %i\n", i1->first, i2->first );
//                }
//            }
//        }
        infile.close();
    } else {
        std::cout << "Unable to open file" << std::endl;
    }

    do {
        //compute D values for all a in A1 and b in B1    
        std::vector<int> dvals = compDVals(a, b, c, set_hashmap, lock_hashmap);
        printf("Initial cutset = %i for iteration %i\n", cutset, count++);
        debug_print_matrices(c, "initial c matrix");
     
//        if (debug) {
//            debug_print_matrices(a, "a"); 
//            debug_print_matrices(b, "b"); 
//            debug_print_matrices(c, "c"); 
//        } 

        
        //clear locked vertices

        int g[cellcount/2];
        for (int n = 0; n < cellcount/2; n++) {
            int ai = 0;
            int bj = 0;
            int gnmax = -1000;
            
            //find a[i] from A1 and b[j] from B1, so g[n] = D[a[i]] + D[b[j]] - 2*c[a[i]][b[j]] is maximal
            for (std::map<int, std::map<int, int> >::iterator i1 = a.begin(); 
                    i1 != a.end(); i1++)
            {
                //check if a[i] was removed from consideration
                if (lock_hashmap.count(i1->first) == UNLOCKED) { 
                    for (std::map<int, std::map<int, int> >::iterator i2 = b.begin(); 
                        i2 != b.end(); i2++)
                    {
                        //check if b[j] was removed from consideration
                        if (lock_hashmap.count(i2->first) == UNLOCKED) { 
                            g[n] = dvals[i1->first] + dvals[i2->first] - (2 * c[i1->first][i2->first]);
                            if (g[n] > gnmax) {
                                gnmax = g[n];
                                ai = i1->first;
                                bj = i2->first;
                            } 
                        }
                    }
                }
            }

            if(debug) printf("found g[n] max of %i @ ai = %i, bi = %i\n", gnmax, ai, bj); 
                    
            //move a[i] to B1 and b[j] to A1
            std::map<int, std::map<int, int> > tmpa(a);
            std::map<int, std::map<int, int> > tmpb(b);
            a.erase(ai);
            b.erase(bj);
            a[bj].insert(tmpb[bj].begin(), tmpb[bj].end());
            b[ai].insert(tmpa[ai].begin(), tmpa[ai].end());
            //update set hashmap
            set_hashmap[ai] = SETB;
            set_hashmap[bj] = SETA;

            if(debug) printf("swapped a[i]=%i with b[j]=%i\n", ai, bj);

            //remove a[i] and b[j] from further consideration in this pass
            lock_hashmap[ai] = LOCKED;
            lock_hashmap[bj] = LOCKED;

            if(debug) printf("locking a[i]=%i and b[j]=%i\n", ai, bj);

            //update D values for the elements of A1 = A1 \ a[i] and B1 = B1 \ b[j]
            //clear matrix C first
            c.clear();
            dvals = compDVals(a, b, c, set_hashmap, lock_hashmap);

            debug_print_matrices(c, "c");
            printf("cutset = %i\n", cutset);
        }

        //find k which maximizes g_max, the sum of g[1],..., g[k]
        int finalk = 0;
        for (int k = 0; k < cellcount/2; k++) {
            int sum = 0;
            for (int i = 0; i < k; i++) {
                sum += g[i];

                if (sum > gmax) {
                    gmax = sum;
                    finalk = k;
                }
            }
        }
        printf("gmax = %i @ k = %i\n", gmax, finalk);

        if (gmax > 0) {
            //TODO: exchange a[1], a[2],...,a[k] with b[1], b[2],...b[k]
            std::map<int, std::map<int, int> > tmpa(a);
            std::map<int, std::map<int, int> > tmpb(b);

            std::map<int, std::map<int, int> >::iterator ia = a.begin();
            std::map<int, std::map<int, int> >::iterator ib = b.begin();
            for (int i = 0; i < finalk; i++) {
                if (ia == a.end() || ib == b.end()) {
                    printf("PROBABLY AN ERROR!\n");
                }

                a.erase(ia->first);
                a[ia->first].insert(tmpb[ib->first].begin(), tmpb[ib->first].end());
                b.erase(ib->first);
                b[ib->first].insert(tmpa[ia->first].begin(), tmpa[ia->first].end());
                ia++;
                ib++;
            }
        }
    } while (gmax > 0);

    //print final partitioning
    int buf_size = cellcount * ((cellcount/2)+1);
    char seta_buf[buf_size];
    memset(seta_buf, 0, buf_size); 
    char setb_buf[cellcount * ((cellcount/2)+1)];
    memset(setb_buf, 0, buf_size); 
    for (std::map<int, int>::iterator i = set_hashmap.begin(); i != set_hashmap.end(); i++) {
        if (set_hashmap[i->first] == SETA) {
            sprintf(seta_buf + std::strlen(seta_buf), "%i ", i->first); 
        } else {
            sprintf(setb_buf + std::strlen(setb_buf), "%i ", i->first); 
        }
    }

    std::cout << "set a: " << seta_buf << std::endl;
    std::cout << "set b: " << setb_buf << std::endl;
    std::cout << "final cutset: " << cutset << std::endl;

    return 0;
}


