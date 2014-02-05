#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cstring>
#include <assert.h>
#include <ctime>

//global variables
bool debug = false;
int cellcount = 0;
int netcount = 0;

enum {
    SETA = 0,  
    SETB = 1
}; 

enum {
    UNLOCKED = 0,
    LOCKED
};

void print_matrix(std::map<int, std::map<int, int> > mat)
{
    int element = 0;
    for (std::map<int, std::map<int, int> >::iterator i1 = mat.begin(); 
            i1 != mat.end(); i1++)
    {
        element++;
        printf("element[%2i]: %2i: ", element, i1->first); 
        for (std::map<int, int>::iterator i2 = i1->second.begin(); 
             i2 != i1->second.end(); i2++) 
        {
            printf("%2i ", i2->first);
        }
        printf("\n");
    }
}

void verify_set_constrs(std::map<int, int> setmap)
{
    int sizea = 0;
    int sizeb = 0;
    int numcells = 0;

    for (std::map<int, int>::iterator i3 = setmap.begin(); i3 != setmap.end(); i3++) {
        numcells++;
        if (setmap[i3->first] == SETA) {
            sizea++;
        } else {
            sizeb++;
        }
    }

    assert(sizea == sizeb);
    assert(numcells == cellcount);
}

std::map<int, int> compDVals(std::map<int, std::map<int, int> > c,
                           std::map<int, int> setmap,
                           std::map<int, int> lockmap) 
{
    std::map<int, int> d;
    int intsum;
    int extsum;
   
    for (std::map<int, int>::iterator i1 = setmap.begin(); i1 != setmap.find(cellcount); i1++) {
        if (lockmap[i1->first] == UNLOCKED) {
            intsum = 0;
            extsum = 0;

            for (std::map<int, int>::iterator i2 = c[i1->first].begin(); 
                 i2 != c[i1->first].end(); i2++)
            {
                if (c[i1->first].count(i2->first) > 0){
                    if (setmap[i1->first] != setmap[i2->first]) 
                    {
                        extsum++;
                    } else {
                        intsum++;
                    }
                }
            }
            if(debug) printf("cell %i has %i internal links and %i external links\n", i1->first, intsum, extsum);
        } 
        d[i1->first] = (extsum - intsum);
    }

    return d;
}

int compute_cutset(std::map<int, int> setmap, std::map<int, std::map<int, int> > c)
{
    int cur_cutset = 0; 

    for (std::map<int, int>::iterator i1 = setmap.begin(); i1 != setmap.end(); i1++) {
        for (std::map<int, int>::iterator i2 = c[i1->first].begin(); 
                 i2 != c[i1->first].end(); i2++)
        {
            if (i2->first > i1->first) { 
                if (setmap[i1->first] != setmap[i2->first]) {
                    cur_cutset++;
                } 
            }
        }
    }
    return cur_cutset;
}


void print_current_sets(std::map<int, int> sets)
{
    int a[cellcount/2];
    int b[cellcount/2];
    int acount = 0;
    int bcount = 0;

    for (std::map<int, int>::iterator i1 = sets.begin(); i1 != sets.end(); i1++) {
        if (sets[i1->first] == SETA) {
            a[acount++] = i1->first;
        } else {
            b[bcount++] = i1->first;
        }
    }

    printf("set a\t\tset b\n");
    for (int i = 0; i < cellcount/2; i++) {
        printf("%5i\t\t%5i\n", a[i], b[i]); 
    }
}

int main(int argc, char* argv[])
{
    int line;
    int linenum = 0;
    std::map<int, std::map<int, int> > c;
    std::map<int, int> set_hashmap, lock_hashmap;
    int cutset = 0;
    int gmax = -1000;
    int count = 0;

    int start_s = clock(); 

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
            } else {
                infile >> tmp2;
                if (set_hashmap.count(tmp1) == 0) {
                    if (set_hashmap.size() < cellcount/2) {
                        if(debug) printf("inserting cell %i into set a\n", tmp1);
                        set_hashmap[tmp1] = SETA;
                    } else { 
                        if(debug) printf("inserting cell %i into set b\n", tmp1);
                        set_hashmap[tmp1] = SETB;
                    }
                }
                
                c[tmp1][tmp2] = 1;
                c[tmp2][tmp1] = 1; 
            } 
            linenum++;
            if(debug) printf("reading file @ line %i\n", linenum);
        }
        infile.close();
    } else {
        std::cout << "Unable to open file" << std::endl;
        return 0;
    }
    
    //fill in unconnected nodes
    if (set_hashmap.size() != cellcount) {
        for (int i = 1; i <= cellcount; i++ ) {
            if (set_hashmap.count(i) == 0) {
                set_hashmap[i] = SETB; 
            }
        }
    }

    do {
        //make sure initial coniditions are satisfied
        verify_set_constrs(set_hashmap); 

        //compute D values for all a in A1 and b in B1    
        std::map<int, int> dvals = compDVals(c, set_hashmap, lock_hashmap);
        if(debug) printf("**********************************************\n");
     
        int gmax = -1000;
        int finalk = 0;
        int sum = 0;
        int g[cellcount/2];
        int ai[cellcount/2];
        int bj[cellcount/2];
        for (int n = 0; n < (cellcount/2); n++) {
            int stop_s = clock();
            std::cout << "elapsed time: " << (stop_s - start_s)/double(CLOCKS_PER_SEC)*1000  << " ms" << std::endl;

            int gnmax = -1000;
            
            if(debug) printf("---- ITER %i ----\n", n);

            if(debug) print_current_sets(set_hashmap);
            if (debug) print_matrix(c);
            cutset = compute_cutset(set_hashmap, c);
            if(debug) printf("cutset = %i\n", cutset);

            //find a[i] from A1 and b[j] from B1, so g[n] = D[a[i]] + D[b[j]] - 2*c[a[i]][b[j]] is maximal
            for (std::map<int, int>::iterator i1 = set_hashmap.begin(); i1 != set_hashmap.end(); i1++) {
                if (set_hashmap[i1->first] == SETA && lock_hashmap[i1->first] != LOCKED) { 
                    for (std::map<int, int>::iterator i2 = set_hashmap.begin(); i2 != set_hashmap.end(); i2++) {
                        if (set_hashmap[i2->first] == SETB && lock_hashmap[i2->first] != LOCKED)
                        {
                            g[n] = dvals[i1->first] + dvals[i2->first] - (2 * c[i1->first].count(i2->first));

                            if(debug) {
                                printf("\tg_%i_%i = D_%i + D_%i - 2c_%i_%i = %i + %i - 2(%i) = %i\n",
                                        i1->first, i2->first, i1->first, i2->first, i1->first, i2->first,
                                        dvals[i1->first], dvals[i2->first], c[i1->first].count(i2->first),
                                        g[n]);
                            }
                            if (g[n] >= gnmax) {
                                gnmax = g[n];
                                ai[n] = i1->first;
                                bj[n] = i2->first;
                            }
                        }
                    }
                }
            }

            assert (ai[n] != 0 && bj[n] != 0);
            g[n] = gnmax;

            if(debug) printf("1) found g[%i]=%i @ ai = %i, bi = %i\n", n, gnmax, ai[n], bj[n]); 
                    
            //move a[i] to B1 and b[j] to A1
            set_hashmap[ai[n]] = SETB;
            set_hashmap[bj[n]] = SETA;

            if(debug) printf("2) swapped a[i]=%i with b[j]=%i\n", ai[n], bj[n]);

            //remove a[i] and b[j] from further consideration in this pass
            lock_hashmap[ai[n]] = LOCKED;
            lock_hashmap[bj[n]] = LOCKED;
            dvals.erase(ai[n]);
            dvals.erase(bj[n]);

            if(debug) printf("3) locking a[i]=%i and b[j]=%i\n", ai[n], bj[n]);
            
            //find k which maximizes gmax, the sum of g[0],...,g[k]
            sum += g[n];
            if (sum > gmax) {
                gmax = sum;
                finalk = n;
            }
            if(debug) printf("sum of g[0...%i] = %i\n", n, sum);

            //update D values for the elements of A1 = A1 \ a[i] and B1 = B1 \ b[j]
            if (n == (cellcount/2)-1) { continue; }
            if (debug) printf("recomputing dvals\n");
            for (std::map<int, int>::iterator i1 = dvals.begin(); i1 != dvals.end(); i1++) {
                //if connected to locked vertices
                if (c[ai[n]].count(i1->first) == 0 && c[bj[n]].count(i1->first) == 0) {
                    continue;
                }

                int old_dval = dvals[i1->first];
                int c_ai_curnode, c_bj_curnode = 0;
                if (c[i1->first].count(ai[n]) > 0) {
                    c_ai_curnode = c[i1->first][ai[n]];
                } else {
                    c_ai_curnode = 0;
                }
                if (c[i1->first].count(bj[n]) > 0) {
                    c_bj_curnode = c[i1->first][bj[n]];
                } else {
                    c_bj_curnode = 0;
                }

                if (set_hashmap[i1->first] == SETA) {
                    dvals[i1->first] = old_dval + 2*(c_ai_curnode - c_bj_curnode);
                    if(debug) {
                        printf("\tD_%i' = D_%i + 2(c_%i_%i - c_%i_%i) = %i + 2(%i - %i) = %i\n",
                               i1->first, i1->first, i1->first, ai[n], i1->first, bj[n], old_dval, c_ai_curnode, 
                               c_bj_curnode, dvals[i1->first]); 
                    }
                } else if (set_hashmap[i1->first] == SETB) {
                    dvals[i1->first] = old_dval + 2*(c_bj_curnode - c_ai_curnode);
                    if(debug) {
                        printf("\tD_%i' = D_%i + 2(c_%i_%i - c_%i_%i) = %i + 2(%i - %i) = %i\n",
                               i1->first, i1->first, i1->first, bj[n], i1->first, ai[n], old_dval, c_bj_curnode, 
                               c_ai_curnode, dvals[i1->first]); 
                    }
                }
            }
        }

        //reset all cell locks
        for (std::map<int, int>::iterator i = lock_hashmap.begin(); i != lock_hashmap.end(); i++) {
            lock_hashmap[i->first] = UNLOCKED;
        }

        if(debug) printf("4) gmax = %i @ k = %i\n", gmax, finalk);

        if (gmax > 0) {
            if(debug) printf("sets before swapping k cells\n");
            if(debug) print_current_sets(set_hashmap);
            if(debug) printf("5) swapping a[0],...a[%i] and b[0],...,b[%i]\n", finalk, finalk);

            //exchange a[1], a[2],...,a[k] with b[1], b[2],...b[k]
            for (int i = 0; i <= finalk; i++) {
                if(debug) printf("swapping %i with %i\n", ai[i], bj[i]);
                //cell was originally in A, so switch it back to A
                set_hashmap[ai[i]] = SETA;
                //cell was originally in B, so switch it back to B
                set_hashmap[bj[i]] = SETB;
            }
        }
    } while (gmax > 0);

    //make sure set constraints are maintained
    verify_set_constrs(set_hashmap);

    //print final partitioning
    print_current_sets(set_hashmap);
    cutset = compute_cutset(set_hashmap, c);
    std::cout << "final cutset: " << cutset << std::endl;

    int stop_s = clock();
    std::cout << "time: " << (stop_s - start_s)/double(CLOCKS_PER_SEC)*1000  << " ms" << std::endl;

    return 0;
}


