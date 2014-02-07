#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <cstring>
#include <assert.h>
#include <pthread.h>
#include <iterator>
#include <vector>
#include <math.h>

//global variables
static bool debug = false;
static int cellcount = 0;
static int netcount = 0;
static int **c;
static int *numlinks;
static std::map<int, int> lock_hashmap;

enum {
    SETA = 1,  
    SETB = 2
}; 

enum {
    UNLOCKED = 0,
    LOCKED
};

struct costInfo{
    int cost;
    int cell;
    int set;
    friend bool operator<(costInfo const& lhs, costInfo const& rhs) { return lhs.cost < rhs.cost; }
};

void print_matrix(int **mat)
{
    for (int i1 = 1; i1 <= cellcount; i1++) {
        printf("%i: ", i1);
        for (int i2 = 0; i2 < numlinks[i1]; i2++) {
            printf("%2i ", c[i1][i2]);
        }
        printf("\n");
    }
}

int searchConnectionsForElement(int srccell, int dstcell, int min, int max) 
{
    if (max < min) {
        return -1;
    } else {
        int mid = float(max+min)/2.0;
        
        if (dstcell < c[srccell][mid]) {
            return searchConnectionsForElement(srccell, dstcell, min, mid-1);
        } else if (dstcell > c[srccell][mid]) {
            return searchConnectionsForElement(srccell, dstcell, mid+1, max);
        } else {
            return mid;
        }
    }
}

void verify_set_constrs(std::set<int> a, std::set<int> b)
{
    if(debug) printf("sizeof(a)=%i\tsizeof(b)=%i\n", int(a.size()), int(b.size()));
    assert(a.size() == b.size());
    assert(int(a.size() + b.size()) == cellcount);
}

std::multiset<struct costInfo> compDVals(int *setarray, int curset) 
{
    std::multiset<struct costInfo> ret;
    int intsum;
    int extsum;
   
    for (int i1 = 1; i1 <= cellcount; i1++) {
        if (setarray[i1] == curset && lock_hashmap[i1] == UNLOCKED) {
            intsum = 0;
            extsum = 0;

            for (int i2 = 0; i2 < numlinks[i1]; i2++) 
            {
                if (setarray[i1] == setarray[c[i1][i2]]) {
                    intsum++;
                } else {
                    extsum++;
                }
            }
            if(debug) printf("cell %i has %i internal links and %i external links\n", i1, intsum, extsum);

            struct costInfo cellinfo;
            cellinfo.cost = extsum - intsum;
            cellinfo.cell = i1;
            cellinfo.set = curset;
            ret.insert(cellinfo);
        } 
    }

    return ret;
}

int compute_cutset(int *setarray)
{
    int cur_cutset = 0; 

    for (int i1 = 1; i1 <= cellcount; i1++) {
        for (int i2 = 0; i2 < numlinks[i1]; i2++) 
        {
            if (c[i1][i2] > i1) { 
                if (setarray[i1] != setarray[c[i1][i2]]) {
                    cur_cutset++;
                } 
            }
        }
    }
    return cur_cutset;
}

void print_current_sets(std::set<int> a, std::set<int> b) {

    for (std::set<int>::iterator i = a.begin(); i != a.end(); i++) {
        printf("%i ", *i);
    }
    printf("\n");
    for (std::set<int>::iterator i = b.begin(); i != b.end(); i++) {
        printf("%i ", *i);
    }
    printf("\n");
}

int main(int argc, char* argv[])
{
    int linenum = 0;
    int *setarray = NULL;
    std::set<int> a, b;
    int cutset = 0;
    int gmax = -1000;
    int generation = 0;

    if (argc < 1) {
        std::cout << "Usage: klalgo [-d]" << std::endl;
        exit(0);
    } else if (argc == 2 && std::string(argv[1]) == "-d") {
        debug = true;
    }

    //get netlist, populate data structs
    int tmp1, tmp2;
    std::vector<std::set<int> > links;
    while (std::cin >> tmp1) {
        if (linenum == 0) {
            cellcount = tmp1;
            setarray = new int[cellcount+1]; //ind 0 will never be used
            memset(setarray, 0, sizeof(int)*(cellcount+1));

            //experimenting with 2d array data struct
            links.resize(cellcount+1);
            numlinks = new int[cellcount+1];
            c = new int*[cellcount+1];
        } else if (linenum == 1) {
            netcount = tmp1;
        } else {
            std::cin >> tmp2;
            if (!setarray[tmp1]) {
                if (int(a.size()) < cellcount/2) {
                    if(debug) printf("inserting cell %i into set a\n", tmp1);
                    a.insert(tmp1); 
                    setarray[tmp1] = SETA;
                } else { 
                    if(debug) printf("inserting cell %i into set b\n", tmp1);
                    b.insert(tmp1);
                    setarray[tmp1] = SETB;
                }
            }
            links[tmp1].insert(tmp2);
            links[tmp2].insert(tmp1);
        } 
        linenum++;
    }

    for (int i1 = 1; i1 <= cellcount; i1++) {
        int ind = 0;
        c[i1] = new int[links[i1].size()];
        numlinks[i1] = links[i1].size();
        for (std::set<int>::iterator i2 = links[i1].begin(); i2 != links[i1].end(); i2++) {
            c[i1][ind++] = *i2; 
        }
    }
    
    //fill in unconnected nodes
    if (b.size() != a.size()) {
        for (int i = 1; i <= cellcount; i++ ) {
            if (!setarray[i]) {
                b.insert(i); 
                setarray[i] = SETB;
            }
        }
    }

    links.clear(); 

    do {
        //make sure initial coniditions are satisfied
        verify_set_constrs(a, b); 

        //compute D values for all a in A1 and b in B1    
        std::multiset<struct costInfo> da = compDVals(setarray, SETA);
        std::multiset<struct costInfo> db = compDVals(setarray, SETB);
        
        if(debug) printf("**********************************************\n");
     
        int gmax = -1000;
        int finalk = 0;
        int sum = 0;
        int g[cellcount/2];
        int ai[cellcount/2];
        int bj[cellcount/2];
        for (int n = 0; n < (cellcount/2); n++) {
            int gnmax = -1000;
            
            printf("---- GENERATION: %i ITERATION: %i ----\n", generation, n);

            if(debug) print_current_sets(a, b);
            if(debug) print_matrix(c);
            cutset = compute_cutset(setarray);
            if(debug) printf("cutset = %i\n", cutset);

            //find a[i] from A1 and b[j] from B1, so g[n] = D[a[i]] + D[b[j]] - 2*c[a[i]][b[j]] is maximal
            //only compute for the first n largest D vals in each set
            int acomps = 0;
            for (std::multiset<struct costInfo>::reverse_iterator ia = da.rbegin(); ia != da.rend(); ia++) {
                if (lock_hashmap[ia->cell] == UNLOCKED) {
                    int bcomps = 0;
                    for (std::multiset<struct costInfo>::reverse_iterator ib = db.rbegin(); ib != db.rend(); ib++) {
                        if (lock_hashmap[ib->cell] == UNLOCKED) {
                            int g = 0;
                            int c_a_b = 0;
                            int cellind = searchConnectionsForElement(ia->cell, ib->cell, 0, numlinks[ia->cell]);
                            if (cellind >= 0) {
                                c_a_b = 1;
                            }
                            
                            g = ia->cost + ib->cost - (2 * c_a_b);

                            if(debug) {
                                printf("\tg_%i_%i = D_%i + D_%i - 2c_%i_%i = %i + %i - 2(%i) = %i\n",
                                        ia->cell, ib->cell, ia->cell, ib->cell, ia->cell, ib->cell,
                                        ia->cost, ib->cost, c_a_b, g);
                            }
                            if (g > gnmax) {
                                gnmax = g;
                                ai[n] = ia->cell;
                                bj[n] = ib->cell;
                            }

                            if (bcomps++ > int(db.size()/10)) { break; }
                        }  
                    }

                    if (acomps++ > int(da.size()/10)) { break; }
                }
            }
          
            assert (ai[n] != 0 && bj[n] != 0);
            g[n] = gnmax;

            if(debug) printf("1) found g[%i]=%i @ ai = %i, bi = %i\n", 
                             n, gnmax, ai[n], bj[n]); 
                    
            //move a[i] to B1 and b[j] to A1
            a.erase(ai[n]);
            a.insert(bj[n]);
            b.erase(bj[n]);
            b.insert(ai[n]);
            setarray[ai[n]] = SETB;
            setarray[bj[n]] = SETA;

            if(debug) printf("2) swapped a[i]=%i with b[j]=%i\n", ai[n], bj[n]);

            //remove a[i] and b[j] from further consideration in this pass
            lock_hashmap[ai[n]] = LOCKED;
            lock_hashmap[bj[n]] = LOCKED;

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
            //recompute set A costs
            std::multiset<struct costInfo> tmpa;
            for (std::multiset<struct costInfo>::iterator ia = da.begin(); ia != da.end(); ia++) {
                int new_dval = 0;
                struct costInfo tmpinfo;
                //if not connected to locked vertices, just copy old dvalue
                if (searchConnectionsForElement(ia->cell, ai[n], 0, numlinks[ia->cell]) < 0 &&
                    searchConnectionsForElement(ia->cell, bj[n], 0, numlinks[ia->cell]) < 0) 
                {
                    new_dval = ia->cost;        
                } else { 
                    int old_dval = ia->cost;
                    int c_ai_curnode, c_bj_curnode = 0;
                    int ai_ind = searchConnectionsForElement(ia->cell, ai[n], 0, numlinks[ia->cell]);
                    if (ai_ind >= 0) {
                        c_ai_curnode = 1;
                    }
                    int bj_ind = searchConnectionsForElement(ia->cell, bj[n], 0, numlinks[ia->cell]);
                    if (bj_ind >= 0) {
                        c_bj_curnode = 1;
                    }
                    
                    new_dval = old_dval + 2*(c_ai_curnode - c_bj_curnode);
                    if(debug) {
                        printf("\tD_%i' = D_%i + 2(c_%i_%i - c_%i_%i) = %i + 2(%i - %i) = %i\n",
                               ia->cell, ia->cell, ia->cell, ai[n], ia->cell, bj[n], old_dval, c_ai_curnode, 
                               c_bj_curnode, ia->cost); 
                    }
                }

                tmpinfo.cost = new_dval;
                tmpinfo.cell = ia->cell;
                tmpinfo.set = SETA;
                tmpa.insert(tmpinfo);
            }
            //repeat for set b
            std::multiset<struct costInfo> tmpb;
            for (std::multiset<struct costInfo>::iterator ib = db.begin(); ib != db.end(); ib++) {
                int new_dval = 0;
                struct costInfo tmpinfo;
                //if connected to locked vertices
                if (searchConnectionsForElement(ib->cell, ai[n], 0, numlinks[ib->cell]) < 0 &&
                    searchConnectionsForElement(ib->cell, bj[n], 0, numlinks[ib->cell]) < 0) 
                {
                    new_dval = ib->cost;
                } else {
                    int old_dval = ib->cost;
                    int c_ai_curnode, c_bj_curnode = 0;
                    int ai_ind = searchConnectionsForElement(ib->cell, ai[n], 0, numlinks[ib->cell]);
                    if (ai_ind >= 0) {
                        c_ai_curnode = 1;
                    }
                    int bj_ind = searchConnectionsForElement(ib->cell, bj[n], 0, numlinks[ib->cell]);
                    if (bj_ind >= 0) {
                        c_bj_curnode = 1;
                    }
                    
                    new_dval = old_dval + 2*(c_bj_curnode - c_ai_curnode);
                    if(debug) {
                        printf("\tD_%i' = D_%i + 2(c_%i_%i - c_%i_%i) = %i + 2(%i - %i) = %i\n",
                               ib->cell, ib->cell, ib->cell, bj[n], ib->cell, ai[n], old_dval, c_bj_curnode, 
                               c_ai_curnode, ib->cost); 
                    }
                }

                tmpinfo.cost = new_dval;
                tmpinfo.cell = ib->cell;
                tmpinfo.set = SETB;
                tmpb.insert(tmpinfo);
            }

            da.clear();
            db.clear();
            da = tmpa;
            db = tmpb;
        }

        //reset all cell locks
        for (std::map<int, int>::iterator i = lock_hashmap.begin(); i != lock_hashmap.end(); i++) {
            lock_hashmap[i->first] = UNLOCKED;
        }

        if(debug) printf("4) gmax = %i @ k = %i\n", gmax, finalk);

        if (gmax > 0) {
            if(debug) printf("sets before swapping k cells\n");
            if(debug) print_current_sets(a, b);
            if(debug) printf("5) swapping a[0],...a[%i] and b[0],...,b[%i]\n", finalk, finalk);

            //exchange a[1], a[2],...,a[k] with b[1], b[2],...b[k]
            for (int i = 0; i <= finalk; i++) {
                if(debug) printf("swapping %i with %i\n", ai[i], bj[i]);
                a.erase(bj[i]);
                a.insert(ai[i]);
                b.erase(ai[i]);
                b.insert(bj[i]);
                setarray[ai[i]] = SETA;
                setarray[bj[i]] = SETB;
            }
        }

        generation++;
    } while (gmax > 0);

    //make sure set constraints are maintained
    verify_set_constrs(a, b);

    //print final partitioning
    cutset = compute_cutset(setarray);
    std::cout << cutset << std::endl;
    print_current_sets(a, b);

    return 0;
}
