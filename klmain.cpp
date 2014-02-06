#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <cstring>
#include <assert.h>
#include <ctime>
#include <pthread.h>
#include <iterator>
#include <math.h>

//global variables
static bool debug = false;
static int cellcount = 0;
static int netcount = 0;
static std::map<int, std::map<int, int> > c;
static std::map<int, int> lock_hashmap;

enum {
    SETA = 0,  
    SETB = 1
}; 

enum {
    UNLOCKED = 0,
    LOCKED
};

typedef struct {
    int thrId;
    std::set<int> a, b;
    std::map<int, int> dvals;
    int numcores;
} threadParamInfo;

typedef struct {
    int ai;
    int bj;
    int gmax;
} gmaxReturnInfo;

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

void verify_set_constrs(std::set<int> a, std::set<int> b)
{
    if(debug) printf("sizeof(a)=%i\tsizeof(b)=%i\n", a.size(), b.size());
    assert(a.size() == b.size());
    assert(a.size() + b.size() == cellcount);
}

std::map<int, int> compDVals(std::set<int> a, std::set<int> b) 
{
    std::map<int, int> d;
    int intsum;
    int extsum;
   
    for (int i1 = 1; i1 < cellcount; i1++) {
        if (lock_hashmap[i1] == UNLOCKED) {
            intsum = 0;
            extsum = 0;

            for (std::map<int, int>::iterator i2 = c[i1].begin(); 
                 i2 != c[i1].end(); i2++)
            {
                if (c[i1].count(i2->first) > 0){
                    if (a.count(i1) && b.count(i2->first) ||
                        b.count(i1) && a.count(i2->first)) 
                    {
                        extsum++;
                    } else {
                        intsum++;
                    }
                }
            }
            if(debug) printf("cell %i has %i internal links and %i external links\n", i1, intsum, extsum);
        } 
        d[i1] = (extsum - intsum);
    }

    return d;
}

int compute_cutset(std::set<int> a, std::set<int> b, std::map<int, std::map<int, int> > c)
{
    int cur_cutset = 0; 

    for (int i1 = 1; i1 <= cellcount; i1++) {
        for (std::map<int, int>::iterator i2 = c[i1].begin(); 
                 i2 != c[i1].end(); i2++)
        {
            if (i2->first > i1) { 
                if (a.count(i1) && b.count(i2->first) ||
                    b.count(i1) && a.count(i2->first)) 
                {
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

gmaxReturnInfo *compGmax(std::set<int> a, std::set<int> b, std::map<int, int> dvals)
{
    gmaxReturnInfo *ret = (gmaxReturnInfo*)malloc(sizeof(gmaxReturnInfo));
    if (ret == NULL) {
        printf("MALLOC ERROR\n");
        exit(1);
    }
    int gmax = -1000;
    int g = 0;

    for (std::set<int>::iterator i1 = a.begin(); i1 != a.end(); i1++) {
        if (lock_hashmap[*i1] != LOCKED) { 
            for (std::set<int>::iterator i2 = b.begin(); i2 != b.end(); i2++) {
                if (lock_hashmap[*i2] != LOCKED)
                {
                    g = dvals[*i1] + dvals[*i2] - (2 * c[*i1].count(*i2));

                    if(debug) {
                        printf("\tg_%i_%i = D_%i + D_%i - 2c_%i_%i = %i + %i - 2(%i) = %i\n",
                                *i1, *i2, *i1, *i2, *i1, *i2,
                                dvals[*i1], dvals[*i2], c[*i1].count(*i2),
                                g);
                    }
                    if (g >= gmax) {
                        gmax = g;
                        ret->gmax = g;
                        ret->ai = *i1;
                        ret->bj = *i2;
                    }
                }
            }
        }
    }
    return ret;
}

void *thr_compG(void *arg) {
    threadParamInfo *info = (threadParamInfo*)arg;
    int count = 0;
    std::set<int> suba;
    int subcellcount = floor(info->a.size()/info->numcores);

    //redefine set A based on input params
    std::set<int>::iterator ia = info->a.begin();
    //advance to starting cell
    std::advance(ia, subcellcount*info->thrId);
    if (debug) printf("creating subset of a with %i elements starting with element %i\n", subcellcount, *ia);
    for (std::set<int>::iterator i = ia; i != info->a.end(); i++) {
        if (count++ < subcellcount) {
            suba.insert(*i);
        } else {
            break;
        }
    }

    gmaxReturnInfo *ret = compGmax(suba, info->b, info->dvals);

    if(debug) printf("thread%i found g=%i @ ai = %i, bi = %i\n", info->thrId, ret->gmax, ret->ai, ret->bj); 

    return (void*) ret; 
}

int main(int argc, char* argv[])
{
    int line;
    int linenum = 0;
    std::set<int> a, b;
    int cutset = 0;
    int gmax = -1000;
    int count = 0;

    int start_s = clock(); 

    if (argc < 1) {
        std::cout << "Usage: klalgo [-d]" << std::endl;
        exit(0);
    } else if (argc == 2 && std::string(argv[1]) == "-d") {
        debug = true;
    }

    //get netlist, populate data structs
    int tmp1, tmp2;
    while (std::cin >> tmp1) {
        if (linenum == 0) {
            cellcount = tmp1;
        } else if (linenum == 1) {
            netcount = tmp1;
        } else {
            std::cin >> tmp2;
            if (a.count(tmp1) == 0 && b.count(tmp1) == 0) {
                if (a.size() < cellcount/2) {
                    if(debug) printf("inserting cell %i into set a\n", tmp1);
                    a.insert(tmp1); 
                } else { 
                    if(debug) printf("inserting cell %i into set b\n", tmp1);
                    b.insert(tmp1);
                }
            }
            
            c[tmp1][tmp2] = 1;
            c[tmp2][tmp1] = 1; 
        } 
        linenum++;
        if(debug) printf("reading file @ line %i\n", linenum);
    }
    
    //fill in unconnected nodes
    if (b.size() != a.size()) {
        for (int i = 1; i <= cellcount; i++ ) {
            if (a.count(i) == 0 && b.count(i) == 0) {
                b.insert(i); 
            }
        }
    }

    print_matrix(c);

    do {
        //make sure initial coniditions are satisfied
        verify_set_constrs(a, b); 

        //compute D values for all a in A1 and b in B1    
        std::map<int, int> dvals = compDVals(a, b);
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

            if(debug) print_current_sets(a, b);
            if (debug) print_matrix(c);
            cutset = compute_cutset(a, b, c);
            if(debug) printf("cutset = %i\n", cutset);

            //find a[i] from A1 and b[j] from B1, so g[n] = D[a[i]] + D[b[j]] - 2*c[a[i]][b[j]] is maximal
            if (CORES > 1) {
                threadParamInfo paramInfo;
                pthread_t thr[CORES];

                std::set<int>::iterator i1 = a.begin();
                paramInfo.a = a;
                paramInfo.b = b;
                paramInfo.dvals = dvals;
                paramInfo.numcores = CORES;
                for (int i = 0; i < CORES; i++) {
                    printf("spawning child thread %i to compute g for first 1/%i of set a cells\n", i, CORES);
                    paramInfo.thrId = i;

                    pthread_create(&thr[i], NULL, &thr_compG, (void*)&paramInfo);
                }
                for (int i = 0; i < CORES; i++) {
                    void *p;
                    pthread_join(thr[i], &p);
                    gmaxReturnInfo *ret = static_cast<gmaxReturnInfo *>(p);

                    //compare gmax returned from threads
                    if (ret->gmax > gnmax) {
                        gnmax = ret->gmax;
                        ai[n] = ret->ai;
                        bj[n] = ret->bj;
                    }
                }
            } else {
                gmaxReturnInfo *ret = compGmax(a, b, dvals);
                ai[n] = ret->ai;
                bj[n] = ret->bj;
                gnmax = ret->gmax;
            }

            assert (ai[n] != 0 && bj[n] != 0);
            g[n] = gnmax;

            if(debug) printf("1) found g[%i]=%i @ ai = %i, bi = %i\n", n, gnmax, ai[n], bj[n]); 
                    
            //move a[i] to B1 and b[j] to A1
            a.erase(ai[n]);
            a.insert(bj[n]);
            b.erase(bj[n]);
            b.insert(ai[n]);

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

                if (a.count(i1->first) > 0) {
                    dvals[i1->first] = old_dval + 2*(c_ai_curnode - c_bj_curnode);
                    if(debug) {
                        printf("\tD_%i' = D_%i + 2(c_%i_%i - c_%i_%i) = %i + 2(%i - %i) = %i\n",
                               i1->first, i1->first, i1->first, ai[n], i1->first, bj[n], old_dval, c_ai_curnode, 
                               c_bj_curnode, dvals[i1->first]); 
                    }
                } else if (b.count(i1->first) > 0) {
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
            if(debug) print_current_sets(a, b);
            if(debug) printf("5) swapping a[0],...a[%i] and b[0],...,b[%i]\n", finalk, finalk);

            //exchange a[1], a[2],...,a[k] with b[1], b[2],...b[k]
            for (int i = 0; i <= finalk; i++) {
                if(debug) printf("swapping %i with %i\n", ai[i], bj[i]);
                a.erase(bj[i]);
                a.insert(ai[i]);
                b.erase(ai[i]);
                b.insert(bj[i]);
            }
        }
    } while (gmax > 0);

    //make sure set constraints are maintained
    verify_set_constrs(a, b);

    //print final partitioning
    cutset = compute_cutset(a, b, c);
    std::cout << "final cutset: " << cutset << std::endl;
    print_current_sets(a, b);

    int stop_s = clock();
    std::cout << "time: " << (stop_s - start_s)/double(CLOCKS_PER_SEC)*1000  << " ms" << std::endl;

    return 0;
}


