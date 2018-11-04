#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>
#include <time.h>
#include <signal.h>
#include <memory>
#include "minisat/core/SolverTypes.h"
#include "minisat/core/Solver.h"


sem_t io_se1, io_se2, io_se3, vc_se1, vc_se2, vc_se3;

//double total_time1, total_time2, total_time3;
bool io_exit = 0;
const int  maxv = 999;
const int max_vertex_num = 20;
int num = -1;
int edgs[max_vertex_num][max_vertex_num];
std::vector <int> v;
std::vector <int> vertex_cover1;
std::vector <int> vertex_cover2;
std::vector <int> vertex_cover3;

bool isnum(char c);

void *io(void *arg);
void *minisat(void *arg);
void *vc1(void *arg);
void *vc2(void *arg);


bool isnum(char c){
    return (c >= '0' && c <= '9');
}

void* io(void *arg){
    while(!std::cin.eof()){
        char comma;
        int err = 0;
        static int err1 = 0;
        static int err2 = 0;
        std::cin.clear();
        std::cin.sync();
        std::cin.get(comma);
        
        switch(comma){
            case 'V' :{
                
                //input
                
                // read the amount of vertices
                if (err2==4){
                    std::cerr << "Error : repeative V command" << std::endl;
                    break;
                }
                
                std::cin >> num;
                
                if (num<0){
                    std::cerr << "Error : no vertex input" << std::endl;
                    break;
                }
                
//                 if (num>20){
//                     std::cerr << "Error : no more than 20 vertices" << std::endl;
//                     break;
//                 }
                
                // reset graph matrix
                
                for (int i=0;i<max_vertex_num;i++){
                    for (int j=0;j<max_vertex_num;j++){
                        edgs[i][j] = maxv;
                    }
                }
                
                
                
                err1 = 1;       //successful V command
                err2 = 4;       //first successful V comand
                break;
            }
                
            case 'E' : {
                
                std::string line;
                std::cin >> line;
                int k = line.size();
                int *e;
                e = new int[k];
                int result;
                int n = 0, m = 0;
                
                if (err1==3){
                    std::cerr << "Error : repeative E command" << std::endl;
                    break;
                }
                
                if (err1!=1){
                    std::cerr << "Error : no V command input" << std::endl;
                    break;
                }
                
                while (n<k){
                    result = 0;
                    
                    if (isnum(line[n])){
                        
                        if (line[n-1]=='-'){
                            err = 2;
                            break;
                        }
                        result = line[n] - '0';
                        
                        while (n<k && isnum(line[++n])){
                            result = result * 10 + line[n] - '0';
                        }
                        
                        if (result >= num || result < 0){
                            err = 1;
                            break;
                        }
                        
                        else
                            e[m++] = result;
                    }
                    ++n;
                }
                
                if (err==1){
                    std::cerr << "Error: vertex doesn't exist" << std::endl;
                    err = 0;
                    if (num==0){
                        err2 = 0;
                    }
                    break;
                }
                
                if (err==2){
                    std::cerr << "Error : negative number" << std::endl;
                    err = 0;
                    break;
                }
                
                // memory edges
                if (m%2){
                    std::cerr << "Error : invalid edges" << std::endl;
                    break;
                }
                
                for (int j=0;j<m;j=j+2){
                    
                    if (e[j] == e[j+1])
                        continue;
                    edgs[e[j]][e[j+1]] = 1;
                    edgs[e[j+1]][e[j]] = 1;
                }
                
                for (int i=0; i<=m;++i)
                {
                    v.push_back(e[i]);
                }
                err1 = 3;       //first successful E command
                err2 = 0;       //can input V again after successful E command
                
                sem_post(&io_se1);
                sem_post(&io_se2);
                sem_post(&io_se3);
                
                // output
                
                sem_wait(&vc_se1);
                sem_wait(&vc_se2);
                sem_wait(&vc_se3);
                
                //                std::cout << "minisat time : " << total_time1 << " ms"<< std::endl;
                //                std::cout << "vc1 time : " << total_time2 << " ms"<< std::endl;
                //                std::cout << "vc2 time : " << total_time3 << " ms"<< std::endl;
                
                ////////////// minisat output ///////////////
                
                if (vertex_cover1.size() > 0){
                    std::sort(vertex_cover1.begin(), vertex_cover1.end());
                    
                    std::cout << "CNF-SAT-VC: ";
                    for (unsigned int i=0;i<vertex_cover1.size()-1;i++){
                        std::cout << vertex_cover1[i] << ",";
                    }
                    std::cout << vertex_cover1[vertex_cover1.size()-1];
                    std::cout << std::endl;
                }
                else{
                    std::cout << "CNF-SAT-VC: " << std::endl;
                }
                
                vertex_cover1.clear();
                
                /////////////////////////////////////////////
                
                
                
                ///////////////// vc1 output ///////////////
                
                if (vertex_cover2.size() > 0){
                    std::sort(vertex_cover2.begin(), vertex_cover2.end());
                    
                    std::cout << "APPROX-VC-1: ";
                    for (unsigned int i=0;i<vertex_cover2.size()-1;i++){
                        std::cout << vertex_cover2[i] << ",";
                    }
                    std::cout << vertex_cover2[vertex_cover2.size()-1];
                    std::cout << std::endl;
                }
                else{
                    std::cout << "APPROX-VC-1: " << std::endl;
                }
                
                vertex_cover2.clear();
                
                /////////////////////////////////////////////
                
                
                
                ///////////////// vc2 output ///////////////
                
                if (vertex_cover3.size() > 0){
                    std::sort(vertex_cover3.begin(), vertex_cover3.end());
                    
                    std::cout << "APPROX-VC-2: ";
                    for (unsigned int i=0;i<vertex_cover3.size()-1;i++){
                        std::cout << vertex_cover3[i] << ",";
                    }
                    std::cout << vertex_cover3[vertex_cover3.size()-1];
                    std::cout << std::endl;
                }
                else{
                    std::cout << "APPROX-VC-2: " << std::endl;
                }
                
                vertex_cover3.clear();
                v.clear();
                
                /////////////////////////////////////////////
            }
                
                // other command
            default:
                
                break;
        }
        
    }
    
    io_exit = 1;
    sem_post(&io_se1);
    sem_post(&io_se2);
    sem_post(&io_se3);
    return 0;
}

void* minisat(void *arg){
    
    while(true){
        
        sem_wait(&io_se1);
        
        if (io_exit==1)
            break;
        
        //        clockid_t cid;
        //        int retcode;
        //        struct timespec start, stop;
        //
        //        retcode = pthread_getcpuclockid(pthread_self(), &cid);
        //
        //        if(retcode){
        //            std::cerr << "Error: could not get thread time(minisat) " << std::endl;
        //        }
        //        if (clock_gettime(cid, &start)==-1){
        //            std::cout << "Error : can not get thread minisat time" << std::endl;
        //        }
        
        /////////////////////////////////////////////////////
        
        int edge_num = 0;
        
        for (int i=0;i<num;i++){
            for (int j=i;j<num;j++){
                if (edgs[i][j] == 1){
                    edge_num++;
                }
            }
        }
        if (edge_num == 0){
            vertex_cover1.clear();
        }
        else{
            
            int vertex_min = 1;
            
            while (vertex_min<num)
            {
                
                // -- allocate on the heap so that we can reset later if needed
                std::unique_ptr<Minisat::Solver> solver(new Minisat::Solver());
                
                Minisat::Lit vertices[num][vertex_min];
                
                
                for (int i=0;i<num;i++){
                    for (int k=0;k<vertex_min;k++){
                        vertices[i][k] = Minisat::mkLit(solver->newVar());
                    }
                }
                
                //1.At least one vertex is the ith vertex in the vertex cover.
                
                for (int k=0;k<vertex_min;k++){
                    Minisat::vec<Minisat::Lit> clauses1;
                    for (int i=0;i<num;i++){
                        clauses1.push(vertices[i][k]);
                    }
                    solver->addClause (clauses1);
                }
                
                //2.No one vertex can appear twice in a vertex cover.
                
                for (int i=0;i<num;i++){
                    for (int m=0;m<vertex_min-1;m++){
                        for (int n=m+1;n<vertex_min;n++){
                            Minisat::vec<Minisat::Lit> clauses2;
                            clauses2.push(~vertices[i][m]);
                            clauses2.push(~vertices[i][n]);
                            solver->addClause (clauses2);
                        }
                    }
                }
                
                //3.No more than one vertex appears in the mth position of the vertex cover.
                
                for (int k=0;k<vertex_min;k++){
                    for (int m=0;m<num-1;m++){
                        for (int n=m+1;n<num;n++){
                            Minisat::vec<Minisat::Lit> clauses3;
                            clauses3.push(~vertices[m][k]);
                            clauses3.push(~vertices[n][k]);
                            solver->addClause (clauses3);
                        }
                    }
                }
                
                //4.Every edge is incident to at least one vertex in the vertex cover.
                
                for (int i=0;i<num;i++){
                    for (int j=i;j<num;j++){
                        if (edgs[i][j] == 1 && i!=j){
                            Minisat::vec<Minisat::Lit> clauses4;
                            for (int k=0;k<vertex_min;k++){
                                clauses4.push(vertices[i][k]);
                                clauses4.push(vertices[j][k]);
                            }
                            solver->addClause (clauses4);
                        }
                    }
                }
                
                bool res = false;
                res = solver->solve();
                
                if (res){
                    //std::cout << vertex_min << " is the number of vertex cover\n";
                    for (int k=0;k<vertex_min;k++){
                        for (int i=0;i<num;i++){
                            if (Minisat::toInt(solver->modelValue(vertices[i][k]))==0){
                                vertex_cover1.push_back(i);
                            }
                        }
                    }
                    
                    break;
                }
                
                solver.reset (new Minisat::Solver());
                
                vertex_min++;
            }
        }
        
        ////////////////////////////////////////////////////////////
        
        //        if (clock_gettime(cid, &stop)==-1){
        //        std::cout << "Error : can not get thread minisat time" << std::endl;
        //        }
        //        total_time1 = (stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_nsec - start.tv_nsec)/1000000;
        
        sem_post(&vc_se1);
    }
    
    return 0;
}

void* vc1(void *arg){
    
    
    while(true)
    {
        
        sem_wait(&io_se2);
        
        if (io_exit==1)
            break;
        
        //        clockid_t cid;
        //        int retcode;
        //        struct timespec start, stop;
        //
        //        retcode = pthread_getcpuclockid(pthread_self(), &cid);
        //
        //        if(retcode){
        //            std::cerr << "Error: could not get thread time(vc1) " << std::endl;
        //        }
        //        if (clock_gettime(cid, &start)==-1){
        //            std::cout << "Error : can not get thread vc1 time" << std::endl;
        //        }
        
        ////////////////////////////////////////////////////////
        int temp_arr[num][num];
        
        for (int i=0;i<num;i++){
            for (int j=0;j<num;j++){
                temp_arr[i][j] = edgs[i][j];
            }
        }
        
        while (true)
        {
            
            int degree_max = 0;
            int max_degree_vertex = 99;     //valid vertex number < 50
            
            for (int i=0;i<num;i++){
                int temp = 0;
                for (int j=0;j<num;j++){
                    if (temp_arr[i][j]==1 && i!=j){
                        temp++;
                    }
                }
                
                if (temp > degree_max){
                    degree_max = temp;
                    max_degree_vertex = i;
                }
            }
            
            if (degree_max == 0)
                break;  // all edges have been removed
            
            for (int i=0;i<num;i++){
                for (int j=0;j<num;j++){
                    if ((max_degree_vertex == i) || (max_degree_vertex == j)){
                        temp_arr[i][j] = maxv;
                    }
                }
            }
            
            vertex_cover2.push_back(max_degree_vertex);
            
        }
        
        /////////////////////////////////////////////////////////
        
        //        if (clock_gettime(cid, &stop)==-1){
        //        std::cout << "Error : can not get thread vc1 time" << std::endl;
        //        }
        //        total_time2 = (stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_nsec - start.tv_nsec)/1000000;
        
        sem_post(&vc_se2);
    }
    
    return 0;
}

void* vc2(void *arg){
    
    while(true)
    {
        
        sem_wait(&io_se3);
        
        if (io_exit==1)
            break;
        
        //        clockid_t cid;
        //        int retcode;
        //        struct timespec start, stop;
        //
        //        retcode = pthread_getcpuclockid(pthread_self(), &cid);
        //
        //        if(retcode){
        //            std::cerr << "Error: could not get thread time(vc2) " << std::endl;
        //        }
        //        if (clock_gettime(cid, &start)==-1){
        //            std::cout << "Error : can not get thread vc2 time" << std::endl;
        //        }
        
        /////////////////////////////////////////////////////
        
        int k1=0;
        int fir = 0;
        for (fir=0;fir<v.size()-1;fir=fir+2){
            if (v[fir]!=v[fir+1]){
                k1=1;
                break;
            }
        }
        if (k1==1){
            
            vertex_cover3.push_back(v[fir]);
            vertex_cover3.push_back(v[fir+1]);
            
            for (unsigned int i = fir+2; i < v.size()-1; i=i+2)
            {
                int k;
                for(unsigned int j=0; j< vertex_cover3.size(); ++j)
                {
                    if((v[i]!=vertex_cover3[j]) && (v[i+1]!=vertex_cover3[j]) && v[i]!=v[i+1])
                    {
                        k=1;
                    }
                    else
                    {
                        k=0;
                        break;
                    }
                }
                if(k ==1)
                {
                    vertex_cover3.push_back(v[i]);
                    vertex_cover3.push_back(v[i+1]);
                }
            }
        }
        
        //////////////////////////////////////////////////////
        
        //        if (clock_gettime(cid, &stop)==-1){
        //        std::cout << "Error : can not get thread vc2 time" << std::endl;
        //        }
        //        total_time3 = (stop.tv_sec - start.tv_sec) * 1000 + (double)(stop.tv_nsec - start.tv_nsec)/1000000;
        
        sem_post(&vc_se3);
    }
    return 0;
}


int main(int argc, char** argv) {
    
    // initialize graph matrix
    for (int i=0;i<max_vertex_num;i++){
        for (int j=0;j<max_vertex_num;j++){
            edgs[i][j] = maxv;
        }
    }
    
    pthread_t thread1, thread2, thread3, thread4;
    int  iret1, iret2, iret3 ,iret4;
    
    sem_init(&vc_se1, 0, 0);
    sem_init(&vc_se2, 0, 0);
    sem_init(&vc_se3, 0, 0);
    sem_init(&io_se1, 0, 0);
    sem_init(&io_se1, 0, 0);
    sem_init(&io_se1, 0, 0);
    
    /* Create independent threads each of which will execute function */
    iret1 = pthread_create( &thread1, NULL, io, NULL);
    if(iret1)
    {
        fprintf(stderr,"Error - pthread_create() return code: %d\n",iret1);
        return -1;
    }
    
    iret2 = pthread_create( &thread2, NULL, minisat, NULL);
    if(iret2)
    {
        fprintf(stderr,"Error - pthread_create() return code: %d\n",iret2);
        return -1;
    }
    
    iret3 = pthread_create( &thread3, NULL, vc1, NULL);
    if(iret3)
    {
        fprintf(stderr,"Error - pthread_create() return code: %d\n",iret3);
        return -1;
    }
    
    iret4 = pthread_create( &thread4, NULL, vc2, NULL);
    if(iret4)
    {
        fprintf(stderr,"Error - pthread_create() return code: %d\n",iret4);
        return -1;
    }
    
    pthread_join( thread1, NULL);
    pthread_join( thread2, NULL);
    pthread_join( thread3, NULL);
    pthread_join( thread4, NULL);
    
    return 0;
}
