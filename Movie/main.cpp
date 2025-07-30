#include "utilityfunction.h"
#include <stdlib.h>
#include "time.h"
#include "DistortedFiltering.h"
#include "RoundRobinSieve.h"
int main(int argc,char *argv[]) {

    //double B=atof(argv[1]);
    //double eps=atof(argv[2]);
    double eps=0.1;
    // cout<<"eps: "<<eps<<endl;

    string filename = "movies.txt";
    dataset *da = new dataset(filename);

    time_t nowtime;
    struct tm* p;;
    time(&nowtime);
    p = localtime(&nowtime);
    string outtext="./result/movie_result_"+to_string((int)ave_num)+"_"+to_string(p->tm_mon+1)+"."+to_string(p->tm_mday)+"_"+to_string(p->tm_hour)+"_"+to_string(p->tm_min)+"_"+to_string(p->tm_sec)+".txt";

    lambda_f=1.0;


    vector<Result> distorted_result;
    vector<Result> rrs_result;

    //vector<Result> icml_result;

    vector<int> ground_set;
    for(int i=0;i<node_num;i++)
        ground_set.push_back(i);

    int m_start=1;
    int m_end=10;
    int m_step=1;

    const int num_runs = 1;
    //for(genres_limit=m_start;genres_limit<=m_end;genres_limit+=m_step)
    for(max_movie = m_start; max_movie <= m_end; max_movie += m_step) {

        double total_rrs_revenue = 0.0;
        long long total_rrs_oracle = 0;
        int total_rrs_round = 0;

        for (int i = 0; i < num_runs; ++i) {
            Result single_run_result = RoundRobinSieve(eps, max_movie);
            total_rrs_revenue += single_run_result.revenue;
            total_rrs_oracle += single_run_result.oracle;
            total_rrs_round += single_run_result.round;
        }

        Result avg_rrs_result;
        avg_rrs_result.revenue = total_rrs_revenue / num_runs;
        avg_rrs_result.oracle = total_rrs_oracle / num_runs;
        avg_rrs_result.round = total_rrs_round / num_runs;

        rrs_result.push_back(avg_rrs_result);


        double total_distorted_revenue = 0.0;
        long long total_distorted_oracle = 0;
        int total_distorted_round = 0;

        for (int i = 0; i < num_runs; ++i) {
            Result single_run_result = DistortedFiltering(eps, max_movie);
            total_distorted_revenue += single_run_result.revenue;
            total_distorted_oracle += single_run_result.oracle;
            total_distorted_round += single_run_result.round;
        }

        Result avg_distorted_result;
        avg_distorted_result.revenue = total_distorted_revenue / num_runs;
        avg_distorted_result.oracle = total_distorted_oracle / num_runs;
        avg_distorted_result.round = total_distorted_round / num_runs;

        distorted_result.push_back(avg_distorted_result);
    }

    ofstream out(outtext);
    out<<"eps: "<<eps<<" lambda_f:"<<lambda_f<<" lambda:"<<da->lambda<<" genres limit:"<<genres_limit<<endl;
    for(int m=m_start;m<=m_end;m+=m_step)
    {
        out<<m<<"\t";
    }
    out<<endl;


    out<<"RoundRobinSieve "<<endl;
    out<<"objective function value: "<<endl;
    for(auto p:rrs_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"query: "<<endl;
    for(auto p:rrs_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;
    out<<"adaptive round: "<<endl;
    for(auto p:rrs_result)
    {
        out<<p.round<<"\t";
    }
    out<<endl;

    out<<"DistortedFiltering "<<endl;
    out<<"objective function value: "<<endl;
    for(auto p:distorted_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"query: "<<endl;
    for(auto p:distorted_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;
    out<<"adaptive round: "<<endl;
    for(auto p:distorted_result)
    {
        out<<p.round<<"\t";
    }
    out<<endl;


    return 0;
}
