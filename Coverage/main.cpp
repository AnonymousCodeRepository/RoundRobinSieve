#include "time.h"
#include "utilityfunction.h"
#include "ParSSP.h"
#include "RepSampling.h"
#include "RoundRobinSieve.h"
#include "RandomMultiGreedy.h"
int main(int argc,char *argv[]) {
    read_data();

    time_t nowtime;
    struct tm *p;;
    time(&nowtime);
    p = localtime(&nowtime);
    string::size_type pos1, pos2, posend;
    pos1 = edge_text.find_last_of("/");
    pos2 = edge_text.rfind("/", pos1 - 1);
    posend = edge_text.find_last_not_of("/");
    string name1 = edge_text.substr(pos2 + 1, pos1 - pos2 - 1);
    string name2 = edge_text.substr(pos1 + 1, posend);
    string result_name = name1 + "_" + name2;
    //cout<<result_name<<endl;
    string outtext =
            "./result/result_" + result_name + "_" + to_string(p->tm_mon + 1) + "." + to_string(p->tm_mday) + "_" +
            to_string(p->tm_hour) + "_" + to_string(p->tm_min) + "_" + to_string(p->tm_sec) + ".txt";

//    S_class test;
//    unordered_set<int> t={476, 539, 378 ,477, 330, 663 ,768 ,649 ,686, 797 ,70 ,970, 600 ,535, 921,1000 ,454, 565 ,616 ,403,956, 652, 958, 800, 841, 957, 878, 912, 43, 274,706 ,743, 369, 1001, 696 ,643 ,954 ,534 ,978, 947 ,789};
//    test.solution=t;
//    cout<<test.f_S()<<endl;
//
//    return 0;



    vector<Result> distorted_result;
    vector<Result> parssp_result;
    vector<Result> repsampling_result;
    vector<Result> rrs_result;
    vector<Result> rmg_result;



    vector<int> ground_set;
    for(int i=0;i<node_num;i++)
        ground_set.push_back(i);

    //double eps = atof(argv[2]);
    //cout << "eps: " << eps << endl;
    double eps=0.1;
    int m_start=100;
    int m_end=1000;
    int m_step=100;

    //int m=999999999;

    const int num_runs = 1; // 定义随机算法的重复运行次数

    for(int m=m_start; m<=m_end; m+=m_step)
    {
        // cout << "Running experiments for m = " << m << "..." << endl;

        rrs_result.push_back(RoundRobinSieve(eps, m));

        double total_parssp_revenue = 0.0;
        long long total_parssp_oracle = 0;
        int total_parssp_round = 0;
        for (int i = 0; i < num_runs; ++i) {
            Result res = ParSSP(eps, m);
            total_parssp_revenue += res.revenue;
            total_parssp_oracle += res.oracle;
            total_parssp_round += res.round;
        }
        Result avg_parssp_result = {
            total_parssp_revenue / num_runs,
            total_parssp_oracle / num_runs,
            total_parssp_round / num_runs
        };
        parssp_result.push_back(avg_parssp_result);

        double total_repsampling_revenue = 0.0;
        long long total_repsampling_oracle = 0;
        int total_repsampling_round = 0;
        for (int i = 0; i < num_runs; ++i) {
            Result res = RepSampling(eps, m);
            total_repsampling_revenue += res.revenue;
            total_repsampling_oracle += res.oracle;
            total_repsampling_round += res.round;
        }
        Result avg_repsampling_result = {
            total_repsampling_revenue / num_runs,
            total_repsampling_oracle / num_runs,
            total_repsampling_round / num_runs
        };
        repsampling_result.push_back(avg_repsampling_result);

        double total_rmg_revenue = 0.0;
        long long total_rmg_oracle = 0;
        int total_rmg_round = 0;
        for (int i = 0; i < num_runs; ++i) {
            Result res = RandomMultiGreedy(eps, m);
            total_rmg_revenue += res.revenue;
            total_rmg_oracle += res.oracle;
            total_rmg_round += res.round;
        }
        Result avg_rmg_result = {
            total_rmg_revenue / num_runs,
            total_rmg_oracle / num_runs,
            total_rmg_round / num_runs
        };
        rmg_result.push_back(avg_rmg_result);

    }

    ofstream out(outtext);
    out<<"group limit: "<<group_limit<<endl;
    out<<"eps: "<<eps<<endl;
    out<<"max node: "<<endl;
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


    out<<"ParSSP "<<endl;
    out<<"objective function value: "<<endl;
    for(auto p:parssp_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"query: "<<endl;
    for(auto p:parssp_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;
    out<<"adaptive round: "<<endl;
    for(auto p:parssp_result)
    {
        out<<p.round<<"\t";
    }
    out<<endl;


    out<<"RepSampling "<<endl;
    out<<"objective function value: "<<endl;
    for(auto p:repsampling_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"query: "<<endl;
    for(auto p:repsampling_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;
    out<<"adaptive round: "<<endl;
    for(auto p:repsampling_result)
    {
        out<<p.round<<"\t";
    }
    out<<endl;



    out<<"RandomMultiGreedy "<<endl;
    out<<"objective function value: "<<endl;
    for(auto p:rmg_result)
    {
        out<<p.revenue<<"\t";
    }
    out<<endl;
    out<<"query: "<<endl;
    for(auto p:rmg_result)
    {
        out<<p.oracle<<"\t";
    }
    out<<endl;
    out<<"adaptive round: "<<endl;
    for(auto p:rmg_result)
    {
        out<<p.round<<"\t";
    }
    out<<endl;


    return 0;
}