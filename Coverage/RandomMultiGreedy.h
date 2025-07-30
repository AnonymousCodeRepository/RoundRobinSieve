#ifndef COVERAGE_MULTIRANDOM_PLUS_H
#define COVERAGE_MULTIRANDOM_PLUS_H

#include "utilityfunction.h"
#include "time.h"
#include "algorithm"
#include "map"

default_random_engine e_for_rmg(12345);
Result RandomMultiGreedy(double eps,int max_node)
{
    cout<<"RandomMultiGreedy & Max node: "<<max_node<<endl;
    long long int oracle_times=0;
    long long int rounds=0;
    vector<S_class> S;

    double p=2.0/(1.0+sqrt((float)k));
    int ell=2;
    bernoulli_distribution u(p);

    for(int i=0;i<ell;i++)
        S.push_back(S_class());

    double delta_f=0.0;
    rounds++;
    for(int iter=0;iter<node_num;iter++)
    {
        oracle_times++;
        double value=f_u(iter);
        if (value > delta_f) {
            delta_f = value;
        }
    }
    vector<int> available(node_num,1);
    for(double tau=delta_f;tau>((eps/(double)max_node)*delta_f);tau=tau*(1.0-eps))
    {
        for(int element=0;element<node_num;element++)
            //for(int element=node_num-1;element>=0;element--)
        {
            if(available[element]==1)//if e\notin S_i,ie., u\in N_i-1
            {
                double max_marginal=-999.0;
                int max_solution=-1;

                bool round_trigger=false;
                for(int j=0;j<S.size();j++)
                {
                    //judge S+{e}\in I or not
                    bool flag=S[j].is_feasible(element,max_node);
                    //if S+{e}\in I
                    if(flag)
                    {
                        round_trigger=true;
                        double marginal=0.0;
                        oracle_times++;
                        marginal = S[j].marginal_gain(element);

                        if(marginal>max_marginal)
                        {
                            max_marginal=marginal;
                            max_solution=j;
                        }
                    }
                }
                if(round_trigger)
                    rounds++;

                if(max_solution==-1)//non element
                    continue;
                if(max_marginal>=tau)//if f(e|S_i)>=tau
                {
                    if(u(e_for_rmg)==1)
                    {
                        S[max_solution].add_element(max_marginal,element);
                    }
                    available[element] = 0;
                }
            }
        }
    }


    S_class *S_star;
    double max_revenue=-999999999.0;

    //    cout<<"MultiRandomAcc & m: "<<max_movie<<endl;
    for(int i=0;i<S.size();i++)
    {
        if(S[i].S_revenue>max_revenue)
        {
            S_star=&S[i];
            max_revenue=S[i].S_revenue;
        }
    }


    cout<<"objective value: "<<(*S_star).S_revenue<<endl;
    cout<<"query: "<<oracle_times<<endl;
    cout<<"round: "<<rounds<<endl;
    return Result((*S_star).S_revenue,(*S_star).S_cost,(*S_star).solution.size(),oracle_times,rounds);
}
#endif //COVERAGE_MULTIRANDOM_PLUS_H
