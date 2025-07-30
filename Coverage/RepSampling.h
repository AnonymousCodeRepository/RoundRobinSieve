//
// Created by CLAKERS on 2023/3/16.
//

#ifndef MOV_REPSAMPLING_H
#define MOV_REPSAMPLING_H
#include "utilityfunction.h"
#include "ParSSP.h"
S_class RandSampling(double eps,double lambda_rs, const vector<int> &V, long long int &rounds, long long int &oracle_times,int max_node) {
    double delta = -999999.0;
    rounds++;
    vector<int> X;
    for (auto e:V) {
        oracle_times++;
        double fu=f_u(e);

        if (fu > delta)
            delta = fu;
        X.push_back(e);
    }
    double delta_0 = lambda_rs * delta;

    int ctr = 0;
    int FLAG = 0;

    vector<int> A_selected(node_num, 0);
    S_class S;
    while (delta >= delta_0)
    {
        while (!X.empty())
        {
            vector<int> A_seq = SampleSeq(S, X, max_node);

            vector<Gi> G;//Add one more set to make G[i]=G_i
            for (int iter = 0; iter <= A_seq.size(); iter++)
                G.push_back(Gi());

            int bl = 1;
            int br = A_seq.size();
            while (bl < br) {
                rounds++;

                int i = floor((bl + br) / 2.0);
                //For simplicity, it can be optimized
                //construct A_i
                fill(A_selected.begin(), A_selected.end(), 0);
                vector<int> Ai;
                for (int iter = 0; iter < i; ++iter) {
                    Ai.push_back(A_seq[iter]);
                    A_selected[A_seq[iter]] = 1;
                }
                //construct S\cup A_i
                S_class S_cup_Ai(S);

                S_cup_Ai.add_batched_element(Ai.begin(), Ai.end());
                //construct Gi, E_i and its value
/*                vector<int> Ei;
                double Ei_value = 0.0;*/
                for (auto &a:X) {
                    //check whether a in Ai
                    if (A_selected[a] == 1) continue;
                    oracle_times++;//the number of queries used to calculate E_i is the most
                    double marginal = S_cup_Ai.marginal_gain(a);

                    //check f(a|S\cup A_i)\geq c(a) and a is feasible for S\cup A_i
                    if (S_cup_Ai.is_feasible(a, max_node) && marginal >= delta) {
                        G[i].set.push_back(a);
                        G[i].value += marginal;
                    }
                    //check f(a|S\cup A_i)< \tau*c(a) for Ei
/*                    if (marginal < 0.0) {
                        Ei.push_back(a);
                        Ei_value += marginal;
                    }*/
                }
                //construct K_i and its value
/*                vector<int> Ki;
                double Ki_value = 0.0;
                //initialize S\cup A_s
                S_gamma S_cup_As;
                S_cup_As.copy(S);
                for (auto &a:Ai) {
                    oracle_times++;
                    double marginal = S_cup_As.marginal(a);
                    if (marginal < 0.0) {
                        Ki.push_back(a);
                        Ki_value += marginal;
                    }
                    S_cup_As.add_element_known(0.0, a);//revenue is useless
                }*/
                bool c1 = false;
                //judge c1
                //////cardinality/ksystem version
                if ((double) G[i].set.size() <= (1.0 - eps) * (double) X.size())
                    c1 = true;

                if (c1) br = i;
                else bl = i + 1;
            }
            int k = br;

            S.add_batched_element(A_seq.begin(), A_seq.begin() + k);
            X.assign(G[k].set.begin(), G[k].set.end());//copy G_k to X;
        }
        delta=delta*(1.0-eps);
        X.clear();
        for(auto e:V)
        {
            if(S.is_feasible(e,max_node)&&S.marginal_gain(e)>=delta)
                X.push_back(e);
        }
    }
    return S;
}
default_random_engine engine2(1);
S_class USM(const unordered_set<int> &A)
{
    S_class T;
    bernoulli_distribution u(0.5);
    for(auto &a:A) {
        if (u(engine2))
            T.add_element(a);
    }
    return T;
}
Result RepSampling(double eps, int max_node)
{
    cout<<"RepSampling & Max node: "<<max_node<<endl;
    double m=1.0+ceil(sqrt((k+1.0)/2.0));
    double lambda_rs=eps*(k+1.0)/m;

    long long int oracle_times=0;
    long long int rounds=0;

    vector<int> node_available(node_num,1);
    vector<int> N;
    for(int i=0;i<node_num;i++)//N has all nodes at first time
        N.push_back(i);

    vector<S_class> U;
    for(int i=0;i<m;i++)
    {
        S_class Si=RandSampling(eps,lambda_rs,N,rounds,oracle_times,max_node);
        S_class Si_=USM(Si.solution);
        U.push_back(Si);
        U.push_back(Si_);

        for(const auto &j:Si.solution) {
            node_available[j]=0;
        }
        N.clear();
        for(int iter=0;iter<node_num;iter++)//rebuild groundset N
        {
            if(node_available[iter]==1)
                N.push_back(iter);
        }
    }
    double max_value=-999999999;
    S_class *max_S;
    for(int i=0;i<U.size();i++)
    {
        if(U[i].S_revenue>max_value)
        {
            max_S = &U[i];
            max_value=U[i].S_revenue;
        }
    }

    Result final((*max_S).S_revenue,(*max_S).S_cost,(*max_S).solution.size(),oracle_times,rounds);

    cout<<"objective value: "<<(*max_S).S_revenue<<endl;
    cout<<"query: "<<oracle_times<<endl;
    cout<<"round: "<<rounds<<endl;

    return final;
}

#endif //MOV_REPSAMPLING_H
