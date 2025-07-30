//
// Created by HP on 25-7-15.
//

#ifndef ALTERNATEFILTER_H
#define ALTERNATEFILTER_H
#include "utilityfunction.h"
double p_af=2.0/(1.0+sqrt(k));
// double p_af=2.0/3.0;
// double p_af=0.7;
default_random_engine e_af(1);
bernoulli_distribution u_af(p_af);

void ProbSampleSieve(const vector<int> &N,vector<int> &V_selected,S_class &T,S_class &T_,double tau,double eps,int max_node,long long int &oracle_times, long long int &rounds) {
    double alpha=1.0+1.0/(k+sqrt(k)-1.0);
    // double alpha=1.0/(1.0+k);
    // cout<<"alpha is setting and p is 1"<<endl;

    double log_val = log(eps/static_cast<double>(max_node)) / log(1.0-eps);
    double delta = 1.0 / ceil(log_val);

    double M = (4 * log(max_node) / eps + 4 * log(max_node / delta) + 4 * sqrt((log(max_node) * log(max_node / delta)) / eps)) / p_af;
    // cout<<"M: "<<M<<endl;

    rounds++;
    vector<int> R;
    for(auto &e:N)
    {
        if(V_selected[e]==0) {
            if(T.is_feasible(e,max_node)) {
                oracle_times++;
                if (T.marginal_gain_scaled(e,alpha)>= tau)
                    R.push_back(e);
            }
        }
    }

    vector<int> last_W;
    bool added_to_T = false;

    for(int t=1;t<=M;t++) {

        if(R.empty())
            break;

        vector<int> W = R;
        int i = 1;
        vector<vector<int>> U;
        U.emplace_back();
        vector<vector<int>> W_sequence;
        W_sequence.emplace_back();
        W_sequence.push_back(W);

        S_class TcupUi(T);
        vector<bool> Ui_selected(node_num, false);

        while(!W.empty()) {

            int idx = rand() % W.size();
            int ui = W[idx];


            vector<int> Ui;
            if(i>1) {
                Ui = U.back();
            }
            Ui.push_back(ui);
            U.push_back(Ui);
            TcupUi.add_element(0.0,ui);
            Ui_selected[ui] = true;


            vector<int> W_next;
            for(auto &e : R) {
                if(!Ui_selected[e] && TcupUi.is_feasible(e, max_node)) {
                    W_next.push_back(e);
                }
            }
            W = W_next;
            W_sequence.push_back(W);
            i++;
        }

        int z = i - 1;
        vector<string> Status(z+1, "Null");


        rounds++;
        for(int i = 1; i <= z; i++) {
            int ui = U[i].back();

            S_class T_cup_Ui_1(T);
            T_cup_Ui_1.add_batched_element(U[i-1].begin(),U[i-1].end());
            double scaled_marginal_gain = T_cup_Ui_1.marginal_gain_scaled(ui,alpha);
            if(scaled_marginal_gain >= tau) {
                Status[i] = "Above_τ";
            }
            if(scaled_marginal_gain <= 0) {
                Status[i] = "NonPositive";
            }
        }


        int i_star = 0;
        for(int i = z; i >= 1; i--) {
            int count = 0;
            for(int j = 1; j <= i; j++) {
                if(Status[j] == "Above_τ") {
                    count++;
                }
            }

            if(count >= (1.0-eps)*i) {
                i_star = i;
                break;
            }
        }


        added_to_T = false;
        if(i_star == 0) {cout<<"i_star wrong！"<<endl;}

        if(i_star > 0) {
            last_W = W_sequence[i_star+1];

            for(auto &u : U[i_star]) {
                V_selected[u] = 1;
            }


            if(u_af(e_af)) {
                added_to_T = true;
                T.add_batched_element(U[i_star].begin(),U[i_star].end());

                for(int i = 1; i <= i_star; i++) {
                    if(Status[i] != "NonPositive") {
                        T_.add_element(U[i].back());
                    }
                }

            }
        }


        rounds++;

        if(added_to_T) {
            R.clear();
            for(auto &e : last_W) {
                if(V_selected[e] == 0) {
                    oracle_times++;
                    if(T.marginal_gain_scaled(e,alpha) >= tau) {
                        R.push_back(e);
                    }
                }
            }
        } else {

            vector<int> R_temp(R);
            R.clear();
            for(auto &e:R_temp) {
                if(V_selected[e]==0) {//only excluded the element that in the V
                    // if(T.is_feasible(e,max_node)) {
                    //     oracle_times++;
                    //     if (T.marginal_gain_scaled(e,alpha)>= tau)
                    //         R.push_back(e);
                    // }

                    R.push_back(e);
                }
            }

        }
    }

}

Result RoundRobinSieve(double eps, int max_node)
{
    int solution_num = 2;
    cout<<"RoundRobinSieve & Max movie: "<<max_node<<endl;
    long long int oracle_times=0;
    long long int rounds=0;
    int u_star=-1;
    double u_value=0.0;
    rounds++;
    vector<int> N;
    for(int e=0;e<node_num;e++) {
        double temp_value = f_u(e);
        oracle_times++;
        if(temp_value<=0)
            continue;
        N.push_back(e);
        if (temp_value > u_value) {
            u_star = e;
            u_value = temp_value;
        }
    }
    double tau_max=u_value;

    double ell=ceil(log(eps/(double)max_node)/log(1.0-eps));
    //    double M=1.0/(pow(eps,3))*log((double)node_num/eps);
    //    double ell=log(eps/(double)node_num)/log(1.0-eps);

    vector<int> V_selected(node_num,0);
    vector<S_class> S(solution_num);
    vector<S_class> S_(solution_num);
    for(int i=0;i<ell;i++)
    {
         int j=i%solution_num;
        //int j=0;
        // cout<<"test always j=0"<<endl;

        double tau_i=tau_max*pow(1.0-eps, i);
        // cout<<"i "<<i<<" tau: "<<tau_i<<endl;
        //cout<<S.Set.size()<<endl;
        ProbSampleSieve(N,V_selected,S[j],S_[j],tau_i,eps,max_node,oracle_times,rounds);
    }

    S_class best_S = S[0];

    if (S[1].S_revenue > best_S.S_revenue) {
        best_S = S[1];
    }
    if (S_[0].S_revenue > best_S.S_revenue) {
        best_S = S_[0];
    }
    if (S_[1].S_revenue > best_S.S_revenue) {
        best_S = S_[1];
    }

    if (u_value > best_S.S_revenue) {

        best_S.replace_with_singleton(u_value, u_star);
    }

    Result final(best_S.S_revenue,best_S.S_cost,best_S.solution.size(),oracle_times,rounds);
    cout<<"objective: "<<best_S.S_revenue<<endl;
    cout<<"query: "<<oracle_times<<endl;
    cout<<"round: "<<rounds<<endl;
    return final;
}

#endif //ALTERNATEFILTER_H
