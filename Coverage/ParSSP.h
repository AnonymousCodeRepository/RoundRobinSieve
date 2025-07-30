
#ifndef MOV_PARSSP_H
#define MOV_PARSSP_H
#include "random"
#include "algorithm"
#include "utilityfunction.h"
double p_2=1.0/(1.0+sqrt(k+1.0));
default_random_engine eng(1);
bernoulli_distribution u_2(p_2);
// random_device rd;
mt19937 g(1);
class Gi
{
public:
    Gi()
    {
        set.clear();
        value=0.0;
    }
    vector<int> set;
    double value;
};
vector<int>SampleSeq(S_class S,const vector<int> &V,int max_node)
{
    vector<int> X;
    for(auto &u:V) {
        if (S.is_feasible(u,max_node))
            X.push_back(u);
    }
    vector<int> A;
    while (!X.empty())
    {
        shuffle(X.begin(), X.end(), g);//random permutation
        for(auto &e:X)
        {
            if(S.is_feasible(e,max_node))
            {
                S.add_element(0.0,e);//revenue is useless, so let marginal gain zero
                A.push_back(e);
            }
            else
                break;
        }
        X.erase(std::remove_if(X.begin(), X.end(), [&S,max_node](int e){
            return !S.is_feasible(e,max_node);
        }), X.end());
    }
    return A;
}
void BatchedSeq2(const vector<int> &N,vector<int> &U_selected,vector<int> &L_selected,S_class &S,double tau,double eps,int M,int max_node,long long int &oracle_times, long long int &rounds)
{
    int ctr=0; int FLAG=0;
    //initialize X
    rounds++;
    vector<int> X;
    for(auto &e:N)
    {
        if(U_selected[e]==0&&L_selected[e]==0) {
            oracle_times++;
            if ((S.is_feasible(e,max_node))&&(S.marginal_gain(e)>= tau))
                X.push_back(e);
        }
    }
    //cout<<X.size()<<endl;
    vector<int> A_selected(node_num,0);
    while(!X.empty())
    {
        vector<int> A_seq=SampleSeq(S,X,max_node);

        vector<Gi> G;//Add one more set to make G[i]=G_i
        for(int iter=0;iter<=A_seq.size();iter++)
            G.push_back(Gi());

        int bl=1;int br=A_seq.size();
        while (bl<br)
        {
            rounds++;

            int i=floor((bl+br)/2.0);
            //For simplicity, it can be optimized
            //construct A_i
            fill(A_selected.begin(), A_selected.end(), 0);
            vector<int> Ai;
            for (int iter = 0; iter < i; ++iter) {
                Ai.push_back(A_seq[iter]);
                A_selected[A_seq[iter]]=1;
            }
            //construct S\cup A_i
            S_class S_cup_Ai(S);

            S_cup_Ai.add_batched_element(Ai.begin(),Ai.end());
            //construct Gi, E_i and its value
            vector<int> Ei; double Ei_value=0.0;
            for (auto &a:X) {
                //check whether a in Ai
                if(A_selected[a]==1) continue;
                oracle_times++;//the number of queries used to calculate E_i is the most
                double marginal=S_cup_Ai.marginal_gain(a);

                //check f(a|S\cup A_i)\geq c(a) and a is feasible for S\cup A_i
                if(S_cup_Ai.is_feasible(a,max_node)&&marginal>=tau) {
                    G[i].set.push_back(a);
                    G[i].value += marginal;
                }
                //check f(a|S\cup A_i)< \tau*c(a) for Ei
                if(marginal<0.0) {
                    Ei.push_back(a);
                    Ei_value += marginal;
                }
            }
            //construct K_i and its value
            vector<int> Ki; double Ki_value=0.0;
            //initialize S\cup A_s
            S_class S_cup_As(S);
            for(auto &a:Ai)
            {
                oracle_times++;
                double marginal=S_cup_As.marginal_gain(a);
                if(marginal<0.0) {
                    Ki.push_back(a); Ki_value += marginal;
                }
                S_cup_As.add_element(0.0,a);//revenue is useless
            }
            bool c1=false; bool c2=false;
            //judge c1
            //////knapsack version

/*            //calculate c(X)
            double X_cost=0.0;
            for(auto &e:X)
                X_cost+=1.0;
            //calculate c(Gi)
            double Gi_cost=0.0;
            for(auto &e:G[i].set)
                Gi_cost+=1.0;
            if(Gi_cost<=(1.0-eps)*X_cost) c1=true;*/

            //////cardinality/ksystem version
            if((double)G[i].set.size()<=(1.0-eps)*(double)X.size())
                c1=true;

            //judge c2
            if(eps*G[i].value<=fabs(Ei_value)+fabs(Ki_value)) c2=true;

            if(c1||c2) br=i;
            else bl=i+1;
            if(c2&&(!c1)) FLAG=1;
            else FLAG=0;
        }
        int k=br;
        //cout<<"k "<<k<<endl;
        //for(int iter=0;iter<k;iter++)
        //    S.add_element(A_seq[iter]);
        for(int iter=0;iter<k;iter++)//U\leftarrow U\cup A_k^*
            U_selected[A_seq[iter]]=1;

        if(u_2(eng))
        {
            S.add_batched_element(A_seq.begin(),A_seq.begin()+k);
            X.assign(G[k].set.begin(), G[k].set.end());//copy G_k to X;
            ctr+=FLAG;
        }
        else {
            X.erase(std::remove_if(X.begin(), X.end(),
                                   [&](auto e){ return U_selected[e] == 1; }),
                    X.end());
/*            for (auto e = X.begin(); e != X.end();) {
                if (U_selected[(*e)] == 1) {
                    e = X.erase(e);
                } else {
                    e++;
                }
            }*/
        }
        if(ctr>=M) {
            for (auto &u:G[k].set) {
                L_selected[u] = 1;
            }
            break;
        }
    }
//    cout<<"S*:"<<endl;
//    cout<<"  revenue: "<<S.S_revenue<<" cost: "<<S.S_cost<<" size: "<<S.Set.size()<<endl;
//    for(int i=0;i<S.Set.size();i++)
//        cout<<S.Set[i]<<" ";
//    cout<<endl;

}
Result ParSSP(double eps, int max_node)
{
    cout<<"ParSSP & Max node: "<<max_node<<endl;
    long long int oracle_times=0;
    long long int rounds=0;
    int u_star=-1;
    double u_value=0.0;
    rounds++;
    vector<int> N;
    for(int e=0;e<node_num;e++) {
        N.push_back(e);

        double temp_value = f_u(e);
        oracle_times++;

        if (temp_value > u_value) {
            u_star = e;
            u_value = temp_value;
        }
    }
    double rho_max=u_value;

    double M=ceil((log(eps/(double)max_node)/log(1.0-eps)+2.0)/(pow(eps,2)));
    double ell=ceil(log(eps/(double)max_node)/log(1.0-eps))+1.0;
//    double M=1.0/(pow(eps,3))*log((double)node_num/eps);
//    double ell=log(eps/(double)node_num)/log(1.0-eps);


    vector<int> U_selected(node_num,0);
    vector<int> L_selected(node_num,0);
    S_class S;
    for(int i=0;i<ell;i++)
    {
        double rho_i=rho_max*pow(1.0-eps, i);
        // cout<<"i "<<i<<" tau: "<<rho_i<<endl;
        //cout<<S.Set.size()<<endl;
        BatchedSeq2(N,U_selected,L_selected,S,rho_i,eps,(int)M,max_node,oracle_times,rounds);
    }
    Result final(S.S_revenue,S.S_cost,S.solution.size(),oracle_times,rounds);
    if(u_value>S.S_revenue)//compare u2_star and S_
    {
        final.revenue=u_value;
        final.size=1;
        S.replace_with_singleton(u_value,u_star);
    }

    cout<<"objective value: "<<S.S_revenue<<endl;
    cout<<"query: "<<oracle_times<<endl;
    cout<<"round: "<<rounds<<endl;
    return final;
}
#endif //MOV_PARSSP_H
