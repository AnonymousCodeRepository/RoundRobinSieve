
#ifndef STREAMING_ALGORITHM_UTILITY_FUNCTIONS_H
#define STREAMING_ALGORITHM_UTILITY_FUNCTIONS_H
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <sys/types.h>
#include <sstream>
#include <climits>
#include <string>
#include <vector>
#include <cmath>
#include <set>
#include <unordered_set>

#include "list"
#include "random"
//#define pi 3.1415926
using namespace std;
//Adventure Fantasy Animation
vector<float> sim_mat_sum;
vector<vector<float>> sim_mat;
vector<vector<int>> genres_index;//Adventure Fantasy Animation (1,1,1)
vector<float> modular_cost;

int max_movie=999999;
int genres_limit=5;
int node_num=-1;
const int d=1;

const double B=999999999.0;
double min_cost=999999999;
double max_cost=-1;
int K;
const int k=3;
int k_for_haba=3;
double lambda_f;
const double ave_num=1.0;

class dataset{
public:
    vector<int> index;
    set<int> v;
    vector<vector<float>> matrix;
    vector<vector<string>> genres;
    vector<string> titles;
    float lambda;
    vector<int> year;
    vector<float> rating;
    vector<float> year_costs;
    vector<float> rate_costs;
    vector<vector<double>> costs;
    //vector<vector<float>> sim_mat;
    dataset(string);
    void cost_normalize();
    void readfile(string filename);
    vector<float> feature_extractor(string features,char delim);
    vector<string> genre_extractor(string genres,char delim);
    void sim_mat_generator();
};
dataset::dataset(string filename)
{
    readfile(filename);
    cost_normalize();
    sim_mat_generator();

    node_num=index.size();
    sim_mat_sum = vector<float>(node_num,0.0);
    //for(int it=0;it < node_num;it++)
        //genres_index.push_back(vector<int>(k,0));
    for(int it=0;it < node_num;it++)
    {
        genres_index.push_back(vector<int>(k,0));
        for(auto p: genres[it])
        {
            if(p=="Adventure")
                genres_index[it][0]=1;
            if(p=="Fantasy")
                genres_index[it][1]=1;
            if(p=="Animation")
                genres_index[it][2]=1;
        }
        /*
        cout<<it<<"\t";
        for(auto p:genres_index[it])
            cout<<p<<"\t";
        cout<<endl;
*/
        for (int itv = 0; itv < node_num; itv++) {
            sim_mat_sum[it] += sim_mat[it][itv];
            //cout<<sim_mat_sum[it]<<endl;
        }
    }

}

void dataset::cost_normalize(){
    float rate = rating.size()/ave_num;
    float rate_coef,year_coef;
    float rate_sum=0.0;
    float year_sum=0.0;
    float rate_base = 10.0,year_base = 1985.0;
    for(int i = 0;i < rating.size();i++){
        rate_sum += rating.at(i);
        year_sum += abs(year.at(i)-year_base);
    }
    rate_sum = rating.size()*rate_base - rate_sum;
    //year_sum = rating.size()*year_base;
    rate_coef = rate/rate_sum;
    year_coef = rate/year_sum;

    for(int i=0;i<d;i++)
        costs.push_back(vector<double>());
    //normalization
    //*
    for(int i = 0;i < rating.size();i++){

        rate_costs.push_back(fabs(rating.at(i)-rate_base)*rate_coef);
        year_costs.push_back(fabs(year.at(i)-year_base)*year_coef);
        /*********not normalized cost*************/
//        rate_costs.push_back(fabs(rating.at(i)-rate_base)*100);
//        year_costs.push_back(fabs(year.at(i)-year_base)*100);
//        rate_costs.push_back(fabs(rating.at(i)-rate_base));
//        year_costs.push_back(fabs(year.at(i)-year_base));

       if(rate_costs[i]<1e-6)
           rate_costs[i]=1e-6;
        if(year_costs[i]<1e-6)
            year_costs[i]=1e-6;

       costs[0].push_back(rate_costs[i]);
       if(d>=2) costs[1].push_back(year_costs[i]);

       modular_cost.push_back(rate_costs[i]);

       //cout<<rate_costs[i]<<endl;
    }
    //*/
//*
    for(int i=0;i<rate_costs.size();i++)
    {
        double temp=min(rate_costs[i],year_costs[i]);
        if((temp<min_cost))
            min_cost=temp;

        double temp2=max(rate_costs[i],year_costs[i]);
        if(temp2>max_cost)
            max_cost=temp2;
    }
    //*/
}


void dataset::readfile(string filename){
    //dataset* data = new dataset();
    ifstream infofile;
    infofile.open(filename.c_str(), ios::in);
    int count = 0;
    char delim = '|';
    if (infofile.is_open()) {
        while (!infofile.eof()) {
            std::string line;
            getline(infofile, line);
            if (line.empty())
                continue;
            std::string::size_type pos = line.find_first_of(delim);
            v.insert(count++);
            int prevpos = 0;
            int length = line.length();
            // id
            string str = line.substr(prevpos, pos - prevpos);
            index.push_back(atoi(str.c_str()));
            // titles
            prevpos = line.find_first_not_of(delim, pos);
            pos = line.find_first_of(delim, prevpos);
            str = line.substr(prevpos, pos - prevpos);
            titles.push_back(str);
            // genres
            prevpos = line.find_first_not_of(delim, pos);
            pos = line.find_first_of(delim, prevpos);
            str = line.substr(prevpos, pos - prevpos);
            genres.push_back(genre_extractor(str,' '));
            // lambda
            lambda = 0.2;
            //features
            prevpos = line.find_first_not_of(delim, pos);
            pos = line.find_first_of(delim, prevpos);
            str = line.substr(prevpos, pos - prevpos);
            matrix.push_back(feature_extractor(str,' '));

            // year
            prevpos = line.find_first_not_of(delim, pos);
            pos = line.find_first_of(delim, prevpos);
            str = line.substr(prevpos, pos - prevpos);
            year.push_back(atoi(str.c_str()));

            // rating
            prevpos = line.find_first_not_of(delim, pos);
            pos = line.find_first_of(delim, prevpos);
            str = line.substr(prevpos, pos - prevpos);
            rating.push_back(atof(str.c_str()));
        }
    }
}

void dataset::sim_mat_generator(){ //计算Mij
    //e^-lambda*dist(v_i,v_j)
    //*
    sim_mat = vector<vector<float>>(index.size(),vector<float>(index.size(),0));
    for(int i = 0;i < index.size();i++)
    {
        for(int j = i;j < index.size();j++)
        {
            float distance = 0;
            float sim_dis = 0;
            for(int k = 0;k < matrix[0].size();k++)
            {
                distance += pow((matrix[i][k]-matrix[j][k]),2);
            }
            distance = sqrt(distance);
            sim_dis = exp(-1.0*lambda*distance);
            //cout<<sim_dis<<endl;
            sim_mat[i][j] = sim_dis;
            sim_mat[j][i] = sim_dis;
        }
    }
    //*/
    //inner product
    /*
    sim_mat = vector<vector<float>>(index.size(),vector<float>(index.size(),0));
    for(int i = 0;i < index.size();i++)
    {
        for(int j = i;j < index.size();j++)
        {
            float product = 0;
            //float sim_dis = 0;
            for(int k = 0;k < matrix[0].size();k++)
            {
                product+=(matrix[i][k]*matrix[j][k]);
                //distance += pow((matrix[i][k]-matrix[j][k]),2);
            }
            //distance = sqrt(distance);
            //sim_dis = exp(-1*lambda*distance);
            //cout<<product<<endl;
            sim_mat[i][j] = product;
            sim_mat[j][i] = product;
        }
    }
    //*/
    //consine
    /*
    sim_mat = vector<vector<float>>(index.size(),vector<float>(index.size(),0));
    for(int i = 0;i < index.size();i++)
    {
        for (int j = i; j < index.size(); j++)
        {
            double product = 0.0;
            double module_a=0.0;
            double module_b=0.0;
            for (int k = 0; k < matrix[0].size(); k++)
            {
                product += (matrix[i][k] * matrix[j][k]);
                module_a+=pow(matrix[i][k],2);
                module_b+=pow(matrix[j][k],2);
            }
            module_a=sqrt(module_a);
            module_b=sqrt(module_b);
            sim_mat[i][j] = product/(module_a*module_b);
            sim_mat[j][i] = product/(module_a*module_b);
        }
    }
    //*/
    //dist(v_i,v_j)
    /*
    sim_mat = vector<vector<float>>(index.size(),vector<float>(index.size(),0));
    for(int i = 0;i < index.size();i++)
    {
        for(int j = i;j < index.size();j++)
        {
            float distance = 0;
            for(int k = 0;k < matrix[0].size();k++)
            {
                distance += pow((matrix[i][k]-matrix[j][k]),2);
            }
            distance = sqrt(distance);
            sim_mat[i][j] = distance;
            sim_mat[j][i] = distance;
        }
    }
*/
}

vector<float> dataset::feature_extractor(string features, char delim){
    vector<float> str;
    std::string::size_type pos = features.find_first_of(delim);
    int prevpos = 0;
    while(string::npos != pos || string::npos != prevpos){
        str.push_back(atof(features.substr(prevpos,pos-prevpos).c_str()));
        prevpos = features.find_first_not_of(delim, pos);
        pos = features.find_first_of(delim, prevpos);
    }
    return str;
}

vector<string> dataset::genre_extractor(string genres, char delim) {
    vector<string> str;
    std::string::size_type pos = genres.find_first_of(delim);
    int prevpos = 0;
    while (string::npos != pos || string::npos != prevpos) {
        str.push_back(genres.substr(prevpos, pos - prevpos).c_str());
        prevpos = genres.find_first_not_of(delim, pos);
        pos = genres.find_first_of(delim, prevpos);
    }
    return str;
}

double g_u(const int &node)
{
    double sum_value=sim_mat_sum[node];
    return (sum_value-lambda_f*sim_mat[node][node])/node_num;
}
double f_u(const int &node)
{
    double sum_value=sim_mat_sum[node];
    double cost=modular_cost[node];

    return (sum_value-lambda_f*sim_mat[node][node])/node_num-cost;
}
double f_u_distorted(const int &node,const double &ratio)
{
    double sum_value=sim_mat_sum[node];
    double cost=modular_cost[node];

    return (sum_value-lambda_f*sim_mat[node][node])/node_num-ratio*cost;
}

/*bool node_budget_feasible(const int &node,double Budget,dataset* data)
{
    *//**********check single node Budget constraint**********//*
    bool feasible=true;
    for(int r=0;r<d;r++)
    {
        if(data->costs[r][node]>Budget)
        {
            feasible=false;
            break;
        }
    }
    return feasible;
}*/
class S_class
{
public:
    S_class()
    {
        S_cost=0.0;
        S_revenue=0.0;
        for(int iter=0;iter<k;iter++)
            genres_sum.push_back(0);
        max_marginal=-999999999.0;
        max_element=-1;
    }
    double S_cost;
    double S_revenue;
    vector<int> genres_sum;
    unordered_set<int> solution;
    //used for offline algorithm
    double max_marginal;
    int max_element;
    /***********new method***************/
    void clear()
    {
        S_cost=0.0;
        for(int iter=0;iter<k;iter++)
            genres_sum[iter]=0;

        solution.clear();
        S_revenue = 0.0;
    }
    void replace_with_singleton(const double &marginal,const int &node)
    {
        clear();
        add_element(marginal,node);
    }
    void add_element(const double &marginal,const int &node)
    {
        solution.emplace(node);
        S_cost+=modular_cost[node];
        S_revenue+=marginal;
        for (int i = 0; i < genres_sum.size(); i++)//S+e\in I
        {
            genres_sum[i] += genres_index[node][i];
        }
    }
    void add_element(const int &node)
    {
        solution.emplace(node);
        S_cost+=modular_cost[node];
        S_revenue+=marginal_gain(node);
        for (int i = 0; i < genres_sum.size(); i++)//S+e\in I
        {
            genres_sum[i] += genres_index[node][i];
        }
    }
    void add_batched_element(vector<int>::iterator seq_begin,vector<int>::iterator seq_end)
    {
        for(vector<int>::iterator iter=seq_begin;iter!=seq_end;++iter)
        {
            // category_sum[lable[(*iter)]]++;
            // Set.push_back((*iter));
            add_element(0.0,*iter);
        }
        S_revenue=f_S();
    }
    double f_S()
    {
        double sum_cost=0.0;
        double sum_value=0.0;
        float M1 = 0,M2 = 0;
        for(int it : solution)
        {
            M1 += sim_mat_sum[it];
            for(int its : solution)
                M2 += sim_mat[it][its];

            //sum cost
            sum_cost+=modular_cost[it];
        }
        sum_value = M1-lambda_f*M2;
        return sum_value/node_num-sum_cost;
    }
    bool is_feasible(int e,int max_movie)
    {
        if(solution.count(e))
            return false;

        if(solution.size()>=max_movie)
            return false;
        bool flag=true;
        for(int iter=0;iter<genres_index[e].size();iter++)
        {
            if(genres_index[e][iter]+genres_sum[iter]>genres_limit)
            {
                flag=false;
                break;
            }
        }
        return flag;
    }

    //used for non-regularized algorithm
    double marginal_gain(int e)
    {
        //Set.push_back(node);
        //selected[node]=1;
        //double f_S_and_u=f_S();
        //Set.pop_back();
        //selected[node]=0;
        float M1=0,M2=0;
        M1 += sim_mat_sum[e];
        for(int it : solution)
        {
            if(e==it) return 0;
            M2 += sim_mat[e][it];
            M2 += sim_mat[it][e];
        }
        M2 += sim_mat[e][e];

        double cost=modular_cost[e];
        return (M1-lambda_f*M2)/node_num-cost;
    }
    double marginal_gain_scaled(int e,const double &ratio)
    {
        //Set.push_back(node);
        //selected[node]=1;
        //double f_S_and_u=f_S();
        //Set.pop_back();
        //selected[node]=0;
        float M1=0,M2=0;
        M1 += sim_mat_sum[e];
        for(int it : solution)
        {
            if(e==it) return 0;
            M2 += sim_mat[e][it];
            M2 += sim_mat[it][e];
        }
        M2 += sim_mat[e][e];

        double cost=modular_cost[e];
        return (M1-lambda_f*M2)/node_num-ratio*cost;
    }

    //used for regularized algorithm, only return submodular marginal gain
    double marginal_g(int e)
    {
        //Set.push_back(node);
        //selected[node]=1;
        //double f_S_and_u=f_S();
        //Set.pop_back();
        //selected[node]=0;
        float M1=0,M2=0;
        M1 += sim_mat_sum[e];
        for(int it : solution)
        {
            if(e==it) return 0;
            M2 += sim_mat[e][it];
            M2 += sim_mat[it][e];
        }
        M2 += sim_mat[e][e];

        return (M1-lambda_f*M2)/node_num;
    }
};
class S_array
{
public:
    vector<S_class> sarray;
    double gamma;
    //S_array(){}
    S_array(const double g,const int h)
    {
        gamma=g;
        for(int i=0;i<h;i++)
            sarray.push_back(S_class());
    }
};
class Result
{
public:
    Result(){}
    Result(double rev,double cos,int siz,long long int ora)
    {
        revenue=rev;
        cost=cos;
        size=siz;
        oracle=ora;
    }
    Result(double rev,double cos,int siz,long long int ora,long long int rou)
    {
        revenue=rev;
        cost=cos;
        size=siz;
        oracle=ora;
        round=rou;
    }

    double revenue;
    long long int oracle;
    long long int round;
    double cost;
    int size;
    long long int time=0;
    long long int max_query=0;
};

#endif //STREAMING_ALGORITHM_UTILITY_FUNCTIONS_H
