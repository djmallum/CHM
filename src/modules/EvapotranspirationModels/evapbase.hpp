class evapT_base
{
public:
    virtual ~evapT_base() = default;

    virtual void CalcEvapT(param_base& baseparam, var_base& basevar, model_output& output) = 0;
    
    double delta(double& t);
    double gamma(double& P_atm, double& t);
}

struct param_base
{

}

struct var_base
{

}

struct model_output
{
    double ET;
}
